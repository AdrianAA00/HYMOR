function [V, D] = EIGENVALUES(A, s, n_modes)
% EIGENVALUES  Compute leading eigenvalues and eigenvectors of the linearized operator.
%
%   [V, D] = EIGENVALUES(A, s, n_modes) dispatches to a CPU or GPU
%   eigenvalue solver based on s.stability_analysis.eigenvalue_solver. Available
%   methods are:
%       "CPU_LU"       - MATLAB eigs with shift-invert (LU factorization).
%       "GPU_TIMESTEPPER_ARNOLDI"  - Arnoldi power method via eigs (GPU or CPU).
%
%   Inputs:
%       A        - Linearized Jacobian matrix (sparse, n x n).
%       s        - Solution structure with solver parameters.
%       n_modes  - Number of eigenvalues/eigenvectors to compute.
%
%   Outputs:
%       V - Matrix of eigenvectors (columns), size n x n_modes.
%       D - Diagonal matrix of eigenvalues, size n_modes x n_modes.
%           Returns V=0, D=0 if linearization is disabled.
%
%   Notes:
%       - For GPU_ARNOLDI, eigenvalues of the time-stepping operator are
%         converted to continuous-time eigenvalues via logarithmic mapping,
%         on the branch given by the Rayleigh quotient of each eigenvector.
%       - The CFL time step is recomputed before GPU analysis.
%       - The time stepper runs on the device chosen by USE_GPU.
%
% Part of: Hypersonics Stability MATLAB Solver - Stability Analysis / Eigenvalues Module

    switch s.stability_analysis.eigenvalue_solver
        case "CPU_LU"
            [V, D] = EIGENVALUES_CPU_LU(A, s, n_modes);
        case "GPU_TIMESTEPPER_ARNOLDI"
            [V, D] = EIGENVALUES_GPU_ARNOLDI(A, s, n_modes);
    end
end

%% ========================================================================
%  CPU Shift-Invert Eigenvalue Solver
%  ========================================================================
function [V, D] = EIGENVALUES_CPU_LU(A, s, n_modes)
% EIGENVALUES_CPU_LU  Compute eigenvalues using MATLAB eigs with
%   shift-invert around a specified sigma in the complex plane.

    fprintf("\n")
    disp("------------------------")
    disp("      eigs MATLAB       ")
    disp("------------------------")

    %% Solver options
    opts.tol = 1e-6;
    opts.maxit = 100;

    %% Frequency sweep parameters
    sigma_imag_min = 0.1;
    sigma_imag_max = 0.50;
    sigma_imag_step = 0.1;
    sigma_real = 0.02;

    N = round((sigma_imag_max - sigma_imag_min) / sigma_imag_step);

    if s.stability_analysis.perturb_shock
        V = zeros(s.mesh.Nchi * s.mesh.Neta * 4 + s.mesh.Nchi, N);
        D = zeros(N, 1);
    else
        V = zeros(s.mesh.Nchi * s.mesh.Neta * 4, N);
        D = zeros(N, 1);
    end

    %% Targeted eigenvalue search
    N = 4;
    tic
    [V, D] = eigs(A, 4, -0.606629 + 0.462148i, opts);
    toc

    %% Display results
    for i = 1:N
        fprintf("eigenvalue ")
        fprintf("%i", i)
        fprintf(" = ")
        fprintf("%f + %fi", real(D(i)), imag(D(i)))
        fprintf("\n")
    end
end

%% ========================================================================
%  GPU/CPU Arnoldi Power Method Eigenvalue Solver
%  ========================================================================
function [V, D] = EIGENVALUES_GPU_ARNOLDI(A, s, n_modes)
% EIGENVALUES_GPU_ARNOLDI  Compute eigenvalues using Arnoldi iteration
%   on the Taylor time stepper exp(A T_ref), on the device chosen by
%   USE_GPU.

    fprintf("\n")
    disp("----------------------------------------------------")
    disp(" Eigenvalue solver: Arnoldi + exponential transform ")
    disp("----------------------------------------------------")

    %% Device of the time stepper
    use_gpu = USE_GPU(s, A);

    %% Time-stepping parameters
    T_ref = 0.5; % Time to evolve
    s = CFL_TIMESTEP(s);
    dt = s.time_integration.dt;
    timesteps = round(T_ref/dt);
    fprintf("Number of timesteps: %f\n", timesteps)
    prop = TAYLOR_PROPAGATOR(A, dt, use_gpu, false);

    %% Arnoldi iteration via eigs
    opts.tol = 1e-4;
    opts.maxit = 40 * n_modes;
    opts.issym = false;
    t_eigs = tic;
    [V, D, flag] = eigs(@(x) ARNOLDI_POWER_GPU(x, prop, timesteps), size(A, 1), n_modes, 'LM', opts);
    fprintf("Eigenvalue solve time: %.6f seconds\n", toc(t_eigs));
    if flag ~= 0
        warning('eigs did not fully converge');
    end

    %% Convert discrete-time eigenvalues to continuous-time
    % log(mu) / T is known up to multiples of 2 pi i / T: frequencies above
    % pi / T alias. The Rayleigh quotient v' A v / v' v of each eigenvector
    % picks the branch (it is accurate to far better than pi / T).
    T = dt * timesteps;
    lambda = log(diag(D)) / T;
    rayleigh = (sum(conj(V) .* (A * V), 1) ./ sum(abs(V).^2, 1)).';
    lambda = lambda + 2i * pi / T * round((imag(rayleigh) - imag(lambda)) * T / (2 * pi));
    D = diag(lambda);

    format long
    for i = 1:n_modes
        fprintf("eigenvalue %i = %f + %fi\n", i, real(lambda(i)), imag(lambda(i)))
    end
end
