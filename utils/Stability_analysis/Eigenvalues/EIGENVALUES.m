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
%         converted to continuous-time eigenvalues via logarithmic mapping.
%       - The CFL time step is recomputed before GPU analysis.
%       - Detects CUDA availability: uses GPU if available, CPU otherwise.
%       - All memory is pre-allocated outside eigs and only copies are done
%         inside the eigs callback. GPU is synchronized before cleanup.
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
%   with a Taylor-expansion time stepper. Automatically detects GPU
%   availability and falls back to CPU if no CUDA device is found.
%   All memory is pre-allocated outside eigs to avoid repeated allocation.

    fprintf("\n")
    disp("----------------------------------------------------")
    disp(" Eigenvalue solver: Arnoldi + exponential transform ")
    disp("----------------------------------------------------")

    %% Detect GPU availability
    use_gpu = (gpuDeviceCount > 0);

    if use_gpu
        gpu = gpuDevice;
        disp(gpu.Name + " selected.")
    else
        disp("No CUDA device found. Running on CPU.")
    end

    %% Time-stepping parameters
    T_ref = 0.5; % Time to evolve
    s = CFL_TIMESTEP(s);
    dt = s.time_integration.dt;
    timesteps = round(T_ref/dt);
    fprintf("Number of timesteps: %f\n", timesteps)
    n = size(A, 1);
    k = n_modes;

    %% Pre-allocate all memory OUTSIDE eigs
    if use_gpu
        Adt = gpuArray(A * dt);
        x_buf    = gpuArray(zeros(n, 1));
        vout_buf = gpuArray(zeros(n, 1));
        term_buf = gpuArray(zeros(n, 1));
    else
        Adt = A * dt;
        x_buf    = zeros(n, 1);
        vout_buf = zeros(n, 1);
        term_buf = zeros(n, 1);
    end

    %% Arnoldi iteration via eigs (no memory allocation inside callback)
    opts.tol = 1e-4;
    opts.maxit = 40 * k;
    opts.issym = false;
    t_eigs = tic;
    [V, D, flag] = eigs(@(x) ARNOLDI_POWER_GPU(x, Adt, timesteps, use_gpu, x_buf, vout_buf, term_buf), n, k, 'LM', opts);
    fprintf("Eigenvalue solve time: %.6f seconds\n", toc(t_eigs));
    if flag ~= 0
        warning('eigs did not fully converge');
    end

    %% GPU synchronization and memory cleanup
    if use_gpu
        wait(gpu);  % Synchronize GPU - catch any async errors
        clear Adt x_buf vout_buf term_buf;
        wait(gpu);  % Verify memory is deallocated
        reset(gpuDevice);
    else
        clear Adt x_buf vout_buf term_buf;
    end

    %% Convert discrete-time eigenvalues to continuous-time
    format long

    for i = 1:k
        D(i, i) = real(log(D(i, i)) / dt / timesteps) + 1j * imag(log(D(i, i)) / dt / timesteps);
        fprintf("eigenvalue ")
        fprintf("%i", i)
        fprintf(" = ")
        fprintf("%f + %fi", real(D(i, i)), imag(D(i, i)))
        fprintf("\n")
    end
end
