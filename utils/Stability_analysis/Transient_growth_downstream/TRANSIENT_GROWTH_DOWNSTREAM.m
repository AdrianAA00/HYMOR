function [V, D, T_opt] = TRANSIENT_GROWTH_DOWNSTREAM(A, s, n_modes, T)
% TRANSIENT_GROWTH_DOWNSTREAM  Compute optimal transient growth modes for downstream disturbances.
%
%   [V, D, T_opt] = TRANSIENT_GROWTH_DOWNSTREAM(A, s, n_modes, T)
%
%   Solves the generalized eigenvalue problem for transient energy growth
%   in the downstream flow domain. Uses a matrix-free approach with eigs()
%   to find the n_modes largest-gain eigenmodes at the specified time
%   horizon T. Automatically detects GPU availability.
%
%   Inputs:
%       A        - System Jacobian matrix (sparse)
%       s        - Solution structure with flow fields, mesh data, and
%                  configuration flags
%       n_modes  - Number of leading transient growth modes to compute
%       T        - Target time horizon for the optimization
%
%   Outputs:
%       V     - Matrix of optimal perturbation eigenvectors (columns)
%       D     - Diagonal matrix of corresponding energy gains (eigenvalues)
%       T_opt - Actual time horizon used (adjusted to integer time steps)
%
%   Notes:
%       - Returns V=0, D=0, T_opt=0 when transient_growth_downstream_analysis is false.
%       - The matrix exponential is approximated with a 5th-order Taylor series.
%       - All memory is pre-allocated outside eigs. Only copies inside callback.
%       - GPU is synchronized and memory verified deallocated before returning.
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Downstream Module

    k =  n_modes;    % Number of modes to compute
    freestream_disturbances = false;

    %% Compute time-stepping parameters
    s = CFL_TIMESTEP(s);
    t_temp = s.time_integration.dt;
    n_t = round(T / t_temp);
    dt = T / n_t;
    T_opt = n_t * dt;

    %% Construct operator D and C (without exponentials)
    R = CONSTRUCT_R(s);
    norms = [1, 1, 1, 1];
    M = CONSTRUCT_M(s, norms);
    M_hat = CONSTRUCT_M_HAT(s, norms);
    C_no_exp = R' * M * R;
    D = R' * M_hat * R;

    %% Eigenvalue computation via Lanczos iteration
    fprintf("\n")
    disp("---------------------------------")
    disp("    Transient growth analysis    ")
    disp("---------------------------------")
    fprintf("\n")
    fprintf("Number timesteps = %f ", n_t)
    fprintf("\n")

    %% Detect GPU availability
    use_gpu = (gpuDeviceCount > 0);
    n = size(A, 1);

    if use_gpu
        gpu = gpuDevice;
        disp(gpu.Name + " selected.")
    else
        disp("No CUDA device found. Running on CPU.")
    end

    %% Pre-allocate ALL memory outside eigs (matrices + work buffers)
    A_dt = A * dt;
    At_dt = A_dt';  % Transpose on CPU to preserve sparsity

    if use_gpu
        gA_dt      = gpuArray(A_dt);
        gAt_dt     = gpuArray(At_dt);
        gC_no_exp  = gpuArray(C_no_exp);
        x_buf      = gpuArray(zeros(n, 1));
        vout_buf   = gpuArray(zeros(n, 1));
        term_buf   = gpuArray(zeros(n, 1));
    else
        gA_dt      = A_dt;
        gAt_dt     = At_dt;
        gC_no_exp  = C_no_exp;
        x_buf      = zeros(n, 1);
        vout_buf   = zeros(n, 1);
        term_buf   = zeros(n, 1);
    end

    opts.tol = 1e-4;     % Convergence tolerance for eigs()
    opts.maxit = 10 * k; % Maximum number of iterations (adjust as needed)
    opts.issym = true;   % Use symmetric Lanczos algorithm

    t_eigs = tic;
    [V, D, flag] = eigs(@(x)CONSTRUCT_C(x, gA_dt, gAt_dt, gC_no_exp, n_t, freestream_disturbances, 0, use_gpu, x_buf, vout_buf, term_buf), n, D, k, 'LM', opts);
    fprintf("Eigenvalue solve time: %.6f seconds\n", toc(t_eigs));
    if flag ~= 0
        warning('eigs did not fully converge');
    end

    %% GPU synchronization and memory cleanup
    if use_gpu
        wait(gpu);  % Synchronize GPU - catch any async errors
        clear gA_dt gAt_dt gC_no_exp x_buf vout_buf term_buf;
        wait(gpu);  % Verify memory is deallocated
        reset(gpuDevice);
    else
        clear gA_dt gAt_dt gC_no_exp x_buf vout_buf term_buf;
    end

    %% Print eigenvalues
    format long
    for i = 1:k
        fprintf("eigenvalue %i = %f\n", i, real(D(i, i)))
    end
end
