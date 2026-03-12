function [V, D, T_opt] = FREESTREAM_RECEPTIVITY(A_, s, n_modes, T, w_infty)
% FREESTREAM_RECEPTIVITY  Compute optimal transient growth modes including freestream disturbances.
%
%   [V, D, T_opt] = FREESTREAM_RECEPTIVITY(A_, s, n_modes, T, w_infty)
%
%   Solves the generalized eigenvalue problem for transient energy growth
%   in the coupled downstream-freestream system. Uses a matrix-free
%   approach with eigs() to find the n_modes largest-gain eigenmodes at
%   the specified time horizon T. Automatically detects GPU availability.
%
%   Inputs:
%       A_                   - Extended system matrix (sparse, from LINEARIZE_L_)
%       s                    - Solution structure with flow fields, mesh data, and
%                              configuration flags
%       n_modes              - Number of leading transient growth modes to compute
%       T                    - Target time horizon for the optimization
%       w_infty              - Vector of freestream disturbance frequencies
%
%   Outputs:
%       V     - Matrix of optimal perturbation eigenvectors (columns)
%       D     - Diagonal matrix of corresponding energy gains (eigenvalues)
%       T_opt - Actual time horizon used (adjusted to integer time steps)
%
%   Notes:
%       - The matrix exponential is approximated with a 5th-order Taylor series.
%       - All memory is pre-allocated outside eigs. Only copies inside callback.
%       - GPU is synchronized and memory verified deallocated before returning.
%       - Downstream DoFs are zeroed at t=0 so that only freestream forcing
%         drives the response (freestream_disturbances = true).
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Freestream Module

    k = n_modes;        % Maximum of modes to compute
    freestream_disturbances = true;

    %% Compute time-stepping parameters
    s = CFL_TIMESTEP(s);
    t_temp = s.time_integration.dt;
    n_t = round(T / t_temp);
    dt = T / n_t;
    T_opt = n_t * dt;

    %% Construct operator D_ and C_ (without exponentials)
    R_ = CONSTRUCT_R_(s, w_infty);
    norms_1 = [1, 1, 1, 1];
    norms_2 = [1, 1, 1, 1];
    M_ = CONSTRUCT_M_(s, norms_1, w_infty);
    scaling_non_temporal = true; % Fixed now
    M_infty_ = CONSTRUCT_M_INFTY_(s, norms_2, T, w_infty, scaling_non_temporal);
    C_no_exp_ = R_' * M_ * R_;
    D_ = R_' * M_infty_ * R_;

    %% Eigenvalue computation
    fprintf("\n")
    disp("-----------------------------------")
    disp("      Freestream receptivity       ")
    disp("-----------------------------------")
    fprintf("\n")
    fprintf("Number timesteps = %f ", n_t)
    fprintf("\n")

    %% Detect GPU availability
    use_gpu = (gpuDeviceCount > 0);
    N_set_0 = 4 * s.mesh.Nchi * s.mesh.Neta;
    n = size(A_, 1);

    if use_gpu
        gpu = gpuDevice;
        disp(gpu.Name + " selected.")
    else
        disp("No CUDA device found. Running on CPU.")
    end

    %% Pre-allocate ALL memory outside eigs (matrices + work buffers)
    A_dt = A_ * dt;
    At_dt = A_dt';  % Transpose on CPU to preserve sparsity

    if use_gpu
        gA_dt      = gpuArray(A_dt);
        gAt_dt     = gpuArray(At_dt);
        gC_no_exp_ = gpuArray(C_no_exp_);
        x_buf      = gpuArray(zeros(n, 1));
        vout_buf   = gpuArray(zeros(n, 1));
        term_buf   = gpuArray(zeros(n, 1));
    else
        gA_dt      = A_dt;
        gAt_dt     = At_dt;
        gC_no_exp_ = C_no_exp_;
        x_buf      = zeros(n, 1);
        vout_buf   = zeros(n, 1);
        term_buf   = zeros(n, 1);
    end

    opts.tol = 1e-4;      % Convergence tolerance for eigs()
    opts.maxit = 10 * k;  % Maximum number of iterations (adjust as needed)
    opts.issym = true;    % Use symmetric Lanczos algorithm

    t_eigs = tic;
    [V, D, flag] = eigs(@(x)CONSTRUCT_C(x, gA_dt, gAt_dt, gC_no_exp_, n_t, freestream_disturbances, N_set_0, use_gpu, x_buf, vout_buf, term_buf), n, D_, k, 'LM', opts);
    fprintf("Eigenvalue solve time: %.6f seconds\n", toc(t_eigs))
    if flag ~= 0
        warning('eigs did not fully converge');
    end


    %% GPU synchronization and memory cleanup
    if use_gpu
        wait(gpu);  % Synchronize GPU - catch any async errors
        clear gA_dt gAt_dt gC_no_exp_ x_buf vout_buf term_buf;
        wait(gpu);  % Verify memory is deallocated
        reset(gpuDevice);
    else
        clear gA_dt gAt_dt gC_no_exp_ x_buf vout_buf term_buf;
    end

    % Numerator gives energy of both real and complex parts of the mode. Each real and imaginary part have approximately same gain,
    % so divide by 2 to get the gain of each mode. Rigorous way would be to compute the numerator gain as real(exp(V*t)')*C*real(exp(V*t)),
    % which is done in LINEAR_INTEGRATION_AND_GAINS, where in depth energy growth analysis are done.
    D = D ./ 2;

    %% Print eigenvalues
    format long
    for i = 1:k
        fprintf("eigenvalue %i = %f\n", i, real(D(i, i)))
    end
end
