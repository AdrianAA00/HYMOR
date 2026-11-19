function [V, D, T_opt] = TRANSIENT_GROWTH_DOWNSTREAM(A, s, n_modes, T)
% TRANSIENT_GROWTH_DOWNSTREAM  Compute optimal transient growth modes for downstream disturbances.
%
%   [V, D, T_opt] = TRANSIENT_GROWTH_DOWNSTREAM(A, s, n_modes, T)
%
%   Solves the generalized eigenvalue problem for transient energy growth
%   in the downstream flow domain. Uses a matrix-free approach with eigs()
%   to find the n_modes largest-gain eigenmodes at the specified time
%   horizon T. Runs on the device chosen by USE_GPU.
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
%       - Largest gains by BLOCK_LANCZOS on the GPU (blocks of 4 vectors per
%         propagation) and by eigs (Lanczos) on the CPU.
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Downstream Module

    k =  n_modes;    % Number of modes to compute

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

    %% Device of the time stepper
    use_gpu = USE_GPU(s, A);
    prop = TAYLOR_PROPAGATOR(A, dt, use_gpu, true);
    C = prop.matrix(C_no_exp);
    apply = @(x) CONSTRUCT_C(x, prop, C, n_t, 0);

    %% Largest gains of the generalized eigenvalue problem
    % On the GPU, block Lanczos propagates BLOCK_SIZE vectors at once: one
    % sparse product per Taylor order for the whole block.
    t_eigs = tic;
    if use_gpu
        BLOCK_SIZE = 4;
        [V, gains, flag] = BLOCK_LANCZOS(apply, D, k, BLOCK_SIZE, 1e-4, 25 * BLOCK_SIZE, true);
        D = diag(gains);
    else
        opts.tol = 1e-4;     % Convergence tolerance for eigs()
        opts.maxit = 10 * k; % Maximum number of iterations (adjust as needed)
        opts.issym = true;   % Use symmetric Lanczos algorithm
        [V, D, flag] = eigs(apply, size(A, 1), D, k, 'LM', opts);
    end
    fprintf("Eigenvalue solve time: %.6f seconds\n", toc(t_eigs));
    if flag ~= 0
        warning('eigs did not fully converge');
    end

    %% Print eigenvalues
    format long
    for i = 1:k
        fprintf("eigenvalue %i = %f\n", i, real(D(i, i)))
    end
end
