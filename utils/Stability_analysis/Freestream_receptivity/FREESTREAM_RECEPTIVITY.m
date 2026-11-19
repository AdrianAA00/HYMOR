function [V, D, T_opt] = FREESTREAM_RECEPTIVITY(A_, s, n_modes, T, w_infty)
% FREESTREAM_RECEPTIVITY  Compute optimal transient growth modes including freestream disturbances.
%
%   [V, D, T_opt] = FREESTREAM_RECEPTIVITY(A_, s, n_modes, T, w_infty)
%
%   Solves the generalized eigenvalue problem for transient energy growth
%   in the coupled downstream-freestream system. Uses a matrix-free
%   approach with eigs() to find the n_modes largest-gain eigenmodes at
%   the specified time horizon T. Runs on the device chosen by USE_GPU.
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
%       - Largest gains by BLOCK_LANCZOS on the GPU (blocks of 4 vectors per
%         propagation) and by eigs (Lanczos) on the CPU.
%       - Downstream DoFs are zeroed at t=0 so that only freestream forcing
%         drives the response (freestream_disturbances = true).
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Freestream Module

    k = n_modes;        % Maximum of modes to compute

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

    %% Device of the time stepper
    use_gpu = USE_GPU(s, A_);
    N_set_0 = numel(PERTURBATION_FIELDS(s)) * s.mesh.Nchi * s.mesh.Neta;   % downstream DoFs, zero at t = 0
    prop = TAYLOR_PROPAGATOR(A_, dt, use_gpu, true);
    C = prop.matrix(C_no_exp_);
    apply = @(x) CONSTRUCT_C(x, prop, C, n_t, N_set_0);

    %% Largest gains of the generalized eigenvalue problem
    % On the GPU, block Lanczos propagates BLOCK_SIZE vectors at once: one
    % sparse product per Taylor order for the whole block.
    t_eigs = tic;
    if use_gpu
        BLOCK_SIZE = 4;
        [V, gains, flag] = BLOCK_LANCZOS(apply, D_, k, BLOCK_SIZE, 1e-4, 25 * BLOCK_SIZE, false);
        D = diag(gains);
    else
        opts.tol = 1e-4;      % Convergence tolerance for eigs()
        opts.maxit = 10 * k;  % Maximum number of iterations (adjust as needed)
        opts.issym = true;    % Use symmetric Lanczos algorithm
        opts.isreal = false;  % Use non real Lanczos algorithm
        [V, D, flag] = eigs(apply, size(A_, 1), D_, k, 'LM', opts);
    end
    fprintf("Eigenvalue solve time: %.6f seconds\n", toc(t_eigs))
    if flag ~= 0
        warning('eigs did not fully converge');
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
