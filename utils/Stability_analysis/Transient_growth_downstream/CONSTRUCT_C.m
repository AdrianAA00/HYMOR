function x = CONSTRUCT_C(x, A_dt, At_dt, C_no_exp, n_t, freestream_disturbances, N_set_0, use_gpu, x_buf, vout_buf, term_buf)
% CONSTRUCT_C  Apply the transient growth operator C to a state vector.
%
%   x = CONSTRUCT_C(x, A_dt, At_dt, C_no_exp, n_t, freestream_disturbances,
%       N_set_0, use_gpu, x_buf, vout_buf, term_buf)
%
%   Computes C*x = exp(A'*t) * R'*M*R * exp(A*t) * x using a Taylor-series
%   approximation for the matrix exponential. Works on both GPU and CPU.
%
%   All buffers (x_buf, vout_buf, term_buf) are pre-allocated OUTSIDE the
%   eigs loop. Inside this function, only memory copies are performed.
%
%   Inputs:
%       x                       - State vector (column vector)
%       A_dt                    - System Jacobian matrix A*dt (gpuArray or CPU sparse)
%       At_dt                   - Transpose (A*dt)' (gpuArray or CPU sparse)
%       C_no_exp                - Energy weight matrix R'*M*R (gpuArray or CPU sparse)
%       n_t                     - Number of time steps (scalar)
%       freestream_disturbances - Flag: true to zero out downstream DoFs at t=0
%       N_set_0                 - Number of downstream DoFs to zero when freestream mode is active
%       use_gpu                 - Logical flag: true for GPU path, false for CPU path
%       x_buf                   - Pre-allocated buffer for input vector
%       vout_buf                - Pre-allocated buffer for output accumulation
%       term_buf                - Pre-allocated buffer for Taylor term
%
%   Outputs:
%       x - Result of applying the operator C to the input vector
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Downstream Module

    %% Zero downstream disturbances if freestream mode
    if freestream_disturbances
        x(1:N_set_0, 1) = 0;
    end

    if use_gpu
        %% GPU path: copy into pre-allocated buffer (no allocation)
        x_buf(:) = x;
        [x_buf] = EXPONENTIAL(x_buf, A_dt, At_dt, C_no_exp, n_t, vout_buf, term_buf);
        x = gather(x_buf);
    else
        %% CPU path: copy into pre-allocated buffer (no allocation)
        x_buf(:) = x;
        [x_buf] = EXPONENTIAL(x_buf, A_dt, At_dt, C_no_exp, n_t, vout_buf, term_buf);
        x = x_buf;
    end
end


function [vout] = EXPONENTIAL(v, A, At, C_no_exp, timesteps, vout, term)
% EXPONENTIAL  Compute exp(A'*t) * C_no_exp * exp(A*t) * v via Taylor series.
%   Uses pre-allocated buffers vout and term to avoid memory allocation.
%
%   Inputs:
%       v         - Input state vector
%       A         - System matrix (A*dt)
%       At        - Transpose of system matrix (A*dt)'
%       C_no_exp  - Energy weight matrix (R'*M*R)
%       timesteps - Number of sub-steps for the exponential
%       vout      - Pre-allocated buffer for output accumulation
%       term      - Pre-allocated buffer for Taylor term
%
%   Outputs:
%       vout - Result vector after applying the full operator

    order = 5; % Fixed order for Taylor expansion

    vout(:) = v;
    for i = 1:timesteps
        term(:) = vout;
        for j = 1:order
            term(:) = (A * term) / j; % Iteratively compute next Taylor term
            vout(:) = vout + term;
        end
    end

    %% Apply energy weight: vout = R'*M*R * exp(A*t) * v
    v(:) = C_no_exp * vout;
    vout(:) = v;

    %% Adjoint propagation: vout = exp(A'*t) * R'*M*R * exp(A*t) * v
    for i = 1:timesteps
        term(:) = vout;
        for j = 1:order
            term(:) = (At * term) / j; % Iteratively compute next Taylor term
            vout(:) = vout + term;
        end
    end
end
