function x = ARNOLDI_POWER_GPU(x, Adt, timesteps, use_gpu, x_buf, vout_buf, term_buf)
% ARNOLDI_POWER_GPU  Matrix exponential action via Taylor expansion.
%
%   x = ARNOLDI_POWER_GPU(x, Adt, timesteps, use_gpu, x_buf, vout_buf, term_buf)
%   computes the action of the matrix exponential exp(A*dt*timesteps) on
%   vector x using a fixed-order Taylor expansion. Works on both GPU and CPU.
%
%   All buffers (x_buf, vout_buf, term_buf) are pre-allocated OUTSIDE the
%   eigs loop to avoid repeated memory allocation. Inside this function,
%   only memory copies are performed (no allocation).
%
%   Inputs:
%       x         - Input vector (CPU array, length n).
%       Adt       - System matrix A*dt (gpuArray or CPU sparse).
%       timesteps - Number of time steps to advance (scalar).
%       use_gpu   - Logical flag: true for GPU path, false for CPU path.
%       x_buf     - Pre-allocated buffer for input vector (gpuArray or CPU).
%       vout_buf  - Pre-allocated buffer for output accumulation (gpuArray or CPU).
%       term_buf  - Pre-allocated buffer for Taylor term (gpuArray or CPU).
%
%   Outputs:
%       x - Result vector after applying the matrix exponential action
%           (CPU array, length n).
%
%   Notes:
%       - Uses a 5th-order Taylor expansion of exp(A*dt) per time step.
%       - When use_gpu=true, the vector is copied to pre-allocated GPU
%         buffers, processed, and gathered back. No GPU allocation occurs.
%       - When use_gpu=false, all operations are performed on CPU.
%
% Part of: Hypersonics Stability MATLAB Solver - Stability Analysis / Eigenvalues Module

    if use_gpu
        x_buf(:) = x;
        [x_buf] = POWER_METHOD(x_buf, Adt, timesteps, vout_buf, term_buf);
        x = gather(x_buf);
    else
        x_buf(:) = x;
        [x_buf] = POWER_METHOD(x_buf, Adt, timesteps, vout_buf, term_buf);
        x = x_buf;
    end
end

%% ========================================================================
%  Taylor Expansion Time Stepper (uses pre-allocated buffers)
%  ========================================================================
function [vout] = POWER_METHOD(v, A, timesteps, vout, term)
% POWER_METHOD  Advance a vector through multiple time steps using a
%   truncated Taylor series expansion of the matrix exponential.
%   Uses pre-allocated buffers vout and term to avoid memory allocation.

    order = 5; % Taylor expansion order
    vout(:) = v;

    for i = 1:timesteps
        factorial_val = 1;
        for j = 1:order
            term(:) = A * v;
            v(:) = term;
            factorial_val = j * factorial_val;
            vout(:) = vout + v / factorial_val;
        end
        v(:) = vout;
    end
end
