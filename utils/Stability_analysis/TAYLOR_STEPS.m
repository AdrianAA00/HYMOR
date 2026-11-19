function x = TAYLOR_STEPS(x, Adt, n_steps)
% TAYLOR_STEPS  Advance x by n_steps order-5 Taylor steps of exp(A dt).
%
%   x = TAYLOR_STEPS(x, Adt, n_steps)
%
%   Each step applies sum_{j=0..5} (A dt)^j / j! in Horner form,
%   w = x + (A dt) w / j for j = 5, ..., 1: five sparse products and no
%   work buffers, so a GPU step is a few kernel launches (latency bound for
%   these operator sizes).
%
%   Inputs:
%       x       - State vector (gpuArray when Adt runs on the GPU).
%       Adt     - Handle computing (A dt) * x, from MATVEC_HANDLE(A * dt, ...)
%                 (or of (A dt)' for the adjoint propagation).
%       n_steps - Number of time steps.
%
%   Outputs:
%       x - exp(A n_steps dt) x to fifth order in dt, on x's device.
%
% Part of: Hypersonics Stability MATLAB Solver - Stability Analysis Module

    order = 5;
    for i = 1:n_steps
        w = x;
        for j = order:-1:1
            w = x + Adt(w) * (1 / j);
        end
        x = w;
    end
end
