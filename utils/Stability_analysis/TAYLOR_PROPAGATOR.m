function prop = TAYLOR_PROPAGATOR(A, dt, use_gpu, adjoint)
% TAYLOR_PROPAGATOR  Device-side operators of the Taylor time stepper.
%
%   prop = TAYLOR_PROPAGATOR(A, dt, use_gpu, adjoint)
%
%   Renumbers the unknowns with reverse Cuthill-McKee (symrcm) and stores
%   A dt in that order on the chosen device (MATVEC_HANDLE). HYMOR numbers
%   the unknowns variable by variable, so a row of A reads four distant
%   parts of the vectors; in RCM order the reads are local. This matters
%   for blocks of vectors that do not fit in the GPU cache: a 16-vector
%   Taylor step is 1.45x faster on the 480 x 160 blunt cone (single
%   vectors are unchanged). States move into and out of that order only
%   through prop.load and prop.unload.
%
%   Inputs:
%       A       - Sparse system matrix (n x n).
%       dt      - Time step.
%       use_gpu - true to run on the GPU (from USE_GPU).
%       adjoint - true to also build the adjoint step (A dt)'.
%
%   Outputs:
%       prop - Struct with
%              .step         - handle, (A dt) * x in RCM order (TAYLOR_STEPS)
%              .step_adjoint - handle, (A dt)' * x (when adjoint)
%              .load         - handle, host vectors (n x b) -> device, RCM order
%              .unload       - handle, the inverse of load
%              .matrix       - handle, M -> MATVEC_HANDLE of M in RCM order
%                              (M may stack n x n blocks vertically)
%
% Part of: Hypersonics Stability MATLAB Solver - Stability Analysis Module

    p = symrcm(spones(A) + spones(A'));
    Ap = A(p, p) * dt;
    prop.step = MATVEC_HANDLE(Ap, use_gpu);
    if adjoint
        prop.step_adjoint = MATVEC_HANDLE(Ap', use_gpu);
    end
    prop.matrix = @(M) MATVEC_HANDLE(M(STACKED_ROWS(p, size(M, 1)), p), use_gpu);
    if use_gpu
        prop.load = @(x) gpuArray(x(p, :));
    else
        prop.load = @(x) x(p, :);
    end
    prop.unload = @(y) UNPERMUTE(gather(y), p);
end


function rows = STACKED_ROWS(p, n_rows)
% STACKED_ROWS  RCM order applied to each n x n block of a vertical stack.
    rows = reshape(p(:) + numel(p) * (0:n_rows / numel(p) - 1), [], 1);
end


function x = UNPERMUTE(y, p)
% UNPERMUTE  x(p, :) = y.
    x = zeros(size(y), 'like', y);
    x(p, :) = y;
end
