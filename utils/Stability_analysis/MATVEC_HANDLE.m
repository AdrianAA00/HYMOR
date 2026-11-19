function mul = MATVEC_HANDLE(M, use_gpu)
% MATVEC_HANDLE  Fastest sparse product x -> M * x on the selected device.
%
%   mul = MATVEC_HANDLE(M, use_gpu)
%
%   On the GPU, M is copied once to a sparse gpuArray. On the CPU, M is
%   stored transposed and the product is formed as (M')' * x: MATLAB
%   multithreads the transposed sparse product but runs M * x on one
%   thread, so this is about ten times faster on a many-core node (1.55 ms
%   -> 0.12 ms for the 240 x 80 blunt-cone operator on 56 threads). Below
%   5e4 nonzeros the plain product is kept: there the threads cost more
%   than they save (about 15%).
%
%   Inputs:
%       M       - Sparse matrix (real or complex).
%       use_gpu - true to keep M on the GPU (x must then be a gpuArray).
%
%   Outputs:
%       mul - Function handle, mul(x) = M * x. It holds only its copy of M.
%
% Part of: Hypersonics Stability MATLAB Solver - Stability Analysis Module

    if use_gpu
        Mg = gpuArray(M);
        mul = @(x) Mg * x;
    elseif nnz(M) < 5e4
        mul = @(x) M * x;   % thread start-up outweighs the work
    else
        Mt = M';
        mul = @(x) Mt' * x;
    end
end
