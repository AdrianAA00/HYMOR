function [V, theta, flag] = BLOCK_LANCZOS(apply, D, k, b, tol, max_basis, is_real)
% BLOCK_LANCZOS  Largest eigenpairs of C v = theta D v by block Lanczos.
%
%   [V, theta, flag] = BLOCK_LANCZOS(apply, D, k, b, tol, max_basis, is_real)
%
%   For a Hermitian positive semi-definite C (given through its action on
%   blocks of b vectors) and a Hermitian positive definite D. With the
%   Cholesky factor D(p, p) = Rc' Rc the problem becomes the standard one
%   H y = theta y, H = Rc^-H P' C P Rc^-1, solved by block Lanczos with full
%   reorthogonalization and Rayleigh-Ritz on the projected matrix.
%   Applying C to b vectors at once is the point: the time stepper then
%   multiplies the sparse operator by an n x b block, about five times
%   cheaper per vector on a GPU than b separate products.
%
%   Inputs:
%       apply     - Handle, apply(X) = C * X for an n x b block X (CPU in and out).
%       D         - Sparse Hermitian positive definite matrix (n x n).
%       k         - Number of eigenpairs wanted.
%       b         - Block size.
%       tol       - Relative residual tolerance, ||H y - theta y|| <= tol * theta.
%       max_basis - Maximum Krylov basis size.
%       is_real   - true for a real problem (real start block).
%
%   Outputs:
%       V     - Eigenvectors (n x k), normalized as V' * D * V = I.
%       theta - Eigenvalues (k x 1), largest first.
%       flag  - 0 if all k pairs converged, 1 otherwise.
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Downstream Module

    n = size(D, 1);
    [Rc, chol_flag, p] = chol(D, 'vector');
    assert(chol_flag == 0, 'BLOCK_LANCZOS: D is not positive definite');
    H = @(Y) H_PRODUCT(Y, apply, Rc, p);

    %% Start block
    X = randn(n, b);
    if ~is_real
        X = complex(X, randn(n, b));
    end
    [Q, ~] = qr(X, 0);
    basis = zeros(n, max_basis + b, 'like', Q);
    basis(:, 1:b) = Q;
    T = zeros(max_basis + b, 'like', Q);
    m = b;

    %% Block Lanczos with full reorthogonalization
    while true
        last = m-b+1:m;
        W = H(basis(:, last));
        Wq = zeros(m, b, 'like', Q);
        for pass = 1:2
            h = basis(:, 1:m)' * W;
            W = W - basis(:, 1:m) * h;
            Wq = Wq + h;
        end
        T(1:m, last) = Wq;
        T(last, 1:m) = Wq';
        Tm = (T(1:m, 1:m) + T(1:m, 1:m)') / 2;
        [Y, th] = eig(Tm, 'vector');
        [th, o] = sort(real(th), 'descend');
        Y = Y(:, o);

        % Residual of each Ritz pair: the part of H basis Y leaving the basis
        [Qn, Bn] = qr(W, 0);
        kk = min(k, m);
        res = vecnorm(Bn * Y(last, 1:kk), 2, 1);
        converged = m >= k && all(res <= tol * abs(th(1:k)).');
        if converged || m + b > max_basis
            break
        end
        basis(:, m+1:m+b) = Qn;
        m = m + b;
    end
    flag = double(~converged);

    %% Ritz vectors back to the original variables
    V = zeros(n, k, 'like', Q);
    V(p, :) = Rc \ (basis(:, 1:m) * Y(:, 1:k));
    theta = th(1:k);
end


function Z = H_PRODUCT(Y, apply, Rc, p)
% H_PRODUCT  H * Y = Rc^-H P' C P Rc^-1 * Y.
    X = zeros(size(Y), 'like', Y);
    X(p, :) = Rc \ Y;
    Z = apply(X);
    Z = Rc' \ Z(p, :);
end
