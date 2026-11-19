function ref = PIPE_FLOW_REFERENCE(Re, alpha, m, N)
% PIPE_FLOW_REFERENCE  Spectral reference for incompressible Hagen-Poiseuille stability.
%
%   ref = PIPE_FLOW_REFERENCE(Re, alpha, m, N)
%
%   Linear stability of U(r) = 1 - r^2 in a pipe of unit radius (velocity
%   scale: centreline velocity) to perturbations
%       [u_r, u_theta, u_x, p](r) * exp(1i*(alpha*x + m*theta) + lambda*t),
%   solved in primitive variables with Chebyshev collocation on r in (0, 1]
%   (even/odd parity across the axis, Trefethen's "Spectral Methods in
%   MATLAB", programs 29/36). It is independent of HYMOR and serves to
%   validate the azimuthal-mode operators of the compressible solver in
%   the low-Mach limit.
%
%   Inputs:
%       Re    - Reynolds number U_c R / nu.
%       alpha - Axial wavenumber.
%       m     - Azimuthal wavenumber.
%       N     - Number of radial collocation points in (0, 1].
%
%   Outputs:
%       ref - Struct with fields
%             .lambda - finite eigenvalues, sorted by decreasing real part
%             .V      - velocity eigenvectors [u_r; u_theta; u_x] at the
%                       interior points (columns ordered as lambda)
%             .w      - quadrature weights of the energy norm
%                       (1/2) int |u|^2 r dr, repeated for the 3 components
%             .r      - interior collocation points
%             .gain   - function handle G = ref.gain(T): optimal energy gain
%                       max |q(T)|^2 / |q(0)|^2 (kinetic energy norm) over
%                       the span of the eigenvectors with real(lambda) > -1
%             .eval   - function handle q = ref.eval(k, r_query) returning the
%                       k-th eigenvector [u_r, u_theta, u_x, p] (columns)
%                       interpolated at the radii r_query (column vector)
%
% Part of: Hypersonics Stability MATLAB Solver - Tutorials (pipe flow)

    %% Chebyshev grid on [-1, 1] with an odd number of intervals (no point at r = 0)
    Nf = 2 * N + 1;
    x = cos(pi * (0:Nf)' / Nf);
    D = CHEB_DIFFERENTIATION(x);
    D2 = D * D;
    r = x(1:N+1);                       % r = 1 (wall) ... smallest positive point
    mirror = Nf + 1 : -1 : Nf + 1 - N;  % columns of the points at -r

    % Fold onto r > 0 for a field of parity sigma: f(-r) = sigma f(r)
    fold = @(A, sigma) A(1:N+1, 1:N+1) + sigma * A(1:N+1, mirror);
    s_scalar = (-1)^m;                  % p, u_x
    s_vector = -(-1)^m;                 % u_r, u_theta

    I = 2:N+1;                          % interior points (velocities vanish at the wall)
    n = numel(I);
    R = diag(r(I));
    Ri = diag(1 ./ r(I));
    U = diag(1 - r(I).^2);
    dU = diag(-2 * r(I));
    Id = eye(n);
    Z = zeros(n);

    %% Operators
    lap = @(sigma) fold(D2, sigma) + diag(1 ./ r) * fold(D, sigma) - m^2 * diag(1 ./ r.^2) - alpha^2 * eye(N+1);
    Lv = lap(s_vector);
    Ls = lap(s_scalar);
    Dv = fold(D, s_vector);
    Dp = fold(D, s_scalar);
    adv = -1i * alpha * U;

    A_rr = adv + (Lv(I, I) - Ri^2) / Re;
    A_rt = -2i * m * Ri^2 / Re;
    A_tr = 2i * m * Ri^2 / Re;
    A_tt = adv + (Lv(I, I) - Ri^2) / Re;
    A_xx = adv + Ls(I, I) / Re;
    A_xr = -dU;

    % Pressure on all N+1 points (wall included); continuity on all N+1 points
    Grad_r = -Dp(I, :);
    Grad_t = -1i * m * Ri * [zeros(n, 1), eye(n)];
    Grad_x = -1i * alpha * [zeros(n, 1), eye(n)];
    Div_r = Dv(:, I) + diag(1 ./ r) * [zeros(1, n); eye(n)];
    Div_t = 1i * m * diag(1 ./ r) * [zeros(1, n); eye(n)];
    Div_x = 1i * alpha * [zeros(1, n); eye(n)];

    A = [A_rr, A_rt, Z,    Grad_r
         A_tr, A_tt, Z,    Grad_t
         A_xr, Z,    A_xx, Grad_x
         Div_r, Div_t, Div_x, zeros(N+1)];
    B = blkdiag(Id, Id, Id, zeros(N+1));

    %% Finite eigenvalues
    [V, L] = eig(A, B, 'vector');
    keep = isfinite(L) & abs(L) < 1e3;
    [~, order] = sort(real(L(keep)), 'descend');
    L = L(keep);
    V = V(:, keep);
    V = V(:, order);
    ref.lambda = L(order);
    ref.V = V(1:3*n, :);
    sigma = [s_vector, s_vector, s_scalar, s_scalar];
    ref.eval = @(k, rq) EVALUATE_MODE(V(:, k), n, N, x, sigma, rq);

    %% Clenshaw-Curtis weights of int_0^1 f r dr (symmetric extension to [-1, 1])
    wcc = CLENSHAW_CURTIS(Nf);
    w = (wcc(I) + wcc(mirror(I))) / 2 .* r(I) / 2;   % (1/2) |u|^2 r
    ref.w = repmat(w(:), 3, 1);
    ref.r = r(I);
    ref.gain = @(T) OPTIMAL_GAIN(ref, T);
end

%% ========================================================================
%  Optimal Transient Growth (eigenvector expansion)
%  ========================================================================
function G = OPTIMAL_GAIN(ref, T)
% OPTIMAL_GAIN  G(T) = || F exp(Lambda T) F^-1 ||_2^2 with F' F = V' W V
%   (Reddy & Henningson), restricted to the eigenvectors with real(lambda) > -1.
    keep = real(ref.lambda) > -1;
    V = ref.V(:, keep);
    F = chol(V' * (ref.w .* V));
    G = zeros(size(T));
    for k = 1:numel(T)
        G(k) = norm(F * diag(exp(ref.lambda(keep) * T(k))) / F)^2;
    end
end

%% ========================================================================
%  Local Helpers
%  ========================================================================
function q = EVALUATE_MODE(v, n, N, x, sigma, rq)
% EVALUATE_MODE  Interpolate one eigenvector at the radii rq (barycentric
%   Chebyshev interpolation of the parity extension on [-1, 1]).
    Nf = numel(x) - 1;
    wb = (-1).^(0:Nf)';
    wb([1, end]) = wb([1, end]) / 2;
    vel = {[0; v(1:n)], [0; v(n+1:2*n)], [0; v(2*n+1:3*n)], v(3*n+1:end)};
    q = zeros(numel(rq), 4);
    for c = 1:4
        f_half = vel{c};                                    % values at r = x(1:N+1)
        f = [f_half; sigma(c) * f_half(end:-1:1)];          % values at x(N+2:end) = -x(N+1:-1:1)
        for i = 1:numel(rq)
            d = rq(i) - x;
            if any(d == 0)
                q(i, c) = f(d == 0);
            else
                q(i, c) = sum(wb .* f ./ d) / sum(wb ./ d);
            end
        end
    end
end

function D = CHEB_DIFFERENTIATION(x)
% CHEB_DIFFERENTIATION  Chebyshev differentiation matrix on the points x.
    Nf = numel(x) - 1;
    c = [2; ones(Nf - 1, 1); 2] .* (-1).^(0:Nf)';
    X = repmat(x, 1, Nf + 1);
    dX = X - X';
    D = (c * (1 ./ c)') ./ (dX + eye(Nf + 1));
    D = D - diag(sum(D, 2));
end

function w = CLENSHAW_CURTIS(Nf)
% CLENSHAW_CURTIS  Clenshaw-Curtis weights on cos(pi*(0:Nf)/Nf).
    theta = pi * (0:Nf)' / Nf;
    w = zeros(Nf + 1, 1);
    ii = 2:Nf;
    v = ones(Nf - 1, 1);
    if mod(Nf, 2) == 0
        w([1, end]) = 1 / (Nf^2 - 1);
        for k = 1:Nf/2 - 1
            v = v - 2 * cos(2 * k * theta(ii)) / (4 * k^2 - 1);
        end
        v = v - cos(Nf * theta(ii)) / (Nf^2 - 1);
    else
        w([1, end]) = 1 / Nf^2;
        for k = 1:(Nf - 1)/2
            v = v - 2 * cos(2 * k * theta(ii)) / (4 * k^2 - 1);
        end
    end
    w(ii) = 2 * v / Nf;
end
