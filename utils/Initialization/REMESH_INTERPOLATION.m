function s = REMESH_INTERPOLATION(s, solution_old)
% REMESH_INTERPOLATION  Transfer the flow fields of solution_old onto the
%   mesh of s by cubic interpolation in the (i, j) index space of the old
%   structured grid.
%
%   s = REMESH_INTERPOLATION(s, solution_old)
%
%   Inputs:
%       s            - Solution struct whose extended mesh (x_Ext, y_Ext)
%                      is the target
%       solution_old - Solution struct holding the fields on the old mesh
%
%   Outputs:
%       s - Same struct with s.var.(rho, rho_u, rho_v, rho_E[, rho_w,
%           gamma_star, cv_star]) filled on the extended (ghost) grid.
%
%   Notes:
%       - Each new cell centre is located in the old grid by inverting the
%         map (i, j) -> (x, y) of the old extended cell centres with Newton
%         iterations (started from the nearest old node). The fields are
%         then interpolated with tensor-product cubic convolution (Keys,
%         a = -1/2) on the uniform (i, j) lattice: local 4 x 4 stencil, C1,
%         third-order accurate.
%       - Interpolating along the grid lines keeps the interpolation error
%         smooth on stretched, curved, high-aspect-ratio grids. Delaunay-
%         linear interpolation in (x, y) instead picks triangle diagonals
%         that vary from cell to cell, and its error is amplified by the
%         discrete derivatives (larger residual jump after a remesh, mostly
%         in the boundary layer and at the stagnation region).
%       - With shock fitting, the old fields from the shocked cell up are
%         replaced by a linear continuation of the last two flow cells, so
%         no stencil straddles the shock or uses the shocked-cell closure.
%         Corner ghost cells (not set by the boundary conditions) are filled
%         by bilinear extrapolation.
%       - New points outside the old grid (ghost and freestream cells) are
%         extrapolated; they are overwritten afterwards by the boundary
%         conditions and the shock treatment.
%
% Part of: Hypersonics Stability MATLAB Solver - Initialization Module

    %% Fields to transfer
    fields = PERTURBATION_FIELDS(solution_old);
    if ~s.chemistry.chemical_equilibrium
        fields = [fields, "gamma_star", "cv_star"];
    end

    %% Locate the new cell centres in the old index space
    [xi, eta] = LOCATE_IN_OLD_GRID(solution_old.mesh.x_Ext, solution_old.mesh.y_Ext, ...
        s.mesh.x_Ext, s.mesh.y_Ext);
    W = KEYS_WEIGHTS(xi, eta, size(solution_old.mesh.x_Ext), false);

    %% Interpolate every field
    for f = fields
        F = PREPARE_OLD_FIELD(solution_old.var.(f), solution_old);
        s.var.(f) = KEYS_EVAL(F, W);
    end
end


function [xi, eta] = LOCATE_IN_OLD_GRID(xo, yo, xn, yn)
% Invert the cubic-convolution map (i, j) -> (x, y) of the old grid at each
% new point with damped Newton iterations, started from the nearest old
% node: nearest old grid line (matched on the first interior row), then
% nearest node along that line.
    [Ni, Nj] = size(xo);
    [Mi, Mj] = size(xn);

    [~, k] = min((xn(:, 2) - xo(:, 2).').^2 + (yn(:, 2) - yo(:, 2).').^2, [], 2);
    xi  = repmat(k, 1, Mj);
    eta = zeros(Mi, Mj);
    for i = 1:Mi
        [~, eta(i, :)] = min((xn(i, :) - xo(k(i), :).').^2 + (yn(i, :) - yo(k(i), :).').^2, [], 1);
    end

    for it = 1:50
        W = KEYS_WEIGHTS(xi, eta, [Ni, Nj], true);
        [X, X_i, X_j] = KEYS_EVAL(xo, W);
        [Y, Y_i, Y_j] = KEYS_EVAL(yo, W);
        rx = xn - X;
        ry = yn - Y;
        det = X_i .* Y_j - X_j .* Y_i;
        dxi  = ( Y_j .* rx - X_j .* ry) ./ det;
        deta = (-Y_i .* rx + X_i .* ry) ./ det;
        step = max(1, max(abs(dxi), abs(deta)));   % at most one cell per iteration
        xi  = xi  + dxi  ./ step;
        eta = eta + deta ./ step;
        if max(abs(dxi(:)) + abs(deta(:))) < 1e-10
            break
        end
    end
    located = isfinite(xi) & isfinite(eta) & abs(dxi) + abs(deta) <= 1e-6;
    if ~all(located(2:end-1, 2:end-1), 'all')   % ghost cells are reset by the BCs
        error('REMESH_INTERPOLATION: new cell centres could not be located in the old grid.');
    end
end


function F = PREPARE_OLD_FIELD(F, solution_old)
% Replace the shocked cell and everything beyond it by a linear
% continuation of the last two flow cells, and fill the four corner ghost
% cells. The shocked cell holds the Rankine-Hugoniot reconstruction, not a
% smooth continuation of the flow field: interpolating through it gives a
% larger residual in the first row below the shock (recomputed anyway by
% UPDATE_SHOCK_BC on the new mesh).
    [Ni, Nj] = size(F);
    if solution_old.shock.enabled
        js = solution_old.shock.cell_indices(:, 1);
        c_col = [js(1); js; js(end)] + 1;   % shocked cell (extended indexing)
        for i = 1:Ni
            m1 = c_col(i) - 1;              % last flow cell below the shock
            k = 1:(Nj - m1);
            F(i, m1 + k) = F(i, m1) + k * (F(i, m1) - F(i, m1 - 1));
        end
    end
    F(1, 1)   = F(2, 1)      + F(1, 2)      - F(2, 2);
    F(Ni, 1)  = F(Ni-1, 1)   + F(Ni, 2)     - F(Ni-1, 2);
    F(1, Nj)  = F(2, Nj)     + F(1, Nj-1)   - F(2, Nj-1);
    F(Ni, Nj) = F(Ni-1, Nj)  + F(Ni, Nj-1)  - F(Ni-1, Nj-1);
end


function W = KEYS_WEIGHTS(xi, eta, sz, derivatives)
% Stencil indices and cubic-convolution weights (Keys 1981, a = -1/2) at
% fractional indices (xi, eta) of a grid of size sz. The grid is padded by
% one node per side with Keys' boundary rule; queries outside the grid are
% extrapolated with the boundary cubic.
    [W.ki, W.ti] = BASE(xi, sz(1));
    [W.kj, W.tj] = BASE(eta, sz(2));
    W.wi = KERNEL(W.ti);  W.wj = KERNEL(W.tj);
    if derivatives
        W.di = KERNEL_DT(W.ti);  W.dj = KERNEL_DT(W.tj);
    end
    W.sz = sz;
    W.derivatives = derivatives;
end


function [k, t] = BASE(q, N)
    k = min(max(floor(q), 1), N - 1);   % left node, 1 <= k <= N-1
    t = q - k;
end


function w = KERNEL(t)
    w = {(-t.^3 + 2*t.^2 - t) / 2, (3*t.^3 - 5*t.^2 + 2) / 2, ...
         (-3*t.^3 + 4*t.^2 + t) / 2, (t.^3 - t.^2) / 2};
end


function w = KERNEL_DT(t)
    w = {(-3*t.^2 + 4*t - 1) / 2, (9*t.^2 - 10*t) / 2, ...
         (-9*t.^2 + 8*t + 1) / 2, (3*t.^2 - 2*t) / 2};
end


function [v, v_i, v_j] = KEYS_EVAL(F, W)
% Evaluate the interpolant of F (and its index derivatives) at the stencils W.
    [Ni, Nj] = size(F);
    P = zeros(Ni + 2, Nj + 2);                 % padded copy, node (i,j) -> (i+1,j+1)
    P(2:end-1, 2:end-1) = F;
    P(1, 2:end-1)   = 3*F(1, :)   - 3*F(2, :)     + F(3, :);
    P(end, 2:end-1) = 3*F(end, :) - 3*F(end-1, :) + F(end-2, :);
    P(:, 1)   = 3*P(:, 2)   - 3*P(:, 3)     + P(:, 4);
    P(:, end) = 3*P(:, end-1) - 3*P(:, end-2) + P(:, end-3);

    v = zeros(size(W.ki));  v_i = v;  v_j = v;
    for a = 1:4
        ii = W.ki + a - 1;                     % padded row of node k-2+a
        for b = 1:4
            jj = W.kj + b - 1;
            Fab = P(ii + (jj - 1) * (Ni + 2));
            v = v + W.wi{a} .* W.wj{b} .* Fab;
            if W.derivatives
                v_i = v_i + W.di{a} .* W.wj{b} .* Fab;
                v_j = v_j + W.wi{a} .* W.dj{b} .* Fab;
            end
        end
    end
end
