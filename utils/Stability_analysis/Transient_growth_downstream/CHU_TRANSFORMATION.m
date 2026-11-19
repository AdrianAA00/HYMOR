function R = CHU_TRANSFORMATION(rho, u, v, p, gamma_star, nv)
% CHU_TRANSFORMATION  Map conservative perturbations to Chu variables.
%
%   R = CHU_TRANSFORMATION(rho, u, v, p, gamma_star, nv)
%
%   Builds the sparse matrix R taking [rho; rho*u; rho*v; rho*E; (rho*w)]'
%   to [p; u; v; S; (w)]' at n points, linearized about a base flow with
%   zero azimuthal velocity (so w' = (rho*w)' / rho and the other rows are
%   unaffected by rho*w). Blocks are ordered by variable, n entries each.
%
%   Inputs:
%       rho, u, v, p, gamma_star - Base-flow fields (any shape, n points).
%       nv - Number of variables (4, or 5 for an azimuthal mode).
%
%   Outputs:
%       R - Sparse (nv*n x nv*n) transformation matrix.
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Downstream Module

    n = numel(rho);
    rho = rho(:);  u = u(:);  v = v(:);  p = p(:);  g = gamma_star(:);

    % {row variable, column variable, values}
    entries = {
        1, 1, (g - 1) .* (u.^2 + v.^2) ./ 2
        1, 2, -u .* (g - 1)
        1, 3, -v .* (g - 1)
        1, 4, g - 1
        2, 1, -u ./ rho
        2, 2, 1 ./ rho
        3, 1, -v ./ rho
        3, 3, 1 ./ rho
        4, 1, (u.^2 + v.^2) ./ (2*p) - g ./ (g - 1) ./ rho
        4, 2, -u ./ p
        4, 3, -v ./ p
        4, 4, 1 ./ p};
    if nv == 5
        entries(end+1, :) = {5, 5, 1 ./ rho};
    end

    idx = (1:n)';
    rows = cellfun(@(r) (r - 1) * n + idx, entries(:, 1), 'UniformOutput', false);
    cols = cellfun(@(c) (c - 1) * n + idx, entries(:, 2), 'UniformOutput', false);
    R = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(entries{:, 3}), nv * n, nv * n);
end
