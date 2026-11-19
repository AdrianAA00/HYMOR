function W = CHU_NORM_WEIGHTS(rho, a, gamma_star, p, norms, nv)
% CHU_NORM_WEIGHTS  Chu energy weights of the variables [p, u, v, S, (w)].
%
%   W = CHU_NORM_WEIGHTS(rho, a, gamma_star, p, norms, nv)
%
%   Returns, per point and per unit measure (volume or area), the diagonal
%   weights of Chu's energy norm in the variables of CHU_TRANSFORMATION:
%       pressure: rho a^2 / (2 (gamma p)^2)
%       velocity: rho / 2                  (u, v and, for nv = 5, w)
%       entropy:  (gamma - 1) p / (2 gamma)
%
%   Inputs:
%       rho, a, gamma_star, p - Base-flow fields (any shape, n points).
%       norms - Logical vector [pressure, u, v, entropy, (w)] selecting the
%               components; the azimuthal velocity follows norms(3) (the
%               radial velocity) when norms has only four entries.
%       nv    - Number of variables (4, or 5 for an azimuthal mode).
%
%   Outputs:
%       W - (n x nv) weights; the columns of deselected components are zero.
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Downstream Module

    if numel(norms) < 5
        norms(5) = norms(3);
    end
    W = zeros(numel(rho), nv);
    if norms(1)
        W(:, 1) = rho(:) .* a(:).^2 ./ (2 .* (gamma_star(:) .* p(:)).^2);
    end
    for k = [2, 3, 5]
        if k <= nv && norms(k)
            W(:, k) = rho(:) / 2;
        end
    end
    if norms(4)
        W(:, 4) = (gamma_star(:) - 1) .* p(:) ./ (2 * gamma_star(:));
    end
end
