function s = COMPUTE_BETA(s)
% COMPUTE_BETA - Compute the local shock wave angle at each sampled point.
%
%   s = COMPUTE_BETA(s)
%
%   Evaluates the first derivative of the shock spline in both x and y
%   directions, then combines with the upstream flow direction to obtain
%   the local oblique-shock angle beta at every shock sample point.
%
%   Inputs:
%       s (struct) - Solution structure containing at minimum:
%                    .shock.spline_func_x - Spline (pp-form) for shock x-coordinates
%                    .shock.spline_func_y - Spline (pp-form) for shock y-coordinates
%                    .mesh.chi            - (Nx x 1) parametric coordinate of shock points
%                    (plus fields required by COMPUTE_UPSTREAM_CONDITIONS)
%
%   Outputs:
%       s (struct) - Updated s with:
%                    .shock.beta - (Nx x 1) local shock wave angle [rad]
%
%   Notes:
%       - The shock angle beta is defined as the angle between the upstream
%         velocity vector and the shock-tangent vector.
%       - Uses MATLAB Curve Fitting Toolbox functions fnder and ppval.
%
%   See also: COMPUTE_UPSTREAM_CONDITIONS, fnder, ppval
%
% Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

    %% Compute shock tangent derivatives
    fderiv_x = fnder(s.shock.spline_func_x, 1);
    deriv_x = ppval(fderiv_x, s.mesh.chi(:, 1));
    fderiv_y = fnder(s.shock.spline_func_y, 1);
    deriv_y = ppval(fderiv_y, s.mesh.chi(:, 1));

    %% Compute shock angle (independent of chi traversal direction)
    [~, u_inf, v_inf, ~, ~] = COMPUTE_UPSTREAM_CONDITIONS(s);

    % Cross and dot products between shock tangent and upstream velocity
    cross_tv = deriv_x .* v_inf - deriv_y .* u_inf;
    dot_tv   = deriv_x .* u_inf + deriv_y .* v_inf;

    % Acute angle between shock tangent line and upstream velocity
    s.shock.beta = atan2(abs(cross_tv), abs(dot_tv));
end
