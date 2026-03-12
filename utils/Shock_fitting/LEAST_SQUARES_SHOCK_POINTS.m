function s = LEAST_SQUARES_SHOCK_POINTS(s)
% LEAST_SQUARES_SHOCK_POINTS - Fit and resample shock points using csaps smoothing spline.
%
%   s = LEAST_SQUARES_SHOCK_POINTS(s)
%
%   Fits a cubic smoothing spline (csaps) through the current shock
%   point coordinates, re-evaluates the spline at the grid parametric
%   locations (chi), and converts the new Cartesian coordinates back
%   to the element-space representation (eta, chi).
%
%   Inputs:
%       s (struct) - Solution structure containing at minimum:
%                    .shock.spline_param   - Smoothing parameter for csaps, p in [0, 1]
%                    .shock.points_chi     - (Nx x 1) parametric coordinates of shock points
%                    .shock.points_x/y     - (Nx x 1) Cartesian shock coordinates
%                    .mesh.chi             - (Nx x Ny) grid parametric coordinates
%
%   Outputs:
%       s (struct) - Updated s with:
%                    .shock.spline_func_x/y - Stored spline functions (pp-form)
%                    .shock.points_x/y      - Resampled Cartesian coordinates
%                    .shock.points_eta      - Radial coordinates in element space
%                    .shock.points_chi      - Updated parametric coordinates
%
%   See also: csaps, ppval, fnder, GO_TO_ELEMENT_SPACE
%
% Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

    %% Fit cubic smoothing spline to shock points
    s.shock.spline_func_x = csaps(s.shock.points_chi, s.shock.points_x, s.shock.spline_param);
    s.shock.spline_func_y = csaps(s.shock.points_chi, s.shock.points_y, s.shock.spline_param);

    %% Resample shock points at grid parametric locations
    s.shock.points_x = ppval(s.shock.spline_func_x, s.mesh.chi(:, 1));
    s.shock.points_y = ppval(s.shock.spline_func_y, s.mesh.chi(:, 1));

    %% Convert back to element space
    [s.shock.points_chi, s.shock.points_eta] = GO_TO_ELEMENT_SPACE(s.shock.points_x, s.shock.points_y, s);
end
