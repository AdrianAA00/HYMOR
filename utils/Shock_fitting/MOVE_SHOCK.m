function s = MOVE_SHOCK(s, dsolution_dt, w)
% MOVE_SHOCK - Advance shock point positions by the computed shock speed.
%
%   s = MOVE_SHOCK(s, dsolution_dt, w)
%
%   Updates the shock point coordinates by adding the product of the
%   shock speed, the Runge-Kutta weight, and the user-defined shock
%   relaxation factor. This is called during the time-integration
%   loop to evolve the shock location.
%
%   Inputs:
%       s     (struct) - Solution structure containing:
%                               .shock.points_x/y    - (Nx x 1) current shock coordinates
%                               .shock.relaxation    - (scalar) under-relaxation factor
%       dsolution_dt (struct) - Time-derivative structure containing:
%                               .shock.speed_x/y     - (Nx x 1) shock velocity components
%       w            (double) - Runge-Kutta stage weight.
%
%   Outputs:
%       s     (struct) - Updated s with new shock coordinates.
%
%   Notes:
%       - The relaxation factor allows under-relaxation of the shock
%         motion for stability (shock_relaxation < 1 damps oscillations).
%
%   See also: UPDATE_SHOCK_JUMP_PROPERTIES
%
% Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

    s.shock.points_x = s.shock.points_x ...
        + dsolution_dt.shock.speed_x * w * s.shock.relaxation;
    s.shock.points_y = s.shock.points_y ...
        + dsolution_dt.shock.speed_y * w * s.shock.relaxation;
end
