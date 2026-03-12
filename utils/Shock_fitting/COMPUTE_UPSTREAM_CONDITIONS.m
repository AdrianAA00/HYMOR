function [p_infty, u_inf, v_inf, rho_inf, e_inf] = COMPUTE_UPSTREAM_CONDITIONS(s)
% COMPUTE_UPSTREAM_CONDITIONS - Evaluate freestream primitive variables at shock points.
%
%   [p_infty, u_inf, v_inf, rho_inf, e_inf] = COMPUTE_UPSTREAM_CONDITIONS(s)
%
%   Retrieves the upstream (pre-shock) conservative variables from
%   UPDATE_SHOCK_UPSTREAM, converts them to primitive form (density,
%   velocity components, internal energy), and computes the upstream
%   pressure using the equation of state.
%
%   Inputs:
%       s (struct) - Solution structure containing upstream flow parameters,
%                    freestream conditions, perturbation amplitudes, and
%                    shock-point coordinates (see UPDATE_SHOCK_UPSTREAM).
%
%   Outputs:
%       p_infty (Nx x 1 double) - Upstream pressure at each shock point.
%       u_inf   (Nx x 1 double) - Upstream x-velocity at each shock point.
%       v_inf   (Nx x 1 double) - Upstream y-velocity at each shock point.
%       rho_inf (Nx x 1 double) - Upstream density at each shock point.
%       e_inf   (Nx x 1 double) - Upstream specific internal energy at each
%                                  shock point.
%
%   Notes:
%       - Pressure is computed from the calorically perfect gas relation:
%         p = rho * e * (gamma_star - 1).
%
%   See also: UPDATE_SHOCK_UPSTREAM (local subfunction)
%
% Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

    [rho, rho_u, rho_v, rho_E] = UPDATE_SHOCK_UPSTREAM(s);

    rho_inf = rho;
    v_inf = rho_v ./ rho_inf;
    u_inf = rho_u ./ rho_inf;
    u_mag = sqrt(v_inf.^2 + u_inf.^2);
    e_inf = (rho_E - rho_inf .* u_mag.^2 / 2) ./ rho_inf;
    p_infty = (e_inf .* rho_inf) * (s.freestream.gamma_star - 1);
end

function [rho, rho_u, rho_v, rho_E] = UPDATE_SHOCK_UPSTREAM(s)
% UPDATE_SHOCK_UPSTREAM - Build upstream conservative state with optional perturbations.
%
%   Computes upstream density, momentum, and total energy at the shock
%   points, either from a prescribed perturbed field (rho_0_upstream_p)
%   or by superimposing a travelling-wave perturbation on the uniform
%   freestream state. Enforces axisymmetric boundary conditions when
%   the dimension flag is set to "3D-axisymmetric".

    u_y = s.freestream.rho_v_0 / s.freestream.rho_0;
    u_x = s.freestream.rho_u_0 / s.freestream.rho_0;
    perturbation = cos(2 * pi * s.freestream.disturbance.k_y * (s.shock.points_y - u_y * s.time_integration.t)) ...
               .* sin(2 * pi * s.freestream.disturbance.k_x * (s.shock.points_x - u_x * s.time_integration.t));

    if isfield(s, 'rho_0_upstream_p')
        rho   = s.freestream.rho_0_p;
        rho_u = s.freestream.rho_u_0_p;
        rho_v = s.freestream.rho_v_0_p;
        rho_E = s.freestream.rho_E_0_p;
    else
        rho   = s.freestream.rho_0   + perturbation * s.freestream.disturbance.amplitude(1);
        rho_u = s.freestream.rho_u_0 + perturbation * s.freestream.disturbance.amplitude(2);
        rho_v = s.freestream.rho_v_0 + perturbation * s.freestream.disturbance.amplitude(3);
        rho_E = s.freestream.rho_E_0 + perturbation * s.freestream.disturbance.amplitude(4);
    end

    %% Apply symmetry boundary conditions for axisymmetric flows
    if s.PDE_dimension == "3D-axisymmetric"
        rho_v(end, 1) = 0;
    end
end
