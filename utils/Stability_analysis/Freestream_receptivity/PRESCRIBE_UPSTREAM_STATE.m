function [s, upstream_fields] = PRESCRIBE_UPSTREAM_STATE(s)
% PRESCRIBE_UPSTREAM_STATE  Hold the upstream state at the shock points.
%
%   [s, upstream_fields] = PRESCRIBE_UPSTREAM_STATE(s)
%
%   Replaces the analytic freestream disturbance by a prescribed upstream
%   state per shock point, initialised to the undisturbed freestream, so
%   that the upstream variables can be perturbed point by point
%   (UPDATE_SHOCK_UPSTREAM, UPDATE_FIELD_UPSTREAM).
%
%   Inputs:
%       s - Solution structure with the freestream state.
%
%   Outputs:
%       s               - Same structure with s.freestream.<field> (one value
%                         per shock point) and s.rho_0_upstream_p set.
%       upstream_fields - Names of the prescribed fields, ordered as
%                         PERTURBATION_FIELDS: "rho_0_p", "rho_u_0_p",
%                         "rho_v_0_p", "rho_E_0_p" and, for an azimuthal
%                         mode, "rho_w_0_p" (zero in the freestream).
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Freestream Module

    upstream_fields = PERTURBATION_FIELDS(s) + "_0_p";
    for f = upstream_fields
        base_field = extractBefore(f, "_p");
        if isfield(s.freestream, base_field)
            s.freestream.(f) = s.freestream.(base_field) * ones(s.mesh.Nchi, 1);
        else
            s.freestream.(f) = zeros(s.mesh.Nchi, 1);   % rho_w: no swirl upstream
        end
    end
    s.rho_0_upstream_p = true;
end
