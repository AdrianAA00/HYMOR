function [s] = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s, chemistry)
% NON_LINEAR_DYNAMICS_NO_DISCONTINUITY  Compute nonlinear dynamics with optional shock masking.
%
%   [s] = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s, chemistry)
%
%   Computes the PDE fluxes and optionally applies a shock-feedback mask
%   to zero out fluxes in cells containing discontinuities. When shock
%   feedback is active, the flux fields are multiplied element-wise by the
%   flow_cells mask so that only smooth-flow regions contribute.
%
%   Inputs:
%       s  - Solution struct containing flow state and shock settings:
%                   .shock.enabled  - Logical flag indicating shock presence
%                   .shock.feedback - Logical flag to enable flux masking
%                   .shock.flow_cells - Mask array (0 at shock, 1 elsewhere)
%                   .flux.rho       - Mass flux field
%                   .flux.rho_u     - x-momentum flux field
%                   .flux.rho_v     - y-momentum flux field
%                   .flux.rho_E     - Total energy flux field
%       chemistry - Chemistry model struct passed to PDE
%
%   Outputs:
%       s  - Updated s struct with computed (and optionally
%                   masked) flux fields
%
%   Notes:
%       - Calls PDE() internally to compute all flux contributions.
%       - When shock_feedback is disabled or no shock is present, fluxes
%         are left unchanged after the PDE call.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% Compute PDE fluxes
    s = PDE(s, chemistry);

    %% Apply shock-feedback mask if active
    if s.shock.enabled == true && s.shock.feedback
        s.flux.rho   = s.flux.rho   .* s.shock.flow_cells;
        s.flux.rho_u = s.flux.rho_u .* s.shock.flow_cells;
        s.flux.rho_v = s.flux.rho_v .* s.shock.flow_cells;
        s.flux.rho_E = s.flux.rho_E .* s.shock.flow_cells;
    end

end
