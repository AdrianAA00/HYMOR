function s = CLEAR_AZIMUTHAL_MODE(s)
% CLEAR_AZIMUTHAL_MODE  Remove an azimuthal-mode set-up from the solution.
%
%   s = CLEAR_AZIMUTHAL_MODE(s)
%
%   Removes s.azimuthal, s.var.rho_w and s.flux.rho_w (SETUP_AZIMUTHAL_MODE,
%   FLUX_AZIMUTHAL_MODE), so that the solution describes the axisymmetric
%   (or 2D) 4-variable problem again.
%
%   Inputs:
%       s - Solution structure.
%
%   Outputs:
%       s - Solution structure without azimuthal-mode fields.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    if isfield(s, 'azimuthal')
        s = rmfield(s, 'azimuthal');
    end
    if isfield(s.var, 'rho_w')
        s.var = rmfield(s.var, 'rho_w');
    end
    if isfield(s.flux, 'rho_w')
        s.flux = rmfield(s.flux, 'rho_w');
    end
end
