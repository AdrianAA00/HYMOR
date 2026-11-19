function fields = PERTURBATION_FIELDS(s)
% PERTURBATION_FIELDS  Conservative variables of the linear state vector, in order.
%
%   fields = PERTURBATION_FIELDS(s)
%
%   Returns ["rho", "rho_u", "rho_v", "rho_E"], followed by "rho_w" (the
%   azimuthal momentum) when an azimuthal mode is set up (s.azimuthal, see
%   SETUP_AZIMUTHAL_MODE). Each block of the state vector holds Nchi*Neta
%   cell values; shock-point DoFs, if any, follow the last block.
%
%   Inputs:
%       s - Solution structure.
%
%   Outputs:
%       fields - (1 x nv) string array of s.var / s.flux field names.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    fields = ["rho", "rho_u", "rho_v", "rho_E"];
    if isfield(s, 'azimuthal')
        fields(end+1) = "rho_w";
    end
end
