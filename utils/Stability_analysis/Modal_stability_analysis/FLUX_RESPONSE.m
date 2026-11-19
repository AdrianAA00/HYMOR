function response = FLUX_RESPONSE(s_pert, s_base, fields, epsilon, chemistry)
% FLUX_RESPONSE  Finite-difference Jacobian column of the discrete operator.
%
%   response = FLUX_RESPONSE(s_pert, s_base, fields, epsilon, chemistry)
%
%   Evaluates the perturbed state and returns (f(x+Dx) - f(x)) / Dx for
%   every flux field and for the shock speed.
%
%   Inputs:
%       s_pert    - Perturbed solution structure (before the evaluation).
%       s_base    - Base solution after NON_LINEAR_DYNAMICS_NO_DISCONTINUITY.
%       fields    - Flux fields to difference (PERTURBATION_FIELDS).
%       epsilon   - Perturbation magnitude.
%       chemistry - Chemistry model structure.
%
%   Outputs:
%       response - Struct with
%                    .flux        - cell array of (Nchi x Neta) differences,
%                                   ordered as fields
%                    .shock_speed - (Nchi x 1) difference of the shock speed
%                                   (shock-fitted cases only)
%
% Part of: Hypersonics Stability MATLAB Solver - Stability Analysis / Linear Stability Module

    s_pert = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s_pert, chemistry);
    response.flux = cell(numel(fields), 1);
    for o = 1:numel(fields)
        response.flux{o} = (s_pert.flux.(fields(o)) - s_base.flux.(fields(o))) / epsilon;
    end
    if s_base.shock.enabled
        response.shock_speed = (s_pert.shock.relative_increase_velocity ...
            - s_base.shock.relative_increase_velocity) / epsilon;
    end
end
