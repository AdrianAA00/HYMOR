function p = UPDATE_PRESSURE(rho, T, chemistry)
% UPDATE_PRESSURE  Compute pressure from density and temperature using chemistry tables.
%
%   p = UPDATE_PRESSURE(rho, T, chemistry)
%
%   Evaluates the pressure using the equation of state p = rho * e * (gamma* - 1),
%   where internal energy and gamma_star are obtained from the chemistry
%   lookup tables. Falls back to zero for non-chemistry cases (placeholder).
%
%   Inputs:
%       rho       - Density field (array)
%       T         - Temperature field (array, same size as rho)
%       chemistry - Chemistry model struct with evaluation functions:
%                   .eval_e          - Internal energy from (T, rho)
%                   .eval_gamma_star - Ratio of specific heats from (rho, e)
%
%   Outputs:
%       p         - Pressure field (array, same size as rho)
%
%   Notes:
%       - When chemistry is inactive, the function returns p = 0 as a
%         placeholder.
%       - This function references s.chemistry.is_chemistry_enabled from
%         the caller's workspace (not passed as a parameter).
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% Compute pressure from equation of state
    if s.chemistry.is_chemistry_enabled
        e = chemistry.eval_e(T, rho);
        gamma_star = chemistry.eval_gamma_star(rho, e);
        p = rho .* e .* (gamma_star - 1);
    else
        p = 0; % Not valid for case without chemistry
    end
end
