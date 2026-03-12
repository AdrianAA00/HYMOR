function s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry)
% UPDATE_THERMODYNAMIC_PROPERTIES  Compute pressure and temperature from conserved variables.
%
%   s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry)
%
%   Evaluates pressure using p = rho * e * (gamma* - 1) and temperature
%   using T = e * energy_factor / cv_star. 
%
%   Inputs:
%       s  - Solution struct containing:
%                   .var.rho, .var.rho_u, .var.rho_v, .var.rho_E - Conserved variables
%                   .var.gamma_star      - Ratio of specific heats field
%                   .var.cv_star         - Specific heat at constant volume
%       chemistry - Chemistry model struct (reserved for future use)
%
%   Outputs:
%       s  - Updated s struct with:
%                   .var.p - Pressure field (non-dimensional, rho * U^2)
%                   .var.T - Temperature field (non-dimensional, U^2 / cv_infty)
%
%   Notes:
%       - Internal energy: e = (rho_E - 0.5*(rho_u^2 + rho_v^2)/rho) / rho
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% Compute internal energy, pressure, and temperature
    rho = s.var.rho;
    e = (s.var.rho_E - (s.var.rho_u.^2 + s.var.rho_v.^2) ./ rho / 2) ./ rho;
    s.var.p = rho .* e .* (s.var.gamma_star - 1);               % Non-dimensional rho * U^2
    s.var.T = e ./ s.var.cv_star;                               % Non-dimensional U^2 / cv_infty

    %% Wall boundary condition at eta0, keep pressure smooth, just update temperature for isothermal walls
    switch s.boundary_conditions.boundary_eta0.name
        case 'no_slip_isothermal'
            s.var.p(:, 1) = 2 * s.var.p(:, 2) - s.var.p(:, 3); % Interpolate pressure
    end
end
