function s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry)
% UPDATE_CHEMISTRY_EQUILIBRIUM  Evaluate equilibrium thermochemical and transport properties.
%
%   s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry)
%
%   Computes the equilibrium values of gamma_star, cv_star, dynamic
%   viscosity (mu), and thermal conductivity (k) from the chemistry lookup
%   tables. When chemical equilibrium is active, the s properties
%   are set to their equilibrium values; otherwise, only the _eq fields
%   are updated. For non-chemistry cases, uniform freestream values are
%   assigned. Derived quantities (Re_flow, Pr_flow) are also updated.
%
%   Inputs:
%       s  - Solution struct containing:
%                   .var.rho, .var.rho_u, .var.rho_v, .var.rho_E - Conserved variables
%                   .chemistry.is_chemistry_enabled - Logical, true if chemistry active
%                   .chemistry.chemical_equilibrium - Logical, true if in equilibrium mode
%                   .freestream.*            - Freestream reference values and
%                                              scaling factors
%       chemistry - Chemistry model struct with evaluation functions:
%                   .eval_gamma_star, .eval_cv_star, .eval_mu, .eval_k
%
%   Outputs:
%       s  - Updated s struct with fields:
%                   .var.gamma_star_eq, .var.gamma_star - Equilibrium ratio of specific heats
%                   .var.cv_star_eq, .var.cv_star       - Equilibrium specific heat at const. vol.
%                   .var.mu_star                        - Dynamic viscosity
%                   .var.Re_flow                        - Local Reynolds number field
%                   .var.k_star                         - Thermal conductivity
%                   .var.Pr_flow                        - Local Prandtl number field
%
%   Notes:
%       - Internal energy is computed as e = (rho_E - KE) / rho.
%       - Dimensional conversions use freestream scaling factors.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% Compute internal energy from conserved variables
    rho = s.var.rho;
    e = (s.var.rho_E - (s.var.rho_u.^2 + s.var.rho_v.^2) ./ rho / 2) ./ rho;

    %% Gamma_star (ratio of specific heats)
    if s.chemistry.is_chemistry_enabled
        s.var.gamma_star_eq = chemistry.eval_gamma_star(s.freestream.rho_factor * rho, s.freestream.energy_factor * e);
        if s.chemistry.chemical_equilibrium
            s.var.gamma_star = s.var.gamma_star_eq;
        end
    else
        s.var.gamma_star_eq = s.freestream.gamma_star * ones(size(rho));
        s.var.gamma_star = s.var.gamma_star_eq;
    end

    %% cv_star (specific heat at constant volume)
    if s.chemistry.is_chemistry_enabled
        cv_star = chemistry.eval_cv_star(s.freestream.rho_factor * rho, s.freestream.energy_factor * e); % Dimensional specific heat at constant volume
        s.var.cv_star_eq = cv_star ./ s.freestream.cv;
        if s.chemistry.chemical_equilibrium
            s.var.cv_star = s.var.cv_star_eq;
        end
    else
        cv = s.freestream.cv * ones(size(rho)); % Dimensional specific heat at constant volume
        s.var.cv_star_eq = cv ./ cv;
        s.var.cv_star = s.var.cv_star_eq;
    end

    %% Dynamic viscosity (mu) and local Reynolds number
    if s.chemistry.is_chemistry_enabled
        s.var.mu_star = chemistry.eval_mu(s.freestream.rho_factor * rho, s.freestream.energy_factor * e) ./ s.freestream.mu;
        s.var.Re_flow = s.freestream.Re ./ s.var.mu_star;
    else
        s.var.mu_star = s.freestream.mu ./ s.freestream.mu * ones(size(rho));
        s.var.Re_flow = s.freestream.Re * ones(size(rho));
    end

    %% Thermal conductivity (k) and local Prandtl number
    if s.chemistry.is_chemistry_enabled
        s.var.k_star = chemistry.eval_k(s.freestream.rho_factor * rho, s.freestream.energy_factor * e) ./ s.freestream.k;
        s.var.Pr_flow =  s.var.mu_star .* s.freestream.mu .* s.var.cv_star .* s.freestream.cv .* s.var.gamma_star ./ s.var.k_star / s.freestream.k;
    else
        s.var.k_star = s.freestream.k / s.freestream.k * ones(size(rho));
        s.var.Pr_flow = s.freestream.Pr * ones(size(rho));
    end
end
