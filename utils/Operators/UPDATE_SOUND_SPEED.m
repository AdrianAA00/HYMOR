function s = UPDATE_SOUND_SPEED(s, chemistry)
% UPDATE_SOUND_SPEED  Compute the local speed of sound from the current flow state.
%
%   s = UPDATE_SOUND_SPEED(s, chemistry)
%
%   Evaluates the speed of sound either from chemistry lookup tables
%   (when chemistry is active) or from the ideal-gas relation
%   a = sqrt(gamma* * p / rho) for a calorically perfect gas.
%
%   Inputs:
%       s  - Solution struct containing:
%                   .var.rho, .var.rho_u, .var.rho_v, .var.rho_E - Conserved variables
%                   .chemistry.is_chemistry_enabled - Logical, true if chemistry active
%                   .var.gamma_star      - Ratio of specific heats field
%                   .var.p               - Pressure field
%                   .freestream.*        - Scaling factors (rho_factor,
%                                          energy_factor, velocity_factor)
%       chemistry - Chemistry model struct with:
%                   .eval_a - Speed of sound evaluation function(rho, e)
%
%   Outputs:
%       s  - Updated s struct with:
%                   .var.a - Speed of sound field (non-dimensional)
%
%   Notes:
%       - Internal energy is computed as e = (rho_E - KE) / rho.
%       - For chemistry cases, the dimensional result from eval_a is
%         non-dimensionalized by velocity_factor.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% Compute internal energy
    rho = s.var.rho;
    e = (s.var.rho_E - (s.var.rho_u.^2 + s.var.rho_v.^2) ./ rho / 2) ./ rho;

    %% Update sound speed
    if s.chemistry.is_chemistry_enabled
        s.var.a = chemistry.eval_a(s.freestream.rho_factor * rho, s.freestream.energy_factor * e) / s.freestream.velocity_factor;
    else
        s.var.a = sqrt(s.var.gamma_star .* s.var.p ./ s.var.rho) / s.freestream.velocity_factor;
    end
end
