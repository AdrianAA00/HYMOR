function [solution_temp, residual] = INTEGRATION_IMPLICIT(s, solution_temp)
% INTEGRATION_IMPLICIT  Perform one implicit Euler update with relaxation.
%
%   [solution_temp, residual] = INTEGRATION_IMPLICIT(s, solution_temp)
%
%   Computes a forward-Euler predictor for each conserved variable (rho,
%   rho_u, rho_v, rho_E), evaluates the maximum absolute difference
%   (residual) between the current iterate and the predictor, and then
%   blends the two using the adaptive relaxation factor stored in
%   s.relax_factor_variable. The shock position is also updated
%   with an explicit Euler step.
%
%   Inputs:
%       s      - struct : Reference (time-level n) flow state
%                                 containing conservative variables (s.var.*),
%                                 dt (s.time_integration.dt), flow_cells mask
%                                 (s.shock.flow_cells), shock data (s.shock.*),
%                                 and the current relaxation factor
%                                 (s.relax_factor_variable).
%       solution_temp - struct : Current implicit iterate whose fluxes
%                                 have already been evaluated by PDE.
%
%   Outputs:
%       solution_temp - struct : Updated iterate after relaxation blending
%                                 for rho, rho_u, rho_v, rho_E and after
%                                 an explicit shock-position update.
%       residual      - double : Maximum absolute residual across all four
%                                 conserved variables (excluding shock and
%                                 ghost cells).
%
%   Notes:
%       - For shock-fitted runs the fluxes are masked by
%         s.shock.flow_cells and the residual is zeroed above the shock
%         row to avoid polluting the convergence measure.
%       - The relaxation update is:
%           q^{k+1} = relax * q^{k} + (1 - relax) * q_predictor
%         where relax = s.relax_factor_variable.
%       - Shock position is advanced with an explicit Euler step scaled
%         by s.shock.relaxation.
%
% Part of: Hypersonics Stability MATLAB Solver - Time Marching Module

    %% Density (rho) update
    if s.shock.enabled == true
        temp = s.var.rho(2:end-1,2:end-1) + s.time_integration.dt * solution_temp.flux.rho .* s.shock.flow_cells;
    else
        temp = s.var.rho(2:end-1,2:end-1) + s.time_integration.dt * solution_temp.flux.rho;
    end
    residual_temp = abs(solution_temp.var.rho(2:end-1,3:end-2) - temp(:,2:end-1));
    if s.shock.enabled == true
        for i = 1:s.mesh.Nchi
            residual_temp(i, s.shock.cell_indices(i,1)-1:end) = 0;
        end
    end
    res_rho = max(residual_temp, [], "all");
    solution_temp.var.rho(2:end-1,2:end-1) = solution_temp.var.rho(2:end-1,2:end-1) * s.relax_factor_variable ...
                                        + (1 - s.relax_factor_variable) * temp;

    %% X-momentum (rho_u) update
    if s.shock.enabled == true
        temp = s.var.rho_u(2:end-1,2:end-1) + s.time_integration.dt * solution_temp.flux.rho_u .* s.shock.flow_cells;
    else
        temp = s.var.rho_u(2:end-1,2:end-1) + s.time_integration.dt * solution_temp.flux.rho_u;
    end
    residual_temp = abs(solution_temp.var.rho_u(2:end-1,3:end-2) - temp(:,2:end-1));
    if s.shock.enabled == true
        for i = 1:s.mesh.Nchi
            residual_temp(i, s.shock.cell_indices(i,1)-1:end) = 0;
        end
    end
    res_rho_u = max(residual_temp, [], "all");
    solution_temp.var.rho_u(2:end-1,2:end-1) = solution_temp.var.rho_u(2:end-1,2:end-1) * s.relax_factor_variable ...
                                           + (1 - s.relax_factor_variable) * temp;

    %% Y-momentum (rho_v) update
    if s.shock.enabled == true
        temp = s.var.rho_v(2:end-1,2:end-1) + s.time_integration.dt * solution_temp.flux.rho_v .* s.shock.flow_cells;
    else
        temp = s.var.rho_v(2:end-1,2:end-1) + s.time_integration.dt * solution_temp.flux.rho_v;
    end
    residual_temp = abs(solution_temp.var.rho_v(2:end-1,3:end-2) - temp(:,2:end-1));
    if s.shock.enabled == true
        for i = 1:s.mesh.Nchi
            residual_temp(i, s.shock.cell_indices(i,1)-1:end) = 0;
        end
    end
    res_rho_v = max(residual_temp, [], "all");
    solution_temp.var.rho_v(2:end-1,2:end-1) = solution_temp.var.rho_v(2:end-1,2:end-1) * s.relax_factor_variable ...
                                           + (1 - s.relax_factor_variable) * temp;

    %% Total energy (rho_E) update
    if s.shock.enabled == true
        temp = s.var.rho_E(2:end-1,2:end-1) + s.time_integration.dt * solution_temp.flux.rho_E .* s.shock.flow_cells;
    else
        temp = s.var.rho_E(2:end-1,2:end-1) + s.time_integration.dt * solution_temp.flux.rho_E;
    end
    residual_temp = abs(solution_temp.var.rho_E(2:end-1,3:end-2) - temp(:,2:end-1));
    if s.shock.enabled == true
        for i = 1:s.mesh.Nchi
            residual_temp(i, s.shock.cell_indices(i,1)-1:end) = 0;
        end
    end
    res_rho_E = max(residual_temp, [], "all");
    solution_temp.var.rho_E(2:end-1,2:end-1) = solution_temp.var.rho_E(2:end-1,2:end-1) * s.relax_factor_variable ...
                                           + (1 - s.relax_factor_variable) * temp;

    %% Overall residual
    residual = max([res_rho, res_rho_u, res_rho_v, res_rho_E]);

    %% Update shock position (explicit Euler step)
    solution_temp.shock.points_x = s.shock.points_x + solution_temp.shock.speed_x * s.time_integration.dt * s.shock.relaxation;
    solution_temp.shock.points_y = s.shock.points_y + solution_temp.shock.speed_y * s.time_integration.dt * s.shock.relaxation;
end
