function [s] = RESTART_SOLUTION(s, solution_old, chemistry, disturbances)
% RESTART_SOLUTION  Restart a s from a previously saved state.
%   Handles re-meshing with field interpolation when the grid dimensions
%   have changed, or directly restores the old s when no re-mesh
%   is required.
%
%   [s] = RESTART_SOLUTION(s, solution_old, chemistry, disturbances)
%
%   Inputs:
%       s      - Current s struct (target grid settings)
%       solution_old  - Previously saved s struct to restart from
%       chemistry     - Chemistry model struct for thermodynamic evaluations
%       disturbances  - Disturbance settings passed to GENERATE_MESH
%
%   Outputs:
%       s      - Updated s struct with restarted fields
%
%   Notes:
%       - If Nx or Ny differ between old and new solutions, remeshing is
%         triggered and all conserved fields are interpolated onto the new
%         grid using linear interpolation.
%       - When a shock is present, only data from inner (sub-shock) cells
%         is used for interpolation to avoid smearing the discontinuity.
%       - After interpolation, boundary conditions, chemistry, and
%         thermodynamic properties are recomputed on the new mesh.
%       - If no remeshing is needed, the old s is copied directly.
%
% Part of: Hypersonics Stability MATLAB Solver - Initialization Module

    %% Check if remeshing is required
    if s.mesh.Nchi ~= solution_old.mesh.Nchi || s.mesh.Neta ~= solution_old.mesh.Neta
        s.remesh = true;
    end

    if s.restart && s.remesh

        %% Initialize variables on new mesh
        s = VARIABLES_INITIALIZATION(s);
        s.linearize = false;

        %% Precompute chi_wall for shock interpolation (same as GENERATE_MESH)
        chi_wall_temp = linspace(0, 1, s.mesh.Nchi+1)';
        if s.curvilinear_mapping.refinement_stagnation.state
            chi_wall_temp = REFINEMENT(chi_wall_temp, ...
                s.curvilinear_mapping.refinement_stagnation.BL_thickness, ...
                0, s.curvilinear_mapping.refinement_stagnation.intensity, 1, 0);
        end
        s.mesh.chi_wall = chi_wall_temp;

        %% Interpolate shock position to new mesh
        if s.shock.enabled == true
            s.shock.points_chi = (s.mesh.chi_wall(1:end-1) + s.mesh.chi_wall(2:end)) / 2;
            shock_points_chi_temp = (s.mesh.chi_wall(1:end-1) + s.mesh.chi_wall(2:end)) / 2;
            s.shock.points_x = ppval(solution_old.shock.spline_func_x, shock_points_chi_temp);
            s.shock.points_y = ppval(solution_old.shock.spline_func_y, shock_points_chi_temp);
            s.shock.spline_func_x = solution_old.shock.spline_func_x;
            s.shock.spline_func_y = solution_old.shock.spline_func_y;
            [s.shock.points_chi, s.shock.points_eta] = GO_TO_ELEMENT_SPACE( ...
                s.shock.points_x, s.shock.points_y, s);
        end

        %% Carry forward time and generate new mesh
        s.time_integration.t = solution_old.time_integration.t;
        s = GENERATE_MESH(s, disturbances);

        interp_method = "linear";

        %% Interpolate fields without shock (full-domain griddata)
        if ~s.shock.enabled
            s.var.rho   = griddata(solution_old.mesh.x_Ext, solution_old.mesh.y_Ext, solution_old.var.rho, ...
                s.mesh.x_Ext, s.mesh.y_Ext, interp_method);
            s.var.rho_u = griddata(solution_old.mesh.x_Ext, solution_old.mesh.y_Ext, solution_old.var.rho_u, ...
                s.mesh.x_Ext, s.mesh.y_Ext, interp_method);
            s.var.rho_v = griddata(solution_old.mesh.x_Ext, solution_old.mesh.y_Ext, solution_old.var.rho_v, ...
                s.mesh.x_Ext, s.mesh.y_Ext, interp_method);
            s.var.rho_E = griddata(solution_old.mesh.x_Ext, solution_old.mesh.y_Ext, solution_old.var.rho_E, ...
                s.mesh.x_Ext, s.mesh.y_Ext, interp_method);
            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star = griddata(solution_old.mesh.x_Ext, solution_old.mesh.y_Ext, solution_old.var.gamma_star, ...
                    s.mesh.x_Ext, s.mesh.y_Ext, interp_method);
                s.var.cv_star = griddata(solution_old.mesh.x_Ext, solution_old.mesh.y_Ext, solution_old.var.cv_star, ...
                    s.mesh.x_Ext, s.mesh.y_Ext, interp_method);
            end
        end

        %% Interpolate fields with shock (sub-shock cells only)
        % Store only data from inner shock cells to avoid building an
        % interpolant across the discontinuity.
        if s.shock.enabled
            count = 1;
            for i = 1:solution_old.mesh.Nchi + 2
                if i > 1 && i < solution_old.mesh.Nchi + 2
                    jmax = solution_old.shock.cell_indices(i - 1, 1) + 2;
                elseif i == 1
                    jmax = solution_old.shock.cell_indices(1, 1) + 2;
                else
                    jmax = solution_old.shock.cell_indices(solution_old.mesh.Nchi, 1) + 2;
                end

                for j = 1:jmax
                    x_temp(count)      = solution_old.mesh.x_Ext(i, j);
                    y_temp(count)      = solution_old.mesh.y_Ext(i, j);
                    temp_rho(count)    = solution_old.var.rho(i, j);
                    temp_rho_u(count)  = solution_old.var.rho_u(i, j);
                    temp_rho_v(count)  = solution_old.var.rho_v(i, j);
                    temp_rho_E(count)  = solution_old.var.rho_E(i, j);
                    if ~s.chemistry.chemical_equilibrium
                        temp_gamma_star(count) = solution_old.var.gamma_star(i, j);
                        temp_cv_star(count)    = solution_old.var.cv_star(i, j);
                    end
                    count = count + 1;
                end
            end

            %% Fit and evaluate interpolants on new mesh
            rho_I   = fit([x_temp(:), y_temp(:)], temp_rho(:), interp_method);
            s.var.rho = rho_I(s.mesh.x_Ext, s.mesh.y_Ext);

            rho_u_I = fit([x_temp(:), y_temp(:)], temp_rho_u(:), interp_method);
            s.var.rho_u = rho_u_I(s.mesh.x_Ext, s.mesh.y_Ext);

            rho_v_I = fit([x_temp(:), y_temp(:)], temp_rho_v(:), interp_method);
            s.var.rho_v = rho_v_I(s.mesh.x_Ext, s.mesh.y_Ext);

            rho_E_I = fit([x_temp(:), y_temp(:)], temp_rho_E(:), interp_method);
            s.var.rho_E = rho_E_I(s.mesh.x_Ext, s.mesh.y_Ext);

            if ~s.chemistry.chemical_equilibrium
                gamma_star_I = fit([x_temp(:), y_temp(:)], temp_gamma_star(:), interp_method);
                s.var.gamma_star = gamma_star_I(s.mesh.x_Ext, s.mesh.y_Ext);

                cv_star_I = fit([x_temp(:), y_temp(:)], temp_cv_star(:), interp_method);
                s.var.cv_star = cv_star_I(s.mesh.x_Ext, s.mesh.y_Ext);
            end
        end

        %% Update freestream and shock boundary conditions
        s = SET_FREESTREAM_PROPERTIES(s, chemistry);

        if s.shock.enabled == true
            s = UPDATE_FLOW_CELLS(s, chemistry);
            s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry);
            s = UPDATE_SOUND_SPEED(s, chemistry);
            s = UPDATE_SHOCK_BC(s, chemistry);
            s = SET_FREESTREAM_PROPERTIES(s, chemistry);
            s = UPDATE_FIELD_UPSTREAM(s);
        end

        %% Extend upstream flow and finalize boundary conditions
        s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry);
        s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry);
        s = APPLY_BOUNDARY_CONDITIONS(s, chemistry);
        s = UPDATE_SOUND_SPEED(s, chemistry);

    elseif s.restart && ~s.remesh
        %% No remesh: copy old s directly
        s = solution_old;
    end
end
