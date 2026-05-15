function updated_handles = PLOT_RUNNING(s, solution_temp_base, handles, chemistry)
% PLOT_RUNNING  Update a live subplot figure with flow-field variables
%   during a time-marching simulation.
%
%   updated_handles = PLOT_RUNNING(s, solution_temp_base, handles,
%       chemistry)
%
%   Creates (on the first qualifying call) or updates (on subsequent calls)
%   a multi-panel subplot figure that displays user-selected flow-field
%   variables.  Pcolor surfaces and optional shock-line overlays are
%   refreshed every s.running_plot.timesteps iterations, providing a
%   live visual monitor of the simulation state.
%
%   Inputs:
%       s           - (struct) current s state containing
%                            conserved variables, grid coordinates,
%                            iteration counter, plotting flags, and
%                            freestream reference values.
%       solution_temp_base - (struct) baseline s used for computing
%                            difference fields (Drho, Du, Dv, De, DP).
%       handles            - (struct) graphics handles with fields: fig,
%                            ax, plot_object, colorbar, shock_line, and
%                            caxis_limits_per_var.  On the first call these
%                            fields may be empty.
%       chemistry          - (struct) chemistry model structure (used by
%                            UPDATE_THERMODYNAMIC_PROPERTIES for pressure
%                            and temperature variables).
%
%   Outputs:
%       updated_handles - (struct) updated graphics-handle struct for reuse
%                         in subsequent iterations.
%
%   Notes:
%       Plotting occurs only when s.running_plot.enabled is true and the
%       current iteration is a multiple of s.running_plot.timesteps.
%       Supported variable names (set via s.running_plot.variable):
%           "grad(rho)/rho", "rho", "u", "v", "vort", "div_U", "e", "P",
%           "T", "gamma*", "cv*", species ("CO2","CO","C","O","N2","N",
%           "NO","H2","H"), and difference fields ("Drho","Du","Dv","De",
%           "DP").
%       Multiple variables are arranged in a 3-column subplot grid.
%       On the first call, pcolor objects, colorbars, and shock lines are
%       created; on subsequent calls only their data properties are updated
%       (no figure recreation) for efficiency.  Colour limits are fixed
%       after the first call to preserve visual consistency.
%
% Part of: Hypersonics Stability MATLAB Solver - Postprocessing Module

    if (s.running_plot.enabled && 0 == mod(s.time_integration.iter, s.running_plot.timesteps))

        %% Prepare temporary s with upstream perturbations
        solution_temp = s;
        if s.shock.enabled
            solution_temp = UPDATE_FIELD_UPSTREAM(solution_temp);
        end
        N_plot_shock_poly = 1000;

        if iscell(solution_temp.running_plot.variable)
            variables_to_plot = solution_temp.running_plot.variable;
        else
            variables_to_plot = {solution_temp.running_plot.variable};
        end
        num_variables = length(variables_to_plot);

        num_subplot_cols = 3;
        num_subplot_rows = ceil(num_variables / 3);

        %% Create or reuse figure
        if isempty(handles.fig) || ~isvalid(handles.fig)
            handles.fig = figure('Name', 'Running Simulation Variables');
            fig_width   = min(1800, 600 * num_subplot_cols);
            fig_height  = min(1600, 800 * num_subplot_rows);
            screen_size = get(groot, 'ScreenSize');
            pos_x = (screen_size(3) - fig_width) / 2;
            pos_y = (screen_size(4) - fig_height) / 2;
            set(handles.fig, 'Position', [pos_x, pos_y, fig_width, fig_height]);
            first_plot_call = true;

            if ~isfield(handles, 'shock_line')
                handles.shock_line = [];
            end
        else
            first_plot_call = false;
        end

        %% Prepare shock line data
        shock_line_exists = false;
        if isfield(solution_temp, 'shock') && solution_temp.shock.enabled == true && ...
                isfield(solution_temp.shock, 'points_chi') && ...
                isfield(solution_temp.shock, 'spline_func_x') && ...
                isfield(solution_temp.shock, 'spline_func_y')

            try
                chi_p = linspace(solution_temp.shock.points_chi(1, 1), ...
                    solution_temp.shock.points_chi(end, 1), N_plot_shock_poly);
                x_p = ppval(solution_temp.shock.spline_func_x, chi_p);
                y_p = ppval(solution_temp.shock.spline_func_y, chi_p);
                shock_line_exists = true;
            catch
                warning('Could not generate shock line data');
                shock_line_exists = false;
            end
        end

        %% Loop over variables to plot
        for v_idx = 1:num_variables
            current_var_name_str = variables_to_plot{v_idx};

            plot_x_coords = [];
            plot_y_coords = [];
            plot_c_data   = [];
            name_latex    = '';

            %% Extract freestream reference values
            U_infty = s.freestream.U;
            p_infty = s.freestream.p;
            e_infty = s.freestream.e;
            if U_infty == 0, U_infty = 1; end
            if p_infty == 0, p_infty = 1; end
            if e_infty == 0, e_infty = 1; end

            %% Variable-specific data extraction
            switch current_var_name_str
                case "grad(rho)/rho"
                    plot_x_coords = solution_temp.mesh.x(2:end-1, 2:end-1);
                    plot_y_coords = solution_temp.mesh.y(2:end-1, 2:end-1);
                    [d_rho_dx, d_rho_dy] = DERIVATIVE(solution_temp.var.rho(2:end-1, 2:end-1), solution_temp);
                    plot_c_data = sqrt(d_rho_dx.^2 + d_rho_dy.^2) ./ solution_temp.var.rho(3:end-2, 3:end-2);
                    name_latex = "$\nabla \rho / \rho $";

                case "rho"
                    plot_x_coords = solution_temp.mesh.x;
                    plot_y_coords = solution_temp.mesh.y;
                    plot_c_data = solution_temp.var.rho(2:end-1, 2:end-1) ./ solution_temp.freestream.rho_0;
                    name_latex = "$\rho/\rho_\infty$";

                case "u"
                    plot_x_coords = solution_temp.mesh.x;
                    plot_y_coords = solution_temp.mesh.y;
                    plot_c_data = (solution_temp.var.rho_u(2:end-1, 2:end-1) ./ solution_temp.var.rho(2:end-1, 2:end-1));
                    name_latex = "$u/U_\infty$";

                case "v"
                    plot_x_coords = solution_temp.mesh.x;
                    plot_y_coords = solution_temp.mesh.y;
                    plot_c_data = (solution_temp.var.rho_v(2:end-1, 2:end-1) ./ solution_temp.var.rho(2:end-1, 2:end-1));
                    name_latex = "$v/U_\infty$";

                case "vort"
                    plot_x_coords = solution_temp.mesh.x(2:end-1, 2:end-1);
                    plot_y_coords = solution_temp.mesh.y(2:end-1, 2:end-1);
                    u_vel_int = solution_temp.var.rho_u(2:end-1, 2:end-1) ./ solution_temp.var.rho(2:end-1, 2:end-1);
                    v_vel_int = solution_temp.var.rho_v(2:end-1, 2:end-1) ./ solution_temp.var.rho(2:end-1, 2:end-1);
                    [~, d_u_dy] = DERIVATIVE(u_vel_int, solution_temp);
                    [d_v_dx, ~] = DERIVATIVE(v_vel_int, solution_temp);
                    plot_c_data = d_v_dx - d_u_dy;
                    plot_c_data = plot_c_data * solution_temp.freestream.L;
                    name_latex = "$\omega L/U_\infty$";

                case "div_U"
                    plot_x_coords = solution_temp.mesh.x(2:end-1, 2:end-1);
                    plot_y_coords = solution_temp.mesh.y(2:end-1, 2:end-1);
                    u_vel_int = solution_temp.var.rho_u(2:end-1, 2:end-1) ./ solution_temp.var.rho(2:end-1, 2:end-1);
                    v_vel_int = solution_temp.var.rho_v(2:end-1, 2:end-1) ./ solution_temp.var.rho(2:end-1, 2:end-1);
                    [d_u_dx, ~] = DERIVATIVE(u_vel_int, solution_temp);
                    [~, d_v_dy] = DERIVATIVE(v_vel_int, solution_temp);
                    plot_c_data = d_u_dx + d_v_dy;
                    plot_c_data = plot_c_data * solution_temp.freestream.L;
                    name_latex = "$\nabla \cdot \mathbf{u} L/U_\infty$";

                case "e"
                    plot_x_coords = solution_temp.mesh.x;
                    plot_y_coords = solution_temp.mesh.y;
                    rho_int    = solution_temp.var.rho(2:end-1, 2:end-1);
                    u_int      = solution_temp.var.rho_u(2:end-1, 2:end-1) ./ rho_int;
                    v_int      = solution_temp.var.rho_v(2:end-1, 2:end-1) ./ rho_int;
                    internal_e = solution_temp.var.rho_E(2:end-1, 2:end-1) - 0.5 * rho_int .* (u_int.^2 + v_int.^2);
                    plot_c_data = internal_e * energy_factor ./ s.freestream.e;
                    name_latex = "$e/e_\infty$";

                case "P"
                    plot_x_coords = solution_temp.mesh.x;
                    plot_y_coords = solution_temp.mesh.y;
                    solution_temp = UPDATE_THERMODYNAMIC_PROPERTIES(solution_temp, chemistry);
                    plot_c_data = solution_temp.var.p(2:end-1, 2:end-1) / s.freestream.p_0;
                    name_latex = "$P/P_\infty$";

                case "T"
                    plot_x_coords = solution_temp.mesh.x;
                    plot_y_coords = solution_temp.mesh.y;
                    plot_c_data = solution_temp.var.T(2:end-1, 2:end-1) .* s.freestream.energy_factor ./ s.freestream.cv;
                    name_latex = "$T[K]$";

                case "gamma*"
                    if s.chemistry.is_chemistry_enabled
                        plot_x_coords = solution_temp.mesh.x;
                        plot_y_coords = solution_temp.mesh.y;
                        plot_c_data = solution_temp.var.gamma_star(2:end-1, 2:end-1);
                        name_latex = "$\gamma^*$";
                    end

                case "cv*"
                    if s.chemistry.is_chemistry_enabled
                        plot_x_coords = solution_temp.mesh.x;
                        plot_y_coords = solution_temp.mesh.y;
                        plot_c_data = solution_temp.var.cv_star(2:end-1, 2:end-1);
                        name_latex = "$c_v^*$";
                    end

                case "CO2"
                    if s.chemistry.is_chemistry_enabled
                        plot_x_coords = solution_temp.mesh.x;
                        plot_y_coords = solution_temp.mesh.y;
                        e = solution_temp.var.rho_E ./ solution_temp.var.rho ...
                            - 1/2 * (solution_temp.var.rho_u.^2 + solution_temp.var.rho_v.^2) ./ solution_temp.var.rho.^2;
                        plot_c_data = chemistry.eval_log10_CO2( ...
                            s.freestream.rho_factor * solution_temp.var.rho(2:end-1, 2:end-1), ...
                            s.freestream.energy_factor * e(2:end-1, 2:end-1));
                        name_latex = "$log_{10}(X_{CO_2})$";
                    end

                case "CO"
                    if s.chemistry.is_chemistry_enabled
                        plot_x_coords = solution_temp.mesh.x;
                        plot_y_coords = solution_temp.mesh.y;
                        e = solution_temp.var.rho_E ./ solution_temp.var.rho ...
                            - 1/2 * (solution_temp.var.rho_u.^2 + solution_temp.var.rho_v.^2) ./ solution_temp.var.rho.^2;
                        plot_c_data = chemistry.eval_log10_CO( ...
                            s.freestream.rho_factor * solution_temp.var.rho(2:end-1, 2:end-1), ...
                            s.freestream.energy_factor * e(2:end-1, 2:end-1));
                        name_latex = "$log_{10}(X_{CO})$";
                    end

                case "C"
                    if s.chemistry.is_chemistry_enabled
                        plot_x_coords = solution_temp.mesh.x;
                        plot_y_coords = solution_temp.mesh.y;
                        e = solution_temp.var.rho_E ./ solution_temp.var.rho ...
                            - 1/2 * (solution_temp.var.rho_u.^2 + solution_temp.var.rho_v.^2) ./ solution_temp.var.rho.^2;
                        plot_c_data = chemistry.eval_log10_C( ...
                            s.freestream.rho_factor * solution_temp.var.rho(2:end-1, 2:end-1), ...
                            s.freestream.energy_factor * e(2:end-1, 2:end-1));
                        name_latex = "$log_{10}(X_{C})$";
                    end

                case "O"
                    if s.chemistry.is_chemistry_enabled
                        plot_x_coords = solution_temp.mesh.x;
                        plot_y_coords = solution_temp.mesh.y;
                        e = solution_temp.var.rho_E ./ solution_temp.var.rho ...
                            - 1/2 * (solution_temp.var.rho_u.^2 + solution_temp.var.rho_v.^2) ./ solution_temp.var.rho.^2;
                        plot_c_data = chemistry.eval_log10_O( ...
                            s.freestream.rho_factor * solution_temp.var.rho(2:end-1, 2:end-1), ...
                            s.freestream.energy_factor * e(2:end-1, 2:end-1));
                        name_latex = "$log_{10}(X_{O})$";
                    end

                case "O2"
                    if s.chemistry.is_chemistry_enabled
                        plot_x_coords = solution_temp.mesh.x;
                        plot_y_coords = solution_temp.mesh.y;
                        e = solution_temp.var.rho_E ./ solution_temp.var.rho ...
                            - 1/2 * (solution_temp.var.rho_u.^2 + solution_temp.var.rho_v.^2) ./ solution_temp.var.rho.^2;
                        plot_c_data = chemistry.eval_log10_O2( ...
                            s.freestream.rho_factor * solution_temp.var.rho(2:end-1, 2:end-1), ...
                            s.freestream.energy_factor * e(2:end-1, 2:end-1));
                        name_latex = "$log_{10}(X_{O_2})$";
                    end

                case "N2"
                    if s.chemistry.is_chemistry_enabled
                        plot_x_coords = solution_temp.mesh.x;
                        plot_y_coords = solution_temp.mesh.y;
                        e = solution_temp.var.rho_E ./ solution_temp.var.rho ...
                            - 1/2 * (solution_temp.var.rho_u.^2 + solution_temp.var.rho_v.^2) ./ solution_temp.var.rho.^2;
                        plot_c_data = chemistry.eval_log10_N2( ...
                            s.freestream.rho_factor * solution_temp.var.rho(2:end-1, 2:end-1), ...
                            s.freestream.energy_factor * e(2:end-1, 2:end-1));
                        name_latex = "$log_{10}(X_{N_2})$";
                    end

                case "N"
                    if s.chemistry.is_chemistry_enabled
                        plot_x_coords = solution_temp.mesh.x;
                        plot_y_coords = solution_temp.mesh.y;
                        e = solution_temp.var.rho_E ./ solution_temp.var.rho ...
                            - 1/2 * (solution_temp.var.rho_u.^2 + solution_temp.var.rho_v.^2) ./ solution_temp.var.rho.^2;
                        plot_c_data = chemistry.eval_log10_N( ...
                            s.freestream.rho_factor * solution_temp.var.rho(2:end-1, 2:end-1), ...
                            s.freestream.energy_factor * e(2:end-1, 2:end-1));
                        name_latex = "$log_{10}(X_{N})$";
                    end

                case "NO"
                    if s.chemistry.is_chemistry_enabled
                        plot_x_coords = solution_temp.mesh.x;
                        plot_y_coords = solution_temp.mesh.y;
                        e = solution_temp.var.rho_E ./ solution_temp.var.rho ...
                            - 1/2 * (solution_temp.var.rho_u.^2 + solution_temp.var.rho_v.^2) ./ solution_temp.var.rho.^2;
                        plot_c_data = chemistry.eval_log10_NO( ...
                            s.freestream.rho_factor * solution_temp.var.rho(2:end-1, 2:end-1), ...
                            s.freestream.energy_factor * e(2:end-1, 2:end-1));
                        name_latex = "$log_{10}(X_{NO})$";
                    end

                case "H2"
                    if s.chemistry.is_chemistry_enabled
                        plot_x_coords = solution_temp.mesh.x;
                        plot_y_coords = solution_temp.mesh.y;
                        e = solution_temp.var.rho_E ./ solution_temp.var.rho ...
                            - 1/2 * (solution_temp.var.rho_u.^2 + solution_temp.var.rho_v.^2) ./ solution_temp.var.rho.^2;
                        plot_c_data = chemistry.eval_log10_H2( ...
                            s.freestream.rho_factor * solution_temp.var.rho(2:end-1, 2:end-1), ...
                            s.freestream.energy_factor * e(2:end-1, 2:end-1));
                        name_latex = "$log_{10}(X_{H_2})$";
                    end

                case "H"
                    if s.chemistry.is_chemistry_enabled
                        plot_x_coords = solution_temp.mesh.x;
                        plot_y_coords = solution_temp.mesh.y;
                        e = solution_temp.var.rho_E ./ solution_temp.var.rho ...
                            - 1/2 * (solution_temp.var.rho_u.^2 + solution_temp.var.rho_v.^2) ./ solution_temp.var.rho.^2;
                        plot_c_data = chemistry.eval_log10_H( ...
                            s.freestream.rho_factor * solution_temp.var.rho(2:end-1, 2:end-1), ...
                            s.freestream.energy_factor * e(2:end-1, 2:end-1));
                        name_latex = "$log_{10}(X_{H})$";
                    end

                %% Difference fields (current minus baseline)
                case "Drho"
                    plot_x_coords = solution_temp.mesh.x;
                    plot_y_coords = solution_temp.mesh.y;
                    rho_curr = solution_temp.var.rho(2:end-1, 2:end-1) ./ solution_temp.freestream.rho_0;
                    rho_base = solution_temp_base.rho(2:end-1, 2:end-1) ./ solution_temp_base.freestream.rho_0;
                    plot_c_data = rho_curr - rho_base;
                    name_latex = "$\Delta (\rho/\rho_\infty)$";

                case "Du"
                    plot_x_coords = solution_temp.mesh.x;
                    plot_y_coords = solution_temp.mesh.y;
                    u_curr = (solution_temp.var.rho_u(2:end-1, 2:end-1) ./ solution_temp.var.rho(2:end-1, 2:end-1)) / U_infty;
                    u_base = (solution_temp_base.rho_u(2:end-1, 2:end-1) ./ solution_temp_base.rho(2:end-1, 2:end-1)) / U_infty;
                    plot_c_data = u_curr - u_base;
                    name_latex = "$\Delta (u/U_\infty)$";

                case "Dv"
                    plot_x_coords = solution_temp.mesh.x;
                    plot_y_coords = solution_temp.mesh.y;
                    v_curr = (solution_temp.var.rho_v(2:end-1, 2:end-1) ./ solution_temp.var.rho(2:end-1, 2:end-1)) / U_infty;
                    v_base = (solution_temp_base.rho_v(2:end-1, 2:end-1) ./ solution_temp_base.rho(2:end-1, 2:end-1)) / U_infty;
                    plot_c_data = v_curr - v_base;
                    name_latex = "$\Delta (v/U_\infty)$";

                case "De"
                    plot_x_coords = solution_temp.mesh.x;
                    plot_y_coords = solution_temp.mesh.y;
                    rho_int_curr    = solution_temp.var.rho(2:end-1, 2:end-1);
                    u_int_curr      = solution_temp.var.rho_u(2:end-1, 2:end-1) ./ rho_int_curr;
                    v_int_curr      = solution_temp.var.rho_v(2:end-1, 2:end-1) ./ rho_int_curr;
                    e_curr_internal = solution_temp.var.rho_E(2:end-1, 2:end-1) ...
                        - 0.5 * rho_int_curr .* (u_int_curr.^2 + v_int_curr.^2);
                    e_curr_norm     = e_curr_internal ./ (solution_temp.freestream.rho_0 * e_infty);

                    rho_int_base    = solution_temp_base.rho(2:end-1, 2:end-1);
                    u_int_base      = solution_temp_base.rho_u(2:end-1, 2:end-1) ./ rho_int_base;
                    v_int_base      = solution_temp_base.rho_v(2:end-1, 2:end-1) ./ rho_int_base;
                    e_base_internal = solution_temp_base.rho_E(2:end-1, 2:end-1) ...
                        - 0.5 * rho_int_base .* (u_int_base.^2 + v_int_base.^2);
                    e_base_norm     = e_base_internal ./ (solution_temp_base.freestream.rho_0 * e_infty);

                    plot_c_data = e_curr_norm - e_base_norm;
                    name_latex = "$\Delta (e/e_\infty)$";

                case "DP"
                    plot_x_coords = solution_temp.mesh.x;
                    plot_y_coords = solution_temp.mesh.y;
                    solution_temp      = UPDATE_THERMODYNAMIC_PROPERTIES(solution_temp, chemistry);
                    solution_temp_base = UPDATE_THERMODYNAMIC_PROPERTIES(solution_temp_base, chemistry);
                    p_curr_norm = solution_temp.var.p(2:end-1, 2:end-1) / p_infty;
                    p_base_norm = solution_temp_base.p(2:end-1, 2:end-1) / p_infty;
                    plot_c_data = p_curr_norm - p_base_norm;
                    name_latex = "$\Delta (P/P_\infty)$";

                case "residuals"
                    if ~s.time_integration.get_residual
                        warning('Data for subplot %s could not be fully generated or is empty. Skipping this subplot.', current_var_name_str);
                        continue
                    end
                    iter = s.time_integration.iter;
                    if first_plot_call
                        ax_r = subplot(num_subplot_rows, num_subplot_cols, v_idx, 'Parent', handles.fig);
                        handles.residual(1) = semilogy(ax_r, iter, s.time_integration.residual.rho, '.-', 'LineWidth', 2, 'DisplayName', '$\rho$');
                        hold(ax_r,'on');
                        handles.residual(2) = semilogy(ax_r, iter, s.time_integration.residual.rho_u, '.-', 'LineWidth', 2, 'DisplayName', '$\rho u$');
                        handles.residual(3) = semilogy(ax_r, iter, s.time_integration.residual.rho_v, '.-', 'LineWidth', 2, 'DisplayName', '$\rho v$');
                        handles.residual(4) = semilogy(ax_r, iter, s.time_integration.residual.rho_E, '.-', 'LineWidth', 2, 'DisplayName', '$\rho E$');
                        hold(ax_r,'off');
                        legend(ax_r,'Interpreter','latex');
                        xlabel(ax_r,'Iteration');
                        title(ax_r,'RMS residuals (post-shock cells)');
                        grid(ax_r,'on');
                        handles.ax(v_idx) = ax_r;
                    else
                        set(handles.residual(1), ...
                            'XData',[handles.residual(1).XData iter], ...
                            'YData',[handles.residual(1).YData s.time_integration.residual.rho]);
                        set(handles.residual(2), ...
                            'XData',[handles.residual(2).XData iter], ...
                            'YData',[handles.residual(2).YData s.time_integration.residual.rho_u]);
                        set(handles.residual(3), ...
                            'XData',[handles.residual(3).XData iter], ...
                            'YData',[handles.residual(3).YData s.time_integration.residual.rho_v]);
                        set(handles.residual(4), ...
                            'XData',[handles.residual(4).XData iter], ...
                            'YData',[handles.residual(4).YData s.time_integration.residual.rho_E]);
                        set(handles.ax(v_idx),'XLimMode','auto','YLimMode','auto');
                    end

                otherwise
                    warning('Unknown variable selected for plotting: %s', current_var_name_str);
                    plot_x_coords = [];
            end

            if ~strcmp(current_var_name_str,'residuals') % Residual plotting is handled separately
                %% Skip if data could not be generated
                if isempty(plot_x_coords) || isempty(plot_y_coords) || isempty(plot_c_data)
                    warning('Data for subplot %s could not be fully generated or is empty. Skipping this subplot.', current_var_name_str);
                    if first_plot_call && num_variables == 1 && isfield(handles, 'fig') && isvalid(handles.fig)
                        try close(handles.fig); catch; end
                        handles.fig = [];
                    end
                    continue;
                end
    
                %% Dimension check before plotting
                if ~isequal(size(plot_x_coords), size(plot_y_coords), size(plot_c_data))
                    warning('Dimension mismatch for %s: X=[%s], Y=[%s], C=[%s]. Skipping subplot.', ...
                        current_var_name_str, mat2str(size(plot_x_coords)), ...
                        mat2str(size(plot_y_coords)), mat2str(size(plot_c_data)));
                    if first_plot_call && num_variables == 1 && isfield(handles, 'fig') && isvalid(handles.fig)
                        try close(handles.fig); catch; end
                        handles.fig = [];
                    end
                    continue;
                end
    
                %% Create or update subplot
                if first_plot_call
                    handles.ax(v_idx) = subplot(num_subplot_rows, num_subplot_cols, v_idx, 'Parent', handles.fig);
                    handles.plot_object(v_idx) = pcolor(handles.ax(v_idx), plot_x_coords, plot_y_coords, plot_c_data);
                    shading(handles.ax(v_idx), 'interp');
                    set(handles.plot_object(v_idx), 'EdgeColor', 'none');
                    colormap(handles.ax(v_idx), "jet");
                    handles.colorbar(v_idx) = colorbar(handles.ax(v_idx));
                    handles.colorbar(v_idx).Label.Interpreter = 'latex';
                    handles.colorbar(v_idx).Label.String      = name_latex;
                    handles.colorbar(v_idx).Label.FontSize    = 10;
    
                    min_val = min(plot_c_data(:));
                    max_val = max(plot_c_data(:));
                    if min_val == max_val
                        min_val = min_val * 2;
                        max_val = max_val * 2;
                        if min_val == 0 && max_val == 0
                            min_val = -0.01;
                            max_val = 0.01;
                        end
                    end
                    if isempty(min_val) || isnan(min_val) || isinf(min_val) || ...
                            isempty(max_val) || isnan(max_val) || isinf(max_val)
                        handles.caxis_limits_per_var{v_idx} = [];
                    else
                        handles.caxis_limits_per_var{v_idx} = [min_val, max_val];
                    end
    
                    if ~isempty(handles.caxis_limits_per_var{v_idx})
                        caxis(handles.ax(v_idx), handles.caxis_limits_per_var{v_idx});
                    end
    
                    title(handles.ax(v_idx), name_latex, 'Interpreter', 'latex', 'FontSize', 11);
                    xlabel(handles.ax(v_idx), '$x/L$', 'Interpreter', 'latex', 'FontSize', 9);
                    ylabel(handles.ax(v_idx), '$y/L$', 'Interpreter', 'latex', 'FontSize', 9);
                    axis(handles.ax(v_idx), 'equal');
                    axis(handles.ax(v_idx), 'tight');
    
                    % Add shock line overlay
                    if shock_line_exists
                        hold(handles.ax(v_idx), 'on');
                        handles.shock_line(v_idx) = plot(handles.ax(v_idx), x_p, y_p, 'black', 'LineWidth', 2);
                        hold(handles.ax(v_idx), 'off');
                    end
    
                else
                    %% Update existing plot data
                    set(handles.plot_object(v_idx), 'XData', plot_x_coords, 'YData', plot_y_coords, 'CData', plot_c_data);
                    if ~isempty(handles.caxis_limits_per_var{v_idx})
                        caxis(handles.ax(v_idx), handles.caxis_limits_per_var{v_idx});
                    else
                        caxis(handles.ax(v_idx), 'auto');
                    end
    
                    % Update shock line
                    if shock_line_exists
                        set(handles.shock_line(v_idx), 'XData', x_p, 'YData', y_p);
                    elseif shock_line_exists && (length(handles.shock_line) < v_idx || ...
                            isempty(handles.shock_line(v_idx)) || ~isvalid(handles.shock_line(v_idx)))
                        hold(handles.ax(v_idx), 'on');
                        handles.shock_line(v_idx) = plot(handles.ax(v_idx), x_p, y_p, 'black', 'LineWidth', 2);
                        hold(handles.ax(v_idx), 'off');
                    end
                end
            end
        end

        %% Update super-title with iteration and time info
        if ~isempty(handles.fig) && isvalid(handles.fig) && any(isgraphics(handles.ax))
            sgtitle(handles.fig, sprintf('Iteration: %d, Time: %.3f', ...
                solution_temp.time_integration.iter, solution_temp.time_integration.t), 'FontSize', 14, 'FontWeight', 'bold');
        end

        drawnow;

        % Option to each generated plot if 'plot_dir' is given
        if isfield(s,'plot_dir')
            savefig(handles.fig, sprintf("%s/runningplt_it%d.fig",s.plot_dir,s.time_integration.iter));
        end
    end
    updated_handles = handles;
end
