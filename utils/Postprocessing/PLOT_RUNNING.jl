using Plots, Printf, LinearAlgebra, LoopVectorization

# PLOT_RUNNING  Update a live subplot figure with flow-field variables
#   during a time-marching simulation.
#
#   updated_handles = PLOT_RUNNING(s, solution_temp_base, handles,
#       chemistry)
#
#   Creates (on the first qualifying call) or updates (on subsequent calls)
#   a multi-panel subplot figure that displays user-selected flow-field
#   variables.  Pcolor surfaces and optional shock-line overlays are
#   refreshed every s["running_plot"]["timesteps"] iterations, providing a
#   live visual monitor of the simulation state.
#
#   Inputs:
#       s                  - (Dict) current solution state containing
#                            conserved variables, grid coordinates,
#                            iteration counter, plotting flags, and
#                            freestream reference values.
#       solution_temp_base - (Dict) baseline solution used for computing
#                            difference fields (Drho, Du, Dv, De, DP).
#       handles            - (Dict) graphics handles with fields: fig,
#                            ax, plot_object, colorbar, shock_line, and
#                            caxis_limits_per_var. On the first call these
#                            fields may be empty/nothing.
#       chemistry          - (Dict) chemistry model structure (used by
#                            UPDATE_THERMODYNAMIC_PROPERTIES for pressure
#                            and temperature variables).
#
#   Outputs:
#       updated_handles - (Dict) updated graphics-handle Dict for reuse
#                         in subsequent iterations.
#
#   Notes:
#       Plotting occurs only when s["running_plot"]["enabled"] is true and the
#       current iteration is a multiple of s["running_plot"]["timesteps"].
#       Supported variable names (set via s["running_plot"]["variable"]):
#           "grad(rho)/rho", "rho", "u", "v", "vort", "div_U", "e", "P",
#           "T", "gamma*", "cv*", species ("CO2","CO","C","O","N2","N",
#           "NO","H2","H"), and difference fields ("Drho","Du","Dv","De",
#           "DP").
#       Multiple variables are arranged in a 3-column subplot grid.
#       On the first call, heatmap objects, colorbars, and shock lines are
#       created; on subsequent calls only their data properties are updated
#       (no figure recreation) for efficiency. Colour limits are fixed
#       after the first call to preserve visual consistency.
#
# Part of: Hypersonics Stability MATLAB Solver - Postprocessing Module

function PLOT_RUNNING(s::Dict{String, Any}, solution_temp_base::Dict{String, Any},
                      handles::Dict{String, Any}, chemistry::Dict{String, Any})

    if s["running_plot"]["enabled"] && mod(s["time_integration"]["iter"], s["running_plot"]["timesteps"]) == 0
        SET_PLOT_DEFAULTS()

        ## Prepare temporary solution with upstream perturbations
        solution_temp = deepcopy(s)
        if s["shock"]["enabled"]
            solution_temp = UPDATE_FIELD_UPSTREAM(solution_temp)
        end
        N_plot_shock_poly = 1000

        if isa(solution_temp["running_plot"]["variable"], Vector)
            variables_to_plot = solution_temp["running_plot"]["variable"]
        else
            variables_to_plot = [solution_temp["running_plot"]["variable"]]
        end
        num_variables = length(variables_to_plot)

        num_subplot_cols = 3
        num_subplot_rows = ceil(Int, num_variables / 3)

        ## Create or reuse figure
        first_plot_call = false
        if !haskey(handles, "fig") || handles["fig"] === nothing
            handles["fig"] = true  # Flag that figure exists
            first_plot_call = true

            if !haskey(handles, "shock_line")
                handles["shock_line"] = Dict{Int, Any}()
            end
            if !haskey(handles, "plots")
                handles["plots"] = Dict{Int, Any}()
            end
            if !haskey(handles, "caxis_limits_per_var")
                handles["caxis_limits_per_var"] = Dict{Int, Any}()
            end
        end

        ## Prepare shock line data
        shock_line_exists = false
        x_p = Float64[]
        y_p = Float64[]
        if haskey(solution_temp, "shock") && solution_temp["shock"]["enabled"] == true &&
                haskey(solution_temp["shock"], "points_chi") &&
                haskey(solution_temp["shock"], "spline_func_x") &&
                haskey(solution_temp["shock"], "spline_func_y")

            try
                chi_p = collect(range(solution_temp["shock"]["points_chi"][1, 1],
                                      solution_temp["shock"]["points_chi"][end, 1],
                                      length=N_plot_shock_poly))
                x_p = ppval(solution_temp["shock"]["spline_func_x"], chi_p)
                y_p = ppval(solution_temp["shock"]["spline_func_y"], chi_p)
                shock_line_exists = true
            catch e
                @warn "Could not generate shock line data"
                shock_line_exists = false
            end
        end

        ## Build subplot array
        subplot_plots = []

        ## Loop over variables to plot
        for v_idx in 1:num_variables
            current_var_name_str = variables_to_plot[v_idx]

            plot_x_coords = Matrix{Float64}(undef, 0, 0)
            plot_y_coords = Matrix{Float64}(undef, 0, 0)
            plot_c_data   = Matrix{Float64}(undef, 0, 0)
            name_latex    = ""

            ## Extract freestream reference values
            U_infty = s["freestream"]["U"]
            p_infty = s["freestream"]["p"]
            e_infty = s["freestream"]["e"]
            if U_infty == 0; U_infty = 1; end
            if p_infty == 0; p_infty = 1; end
            if e_infty == 0; e_infty = 1; end

            ## Variable-specific data extraction
            if current_var_name_str == "grad(rho)/rho"
                plot_x_coords = solution_temp["mesh"]["x"][2:end-1, 2:end-1]
                plot_y_coords = solution_temp["mesh"]["y"][2:end-1, 2:end-1]
                (d_rho_dx, d_rho_dy) = DERIVATIVE(solution_temp["var"]["rho"][2:end-1, 2:end-1], solution_temp)
                plot_c_data = sqrt.(d_rho_dx.^2 .+ d_rho_dy.^2) ./ solution_temp["var"]["rho"][3:end-2, 3:end-2]
                name_latex = L"\nabla \rho / \rho"

            elseif current_var_name_str == "rho"
                plot_x_coords = solution_temp["mesh"]["x"]
                plot_y_coords = solution_temp["mesh"]["y"]
                plot_c_data = solution_temp["var"]["rho"][2:end-1, 2:end-1] ./ solution_temp["freestream"]["rho_0"]
                name_latex = L"\rho/\rho_\infty"

            elseif current_var_name_str == "u"
                plot_x_coords = solution_temp["mesh"]["x"]
                plot_y_coords = solution_temp["mesh"]["y"]
                plot_c_data = solution_temp["var"]["rho_u"][2:end-1, 2:end-1] ./ solution_temp["var"]["rho"][2:end-1, 2:end-1]
                name_latex = L"u/U_\infty"

            elseif current_var_name_str == "v"
                plot_x_coords = solution_temp["mesh"]["x"]
                plot_y_coords = solution_temp["mesh"]["y"]
                plot_c_data = solution_temp["var"]["rho_v"][2:end-1, 2:end-1] ./ solution_temp["var"]["rho"][2:end-1, 2:end-1]
                name_latex = L"v/U_\infty"

            elseif current_var_name_str == "vort"
                plot_x_coords = solution_temp["mesh"]["x"][2:end-1, 2:end-1]
                plot_y_coords = solution_temp["mesh"]["y"][2:end-1, 2:end-1]
                u_vel_int = solution_temp["var"]["rho_u"][2:end-1, 2:end-1] ./ solution_temp["var"]["rho"][2:end-1, 2:end-1]
                v_vel_int = solution_temp["var"]["rho_v"][2:end-1, 2:end-1] ./ solution_temp["var"]["rho"][2:end-1, 2:end-1]
                (_, d_u_dy) = DERIVATIVE(u_vel_int, solution_temp)
                (d_v_dx, _) = DERIVATIVE(v_vel_int, solution_temp)
                plot_c_data = d_v_dx .- d_u_dy
                plot_c_data = plot_c_data .* solution_temp["freestream"]["L"]
                name_latex = L"\omega L/U_\infty"

            elseif current_var_name_str == "div_U"
                plot_x_coords = solution_temp["mesh"]["x"][2:end-1, 2:end-1]
                plot_y_coords = solution_temp["mesh"]["y"][2:end-1, 2:end-1]
                u_vel_int = solution_temp["var"]["rho_u"][2:end-1, 2:end-1] ./ solution_temp["var"]["rho"][2:end-1, 2:end-1]
                v_vel_int = solution_temp["var"]["rho_v"][2:end-1, 2:end-1] ./ solution_temp["var"]["rho"][2:end-1, 2:end-1]
                (d_u_dx, _) = DERIVATIVE(u_vel_int, solution_temp)
                (_, d_v_dy) = DERIVATIVE(v_vel_int, solution_temp)
                plot_c_data = d_u_dx .+ d_v_dy
                plot_c_data = plot_c_data .* solution_temp["freestream"]["L"]
                name_latex = L"\nabla \cdot \mathbf{u} L/U_\infty"

            elseif current_var_name_str == "e"
                plot_x_coords = solution_temp["mesh"]["x"]
                plot_y_coords = solution_temp["mesh"]["y"]
                rho_int    = solution_temp["var"]["rho"][2:end-1, 2:end-1]
                u_int      = solution_temp["var"]["rho_u"][2:end-1, 2:end-1] ./ rho_int
                v_int      = solution_temp["var"]["rho_v"][2:end-1, 2:end-1] ./ rho_int
                internal_e = solution_temp["var"]["rho_E"][2:end-1, 2:end-1] .- 0.5 .* rho_int .* (u_int.^2 .+ v_int.^2)
                energy_factor = haskey(s["freestream"], "energy_factor") ? s["freestream"]["energy_factor"] : 1.0
                plot_c_data = internal_e .* energy_factor ./ s["freestream"]["e"]
                name_latex = L"e/e_\infty"

            elseif current_var_name_str == "P"
                plot_x_coords = solution_temp["mesh"]["x"]
                plot_y_coords = solution_temp["mesh"]["y"]
                solution_temp = UPDATE_THERMODYNAMIC_PROPERTIES(solution_temp, chemistry)
                plot_c_data = solution_temp["var"]["p"][2:end-1, 2:end-1] ./ s["freestream"]["p_0"]
                name_latex = L"P/P_\infty"

            elseif current_var_name_str == "T"
                plot_x_coords = solution_temp["mesh"]["x"]
                plot_y_coords = solution_temp["mesh"]["y"]
                plot_c_data = solution_temp["var"]["T"][2:end-1, 2:end-1] .* s["freestream"]["energy_factor"] ./ s["freestream"]["cv"]
                name_latex = L"T[K]"

            elseif current_var_name_str == "gamma*"
                if s["chemistry"]["is_chemistry_enabled"]
                    plot_x_coords = solution_temp["mesh"]["x"]
                    plot_y_coords = solution_temp["mesh"]["y"]
                    plot_c_data = solution_temp["var"]["gamma_star"][2:end-1, 2:end-1]
                    name_latex = L"\gamma^*"
                end

            elseif current_var_name_str == "cv*"
                if s["chemistry"]["is_chemistry_enabled"]
                    plot_x_coords = solution_temp["mesh"]["x"]
                    plot_y_coords = solution_temp["mesh"]["y"]
                    plot_c_data = solution_temp["var"]["cv_star"][2:end-1, 2:end-1]
                    name_latex = L"c_v^*"
                end

            elseif current_var_name_str == "CO2"
                if s["chemistry"]["is_chemistry_enabled"]
                    plot_x_coords = solution_temp["mesh"]["x"]
                    plot_y_coords = solution_temp["mesh"]["y"]
                    e_field = solution_temp["var"]["rho_E"] ./ solution_temp["var"]["rho"] .-
                        1/2 .* (solution_temp["var"]["rho_u"].^2 .+ solution_temp["var"]["rho_v"].^2) ./ solution_temp["var"]["rho"].^2
                    plot_c_data = chemistry["eval_log10_CO2"](
                        s["freestream"]["rho_factor"] .* solution_temp["var"]["rho"][2:end-1, 2:end-1],
                        s["freestream"]["energy_factor"] .* e_field[2:end-1, 2:end-1])
                    name_latex = L"log_{10}(X_{CO_2})"
                end

            elseif current_var_name_str == "CO"
                if s["chemistry"]["is_chemistry_enabled"]
                    plot_x_coords = solution_temp["mesh"]["x"]
                    plot_y_coords = solution_temp["mesh"]["y"]
                    e_field = solution_temp["var"]["rho_E"] ./ solution_temp["var"]["rho"] .-
                        1/2 .* (solution_temp["var"]["rho_u"].^2 .+ solution_temp["var"]["rho_v"].^2) ./ solution_temp["var"]["rho"].^2
                    plot_c_data = chemistry["eval_log10_CO"](
                        s["freestream"]["rho_factor"] .* solution_temp["var"]["rho"][2:end-1, 2:end-1],
                        s["freestream"]["energy_factor"] .* e_field[2:end-1, 2:end-1])
                    name_latex = L"log_{10}(X_{CO})"
                end

            elseif current_var_name_str == "C"
                if s["chemistry"]["is_chemistry_enabled"]
                    plot_x_coords = solution_temp["mesh"]["x"]
                    plot_y_coords = solution_temp["mesh"]["y"]
                    e_field = solution_temp["var"]["rho_E"] ./ solution_temp["var"]["rho"] .-
                        1/2 .* (solution_temp["var"]["rho_u"].^2 .+ solution_temp["var"]["rho_v"].^2) ./ solution_temp["var"]["rho"].^2
                    plot_c_data = chemistry["eval_log10_C"](
                        s["freestream"]["rho_factor"] .* solution_temp["var"]["rho"][2:end-1, 2:end-1],
                        s["freestream"]["energy_factor"] .* e_field[2:end-1, 2:end-1])
                    name_latex = L"log_{10}(X_{C})"
                end

            elseif current_var_name_str == "O"
                if s["chemistry"]["is_chemistry_enabled"]
                    plot_x_coords = solution_temp["mesh"]["x"]
                    plot_y_coords = solution_temp["mesh"]["y"]
                    e_field = solution_temp["var"]["rho_E"] ./ solution_temp["var"]["rho"] .-
                        1/2 .* (solution_temp["var"]["rho_u"].^2 .+ solution_temp["var"]["rho_v"].^2) ./ solution_temp["var"]["rho"].^2
                    plot_c_data = chemistry["eval_log10_O"](
                        s["freestream"]["rho_factor"] .* solution_temp["var"]["rho"][2:end-1, 2:end-1],
                        s["freestream"]["energy_factor"] .* e_field[2:end-1, 2:end-1])
                    name_latex = L"log_{10}(X_{O})"
                end

            elseif current_var_name_str == "O2"
                if s["chemistry"]["is_chemistry_enabled"]
                    plot_x_coords = solution_temp["mesh"]["x"]
                    plot_y_coords = solution_temp["mesh"]["y"]
                    e_field = solution_temp["var"]["rho_E"] ./ solution_temp["var"]["rho"] .-
                        1/2 .* (solution_temp["var"]["rho_u"].^2 .+ solution_temp["var"]["rho_v"].^2) ./ solution_temp["var"]["rho"].^2
                    plot_c_data = chemistry["eval_log10_O2"](
                        s["freestream"]["rho_factor"] .* solution_temp["var"]["rho"][2:end-1, 2:end-1],
                        s["freestream"]["energy_factor"] .* e_field[2:end-1, 2:end-1])
                    name_latex = L"log_{10}(X_{O_2})"
                end

            elseif current_var_name_str == "N2"
                if s["chemistry"]["is_chemistry_enabled"]
                    plot_x_coords = solution_temp["mesh"]["x"]
                    plot_y_coords = solution_temp["mesh"]["y"]
                    e_field = solution_temp["var"]["rho_E"] ./ solution_temp["var"]["rho"] .-
                        1/2 .* (solution_temp["var"]["rho_u"].^2 .+ solution_temp["var"]["rho_v"].^2) ./ solution_temp["var"]["rho"].^2
                    plot_c_data = chemistry["eval_log10_N2"](
                        s["freestream"]["rho_factor"] .* solution_temp["var"]["rho"][2:end-1, 2:end-1],
                        s["freestream"]["energy_factor"] .* e_field[2:end-1, 2:end-1])
                    name_latex = L"log_{10}(X_{N_2})"
                end

            elseif current_var_name_str == "N"
                if s["chemistry"]["is_chemistry_enabled"]
                    plot_x_coords = solution_temp["mesh"]["x"]
                    plot_y_coords = solution_temp["mesh"]["y"]
                    e_field = solution_temp["var"]["rho_E"] ./ solution_temp["var"]["rho"] .-
                        1/2 .* (solution_temp["var"]["rho_u"].^2 .+ solution_temp["var"]["rho_v"].^2) ./ solution_temp["var"]["rho"].^2
                    plot_c_data = chemistry["eval_log10_N"](
                        s["freestream"]["rho_factor"] .* solution_temp["var"]["rho"][2:end-1, 2:end-1],
                        s["freestream"]["energy_factor"] .* e_field[2:end-1, 2:end-1])
                    name_latex = L"log_{10}(X_{N})"
                end

            elseif current_var_name_str == "NO"
                if s["chemistry"]["is_chemistry_enabled"]
                    plot_x_coords = solution_temp["mesh"]["x"]
                    plot_y_coords = solution_temp["mesh"]["y"]
                    e_field = solution_temp["var"]["rho_E"] ./ solution_temp["var"]["rho"] .-
                        1/2 .* (solution_temp["var"]["rho_u"].^2 .+ solution_temp["var"]["rho_v"].^2) ./ solution_temp["var"]["rho"].^2
                    plot_c_data = chemistry["eval_log10_NO"](
                        s["freestream"]["rho_factor"] .* solution_temp["var"]["rho"][2:end-1, 2:end-1],
                        s["freestream"]["energy_factor"] .* e_field[2:end-1, 2:end-1])
                    name_latex = L"log_{10}(X_{NO})"
                end

            elseif current_var_name_str == "H2"
                if s["chemistry"]["is_chemistry_enabled"]
                    plot_x_coords = solution_temp["mesh"]["x"]
                    plot_y_coords = solution_temp["mesh"]["y"]
                    e_field = solution_temp["var"]["rho_E"] ./ solution_temp["var"]["rho"] .-
                        1/2 .* (solution_temp["var"]["rho_u"].^2 .+ solution_temp["var"]["rho_v"].^2) ./ solution_temp["var"]["rho"].^2
                    plot_c_data = chemistry["eval_log10_H2"](
                        s["freestream"]["rho_factor"] .* solution_temp["var"]["rho"][2:end-1, 2:end-1],
                        s["freestream"]["energy_factor"] .* e_field[2:end-1, 2:end-1])
                    name_latex = L"log_{10}(X_{H_2})"
                end

            elseif current_var_name_str == "H"
                if s["chemistry"]["is_chemistry_enabled"]
                    plot_x_coords = solution_temp["mesh"]["x"]
                    plot_y_coords = solution_temp["mesh"]["y"]
                    e_field = solution_temp["var"]["rho_E"] ./ solution_temp["var"]["rho"] .-
                        1/2 .* (solution_temp["var"]["rho_u"].^2 .+ solution_temp["var"]["rho_v"].^2) ./ solution_temp["var"]["rho"].^2
                    plot_c_data = chemistry["eval_log10_H"](
                        s["freestream"]["rho_factor"] .* solution_temp["var"]["rho"][2:end-1, 2:end-1],
                        s["freestream"]["energy_factor"] .* e_field[2:end-1, 2:end-1])
                    name_latex = L"log_{10}(X_{H})"
                end

            ## Difference fields (current minus baseline)
            elseif current_var_name_str == "Drho"
                plot_x_coords = solution_temp["mesh"]["x"]
                plot_y_coords = solution_temp["mesh"]["y"]
                rho_curr = solution_temp["var"]["rho"][2:end-1, 2:end-1] ./ solution_temp["freestream"]["rho_0"]
                rho_base = solution_temp_base["rho"][2:end-1, 2:end-1] ./ solution_temp_base["freestream"]["rho_0"]
                plot_c_data = rho_curr .- rho_base
                name_latex = L"\Delta (\rho/\rho_\infty)"

            elseif current_var_name_str == "Du"
                plot_x_coords = solution_temp["mesh"]["x"]
                plot_y_coords = solution_temp["mesh"]["y"]
                u_curr = (solution_temp["var"]["rho_u"][2:end-1, 2:end-1] ./ solution_temp["var"]["rho"][2:end-1, 2:end-1]) ./ U_infty
                u_base = (solution_temp_base["rho_u"][2:end-1, 2:end-1] ./ solution_temp_base["rho"][2:end-1, 2:end-1]) ./ U_infty
                plot_c_data = u_curr .- u_base
                name_latex = L"\Delta (u/U_\infty)"

            elseif current_var_name_str == "Dv"
                plot_x_coords = solution_temp["mesh"]["x"]
                plot_y_coords = solution_temp["mesh"]["y"]
                v_curr = (solution_temp["var"]["rho_v"][2:end-1, 2:end-1] ./ solution_temp["var"]["rho"][2:end-1, 2:end-1]) ./ U_infty
                v_base = (solution_temp_base["rho_v"][2:end-1, 2:end-1] ./ solution_temp_base["rho"][2:end-1, 2:end-1]) ./ U_infty
                plot_c_data = v_curr .- v_base
                name_latex = L"\Delta (v/U_\infty)"

            elseif current_var_name_str == "De"
                plot_x_coords = solution_temp["mesh"]["x"]
                plot_y_coords = solution_temp["mesh"]["y"]
                rho_int_curr    = solution_temp["var"]["rho"][2:end-1, 2:end-1]
                u_int_curr      = solution_temp["var"]["rho_u"][2:end-1, 2:end-1] ./ rho_int_curr
                v_int_curr      = solution_temp["var"]["rho_v"][2:end-1, 2:end-1] ./ rho_int_curr
                e_curr_internal = solution_temp["var"]["rho_E"][2:end-1, 2:end-1] .-
                    0.5 .* rho_int_curr .* (u_int_curr.^2 .+ v_int_curr.^2)
                e_curr_norm     = e_curr_internal ./ (solution_temp["freestream"]["rho_0"] * e_infty)

                rho_int_base    = solution_temp_base["rho"][2:end-1, 2:end-1]
                u_int_base      = solution_temp_base["rho_u"][2:end-1, 2:end-1] ./ rho_int_base
                v_int_base      = solution_temp_base["rho_v"][2:end-1, 2:end-1] ./ rho_int_base
                e_base_internal = solution_temp_base["rho_E"][2:end-1, 2:end-1] .-
                    0.5 .* rho_int_base .* (u_int_base.^2 .+ v_int_base.^2)
                e_base_norm     = e_base_internal ./ (solution_temp_base["freestream"]["rho_0"] * e_infty)

                plot_c_data = e_curr_norm .- e_base_norm
                name_latex = L"\Delta (e/e_\infty)"

            elseif current_var_name_str == "DP"
                plot_x_coords = solution_temp["mesh"]["x"]
                plot_y_coords = solution_temp["mesh"]["y"]
                solution_temp      = UPDATE_THERMODYNAMIC_PROPERTIES(solution_temp, chemistry)
                solution_temp_base = UPDATE_THERMODYNAMIC_PROPERTIES(solution_temp_base, chemistry)
                p_curr_norm = solution_temp["var"]["p"][2:end-1, 2:end-1] ./ p_infty
                p_base_norm = solution_temp_base["p"][2:end-1, 2:end-1] ./ p_infty
                plot_c_data = p_curr_norm .- p_base_norm
                name_latex = L"\Delta (P/P_\infty)"

            else
                @warn "Unknown variable selected for plotting: $current_var_name_str"
                plot_x_coords = Matrix{Float64}(undef, 0, 0)
            end

            ## Skip if data could not be generated
            if length(plot_x_coords) == 0 || length(plot_y_coords) == 0 || length(plot_c_data) == 0
                @warn "Data for subplot $current_var_name_str could not be fully generated or is empty. Skipping this subplot."
                continue
            end

            ## Dimension check before plotting
            if size(plot_x_coords) != size(plot_y_coords) || size(plot_x_coords) != size(plot_c_data)
                @warn "Dimension mismatch for $current_var_name_str: X=$(size(plot_x_coords)), Y=$(size(plot_y_coords)), C=$(size(plot_c_data)). Skipping subplot."
                continue
            end

            ## Create subplot
            # Compute colour limits
            min_val = minimum(plot_c_data)
            max_val = maximum(plot_c_data)
            if min_val == max_val
                min_val = min_val * 2
                max_val = max_val * 2
                if min_val == 0.0 && max_val == 0.0
                    min_val = -0.01
                    max_val = 0.01
                end
            end

            clims_val = (min_val, max_val)
            if isnan(min_val) || isinf(min_val) || isnan(max_val) || isinf(max_val)
                clims_val = :auto
            end

            # Store caxis limits for reuse
            if first_plot_call
                handles["caxis_limits_per_var"][v_idx] = clims_val
            else
                if haskey(handles["caxis_limits_per_var"], v_idx)
                    clims_val = handles["caxis_limits_per_var"][v_idx]
                end
            end

            # Use PyPlot pcolormesh with full 2-D coordinate arrays
            # for correct rendering on curvilinear grids
            (fig_sp, ax_sp) = curvilinear_heatmap(plot_x_coords, plot_y_coords, plot_c_data,
                                clims=clims_val,
                                colorbar_label=string(name_latex),
                                title_str=string(name_latex),
                                xlabel_str="x/L",
                                ylabel_str="y/L")

            # Add shock line overlay
            if shock_line_exists
                ax_sp.plot(x_p, y_p, color="black", linewidth=2)
            end

            push!(subplot_plots, (fig_sp, ax_sp))
        end

        ## Display all subplot figures
        if !isempty(subplot_plots)
            for (fig_sp, ax_sp) in subplot_plots
                fig_sp.suptitle(@sprintf("Iteration: %d, Time: %.3f",
                                        solution_temp["time_integration"]["iter"],
                                        solution_temp["time_integration"]["t"]))
                fig_sp.tight_layout()
                save_or_display_figure(fig_sp, name=@sprintf("running_iter_%06d", solution_temp["time_integration"]["iter"]))
            end
            handles["combined_plot"] = subplot_plots
        end

    end

    return handles
end
