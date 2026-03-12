using Plots, Printf, LinearAlgebra, LoopVectorization

# PLOT_MODES_CHU_NORM  Visualise energy-budget fields.
#
#   PLOT_MODES_CHU_NORM(freestream_disturbances, A, s, chemistry, V,
#        T_plot, w_infty)
#
#   Generates publication-quality 2-D pseudocolor plots of perturbation
#   quantities (pressure, density, velocity components, vorticity,
#   divergence, entropy, energy flux) and, when applicable, kinetic and
#   entropic energy-budget maps.  Optionally evolves the perturbation state
#   forward in time on the GPU using a Taylor-series approximation to the
#   matrix exponential.
#
#   Inputs:
#       freestream_disturbances - (Bool) true when upstream disturbance
#                                 modes are included in the state vector V.
#       A                       - (sparse matrix) linearised operator used
#                                 for GPU time integration.
#       s                       - (Dict) base-flow solution containing
#                                 grid coordinates, thermodynamic fields,
#                                 shock geometry, and solver parameters.
#       chemistry               - (Dict) chemistry model data (needed for
#                                 entropy budgets and remeshing).
#       V                       - (column vector) perturbation state vector
#                                 [rho'; rho_u'; rho_v'; rho_E'; ...]; may
#                                 include shock displacement and freestream
#                                 mode amplitudes.
#       T_plot                  - (scalar) physical time at which to
#                                 evaluate. T_plot = 0 plots the initial
#                                 perturbation; T_plot > 0 triggers GPU
#                                 time integration.
#       w_infty                 - (vector) angular frequencies of the
#                                 freestream disturbance modes (empty when
#                                 not used).
#
#   Outputs:
#       (none) - all output is graphical (Plots.jl figure windows).
#
#   Notes:
#       - Auxiliary ghost / shock cells are stripped via
#         REMOVE_AUX_SHOCK_CELLS before any processing.
#       - When T_plot == 0 and freestream_disturbances is true the solution
#         is remeshed with a wider wall-normal extent so the freestream
#         region is visible.
#       - GPU time integration uses LINEAR_INTEGRATION_GPU (local
#         subfunction) with a user-specified Taylor-series order.
#       - Scaling factors are derived from either the energy-flux inflow
#         rate or the maximum total energy rate of change at t = 0.
#       - All contour plots are created by CREATE_PUBLICATION_PLOT.
#       - Figure export is available through EXPORT_ALL_FIGURES.
#
# Part of: Hypersonics Stability MATLAB Solver - Postprocessing Module

function PLOT_MODES_CHU_NORM(freestream_disturbances::Bool, A,
                             s::Dict{String, Any}, chemistry::Dict{String, Any},
                             V::AbstractVector, T_plot::Real,
                             w_infty::AbstractVector)

    SET_PLOT_DEFAULTS()

    if s["shock"]["enabled"]
        s = REMOVE_AUX_SHOCK_CELLS(s)
    end

    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]
    N_v = Nchi * Neta
    if freestream_disturbances
        N_f = length(w_infty)
    else
        N_f = 0
    end

    ## Extract perturbation fields from eigenvector
    pert_rho   = real.(reshape(V[(0)*N_v+1:1*N_v], Nchi, Neta))
    pert_rho_u = real.(reshape(V[(1)*N_v+1:2*N_v], Nchi, Neta))
    pert_rho_v = real.(reshape(V[(2)*N_v+1:3*N_v], Nchi, Neta))
    pert_rho_E = real.(reshape(V[(3)*N_v+1:4*N_v], Nchi, Neta))

    ## Extract freestream disturbance amplitudes (if present)
    size_freestream_problem_1 = 4 * Nchi * Neta + 4 * Nchi * N_f + Nchi
    size_freestream_problem_2 = 4 * Nchi * Neta + 4 * Nchi * N_f

    pert_infty_all = Vector{eltype(V)}()
    if length(V) == size_freestream_problem_1
        idx_start_infty = 4 * Nchi * Neta + Nchi
        pert_infty_all = V[idx_start_infty+1:size_freestream_problem_1]
    elseif length(V) == size_freestream_problem_2
        idx_start_infty = 4 * Nchi * Neta
        pert_infty_all = V[idx_start_infty+1:size_freestream_problem_2]
    end

    ## Compute non-dimensional scaling factors
    scale_kinetic = 0.0
    scale_entropic = 0.0

    if s["shock"]["enabled"] && freestream_disturbances && T_plot > 0
        solution_remesh = false
        (pert_rho, pert_rho_u, pert_rho_v, pert_rho_E) = ADD_FREESTREAM_DISTURBANCES(
            pert_rho, pert_rho_u, pert_rho_v, pert_rho_E,
            pert_infty_all, s, w_infty, T_plot, solution_remesh)
        (_, dE_dt_V_inflow) = ENERGY_INFLOW(V, s, w_infty)
        scale_kinetic  = dE_dt_V_inflow
        scale_entropic = dE_dt_V_inflow

    elseif s["shock"]["enabled"] && freestream_disturbances && T_plot == 0
        # Remesh to visualise freestream region
        s["shock"]["remesh_shock_distance"] = 8
        s["mesh"]["Neta"] = s["mesh"]["Neta"] * 5
        solution_old = deepcopy(s)
        disturbances = true
        s = RESTART_SOLUTION(s, solution_old, chemistry, disturbances)
        s = REMOVE_AUX_SHOCK_CELLS(s)
        solution_remesh = true
        (pert_rho, pert_rho_u, pert_rho_v, pert_rho_E) = ADD_FREESTREAM_DISTURBANCES(
            pert_rho, pert_rho_u, pert_rho_v, pert_rho_E,
            pert_infty_all, s, w_infty, 0, solution_remesh)
    else
        # No freestream disturbances -- scale by max total energy rate
        output_flow = true
        budgets_kinetic = COMPUTE_BUDGETS_KINETIC(V, s, output_flow)
        temp_kinetic = budgets_kinetic["A_adv"] .+ budgets_kinetic["P_mom"] .+
            budgets_kinetic["P_mass"] .+ budgets_kinetic["Dilat_P"] .+
            budgets_kinetic["Transport_tau"] .+ budgets_kinetic["Transport_p"] .+
            budgets_kinetic["Dissipation"]
        budgets_entropy = COMPUTE_BUDGETS_ENTROPY(V, s, chemistry, output_flow)
        temp_entropic = budgets_entropy["A_adv"] .+ budgets_entropy["P_mom"] .+
            budgets_entropy["P_mass"] .+ budgets_entropy["P_T"] .+
            budgets_entropy["Transport"] .+ budgets_entropy["Dissipation"] .+
            budgets_entropy["Source"]
        temp = temp_entropic .+ temp_kinetic
        scale_kinetic  = maximum(abs.(temp))
        scale_entropic = maximum(abs.(temp))
    end

    ## Time-advance perturbation on GPU (optional)
    if T_plot > 0
        order = 5

        s = CFL_TIMESTEP(s)
        t_temp = s["time_integration"]["dt"] / s["time_integration"]["CFL"]
        n_t = round(Int, T_plot / t_temp)
        dt  = T_plot / n_t
        @printf("\n")
        @printf("Number timesteps = %f ", n_t)
        @printf("\n")

        gv     = V
        gA     = A
        gorder = order
        gn_t   = n_t
        gdt    = dt

        # Note: GPU device selection in Julia uses CUDA.jl
        # gpu = CUDA.device()
        # println(CUDA.name(gpu) * " GPU selected.")

        t_start = time()
        V = LINEAR_INTEGRATION_GPU(gv, gn_t, gdt, gorder, gA)
        t_elapsed = time() - t_start
        @printf("Elapsed time: %.6f seconds\n", t_elapsed)

        pert_rho   = real.(reshape(V[(0)*N_v+1:1*N_v], Nchi, Neta))
        pert_rho_u = real.(reshape(V[(1)*N_v+1:2*N_v], Nchi, Neta))
        pert_rho_v = real.(reshape(V[(2)*N_v+1:3*N_v], Nchi, Neta))
        pert_rho_E = real.(reshape(V[(3)*N_v+1:4*N_v], Nchi, Neta))

        solution_remesh = false
        if freestream_disturbances
            (pert_rho, pert_rho_u, pert_rho_v, pert_rho_E) = ADD_FREESTREAM_DISTURBANCES(
                pert_rho, pert_rho_u, pert_rho_v, pert_rho_E,
                pert_infty_all, s, w_infty, T_plot, solution_remesh)
        end
    end

    ## Derived perturbation quantities at time T
    p_0   = s["var"]["p"][2:end-1, 2:end-1]
    rho_0 = s["var"]["rho"][2:end-1, 2:end-1]
    u_0   = s["var"]["rho_u"][2:end-1, 2:end-1] ./ s["var"]["rho"][2:end-1, 2:end-1]
    v_0   = s["var"]["rho_v"][2:end-1, 2:end-1] ./ s["var"]["rho"][2:end-1, 2:end-1]
    pert_u = (pert_rho_u .- u_0 .* pert_rho) ./ rho_0
    pert_v = (pert_rho_v .- v_0 .* pert_rho) ./ rho_0

    cell_x_wall_normal = (s["mesh"]["x_wall_normal"][1:end-1, 1] .+ s["mesh"]["x_wall_normal"][2:end, 1]) ./ 2
    cell_y_wall_normal = (s["mesh"]["y_wall_normal"][1:end-1, 1] .+ s["mesh"]["y_wall_normal"][2:end, 1]) ./ 2
    cell_x_wall_tan    =  cell_y_wall_normal
    cell_y_wall_tan    = -cell_x_wall_normal

    pert_u_tan  = pert_u .* cell_x_wall_tan .+ pert_v .* cell_y_wall_tan
    pert_u_norm = pert_u .* cell_x_wall_normal .+ pert_v .* cell_y_wall_normal
    pert_u_mag  = (pert_u .* u_0 .+ pert_v .* v_0) ./ sqrt.(u_0.^2 .+ v_0.^2)

    t1 = -u_0 .* pert_rho_u
    t2 = -v_0 .* pert_rho_v
    t3 = (u_0.^2 .+ v_0.^2) .* pert_rho ./ 2
    pert_p = (s["var"]["gamma_star"][2:end-1, 2:end-1] .- 1) .* (pert_rho_E .+ t1 .+ t2 .+ t3)

    pert_u_Ext = EXTEND_TO_GHOST_POINTS(pert_u, s)
    pert_v_Ext = EXTEND_TO_GHOST_POINTS(pert_v, s)
    (dpert_u_dx, dpert_u_dy) = DERIVATIVE_EXT(pert_u_Ext, s)
    (dpert_v_dx, dpert_v_dy) = DERIVATIVE_EXT(pert_v_Ext, s)
    pert_div  = dpert_u_dx .+ dpert_v_dy
    pert_vort = dpert_v_dx .- dpert_u_dy

    pert_entropy = pert_p ./ p_0 ./ (s["var"]["gamma_star"][2:end-1, 2:end-1] .- 1) .-
        pert_rho ./ rho_0 .* s["var"]["gamma_star"][2:end-1, 2:end-1] ./
        (s["var"]["gamma_star"][2:end-1, 2:end-1] .- 1)

    ## Shock spline (possibly perturbed)
    N_plot_shock_poly = 1000
    x_p = zeros(0)
    y_p = zeros(0)
    if s["shock"]["enabled"] == true
        solution2 = deepcopy(s)
        if s["stability_analysis"]["perturb_shock"]
            dr_shock = real.(V[(4)*N_v+1:4*N_v+Nchi])
            ang = s["shock"]["beta"]
            solution2["shock"]["points_x"] = s["shock"]["points_x"] .-
                dr_shock .* sin.(ang) ./ s["stability_analysis"]["perturbation_magnitude"] ./ 2000000
            solution2["shock"]["points_y"] = s["shock"]["points_y"] .+
                dr_shock .* cos.(ang) ./ s["stability_analysis"]["perturbation_magnitude"] ./ 2000000
            fig_shock, ax_shock = PyPlot.subplots()
            ax_shock.plot(solution2["shock"]["points_x"], solution2["shock"]["points_y"], label="original")
            solution2 = LEAST_SQUARES_SHOCK_POINTS(solution2)
            ax_shock.plot(solution2["shock"]["points_x"], solution2["shock"]["points_y"], label="filtered")
            ax_shock.legend()
            fig_shock.tight_layout()
            save_or_display_figure(fig_shock, name="shock_spline_chu")
        elseif T_plot > 0
            valid_ix      = collect(1:Nchi)
            sc_idy        = s["shock"]["cell_indices"][valid_ix, 1]
            lin_idx_shock = [CartesianIndex(valid_ix[k], sc_idy[k] - 1) for k in 1:Nchi]
            p_shock       = real.([pert_p[idx] for idx in lin_idx_shock])
            max_p_shock   = maximum(p_shock)
            min_p_shock   = minimum(p_shock)
            p_shock       = (p_shock .- min_p_shock) ./ (max_p_shock - min_p_shock) .- 1/2
            displace      = -sin.(pi .* p_shock)
            solution2["shock"]["points_x"] = s["shock"]["points_x"] .+ displace ./ 500
            solution2["shock"]["points_y"] = s["shock"]["points_y"] .+ displace ./ 500
            solution2 = LEAST_SQUARES_SHOCK_POINTS(solution2)
        end

        chi_p = collect(range(solution2["shock"]["points_chi"][1, 1],
                              solution2["shock"]["points_chi"][end, 1],
                              length=N_plot_shock_poly))
        x_p = ppval(solution2["shock"]["spline_func_x"], chi_p)
        y_p = ppval(solution2["shock"]["spline_func_y"], chi_p)
    end

    ## Configure plot parameters (font sizes and figure dimensions are
    ## handled automatically by SET_PLOT_DEFAULTS / auto_figure_size)
    plot_params = Dict{String, Any}()
    plot_params["freestream_disturbance"] = (T_plot == 0 && freestream_disturbances)
    plot_params["line_width"]      = 1.5
    plot_params["colormap_name"]   = :turbo
    plot_params["use_common_clim"] = false

    ## Define energy-budget plot lists
    if T_plot > 0 || !freestream_disturbances
        output_flow = true
        budgets_kinetic = COMPUTE_BUDGETS_KINETIC(V, s, output_flow)

        kinetic_plots = [
            (budgets_kinetic["A_adv"] ./ scale_kinetic,          L"\mathcal{A}^k"),
            (budgets_kinetic["P_mom"] ./ scale_kinetic,          L"\mathcal{P}_u^k"),
            (budgets_kinetic["P_mass"] ./ scale_kinetic,         L"\mathcal{P}_\rho^k"),
            (budgets_kinetic["Dilat_P"] ./ scale_kinetic,        L"\Pi_d^k"),
            (budgets_kinetic["Transport_p"] ./ scale_kinetic,    L"\mathcal{T}^k_p"),
            (budgets_kinetic["Transport_tau"] ./ scale_kinetic,  L"\mathcal{T}^k_\tau"),
            (budgets_kinetic["Dissipation"] ./ scale_kinetic,    L"\mathcal{D}^k")
        ]

        output_flow = true
        budgets_entropy = COMPUTE_BUDGETS_ENTROPY(V, s, chemistry, output_flow)

        entropy_plots = [
            (budgets_entropy["A_adv"] ./ scale_entropic,         L"\mathcal{A}^s"),
            (budgets_entropy["P_mom"] ./ scale_entropic,         L"\mathcal{P}_u^s"),
            (budgets_entropy["P_mass"] ./ scale_entropic,        L"\mathcal{P}_\rho^s"),
            (budgets_entropy["P_T"] ./ scale_entropic,           L"\mathcal{P}_T^s"),
            (budgets_entropy["Transport"] ./ scale_entropic,     L"\mathcal{T}^s"),
            (budgets_entropy["Dissipation"] ./ scale_entropic,   L"\mathcal{D}^s"),
            (budgets_entropy["Source"] ./ scale_entropic,        L"\mathcal{S}^s")
        ]
    end

    ## Build coordinate arrays and create all plots
    @printf("Creating perturbation plots...\n")

    y_coord = s["mesh"]["y"]
    x_coord = s["mesh"]["x"]
    y_shock = y_p
    x_shock = x_p

    ## Render kinetic and entropic budget plots
    if T_plot > 0 || !freestream_disturbances
        @printf("Creating kinetic budget plots...\n")
        for i in 1:length(kinetic_plots)
            data      = kinetic_plots[i][1]
            cb_label  = kinetic_plots[i][2]
            data_plot = real.(data)
            CREATE_PUBLICATION_PLOT_CHU(y_coord, x_coord, data_plot, cb_label,
                s["shock"]["enabled"], y_shock, x_shock, plot_params, s)
        end

        @printf("Creating entropy budget plots...\n")
        for i in 1:length(entropy_plots)
            data      = entropy_plots[i][1]
            cb_label  = entropy_plots[i][2]
            data_plot = real.(data)
            CREATE_PUBLICATION_PLOT_CHU(y_coord, x_coord, data_plot, cb_label,
                s["shock"]["enabled"], y_shock, x_shock, plot_params, s)
        end
    end

    @printf("All plots created successfully!\n")

    return nothing
end


# Note: LINEAR_INTEGRATION_GPU is defined in PLOT_MODES.jl


# ========================================================================
#  LOCAL FUNCTION: CREATE_PUBLICATION_PLOT_CHU
# ========================================================================
# CREATE_PUBLICATION_PLOT_CHU  Generate a single publication-quality pcolor
#   figure with a symmetric colorbar centred at zero.
#   (Note: arguments are (y, x, data) matching the MATLAB call order in
#    PLOT_MODES_CHU_NORM where x/y are swapped relative to PLOT_MODES.)
#
#   Inputs:
#       x, y      - (2-D arrays) coordinate grids for heatmap.
#       data      - (2-D array) field data to plot.
#       cb_label  - (String/LaTeXString) label for the colorbar.
#       has_shock - (Bool) overlay shock line when true.
#       x_shock, y_shock - (vectors) shock-line coordinates.
#       params    - (Dict) figure, font, and colormap parameters.
#       s         - (Dict) solution struct.
#
#   Outputs:
#       (none) - displays a new figure window.

function CREATE_PUBLICATION_PLOT_CHU(x::AbstractMatrix, y::AbstractMatrix,
                                     data::AbstractMatrix, cb_label,
                                     has_shock::Bool, x_shock::AbstractVector,
                                     y_shock::AbstractVector,
                                     params::Dict{String, Any},
                                     s::Dict{String, Any})

    ## Symmetric colour axis truncated to 1 significant figure
    max_abs_val = maximum(abs.(data))
    function truncate_to_1sig(x_val)
        if x_val == 0.0
            return 0.0
        end
        return floor(x_val / 10.0^floor(log10(abs(x_val)))) * 10.0^floor(log10(abs(x_val)))
    end
    max_abs_val = truncate_to_1sig(max_abs_val)

    scaling_range = 1 / 2
    clim_lo = -max_abs_val * scaling_range
    clim_hi =  max_abs_val * scaling_range

    ## Automatic figure size from data aspect ratio
    fig_size = auto_figure_size(x, y)

    ## Pseudocolor plot using PyPlot pcolormesh for curvilinear grids
    # Note: caller passes (y_coord, x_coord) so x=mesh.y and y=mesh.x here.
    # For pcolormesh we swap back: y (=mesh.x) on horizontal, x (=mesh.y) on vertical
    (fig, ax) = curvilinear_heatmap(y, x, data,
                       cmap=params["colormap_name"],
                       clims=(clim_lo, clim_hi),
                       colorbar_label=string(cb_label),
                       xlabel_str="x",
                       ylabel_str="y",
                       fig_size=fig_size)

    if has_shock
        ax.plot(x_shock, y_shock, color="black", linewidth=params["line_width"])
    end

    fig.tight_layout()
    save_or_display_figure(fig, name=string("mode_chu_", replace(string(cb_label), r"[^\w]" => "_")))

    return nothing
end

# Note: EXPORT_ALL_FIGURES, SCALINGS, and ADD_FREESTREAM_DISTURBANCES
# are defined in PLOT_MODES.jl and shared across postprocessing files.
