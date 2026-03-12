using Plots, Printf, LinearAlgebra, LoopVectorization

# PLOT_MODES  Visualise perturbation eigenmodes and energy-budget fields.
#
#   PLOT_MODES(freestream_disturbances, A, s, chemistry, V,
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

function PLOT_MODES(freestream_disturbances::Bool, A, s::Dict{String, Any},
                    chemistry::Dict{String, Any}, V::AbstractVector,
                    T_plot::Real, w_infty::AbstractVector=Float64[])

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
        (scale_u, scale_vort, scale_div, scale_rho, scale_p, scale_s) = SCALINGS(
            s, pert_rho, pert_rho_u, pert_rho_v, pert_rho_E)
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
        (scale_u, scale_vort, scale_div, scale_rho, scale_p, scale_s) = SCALINGS(
            s, pert_rho, pert_rho_u, pert_rho_v, pert_rho_E)

    else
        # No freestream disturbances -- scale by max total energy rate
        (scale_u, scale_vort, scale_div, scale_rho, scale_p, scale_s) = SCALINGS(
            s, pert_rho, pert_rho_u, pert_rho_v, pert_rho_E)
    end

    ## Time-advance perturbation on GPU (optional)
    if T_plot > 0
        order = 5

        s = CFL_TIMESTEP(s)
        n_t = round(Int, T_plot / s["time_integration"]["dt"])
        dt = T_plot / n_t
        @printf("\n")
        @printf("Number timesteps = %f ", n_t)
        @printf("\n")

        ## Detect GPU availability
        use_gpu = CUDA.functional()
        n_vec = size(A, 1)

        if use_gpu
            gpu = CUDA.device()
            println("$(CUDA.name(gpu)) GPU selected.")
            gA       = CuSparseMatrixCSC(A * dt)
            v_buf    = CUDA.zeros(ComplexF64, n_vec) ## Can be complex in freestream receptivity case
            vout_buf = CUDA.zeros(ComplexF64, n_vec)
            term_buf = CUDA.zeros(ComplexF64, n_vec)
        else
            println("No CUDA device found. Running on CPU.")
            gA       = A * dt
            v_buf    = zeros(ComplexF64, n_vec)
            vout_buf = zeros(ComplexF64, n_vec)
            term_buf = zeros(ComplexF64, n_vec)
        end

        t_start = time()
        V = LINEAR_INTEGRATION_GPU(V, n_t, order, gA, use_gpu, v_buf, vout_buf, term_buf)
        t_elapsed = time() - t_start
        @printf("Elapsed time: %.6f seconds\n", t_elapsed)

        ## GPU synchronization and memory cleanup
        if use_gpu
            CUDA.synchronize()
            CUDA.unsafe_free!(gA)
            CUDA.unsafe_free!(v_buf)
            CUDA.unsafe_free!(vout_buf)
            CUDA.unsafe_free!(term_buf)
            CUDA.synchronize()
            CUDA.reclaim()
        end

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
            save_or_display_figure(fig_shock, name="shock_spline")
        elseif T_plot > 0
            valid_ix      = collect(1:Nchi)
            sc_idy        = s["shock"]["cell_indices"][valid_ix, 1]
            lin_idx_shock = [CartesianIndex(valid_ix[k], sc_idy[k] - 1) for k in 1:Nchi]
            p_shock       = real.([pert_p[idx] for idx in lin_idx_shock])
            max_p_shock   = maximum(p_shock)
            min_p_shock   = minimum(p_shock)
            if max_p_shock ≈ min_p_shock
                displace = zeros(length(p_shock))
            else
                p_shock   = (p_shock .- min_p_shock) ./ (max_p_shock - min_p_shock) .- 1/2
                displace  = -sin.(pi .* p_shock)
            end
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

    ## Define perturbation plot list
    pert_plots = [
        (pert_p ./ scale_p,           L"\frac{P'(T)}{max(|P'(0)|)}"),
        (pert_rho ./ scale_rho,       L"\frac{\rho'(T)}{max(|\rho'(0)|)}"),
        (pert_u_mag ./ scale_u,       L"\frac{\|\vec{u}\|_2'(T)}{max(|\|\vec{u}\|_2'(0)|)}"),
        (pert_vort ./ scale_vort,     L"\frac{\omega_z'(T)}{max(|\omega_z'(0)|)}"),
        (pert_entropy ./ scale_s,     L"\frac{s'(T)}{max(|s'(0)|)}")
    ]

    ## Build coordinate arrays and create all plots
    @printf("Creating perturbation plots...\n")

    y_coord = s["mesh"]["y"]
    x_coord = s["mesh"]["x"]
    if s["shock"]["enabled"]
        y_shock = y_p
        x_shock = x_p
    else
        y_shock = [0.0]
        x_shock = [0.0]
    end

    ## Render perturbation field plots
    for i in 1:length(pert_plots)
        data      = pert_plots[i][1]
        cb_label  = pert_plots[i][2]
        data_plot = real.(data)
        CREATE_PUBLICATION_PLOT(x_coord, y_coord, data_plot, cb_label,
            s["shock"]["enabled"], x_shock, y_shock, plot_params, s)
    end

    @printf("All plots created successfully!\n")

    return nothing
end


# ========================================================================
#  LOCAL FUNCTION: LINEAR_INTEGRATION_GPU
# ========================================================================
# LINEAR_INTEGRATION_GPU  Approximate exp(A*t)*x via truncated Taylor series.
#   Works on both GPU and CPU using pre-allocated buffers.
#
#   Inputs:
#       x        - (Vector) initial state vector (CPU).
#       n_t      - (Int) number of time steps.
#       order    - (Int) Taylor-series truncation order.
#       Adt      - system matrix A*dt (CuSparse or CPU sparse).
#       use_gpu  - (Bool) true for GPU path, false for CPU path.
#       v_buf    - Pre-allocated buffer for working vector.
#       vout_buf - Pre-allocated buffer for output accumulation.
#       term_buf - Pre-allocated buffer for Taylor term.
#
#   Outputs:
#       voutCPU - (Vector) final state vector (CPU).

function LINEAR_INTEGRATION_GPU(x::AbstractVector, n_t::Int, order::Int,
                                Adt, use_gpu::Bool, v_buf, vout_buf, term_buf)
    copyto!(v_buf, x)
    copyto!(vout_buf, v_buf)

    for i in 1:n_t
        factorial_val = 1
        for j in 1:order
            mul!(term_buf, Adt, v_buf, 1.0, 0.0)
            copyto!(v_buf, term_buf)
            factorial_val = j * factorial_val
            vout_buf .+= v_buf .* (1.0 / factorial_val)
        end
        copyto!(v_buf, vout_buf)
    end

    if use_gpu
        return Array(vout_buf)
    else
        return copy(vout_buf)
    end
end


# ========================================================================
#  LOCAL FUNCTION: CREATE_PUBLICATION_PLOT
# ========================================================================
# CREATE_PUBLICATION_PLOT  Generate a single publication-quality pcolor
#   figure with a symmetric colorbar centred at zero.
#
#   Inputs:
#       x, y      - (2-D arrays) coordinate grids for heatmap.
#       data      - (2-D array) field data to plot.
#       cb_label  - (String/LaTeXString) label for the colorbar.
#       has_shock - (Bool) overlay shock line when true.
#       x_shock, y_shock - (vectors) shock-line coordinates.
#       params    - (Dict) figure, font, and colormap parameters.
#       s         - (Dict) solution struct (boundary_type, L, etc.).
#
#   Outputs:
#       (none) - displays a new figure window.

function CREATE_PUBLICATION_PLOT(x::AbstractMatrix, y::AbstractMatrix,
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

    if haskey(s, "running_plot") && haskey(s["running_plot"], "scaling_range")
        scaling_range = s["running_plot"]["scaling_range"]
    else
        scaling_range = 1 / 2
    end

    clim_lo = -max_abs_val * scaling_range
    clim_hi =  max_abs_val * scaling_range

    ## Automatic figure size from data aspect ratio
    fig_size = auto_figure_size(x, y)

    ## Pseudocolor plot using PyPlot pcolormesh for curvilinear grids
    (fig, ax) = curvilinear_heatmap(x, y, data,
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
    save_or_display_figure(fig, name=string("mode_", replace(string(cb_label), r"[^\w]" => "_")))

    return nothing
end


# ========================================================================
#  LOCAL FUNCTION: EXPORT_ALL_FIGURES
# ========================================================================
# EXPORT_ALL_FIGURES  Export a list of PyPlot figure objects to files.
#
#   EXPORT_ALL_FIGURES(figs, format, dpi)
#
#   Inputs:
#       figs   - Vector of PyPlot figure objects to export.
#       format - (String) "png", "pdf", "eps", or "tiff".
#       dpi    - (Int, optional) resolution for raster formats (default 600).
#
#   Outputs:
#       (none) - writes figure_001, figure_002, ... files to disk.

function EXPORT_ALL_FIGURES(figs::Vector, format::String, dpi::Int=600)
    fmt = lowercase(format)
    if !(fmt in ("png", "pdf", "eps", "tiff"))
        @warn "Unknown format: $format"
        return nothing
    end

    for i in 1:length(figs)
        filename = @sprintf("figure_%03d.%s", i, fmt)
        figs[i].savefig(filename, dpi=dpi, bbox_inches="tight")
        @printf("Exported: %s\n", filename)
    end
    return nothing
end


# ========================================================================
#  LOCAL FUNCTION: SCALINGS
# ========================================================================
# SCALINGS  Compute maximum-absolute-value scaling factors for each
#   perturbation quantity, used to normalise contour plots.
#
#   Inputs:
#       s               - (Dict) base-flow solution.
#       pert_rho, pert_rho_u,
#       pert_rho_v, pert_rho_E - (2-D arrays) perturbation conservative
#                                 variables.
#
#   Outputs:
#       scale_u    - max |u'_mag|
#       scale_vort - max |vorticity'|
#       scale_div  - max |divergence'|
#       scale_rho  - max |rho'|
#       scale_p    - max |p'|
#       scale_s    - max |entropy'|

function SCALINGS(s::Dict{String, Any}, pert_rho::AbstractMatrix,
                  pert_rho_u::AbstractMatrix, pert_rho_v::AbstractMatrix,
                  pert_rho_E::AbstractMatrix)

    p_0   = s["var"]["p"][2:end-1, 2:end-1]
    rho_0 = s["var"]["rho"][2:end-1, 2:end-1]
    u_0   = s["var"]["rho_u"][2:end-1, 2:end-1] ./ s["var"]["rho"][2:end-1, 2:end-1]
    v_0   = s["var"]["rho_v"][2:end-1, 2:end-1] ./ s["var"]["rho"][2:end-1, 2:end-1]

    pert_u     = (pert_rho_u .- u_0 .* pert_rho) ./ rho_0
    pert_v     = (pert_rho_v .- v_0 .* pert_rho) ./ rho_0
    pert_u_mag = (pert_u .* u_0 .+ pert_v .* v_0) ./ sqrt.(u_0.^2 .+ v_0.^2)

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

    scale_u    = maximum(abs.(pert_u_mag))
    scale_vort = maximum(abs.(pert_vort))
    scale_div  = maximum(abs.(pert_div))
    scale_rho  = maximum(abs.(pert_rho))
    scale_p    = maximum(abs.(pert_p))
    scale_s    = maximum(abs.(pert_entropy))

    return (scale_u, scale_vort, scale_div, scale_rho, scale_p, scale_s)
end


# ========================================================================
#  LOCAL FUNCTION: ADD_FREESTREAM_DISTURBANCES
# ========================================================================
# ADD_FREESTREAM_DISTURBANCES  Superpose freestream Fourier modes onto the
#   interior perturbation fields beyond the shock.
#
#   Inputs:
#       pert_rho, pert_rho_u,
#       pert_rho_v, pert_rho_E - (Nx x Ny) interior perturbation arrays.
#       pert_infty_all         - (vector) packed freestream amplitudes.
#       s                      - (Dict) grid and shock data.
#       w_infty                - (vector) freestream frequencies.
#       T_plot                 - (scalar) evaluation time.
#       solution_remesh        - (Bool) if true, zero interior first.
#
#   Outputs:
#       pert_rho, pert_rho_u,
#       pert_rho_v, pert_rho_E - Updated perturbation arrays with
#                                 freestream contributions added.

function ADD_FREESTREAM_DISTURBANCES(pert_rho::AbstractMatrix, pert_rho_u::AbstractMatrix,
                                     pert_rho_v::AbstractMatrix, pert_rho_E::AbstractMatrix,
                                     pert_infty_all::AbstractVector,
                                     s::Dict{String, Any}, w_infty::AbstractVector,
                                     T_plot::Real, solution_remesh::Bool)

    ## Unpack freestream disturbance coefficients
    N_f        = length(w_infty)
    Nx         = s["mesh"]["Nchi"]
    Ny         = s["mesh"]["Neta"]
    size_infty = Nx * N_f

    rho_infty_all   = zeros(ComplexF64, size_infty)
    rho_u_infty_all = zeros(ComplexF64, size_infty)
    rho_v_infty_all = zeros(ComplexF64, size_infty)
    rho_E_infty_all = zeros(ComplexF64, size_infty)

    for i in 0:(N_f - 1)
        rho_infty_all[i*Nx+1:i*Nx+Nx]   = pert_infty_all[4*i*Nx+1       : 4*i*Nx+Nx]
        rho_u_infty_all[i*Nx+1:i*Nx+Nx] = pert_infty_all[4*i*Nx+Nx+1    : 4*i*Nx+2*Nx]
        rho_v_infty_all[i*Nx+1:i*Nx+Nx] = pert_infty_all[4*i*Nx+2*Nx+1  : 4*i*Nx+3*Nx]
        rho_E_infty_all[i*Nx+1:i*Nx+Nx] = pert_infty_all[4*i*Nx+3*Nx+1  : 4*i*Nx+4*Nx]
    end

    ## Zero freestream region and superpose Fourier modes
    indices   = collect(1:Nx)
    shock_idx = s["shock"]["cell_indices"][indices, 1]
    shock_x   = s["shock"]["points_x"][indices, 1]

    # Make mutable copies promoted to ComplexF64 (freestream contributions are complex)
    pert_rho   = ComplexF64.(pert_rho)
    pert_rho_u = ComplexF64.(pert_rho_u)
    pert_rho_v = ComplexF64.(pert_rho_v)
    pert_rho_E = ComplexF64.(pert_rho_E)

    if solution_remesh
        pert_rho   = zeros(ComplexF64, Nx, Ny)
        pert_rho_u = zeros(ComplexF64, Nx, Ny)
        pert_rho_v = zeros(ComplexF64, Nx, Ny)
        pert_rho_E = zeros(ComplexF64, Nx, Ny)
    else
        for k in 1:Nx
            pert_rho[k, shock_idx[k]:end]   .= 0.0
            pert_rho_u[k, shock_idx[k]:end] .= 0.0
            pert_rho_v[k, shock_idx[k]:end] .= 0.0
            pert_rho_E[k, shock_idx[k]:end] .= 0.0
        end
    end

    factor = 1
    for i in 1:N_f
        t = T_plot
        for k in 1:Nx
            y_row = s["mesh"]["x"][k, shock_idx[k]:end] .- shock_x[k]
            phase = exp.(1im .* w_infty[i] .* (y_row .- t)) ./ factor

            pert_rho[k, shock_idx[k]:end]   .= pert_rho[k, shock_idx[k]:end] .+
                rho_infty_all[k + (i-1)*Nx] .* phase
            pert_rho_u[k, shock_idx[k]:end] .= pert_rho_u[k, shock_idx[k]:end] .+
                rho_u_infty_all[k + (i-1)*Nx] .* phase
            pert_rho_v[k, shock_idx[k]:end] .= pert_rho_v[k, shock_idx[k]:end] .+
                rho_v_infty_all[k + (i-1)*Nx] .* phase
            pert_rho_E[k, shock_idx[k]:end] .= pert_rho_E[k, shock_idx[k]:end] .+
                rho_E_infty_all[k + (i-1)*Nx] .* phase
        end
    end

    return (pert_rho, pert_rho_u, pert_rho_v, pert_rho_E)
end
