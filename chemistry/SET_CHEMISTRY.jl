using Printf
using Dates
using Interpolations
using NaturalNeighbours
using Statistics
using Plots
using LinearAlgebra

"""
    SET_CHEMISTRY(s)

Initialize and configure chemistry module for the solver.

Master setup function for the chemistry module. Loads gas property data
from equilibrium model files, creates interpolation fits for all
thermodynamic and transport properties, validates the fits, generates
diagnostic plots, and optionally loads non-equilibrium relaxation models.

# Arguments
- `s::Dict{String,Any}`: Solution structure containing configuration fields

# Returns
- `chemistry::Dict{String,Any}`: Fully initialized chemistry structure

# Part of: Hypersonics Stability Julia Solver - Chemistry Module
"""
function SET_CHEMISTRY(s::Dict{String,Any})

    ## Early return if chemistry is disabled
    if !s["chemistry"]["is_chemistry_enabled"]
        println("No chemistry enabled...")
        s["chemistry"]["is_chemistry_enabled"] = false
        s["chemistry"]["chemistry_type"] = "None"
        s["chemistry"]["chemical_equilibrium"] = false
        s["chemistry"]["chemistry_composition"] = "None"
        chemistry = Dict{String,Any}()
        chemistry["planet"] = "NONE"
        chemistry["gas_model"] = "NONE"
        return chemistry
    end

    ## Load equilibrium gas property data
    println("Loading chemistry data...")
    name = s["solver_dir"] * "chemistry/equilibrium_models/gas_properties_" * s["chemistry"]["chemistry_type"] * "_" * s["chemistry"]["chemistry_composition"] * ".txt"

    chemistry = Dict{String,Any}()
    chemistry = READ_GAS_PROPERTIES(chemistry, name)
    chemistry["detailed_output"] = false
    @printf("Initializing chemistry module with detailed output: %d\n", chemistry["detailed_output"])

    ## Load frozen gas property data for non-equilibrium shock jumps
    if !s["chemistry"]["chemical_equilibrium"]
        if s["chemistry"]["chemistry_type"] == "Chemical-RTV" || s["chemistry"]["chemistry_type"] == "Frozen-RTV"
            chemistry_frozen = "Frozen-RTV"
        elseif s["chemistry"]["chemistry_type"] == "Chemical-RTVE" || s["chemistry"]["chemistry_type"] == "Frozen-RTVE"
            chemistry_frozen = "Frozen-RTVE"
        end

        chemistry["frozen"] = Dict{String,Any}()
        name_frozen = s["solver_dir"] * "chemistry/equilibrium_models/gas_properties_" * chemistry_frozen * "_" * s["chemistry"]["chemistry_composition"] * ".txt"
        chemistry["frozen"] = READ_GAS_PROPERTIES(chemistry["frozen"], name_frozen)
    end

    ## Fit interpolation models
    println("\nFitting chemistry models...")
    chemistry = FIT_CHEMISTRY(chemistry)
    if !s["chemistry"]["chemical_equilibrium"]
        chemistry["frozen"] = FIT_CHEMISTRY(chemistry["frozen"])
    end

    ## Run diagnostics
    if chemistry["detailed_output"]
        TEST_FIT_PERFORMANCE(chemistry)
        PLOT_FIT_EVALUATION(chemistry)

        println("\nCreating 3D surface plots...")
        PLOT_3D_CHEMISTRY_SURFACES(chemistry)
        TEST_RANKINE_HUGONIOT_SOLVER(chemistry, s)
    end

    ## Load non-equilibrium relaxation model
    if !s["chemistry"]["chemical_equilibrium"]
        if s["chemistry"]["chemistry_type"] == "Frozen-RTV"
            error("Frozen-RTV: frozen chemistry cannot have non equilibrium model")
        end

        println()
        println("Loading non-equilibrium model...")
        name = s["solver_dir"] * "chemistry/non_equilibrium_models/non_equi_model_" * s["chemistry"]["chemistry_composition"] * ".txt"
        chemistry = READ_NON_EQUILIBRIUM_MODEL(chemistry, name)
        println()
    end

    return chemistry
end


# ========================================================================
#  FIT_CHEMISTRY - Create interpolation models for all properties
# ========================================================================
"""
    FIT_CHEMISTRY(chemistry)

Creates fast interpolation models with log transforms.

Builds interpolation models for T, gamma_star, cv_star, mu, k, a, s and
log10 species concentrations as functions of (rho, e). Uses log10
transformations for density and energy, with all variables scaled to
[-1, 1] for numerical accuracy.
"""
function FIT_CHEMISTRY(chemistry::Dict{String,Any})
    println("Starting chemistry data fitting with log(rho) + scaled interpolation...")
    t_start = time()

    ## Validate input data
    if !haskey(chemistry, "rho") || !haskey(chemistry, "e")
        error("Chemistry structure must contain rho and E fields")
    end

    required_fields = ["T", "gamma_star", "cv_star", "mu", "k"]
    for field in required_fields
        if !haskey(chemistry, field)
            error("Chemistry structure must contain $field field")
        end
    end

    ## Extract data vectors
    rho = vec(chemistry["rho"])
    e = vec(chemistry["e"])
    T = vec(chemistry["T"])
    gamma_star = vec(chemistry["gamma_star"])
    cv_star = vec(chemistry["cv_star"])
    mu = vec(chemistry["mu"])
    k = vec(chemistry["k"])

    # Check for optional sound speed data
    if chemistry["has_sound_speed"]
        a = vec(chemistry["a"])
        has_a = true
        println("Sound speed data detected - will be included in fits.")
    else
        a = Float64[]
        has_a = false
        println("Sound speed data not available - skipping sound speed fits.")
    end

    # Check for optional entropy data
    if chemistry["has_entropy"]
        s = vec(chemistry["s"])
        has_s = true
        println("Entropy data detected - will be included in fits.")
    else
        s = Float64[]
        has_s = false
        println("Entropy data not available - skipping entropy fits.")
    end

    n_points = length(rho)
    @printf("Fitting %d data points...\n", n_points)

    ## Apply logarithmic transformations
    println("Applying logarithmic transformations: log10(rho), log10(e), log10(T)...")
    log_rho = log10.(rho)
    log_e = log10.(e)
    log_T = log10.(T)

    ## Remove invalid data points
    valid_idx = .!(isnan.(log_rho) .| isnan.(log_e) .| isnan.(log_T) .| isnan.(gamma_star) .|
                    isnan.(cv_star) .| isnan.(mu) .| isnan.(k) .|
                    isinf.(log_rho) .| isinf.(log_e) .| isinf.(log_T) .| isinf.(gamma_star) .|
                    isinf.(cv_star) .| isinf.(mu) .| isinf.(k))

    if has_a
        valid_idx .&= .!(isnan.(a) .| isinf.(a))
    end
    if has_s
        valid_idx .&= .!(isnan.(s) .| isinf.(s))
    end

    if sum(valid_idx) < n_points
        @printf("Warning: Removed %d invalid data points for basic properties\n", n_points - sum(valid_idx))
        rho = rho[valid_idx]
        log_rho = log_rho[valid_idx]
        e = e[valid_idx]
        log_e = log_e[valid_idx]
        T = T[valid_idx]
        log_T = log_T[valid_idx]
        gamma_star = gamma_star[valid_idx]
        cv_star = cv_star[valid_idx]
        mu = mu[valid_idx]
        k = k[valid_idx]
        if has_a
            a = a[valid_idx]
        end
        if has_s
            s = s[valid_idx]
        end
        n_points = length(rho)
    end

    ## Store data ranges for scaling
    if !haskey(chemistry, "fit_info")
        chemistry["fit_info"] = Dict{String,Any}()
    end
    chemistry["fit_info"]["rho_range"] = [minimum(rho), maximum(rho)]
    chemistry["fit_info"]["log_rho_range"] = [minimum(log_rho), maximum(log_rho)]
    chemistry["fit_info"]["e_range"] = [minimum(e), maximum(e)]
    chemistry["fit_info"]["log_e_range"] = [minimum(log_e), maximum(log_e)]
    chemistry["fit_info"]["T_range"] = [minimum(T), maximum(T)]
    chemistry["fit_info"]["log_T_range"] = [minimum(log_T), maximum(log_T)]
    chemistry["fit_info"]["gamma_star_range"] = [minimum(gamma_star), maximum(gamma_star)]
    chemistry["fit_info"]["cv_star_range"] = [minimum(cv_star), maximum(cv_star)]
    chemistry["fit_info"]["mu_range"] = [minimum(mu), maximum(mu)]
    chemistry["fit_info"]["k_range"] = [minimum(k), maximum(k)]
    if has_a
        chemistry["fit_info"]["a_range"] = [minimum(a), maximum(a)]
    end
    if has_s
        chemistry["fit_info"]["s_range"] = [minimum(s), maximum(s)]
    end
    chemistry["fit_info"]["n_points"] = n_points

    ## Print data ranges
    println("Original data ranges:")
    @printf("  rho:        [%.3e, %.3e] kg/m^3\n", chemistry["fit_info"]["rho_range"]...)
    @printf("  log10(rho): [%.3f, %.3f]\n", chemistry["fit_info"]["log_rho_range"]...)
    @printf("  e:          [%.3e, %.3e] J/kg\n", chemistry["fit_info"]["e_range"]...)
    @printf("  log10(e):   [%.3f, %.3f]\n", chemistry["fit_info"]["log_e_range"]...)
    @printf("  T:          [%.1f, %.1f] K\n", chemistry["fit_info"]["T_range"]...)
    @printf("  log10(T):   [%.3f, %.3f]\n", chemistry["fit_info"]["log_T_range"]...)
    @printf("  gamma*:     [%.4f, %.4f]\n", chemistry["fit_info"]["gamma_star_range"]...)
    @printf("  cv*:        [%.1f, %.1f] J/kg/K\n", chemistry["fit_info"]["cv_star_range"]...)
    @printf("  mu:         [%.2e, %.2e] Pa*s\n", chemistry["fit_info"]["mu_range"]...)
    @printf("  k:          [%.4f, %.4f] W/m/K\n", chemistry["fit_info"]["k_range"]...)
    if has_a
        @printf("  a:          [%.1f, %.1f] m/s\n", chemistry["fit_info"]["a_range"]...)
    end
    if has_s
        @printf("  s:          [%.1f, %.1f] J/kg/K\n", chemistry["fit_info"]["s_range"]...)
    end

    ## Scale all data to [-1, 1]
    if chemistry["detailed_output"]
        println("\nScaling all data to [-1, 1] for improved interpolation...")
    end

    # Scale input variables
    log_rho_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_rho, chemistry["fit_info"]["log_rho_range"])
    log_e_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_e, chemistry["fit_info"]["log_e_range"])

    # Scale output variables
    T_scaled = SCALE_TO_MINUS_ONE_TO_ONE(T, chemistry["fit_info"]["T_range"])
    gamma_star_scaled = SCALE_TO_MINUS_ONE_TO_ONE(gamma_star, chemistry["fit_info"]["gamma_star_range"])
    cv_star_scaled = SCALE_TO_MINUS_ONE_TO_ONE(cv_star, chemistry["fit_info"]["cv_star_range"])
    mu_scaled = SCALE_TO_MINUS_ONE_TO_ONE(mu, chemistry["fit_info"]["mu_range"])
    k_scaled = SCALE_TO_MINUS_ONE_TO_ONE(k, chemistry["fit_info"]["k_range"])
    if has_a
        a_scaled = SCALE_TO_MINUS_ONE_TO_ONE(a, chemistry["fit_info"]["a_range"])
    end
    if has_s
        s_scaled = SCALE_TO_MINUS_ONE_TO_ONE(s, chemistry["fit_info"]["s_range"])
    end

    ## Create regular grid for interpolants
    n_grid_log_rho = 100
    n_grid_log_e = 100

    if chemistry["detailed_output"]
        @printf("\nCreating regular grid (%dx%d = %d points)...\n", n_grid_log_rho, n_grid_log_e, n_grid_log_rho * n_grid_log_e)
    end

    log_rho_grid_vec = range(-1, 1, length=n_grid_log_rho)
    log_e_grid_vec = range(-1, 1, length=n_grid_log_e)

    ## Interpolate scattered data onto regular grid using NaturalNeighbours Farin C1
    # Analogous to MATLAB's griddata(..., 'cubic') - C1 Bézier on Delaunay triangles
    println("\nCreating interpolation objects for each property...")

    # Build query grid in ndgrid ordering (first dim = log_e, second dim = log_rho)
    log_e_query = vec(Float64[log_e_grid_vec[i] for i in 1:n_grid_log_e, j in 1:n_grid_log_rho])
    log_rho_query = vec(Float64[log_rho_grid_vec[j] for i in 1:n_grid_log_e, j in 1:n_grid_log_rho])

    T_grid = reshape(NaturalNeighbours.interpolate(log_e_scaled, log_rho_scaled, T_scaled; derivatives=true)(
        log_e_query, log_rho_query; method=Farin()), n_grid_log_e, n_grid_log_rho)
    println("  T: done")

    gamma_star_grid = reshape(NaturalNeighbours.interpolate(log_e_scaled, log_rho_scaled, gamma_star_scaled; derivatives=true)(
        log_e_query, log_rho_query; method=Farin()), n_grid_log_e, n_grid_log_rho)
    println("  gamma_star: done")

    cv_star_grid = reshape(NaturalNeighbours.interpolate(log_e_scaled, log_rho_scaled, cv_star_scaled; derivatives=true)(
        log_e_query, log_rho_query; method=Farin()), n_grid_log_e, n_grid_log_rho)
    println("  cv_star: done")

    mu_grid = reshape(NaturalNeighbours.interpolate(log_e_scaled, log_rho_scaled, mu_scaled; derivatives=true)(
        log_e_query, log_rho_query; method=Farin()), n_grid_log_e, n_grid_log_rho)
    println("  mu: done")

    k_grid = reshape(NaturalNeighbours.interpolate(log_e_scaled, log_rho_scaled, k_scaled; derivatives=true)(
        log_e_query, log_rho_query; method=Farin()), n_grid_log_e, n_grid_log_rho)
    println("  k: done")

    if has_a
        a_grid = reshape(NaturalNeighbours.interpolate(log_e_scaled, log_rho_scaled, a_scaled; derivatives=true)(
            log_e_query, log_rho_query; method=Farin()), n_grid_log_e, n_grid_log_rho)
        println("  a (sound speed): done")
    end

    if has_s
        s_grid = reshape(NaturalNeighbours.interpolate(log_e_scaled, log_rho_scaled, s_scaled; derivatives=true)(
            log_e_query, log_rho_query; method=Farin()), n_grid_log_e, n_grid_log_rho)
        println("  s (entropy): done")
    end

    ## Create interpolation objects using Interpolations.jl
    chemistry["fit_T"] = Interpolations.extrapolate(
        Interpolations.scale(Interpolations.interpolate(T_grid, BSpline(Cubic(Line(OnGrid())))), log_e_grid_vec, log_rho_grid_vec),
        Line())
    chemistry["fit_gamma_star"] = Interpolations.extrapolate(
        Interpolations.scale(Interpolations.interpolate(gamma_star_grid, BSpline(Cubic(Line(OnGrid())))), log_e_grid_vec, log_rho_grid_vec),
        Line())
    chemistry["fit_cv_star"] = Interpolations.extrapolate(
        Interpolations.scale(Interpolations.interpolate(cv_star_grid, BSpline(Cubic(Line(OnGrid())))), log_e_grid_vec, log_rho_grid_vec),
        Line())
    chemistry["fit_mu"] = Interpolations.extrapolate(
        Interpolations.scale(Interpolations.interpolate(mu_grid, BSpline(Cubic(Line(OnGrid())))), log_e_grid_vec, log_rho_grid_vec),
        Line())
    chemistry["fit_k"] = Interpolations.extrapolate(
        Interpolations.scale(Interpolations.interpolate(k_grid, BSpline(Cubic(Line(OnGrid())))), log_e_grid_vec, log_rho_grid_vec),
        Line())
    if has_a
        chemistry["fit_a"] = Interpolations.extrapolate(
            Interpolations.scale(Interpolations.interpolate(a_grid, BSpline(Cubic(Line(OnGrid())))), log_e_grid_vec, log_rho_grid_vec),
            Line())
    end
    if has_s
        chemistry["fit_s"] = Interpolations.extrapolate(
            Interpolations.scale(Interpolations.interpolate(s_grid, BSpline(Cubic(Line(OnGrid())))), log_e_grid_vec, log_rho_grid_vec),
            Line())
    end

    ## Create inverse interpolant: e(T, rho)
    println("  e: done")

    log_T_grid_vec = range(-1, 1, length=n_grid_log_e)
    log_rho_grid_vec_inv = range(-1, 1, length=n_grid_log_rho)

    log_T_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_T, chemistry["fit_info"]["log_T_range"])
    log_e_scaled_data = SCALE_TO_MINUS_ONE_TO_ONE(log_e, chemistry["fit_info"]["log_e_range"])

    log_T_query = vec(Float64[log_T_grid_vec[i] for i in 1:length(log_T_grid_vec), j in 1:length(log_rho_grid_vec_inv)])
    log_rho_inv_query = vec(Float64[log_rho_grid_vec_inv[j] for i in 1:length(log_T_grid_vec), j in 1:length(log_rho_grid_vec_inv)])

    log_e_grid_inv = reshape(NaturalNeighbours.interpolate(log_T_scaled, log_rho_scaled, log_e_scaled_data; derivatives=true)(
        log_T_query, log_rho_inv_query; method=Farin()), length(log_T_grid_vec), length(log_rho_grid_vec_inv))

    chemistry["fit_log_e_from_T_rho"] = Interpolations.extrapolate(
        Interpolations.scale(Interpolations.interpolate(log_e_grid_inv, BSpline(Cubic(Line(OnGrid())))), log_T_grid_vec, log_rho_grid_vec_inv),
        Line())

    ## Create species interpolants
    if chemistry["has_molar_data"]
        if chemistry["detailed_output"]
            println("\nSpecies data detected - creating interpolants for species concentrations...")
        end

        chemistry["species_fits"] = Dict{String,Any}()
        chemistry["species_ranges"] = Dict{String,Any}()

        for species in chemistry["species_list"]
            log10_field = "log10_" * species

            if haskey(chemistry, log10_field)
                log10_conc = chemistry[log10_field][valid_idx]

                # Additional filtering for species data
                species_valid_idx = .!(isnan.(log10_conc) .| isinf.(log10_conc) .| (log10_conc .< -50))

                if sum(species_valid_idx) > 10
                    range_name = "range_log10_" * species
                    chemistry["species_ranges"][range_name] = [minimum(log10_conc[species_valid_idx]), maximum(log10_conc[species_valid_idx])]

                    log10_conc_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log10_conc[species_valid_idx], chemistry["species_ranges"][range_name])

                    log10_conc_grid = reshape(NaturalNeighbours.interpolate(
                        log_e_scaled[species_valid_idx], log_rho_scaled[species_valid_idx], log10_conc_scaled;
                        derivatives=true)(log_e_query, log_rho_query; method=Farin()), n_grid_log_e, n_grid_log_rho)

                    fit_name = "fit_log10_" * species
                    chemistry["species_fits"][fit_name] = Interpolations.extrapolate(
                        Interpolations.scale(Interpolations.interpolate(log10_conc_grid, BSpline(Cubic(Line(OnGrid())))), log_e_grid_vec, log_rho_grid_vec),
                        Line())
                    @printf("  X_%s: done\n", species)

                    if chemistry["detailed_output"]
                        @printf("  X_%s: %d valid points, range [%.2f, %.2f]\n",
                                species, sum(species_valid_idx), chemistry["species_ranges"][range_name]...)
                    end
                else
                    @printf("    %s: insufficient data (%d points), skipping\n", species, sum(species_valid_idx))
                end
            end
        end
    end

    ## Store grid information
    chemistry["fit_info"]["grid_size"] = [n_grid_log_e, n_grid_log_rho]
    chemistry["fit_info"]["log_rho_grid_vec"] = log_rho_grid_vec
    chemistry["fit_info"]["log_e_grid_vec"] = log_e_grid_vec

    ## Create evaluation function handles (closures)
    chemistry["eval_T"] = (rho_new, e_new) -> EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(
        rho_new, e_new, chemistry["fit_T"], chemistry["fit_info"]["log_rho_range"],
        chemistry["fit_info"]["log_e_range"], chemistry["fit_info"]["T_range"])

    chemistry["eval_gamma_star"] = (rho_new, e_new) -> EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(
        rho_new, e_new, chemistry["fit_gamma_star"], chemistry["fit_info"]["log_rho_range"],
        chemistry["fit_info"]["log_e_range"], chemistry["fit_info"]["gamma_star_range"])

    chemistry["eval_cv_star"] = (rho_new, e_new) -> EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(
        rho_new, e_new, chemistry["fit_cv_star"], chemistry["fit_info"]["log_rho_range"],
        chemistry["fit_info"]["log_e_range"], chemistry["fit_info"]["cv_star_range"])

    chemistry["eval_mu"] = (rho_new, e_new) -> EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(
        rho_new, e_new, chemistry["fit_mu"], chemistry["fit_info"]["log_rho_range"],
        chemistry["fit_info"]["log_e_range"], chemistry["fit_info"]["mu_range"])

    chemistry["eval_k"] = (rho_new, e_new) -> EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(
        rho_new, e_new, chemistry["fit_k"], chemistry["fit_info"]["log_rho_range"],
        chemistry["fit_info"]["log_e_range"], chemistry["fit_info"]["k_range"])

    if has_a
        chemistry["eval_a"] = (rho_new, e_new) -> EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(
            rho_new, e_new, chemistry["fit_a"], chemistry["fit_info"]["log_rho_range"],
            chemistry["fit_info"]["log_e_range"], chemistry["fit_info"]["a_range"])
    end
    if has_s
        chemistry["eval_s"] = (rho_new, e_new) -> EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(
            rho_new, e_new, chemistry["fit_s"], chemistry["fit_info"]["log_rho_range"],
            chemistry["fit_info"]["log_e_range"], chemistry["fit_info"]["s_range"])
    end

    # Inverse evaluation function: E(T, rho)
    chemistry["eval_e"] = (T_new, rho_new) -> EVAL_GRIDDED_INTERPOLANT_INVERSE_LOG_INPUTS(
        T_new, rho_new, chemistry["fit_log_e_from_T_rho"], chemistry["fit_info"]["log_T_range"],
        chemistry["fit_info"]["log_rho_range"], chemistry["fit_info"]["log_e_range"],
        chemistry["fit_info"]["e_range"])

    ## Create species evaluation functions
    if chemistry["has_molar_data"]
        chemistry["species_eval_info"] = []

        for species in chemistry["species_list"]
            fit_name = "fit_log10_" * species
            eval_name = "eval_log10_" * species
            range_name = "range_log10_" * species

            if haskey(chemistry["species_fits"], fit_name) && haskey(chemistry["species_ranges"], range_name)
                eval_info = Dict{String,Any}()
                eval_info["species"] = species
                eval_info["fit_object"] = chemistry["species_fits"][fit_name]
                eval_info["log_rho_range"] = chemistry["fit_info"]["log_rho_range"]
                eval_info["log_e_range"] = chemistry["fit_info"]["log_e_range"]
                eval_info["species_range"] = chemistry["species_ranges"][range_name]

                push!(chemistry["species_eval_info"], eval_info)

                # Create closure for this species
                local sp = species
                chemistry[eval_name] = (rho_new, e_new) -> EVAL_SPECIES_GRIDDED_LOG_INPUTS(chemistry, sp, rho_new, e_new)
            end
        end

        chemistry["eval_all_species"] = (rho_new, e_new) -> EVALUATE_ALL_SPECIES(chemistry, rho_new, e_new)
    end

    # Vectorized evaluation for all basic properties
    chemistry["eval_all"] = (rho_new, e_new) -> EVALUATE_ALL_PROPERTIES(chemistry, rho_new, e_new, has_a)

    ## Validate fits on original data
    if chemistry["detailed_output"]
        println("Validating fits on original data...")
        chemistry["fit_validation"] = VALIDATE_FITS(chemistry, rho, e, T, gamma_star, cv_star, mu, k, a, s)
    end

    ## Store fitting metadata
    chemistry["fit_info"]["griddata_method"] = "NaturalNeighbours Farin C1 (Delaunay)"
    chemistry["fit_info"]["interp_method"] = "Cubic BSpline"
    chemistry["fit_info"]["extrapolation"] = "Linear"
    chemistry["fit_info"]["fit_time"] = time() - t_start
    chemistry["fit_info"]["date_fitted"] = string(Dates.now())
    chemistry["fit_info"]["has_sound_speed_fit"] = has_a
    chemistry["fit_info"]["has_entropy_fit"] = has_s
    chemistry["fit_info"]["interpolant_type"] = "Interpolations.jl_with_log_rho_log_e_log_T"
    chemistry["fit_info"]["uses_log_rho"] = true
    chemistry["fit_info"]["uses_log_e"] = true
    chemistry["fit_info"]["uses_log_T"] = true

    ## Print validation summary
    if chemistry["detailed_output"]
        @printf("Fitting completed in %.2f seconds\n", chemistry["fit_info"]["fit_time"])
        println("Validation results for basic properties:")
        for field in keys(chemistry["fit_validation"])
            if occursin("rmse", field)
                @printf("  %s: %.6e\n", field, chemistry["fit_validation"][field])
            elseif occursin("r2", field)
                @printf("  %s: %.6f\n", field, chemistry["fit_validation"][field])
            end
        end

        if chemistry["has_molar_data"]
            @printf("Species fitting completed for: %s\n",
                    join(keys(chemistry["species_fits"]), ", "))
        end

        @printf("\nData approach: Scattered data interpolated onto %dx%d regular grid.\n", n_grid_log_e, n_grid_log_rho)
        println("Data transformation: log10(rho), log10(e), log10(T) applied before [-1,1] scaling.")
        println("Data scaling: All variables scaled to [-1,1] for improved numerical accuracy.")
        @printf("Using Interpolations.jl with %s interpolation for fast, smooth evaluation.\n", "Cubic BSpline")
        println("Inverse fit E(T,rho) created in log-space for Rankine-Hugoniot solver integration.")
        if has_s
            println("Entropy s(rho,E) interpolation included.")
        end
        println("\nNote: log transformations improve interpolation over wide ranges.")
    end

    return chemistry
end


# ========================================================================
#  EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS - Evaluate with log transforms
# ========================================================================
"""
    EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(rho_new, e_new, interpolant, log_rho_range, log_e_range, values_range)

Evaluate property using log-space interpolation. Transforms query points to
log10 space, scales to [-1, 1], evaluates the interpolant, and unscales
the result to original units.
"""
function EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(rho_new, e_new, interpolant, log_rho_range, log_e_range, values_range)
    # Wrap scalars into 1-element arrays so boolean indexing works uniformly
    scalar_input = isa(rho_new, Number)
    if scalar_input
        rho_new = [rho_new]
        e_new   = [e_new]
    end

    # Handle non-positive values with extrapolation instead of NaN
    result = zeros(size(rho_new))
    valid_mask = (rho_new .> 0) .& (e_new .> 0)

    if !any(valid_mask)
        # All values are non-positive - use minimum log ranges
        log_rho_min = log_rho_range[1]
        log_e_min = log_e_range[1]
        log_rho_new_extrap = fill(log_rho_min, size(rho_new))
        log_e_new_extrap = fill(log_e_min, size(e_new))
        log_rho_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_rho_new_extrap, log_rho_range)
        log_e_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_e_new_extrap, log_e_range)
        result_scaled = interpolant.(log_e_new_scaled, log_rho_new_scaled)
        result = SCALE_FROM_MINUS_ONE_TO_ONE(result_scaled, values_range)
        return scalar_input ? result[1] : result
    end

    # Process valid positive rho and e values
    if any(valid_mask)
        log_rho_new = log10.(rho_new[valid_mask])
        log_e_new = log10.(e_new[valid_mask])

        log_rho_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_rho_new, log_rho_range)
        log_e_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_e_new, log_e_range)

        result_scaled = interpolant.(log_e_new_scaled, log_rho_new_scaled)
        result[valid_mask] .= SCALE_FROM_MINUS_ONE_TO_ONE(result_scaled, values_range)
    end

    # For invalid values, extrapolate using minimum log values
    if any(.!valid_mask)
        log_rho_min = log_rho_range[1]
        log_e_min = log_e_range[1]

        n_invalid = sum(.!valid_mask)
        log_rho_extrap = fill(log_rho_min, n_invalid)
        log_e_extrap = fill(log_e_min, n_invalid)

        # Override with actual log values where possible
        invalid_rho = rho_new[.!valid_mask]
        invalid_e = e_new[.!valid_mask]
        valid_rho_only = invalid_rho .> 0
        valid_e_only = invalid_e .> 0

        if any(valid_rho_only)
            log_rho_extrap[valid_rho_only] .= log10.(invalid_rho[valid_rho_only])
        end
        if any(valid_e_only)
            log_e_extrap[valid_e_only] .= log10.(invalid_e[valid_e_only])
        end

        log_rho_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_rho_extrap, log_rho_range)
        log_e_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_e_extrap, log_e_range)

        result_scaled = interpolant.(log_e_new_scaled, log_rho_new_scaled)
        result[.!valid_mask] .= SCALE_FROM_MINUS_ONE_TO_ONE(result_scaled, values_range)
    end

    return scalar_input ? result[1] : result
end


# ========================================================================
#  EVAL_GRIDDED_INTERPOLANT_INVERSE_LOG_INPUTS - Inverse evaluation e(T, rho)
# ========================================================================
"""
    EVAL_GRIDDED_INTERPOLANT_INVERSE_LOG_INPUTS(T_new, rho_new, interpolant, log_T_range, log_rho_range, log_e_range, e_range)

Evaluate e(T, rho) in log-space. Returns energy in original units.
"""
function EVAL_GRIDDED_INTERPOLANT_INVERSE_LOG_INPUTS(T_new, rho_new, interpolant, log_T_range, log_rho_range, log_e_range, e_range)
    # Wrap scalars into 1-element arrays so boolean indexing works uniformly
    scalar_input = isa(rho_new, Number)
    if scalar_input
        T_new   = [T_new]
        rho_new = [rho_new]
    end

    # Handle non-positive values with extrapolation
    result = zeros(size(rho_new))
    valid_mask = (T_new .> 0) .& (rho_new .> 0)

    if !any(valid_mask)
        log_T_min = log_T_range[1]
        log_rho_min = log_rho_range[1]
        log_T_new_extrap = fill(log_T_min, size(T_new))
        log_rho_new_extrap = fill(log_rho_min, size(rho_new))
        log_T_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_T_new_extrap, log_T_range)
        log_rho_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_rho_new_extrap, log_rho_range)
        log_e_result_scaled = interpolant.(log_T_new_scaled, log_rho_new_scaled)
        log_e_result = SCALE_FROM_MINUS_ONE_TO_ONE(log_e_result_scaled, log_e_range)
        result = 10.0 .^ log_e_result
        return scalar_input ? result[1] : result
    end
    if any(valid_mask)
        log_T_new = log10.(T_new[valid_mask])
        log_rho_new = log10.(rho_new[valid_mask])

        log_T_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_T_new, log_T_range)
        log_rho_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_rho_new, log_rho_range)

        log_e_result_scaled = interpolant.(log_T_new_scaled, log_rho_new_scaled)
        log_e_result = SCALE_FROM_MINUS_ONE_TO_ONE(log_e_result_scaled, log_e_range)
        result[valid_mask] .= 10.0 .^ log_e_result
    end

    # For invalid values, extrapolate using minimum log values
    if any(.!valid_mask)
        log_T_min = log_T_range[1]
        log_rho_min = log_rho_range[1]

        n_invalid = sum(.!valid_mask)
        log_T_extrap = fill(log_T_min, n_invalid)
        log_rho_extrap = fill(log_rho_min, n_invalid)

        invalid_T = T_new[.!valid_mask]
        invalid_rho = rho_new[.!valid_mask]
        valid_T_only = invalid_T .> 0
        valid_rho_only = invalid_rho .> 0

        if any(valid_T_only)
            log_T_extrap[valid_T_only] .= log10.(invalid_T[valid_T_only])
        end
        if any(valid_rho_only)
            log_rho_extrap[valid_rho_only] .= log10.(invalid_rho[valid_rho_only])
        end

        log_T_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_T_extrap, log_T_range)
        log_rho_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_rho_extrap, log_rho_range)

        log_e_result_scaled = interpolant.(log_T_new_scaled, log_rho_new_scaled)
        log_e_result = SCALE_FROM_MINUS_ONE_TO_ONE(log_e_result_scaled, log_e_range)
        result[.!valid_mask] .= 10.0 .^ log_e_result
    end

    return scalar_input ? result[1] : result
end


# ========================================================================
#  EVAL_SPECIES_GRIDDED_LOG_INPUTS - Evaluate species concentration
# ========================================================================
"""
    EVAL_SPECIES_GRIDDED_LOG_INPUTS(chemistry, species_name, rho_new, e_new)

Evaluate log10 species concentration.
"""
function EVAL_SPECIES_GRIDDED_LOG_INPUTS(chemistry::Dict{String,Any}, species_name::String, rho_new, e_new)
    for eval_info in chemistry["species_eval_info"]
        if eval_info["species"] == species_name
            result = EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(rho_new, e_new, eval_info["fit_object"],
                eval_info["log_rho_range"], eval_info["log_e_range"], eval_info["species_range"])
            return result
        end
    end
    # Species not found
    return fill(NaN, size(rho_new))
end


# ========================================================================
#  SCALE_TO_MINUS_ONE_TO_ONE - Scale data to [-1, 1]
# ========================================================================
"""
    SCALE_TO_MINUS_ONE_TO_ONE(data, data_range)

Linear scaling from original range to [-1, 1].
"""
function SCALE_TO_MINUS_ONE_TO_ONE(data, data_range)
    if data_range[2] == data_range[1]
        return zeros(size(data))
    else
        return 2.0 .* (data .- data_range[1]) ./ (data_range[2] - data_range[1]) .- 1.0
    end
end


# ========================================================================
#  SCALE_FROM_MINUS_ONE_TO_ONE - Unscale data from [-1, 1]
# ========================================================================
"""
    SCALE_FROM_MINUS_ONE_TO_ONE(scaled_data, data_range)

Linear unscaling from [-1, 1] to original range.
"""
function SCALE_FROM_MINUS_ONE_TO_ONE(scaled_data, data_range)
    if data_range[2] == data_range[1]
        return fill(data_range[1], size(scaled_data))
    else
        return (scaled_data .+ 1.0) .* (data_range[2] - data_range[1]) ./ 2.0 .+ data_range[1]
    end
end


# ========================================================================
#  EVALUATE_ALL_PROPERTIES - Batch evaluate all basic properties
# ========================================================================
"""
    EVALUATE_ALL_PROPERTIES(chemistry, rho_new, e_new, has_a)

Evaluate T, gamma*, cv*, mu, k (and a) at once.
"""
function EVALUATE_ALL_PROPERTIES(chemistry::Dict{String,Any}, rho_new, e_new, has_a::Bool)
    result = Dict{String,Any}()
    result["T"] = chemistry["eval_T"](rho_new, e_new)
    result["gamma_star"] = chemistry["eval_gamma_star"](rho_new, e_new)
    result["cv_star"] = chemistry["eval_cv_star"](rho_new, e_new)
    result["mu"] = chemistry["eval_mu"](rho_new, e_new)
    result["k"] = chemistry["eval_k"](rho_new, e_new)
    if has_a
        result["a"] = chemistry["eval_a"](rho_new, e_new)
    end
    return result
end


# ========================================================================
#  EVALUATE_ALL_SPECIES - Batch evaluate all species concentrations
# ========================================================================
"""
    EVALUATE_ALL_SPECIES(chemistry, rho_new, e_new)

Evaluate all species log10 concentrations at once.
"""
function EVALUATE_ALL_SPECIES(chemistry::Dict{String,Any}, rho_new, e_new)
    result = Dict{String,Any}()
    for species in chemistry["species_list"]
        eval_name = "eval_log10_" * species
        if haskey(chemistry, eval_name)
            result[species] = chemistry[eval_name](rho_new, e_new)
        else
            result[species] = NaN
        end
    end
    return result
end


# ========================================================================
#  VALIDATE_FITS - Compute error metrics for all fitted properties
# ========================================================================
"""
    VALIDATE_FITS(chemistry, rho, e, T_true, gamma_star_true, cv_star_true, mu_true, k_true, a_true, s_true)

Validate interpolation fits against original data. Computes RMSE, MAE, and R^2.
"""
function VALIDATE_FITS(chemistry::Dict{String,Any}, rho, e, T_true, gamma_star_true, cv_star_true, mu_true, k_true, a_true, s_true)
    ## Predict using fits
    T_pred = chemistry["eval_T"](rho, e)
    gamma_star_pred = chemistry["eval_gamma_star"](rho, e)
    cv_star_pred = chemistry["eval_cv_star"](rho, e)
    mu_pred = chemistry["eval_mu"](rho, e)
    k_pred = chemistry["eval_k"](rho, e)
    e_pred = chemistry["eval_e"](T_true, rho)

    ## Create valid-data masks
    valid_T = .!isnan.(T_pred)
    valid_gamma_star = .!isnan.(gamma_star_pred)
    valid_cv_star = .!isnan.(cv_star_pred)
    valid_mu = .!isnan.(mu_pred)
    valid_k = .!isnan.(k_pred)
    valid_e = .!isnan.(e_pred)

    # Helper for safe total sum of squares
    safe_sst(true_data, valid_mask) = sum((true_data[valid_mask] .- mean(true_data[valid_mask])).^2)
    sst_threshold = 1e-12

    validation = Dict{String,Any}()

    ## Temperature validation
    validation["T_rmse"] = sqrt(mean((T_pred[valid_T] .- T_true[valid_T]).^2))
    validation["T_mae"] = mean(abs.(T_pred[valid_T] .- T_true[valid_T]))
    sst_T = safe_sst(T_true, valid_T)
    validation["T_r2"] = sst_T < sst_threshold ? 1.0 : 1.0 - sum((T_pred[valid_T] .- T_true[valid_T]).^2) / sst_T

    ## gamma_star validation
    validation["gamma_star_rmse"] = sqrt(mean((gamma_star_pred[valid_gamma_star] .- gamma_star_true[valid_gamma_star]).^2))
    validation["gamma_star_mae"] = mean(abs.(gamma_star_pred[valid_gamma_star] .- gamma_star_true[valid_gamma_star]))
    sst_gs = safe_sst(gamma_star_true, valid_gamma_star)
    validation["gamma_star_r2"] = sst_gs < sst_threshold ? 1.0 : 1.0 - sum((gamma_star_pred[valid_gamma_star] .- gamma_star_true[valid_gamma_star]).^2) / sst_gs

    ## cv_star validation
    validation["cv_star_rmse"] = sqrt(mean((cv_star_pred[valid_cv_star] .- cv_star_true[valid_cv_star]).^2))
    validation["cv_star_mae"] = mean(abs.(cv_star_pred[valid_cv_star] .- cv_star_true[valid_cv_star]))
    sst_cv = safe_sst(cv_star_true, valid_cv_star)
    validation["cv_star_r2"] = sst_cv < sst_threshold ? 1.0 : 1.0 - sum((cv_star_pred[valid_cv_star] .- cv_star_true[valid_cv_star]).^2) / sst_cv

    ## Viscosity validation
    validation["mu_rmse"] = sqrt(mean((mu_pred[valid_mu] .- mu_true[valid_mu]).^2))
    validation["mu_mae"] = mean(abs.(mu_pred[valid_mu] .- mu_true[valid_mu]))
    sst_mu = safe_sst(mu_true, valid_mu)
    validation["mu_r2"] = sst_mu < sst_threshold ? 1.0 : 1.0 - sum((mu_pred[valid_mu] .- mu_true[valid_mu]).^2) / sst_mu

    ## Thermal conductivity validation
    validation["k_rmse"] = sqrt(mean((k_pred[valid_k] .- k_true[valid_k]).^2))
    validation["k_mae"] = mean(abs.(k_pred[valid_k] .- k_true[valid_k]))
    sst_k = safe_sst(k_true, valid_k)
    validation["k_r2"] = sst_k < sst_threshold ? 1.0 : 1.0 - sum((k_pred[valid_k] .- k_true[valid_k]).^2) / sst_k

    ## Inverse fit E(T, rho) validation
    validation["e_rmse"] = sqrt(mean((e_pred[valid_e] .- e[valid_e]).^2))
    validation["e_mae"] = mean(abs.(e_pred[valid_e] .- e[valid_e]))
    sst_e = safe_sst(e, valid_e)
    validation["e_r2"] = sst_e < sst_threshold ? 1.0 : 1.0 - sum((e_pred[valid_e] .- e[valid_e]).^2) / sst_e

    ## Sound speed validation (if available)
    if !isempty(a_true)
        a_pred = chemistry["eval_a"](rho, e)
        valid_a = .!isnan.(a_pred)
        validation["a_rmse"] = sqrt(mean((a_pred[valid_a] .- a_true[valid_a]).^2))
        validation["a_mae"] = mean(abs.(a_pred[valid_a] .- a_true[valid_a]))
        sst_a = safe_sst(a_true, valid_a)
        validation["a_r2"] = sst_a < sst_threshold ? 1.0 : 1.0 - sum((a_pred[valid_a] .- a_true[valid_a]).^2) / sst_a
        if !haskey(validation, "predictions")
            validation["predictions"] = Dict{String,Any}()
        end
        validation["predictions"]["a"] = a_pred
    end

    ## Entropy validation (if available)
    if !isempty(s_true)
        println("\n=== Entropy Validation ===")

        # Method 1: Direct interpolation
        println("Method 1: Direct interpolation s(rho,e)...")
        s_pred_interp = chemistry["eval_s"](rho, e)
        valid_s_interp = .!isnan.(s_pred_interp)

        validation["s_interp_rmse"] = sqrt(mean((s_pred_interp[valid_s_interp] .- s_true[valid_s_interp]).^2))
        validation["s_interp_mae"] = mean(abs.(s_pred_interp[valid_s_interp] .- s_true[valid_s_interp]))
        sst_s_interp = safe_sst(s_true, valid_s_interp)
        validation["s_interp_r2"] = sst_s_interp < sst_threshold ? 1.0 : 1.0 - sum((s_pred_interp[valid_s_interp] .- s_true[valid_s_interp]).^2) / sst_s_interp

        @printf("  Interpolation R^2: %.6f\n", validation["s_interp_r2"])
        @printf("  Interpolation RMSE: %.4f J/kg/K\n", validation["s_interp_rmse"])

        # Use interpolation as primary
        validation["s_rmse"] = validation["s_interp_rmse"]
        validation["s_mae"] = validation["s_interp_mae"]
        validation["s_r2"] = validation["s_interp_r2"]

        println("===========================")

        if !haskey(validation, "predictions")
            validation["predictions"] = Dict{String,Any}()
        end
        validation["predictions"]["s"] = s_pred_interp
    end

    ## Store predicted values
    if !haskey(validation, "predictions")
        validation["predictions"] = Dict{String,Any}()
    end
    validation["predictions"]["T"] = T_pred
    validation["predictions"]["gamma_star"] = gamma_star_pred
    validation["predictions"]["cv_star"] = cv_star_pred
    validation["predictions"]["mu"] = mu_pred
    validation["predictions"]["k"] = k_pred
    validation["predictions"]["e"] = e_pred

    return validation
end


# ========================================================================
#  TEST_FIT_PERFORMANCE - Benchmark interpolation speed
# ========================================================================
"""
    TEST_FIT_PERFORMANCE(chemistry)

Test evaluation speed of the fitted models.
"""
function TEST_FIT_PERFORMANCE(chemistry::Dict{String,Any})
    if !haskey(chemistry, "fit_T")
        error("Chemistry data must be fitted first. Run FIT_CHEMISTRY(chemistry)")
    end

    println("\nTesting fit performance with scaled interpolation...")

    ## Generate test points
    n_test = 10000
    rho_range = chemistry["fit_info"]["rho_range"]
    e_range = chemistry["fit_info"]["e_range"]
    T_range = chemistry["fit_info"]["T_range"]

    rho_test = rho_range[1] .+ (rho_range[2] - rho_range[1]) .* rand(n_test)
    e_test = e_range[1] .+ (e_range[2] - e_range[1]) .* rand(n_test)
    T_test = T_range[1] .+ (T_range[2] - T_range[1]) .* rand(n_test)

    ## Time vectorized evaluations
    @printf("Timing %d evaluations:\n", n_test)

    t = @elapsed T_eval = chemistry["eval_T"](rho_test, e_test)
    @printf("  T evaluation: %.4f seconds (%.2e sec/point)\n", t, t / n_test)

    t = @elapsed gamma_star_test = chemistry["eval_gamma_star"](rho_test, e_test)
    @printf("  gamma_star evaluation: %.4f seconds (%.2e sec/point)\n", t, t / n_test)

    if chemistry["fit_info"]["has_sound_speed_fit"]
        t = @elapsed a_test = chemistry["eval_a"](rho_test, e_test)
        @printf("  Sound speed evaluation: %.4f seconds (%.2e sec/point)\n", t, t / n_test)
    end

    t = @elapsed all_props = chemistry["eval_all"](rho_test, e_test)
    @printf("  All basic properties: %.4f seconds (%.2e sec/point)\n", t, t / n_test)

    t = @elapsed e_inverse = chemistry["eval_e"](T_test, rho_test)
    @printf("  E(T,rho) inverse evaluation: %.4f seconds (%.2e sec/point)\n", t, t / n_test)

    if chemistry["has_molar_data"]
        t = @elapsed all_species = chemistry["eval_all_species"](rho_test, e_test)
        @printf("  All species: %.4f seconds (%.2e sec/point)\n", t, t / n_test)
    end

    ## Time single-point evaluations
    rho_single = mean(rho_range)
    e_single = mean(e_range)
    T_single = mean(T_range)

    t = @elapsed begin
        for i in 1:1000
            T_eval_single = chemistry["eval_T"](rho_single, e_single)
        end
    end
    @printf("  Single point (1000 calls): %.4f seconds (%.2e sec/call)\n", t, t / 1000)

    t = @elapsed begin
        for i in 1:1000
            e_inverse_single = chemistry["eval_e"](T_single, rho_single)
        end
    end
    @printf("  Single inverse E(T,rho) (1000 calls): %.4f seconds (%.2e sec/call)\n", t, t / 1000)

    ## Print summary
    println("\nInterpolations.jl with [-1,1] scaling provides fast, smooth interpolation.")
    @printf("Interpolant type: %s\n", chemistry["fit_info"]["interpolant_type"])
    @printf("Grid resolution: %dx%d points\n", chemistry["fit_info"]["grid_size"]...)
    @printf("Griddata method: %s\n", chemistry["fit_info"]["griddata_method"])
    @printf("Interpolation method: %s\n", chemistry["fit_info"]["interp_method"])
    println("Inverse E(T,rho) function ready for Rankine-Hugoniot solver.")
    if chemistry["fit_info"]["has_sound_speed_fit"]
        println("Sound speed interpolation included.")
    end
end


# ========================================================================
#  PLOT_FIT_EVALUATION - Validation scatter plots
# ========================================================================
"""
    PLOT_FIT_EVALUATION(chemistry)

Plot predicted vs. true values for all fitted properties.
"""
function PLOT_FIT_EVALUATION(chemistry::Dict{String,Any})
    if !haskey(chemistry, "fit_validation")
        error("Chemistry data must be fitted first. Run FIT_CHEMISTRY(chemistry)")
    end

    val = chemistry["fit_validation"]
    has_a = haskey(val, "a_rmse")
    has_s = haskey(val, "s_rmse")

    plots_list = []

    # Temperature
    p1 = scatter(chemistry["T"], val["predictions"]["T"], ms=1, label="", xlabel="True T (K)", ylabel="Predicted T (K)",
                 title=@sprintf("Temperature (R² = %.4f)", val["T_r2"]))
    plot!(p1, [minimum(chemistry["T"]), maximum(chemistry["T"])], [minimum(chemistry["T"]), maximum(chemistry["T"])], lw=2, color=:red, label="")
    push!(plots_list, p1)

    # gamma_star
    p2 = scatter(chemistry["gamma_star"], val["predictions"]["gamma_star"], ms=1, label="", xlabel="True γ*", ylabel="Predicted γ*",
                 title=@sprintf("Heat Cap. Ratio (R² = %.4f)", val["gamma_star_r2"]))
    plot!(p2, [minimum(chemistry["gamma_star"]), maximum(chemistry["gamma_star"])], [minimum(chemistry["gamma_star"]), maximum(chemistry["gamma_star"])], lw=2, color=:red, label="")
    push!(plots_list, p2)

    # cv_star
    p3 = scatter(chemistry["cv_star"], val["predictions"]["cv_star"], ms=1, label="", xlabel="True cv* (J/kg/K)", ylabel="Predicted cv*",
                 title=@sprintf("Specific Heat (R² = %.4f)", val["cv_star_r2"]))
    plot!(p3, [minimum(chemistry["cv_star"]), maximum(chemistry["cv_star"])], [minimum(chemistry["cv_star"]), maximum(chemistry["cv_star"])], lw=2, color=:red, label="")
    push!(plots_list, p3)

    # Viscosity
    p4 = scatter(chemistry["mu"], val["predictions"]["mu"], ms=1, label="", xlabel="True μ (Pa·s)", ylabel="Predicted μ",
                 title=@sprintf("Viscosity (R² = %.4f)", val["mu_r2"]))
    plot!(p4, [minimum(chemistry["mu"]), maximum(chemistry["mu"])], [minimum(chemistry["mu"]), maximum(chemistry["mu"])], lw=2, color=:red, label="")
    push!(plots_list, p4)

    # Thermal conductivity
    p5 = scatter(chemistry["k"], val["predictions"]["k"], ms=1, label="", xlabel="True k (W/m/K)", ylabel="Predicted k",
                 title=@sprintf("Thermal Cond. (R² = %.4f)", val["k_r2"]))
    plot!(p5, [minimum(chemistry["k"]), maximum(chemistry["k"])], [minimum(chemistry["k"]), maximum(chemistry["k"])], lw=2, color=:red, label="")
    push!(plots_list, p5)

    # Energy inverse fit
    p6 = scatter(chemistry["e"], val["predictions"]["e"], ms=1, label="", xlabel="True E (J/kg)", ylabel="Predicted E",
                 title=@sprintf("Energy E(T,ρ) (R² = %.4f)", val["e_r2"]))
    plot!(p6, [minimum(chemistry["e"]), maximum(chemistry["e"])], [minimum(chemistry["e"]), maximum(chemistry["e"])], lw=2, color=:red, label="")
    push!(plots_list, p6)

    if has_a
        p7 = scatter(chemistry["a"], val["predictions"]["a"], ms=1, label="", xlabel="True a (m/s)", ylabel="Predicted a",
                     title=@sprintf("Sound Speed (R² = %.4f)", val["a_r2"]))
        plot!(p7, [minimum(chemistry["a"]), maximum(chemistry["a"])], [minimum(chemistry["a"]), maximum(chemistry["a"])], lw=2, color=:red, label="")
        push!(plots_list, p7)
    end

    if has_s
        p8 = scatter(chemistry["s"], val["predictions"]["s"], ms=1, label="", xlabel="True s (J/kg/K)", ylabel="Predicted s",
                     title=@sprintf("Entropy (R² = %.4f)", val["s_r2"]))
        plot!(p8, [minimum(chemistry["s"]), maximum(chemistry["s"])], [minimum(chemistry["s"]), maximum(chemistry["s"])], lw=2, color=:red, label="")
        push!(plots_list, p8)
    end

    n_plots = length(plots_list)
    ncols = 3
    nrows = ceil(Int, n_plots / ncols)
    p_all = plot(plots_list..., layout=(nrows, ncols), size=(1200, 300*nrows))
    display(p_all)
end


# ========================================================================
#  PLOT_3D_CHEMISTRY_SURFACES - 3D surface visualization
# ========================================================================
"""
    PLOT_3D_CHEMISTRY_SURFACES(chemistry)

Create 3D surface plots of all chemistry properties.
"""
function PLOT_3D_CHEMISTRY_SURFACES(chemistry::Dict{String,Any})
    if !haskey(chemistry, "fit_T")
        error("Chemistry data must be fitted first. Run FIT_CHEMISTRY(chemistry)")
    end

    has_a = chemistry["fit_info"]["has_sound_speed_fit"]
    has_s = haskey(chemistry["fit_info"], "has_entropy_fit") && chemistry["fit_info"]["has_entropy_fit"]

    println("Creating 3D surface plots with original data points...")

    ## Get data ranges and create evaluation grids
    rho_range = chemistry["fit_info"]["rho_range"]
    e_range = chemistry["fit_info"]["e_range"]
    T_range = chemistry["fit_info"]["T_range"]

    n_points = 200

    e_grid = collect(range(e_range[1], e_range[2], length=n_points))
    rho_grid = collect(range(rho_range[1], rho_range[2], length=n_points))
    e_mesh = [e_grid[i] for i in 1:n_points, j in 1:n_points]
    rho_mesh = [rho_grid[j] for i in 1:n_points, j in 1:n_points]

    ## Evaluate properties on (E, rho) grid
    gamma_star_e_rho = reshape(chemistry["eval_gamma_star"](vec(rho_mesh), vec(e_mesh)), size(rho_mesh))
    T_e_rho = reshape(chemistry["eval_T"](vec(rho_mesh), vec(e_mesh)), size(rho_mesh))

    ## Figure: T(E, rho)
    p = surface(e_grid ./ 1e6, rho_grid, T_e_rho' ./ 1000,
                xlabel="E (MJ/kg)", ylabel="ρ (kg/m³)", zlabel="T (kK)",
                title="Temperature T(E,ρ)", colorbar=true)
    display(p)

    println("3D surface plots created successfully.")
    @printf("\nData coverage statistics:\n")
    @printf("  Temperature range: %.1f - %.1f K\n", T_range...)
    @printf("  Energy range: %.2e - %.2e J/kg\n", e_range...)
    @printf("  Density range: %.2e - %.2e kg/m^3\n", rho_range...)
end
