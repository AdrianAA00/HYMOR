# UPDATE_SHOCK_JUMP_PROPERTIES  Update post-shock flow properties via Rankine-Hugoniot relations.
#
#   s = UPDATE_SHOCK_JUMP_PROPERTIES(s, chemistry) solves the Rankine-Hugoniot
#   jump conditions across the fitted shock using a Newton-Raphson iteration
#   on the Riemann invariant formulation. The routine interpolates the
#   downstream flow field to the shock location, iterates for the post-shock
#   pressure, and back-computes all conservative variables and shock speeds.
#
#   Inputs:
#       s         - Solution structure containing the flow field, grid, shock
#                   geometry, and solver parameters.
#       chemistry - Chemistry model structure providing thermodynamic property
#                   evaluators (gamma_star, cv_star) for equilibrium or frozen
#                   chemistry.
#
#   Outputs:
#       s - Updated s structure with:
#             s["shock"]["properties"]  (rho, rho_u, rho_v, rho_E, p, gamma_star, cv_star)
#             s["shock"]["speed_x"], s["shock"]["speed_y"]
#             s["shock"]["relative_increase_velocity"]
#           and interpolated shocked-cell values in the flow field arrays.
#
#   Notes:
#       - Dispatches to UPDATE_SHOCK_JUMP_PROPERTIES_PRESSURE_3 (active method).
#       - Supports both Lagrangian and Eulerian shock-speed formulations.
#       - Handles calorically perfect gas, frozen chemistry, and equilibrium
#         chemistry via the chemistry argument.
#       - Ghost-cell extrapolation is performed via EXTRAPOLATE_CELLS_SHOCK.
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

function UPDATE_SHOCK_JUMP_PROPERTIES(s::Dict{String, Any}, chemistry::Dict{String, Any})
    s = UPDATE_SHOCK_JUMP_PROPERTIES_PRESSURE_3(s, chemistry)
    return s
end

# ========================================================================
#  Riemann-Invariant Rankine-Hugoniot Solver
# ========================================================================
function UPDATE_SHOCK_JUMP_PROPERTIES_RIEMANN(s::Dict{String, Any}, chemistry::Dict{String, Any})
# UPDATE_SHOCK_JUMP_PROPERTIES_RIEMANN  Solve R-H with Newton-Raphson on
#   the Riemann invariant (C- characteristic) to obtain post-shock pressure.

    ## Extract dict refs
    var   = s["var"]::Dict{String, Any}
    mesh  = s["mesh"]::Dict{String, Any}
    shock = s["shock"]::Dict{String, Any}
    chem  = s["chemistry"]::Dict{String, Any}
    fs    = s["freestream"]::Dict{String, Any}
    Nchi  = mesh["Nchi"]::Int

    rho_arr       = var["rho"]::Matrix{Float64}
    rho_u_arr     = var["rho_u"]::Matrix{Float64}
    rho_v_arr     = var["rho_v"]::Matrix{Float64}
    rho_E_arr     = var["rho_E"]::Matrix{Float64}
    p_arr         = var["p"]::Matrix{Float64}
    a_arr         = var["a"]::Matrix{Float64}
    x_mesh        = mesh["x"]::Matrix{Float64}
    y_mesh        = mesh["y"]::Matrix{Float64}

    chem_eq = chem["chemical_equilibrium"]::Bool
    is_chem = chem["is_chemistry_enabled"]::Bool
    rho_factor    = fs["rho_factor"]::Float64
    energy_factor = fs["energy_factor"]::Float64

    ## Upstream conditions
    p_infty, u_inf, v_inf, rho_inf, e_inf = COMPUTE_UPSTREAM_CONDITIONS(s)
    s = UPDATE_FIELD_UPSTREAM(s)

    ## Index computation for shocked cells
    valid_ix = collect(1:Nchi)
    sc_idy = shock["cell_indices"][valid_ix, 1]

    lin_idx_0 = CartesianIndex.(valid_ix, sc_idy)
    lin_idx_1 = CartesianIndex.(valid_ix, sc_idy .- 1)
    lin_idx_2 = CartesianIndex.(valid_ix, sc_idy .- 2)

    x_0 = x_mesh[lin_idx_0]
    y_0 = y_mesh[lin_idx_0]
    x_1 = x_mesh[lin_idx_1]
    y_1 = y_mesh[lin_idx_1]
    x_2 = x_mesh[lin_idx_2]
    y_2 = y_mesh[lin_idx_2]
    x_s = shock["points_x"][valid_ix]
    y_s = shock["points_y"][valid_ix]

    ## Distance computation for interpolation
    dist01 = @. sqrt((x_0 - x_1)^2 + (y_0 - y_1)^2)
    dist12 = @. sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2)
    dist1s = @. ((x_s - x_1)*(x_1 - x_2) + (y_s - y_1)*(y_1 - y_2)) / dist12
    dist0s = @. ((x_s - x_0)*(x_0 - x_1) + (y_s - y_0)*(y_0 - y_1)) / dist01

    ## Flow field gradient computation for linear interpolation
    p_1 = p_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
    p_2 = p_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
    slope_p = @. (p_1 - p_2) / dist12

    rho_1 = rho_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
    rho_2 = rho_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
    slope_rho = @. (rho_1 - rho_2) / dist12

    rho_u_1 = rho_u_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
    rho_u_2 = rho_u_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
    slope_rho_u = @. (rho_u_1 - rho_u_2) / dist12

    rho_v_1 = rho_v_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
    rho_v_2 = rho_v_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
    slope_rho_v = @. (rho_v_1 - rho_v_2) / dist12

    rho_E_1 = rho_E_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
    rho_E_2 = rho_E_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
    slope_rho_E = @. (rho_E_1 - rho_E_2) / dist12

    a_1 = a_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
    a_2 = a_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
    slope_a = @. (a_1 - a_2) / dist12

    ## Internal energy gradient (use fused broadcast to avoid temporaries)
    rho_velocity_squared = @. (rho_u_arr^2 + rho_v_arr^2) / rho_arr
    e = @. (rho_E_arr - rho_velocity_squared / 2) / rho_arr
    e_1 = e[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
    e_2 = e[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
    slope_e = @. (e_1 - e_2) / dist12

    ## Non-equilibrium chemistry gradients (gamma_star, cv_star)
    slope_gamma_star = nothing
    slope_cv_star = nothing
    if !chem_eq && is_chem
        gamma_star_arr = var["gamma_star"]::Matrix{Float64}
        cv_star_arr    = var["cv_star"]::Matrix{Float64}
        gamma_star_1 = gamma_star_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
        gamma_star_2 = gamma_star_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
        slope_gamma_star = @. (gamma_star_1 - gamma_star_2) / dist12
        cv_star_1 = cv_star_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
        cv_star_2 = cv_star_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
        slope_cv_star = @. (cv_star_1 - cv_star_2) / dist12
    end

    ## Initial guess from flow field interpolation
    p_s = @. slope_p * dist1s + p_1
    e_s = @. slope_e * dist1s + e_1
    rho_s = @. slope_rho * dist1s + rho_1
    a_s = @. slope_a * dist1s + a_1
    rho_u_s = @. slope_rho_u * dist1s + rho_u_1
    rho_v_s = @. slope_rho_v * dist1s + rho_v_1
    u_s = rho_u_s ./ rho_s
    v_s = rho_v_s ./ rho_s
    beta_arr = shock["beta"]
    ang_inf = @views beta_arr[valid_ix, 1] .- atan.(v_inf, u_inf)
    ang_s = @views beta_arr[valid_ix, 1] .- atan.(v_s, u_s)

    ## Thermodynamic property evaluation (gamma, cv)
    cv_s = nothing
    gamma_s = nothing
    gamma_inf = nothing
    if is_chem
        if chem_eq
            gamma_s = chemistry["eval_gamma_star"](rho_s * rho_factor, e_s * energy_factor)
            cv_s = chemistry["eval_cv_star"](rho_s * rho_factor, e_s * energy_factor)
        else
            gamma_s = chemistry["frozen"]["eval_gamma_star"](rho_s * rho_factor, e_s * energy_factor)
            cv_s = chemistry["frozen"]["eval_cv_star"](rho_s * rho_factor, e_s * energy_factor)
        end
        gamma_inf = fs["gamma_star"]
    else
        cv_s = 1
        gamma_s = fs["gamma"]
        gamma_inf = fs["gamma"]
    end

    ## Riemann invariant along C- characteristic
    @views normal_velocity_inf = sin.(ang_inf[valid_ix, 1]) .* sqrt.(u_inf.^2 .+ v_inf.^2)
    @views normal_velocity_s = sin.(ang_s[valid_ix, 1]) .* sqrt.(u_s.^2 .+ v_s.^2)
    J_minus = @. p_s - rho_s * a_s * normal_velocity_s

    ## Newton-Raphson iteration for post-shock pressure
    max_iter = 1000
    tol = 1e-12
    epsilon = 1e-4 * minimum(abs.(p_s))

    n_points = length(p_s)
    converged = falses(n_points)

    iter = 0
    for it in 1:max_iter
        iter = it
        # Evaluate residual at current pressure
        residual = COMPUTE_RESIDUAL_RIEMANN(p_infty, normal_velocity_inf,
            rho_inf, gamma_inf, p_s, a_s, gamma_s, J_minus)

        # Check convergence
        not_converged = .!converged .& (abs.(residual) .> tol)
        if !any(not_converged)
            break
        end

        # Numerical Jacobian via forward difference
        p_s_perturbed = copy(p_s)
        p_s_perturbed[not_converged] .= p_s[not_converged] .+ epsilon

        residual_perturbed = COMPUTE_RESIDUAL_RIEMANN(p_infty, normal_velocity_inf,
            rho_inf, gamma_inf, p_s_perturbed, a_s, gamma_s, J_minus)

        dR_dP = (residual_perturbed .- residual) / epsilon

        # Newton update
        delta_p = .-residual ./ dR_dP
        p_s[not_converged] .= p_s[not_converged] .+ delta_p[not_converged]

        # Update convergence flags
        converged .= converged .| (abs.(residual) .<= tol)
    end

    if iter == max_iter && any(.!converged)
        @warn "$(sum(.!converged)) points did not converge after $max_iter iterations"
    end

    ## Post-shock state from converged pressure
    y = p_s ./ p_infty

    # Density ratio from Rankine-Hugoniot
    A = @. (gamma_s + 1) / (gamma_s - 1)
    B = @. (gamma_inf + 1) / (gamma_inf - 1)
    rho_s = @. rho_inf * (y * A + 1) / (y + B)

    # Internal energy
    e_s = @. p_s / ((gamma_s - 1) * rho_s)

    ## Shock velocity computation
    @views normal_velocity_inf = sin.(ang_inf[valid_ix, 1]) .* sqrt.(u_inf.^2 .+ v_inf.^2)
    temp1 = @. 2 * y / (gamma_s - 1) - 2 / (gamma_inf - 1)
    temp2 = @. y * A + 1
    temp3 = @. p_infty * (y - 1) / rho_inf
    shock["relative_increase_velocity"] = @. -normal_velocity_inf + sqrt(temp3 * temp2 / temp1)

    ## Lab-frame velocity decomposition
    @views tangential_velocity_s = cos.(ang_inf[valid_ix, 1]) .* sqrt.(u_inf.^2 .+ v_inf.^2)
    rel_inc_vel = shock["relative_increase_velocity"]
    normal_velocity_s = @. -rel_inc_vel + rho_inf / rho_s * (rel_inc_vel + normal_velocity_inf)

    ## Shock speed in grid coordinates
    formulation = shock["formulation"]::String
    if formulation == "Lagrangian"
        shock["speed_x"] = @. -rel_inc_vel * sin(beta_arr)
        shock["speed_y"] = @. rel_inc_vel * cos(beta_arr)
    elseif formulation == "Eulerian"
        shock["speed_x"] = @. -rel_inc_vel * sin(beta_arr)
        shock["speed_y"] = @. rel_inc_vel * cos(beta_arr)
        shock["speed_x"] = @. shock["speed_x"] - shock["speed_y"] * tan(pi/2 - beta_arr)
        shock["speed_y"] = zeros(size(rel_inc_vel))
    end

    ## Conservative variable reconstruction at the shock
    @views rho_u_s = @. rho_s * (normal_velocity_s * sin(beta_arr[valid_ix, 1]) +
        tangential_velocity_s * cos(beta_arr[valid_ix, 1]))
    @views rho_v_s = @. rho_s * (-normal_velocity_s * cos(beta_arr[valid_ix, 1]) +
        tangential_velocity_s * sin(beta_arr[valid_ix, 1]))
    rho_E_s = @. e_s * rho_s + (rho_u_s^2 + rho_v_s^2) / (2 * rho_s)

    ## Store shock properties
    props = shock["properties"]::Dict{String, Any}
    props["rho"] = rho_s
    props["rho_u"] = rho_u_s
    props["rho_v"] = rho_v_s
    props["rho_E"] = rho_E_s
    props["p"] = p_s
    props["gamma_star"] = gamma_s
    props["cv_star"] = cv_s

    ## Update shocked cells via linear interpolation
    idx_0 = CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1)

    rho_arr[idx_0] .= @. -slope_rho * dist0s + rho_s
    rho_u_arr[idx_0] .= @. -slope_rho_u * dist0s + rho_u_s
    rho_v_arr[idx_0] .= @. -slope_rho_v * dist0s + rho_v_s
    rho_E_arr[idx_0] .= @. -slope_rho_E * dist0s + rho_E_s
    p_arr[idx_0] .= @. -slope_p * dist0s + p_s

    if !chem_eq && is_chem
        gamma_star_arr = var["gamma_star"]::Matrix{Float64}
        cv_star_arr    = var["cv_star"]::Matrix{Float64}
        gamma_star_arr[idx_0] .= @. -slope_gamma_star * dist0s + gamma_s
        cv_star_arr[idx_0] .= @. -slope_cv_star * dist0s + cv_s
    end

    ## Extrapolate to ghost cells
    s = EXTRAPOLATE_CELLS_SHOCK(s)

    return s
end

# ========================================================================
#  Riemann Invariant Residual
# ========================================================================
function COMPUTE_RESIDUAL_RIEMANN(p_infty, normal_velocity_inf,
        rho_inf, gamma_inf, p_s, a_s, gamma_s, J_minus)
# COMPUTE_RESIDUAL_RIEMANN  Evaluate the Riemann-invariant residual for
#   the Newton-Raphson pressure iteration.
#
#   residual = P_s - rho_s * a_s * U_n,s - J^-

    # Density ratio from Rankine-Hugoniot
    y = p_s ./ p_infty
    A = @. (gamma_s + 1) / (gamma_s - 1)
    B = @. (gamma_inf + 1) / (gamma_inf - 1)
    rho_s = @. rho_inf * (y * A + 1) / (y + B)
    temp1 = @. 2 * y / (gamma_s - 1) - 2 / (gamma_inf - 1)
    temp2 = @. y * A + 1
    temp3 = @. p_infty * (y - 1) / rho_inf
    relative_increase_velocity = @. -normal_velocity_inf + sqrt(temp3 * temp2 / temp1)

    # Post-shock normal velocity in lab frame
    normal_velocity_s = @. -relative_increase_velocity + rho_inf / rho_s * (relative_increase_velocity + normal_velocity_inf)

    # Residual: R = P_s - rho_s * a_s * U_n,s - J^-
    residual = @. p_s - rho_s * a_s * normal_velocity_s - J_minus

    return residual
end

# ========================================================================
#  Pressure-Based Rankine-Hugoniot Solver
# ========================================================================
function UPDATE_SHOCK_JUMP_PROPERTIES_PRESSURE_3(s::Dict{String, Any}, chemistry::Dict{String, Any})
# UPDATE_SHOCK_JUMP_PROPERTIES_PRESSURE_3  Solve R-H based on direct
#   pressure-jump interpolation (alternative to Riemann invariant method).

    ## Extract dict refs
    var   = s["var"]::Dict{String, Any}
    mesh  = s["mesh"]::Dict{String, Any}
    shock = s["shock"]::Dict{String, Any}
    chem  = s["chemistry"]::Dict{String, Any}
    fs    = s["freestream"]::Dict{String, Any}
    Nchi  = mesh["Nchi"]::Int

    rho_arr       = var["rho"]::Matrix{Float64}
    rho_u_arr     = var["rho_u"]::Matrix{Float64}
    rho_v_arr     = var["rho_v"]::Matrix{Float64}
    rho_E_arr     = var["rho_E"]::Matrix{Float64}
    p_arr         = var["p"]::Matrix{Float64}
    x_mesh        = mesh["x"]::Matrix{Float64}
    y_mesh        = mesh["y"]::Matrix{Float64}

    chem_eq = chem["chemical_equilibrium"]::Bool
    is_chem = chem["is_chemistry_enabled"]::Bool
    rho_factor    = fs["rho_factor"]::Float64
    energy_factor = fs["energy_factor"]::Float64

    ## Upstream conditions
    p_infty, u_inf, v_inf, rho_inf, e_inf = COMPUTE_UPSTREAM_CONDITIONS(s)
    s = UPDATE_FIELD_UPSTREAM(s)

    ## Index computation for shocked cells
    valid_ix = collect(1:Nchi)
    sc_idy = shock["cell_indices"][valid_ix, 1]

    lin_idx_0 = CartesianIndex.(valid_ix, sc_idy)
    lin_idx_1 = CartesianIndex.(valid_ix, sc_idy .- 1)
    lin_idx_2 = CartesianIndex.(valid_ix, sc_idy .- 2)

    x_0 = x_mesh[lin_idx_0]
    y_0 = y_mesh[lin_idx_0]
    x_1 = x_mesh[lin_idx_1]
    y_1 = y_mesh[lin_idx_1]
    x_2 = x_mesh[lin_idx_2]
    y_2 = y_mesh[lin_idx_2]
    x_s = shock["points_x"][valid_ix]
    y_s = shock["points_y"][valid_ix]

    ## Distance computation for interpolation
    dist01 = @. sqrt((x_0 - x_1)^2 + (y_0 - y_1)^2)
    dist12 = @. sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2)
    dist1s = @. ((x_s - x_1)*(x_1 - x_2) + (y_s - y_1)*(y_1 - y_2)) / dist12
    dist0s = @. ((x_s - x_0)*(x_0 - x_1) + (y_s - y_0)*(y_0 - y_1)) / dist01

    ## Flow field gradient computation
    p_0 = p_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 0)]
    p_1 = p_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
    p_2 = p_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
    slope_p = @. (p_1 - p_2) / dist12

    rho_1 = rho_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
    rho_2 = rho_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
    slope_rho = @. (rho_1 - rho_2) / dist12

    rho_u_1 = rho_u_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
    rho_u_2 = rho_u_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
    slope_rho_u = @. (rho_u_1 - rho_u_2) / dist12

    rho_v_1 = rho_v_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
    rho_v_2 = rho_v_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
    slope_rho_v = @. (rho_v_1 - rho_v_2) / dist12

    rho_E_1 = rho_E_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
    rho_E_2 = rho_E_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
    slope_rho_E = @. (rho_E_1 - rho_E_2) / dist12

    ## Internal energy gradient (use fused broadcast to avoid temporaries)
    rho_velocity_squared = @. (rho_u_arr^2 + rho_v_arr^2) / rho_arr
    e = @. (rho_E_arr - rho_velocity_squared / 2) / rho_arr
    e_1 = e[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
    e_2 = e[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
    slope_e = @. (e_1 - e_2) / dist12

    ## Non-equilibrium chemistry gradients
    slope_gamma_star = nothing
    slope_cv_star = nothing
    if !chem_eq && is_chem
        gamma_star_arr = var["gamma_star"]::Matrix{Float64}
        cv_star_arr    = var["cv_star"]::Matrix{Float64}
        gamma_star_1 = gamma_star_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
        gamma_star_2 = gamma_star_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
        slope_gamma_star = @. (gamma_star_1 - gamma_star_2) / dist12
        cv_star_1 = cv_star_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
        cv_star_2 = cv_star_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
        slope_cv_star = @. (cv_star_1 - cv_star_2) / dist12
    end

    ## Interpolated state at shock location (initial values for gamma evaluation)
    e_s = @. slope_e * dist1s + e_1
    rho_s = @. slope_rho * dist1s + rho_1

    ## Thermodynamic property evaluation (gamma, cv)
    cv_s = nothing
    gamma_s = nothing
    gamma_inf = nothing
    if is_chem
        if chem_eq
            gamma_s = chemistry["eval_gamma_star"](rho_s * rho_factor, e_s * energy_factor)
            cv_s = chemistry["eval_cv_star"](rho_s * rho_factor, e_s * energy_factor)
        else
            gamma_s = chemistry["frozen"]["eval_gamma_star"](rho_s * rho_factor, e_s * energy_factor)
            cv_s = chemistry["frozen"]["eval_cv_star"](rho_s * rho_factor, e_s * energy_factor)
        end
        gamma_inf = fs["gamma_star"]
    else
        gamma_s = fs["gamma"]
        gamma_inf = fs["gamma"]
        cv_s = 1
    end

    ## Pressure interpolation at shock
    p_s = @. slope_p * dist1s + p_1

    ## Pressure ratio and validity check
    y = p_s ./ p_infty

    if any(y .< 1)
        @warn "ERROR: Mach number less than 1"
        @warn "Add smoothing to shock and increase viscosity_scaling"
    end

    ## Density and energy from Rankine-Hugoniot
    A = @. (gamma_s + 1) / (gamma_s - 1)
    B = @. (gamma_inf + 1) / (gamma_inf - 1)
    rho_s = @. rho_inf * (y * A + 1) / (y + B)
    e_s = @. p_s / ((gamma_s - 1) * rho_s)

    ## Shock velocity computation
    beta_arr = shock["beta"]
    @views ang_inf = beta_arr[valid_ix, 1] .- atan.(v_inf, u_inf)
    @views normal_velocity_inf = sin.(ang_inf[valid_ix, 1]) .* sqrt.(u_inf.^2 .+ v_inf.^2)
    temp1 = @. 2 * y / (gamma_s - 1) - 2 / (gamma_inf - 1)
    temp2 = @. y * A + 1
    temp3 = @. p_infty * (y - 1) / rho_inf
    shock["relative_increase_velocity"] = @. -normal_velocity_inf + sqrt(temp3 * temp2 / temp1)

    ## Lab-frame velocity decomposition
    @views tangential_velocity_s = cos.(ang_inf[valid_ix, 1]) .* sqrt.(u_inf.^2 .+ v_inf.^2)
    rel_inc_vel = shock["relative_increase_velocity"]
    normal_velocity_s = @. -rel_inc_vel + rho_inf / rho_s * (rel_inc_vel + normal_velocity_inf)

    ## Shock speed in grid coordinates
    formulation = shock["formulation"]::String
    if formulation == "Lagrangian"
        shock["speed_x"] = @. -rel_inc_vel * sin(beta_arr)
        shock["speed_y"] = @. rel_inc_vel * cos(beta_arr)
    elseif formulation == "Eulerian"
        shock["speed_x"] = @. -rel_inc_vel * sin(beta_arr)
        shock["speed_y"] = @. rel_inc_vel * cos(beta_arr)
        shock["speed_x"] = @. shock["speed_x"] - shock["speed_y"] * tan(pi/2 - beta_arr)
        shock["speed_y"] = zeros(size(rel_inc_vel))
    end

    ## Conservative variable reconstruction at the shock
    @views rho_u_s = @. rho_s * (normal_velocity_s * sin(beta_arr[valid_ix, 1]) +
        tangential_velocity_s * cos(beta_arr[valid_ix, 1]))
    @views rho_v_s = @. rho_s * (-normal_velocity_s * cos(beta_arr[valid_ix, 1]) +
        tangential_velocity_s * sin(beta_arr[valid_ix, 1]))
    rho_E_s = @. e_s * rho_s + (rho_u_s^2 + rho_v_s^2) / (2 * rho_s)

    ## Store shock properties
    props = shock["properties"]::Dict{String, Any}
    props["rho"] = rho_s
    props["rho_u"] = rho_u_s
    props["rho_v"] = rho_v_s
    props["rho_E"] = rho_E_s
    props["p"] = p_s
    props["gamma_star"] = gamma_s
    props["cv_star"] = cv_s

    ## Update shocked cells via linear interpolation
    idx_0 = CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1)

    rho_arr[idx_0] .= @. -slope_rho * dist0s + rho_s
    rho_u_arr[idx_0] .= @. -slope_rho_u * dist0s + rho_u_s
    rho_v_arr[idx_0] .= @. -slope_rho_v * dist0s + rho_v_s
    rho_E_arr[idx_0] .= @. -slope_rho_E * dist0s + rho_E_s
    p_arr[idx_0] .= @. -slope_p * dist0s + p_s

    if !chem_eq && is_chem
        gamma_star_arr = var["gamma_star"]::Matrix{Float64}
        cv_star_arr    = var["cv_star"]::Matrix{Float64}
        gamma_star_arr[idx_0] .= @. -slope_gamma_star * dist0s + gamma_s
        cv_star_arr[idx_0] .= @. -slope_cv_star * dist0s + cv_s
    end

    ## Extrapolate to ghost cells
    s = EXTRAPOLATE_CELLS_SHOCK(s)

    return s
end

# ========================================================================
#  Mach-Number-Based Pressure Solver (Legacy)
# ========================================================================
function UPDATE_SHOCK_JUMP_PROPERTIES_PRESSURE(s::Dict{String, Any})
# UPDATE_SHOCK_JUMP_PROPERTIES_PRESSURE  Legacy solver using Mach-number-
#   based pressure jump for shock fitting with or without shock feedback.

    ## Extract dict refs
    var   = s["var"]::Dict{String, Any}
    mesh  = s["mesh"]::Dict{String, Any}
    shock = s["shock"]::Dict{String, Any}
    fs    = s["freestream"]::Dict{String, Any}
    Nchi  = mesh["Nchi"]::Int

    rho_arr       = var["rho"]::Matrix{Float64}
    rho_u_arr     = var["rho_u"]::Matrix{Float64}
    rho_v_arr     = var["rho_v"]::Matrix{Float64}
    rho_E_arr     = var["rho_E"]::Matrix{Float64}
    p_arr         = var["p"]::Matrix{Float64}
    x_mesh        = mesh["x"]::Matrix{Float64}
    y_mesh        = mesh["y"]::Matrix{Float64}

    if shock["feedback"]
        ## Precompute factors
        gamma = fs["gamma"]
        factor_M = (gamma + 1) / (2 * gamma)
        factor_add = (gamma - 1) / (2 * gamma)
        shock_diff_factor = Nchi^2 / s["Re_shock"]

        p_infty, u_inf, v_inf, rho_inf, e_inf = COMPUTE_UPSTREAM_CONDITIONS(s)
        a_s_inf = @. sqrt(fs["gamma"] * p_infty / rho_inf)
        u_mag = @. sqrt(v_inf^2 + u_inf^2)
        Mach_inf = u_mag ./ a_s_inf

        ## Index computation for shocked cells
        valid_ix = collect(1:Nchi)
        sc_idy = shock["cell_indices"][valid_ix, 1]

        lin_idx_0 = CartesianIndex.(valid_ix, sc_idy)
        lin_idx_1 = CartesianIndex.(valid_ix, sc_idy .- 1)
        lin_idx_2 = CartesianIndex.(valid_ix, sc_idy .- 2)

        x_0 = x_mesh[lin_idx_0]
        y_0 = y_mesh[lin_idx_0]
        x_1 = x_mesh[lin_idx_1]
        y_1 = y_mesh[lin_idx_1]
        x_2 = x_mesh[lin_idx_2]
        y_2 = y_mesh[lin_idx_2]
        x_s = shock["points_x"][valid_ix]
        y_s = shock["points_y"][valid_ix]

        ## Distance computation
        dist01 = @. sqrt((x_0 - x_1)^2 + (y_0 - y_1)^2)
        dist12 = @. sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2)
        dist1s = @. ((x_s - x_1)*(x_1 - x_2) + (y_s - y_1)*(y_1 - y_2)) / dist12
        dist0s = @. ((x_s - x_0)*(x_0 - x_1) + (y_s - y_0)*(y_0 - y_1)) / dist01

        ## Pressure interpolation
        p_1 = p_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
        p_2 = p_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
        slope_p = @. (p_1 - p_2) / dist12
        p_s = @. slope_p * dist1s + p_1

        ## Flow field gradients
        rho_1 = rho_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
        rho_2 = rho_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
        slope_rho = @. (rho_1 - rho_2) / dist12

        rho_u_1 = rho_u_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
        rho_u_2 = rho_u_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
        slope_rho_u = @. (rho_u_1 - rho_u_2) / dist12

        rho_v_1 = rho_v_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
        rho_v_2 = rho_v_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
        slope_rho_v = @. (rho_v_1 - rho_v_2) / dist12

        rho_E_1 = rho_E_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)]
        rho_E_2 = rho_E_arr[CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)]
        slope_rho_E = @. (rho_E_1 - rho_E_2) / dist12

        ## Shock Mach number and velocity
        beta_arr = shock["beta"]
        shock["M"] = @. real(sqrt(p_s / p_infty * factor_M + factor_add))
        if any(shock["M"] .< 1)
            @warn "ERROR: Mach number less than 1"
            @warn "Add smoothing to shock and increase viscosity_scaling"
        end

        @views ang = beta_arr[valid_ix, 1] .- atan.(v_inf, u_inf)
        shock["relative_increase_velocity"] = @. (shock["M"] - sin(ang) * Mach_inf) * a_s_inf

        ## Shock speed in grid coordinates
        rel_inc_vel = shock["relative_increase_velocity"]
        formulation = shock["formulation"]::String
        if formulation == "Lagrangian"
            shock["speed_x"] = @. -rel_inc_vel * sin(beta_arr)
            shock["speed_y"] = @. rel_inc_vel * cos(beta_arr)
        elseif formulation == "Eulerian"
            shock["speed_x"] = @. -rel_inc_vel * sin(beta_arr)
            shock["speed_y"] = @. rel_inc_vel * cos(beta_arr)
            shock["speed_x"] = @. shock["speed_x"] - shock["speed_y"] * tan(pi/2 - beta_arr)
            shock["speed_y"] = zeros(size(rel_inc_vel))
        end

        ## Conservative variable reconstruction (supersonic cells only)
        valid_ix = findall(shock["M"][:, 1] .> 1)

        if !isempty(valid_ix)
            M = shock["M"][valid_ix, 1]
            rho_s = @. (gamma + 1) * M^2 / ((gamma - 1) * M^2 + 2) * rho_inf

            @views tangential_velocity = cos.(ang[valid_ix, 1]) .* sqrt.(u_inf.^2 .+ v_inf.^2)
            normal_velocity = @. (M * a_s_inf) * (rho_inf / rho_s)

            shock_speed_x = shock["speed_x"][valid_ix, 1]
            shock_speed_y = shock["speed_y"][valid_ix, 1]
            @views rho_u_s = @. rho_s * (normal_velocity * sin(beta_arr[valid_ix, 1]) +
                tangential_velocity * cos(beta_arr[valid_ix, 1]) + shock_speed_x)
            @views rho_v_s = @. rho_s * (-normal_velocity * cos(beta_arr[valid_ix, 1]) +
                tangential_velocity * sin(beta_arr[valid_ix, 1]) + shock_speed_y)

            rho_E_s = @. p_s / (gamma - 1) + (rho_u_s^2 + rho_v_s^2) / (2 * rho_s)

            # Store properties
            props = shock["properties"]::Dict{String, Any}
            props["rho"] = rho_s
            props["rho_u"] = rho_u_s
            props["rho_v"] = rho_v_s
            props["rho_E"] = rho_E_s
            props["p"] = p_s
            props["gamma_star"] = gamma_s
            props["cv_star"] = cv_s

            # Update shocked cells via linear interpolation
            idx_0 = CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1)

            rho_arr[idx_0] .= @. -slope_rho * dist0s + rho_s
            rho_u_arr[idx_0] .= @. -slope_rho_u * dist0s + rho_u_s
            rho_v_arr[idx_0] .= @. -slope_rho_v * dist0s + rho_v_s
            rho_E_arr[idx_0] .= @. -slope_rho_E * dist0s + rho_E_s
            p_arr[idx_0] .= @. -slope_p * dist0s + p_s
        end
    else
        ## Passive shock (no feedback) -- simple extrapolation
        valid_ix = collect(1:Nchi)
        sc_idy = shock["cell_indices"][valid_ix, 1]

        idx_0 = CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1)
        idx_1 = CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 1)
        idx_2 = CartesianIndex.(valid_ix .+ 1, sc_idy .+ 1 .- 2)

        rho_E_arr[idx_0] .= 2 .* rho_E_arr[idx_1] .- rho_E_arr[idx_2]
    end

    ## Extrapolate to ghost cells
    s = EXTRAPOLATE_CELLS_SHOCK(s)

    return s
end
