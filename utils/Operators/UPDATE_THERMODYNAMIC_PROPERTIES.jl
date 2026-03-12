function UPDATE_THERMODYNAMIC_PROPERTIES(s::Dict{String, Any}, chemistry::Dict{String, Any})
# UPDATE_THERMODYNAMIC_PROPERTIES - Compute pressure and temperature from conserved variables.
#
# Evaluates pressure using p = rho * e * (gamma* - 1) and temperature
# using T = e / cv_star.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract dict refs
    var = s["var"]::Dict{String, Any}
    bc  = s["boundary_conditions"]::Dict{String, Any}

    ## Compute internal energy, pressure, and temperature
    rho         = var["rho"]::Matrix{Float64}
    rho_E       = var["rho_E"]::Matrix{Float64}
    rho_u_arr   = var["rho_u"]::Matrix{Float64}
    rho_v_arr   = var["rho_v"]::Matrix{Float64}
    gamma_star  = var["gamma_star"]::Matrix{Float64}
    cv_star     = var["cv_star"]::Matrix{Float64}

    ## Compute internal energy, pressure, and temperature via fused broadcasting
    e = @. (rho_E - (rho_u_arr^2 + rho_v_arr^2) / rho / 2) / rho
    var["p"] = @. rho * e * (gamma_star - 1)               # Non-dimensional rho * U^2
    var["T"] = @. e / cv_star                                # Non-dimensional U^2 / cv_infty

    ## Wall boundary condition at eta0, keep pressure smooth, just update temperature for isothermal walls
    bc_eta0 = bc["boundary_eta0"]["name"]
    if bc_eta0 == "no_slip_isothermal"
        p = var["p"]::Matrix{Float64}
        @views @. p[:, 1] = 2 * p[:, 2] - p[:, 3]  # Interpolate pressure
    end

    return s
end
