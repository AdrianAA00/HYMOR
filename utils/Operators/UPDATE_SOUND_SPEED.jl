function UPDATE_SOUND_SPEED(s::Dict{String, Any}, chemistry::Dict{String, Any})
# UPDATE_SOUND_SPEED - Compute the local speed of sound from the current flow state.
#
# Evaluates the speed of sound either from chemistry lookup tables
# or from the ideal-gas relation a = sqrt(gamma* * p / rho).
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract dict refs
    var  = s["var"]::Dict{String, Any}
    chem = s["chemistry"]::Dict{String, Any}
    fs   = s["freestream"]::Dict{String, Any}

    ## Compute internal energy
    rho       = var["rho"]::Matrix{Float64}
    rho_E     = var["rho_E"]::Matrix{Float64}
    rho_u_arr = var["rho_u"]::Matrix{Float64}
    rho_v_arr = var["rho_v"]::Matrix{Float64}
    e = @. (rho_E - (rho_u_arr^2 + rho_v_arr^2) / rho / 2) / rho

    ## Update sound speed
    rho_factor    = fs["rho_factor"]::Float64
    energy_factor = fs["energy_factor"]::Float64
    vel_factor    = fs["velocity_factor"]::Float64

    if chem["is_chemistry_enabled"]::Bool
        eval_a = chemistry["eval_a"]
        rho_dim = @. rho_factor * rho
        e_dim   = @. energy_factor * e
        var["a"] = @. eval_a(rho_dim, e_dim) / vel_factor
    else
        p          = var["p"]::Matrix{Float64}
        gamma_star = var["gamma_star"]::Matrix{Float64}
        var["a"] = @. sqrt(gamma_star * p / rho) / vel_factor
    end

    return s
end
