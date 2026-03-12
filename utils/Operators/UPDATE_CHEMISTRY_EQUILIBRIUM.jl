function UPDATE_CHEMISTRY_EQUILIBRIUM(s::Dict{String, Any}, chemistry::Dict{String, Any})
# UPDATE_CHEMISTRY_EQUILIBRIUM - Evaluate equilibrium thermochemical and transport properties.
#
# Computes the equilibrium values of gamma_star, cv_star, dynamic
# viscosity (mu), and thermal conductivity (k) from the chemistry lookup
# tables. Derived quantities (Re_flow, Pr_flow) are also updated.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract dict refs
    var  = s["var"]::Dict{String, Any}
    chem = s["chemistry"]::Dict{String, Any}
    fs   = s["freestream"]::Dict{String, Any}

    ## Compute internal energy from conserved variables
    rho       = var["rho"]::Matrix{Float64}
    rho_E     = var["rho_E"]::Matrix{Float64}
    rho_u_arr = var["rho_u"]::Matrix{Float64}
    rho_v_arr = var["rho_v"]::Matrix{Float64}
    e = @. (rho_E - (rho_u_arr^2 + rho_v_arr^2) / rho / 2) / rho

    rho_factor    = fs["rho_factor"]::Float64
    energy_factor = fs["energy_factor"]::Float64
    is_chem       = chem["is_chemistry_enabled"]::Bool
    is_eq         = chem["chemical_equilibrium"]::Bool

    ## Pre-compute dimensional state once (avoids 4x redundant allocations)
    if is_chem
        rho_dim = @. rho_factor * rho
        e_dim   = @. energy_factor * e
    end

    ## Gamma_star (ratio of specific heats)
    if is_chem
        eval_gamma_star = chemistry["eval_gamma_star"]
        var["gamma_star_eq"] = eval_gamma_star(rho_dim, e_dim)
        if is_eq
            var["gamma_star"] = var["gamma_star_eq"]
        end
    else
        gs_val = fs["gamma_star"]::Float64
        gamma_star_arr = var["gamma_star"]::Matrix{Float64}
        fill!(gamma_star_arr, gs_val)
        var["gamma_star_eq"] = gamma_star_arr
    end

    ## cv_star (specific heat at constant volume)
    cv_fs = fs["cv"]::Float64
    if is_chem
        eval_cv_star = chemistry["eval_cv_star"]
        cv_star_dim = eval_cv_star(rho_dim, e_dim)
        var["cv_star_eq"] = @. cv_star_dim / cv_fs
        if is_eq
            var["cv_star"] = var["cv_star_eq"]
        end
    else
        cv_star_arr = var["cv_star"]::Matrix{Float64}
        fill!(cv_star_arr, 1.0)
        var["cv_star_eq"] = cv_star_arr
    end

    ## Dynamic viscosity (mu) and local Reynolds number
    Re_fs = Float64(fs["Re"])
    mu_fs = fs["mu"]::Float64
    if is_chem
        eval_mu = chemistry["eval_mu"]
        var["mu_star"] = @. eval_mu(rho_dim, e_dim) / mu_fs
        mu_star_arr = var["mu_star"]::Matrix{Float64}
        var["Re_flow"] = @. Re_fs / mu_star_arr
    else
        mu_star_arr = var["mu_star"]::Matrix{Float64}
        fill!(mu_star_arr, 1.0)
        Re_flow_arr = var["Re_flow"]::Matrix{Float64}
        fill!(Re_flow_arr, Re_fs)
    end

    ## Thermal conductivity (k) and local Prandtl number
    k_fs = fs["k"]::Float64
    if is_chem
        eval_k = chemistry["eval_k"]
        var["k_star"] = @. eval_k(rho_dim, e_dim) / k_fs
        gamma_star_arr = var["gamma_star"]::Matrix{Float64}
        cv_star_arr    = var["cv_star"]::Matrix{Float64}
        mu_star_arr    = var["mu_star"]::Matrix{Float64}
        k_star_arr     = var["k_star"]::Matrix{Float64}
        var["Pr_flow"] = @. mu_star_arr * mu_fs * cv_star_arr * cv_fs * gamma_star_arr / k_star_arr / k_fs
    else
        Pr_fs = fs["Pr"]::Float64
        k_star_arr = var["k_star"]::Matrix{Float64}
        fill!(k_star_arr, 1.0)
        Pr_flow_arr = var["Pr_flow"]::Matrix{Float64}
        fill!(Pr_flow_arr, Pr_fs)
    end

    return s
end
