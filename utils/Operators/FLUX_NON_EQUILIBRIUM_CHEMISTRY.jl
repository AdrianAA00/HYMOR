function FLUX_NON_EQUILIBRIUM_CHEMISTRY(s::Dict{String, Any}, chemistry::Dict{String, Any})
# FLUX_NON_EQUILIBRIUM_CHEMISTRY - Compute fluxes for non-equilibrium chemistry transport.
#
# Handles advection, relaxation towards equilibrium (analytic exponential
# integrator), and optional viscous diffusion of gamma_star and cv_star.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract dict refs to avoid repeated lookups
    chem        = s["chemistry"]::Dict{String, Any}
    var         = s["var"]::Dict{String, Any}
    flux        = s["flux"]::Dict{String, Any}
    mesh        = s["mesh"]::Dict{String, Any}

    if chem["is_chemistry_enabled"]::Bool
        if chem["chemical_equilibrium"]::Bool
            ## Chemical equilibrium mode
            var["gamma_star"] = var["gamma_star_eq"]
            var["cv_star"]    = var["cv_star_eq"]
            flux_gs = flux["gamma_star"]::Matrix{Float64}
            flux_cv = flux["cv_star"]::Matrix{Float64}
            fill!(flux_gs, 0.0)
            fill!(flux_cv, 0.0)
        else
            ## Advection equations
            g  = var["gamma_star"]::Matrix{Float64}
            cv = var["cv_star"]::Matrix{Float64}
            dgamma_dx, dgamma_dy = DERIVATIVE_EXT(g, s)
            dcv_dx, dcv_dy       = DERIVATIVE_EXT(cv, s)

            rho_arr   = var["rho"]::Matrix{Float64}
            rho_u_arr = var["rho_u"]::Matrix{Float64}
            rho_v_arr = var["rho_v"]::Matrix{Float64}

            @views u = rho_u_arr[2:end-1, 2:end-1] ./ rho_arr[2:end-1, 2:end-1]
            @views v = rho_v_arr[2:end-1, 2:end-1] ./ rho_arr[2:end-1, 2:end-1]
            flux["gamma_star"] = @. -u * dgamma_dx - v * dgamma_dy
            flux["cv_star"]    = @. -u * dcv_dx    - v * dcv_dy

            ## Relaxation towards chemical equilibrium
            s = UPDATE_TAU_CHEMISTRY(s, chemistry)

            # Analytic exponential integration to avoid stiffness
            dt               = s["time_integration"]["dt"]::Float64
            gamma_star_eq    = var["gamma_star_eq"]::Matrix{Float64}
            cv_star_eq       = var["cv_star_eq"]::Matrix{Float64}
            tau_gamma        = s["tau_gamma_star"]::Matrix{Float64}
            tau_cv           = s["tau_cv_star"]::Matrix{Float64}
            flux_gamma       = flux["gamma_star"]::Matrix{Float64}
            flux_cv          = flux["cv_star"]::Matrix{Float64}

            @views @. flux_gamma = flux_gamma +
                (gamma_star_eq[2:end-1, 2:end-1] - g[2:end-1, 2:end-1]) / dt *
                (1 - exp(-dt / tau_gamma[2:end-1, 2:end-1]))
            @views @. flux_cv = flux_cv +
                (cv_star_eq[2:end-1, 2:end-1] - cv[2:end-1, 2:end-1]) / dt *
                (1 - exp(-dt / tau_cv[2:end-1, 2:end-1]))

            ## Viscous diffusion term (species diffusion analogue)
            x_Ext = mesh["x_Ext"]::Matrix{Float64}
            y_Ext = mesh["y_Ext"]::Matrix{Float64}

            @views centroids_bt_distance = @. sqrt((y_Ext[2:end-1, 2:end] - y_Ext[2:end-1, 1:end-1])^2 +
                                                   (x_Ext[2:end-1, 2:end] - x_Ext[2:end-1, 1:end-1])^2)
            @views centroids_lr_distance = @. sqrt((y_Ext[2:end, 2:end-1] - y_Ext[1:end-1, 2:end-1])^2 +
                                                   (x_Ext[2:end, 2:end-1] - x_Ext[1:end-1, 2:end-1])^2)

            @views dissip_gamma = @. (g[1:end-2, 2:end-1] - 2 * g[2:end-1, 2:end-1] + g[3:end, 2:end-1]) / centroids_lr_distance[1:end-1, :]^2 +
                                     (g[2:end-1, 1:end-2] - 2 * g[2:end-1, 2:end-1] + g[2:end-1, 3:end]) / centroids_bt_distance[:, 1:end-1]^2
            @views dissip_cv = @. (cv[1:end-2, 2:end-1] - 2 * cv[2:end-1, 2:end-1] + cv[3:end, 2:end-1]) / centroids_lr_distance[1:end-1, :]^2 +
                                   (cv[2:end-1, 1:end-2] - 2 * cv[2:end-1, 2:end-1] + cv[2:end-1, 3:end]) / centroids_bt_distance[:, 1:end-1]^2

            Re_inv = 1.0 / (s["freestream"]["Re"]::Float64)
            @. flux_gamma = flux_gamma + Re_inv * dissip_gamma / 2
            @. flux_cv    = flux_cv    + Re_inv * dissip_cv / 2
        end
    end

    return s
end

function UPDATE_TAU_CHEMISTRY(s::Dict{String, Any}, chemistry::Dict{String, Any})
# UPDATE_TAU_CHEMISTRY - Update relaxation times for non-equilibrium chemistry.
#
# Computes tau_gamma_star and tau_cv_star from Arrhenius-type fits.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract dict refs
    chem = s["chemistry"]::Dict{String, Any}
    var  = s["var"]::Dict{String, Any}
    fs   = s["freestream"]::Dict{String, Any}
    neq  = chemistry["neq"]::Dict{String, Any}

    if chem["chemical_equilibrium"]
        ## Chemical equilibrium: infinite relaxation times
        s["tau_gamma_star"] = Inf
        s["tau_cv_star"] = Inf
    else
        ## Extract dimensional temperature and pressure
        T = var["T"]::Matrix{Float64}
        vel_factor = fs["velocity_factor"]::Float64
        rho_factor = fs["rho_factor"]::Float64
        p_arr = var["p"]::Matrix{Float64}
        vel2_rho = vel_factor^2 * rho_factor
        p = @. p_arr * vel2_rho

        # Determine model type (default to linear for backward compatibility)
        if haskey(chem, "non_equilibrium_model")
            model_type = chem["non_equilibrium_model"]::String
        else
            model_type = "linear"
        end

        ## Compute relaxation times from Arrhenius fits
        m_gamma = neq["m_gamma"]::Float64
        m_cv    = neq["m_cv"]::Float64

        if lowercase(model_type) == "linear"
            # Linear model: ln(tau*p^m) = ln(T) + a0 + a1/T
            gamma_lin = neq["gamma"]["linear"]::Dict{String, Any}
            cv_lin    = neq["cv"]["linear"]::Dict{String, Any}
            g_a0 = gamma_lin["a0"]::Float64;  g_a1 = gamma_lin["a1"]::Float64
            c_a0 = cv_lin["a0"]::Float64;     c_a1 = cv_lin["a1"]::Float64
            s["tau_gamma_star"] = @. exp(g_a0 + g_a1 / T + log(T)) / p^m_gamma
            s["tau_cv_star"]    = @. exp(c_a0 + c_a1 / T + log(T)) / p^m_cv
        else
            # Quadratic model: ln(tau*p^m) = ln(T) + a0 + a1/T + a2/T^2
            gamma_quad = neq["gamma"]["quadratic"]::Dict{String, Any}
            cv_quad    = neq["cv"]["quadratic"]::Dict{String, Any}
            g_a0 = gamma_quad["a0"]::Float64;  g_a1 = gamma_quad["a1"]::Float64;  g_a2 = gamma_quad["a2"]::Float64
            c_a0 = cv_quad["a0"]::Float64;     c_a1 = cv_quad["a1"]::Float64;     c_a2 = cv_quad["a2"]::Float64
            s["tau_gamma_star"] = @. exp(g_a0 + g_a1 / T + g_a2 / T^2 + log(T)) / p^m_gamma
            s["tau_cv_star"]    = @. exp(c_a0 + c_a1 / T + c_a2 / T^2 + log(T)) / p^m_cv
        end

        ## Non-dimensionalise relaxation times
        U_over_L = fs["U"]::Float64 / fs["L"]::Float64
        tau_gs = s["tau_gamma_star"]::Matrix{Float64}
        tau_cv = s["tau_cv_star"]::Matrix{Float64}
        @. tau_gs *= U_over_L
        @. tau_cv *= U_over_L
    end

    return s
end
