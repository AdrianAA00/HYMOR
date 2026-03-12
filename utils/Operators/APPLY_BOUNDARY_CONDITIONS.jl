function APPLY_BOUNDARY_CONDITIONS(s::Dict{String, Any}, chemistry::Dict{String, Any})
# APPLY_BOUNDARY_CONDITIONS - Enforce boundary conditions on all domain edges
#
# Applies ghost-cell boundary conditions on the four edges of a 2D
# structured grid (eta0, eta1, chi0, chi1). Each edge can be independently
# assigned a boundary condition type.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract dict refs once at top to avoid repeated lookups
    var  = s["var"]::Dict{String, Any}
    bc   = s["boundary_conditions"]::Dict{String, Any}
    fs   = s["freestream"]::Dict{String, Any}
    mesh = s["mesh"]::Dict{String, Any}
    chem = s["chemistry"]::Dict{String, Any}

    rho         = var["rho"]::Matrix{Float64}
    rho_u_arr   = var["rho_u"]::Matrix{Float64}
    rho_v_arr   = var["rho_v"]::Matrix{Float64}
    rho_E_arr   = var["rho_E"]::Matrix{Float64}
    p_arr       = var["p"]::Matrix{Float64}
    T_arr       = var["T"]::Matrix{Float64}

    chem_eq = chem["chemical_equilibrium"]::Bool

    # Only extract these if chemistry is non-equilibrium
    gamma_star_arr = var["gamma_star"]::Matrix{Float64}
    cv_star_arr    = var["cv_star"]::Matrix{Float64}

    # Freestream values
    rho_0   = fs["rho_0"]
    rho_u_0 = fs["rho_u_0"]
    rho_v_0 = fs["rho_v_0"]
    rho_E_0 = haskey(fs, "rho_E_0") ? fs["rho_E_0"] : nothing
    p_0     = haskey(fs, "p_0") ? fs["p_0"] : nothing

    ##########################################################################
    ## Left boundary (chi0)
    bc_chi0 = bc["boundary_chi0"]["name"]

    if bc_chi0 == "inflow_subsonic"
        @views rho[1,:] .= rho_0
        @views rho_u_arr[1,:] .= rho_u_0
        @views rho_v_arr[1,:] .= rho_v_0
        @views @. rho_E_arr[1,:] = 2 * rho_E_arr[2,:] - rho_E_arr[3,:]
        @views @. p_arr[1,:] = 2 * p_arr[2,:] - p_arr[3,:]

        if !chem_eq
            gamma_star_eq = var["gamma_star_eq"]::Matrix{Float64}
            cv_star_eq    = var["cv_star_eq"]::Matrix{Float64}
            @views gamma_star_arr[1,:] .= gamma_star_eq[1,:]
            @views cv_star_arr[1,:] .= cv_star_eq[1,:]
        end

    elseif bc_chi0 == "inflow_supersonic"
        @views rho[1,:] .= rho_0
        @views rho_u_arr[1,:] .= rho_u_0
        @views rho_v_arr[1,:] .= rho_v_0
        @views rho_E_arr[1,:] .= rho_E_0
        @views p_arr[1,:] .= p_0

        if !chem_eq
            gamma_star_eq = var["gamma_star_eq"]::Matrix{Float64}
            cv_star_eq    = var["cv_star_eq"]::Matrix{Float64}
            @views gamma_star_arr[1,:] .= gamma_star_eq[1,:]
            @views cv_star_arr[1,:] .= cv_star_eq[1,:]
        end

    elseif bc_chi0 == "no_slip_adiabatic"
        @views rho[1,:] .= rho[2,:]
        @views rho_u_arr[1,:] .= .-rho_u_arr[2,:]
        @views rho_v_arr[1,:] .= .-rho_v_arr[2,:]
        @views @. rho_E_arr[1,:] = 2 * rho_E_arr[2,:] - rho_E_arr[3,:]
        @views @. p_arr[1,:] = 2 * p_arr[2,:] - p_arr[3,:]
        @views T_arr[1,:] .= T_arr[2,:]  # No temperature gradient for adiabatic wall

        if !chem_eq
            @views @. gamma_star_arr[1,:] = 2 * gamma_star_arr[2,:] - gamma_star_arr[3,:]
            @views @. cv_star_arr[1,:] = 2 * cv_star_arr[2,:] - cv_star_arr[3,:]
        end

    elseif bc_chi0 == "no_slip_isothermal"
        @views rho[1,:] .= rho[2,:]
        @views rho_u_arr[1,:] .= .-rho_u_arr[2,:]
        @views rho_v_arr[1,:] .= .-rho_v_arr[2,:]
        @views @. rho_E_arr[1,:] = 2 * rho_E_arr[2,:] - rho_E_arr[3,:]
        @views @. p_arr[1,:] = 2 * p_arr[2,:] - p_arr[3,:]
        T_w = bc["boundary_chi0"]["Tw"] * fs["cv"] / fs["energy_factor"]  # Non-dimensional wall temperature
        @views @. T_arr[1,:] = 2 * T_w - T_arr[2,:]  # Set wall temperature

        if !chem_eq
            @views gamma_star_arr[1,:] .= gamma_star_arr[2,:]
            @views cv_star_arr[1,:] .= cv_star_arr[2,:]
        end

    elseif bc_chi0 == "outflow_subsonic"
        @views @. rho[1,:] = 2 * rho[2,:] - rho[3,:]
        @views @. rho_E_arr[1,:] = 2 * rho_E_arr[2,:] - rho_E_arr[3,:]
        @views @. rho_u_arr[1,:] = 2 * rho_u_arr[2,:] - rho_u_arr[3,:]
        @views @. rho_v_arr[1,:] = 2 * rho_v_arr[2,:] - rho_v_arr[3,:]
        @views p_arr[1,:] .= p_arr[2,:]  # zero pressure gradient

        if !chem_eq
            @views gamma_star_arr[1,:] .= gamma_star_arr[2,:]
            @views cv_star_arr[1,:] .= cv_star_arr[2,:]
        end

    elseif bc_chi0 == "outflow_supersonic"
        @views @. rho[1,:] = 2 * rho[2,:] - rho[3,:]
        @views @. rho_E_arr[1,:] = 2 * rho_E_arr[2,:] - rho_E_arr[3,:]
        @views @. rho_u_arr[1,:] = 2 * rho_u_arr[2,:] - rho_u_arr[3,:]
        @views @. rho_v_arr[1,:] = 2 * rho_v_arr[2,:] - rho_v_arr[3,:]
        @views @. p_arr[1,:] = 2 * p_arr[2,:] - p_arr[3,:]

        if !chem_eq
            @views @. gamma_star_arr[1,:] = 2 * gamma_star_arr[2,:] - gamma_star_arr[3,:]
            @views @. cv_star_arr[1,:] = 2 * cv_star_arr[2,:] - cv_star_arr[3,:]
        end

    elseif bc_chi0 == "outflow_supersonic_1st"
        @views rho_u_arr[1,:] .= rho_u_arr[2,:]
        @views rho_v_arr[1,:] .= rho_v_arr[2,:]
        @views rho[1,:] .= rho[2,:]
        @views rho_E_arr[1,:] .= rho_E_arr[2,:]
        @views p_arr[1,:] .= p_arr[2,:]

        if !chem_eq
            @views gamma_star_arr[1,:] .= gamma_star_arr[2,:]
            @views cv_star_arr[1,:] .= cv_star_arr[2,:]
        end

    elseif bc_chi0 == "outflow_NRCBC"
        # LODI boundary conditions (Poinsot & Lele): L1 = 0
        # Extrapolate all variables (supersonic default)
        @views @. rho[1,:] = 2 * rho[2,:] - rho[3,:]
        @views @. rho_u_arr[1,:] = 2 * rho_u_arr[2,:] - rho_u_arr[3,:]
        @views @. rho_v_arr[1,:] = 2 * rho_v_arr[2,:] - rho_v_arr[3,:]
        @views @. rho_E_arr[1,:] = 2 * rho_E_arr[2,:] - rho_E_arr[3,:]
        @views @. p_arr[1,:] = 2 * p_arr[2,:] - p_arr[3,:]

        # Sound speed and density at first interior cell
        a_arr = var["a"]::Matrix{Float64}
        lr_x_normal = mesh["lr_x_normal"]::Matrix{Float64}
        lr_y_normal = mesh["lr_y_normal"]::Matrix{Float64}

        @views begin
            a = a_arr[2, 2:end-1]
            rho_int = rho[2, 2:end-1]

            # Normal velocity at ghost cell, positive exiting domain
            u_n_ghost = @. -(rho_u_arr[1, 2:end-1] * lr_x_normal[1,:] +
                          rho_v_arr[1, 2:end-1] * lr_y_normal[1,:]) / rho[1, 2:end-1]

            # Normal velocity at first interior cell
            u_n_int = @. -(rho_u_arr[2, 2:end-1] * lr_x_normal[1,:] +
                        rho_v_arr[2, 2:end-1] * lr_y_normal[1,:]) / rho_int

            # LODI correction for subsonic cells (L1 = 0 => dp/dx_n = rho*c*du_n/dx_n)
            subsonic = u_n_ghost .< a
            if any(subsonic)
                p_int = p_arr[2, 2:end-1]
                p_lodi = @. p_int + rho_int * a * (u_n_ghost - u_n_int)
                @. p_arr[1, 2:end-1] = subsonic * p_lodi + (!subsonic) * p_arr[1, 2:end-1]
                rho_e = @. p_arr[1, 2:end-1] / (gamma_star_arr[1, 2:end-1] - 1)
                rho_E_lodi = @. rho_e + 0.5 * (rho_u_arr[1, 2:end-1]^2 + rho_v_arr[1, 2:end-1]^2) / rho[1, 2:end-1]
                @. rho_E_arr[1, 2:end-1] = subsonic * rho_E_lodi + (!subsonic) * rho_E_arr[1, 2:end-1]
                s["p_int"] = p_int
                s["p_lodi"] = p_lodi
            end
        end

        if !chem_eq
            gamma_star_eq = var["gamma_star_eq"]::Matrix{Float64}
            cv_star_eq    = var["cv_star_eq"]::Matrix{Float64}
            @views gamma_star_arr[1,:] .= gamma_star_eq[1,:]
            @views cv_star_arr[1,:] .= cv_star_eq[1,:]
        end

    elseif bc_chi0 == "periodic"
        @views rho_u_arr[1,:] .= rho_u_arr[end-1,:]
        @views rho_v_arr[1,:] .= rho_v_arr[end-1,:]
        @views rho[1,:] .= rho[end-1,:]
        @views rho_E_arr[1,:] .= rho_E_arr[end-1,:]
        @views p_arr[1,:] .= p_arr[end-1,:]
        @views T_arr[1,:] .= T_arr[end-1,:]

        if !chem_eq
            @views gamma_star_arr[1,:] .= gamma_star_arr[end-1,:]
            @views cv_star_arr[1,:] .= cv_star_arr[end-1,:]
        end

    elseif bc_chi0 == "symmetry"
        @views rho[1,:] .= rho[2,:]
        @views p_arr[1,:] .= p_arr[2,:]

        lr_x_normal = mesh["lr_x_normal"]::Matrix{Float64}
        lr_y_normal = mesh["lr_y_normal"]::Matrix{Float64}

        # Remove wall-normal velocity component via reflection
        @views begin
            flux_vel = @. rho_u_arr[2, 2:end-1] * lr_x_normal[1,:] +
                       rho_v_arr[2, 2:end-1] * lr_y_normal[1,:]
            @. rho_u_arr[1, 2:end-1] = rho_u_arr[2, 2:end-1] - 2 * flux_vel * lr_x_normal[1,:]
            @. rho_v_arr[1, 2:end-1] = rho_v_arr[2, 2:end-1] - 2 * flux_vel * lr_y_normal[1,:]
        end
        rho_u_arr[1,1] = rho_u_arr[1,2]
        rho_v_arr[1,1] = rho_v_arr[1,2]
        rho_u_arr[1,end] = rho_u_arr[1,end-1]
        rho_v_arr[1,end] = rho_v_arr[1,end-1]

        # Remove heat flux (adiabatic symmetry)
        @views rho_E_arr[1,:] .= rho_E_arr[2,:]
        @views T_arr[1,:] .= T_arr[2,:]

        if !chem_eq
            @views gamma_star_arr[1,:] .= gamma_star_arr[2,:]
            @views cv_star_arr[1,:] .= cv_star_arr[2,:]
        end

    else
        error("Non defined Boundary-condition, chi0")
    end

    ##########################################################################
    ## Right boundary (chi1)
    bc_chi1 = bc["boundary_chi1"]["name"]

    if bc_chi1 == "inflow_subsonic"
        @views rho_u_arr[end,:] .= rho_u_0
        @views rho_v_arr[end,:] .= rho_v_0
        @views rho[end,:] .= rho_0
        @views @. rho_E_arr[end,:] = 2 * rho_E_arr[end-1,:] - rho_E_arr[end-2,:]
        @views @. p_arr[end,:] = 2 * p_arr[end-1,:] - p_arr[end-2,:]

        if !chem_eq
            gamma_star_eq = var["gamma_star_eq"]::Matrix{Float64}
            cv_star_eq    = var["cv_star_eq"]::Matrix{Float64}
            @views gamma_star_arr[end,:] .= gamma_star_eq[end,:]
            @views cv_star_arr[end,:] .= cv_star_eq[end,:]
        end

    elseif bc_chi1 == "inflow_supersonic"
        @views rho_u_arr[end,:] .= rho_u_0
        @views rho_v_arr[end,:] .= rho_v_0
        @views rho[end,:] .= rho_0
        @views rho_E_arr[end,:] .= rho_E_0
        @views p_arr[end,:] .= p_0

        if !chem_eq
            gamma_star_eq = var["gamma_star_eq"]::Matrix{Float64}
            cv_star_eq    = var["cv_star_eq"]::Matrix{Float64}
            @views gamma_star_arr[end,:] .= gamma_star_eq[end,:]
            @views cv_star_arr[end,:] .= cv_star_eq[end,:]
        end

    elseif bc_chi1 == "no_slip_adiabatic"
        @views rho[end,:] .= rho[end-1,:]
        @views rho_u_arr[end,:] .= .-rho_u_arr[end-1,:]
        @views rho_v_arr[end,:] .= .-rho_v_arr[end-1,:]
        @views @. rho_E_arr[end,:] = 2 * rho_E_arr[end-1,:] - rho_E_arr[end-2,:]
        @views @. p_arr[end,:] = 2 * p_arr[end-1,:] - p_arr[end-2,:]
        @views T_arr[end,:] .= T_arr[end-1,:]  # No temperature gradient for adiabatic wall

        if !chem_eq
            @views @. gamma_star_arr[end,:] = 2 * gamma_star_arr[end-1,:] - gamma_star_arr[end-2,:]
            @views @. cv_star_arr[end,:] = 2 * cv_star_arr[end-1,:] - cv_star_arr[end-2,:]
        end

    elseif bc_chi1 == "no_slip_isothermal"
        @views rho[end,:] .= rho[end-1,:]
        @views rho_u_arr[end,:] .= .-rho_u_arr[end-1,:]
        @views rho_v_arr[end,:] .= .-rho_v_arr[end-1,:]
        @views @. rho_E_arr[end,:] = 2 * rho_E_arr[end-1,:] - rho_E_arr[end-2,:]
        @views @. p_arr[end,:] = 2 * p_arr[end-1,:] - p_arr[end-2,:]
        T_w = bc["boundary_chi0"]["Tw"] * fs["cv"] / fs["energy_factor"]  # Non-dimensional wall temperature
        @views @. T_arr[end,:] = 2 * T_w - T_arr[end-1,:]  # Set wall temperature

        if !chem_eq
            @views gamma_star_arr[end,:] .= gamma_star_arr[end-1,:]
            @views cv_star_arr[end,:] .= cv_star_arr[end-1,:]
        end

    elseif bc_chi1 == "outflow_supersonic"
        @views @. rho[end,:] = 2 * rho[end-1,:] - rho[end-2,:]
        @views @. rho_E_arr[end,:] = 2 * rho_E_arr[end-1,:] - rho_E_arr[end-2,:]
        @views @. rho_u_arr[end,:] = 2 * rho_u_arr[end-1,:] - rho_u_arr[end-2,:]
        @views @. rho_v_arr[end,:] = 2 * rho_v_arr[end-1,:] - rho_v_arr[end-2,:]
        @views @. p_arr[end,:] = 2 * p_arr[end-1,:] - p_arr[end-2,:]

        if !chem_eq
            @views @. gamma_star_arr[end,:] = 2 * gamma_star_arr[end-1,:] - gamma_star_arr[end-2,:]
            @views @. cv_star_arr[end,:] = 2 * cv_star_arr[end-1,:] - cv_star_arr[end-2,:]
        end

    elseif bc_chi1 == "outflow_supersonic_1st"
        @views rho_u_arr[end,:] .= rho_u_arr[end-1,:]
        @views rho_v_arr[end,:] .= rho_v_arr[end-1,:]
        @views rho[end,:] .= rho[end-1,:]
        @views rho_E_arr[end,:] .= rho_E_arr[end-1,:]
        @views p_arr[end,:] .= p_arr[end-1,:]

        if !chem_eq
            @views gamma_star_arr[end,:] .= gamma_star_arr[end-1,:]
            @views cv_star_arr[end,:] .= cv_star_arr[end-1,:]
        end

    elseif bc_chi1 == "outflow_subsonic"
        @views @. rho_u_arr[end,:] = 2 * rho_u_arr[end-1,:] - rho_u_arr[end-2,:]
        @views @. rho_v_arr[end,:] = 2 * rho_v_arr[end-1,:] - rho_v_arr[end-2,:]
        @views @. rho[end,:] = 2 * rho[end-1,:] - rho[end-2,:]
        @views @. rho_E_arr[end,:] = 2 * rho_E_arr[end-1,:] - rho_E_arr[end-2,:]
        @views p_arr[end,:] .= p_arr[end-1,:]  # zero pressure gradient

        if !chem_eq
            @views gamma_star_arr[end,:] .= gamma_star_arr[end-1,:]
            @views cv_star_arr[end,:] .= cv_star_arr[end-1,:]
        end

    elseif bc_chi1 == "outflow_NRCBC"
        # LODI boundary conditions (Poinsot & Lele): L1 = 0
        # Extrapolate all variables (supersonic default)
        @views @. rho[end,:] = 2 * rho[end-1,:] - rho[end-2,:]
        @views @. rho_u_arr[end,:] = 2 * rho_u_arr[end-1,:] - rho_u_arr[end-2,:]
        @views @. rho_v_arr[end,:] = 2 * rho_v_arr[end-1,:] - rho_v_arr[end-2,:]
        @views @. rho_E_arr[end,:] = 2 * rho_E_arr[end-1,:] - rho_E_arr[end-2,:]
        @views @. p_arr[end,:] = 2 * p_arr[end-1,:] - p_arr[end-2,:]

        # Sound speed and density at last interior cell
        a_arr = var["a"]::Matrix{Float64}
        lr_x_normal = mesh["lr_x_normal"]::Matrix{Float64}
        lr_y_normal = mesh["lr_y_normal"]::Matrix{Float64}

        @views begin
            a = a_arr[end-1, 2:end-1]
            rho_int = rho[end-1, 2:end-1]

            # Normal velocity at ghost cell, positive exiting domain
            u_n_ghost = @. (rho_u_arr[end, 2:end-1] * lr_x_normal[end,:] +
                        rho_v_arr[end, 2:end-1] * lr_y_normal[end,:]) / rho[end, 2:end-1]

            # Normal velocity at last interior cell
            u_n_int = @. (rho_u_arr[end-1, 2:end-1] * lr_x_normal[end,:] +
                       rho_v_arr[end-1, 2:end-1] * lr_y_normal[end,:]) / rho_int

            # LODI correction for subsonic cells (L1 = 0 => dp/dx_n = rho*c*du_n/dx_n)
            subsonic = u_n_ghost .< a
            if any(subsonic)
                p_int = p_arr[end-1, 2:end-1]
                p_lodi = @. p_int + rho_int * a * (u_n_ghost - u_n_int)
                @. p_arr[end, 2:end-1] = subsonic * p_lodi + (!subsonic) * p_arr[end, 2:end-1]
                rho_e = @. p_arr[end, 2:end-1] / (gamma_star_arr[end, 2:end-1] - 1)
                rho_E_lodi = @. rho_e + 0.5 * (rho_u_arr[end, 2:end-1]^2 + rho_v_arr[end, 2:end-1]^2) / rho[end, 2:end-1]
                @. rho_E_arr[end, 2:end-1] = subsonic * rho_E_lodi + (!subsonic) * rho_E_arr[end, 2:end-1]
            end
        end

        if !chem_eq
            gamma_star_eq = var["gamma_star_eq"]::Matrix{Float64}
            cv_star_eq    = var["cv_star_eq"]::Matrix{Float64}
            @views gamma_star_arr[end,:] .= gamma_star_eq[end,:]
            @views cv_star_arr[end,:] .= cv_star_eq[end,:]
        end

    elseif bc_chi1 == "periodic"
        @views rho_u_arr[end,:] .= rho_u_arr[2,:]
        @views rho_v_arr[end,:] .= rho_v_arr[2,:]
        @views rho[end,:] .= rho[2,:]
        @views rho_E_arr[end,:] .= rho_E_arr[2,:]
        @views p_arr[end,:] .= p_arr[2,:]
        @views T_arr[end,:] .= T_arr[2,:]

        if !chem_eq
            @views gamma_star_arr[end,:] .= gamma_star_arr[2,:]
            @views cv_star_arr[end,:] .= cv_star_arr[2,:]
        end

    elseif bc_chi1 == "symmetry"
        @views rho[end,:] .= rho[end-1,:]
        @views p_arr[end,:] .= p_arr[end-1,:]

        lr_x_normal = mesh["lr_x_normal"]::Matrix{Float64}
        lr_y_normal = mesh["lr_y_normal"]::Matrix{Float64}

        # Remove wall-normal velocity component via reflection
        @views begin
            flux_vel = @. rho_u_arr[end-1, 2:end-1] * lr_x_normal[end,:] +
                       rho_v_arr[end-1, 2:end-1] * lr_y_normal[end,:]
            @. rho_u_arr[end, 2:end-1] = rho_u_arr[end-1, 2:end-1] - 2 * flux_vel * lr_x_normal[end,:]
            @. rho_v_arr[end, 2:end-1] = rho_v_arr[end-1, 2:end-1] - 2 * flux_vel * lr_y_normal[end,:]
        end
        rho_u_arr[end,1] = rho_u_arr[end,2]
        rho_v_arr[end,1] = rho_v_arr[end,2]
        rho_u_arr[end,end] = rho_u_arr[end,end-1]
        rho_v_arr[end,end] = rho_v_arr[end,end-1]

        # Remove heat flux (adiabatic symmetry)
        @views rho_E_arr[end,:] .= rho_E_arr[end-1,:]
        @views T_arr[end,:] .= T_arr[end-1,:]

        if !chem_eq
            @views gamma_star_arr[end,:] .= gamma_star_arr[end-1,:]
            @views cv_star_arr[end,:] .= cv_star_arr[end-1,:]
        end

    else
        error("Non defined Boundary-condition, chi1")
    end

    ##########################################################################
    ## Bottom boundary (eta0)
    bc_eta0 = bc["boundary_eta0"]["name"]

    if bc_eta0 == "symmetry"
        @views p_arr[:,1] .= p_arr[:,2]
        @views rho[:,1] .= rho[:,2]

        bt_x_normal = mesh["bt_x_normal"]::Matrix{Float64}
        bt_y_normal = mesh["bt_y_normal"]::Matrix{Float64}

        # Remove wall-normal velocity component via reflection
        @views begin
            flux_vel = @. rho_u_arr[2:end-1,2] * bt_x_normal[:,1] +
                       rho_v_arr[2:end-1,2] * bt_y_normal[:,1]
            @. rho_u_arr[2:end-1,1] = rho_u_arr[2:end-1,2] - 2 * flux_vel * bt_x_normal[:,1]
            @. rho_v_arr[2:end-1,1] = rho_v_arr[2:end-1,2] - 2 * flux_vel * bt_y_normal[:,1]
        end
        rho_u_arr[1,1] = rho_u_arr[2,1]
        rho_v_arr[1,1] = rho_v_arr[2,1]
        rho_u_arr[end,1] = rho_u_arr[end-1,1]
        rho_v_arr[end,1] = rho_v_arr[end-1,1]

        # Remove heat flux (adiabatic symmetry)
        @views rho_E_arr[:,1] .= rho_E_arr[:,2]
        @views T_arr[:,1] .= T_arr[:,2]

        if !chem_eq
            @views gamma_star_arr[:,1] .= gamma_star_arr[:,2]
            @views cv_star_arr[:,1] .= cv_star_arr[:,2]
        end

    elseif bc_eta0 == "no_slip_adiabatic"
        @views rho[:,1] .= rho[:,2]
        @views rho_u_arr[:,1] .= .-rho_u_arr[:,2]
        @views rho_v_arr[:,1] .= .-rho_v_arr[:,2]
        @views @. rho_E_arr[:,1] = 2 * rho_E_arr[:,2] - rho_E_arr[:,3]
        @views @. p_arr[:,1] = 2 * p_arr[:,2] - p_arr[:,3]
        @views T_arr[:,1] .= T_arr[:,2]  # No temperature gradient for adiabatic wall

        if !chem_eq
            @views @. gamma_star_arr[:,1] = 2 * gamma_star_arr[:,2] - gamma_star_arr[:,3]
            @views @. cv_star_arr[:,1] = 2 * cv_star_arr[:,2] - cv_star_arr[:,3]
        end

    elseif bc_eta0 == "no_slip_isothermal"
        @views rho[:,1] .= rho[:,2]
        @views rho_u_arr[:,1] .= .-rho_u_arr[:,2]
        @views rho_v_arr[:,1] .= .-rho_v_arr[:,2]
        @views @. rho_E_arr[:,1] = 2 * rho_E_arr[:,2] - rho_E_arr[:,3]
        @views @. p_arr[:,1] = 2 * p_arr[:,2] - p_arr[:,3]
        T_w = bc["boundary_eta0"]["Tw"] * fs["cv"] / fs["energy_factor"]  # Non-dimensional wall temperature
        @views @. T_arr[:,1] = 2 * T_w - T_arr[:,2]  # Set wall temperature

        if !chem_eq
            @views gamma_star_arr[:,1] .= gamma_star_arr[:,2]
            @views cv_star_arr[:,1] .= cv_star_arr[:,2]
        end

    elseif bc_eta0 == "outflow_supersonic"
        @views @. rho[:,1] = 2 * rho[:,2] - rho[:,3]
        @views @. rho_E_arr[:,1] = 2 * rho_E_arr[:,2] - rho_E_arr[:,3]
        @views @. rho_u_arr[:,1] = 2 * rho_u_arr[:,2] - rho_u_arr[:,3]
        @views @. rho_v_arr[:,1] = 2 * rho_v_arr[:,2] - rho_v_arr[:,3]
        @views @. p_arr[:,1] = 2 * p_arr[:,2] - p_arr[:,3]

        if !chem_eq
            @views @. gamma_star_arr[:,1] = 2 * gamma_star_arr[:,2] - gamma_star_arr[:,3]
            @views @. cv_star_arr[:,1] = 2 * cv_star_arr[:,2] - cv_star_arr[:,3]
        end

    elseif bc_eta0 == "outflow_supersonic_1st"
        @views rho[:,1] .= rho[:,2]
        @views rho_u_arr[:,1] .= rho_u_arr[:,2]
        @views rho_v_arr[:,1] .= rho_v_arr[:,2]
        @views rho_E_arr[:,1] .= rho_E_arr[:,2]
        @views p_arr[:,1] .= p_arr[:,2]

        if !chem_eq
            @views gamma_star_arr[:,1] .= gamma_star_arr[:,2]
            @views cv_star_arr[:,1] .= cv_star_arr[:,2]
        end

    elseif bc_eta0 == "outflow_subsonic"
        @views @. rho[:,1] = 2 * rho[:,2] - rho[:,3]
        @views @. rho_E_arr[:,1] = 2 * rho_E_arr[:,2] - rho_E_arr[:,3]
        @views @. rho_u_arr[:,1] = 2 * rho_u_arr[:,2] - rho_u_arr[:,3]
        @views @. rho_v_arr[:,1] = 2 * rho_v_arr[:,2] - rho_v_arr[:,3]
        @views p_arr[:,1] .= p_arr[:,2]  # zero pressure gradient

        if !chem_eq
            @views gamma_star_arr[:,1] .= gamma_star_arr[:,2]
            @views cv_star_arr[:,1] .= cv_star_arr[:,2]
        end

    elseif bc_eta0 == "outflow_NRCBC"
        # LODI boundary conditions (Poinsot & Lele): L1 = 0
        # Extrapolate all variables (supersonic default)
        @views @. rho[:,1] = 2 * rho[:,2] - rho[:,3]
        @views @. rho_u_arr[:,1] = 2 * rho_u_arr[:,2] - rho_u_arr[:,3]
        @views @. rho_v_arr[:,1] = 2 * rho_v_arr[:,2] - rho_v_arr[:,3]
        @views @. rho_E_arr[:,1] = 2 * rho_E_arr[:,2] - rho_E_arr[:,3]
        @views @. p_arr[:,1] = 2 * p_arr[:,2] - p_arr[:,3]

        # Sound speed and density at first interior cell
        a_arr = var["a"]::Matrix{Float64}
        bt_x_normal = mesh["bt_x_normal"]::Matrix{Float64}
        bt_y_normal = mesh["bt_y_normal"]::Matrix{Float64}

        @views begin
            a = a_arr[2:end-1, 2]
            rho_int = rho[2:end-1, 2]

            # Normal velocity at ghost cell, positive exiting domain
            u_n_ghost = @. -(rho_u_arr[2:end-1, 1] * bt_x_normal[:,1] +
                          rho_v_arr[2:end-1, 1] * bt_y_normal[:,1]) / rho[2:end-1, 1]

            # Normal velocity at first interior cell
            u_n_int = @. -(rho_u_arr[2:end-1, 2] * bt_x_normal[:,1] +
                        rho_v_arr[2:end-1, 2] * bt_y_normal[:,1]) / rho_int

            # LODI correction for subsonic cells (L1 = 0 => dp/dx_n = rho*c*du_n/dx_n)
            subsonic = u_n_ghost .< a
            if any(subsonic)
                p_int = p_arr[2:end-1, 2]
                p_lodi = @. p_int + rho_int * a * (u_n_ghost - u_n_int)
                # Apply LODI pressure only to subsonic cells, keep extrapolated for supersonic
                @. p_arr[2:end-1, 1] = subsonic * p_lodi + (!subsonic) * p_arr[2:end-1, 1]
                # Recompute total energy from corrected pressure
                rho_e = @. p_arr[2:end-1, 1] / (gamma_star_arr[2:end-1, 1] - 1)
                rho_E_lodi = @. rho_e + 0.5 * (rho_u_arr[2:end-1, 1]^2 + rho_v_arr[2:end-1, 1]^2) / rho[2:end-1, 1]
                @. rho_E_arr[2:end-1, 1] = subsonic * rho_E_lodi + (!subsonic) * rho_E_arr[2:end-1, 1]
            end
        end

        if !chem_eq
            gamma_star_eq = var["gamma_star_eq"]::Matrix{Float64}
            cv_star_eq    = var["cv_star_eq"]::Matrix{Float64}
            @views gamma_star_arr[:,1] .= gamma_star_eq[:,1]
            @views cv_star_arr[:,1] .= cv_star_eq[:,1]
        end

    elseif bc_eta0 == "inflow_supersonic"
        @views rho[:,1] .= rho_0
        @views rho_u_arr[:,1] .= rho_u_0
        @views rho_v_arr[:,1] .= rho_v_0
        @views rho_E_arr[:,1] .= rho_E_0
        @views p_arr[:,1] .= p_0

        if !chem_eq
            gamma_star_eq = var["gamma_star_eq"]::Matrix{Float64}
            cv_star_eq    = var["cv_star_eq"]::Matrix{Float64}
            @views gamma_star_arr[:,1] .= gamma_star_eq[:,1]
            @views cv_star_arr[:,1] .= cv_star_eq[:,1]
        end

    elseif bc_eta0 == "inflow_subsonic"
        @views rho[:,1] .= rho_0
        @views rho_u_arr[:,1] .= rho_u_0
        @views rho_v_arr[:,1] .= rho_v_0
        @views @. rho_E_arr[:,1] = 2 * rho_E_arr[:,2] - rho_E_arr[:,3]
        @views @. p_arr[:,1] = 2 * p_arr[:,2] - p_arr[:,3]  # Release acoustic characteristic for subsonic inflow

        if !chem_eq
            gamma_star_eq = var["gamma_star_eq"]::Matrix{Float64}
            cv_star_eq    = var["cv_star_eq"]::Matrix{Float64}
            @views gamma_star_arr[:,1] .= gamma_star_eq[:,1]
            @views cv_star_arr[:,1] .= cv_star_eq[:,1]
        end

    elseif bc_eta0 == "periodic"
        @views rho[:,1] .= rho[:,end-1]
        @views rho_u_arr[:,1] .= rho_u_arr[:,end-1]
        @views rho_v_arr[:,1] .= rho_v_arr[:,end-1]
        @views rho_E_arr[:,1] .= rho_E_arr[:,end-1]
        @views p_arr[:,1] .= p_arr[:,end-1]
        @views T_arr[:,1] .= T_arr[:,end-1]

        if !chem_eq
            @views gamma_star_arr[:,1] .= gamma_star_arr[:,end-1]
            @views cv_star_arr[:,1] .= cv_star_arr[:,end-1]
        end

    else
        error("Non defined Boundary-condition, eta0")
    end

    ##########################################################################
    ## Top boundary (eta1)
    bc_eta1 = bc["boundary_eta1"]["name"]

    if bc_eta1 == "symmetry"
        @views p_arr[:,end] .= p_arr[:,end-1]
        @views rho[:,end] .= rho[:,end-1]

        bt_x_normal = mesh["bt_x_normal"]::Matrix{Float64}
        bt_y_normal = mesh["bt_y_normal"]::Matrix{Float64}

        # Remove wall-normal velocity component
        @views begin
            flux_vel = @. rho_u_arr[2:end-1, end-1] * bt_x_normal[:,end] +
                       rho_v_arr[2:end-1, end-1] * bt_y_normal[:,end]
            @. rho_u_arr[2:end-1, end] = rho_u_arr[2:end-1, end-1] - 2 * flux_vel * bt_x_normal[:,end]
            @. rho_v_arr[2:end-1, end] = rho_v_arr[2:end-1, end-1] - 2 * flux_vel * bt_y_normal[:,end]
        end
        rho_u_arr[1,end] = rho_u_arr[2,end]
        rho_v_arr[1,end] = rho_v_arr[2,end]
        rho_u_arr[end,end] = rho_u_arr[end-1,end]
        rho_v_arr[end,end] = rho_v_arr[end-1,end]

        # Remove heat flux (adiabatic symmetry)
        @views rho_E_arr[:,end] .= rho_E_arr[:,end-1]
        @views T_arr[:,end] .= T_arr[:,end-1]

        if !chem_eq
            @views gamma_star_arr[:,end] .= gamma_star_arr[:,end-1]
            @views cv_star_arr[:,end] .= cv_star_arr[:,end-1]
        end

    elseif bc_eta1 == "no_slip_adiabatic"
        @views rho[:,end] .= rho[:,end-1]
        @views rho_u_arr[:,end] .= .-rho_u_arr[:,end-1]
        @views rho_v_arr[:,end] .= .-rho_v_arr[:,end-1]
        @views @. rho_E_arr[:,end] = 2 * rho_E_arr[:,end-1] - rho_E_arr[:,end-2]
        @views @. p_arr[:,end] = 2 * p_arr[:,end-1] - p_arr[:,end-2]
        @views T_arr[:,end] .= T_arr[:,end-1]  # No temperature gradient for adiabatic wall

        if !chem_eq
            @views gamma_star_arr[:,end] .= gamma_star_arr[:,end-1]
            @views cv_star_arr[:,end] .= cv_star_arr[:,end-1]
        end

    elseif bc_eta1 == "no_slip_isothermal"
        @views rho[:,end] .= rho[:,end-1]
        @views rho_u_arr[:,end] .= .-rho_u_arr[:,end-1]
        @views rho_v_arr[:,end] .= .-rho_v_arr[:,end-1]
        @views @. rho_E_arr[:,end] = 2 * rho_E_arr[:,end-1] - rho_E_arr[:,end-2]
        @views @. p_arr[:,end] = 2 * p_arr[:,end-1] - p_arr[:,end-2]
        T_w = bc["boundary_eta1"]["Tw"] * fs["cv"] / fs["energy_factor"]  # Non-dimensional wall temperature
        @views @. T_arr[:,end] = 2 * T_w - T_arr[:,end-1]  # Set wall temperature

        if !chem_eq
            @views gamma_star_arr[:,end] .= gamma_star_arr[:,end-1]
            @views cv_star_arr[:,end] .= cv_star_arr[:,end-1]
        end

    elseif bc_eta1 == "outflow_supersonic"
        @views @. rho[:,end] = 2 * rho[:,end-1] - rho[:,end-2]
        @views @. rho_E_arr[:,end] = 2 * rho_E_arr[:,end-1] - rho_E_arr[:,end-2]
        @views @. rho_u_arr[:,end] = 2 * rho_u_arr[:,end-1] - rho_u_arr[:,end-2]
        @views @. rho_v_arr[:,end] = 2 * rho_v_arr[:,end-1] - rho_v_arr[:,end-2]
        @views @. p_arr[:,end] = 2 * p_arr[:,end-1] - p_arr[:,end-2]

        if !chem_eq
            @views @. gamma_star_arr[:,end] = 2 * gamma_star_arr[:,end-1] - gamma_star_arr[:,end-2]
            @views @. cv_star_arr[:,end] = 2 * cv_star_arr[:,end-1] - cv_star_arr[:,end-2]
        end

    elseif bc_eta1 == "outflow_supersonic_1st"
        @views rho[:,end] .= rho[:,end-1]
        @views rho_u_arr[:,end] .= rho_u_arr[:,end-1]
        @views rho_v_arr[:,end] .= rho_v_arr[:,end-1]
        @views rho_E_arr[:,end] .= rho_E_arr[:,end-1]
        @views p_arr[:,end] .= p_arr[:,end-1]

        if !chem_eq
            @views gamma_star_arr[:,end] .= gamma_star_arr[:,end-1]
            @views cv_star_arr[:,end] .= cv_star_arr[:,end-1]
        end

    elseif bc_eta1 == "outflow_subsonic"
        @views @. rho[:,end] = 2 * rho[:,end-1] - rho[:,end-2]
        @views @. rho_E_arr[:,end] = 2 * rho_E_arr[:,end-1] - rho_E_arr[:,end-2]
        @views @. rho_u_arr[:,end] = 2 * rho_u_arr[:,end-1] - rho_u_arr[:,end-2]
        @views @. rho_v_arr[:,end] = 2 * rho_v_arr[:,end-1] - rho_v_arr[:,end-2]
        @views p_arr[:,end] .= p_arr[:,end-1]  # zero pressure gradient

        if !chem_eq
            @views gamma_star_arr[:,end] .= gamma_star_arr[:,end-1]
            @views cv_star_arr[:,end] .= cv_star_arr[:,end-1]
        end

    elseif bc_eta1 == "outflow_NRCBC"
        # LODI boundary conditions (Poinsot & Lele): L1 = 0
        # Extrapolate all variables (supersonic default)
        @views @. rho[:,end] = 2 * rho[:,end-1] - rho[:,end-2]
        @views @. rho_u_arr[:,end] = 2 * rho_u_arr[:,end-1] - rho_u_arr[:,end-2]
        @views @. rho_v_arr[:,end] = 2 * rho_v_arr[:,end-1] - rho_v_arr[:,end-2]
        @views @. rho_E_arr[:,end] = 2 * rho_E_arr[:,end-1] - rho_E_arr[:,end-2]
        @views @. p_arr[:,end] = 2 * p_arr[:,end-1] - p_arr[:,end-2]

        # Sound speed and density at last interior cell
        a_arr = var["a"]::Matrix{Float64}
        bt_x_normal = mesh["bt_x_normal"]::Matrix{Float64}
        bt_y_normal = mesh["bt_y_normal"]::Matrix{Float64}

        @views begin
            a = a_arr[2:end-1, end-1]
            rho_int = rho[2:end-1, end-1]

            # Normal velocity at ghost cell, positive exiting domain (bt_normal points +eta = outward at eta1)
            u_n_ghost = @. (rho_u_arr[2:end-1, end] * bt_x_normal[:,end] +
                         rho_v_arr[2:end-1, end] * bt_y_normal[:,end]) / rho[2:end-1, end]

            # Normal velocity at last interior cell
            u_n_int = @. (rho_u_arr[2:end-1, end-1] * bt_x_normal[:,end] +
                       rho_v_arr[2:end-1, end-1] * bt_y_normal[:,end]) / rho_int

            # LODI correction for subsonic cells (L1 = 0 => dp/dx_n = rho*c*du_n/dx_n)
            subsonic = u_n_ghost .< a
            if any(subsonic)
                p_int = p_arr[2:end-1, end-1]
                p_lodi = @. p_int + rho_int * a * (u_n_ghost - u_n_int)
                @. p_arr[2:end-1, end] = subsonic * p_lodi + (!subsonic) * p_arr[2:end-1, end]
                rho_e = @. p_arr[2:end-1, end] / (gamma_star_arr[2:end-1, end] - 1)
                rho_E_lodi = @. rho_e + 0.5 * (rho_u_arr[2:end-1, end]^2 + rho_v_arr[2:end-1, end]^2) / rho[2:end-1, end]
                @. rho_E_arr[2:end-1, end] = subsonic * rho_E_lodi + (!subsonic) * rho_E_arr[2:end-1, end]
            end
        end

        if !chem_eq
            gamma_star_eq = var["gamma_star_eq"]::Matrix{Float64}
            cv_star_eq    = var["cv_star_eq"]::Matrix{Float64}
            @views gamma_star_arr[:,end] .= gamma_star_eq[:,end]
            @views cv_star_arr[:,end] .= cv_star_eq[:,end]
        end

    elseif bc_eta1 == "inflow_supersonic"
        @views rho_u_arr[:,end] .= rho_u_0
        @views rho_v_arr[:,end] .= rho_v_0
        @views rho[:,end] .= rho_0
        @views rho_E_arr[:,end] .= rho_E_0
        @views p_arr[:,end] .= p_0

        if !chem_eq
            gamma_star_eq = var["gamma_star_eq"]::Matrix{Float64}
            cv_star_eq    = var["cv_star_eq"]::Matrix{Float64}
            @views gamma_star_arr[:,end] .= gamma_star_eq[:,end]
            @views cv_star_arr[:,end] .= cv_star_eq[:,end]
        end

    elseif bc_eta1 == "inflow_subsonic"
        @views rho_u_arr[:,end] .= rho_u_0
        @views rho_v_arr[:,end] .= rho_v_0
        @views rho[:,end] .= rho_0
        @views @. rho_E_arr[:,end] = 2 * rho_E_arr[:,end-1] - rho_E_arr[:,end-2]
        @views @. p_arr[:,end] = 2 * p_arr[:,end-1] - p_arr[:,end-2]

        if !chem_eq
            gamma_star_eq = var["gamma_star_eq"]::Matrix{Float64}
            cv_star_eq    = var["cv_star_eq"]::Matrix{Float64}
            @views gamma_star_arr[:,end] .= gamma_star_eq[:,end]
            @views cv_star_arr[:,end] .= cv_star_eq[:,end]
        end

    elseif bc_eta1 == "periodic"
        @views rho[:,end] .= rho[:,2]
        @views rho_u_arr[:,end] .= rho_u_arr[:,2]
        @views rho_v_arr[:,end] .= rho_v_arr[:,2]
        @views rho_E_arr[:,end] .= rho_E_arr[:,2]
        @views p_arr[:,end] .= p_arr[:,2]
        @views T_arr[:,end] .= T_arr[:,2]

        if !chem_eq
            @views gamma_star_arr[:,end] .= gamma_star_arr[:,2]
            @views cv_star_arr[:,end] .= cv_star_arr[:,2]
        end

    elseif bc_eta1 == "shock"  # Only for eta1
        # Shock update handled within UPDATE_SHOCK_BC

    elseif bc_eta1 == "lid_driven_cavity"  # Only for eta1
        # Zero pressure gradient at lid
        @views p_arr[:,end] .= p_arr[:,end-1]
        x_Ext = mesh["x_Ext"]::Matrix{Float64}

        @views begin
            e2 = @. (rho_E_arr[:,end-1] - (rho_u_arr[:,end-1]^2 + rho_v_arr[:,end-1]^2) /
                rho[:,end-1] / 2) / rho[:,end-1]
            @. rho[:,end] = rho[:,end-1] * e2 / s["e_ref"]

            # Impose lid velocity profile (quartic bump)
            u_tan = @. 16 * (x_Ext[:,end]^2) * (1 - x_Ext[:,end])^2
            u_ghost = @. 2 * u_tan - rho_u_arr[:,end-1] / rho[:,end-1]
            rho_v_arr[:,end] .= .-rho_v_arr[:,end-1]
            @. rho_u_arr[:,end] = u_ghost * rho[:,end]

            # Isothermal wall condition
            eghost = @. 2 * s["e_ref"] - e2
            @. rho_E_arr[:,end] = eghost * rho[:,end] +
                (rho_u_arr[:,end]^2 + rho_v_arr[:,end]^2) / rho[:,end] / 2
        end

        if !chem_eq
            @views @. gamma_star_arr[:,end] = 2 * gamma_star_arr[:,end-1] - gamma_star_arr[:,end-2]
            @views @. cv_star_arr[:,end] = 2 * cv_star_arr[:,end-1] - cv_star_arr[:,end-2]
        end

    else
        error("Non defined Boundary-condition, eta1")
    end

    return s
end
