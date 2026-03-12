using LinearAlgebra

function COMPUTE_BUDGETS_KINETIC(V, s, output_flow)
    # COMPUTE_BUDGETS_KINETIC  Compute kinetic energy budget terms for a perturbation field.
    #
    #   budgets_kinetic = COMPUTE_BUDGETS_KINETIC(V, s, output_flow)
    #
    #   Computes the individual terms in the kinetic energy component of the Chu
    #   energy budget equation. The kinetic budget decomposes the time rate of
    #   change of perturbation kinetic energy into advection, momentum and mass
    #   production, pressure-dilatation, viscous and pressure transport, and
    #   viscous dissipation. Both 2-D Cartesian and 3-D axisymmetric formulations
    #   are supported.
    #
    #   Inputs:
    #       V           - (Vector) Global perturbation state vector
    #                     containing [rho'; rho_u'; rho_v'; rho_E'] stacked over
    #                     the Nx-by-Ny grid (length >= 4*Nx*Ny).
    #       s           - (Dict) Solution structure with base-flow fields.
    #       output_flow - (Bool) If true, the full 2-D field of every budget
    #                     term is stored in the output dict.
    #
    #   Outputs:
    #       budgets_kinetic - (Dict{String,Any}) containing volume-integrated
    #                         budget terms and optionally the spatially resolved 2-D fields.

    N_v = s["mesh"]["Nchi"] * s["mesh"]["Neta"]

    ## Extract perturbations from state vector
    pert_rho   = real.(reshape(V[0*N_v+1 : 1*N_v, 1], s["mesh"]["Nchi"], s["mesh"]["Neta"]))
    pert_rho_u = real.(reshape(V[1*N_v+1 : 2*N_v, 1], s["mesh"]["Nchi"], s["mesh"]["Neta"]))
    pert_rho_v = real.(reshape(V[2*N_v+1 : 3*N_v, 1], s["mesh"]["Nchi"], s["mesh"]["Neta"]))
    pert_rho_E = real.(reshape(V[3*N_v+1 : 4*N_v, 1], s["mesh"]["Nchi"], s["mesh"]["Neta"]))

    ## Extract base-flow variables
    rho_0   = s["var"]["rho"][2:end-1, 2:end-1]
    rho_u_0 = s["var"]["rho_u"][2:end-1, 2:end-1]
    rho_v_0 = s["var"]["rho_v"][2:end-1, 2:end-1]
    u_0 = rho_u_0 ./ rho_0
    v_0 = s["var"]["rho_v"][2:end-1, 2:end-1] ./ rho_0
    p_0 = s["var"]["p"][2:end-1, 2:end-1]
    a_0 = s["var"]["a"][2:end-1, 2:end-1]
    gamma_star_0 = s["var"]["gamma_star"][2:end-1, 2:end-1]

    u_0_Ext = EXTEND_TO_GHOST_POINTS(s["var"]["rho_u"][2:end-1, 2:end-1] ./ s["var"]["rho"][2:end-1, 2:end-1], s)
    v_0_Ext = EXTEND_TO_GHOST_POINTS(s["var"]["rho_v"][2:end-1, 2:end-1] ./ s["var"]["rho"][2:end-1, 2:end-1], s)

    ## Non-dimensional viscosity
    mu_0 = s["var"]["mu_star"][2:end-1, 2:end-1] ./ s["freestream"]["Re"]

    volume     = s["mesh"]["volume"]
    flow_cells = s["shock"]["flow_cells"]

    ## Compute base-flow derivatives
    (du0_dx, du0_dy) = DERIVATIVE_EXT(u_0_Ext, s)
    (dv0_dx, dv0_dy) = DERIVATIVE_EXT(v_0_Ext, s)

    ## Transform perturbations to Chu variables
    (pert_p, pert_u, pert_v, pert_entropy) = GET_CHU_COMPONENTS(
        rho_0, rho_u_0, rho_v_0, p_0, a_0, gamma_star_0,
        pert_rho, pert_rho_u, pert_rho_v, pert_rho_E)

    pert_u_ext = EXTEND_TO_GHOST_POINTS(pert_u, s)
    pert_v_ext = EXTEND_TO_GHOST_POINTS(pert_v, s)
    pert_p_ext = EXTEND_TO_GHOST_POINTS(pert_p, s)

    (dpert_u_dx, dpert_u_dy) = DERIVATIVE_EXT(pert_u_ext, s)
    (dpert_v_dx, dpert_v_dy) = DERIVATIVE_EXT(pert_v_ext, s)

    x = s["mesh"]["x"]  # Radial coordinate at cell centers

    ## Compute budget terms
    if s["PDE_dimension"] == "3D-axisymmetric"
        # ========== ADVECTION TERM ==========
        # A_adv = -rho_0 * u_{j,0} * d/dx_j(u'_i * u'_i / 2)
        A_adv = -rho_0 .* u_0 .* pert_u .* dpert_u_dx .-
                 rho_0 .* u_0 .* pert_v .* dpert_v_dx .-
                 rho_0 .* v_0 .* pert_u .* dpert_u_dy .-
                 rho_0 .* v_0 .* pert_v .* dpert_v_dy

        # ========== PRODUCTION TERMS ==========
        # Axisymmetric: x = radial (r), y = axial (z)
        P_mom = -rho_0 .* (pert_u .* pert_u .* du0_dx .+
                           pert_u .* pert_v .* du0_dy .+
                           pert_v .* pert_u .* dv0_dx .+
                           pert_v .* pert_v .* dv0_dy)

        P_mass = -pert_rho .* (pert_u .* u_0 .* du0_dx .+
                               pert_u .* v_0 .* du0_dy .+
                               pert_v .* u_0 .* dv0_dx .+
                               pert_v .* v_0 .* dv0_dy)

        # ========== PRESSURE-DILATATION ==========
        # Pi_d = p' * (du'/dx + u'/x + dv'/dy)
        Dilat_P = pert_p .* (dpert_u_dx .+ pert_u ./ x .+ dpert_v_dy)

        # ========== VISCOUS STRESS TENSOR ==========
        div_pert = dpert_u_dx .+ pert_u ./ x .+ dpert_v_dy

        tau_xx          = mu_0 .* (2 .* dpert_u_dx   .- 2/3 .* div_pert)
        tau_yy          = mu_0 .* (2 .* dpert_v_dy   .- 2/3 .* div_pert)
        tau_theta_theta = mu_0 .* (2 .* pert_u ./ x  .- 2/3 .* div_pert)
        tau_xy          = mu_0 .* (dpert_u_dy .+ dpert_v_dx)

        # ========== TRANSPORT TERM ==========
        # Pressure transport fluxes (with axisymmetric geometric factor)
        Flux_x_p   = x .* (-pert_u .* pert_p)
        Flux_x_tau = x .* (pert_u .* tau_xx .+ pert_v .* tau_xy)
        Flux_y_p   = -pert_v .* pert_p
        Flux_y_tau = pert_u .* tau_xy .+ pert_v .* tau_yy

        Flux_x_p_Ext   = EXTEND_TO_GHOST_POINTS(Flux_x_p, s)
        Flux_x_tau_Ext = EXTEND_TO_GHOST_POINTS(Flux_x_tau, s)
        Flux_y_p_Ext   = EXTEND_TO_GHOST_POINTS(Flux_y_p, s)
        Flux_y_tau_Ext = EXTEND_TO_GHOST_POINTS(Flux_y_tau, s)

        (dFlux_x_p_dx, _)   = DERIVATIVE_EXT(Flux_x_p_Ext, s)
        (dFlux_x_tau_dx, _)  = DERIVATIVE_EXT(Flux_x_tau_Ext, s)
        (_, dFlux_y_p_dy)    = DERIVATIVE_EXT(Flux_y_p_Ext, s)
        (_, dFlux_y_tau_dy)  = DERIVATIVE_EXT(Flux_y_tau_Ext, s)

        Transport_p   = (1 ./ x) .* dFlux_x_p_dx   .+ dFlux_y_p_dy
        Transport_tau = (1 ./ x) .* dFlux_x_tau_dx .+ dFlux_y_tau_dy

        # ========== DISSIPATION TERM ==========
        Dissipation = -(tau_xx          .* dpert_u_dx                   .+
                        tau_theta_theta .* (pert_u ./ x)               .+
                        tau_yy          .* dpert_v_dy                   .+
                        tau_xy          .* (dpert_u_dy .+ dpert_v_dx))

    else  # Cartesian 2D
        # ========== ADVECTION TERM ==========
        A_adv = -rho_0 .* u_0 .* pert_u .* dpert_u_dx .-
                 rho_0 .* u_0 .* pert_v .* dpert_v_dx .-
                 rho_0 .* v_0 .* pert_u .* dpert_u_dy .-
                 rho_0 .* v_0 .* pert_v .* dpert_v_dy

        # ========== PRODUCTION TERMS ==========
        P_mom = -rho_0 .* (pert_u .* pert_u .* du0_dx .+
                           pert_u .* pert_v .* du0_dy .+
                           pert_v .* pert_u .* dv0_dx .+
                           pert_v .* pert_v .* dv0_dy)

        P_mass = -pert_rho .* (pert_u .* u_0 .* du0_dx .+
                               pert_u .* v_0 .* du0_dy .+
                               pert_v .* u_0 .* dv0_dx .+
                               pert_v .* v_0 .* dv0_dy)

        # ========== PRESSURE-DILATATION ==========
        Dilat_P = pert_p .* (dpert_u_dx .+ dpert_v_dy)

        # ========== VISCOUS STRESS TENSOR ==========
        div_pert = dpert_u_dx .+ dpert_v_dy

        tau_xx = mu_0 .* (2 .* dpert_u_dx .- 2/3 .* div_pert)
        tau_yy = mu_0 .* (2 .* dpert_v_dy .- 2/3 .* div_pert)
        tau_xy = mu_0 .* (dpert_u_dy .+ dpert_v_dx)

        # ========== TRANSPORT TERM ==========
        Flux_x_p   = -pert_u .* pert_p
        Flux_x_tau = pert_u .* tau_xx .+ pert_v .* tau_xy
        Flux_y_p   = -pert_v .* pert_p
        Flux_y_tau = pert_u .* tau_xy .+ pert_v .* tau_yy

        Flux_x_p_Ext   = EXTEND_TO_GHOST_POINTS(Flux_x_p, s)
        Flux_x_tau_Ext = EXTEND_TO_GHOST_POINTS(Flux_x_tau, s)
        Flux_y_p_Ext   = EXTEND_TO_GHOST_POINTS(Flux_y_p, s)
        Flux_y_tau_Ext = EXTEND_TO_GHOST_POINTS(Flux_y_tau, s)

        (dFlux_x_p_dx, _)   = DERIVATIVE_EXT(Flux_x_p_Ext, s)
        (dFlux_x_tau_dx, _)  = DERIVATIVE_EXT(Flux_x_tau_Ext, s)
        (_, dFlux_y_p_dy)    = DERIVATIVE_EXT(Flux_y_p_Ext, s)
        (_, dFlux_y_tau_dy)  = DERIVATIVE_EXT(Flux_y_tau_Ext, s)

        Transport_p   = dFlux_x_p_dx   .+ dFlux_y_p_dy
        Transport_tau = dFlux_x_tau_dx .+ dFlux_y_tau_dy

        # ========== DISSIPATION TERM ==========
        Dissipation = -(tau_xx .* dpert_u_dx                   .+
                        tau_yy .* dpert_v_dy                   .+
                        tau_xy .* (dpert_u_dy .+ dpert_v_dx))
    end

    ## Compute volume-integrated budgets
    budgets_kinetic = Dict{String,Any}()
    budgets_kinetic["A_adv_sum"]         = sum(A_adv         .* volume .* flow_cells)
    budgets_kinetic["P_mom_sum"]         = sum(P_mom         .* volume .* flow_cells)
    budgets_kinetic["P_mass_sum"]        = sum(P_mass        .* volume .* flow_cells)
    budgets_kinetic["Dilat_P_sum"]       = sum(Dilat_P       .* volume .* flow_cells)
    budgets_kinetic["Transport_p_sum"]   = sum(Transport_p   .* volume .* flow_cells)
    budgets_kinetic["Transport_tau_sum"] = sum(Transport_tau  .* volume .* flow_cells)
    budgets_kinetic["Dissipation_sum"]   = sum(Dissipation   .* volume .* flow_cells)

    ## Optionally store full spatial fields
    if output_flow
        budgets_kinetic["A_adv"]         = A_adv
        budgets_kinetic["P_mom"]         = P_mom
        budgets_kinetic["P_mass"]        = P_mass
        budgets_kinetic["Dilat_P"]       = Dilat_P
        budgets_kinetic["Transport_p"]   = Transport_p
        budgets_kinetic["Transport_tau"] = Transport_tau
        budgets_kinetic["Dissipation"]   = Dissipation
    end

    return budgets_kinetic
end
