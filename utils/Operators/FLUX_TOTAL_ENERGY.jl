function FLUX_TOTAL_ENERGY(s::Dict{String, Any})
# FLUX_TOTAL_ENERGY - Compute the inviscid total energy flux for the 2D Euler equations.
#
# Evaluates the net total energy flux through all four faces using
# midpoint interpolation and midpoint quadrature. The energy flux is
# (rho_E + p) * V_n.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract arrays to avoid repeated dictionary lookups
    rho    = s["var"]["rho"]::Matrix{Float64}
    rho_u  = s["var"]["rho_u"]::Matrix{Float64}
    rho_v  = s["var"]["rho_v"]::Matrix{Float64}
    rho_E  = s["var"]["rho_E"]::Matrix{Float64}
    p_arr  = s["var"]["p"]::Matrix{Float64}
    lr_xn  = s["mesh"]["lr_x_normal"]::Matrix{Float64}
    lr_yn  = s["mesh"]["lr_y_normal"]::Matrix{Float64}
    lr_a   = s["mesh"]["lr_area"]::Matrix{Float64}
    bt_xn  = s["mesh"]["bt_x_normal"]::Matrix{Float64}
    bt_yn  = s["mesh"]["bt_y_normal"]::Matrix{Float64}
    bt_a   = s["mesh"]["bt_area"]::Matrix{Float64}
    vol    = s["mesh"]["volume"]::Matrix{Float64}
    flux_rhoE = s["flux"]["rho_E"]::Matrix{Float64}

    ## East and west faces
    # Midpoint interpolation (2nd order on rectangles)
    @views temp_rho_E_face = @. (rho_E[1:end-1, 2:end-1] + rho_E[2:end, 2:end-1]) / 2
    @views temp_P_face = @. (p_arr[1:end-1, 2:end-1] + p_arr[2:end, 2:end-1]) / 2
    @views temp_vel_face = @. (rho_u[1:end-1, 2:end-1] + rho_u[2:end, 2:end-1]) /
                               (rho[2:end, 2:end-1] + rho[1:end-1, 2:end-1]) * lr_xn +
                              (rho_v[1:end-1, 2:end-1] + rho_v[2:end, 2:end-1]) /
                               (rho[2:end, 2:end-1] + rho[1:end-1, 2:end-1]) * lr_yn

    # Midpoint quadrature (2nd order on rectangles)
    @views @. flux_rhoE = -((temp_rho_E_face[2:end, :] + temp_P_face[2:end, :]) * temp_vel_face[2:end, :]) * lr_a[2:end, :]            # Outgoing fluxes
    @views @. flux_rhoE = flux_rhoE + ((temp_rho_E_face[1:end-1, :] + temp_P_face[1:end-1, :]) * temp_vel_face[1:end-1, :]) * lr_a[1:end-1, :]  # Incoming fluxes

    ## South and north faces
    # Midpoint interpolation (2nd order on rectangles)
    @views temp_rho_E_face2 = @. (rho_E[2:end-1, 1:end-1] + rho_E[2:end-1, 2:end]) / 2
    @views temp_P_face2 = @. (p_arr[2:end-1, 1:end-1] + p_arr[2:end-1, 2:end]) / 2
    @views temp_vel_face2 = @. (rho_u[2:end-1, 1:end-1] + rho_u[2:end-1, 2:end]) /
                                (rho[2:end-1, 2:end] + rho[2:end-1, 1:end-1]) * bt_xn +
                               (rho_v[2:end-1, 1:end-1] + rho_v[2:end-1, 2:end]) /
                                (rho[2:end-1, 2:end] + rho[2:end-1, 1:end-1]) * bt_yn

    # Midpoint quadrature (2nd order on rectangles)
    @views @. flux_rhoE = flux_rhoE - ((temp_rho_E_face2[:, 2:end] + temp_P_face2[:, 2:end]) * temp_vel_face2[:, 2:end]) * bt_a[:, 2:end]            # Outgoing fluxes
    @views @. flux_rhoE = flux_rhoE + ((temp_rho_E_face2[:, 1:end-1] + temp_P_face2[:, 1:end-1]) * temp_vel_face2[:, 1:end-1]) * bt_a[:, 1:end-1]  # Incoming fluxes

    ## Normalise by cell volume
    @. flux_rhoE = flux_rhoE / vol

    return s
end
