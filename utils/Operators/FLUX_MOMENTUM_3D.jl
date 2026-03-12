function FLUX_MOMENTUM_3D(s::Dict{String, Any})
# FLUX_MOMENTUM_3D - Compute inviscid momentum fluxes for the 3D axisymmetric Euler equations.
#
# Includes east/west, south/north faces plus front/back pressure term for
# x-momentum arising from the axisymmetric geometry.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract arrays to avoid repeated dictionary lookups
    rho    = s["var"]["rho"]::Matrix{Float64}
    rho_u  = s["var"]["rho_u"]::Matrix{Float64}
    rho_v  = s["var"]["rho_v"]::Matrix{Float64}
    p_arr  = s["var"]["p"]::Matrix{Float64}
    lr_xn  = s["mesh"]["lr_x_normal"]::Matrix{Float64}
    lr_yn  = s["mesh"]["lr_y_normal"]::Matrix{Float64}
    lr_a   = s["mesh"]["lr_area"]::Matrix{Float64}
    bt_xn  = s["mesh"]["bt_x_normal"]::Matrix{Float64}
    bt_yn  = s["mesh"]["bt_y_normal"]::Matrix{Float64}
    bt_a   = s["mesh"]["bt_area"]::Matrix{Float64}
    vol    = s["mesh"]["volume"]::Matrix{Float64}
    fb_yn_front = s["mesh"]["fb_y_normal_front"]::Matrix{Float64}
    fb_yn_back  = s["mesh"]["fb_y_normal_back"]::Matrix{Float64}
    fb_a   = s["mesh"]["fb_area"]::Matrix{Float64}
    y_mesh = s["mesh"]["y"]::Matrix{Float64}
    flux_rho_u = s["flux"]["rho_u"]::Matrix{Float64}
    flux_rho_v = s["flux"]["rho_v"]::Matrix{Float64}

    ## X-momentum (u-component)

    # East and west faces - midpoint interpolation (2nd order on rectangles)
    @views temp_rho_u_face = @. (rho_u[1:end-1, 2:end-1] + rho_u[2:end, 2:end-1]) / 2
    @views temp_P_face = @. (p_arr[1:end-1, 2:end-1] + p_arr[2:end, 2:end-1]) * lr_xn / 2
    @views temp_vel_face = @. (rho_u[1:end-1, 2:end-1] + rho_u[2:end, 2:end-1]) /
                               (rho[2:end, 2:end-1] + rho[1:end-1, 2:end-1]) * lr_xn +
                              (rho_v[1:end-1, 2:end-1] + rho_v[2:end, 2:end-1]) /
                               (rho[2:end, 2:end-1] + rho[1:end-1, 2:end-1]) * lr_yn

    # East and west faces - midpoint quadrature (2nd order on rectangles)
    @views @. flux_rho_u = -(temp_rho_u_face[2:end, :] * temp_vel_face[2:end, :] + temp_P_face[2:end, :]) * lr_a[2:end, :]            # Outgoing fluxes
    @views @. flux_rho_u = flux_rho_u + (temp_rho_u_face[1:end-1, :] * temp_vel_face[1:end-1, :] + temp_P_face[1:end-1, :]) * lr_a[1:end-1, :]  # Incoming fluxes

    # South and north faces - midpoint interpolation (2nd order on rectangles)
    @views temp_rho_u_face2 = @. (rho_u[2:end-1, 1:end-1] + rho_u[2:end-1, 2:end]) / 2
    @views temp_P_face2 = @. (p_arr[2:end-1, 1:end-1] + p_arr[2:end-1, 2:end]) * bt_xn / 2
    @views temp_vel_face2 = @. (rho_u[2:end-1, 1:end-1] + rho_u[2:end-1, 2:end]) /
                                (rho[2:end-1, 2:end] + rho[2:end-1, 1:end-1]) * bt_xn +
                               (rho_v[2:end-1, 1:end-1] + rho_v[2:end-1, 2:end]) /
                                (rho[2:end-1, 2:end] + rho[2:end-1, 1:end-1]) * bt_yn

    # South and north faces - midpoint quadrature (2nd order on rectangles)
    @views @. flux_rho_u = flux_rho_u - (temp_rho_u_face2[:, 2:end] * temp_vel_face2[:, 2:end] + temp_P_face2[:, 2:end]) * bt_a[:, 2:end]            # Outgoing fluxes
    @views @. flux_rho_u = flux_rho_u + (temp_rho_u_face2[:, 1:end-1] * temp_vel_face2[:, 1:end-1] + temp_P_face2[:, 1:end-1]) * bt_a[:, 1:end-1]  # Incoming fluxes

    # Normalise by cell volume
    @. flux_rho_u = flux_rho_u / vol

    ## Y-momentum (v-component)

    # East and west faces - midpoint interpolation (2nd order on rectangles)
    @views temp_rho_v_face = @. (rho_v[1:end-1, 2:end-1] + rho_v[2:end, 2:end-1]) / 2
    @views temp_P_face = @. (p_arr[1:end-1, 2:end-1] + p_arr[2:end, 2:end-1]) * lr_yn / 2

    # East and west faces - midpoint quadrature (2nd order on rectangles)
    @views @. flux_rho_v = -(temp_rho_v_face[2:end, :] * temp_vel_face[2:end, :] + temp_P_face[2:end, :]) * lr_a[2:end, :]            # Outgoing fluxes
    @views @. flux_rho_v = flux_rho_v + (temp_rho_v_face[1:end-1, :] * temp_vel_face[1:end-1, :] + temp_P_face[1:end-1, :]) * lr_a[1:end-1, :]  # Incoming fluxes

    # South and north faces - midpoint interpolation (2nd order on rectangles)
    @views temp_rho_v_face2 = @. (rho_v[2:end-1, 1:end-1] + rho_v[2:end-1, 2:end]) / 2
    @views temp_P_face2 = @. (p_arr[2:end-1, 1:end-1] + p_arr[2:end-1, 2:end]) * bt_yn / 2

    # South and north faces - midpoint quadrature (2nd order on rectangles)
    @views @. flux_rho_v = flux_rho_v - (temp_rho_v_face2[:, 2:end] * temp_vel_face2[:, 2:end] + temp_P_face2[:, 2:end]) * bt_a[:, 2:end]            # Outgoing fluxes
    @views @. flux_rho_v = flux_rho_v + (temp_rho_v_face2[:, 1:end-1] * temp_vel_face2[:, 1:end-1] + temp_P_face2[:, 1:end-1]) * bt_a[:, 1:end-1]  # Incoming fluxes

    # Front and back faces (3D axisymmetric pressure contribution)
    @views temp_P_face_front = @. p_arr[2:end-1, 2:end-1] * fb_yn_front
    @views temp_P_face_back = @. p_arr[2:end-1, 2:end-1] * fb_yn_back
    @. flux_rho_v = flux_rho_v + (temp_P_face_front - temp_P_face_back) * fb_a * sign(y_mesh)

    # Normalise by cell volume
    @. flux_rho_v = flux_rho_v / vol

    return s
end
