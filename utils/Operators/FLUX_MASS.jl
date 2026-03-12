function FLUX_MASS(s::Dict{String, Any})
# FLUX_MASS - Compute the inviscid mass flux for the 2D continuity equation.
#
# Evaluates the net mass flux through all four faces (east, west, south,
# north) of each interior cell using midpoint interpolation and midpoint
# quadrature, both second-order accurate on rectangular cells.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract arrays to avoid repeated dictionary lookups
    rho_u  = s["var"]["rho_u"]::Matrix{Float64}
    rho_v  = s["var"]["rho_v"]::Matrix{Float64}
    lr_xn  = s["mesh"]["lr_x_normal"]::Matrix{Float64}
    lr_yn  = s["mesh"]["lr_y_normal"]::Matrix{Float64}
    lr_a   = s["mesh"]["lr_area"]::Matrix{Float64}
    bt_xn  = s["mesh"]["bt_x_normal"]::Matrix{Float64}
    bt_yn  = s["mesh"]["bt_y_normal"]::Matrix{Float64}
    bt_a   = s["mesh"]["bt_area"]::Matrix{Float64}
    vol    = s["mesh"]["volume"]::Matrix{Float64}
    flux_rho = s["flux"]["rho"]::Matrix{Float64}

    ## East and west faces
    # Midpoint interpolation (2nd order on rectangles) - fused single pass
    @views temp_rho_u_face = @. (rho_u[1:end-1, 2:end-1] + rho_u[2:end, 2:end-1]) * lr_xn / 2 +
                                 (rho_v[1:end-1, 2:end-1] + rho_v[2:end, 2:end-1]) * lr_yn / 2

    # Midpoint quadrature (2nd order on rectangles)
    @views @. flux_rho = -temp_rho_u_face[2:end, :] * lr_a[2:end, :]            # Outgoing fluxes
    @views @. flux_rho = flux_rho + temp_rho_u_face[1:end-1, :] * lr_a[1:end-1, :]  # Incoming fluxes

    ## South and north faces
    # Midpoint interpolation (2nd order on rectangles) - fused single pass
    @views temp_rho_u_face2 = @. (rho_u[2:end-1, 1:end-1] + rho_u[2:end-1, 2:end]) * bt_xn / 2 +
                                  (rho_v[2:end-1, 1:end-1] + rho_v[2:end-1, 2:end]) * bt_yn / 2

    # Midpoint quadrature (2nd order on rectangles)
    @views @. flux_rho = flux_rho - temp_rho_u_face2[:, 2:end] * bt_a[:, 2:end]            # Outgoing fluxes
    @views @. flux_rho = flux_rho + temp_rho_u_face2[:, 1:end-1] * bt_a[:, 1:end-1]  # Incoming fluxes

    ## Normalise by cell volume
    @. flux_rho = flux_rho / vol

    return s
end
