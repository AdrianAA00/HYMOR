function FILTERING(s::Dict{String, Any})
# FILTERING - Apply numerical dissipation filtering to conserved variables.
#
# Blends the computed flux with a 5-point Laplacian smoothing stencil to
# provide numerical dissipation. Each conserved variable (rho, rho*u,
# rho*v, rho*E) is filtered independently using its own numerical
# viscosity coefficient (mu_numerical_*).
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract variable arrays (ghost-padded)
    var_rho   = s["var"]["rho"]::Matrix{Float64}
    var_rho_u = s["var"]["rho_u"]::Matrix{Float64}
    var_rho_v = s["var"]["rho_v"]::Matrix{Float64}
    var_rho_E = s["var"]["rho_E"]::Matrix{Float64}

    ## Extract flux arrays (interior-sized)
    flux_rho   = s["flux"]["rho"]::Matrix{Float64}
    flux_rho_u = s["flux"]["rho_u"]::Matrix{Float64}
    flux_rho_v = s["flux"]["rho_v"]::Matrix{Float64}
    flux_rho_E = s["flux"]["rho_E"]::Matrix{Float64}

    ## Extract numerical dissipation coefficients
    mu_rho   = s["numerical_dissipation"]["mu_rho"]::Float64
    mu_rho_u = s["numerical_dissipation"]["mu_rho_u"]::Float64
    mu_rho_v = s["numerical_dissipation"]["mu_rho_v"]::Float64
    mu_rho_E = s["numerical_dissipation"]["mu_rho_E"]::Float64

    ## Filter density flux
    @views @. flux_rho = flux_rho * (1 - mu_rho) + mu_rho * (4 * var_rho[2:end-1, 2:end-1] + var_rho[1:end-2, 2:end-1] + var_rho[2:end-1, 1:end-2] + var_rho[3:end, 2:end-1] + var_rho[2:end-1, 3:end]) / 8

    ## Filter x-momentum flux
    @views @. flux_rho_u = flux_rho_u * (1 - mu_rho_u) + mu_rho_u * (4 * var_rho_u[2:end-1, 2:end-1] + var_rho_u[1:end-2, 2:end-1] + var_rho_u[2:end-1, 1:end-2] + var_rho_u[3:end, 2:end-1] + var_rho_u[2:end-1, 3:end]) / 8

    ## Filter y-momentum flux
    @views @. flux_rho_v = flux_rho_v * (1 - mu_rho_v) + mu_rho_v * (4 * var_rho_v[2:end-1, 2:end-1] + var_rho_v[1:end-2, 2:end-1] + var_rho_v[2:end-1, 1:end-2] + var_rho_v[3:end, 2:end-1] + var_rho_v[2:end-1, 3:end]) / 8

    ## Filter total energy flux
    @views @. flux_rho_E = flux_rho_E * (1 - mu_rho_E) + mu_rho_E * (4 * var_rho_E[2:end-1, 2:end-1] + var_rho_E[1:end-2, 2:end-1] + var_rho_E[2:end-1, 1:end-2] + var_rho_E[3:end, 2:end-1] + var_rho_E[2:end-1, 3:end]) / 8

    return s
end
