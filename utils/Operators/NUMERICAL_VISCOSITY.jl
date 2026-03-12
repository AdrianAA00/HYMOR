function NUMERICAL_VISCOSITY(s::Dict{String, Any})
# NUMERICAL_VISCOSITY - Add artificial numerical viscosity to all conserved variable fluxes.
#
# Applies a Laplacian-type artificial dissipation to the density,
# momentum, and total energy flux fields.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract arrays to avoid repeated dictionary lookups
    var_rho   = s["var"]["rho"]::Matrix{Float64}
    var_rho_u = s["var"]["rho_u"]::Matrix{Float64}
    var_rho_v = s["var"]["rho_v"]::Matrix{Float64}
    var_rho_E = s["var"]["rho_E"]::Matrix{Float64}
    flux_rho   = s["flux"]["rho"]::Matrix{Float64}
    flux_rho_u = s["flux"]["rho_u"]::Matrix{Float64}
    flux_rho_v = s["flux"]["rho_v"]::Matrix{Float64}
    flux_rho_E = s["flux"]["rho_E"]::Matrix{Float64}
    bt_l  = s["bt_length"]::Matrix{Float64}
    lr_l  = s["lr_length"]::Matrix{Float64}
    vol   = s["mesh"]["volume"]::Matrix{Float64}
    mu_rho   = s["numerical_dissipation"]["mu_rho"]::Float64
    mu_rho_u = s["numerical_dissipation"]["mu_rho_u"]::Float64
    mu_rho_v = s["numerical_dissipation"]["mu_rho_v"]::Float64
    mu_rho_E = s["numerical_dissipation"]["mu_rho_E"]::Float64

    ## Density numerical viscosity
    @views @. flux_rho = flux_rho + mu_rho * (
        -var_rho[2:end-1, 2:end-1] * (bt_l[:, 1:end-1] + bt_l[:, 2:end] + lr_l[1:end-1, :] + lr_l[2:end, :]) +
         var_rho[1:end-2, 2:end-1] * lr_l[1:end-1, :] +
         var_rho[2:end-1, 1:end-2] * bt_l[:, 1:end-1] +
         var_rho[3:end, 2:end-1]   * lr_l[2:end, :] +
         var_rho[2:end-1, 3:end]   * bt_l[:, 2:end]
    ) / vol

    ## x-momentum numerical viscosity
    @views @. flux_rho_u = flux_rho_u + mu_rho_u * (
        -var_rho_u[2:end-1, 2:end-1] * (bt_l[:, 1:end-1] + bt_l[:, 2:end] + lr_l[1:end-1, :] + lr_l[2:end, :]) +
         var_rho_u[1:end-2, 2:end-1] * lr_l[1:end-1, :] +
         var_rho_u[2:end-1, 1:end-2] * bt_l[:, 1:end-1] +
         var_rho_u[3:end, 2:end-1]   * lr_l[2:end, :] +
         var_rho_u[2:end-1, 3:end]   * bt_l[:, 2:end]
    ) / vol

    ## y-momentum numerical viscosity
    @views @. flux_rho_v = flux_rho_v + mu_rho_v * (
        -var_rho_v[2:end-1, 2:end-1] * (bt_l[:, 1:end-1] + bt_l[:, 2:end] + lr_l[1:end-1, :] + lr_l[2:end, :]) +
         var_rho_v[1:end-2, 2:end-1] * lr_l[1:end-1, :] +
         var_rho_v[2:end-1, 1:end-2] * bt_l[:, 1:end-1] +
         var_rho_v[3:end, 2:end-1]   * lr_l[2:end, :] +
         var_rho_v[2:end-1, 3:end]   * bt_l[:, 2:end]
    ) / vol

    ## Total energy numerical viscosity
    @views @. flux_rho_E = flux_rho_E + mu_rho_E * (
        -var_rho_E[2:end-1, 2:end-1] * (bt_l[:, 1:end-1] + bt_l[:, 2:end] + lr_l[1:end-1, :] + lr_l[2:end, :]) +
         var_rho_E[1:end-2, 2:end-1] * lr_l[1:end-1, :] +
         var_rho_E[2:end-1, 1:end-2] * bt_l[:, 1:end-1] +
         var_rho_E[3:end, 2:end-1]   * lr_l[2:end, :] +
         var_rho_E[2:end-1, 3:end]   * bt_l[:, 2:end]
    ) / vol

    return s
end
