function NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s::Dict{String, Any}, chemistry::Dict{String, Any})
# NON_LINEAR_DYNAMICS_NO_DISCONTINUITY - Compute nonlinear dynamics with optional shock masking.
#
# Computes the PDE fluxes and optionally applies a shock-feedback mask
# to zero out fluxes in cells containing discontinuities.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Compute PDE fluxes
    s = PDE(s, chemistry)

    ## Apply shock-feedback mask if active
    shock_dict = s["shock"]::Dict{String, Any}
    if shock_dict["enabled"] == true && shock_dict["feedback"]
        flow_cells = shock_dict["flow_cells"]::Matrix{Float64}
        flux       = s["flux"]::Dict{String, Any}
        flux_rho   = flux["rho"]::Matrix{Float64}
        flux_rho_u = flux["rho_u"]::Matrix{Float64}
        flux_rho_v = flux["rho_v"]::Matrix{Float64}
        flux_rho_E = flux["rho_E"]::Matrix{Float64}
        @. flux_rho   *= flow_cells
        @. flux_rho_u *= flow_cells
        @. flux_rho_v *= flow_cells
        @. flux_rho_E *= flow_cells
    end

    return s
end
