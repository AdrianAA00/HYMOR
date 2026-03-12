# INTEGRATION  Apply a weighted Euler update to the conserved variables.
#
#   s = INTEGRATION(s, dsolution_dt, w)
#
#   Advances each conserved field (rho, rho_u, rho_v, rho_E, and
#   optionally gamma_star and cv_star for non-equilibrium chemistry) by
#   the increment  w * flux.  When shock fitting is active the fluxes
#   are masked by s["shock"]["flow_cells"] and the shock position is updated
#   via MOVE_SHOCK.
#
#   Inputs:
#       s              (Dict{String,Any}) - Current flow s with conservative
#                                           variables and grid / shock data.
#       dsolution_dt   (Dict{String,Any}) - Time-derivative (flux) data containing
#                                           ["flux"]["rho"], ["flux"]["rho_u"],
#                                           ["flux"]["rho_v"], ["flux"]["rho_E"],
#                                           and (when applicable)
#                                           ["flux"]["gamma_star"], ["flux"]["cv_star"].
#       w              (Float64)          - Integration weight (e.g. dt for forward
#                                           Euler, or a Runge-Kutta stage weight).
#
#   Outputs:
#       s              (Dict{String,Any}) - Updated s with incremented
#                                           conservative variables (interior cells
#                                           only, indices 2:end-1).
#
#   Notes:
#       - For shock-fitted runs (s["shock"]["enabled"] == true), fluxes are
#         multiplied by s["shock"]["flow_cells"] to zero out non-flow regions,
#         and MOVE_SHOCK is called to update the shock geometry.
#       - Non-equilibrium chemistry fields (gamma_star, cv_star) are
#         integrated only when s["chemistry"]["chemical_equilibrium"] is false.
#
# Part of: Hypersonics Stability Julia Solver - Time Marching Module

function INTEGRATION(s::Dict{String,Any}, dsolution_dt::Dict{String,Any}, w::Float64)

    ## Extract dict refs with type assertions for compiler optimization
    var  = s["var"]::Dict{String,Any}
    flux = dsolution_dt["flux"]::Dict{String,Any}
    chem_eq = s["chemistry"]["chemical_equilibrium"]::Bool

    rho_arr   = var["rho"]::Matrix{Float64}
    rho_u_arr = var["rho_u"]::Matrix{Float64}
    rho_v_arr = var["rho_v"]::Matrix{Float64}
    rho_E_arr = var["rho_E"]::Matrix{Float64}

    fl_rho   = flux["rho"]::Matrix{Float64}
    fl_rho_u = flux["rho_u"]::Matrix{Float64}
    fl_rho_v = flux["rho_v"]::Matrix{Float64}
    fl_rho_E = flux["rho_E"]::Matrix{Float64}

    @views if s["shock"]["enabled"]::Bool
        ## Shock-fitted integration (masked by flow_cells)
        flow_cells = s["shock"]["flow_cells"]::Matrix{Float64}
        @. rho_arr[2:end-1, 2:end-1]   = rho_arr[2:end-1, 2:end-1]   + w * fl_rho   * flow_cells
        @. rho_u_arr[2:end-1, 2:end-1] = rho_u_arr[2:end-1, 2:end-1] + w * fl_rho_u * flow_cells
        @. rho_v_arr[2:end-1, 2:end-1] = rho_v_arr[2:end-1, 2:end-1] + w * fl_rho_v * flow_cells
        @. rho_E_arr[2:end-1, 2:end-1] = rho_E_arr[2:end-1, 2:end-1] + w * fl_rho_E * flow_cells

        if !chem_eq
            gs_arr = var["gamma_star"]::Matrix{Float64}
            cv_arr = var["cv_star"]::Matrix{Float64}
            fl_gs  = flux["gamma_star"]::Matrix{Float64}
            fl_cv  = flux["cv_star"]::Matrix{Float64}
            @. gs_arr[2:end-1, 2:end-1] = gs_arr[2:end-1, 2:end-1] + w * fl_gs * flow_cells
            @. cv_arr[2:end-1, 2:end-1] = cv_arr[2:end-1, 2:end-1] + w * fl_cv * flow_cells
        end

        # Update shock geometry
        s = MOVE_SHOCK(s, dsolution_dt, w)
    else
        ## Standard integration (no shock fitting)
        @. rho_arr[2:end-1, 2:end-1]   = rho_arr[2:end-1, 2:end-1]   + w * fl_rho
        @. rho_u_arr[2:end-1, 2:end-1] = rho_u_arr[2:end-1, 2:end-1] + w * fl_rho_u
        @. rho_v_arr[2:end-1, 2:end-1] = rho_v_arr[2:end-1, 2:end-1] + w * fl_rho_v
        @. rho_E_arr[2:end-1, 2:end-1] = rho_E_arr[2:end-1, 2:end-1] + w * fl_rho_E

        if !chem_eq
            gs_arr = var["gamma_star"]::Matrix{Float64}
            cv_arr = var["cv_star"]::Matrix{Float64}
            fl_gs  = flux["gamma_star"]::Matrix{Float64}
            fl_cv  = flux["cv_star"]::Matrix{Float64}
            @. gs_arr[2:end-1, 2:end-1] = gs_arr[2:end-1, 2:end-1] + w * fl_gs
            @. cv_arr[2:end-1, 2:end-1] = cv_arr[2:end-1, 2:end-1] + w * fl_cv
        end
    end

    return s
end
