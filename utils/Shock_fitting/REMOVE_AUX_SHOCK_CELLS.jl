# REMOVE_AUX_SHOCK_CELLS - Reset auxiliary upstream cells to freestream values.
#
#   s = REMOVE_AUX_SHOCK_CELLS(s)
#
#   Replaces all field values in cells upstream of the shocked cell
#   (including the shocked cell itself) with uniform freestream
#   conditions. This cleans up auxiliary/ghost cells that were
#   populated during the shock-fitting extrapolation step.
#
#   Inputs:
#       s (Dict) - Solution structure containing at minimum:
#                           s["shock"]["enabled"]         - (Bool) shock-fitting flag
#                           s["mesh"]["Nchi"]             - Number of streamwise cells
#                           s["shock"]["cell_indices"]    - (Nx x 1) shocked-cell indices
#                           s["freestream"]["rho_0"], s["freestream"]["rho_u_0"],
#                           s["freestream"]["rho_v_0"], s["freestream"]["rho_E_0"],
#                           s["freestream"]["p_0"]        - Scalar freestream conservative values
#                           s["freestream"]["gamma_star"]
#                           s["freestream"]["T"]          - Freestream thermodynamic properties
#
#   Outputs:
#       s (Dict) - Solution with upstream cells set to freestream values.
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

function REMOVE_AUX_SHOCK_CELLS(s::Dict{String, Any})
    shock = s["shock"]::Dict{String, Any}

    if shock["enabled"]
        ## Extract dict refs
        var   = s["var"]::Dict{String, Any}
        fs    = s["freestream"]::Dict{String, Any}
        mesh  = s["mesh"]::Dict{String, Any}
        Nchi  = mesh["Nchi"]::Int

        rho_arr       = var["rho"]::Matrix{Float64}
        rho_u_arr     = var["rho_u"]::Matrix{Float64}
        rho_v_arr     = var["rho_v"]::Matrix{Float64}
        rho_E_arr     = var["rho_E"]::Matrix{Float64}
        p_arr         = var["p"]::Matrix{Float64}
        gamma_star_arr = var["gamma_star"]::Matrix{Float64}
        cv_star_arr   = var["cv_star"]::Matrix{Float64}
        T_arr         = var["T"]::Matrix{Float64}
        cell_indices  = shock["cell_indices"]

        rho_0       = fs["rho_0"]::Float64
        rho_u_0     = fs["rho_u_0"]::Float64
        rho_v_0     = fs["rho_v_0"]::Float64
        rho_E_0     = fs["rho_E_0"]::Float64
        p_0         = fs["p_0"]::Float64
        gamma_star_fs = fs["gamma_star"]::Float64
        T_fs        = fs["T"]::Float64

        @inbounds for i in 1:Nchi
            sc_i = cell_indices[i, 1]
            @views rho_arr[i + 1, sc_i + 1:end]        .= rho_0
            @views rho_u_arr[i + 1, sc_i + 1:end]      .= rho_u_0
            @views rho_v_arr[i + 1, sc_i + 1:end]      .= rho_v_0
            @views rho_E_arr[i + 1, sc_i + 1:end]      .= rho_E_0
            @views p_arr[i + 1, sc_i + 1:end]          .= p_0
            @views gamma_star_arr[i + 1, sc_i + 1:end] .= gamma_star_fs
            @views cv_star_arr[i + 1, sc_i + 1:end]    .= 1  # Non-dimensional cv_star is 1 for freestream
            @views T_arr[i + 1, sc_i + 1:end]          .= T_fs
        end
    end

    return s
end
