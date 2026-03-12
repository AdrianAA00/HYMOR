# UPDATE_FIELD_UPSTREAM - Set upstream field cells to freestream values with perturbations.
#
#   s = UPDATE_FIELD_UPSTREAM(s)
#
#   Fills all grid cells upstream of the shock (beyond the two ghost-cell
#   layers used for extrapolation) with either a prescribed perturbed
#   upstream state or a uniform freestream plus travelling-wave
#   perturbation. Also sets gamma_star and cv_star to freestream values
#   for non-equilibrium chemistry cases.
#
#   Inputs:
#       s (Dict) - Solution structure containing at minimum:
#                    s["mesh"]["Nchi"]              - Number of streamwise cells
#                    s["shock"]["cell_indices"]     - (Nx x 1) shocked-cell indices
#                    s["freestream"]["rho_0"], s["freestream"]["rho_u_0"],
#                    s["freestream"]["rho_v_0"], s["freestream"]["rho_E_0"]
#                                                  - Scalar uniform freestream values
#                    s["freestream"]["disturbance"]["k_y"], s["freestream"]["disturbance"]["k_x"]
#                                                  - Perturbation wavenumbers (tangential, normal)
#                    s["freestream"]["disturbance"]["amplitude"]     - (4 x 1) perturbation amplitudes
#                    s["mesh"]["x_Ext"], s["mesh"]["y_Ext"] - Extended grid coordinates
#                    s["time_integration"]["t"]      - Current simulation time
#                    s["chemistry"]["chemical_equilibrium"] - (Bool) equilibrium chemistry flag
#                    s["freestream"]["gamma_star"]   - Freestream thermodynamic properties
#                    (optionally) s["freestream"]["rho_0_p"], etc. - Prescribed perturbed fields
#
#   Outputs:
#       s (Dict) - Solution with upstream cells populated.
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

function UPDATE_FIELD_UPSTREAM(s::Dict{String, Any})
    ## Extract dict refs
    var   = s["var"]::Dict{String, Any}
    fs    = s["freestream"]::Dict{String, Any}
    mesh  = s["mesh"]::Dict{String, Any}
    shock = s["shock"]::Dict{String, Any}
    chem  = s["chemistry"]::Dict{String, Any}
    ti    = s["time_integration"]::Dict{String, Any}
    Nchi  = mesh["Nchi"]::Int

    rho_arr       = var["rho"]::Matrix{Float64}
    rho_u_arr     = var["rho_u"]::Matrix{Float64}
    rho_v_arr     = var["rho_v"]::Matrix{Float64}
    rho_E_arr     = var["rho_E"]::Matrix{Float64}
    gamma_star_arr = var["gamma_star"]::Matrix{Float64}
    cv_star_arr   = var["cv_star"]::Matrix{Float64}
    cell_indices  = shock["cell_indices"]
    chem_eq       = chem["chemical_equilibrium"]::Bool

    x_Ext = mesh["x_Ext"]::Matrix{Float64}
    y_Ext = mesh["y_Ext"]::Matrix{Float64}
    t     = ti["t"]::Float64

    rho_0   = fs["rho_0"]::Float64
    rho_u_0 = fs["rho_u_0"]::Float64
    rho_v_0 = fs["rho_v_0"]::Float64
    rho_E_0 = fs["rho_E_0"]::Float64
    gamma_star_fs = fs["gamma_star"]::Float64

    disturbance = fs["disturbance"]::Dict{String, Any}
    k_y = disturbance["k_y"]::Float64
    k_x = disturbance["k_x"]::Float64
    amplitude = disturbance["amplitude"]

    u_y = rho_v_0 / rho_0
    u_x = rho_u_0 / rho_0
    two_pi_ky = 2.0 * pi * k_y
    two_pi_kx = 2.0 * pi * k_x
    u_x_t = u_x * t
    u_y_t = u_y * t
    @inbounds perturbation = @. cos(two_pi_ky * (x_Ext - u_x_t)) *
                                sin(two_pi_kx * (y_Ext - u_y_t))

    if haskey(s, "rho_0_upstream_p")
        ## Prescribed perturbed upstream state
        rho_0_p   = fs["rho_0_p"]
        rho_u_0_p = fs["rho_u_0_p"]
        rho_v_0_p = fs["rho_v_0_p"]
        rho_E_0_p = fs["rho_E_0_p"]

        @inbounds for i in 1:Nchi
            sc_i = cell_indices[i, 1]
            @views rho_arr[i + 1, sc_i + 3:end]       .= rho_0_p[i, 1]
            @views rho_u_arr[i + 1, sc_i + 3:end]     .= rho_u_0_p[i, 1]
            @views rho_v_arr[i + 1, sc_i + 3:end]     .= rho_v_0_p[i, 1]
            @views rho_E_arr[i + 1, sc_i + 3:end]     .= rho_E_0_p[i, 1]
            if !chem_eq
                @views gamma_star_arr[i + 1, sc_i + 3:end] .= gamma_star_fs
                @views cv_star_arr[i + 1, sc_i + 3:end]    .= 1
            end
        end
    else
        ## Uniform freestream plus travelling-wave perturbation
        amp1 = amplitude[1]
        amp2 = amplitude[2]
        amp3 = amplitude[3]
        amp4 = amplitude[4]
        @inbounds for i in 1:Nchi
            sc_i = cell_indices[i, 1]
            @views rho_arr[i + 1, sc_i + 3:end]   .= rho_0   .+ perturbation[i + 1, sc_i + 3:end] .* amp1
            @views rho_u_arr[i + 1, sc_i + 3:end] .= rho_u_0 .+ perturbation[i + 1, sc_i + 3:end] .* amp2
            @views rho_v_arr[i + 1, sc_i + 3:end] .= rho_v_0 .+ perturbation[i + 1, sc_i + 3:end] .* amp3
            @views rho_E_arr[i + 1, sc_i + 3:end] .= rho_E_0 .+ perturbation[i + 1, sc_i + 3:end] .* amp4
            if !chem_eq
                @views gamma_star_arr[i + 1, sc_i + 3:end] .= gamma_star_fs
                @views cv_star_arr[i + 1, sc_i + 3:end]    .= 1
            end
        end
    end

    return s
end
