# INITIAL_SOLUTION_POST_SHOCK - Initialize post-shock flow field using Rankine-Hugoniot relations.
#
#   s = INITIAL_SOLUTION_POST_SHOCK(s, chemistry)
#
#   Computes the initial post-shock conservative variables (rho, rho_u,
#   rho_v, rho_E) at every streamwise station by applying oblique-shock
#   Rankine-Hugoniot jump conditions. Supports both calorically perfect
#   gas and finite-rate / equilibrium chemistry models.
#
#   Inputs:
#       s  (Dict) - Solution structure containing at minimum:
#                            s["mesh"]["Nchi"], s["mesh"]["Neta"]  - Grid dimensions
#                            s["freestream"]["Mach"]        - Freestream Mach number
#                            s["freestream"]["gamma"]       - Ratio of specific heats (perfect gas)
#                            s["shock"]["beta"]             - (Nx x 1) local shock angles [rad]
#                            s["shock"]["cell_indices"]     - (Nx x 1) shocked-cell column indices
#                            s["chemistry"]["is_chemistry_enabled"] - (Bool) chemistry enabled flag
#                            s["chemistry"]["chemical_equilibrium"] - (Bool) equilibrium chemistry flag
#                            s["freestream"]                - Freestream reference quantities
#       chemistry (Dict) - Chemistry model structure; used when
#                            chemistry_state is true. Contains either
#                            equilibrium or frozen sub-models with
#                            SOLVE_RANKINE_HUGONIOT_CHEMISTRY and
#                            eval_gamma_star / eval_cv_star methods.
#
#   Outputs:
#       s  (Dict) - Updated s with post-shock fields populated
#                            from the body surface up to the shock location.
#
#   Notes:
#       - For perfect gas, the standard oblique-shock relations are used.
#         If M*sin(beta) <= 1 at a station, that station is treated as
#         unshocked and set to freestream conditions.
#       - For reacting flows, SOLVE_RANKINE_HUGONIOT_CHEMISTRY is called
#         to obtain the post-shock state, then results are non-dimensionalized.
#       - After the main loop, gamma_star and cv_star are evaluated for
#         frozen chemistry cases.
#
#   See also: SOLVE_RANKINE_HUGONIOT_CHEMISTRY, COMPUTE_BETA
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

function INITIAL_SOLUTION_POST_SHOCK(s::Dict{String, Any}, chemistry::Dict{String, Any})
    ## Extract dict refs
    var   = s["var"]::Dict{String, Any}
    fs    = s["freestream"]::Dict{String, Any}
    chem  = s["chemistry"]::Dict{String, Any}
    shock = s["shock"]::Dict{String, Any}
    mesh  = s["mesh"]::Dict{String, Any}
    Nchi  = mesh["Nchi"]::Int
    Neta  = mesh["Neta"]::Int

    rho_arr   = var["rho"]::Matrix{Float64}
    rho_u_arr = var["rho_u"]::Matrix{Float64}
    rho_v_arr = var["rho_v"]::Matrix{Float64}
    rho_E_arr = var["rho_E"]::Matrix{Float64}
    cell_indices = shock["cell_indices"]
    beta_arr     = shock["beta"]

    is_chem = chem["is_chemistry_enabled"]::Bool
    chem_eq = chem["chemical_equilibrium"]::Bool

    if !is_chem
        ## Perfect gas oblique-shock initialization
        M     = fs["Mach"]
        gamma = fs["gamma"]
        rho_u_0_fs = fs["rho_u_0"]
        rho_v_0_fs = fs["rho_v_0"]

        for i in 1:Nchi
            beta = beta_arr[i, 1]

            if (M * sin(beta) > 1)
                ang = beta - atan(rho_v_0_fs, rho_u_0_fs)
                temp = (gamma + 1) * (M * sin(beta))^2 /
                     ((gamma - 1) * (M * sin(beta))^2 + 2)
                @views rho_arr[i + 1, 1:cell_indices[i, 1] + 2] .= temp

                normal_velocity = sin(ang) / temp
                tangential_velocity = cos(ang)
                @views rho_u_arr[i + 1, 1:cell_indices[i, 1] + 2] .=
                    temp * (normal_velocity * sin(beta) + tangential_velocity * cos(beta))
                @views rho_v_arr[i + 1, 1:cell_indices[i, 1] + 2] .=
                    temp * (-normal_velocity * cos(beta) + tangential_velocity * sin(beta))

                # Total enthalpy is preserved across the shock
                p_inf = 1 / gamma / M^2
                k = (rho_u_arr[i + 1, 1]^2 + rho_v_arr[i + 1, 1]^2) / 2 / temp^2
                p = p_inf * (2 * gamma * (M * sin(beta))^2 - (gamma - 1)) /
                  (gamma + 1)
                @views rho_E_arr[i + 1, 1:cell_indices[i, 1] + 2] .=
                    p / (gamma - 1) + k * temp
            else
                # No shock at this station: set to freestream
                @views rho_arr[i + 1, 1:cell_indices[i, 1] + 2] .= 1
                @views rho_u_arr[i + 1, 1:cell_indices[i, 1] + 2] .=
                    rho_u_arr[Nchi + 1, Neta + 1]
                @views rho_v_arr[i + 1, 1:cell_indices[i, 1] + 2] .=
                    rho_v_arr[Nchi + 1, Neta + 1]
                @views rho_E_arr[i + 1, 1:cell_indices[i, 1] + 2] .=
                    (1/2 + 1 / (gamma - 1) / gamma / M^2)
            end
        end
    else
        ## Reacting-gas oblique-shock initialization
        rho_u_0_fs = fs["rho_u_0"]
        rho_v_0_fs = fs["rho_v_0"]
        U_fs       = fs["U"]
        e_fs       = fs["e"]
        rho_fs     = fs["rho"]
        rho_factor    = fs["rho_factor"]::Float64
        energy_factor = fs["energy_factor"]::Float64
        vel_factor    = fs["velocity_factor"]::Float64

        for i in 1:Nchi
            beta = beta_arr[i, 1]
            ang = beta - atan(rho_v_0_fs, rho_u_0_fs)
            w_1 = sin(ang) * U_fs
            e_1 = e_fs
            rho_1 = rho_fs

            if chem_eq
                rho_2, e_2, w_2 = SOLVE_RANKINE_HUGONIOT_CHEMISTRY(chemistry, w_1, rho_1, e_1)
            else
                rho_2, e_2, w_2 = SOLVE_RANKINE_HUGONIOT_CHEMISTRY(chemistry["frozen"], w_1, rho_1, e_1)
            end
            # Extract scalars from 1-element vectors returned by the solver
            rho_2 = rho_2[1]; e_2 = e_2[1]; w_2 = w_2[1]

            # Non-dimensionalize post-shock state
            rho_2_nd = rho_2 / rho_factor
            e_2_nd   = e_2   / energy_factor
            w_2_nd   = w_2   / vel_factor

            normal_velocity     = w_2_nd
            tangential_velocity = cos(ang) * U_fs / vel_factor

            @views rho_arr[i + 1, 1:cell_indices[i, 1] + 2] .= rho_2_nd
            @views rho_u_arr[i + 1, 1:cell_indices[i, 1] + 2] .=
                rho_2_nd * (normal_velocity * sin(beta) + tangential_velocity * cos(beta))
            @views rho_v_arr[i + 1, 1:cell_indices[i, 1] + 2] .=
                rho_2_nd * (-normal_velocity * cos(beta) + tangential_velocity * sin(beta))

            k = (normal_velocity^2 + tangential_velocity^2) / 2
            @views rho_E_arr[i + 1, 1:cell_indices[i, 1] + 2] .=
                rho_2_nd * e_2_nd + k * rho_2_nd
        end
    end

    ## Evaluate thermodynamic properties for frozen chemistry
    if is_chem
        if !chem_eq
            rho_factor    = fs["rho_factor"]::Float64
            energy_factor = fs["energy_factor"]::Float64
            cv_fs         = fs["cv"]::Float64

            e = @. (rho_E_arr - (rho_u_arr^2 + rho_v_arr^2) / rho_arr / 2) / rho_arr
            var["gamma_star"] = chemistry["frozen"]["eval_gamma_star"](
                rho_factor * rho_arr, energy_factor * e)
            var["cv_star"] = chemistry["frozen"]["eval_cv_star"](
                rho_factor * rho_arr, energy_factor * e) ./ cv_fs
        end
    end

    return s
end
