# SET_FREESTREAM_PROPERTIES  Compute freestream thermodynamic and transport properties.
#   Evaluates velocity magnitude, internal energy, sound speed, specific
#   heat ratio, viscosity, thermal conductivity, and pressure at freestream
#   conditions.  Also sets the non-dimensional upstream conserved-variable
#   reference values used by the solver.
#
#   s = SET_FREESTREAM_PROPERTIES(s, chemistry)
#
#   Inputs:
#       s         - (Dict) Solution struct containing freestream state.
#       chemistry - (Dict) Chemistry model struct providing evaluation functions
#                   when chemistry is enabled.
#
#   Outputs:
#       s         - (Dict) Updated struct with computed freestream properties.
#
#   Notes:
#       - When s["chemistry"]["is_chemistry_enabled"] is true, thermodynamic
#         properties are computed from the chemistry model; otherwise,
#         ideal-gas relations based on Mach number and gamma are used.
#       - Non-dimensionalization uses freestream density and velocity
#         magnitude as reference scales.
#
# Part of: Hypersonics Stability MATLAB Solver - Initialization Module

function SET_FREESTREAM_PROPERTIES(s::Dict{String,Any}, chemistry::Dict{String,Any})

    ## Velocity magnitude and kinetic energy
    s["freestream"]["U"] = sqrt(s["freestream"]["u"]^2 + s["freestream"]["v"]^2)
    s["freestream"]["kin"] = 0.5 * s["freestream"]["U"]^2

    ## Dimensional reference factors for non-dimensionalization (based on rho, U_mag, L)
    s["freestream"]["velocity_factor"] = s["freestream"]["U"]
    s["freestream"]["energy_factor"]   = s["freestream"]["U"]^2
    s["freestream"]["rho_factor"]      = s["freestream"]["rho"]

    ## Internal energy [J/kg]
    if s["chemistry"]["is_chemistry_enabled"]
        s["freestream"]["e"] = chemistry["eval_e"](s["freestream"]["T"], s["freestream"]["rho"])
    else
        s["freestream"]["e"] = s["freestream"]["energy_factor"] *
            1 / (s["freestream"]["gamma"] - 1) / s["freestream"]["gamma"] / s["freestream"]["Mach"]^2
    end

    ## Sound speed [m/s]
    if s["chemistry"]["is_chemistry_enabled"]
        s["freestream"]["a"] = chemistry["eval_a"](s["freestream"]["rho"], s["freestream"]["e"])
    else
        s["freestream"]["a"] = s["freestream"]["U"] / s["freestream"]["Mach"]
    end

    ## Total energy [J/kg]
    s["freestream"]["E"] = s["freestream"]["e"] + s["freestream"]["kin"]

    ## Effective specific heat ratio gamma* [non-dimensional]
    if s["chemistry"]["is_chemistry_enabled"]
        s["freestream"]["gamma_star"] = chemistry["eval_gamma_star"](
            s["freestream"]["rho"], s["freestream"]["e"])
    else
        s["freestream"]["gamma_star"] = s["freestream"]["gamma"]
    end

    ## Effective specific heat cv [J/(K*kg)]
    if s["chemistry"]["is_chemistry_enabled"]
        s["freestream"]["cv"] = chemistry["eval_cv_star"](s["freestream"]["rho"], s["freestream"]["e"])
    else
        s["freestream"]["cv"] = s["freestream"]["e"] / s["freestream"]["T"]
    end

    ## Dynamic viscosity [kg/(m*s)]
    # Set length factor to satisfy Reynolds number definition: Re = rho * U * L / mu
    if haskey(s["freestream"], "Re")
        if s["chemistry"]["is_chemistry_enabled"]
            s["freestream"]["mu"] = chemistry["eval_mu"](s["freestream"]["rho"], s["freestream"]["e"])
            s["freestream"]["L"] = s["freestream"]["mu"] * s["freestream"]["Re"] / s["freestream"]["rho"] / s["freestream"]["U"]
        else
            s["freestream"]["L"] = 1.0
            s["freestream"]["mu"] = s["freestream"]["rho"] * s["freestream"]["U"] * s["freestream"]["L"] / s["freestream"]["Re"]
        end
    elseif haskey(s["freestream"], "L_ref")
        if s["chemistry"]["is_chemistry_enabled"]
            s["freestream"]["mu"] = chemistry["eval_mu"](s["freestream"]["rho"], s["freestream"]["e"])
            s["freestream"]["Re"] = s["freestream"]["rho"] * s["freestream"]["U"] * s["freestream"]["L_ref"] / s["freestream"]["mu"]
            s["freestream"]["L"] = s["freestream"]["L_ref"]
        else
            error("Reynolds number not specified. \n")
        end
    end

    ## Thermal conductivity [W/(m*K)]
    if s["chemistry"]["is_chemistry_enabled"]
        s["freestream"]["k"] = chemistry["eval_k"](s["freestream"]["rho"], s["freestream"]["e"])
        s["freestream"]["Pr"] = s["freestream"]["mu"] * s["freestream"]["cv"] * s["freestream"]["gamma_star"] / s["freestream"]["k"]
    else
        s["freestream"]["k"] = s["freestream"]["gamma"] * s["freestream"]["cv"] *
            s["freestream"]["mu"] / s["freestream"]["Pr"]
    end

    ## Freestream pressure [Pa]
    s["freestream"]["p"] = s["freestream"]["rho"] * s["freestream"]["e"] *
        (s["freestream"]["gamma_star"] - 1)

    ## Non-dimensional upstream reference state
    s["freestream"]["rho_u_0"] = s["freestream"]["rho"] * s["freestream"]["u"] / (s["freestream"]["U"] * s["freestream"]["rho_factor"])
    s["freestream"]["rho_v_0"] = s["freestream"]["rho"] * s["freestream"]["v"] / (s["freestream"]["U"] * s["freestream"]["rho_factor"])
    s["freestream"]["rho_0"]   = s["freestream"]["rho"] / s["freestream"]["rho_factor"]
    s["freestream"]["rho_E_0"] = s["freestream"]["rho"] * s["freestream"]["E"] /
        (s["freestream"]["energy_factor"] * s["freestream"]["rho_factor"])
    s["freestream"]["p_0"]     = s["freestream"]["p"] /
        (s["freestream"]["rho_factor"] * s["freestream"]["velocity_factor"]^2)
    s["freestream"]["a_0"] = s["freestream"]["a"] / s["freestream"]["velocity_factor"]

    return s
end
