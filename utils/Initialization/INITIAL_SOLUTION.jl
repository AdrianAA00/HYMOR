# INITIAL_SOLUTION  Set initial flow field to uniform freestream conditions.
#
#   s = INITIAL_SOLUTION(s, chemistry)
#
#   Initializes all conservative flow variables (density, momentum, and
#   total energy) to the upstream freestream values across the entire
#   computational domain. After setting the uniform field, the chemistry
#   equilibrium state is updated. For frozen (non-equilibrium) chemistry,
#   gamma_star and cv_star are also evaluated from the thermodynamic tables.
#
#   Inputs:
#       s         - (Dict) Solver data structure.
#       chemistry - (Dict) Chemistry model with evaluation functions.
#
#   Outputs:
#       s         - (Dict) Updated solver structure with uniform initial
#                   flow field and chemistry properties.
#
#   See also: SET_FREESTREAM_PROPERTIES, UPDATE_CHEMISTRY_EQUILIBRIUM
#
#   Part of: Hypersonics Stability MATLAB Solver - Initialization Module

function INITIAL_SOLUTION(s::Dict{String,Any}, chemistry::Dict{String,Any})

    ## Set freestream reference properties
    s = SET_FREESTREAM_PROPERTIES(s, chemistry)

    ## Initialize conservative variables to uniform upstream state
    field_size = size(s["var"]["rho_E"])

    boundary_type = s["curvilinear_mapping"]["boundary_type"]

    if boundary_type == "MSL"           ||
       boundary_type == "blunt_cone"    ||
       boundary_type == "circle"        ||
       boundary_type == "lid_driven_cavity"

        s["var"]["rho_u"] = fill(s["freestream"]["rho_u_0"], field_size)
        s["var"]["rho_v"] = fill(s["freestream"]["rho_v_0"], field_size)
        s["var"]["rho"]   = fill(s["freestream"]["rho_0"],   field_size)
        s["var"]["rho_E"] = fill(s["freestream"]["rho_E_0"], field_size)

    elseif boundary_type == "channel"
        # Parabolic velocity profile: u(y) = u_max * y*(2-y)
        #   zero at y = 0 and y = 2, maximum at y = 1
        y = 2 .* s["mesh"]["y_Ext"] ./ s["curvilinear_mapping"]["Ly"]  # Normalize to [0,2] for the profile
        u_profile = y .* (2 .- y)  # ranges from 0 to 1

        s["var"]["rho"]   = fill(s["freestream"]["rho_0"], field_size)
        s["var"]["rho_v"] = zeros(field_size)  # v = 0
        s["var"]["rho_u"] = s["freestream"]["rho_u_0"] .* u_profile

        # Internal energy from freestream, then add local kinetic energy
        rho_e_internal = s["freestream"]["rho_E_0"] -
            0.5 * (s["freestream"]["rho_u_0"]^2 + s["freestream"]["rho_v_0"]^2) / s["freestream"]["rho_0"]
        s["var"]["rho_E"] = rho_e_internal .+ 0.5 .* s["var"]["rho_u"].^2 ./ s["var"]["rho"]

    else
        error("Unsupported boundary type: $(boundary_type)")
    end

    ## Update chemistry equilibrium state
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry)

    ## Evaluate frozen chemistry properties if applicable
    if s["chemistry"]["is_chemistry_enabled"]
        if !s["chemistry"]["chemical_equilibrium"]
            rho = s["var"]["rho"]
            e = (s["var"]["rho_E"] .- (s["var"]["rho_u"].^2 .+ s["var"]["rho_v"].^2) ./ rho ./ 2) ./ rho
            s["var"]["gamma_star"] = chemistry["frozen"]["eval_gamma_star"](
                s["freestream"]["rho_factor"] .* rho,
                s["freestream"]["energy_factor"] .* e)
            s["var"]["cv_star"] = chemistry["frozen"]["eval_cv_star"](
                s["freestream"]["rho_factor"] .* rho,
                s["freestream"]["energy_factor"] .* e) ./ s["freestream"]["cv"]
        end
    end

    return s
end
