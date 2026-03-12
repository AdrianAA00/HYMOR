# UPDATE_FLOW_CELLS - Detect shocked cells and update thermodynamic properties.
#
#   s = UPDATE_FLOW_CELLS(s, chemistry)
#
#   When shock fitting is active, this function orchestrates the sequence
#   of operations needed after the flow field has been advanced in time:
#     1. Detect which cells are crossed by the shock.
#     2. Extrapolate flow quantities into the ghost cells beyond the shock.
#     3. Update chemical equilibrium composition (if applicable).
#     4. Recompute thermodynamic properties (pressure, sound speed, etc.).
#
#   Inputs:
#       s  (Dict) - Solution structure with s["shock"] flag and flow fields.
#       chemistry (Dict) - Chemistry model structure.
#
#   Outputs:
#       s  (Dict) - Updated s with consistent cell flags,
#                            ghost-cell values, and thermodynamic properties.
#
#   See also: DETECT_CELLS_SHOCKED, EXTRAPOLATE_CELLS_SHOCK,
#             UPDATE_CHEMISTRY_EQUILIBRIUM, UPDATE_THERMODYNAMIC_PROPERTIES
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

function UPDATE_FLOW_CELLS(s::Dict{String, Any}, chemistry::Dict{String, Any})
    shock = s["shock"]::Dict{String, Any}
    if shock["enabled"]
        s = DETECT_CELLS_SHOCKED(s)
        s = EXTRAPOLATE_CELLS_SHOCK(s)
        s = UPDATE_FIELD_UPSTREAM(s)
        s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry)
        s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry)
        s = UPDATE_SOUND_SPEED(s, chemistry)
    end

    return s
end
