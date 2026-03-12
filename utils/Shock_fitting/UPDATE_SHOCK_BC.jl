# UPDATE_SHOCK_BC - Apply shock boundary conditions via Rankine-Hugoniot jump relations.
#
#   s = UPDATE_SHOCK_BC(s, chemistry)
#
#   When shock fitting is active, this function performs the full shock
#   boundary-condition update sequence:
#     1. Fit/resample the shock points (LEAST_SQUARES_SHOCK_POINTS).
#     2. Compute the local shock angle beta (COMPUTE_BETA).
#     3. Enforce symmetry at the last station (beta = pi/2).
#     4. Update shocked-cell values from Rankine-Hugoniot relations
#        (UPDATE_SHOCK_JUMP_PROPERTIES).
#     5. Extrapolate into ghost cells (EXTRAPOLATE_CELLS_SHOCK).
#
#   Inputs:
#       s  (Dict) - Solution structure with s["shock"]["enabled"] flag and flow fields.
#       chemistry (Dict) - Chemistry model structure.
#
#   Outputs:
#       s  (Dict) - Updated s with shock boundary conditions applied.
#
#   See also: LEAST_SQUARES_SHOCK_POINTS, COMPUTE_BETA,
#             UPDATE_SHOCK_JUMP_PROPERTIES, EXTRAPOLATE_CELLS_SHOCK
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

function UPDATE_SHOCK_BC(s::Dict{String, Any}, chemistry::Dict{String, Any})
    shock = s["shock"]::Dict{String, Any}
    if shock["enabled"]
        s = LEAST_SQUARES_SHOCK_POINTS(s)
        s = COMPUTE_BETA(s)
        bc = s["boundary_conditions"]::Dict{String, Any}
        if bc["boundary_chi0"]["name"] == "symmetry"
            shock["beta"][1, 1] = pi / 2  # Symmetry condition
        end
        s = UPDATE_SHOCK_JUMP_PROPERTIES(s, chemistry)
        s = EXTRAPOLATE_CELLS_SHOCK(s)
    end

    return s
end
