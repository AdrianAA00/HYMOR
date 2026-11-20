# RESTART_SOLUTION  Restart a solution from a previously saved state.
#   Handles re-meshing with field interpolation when the grid dimensions
#   have changed, or directly restores the old solution when no re-mesh
#   is required.
#
#   s = RESTART_SOLUTION(s, solution_old, chemistry, disturbances)
#
#   Inputs:
#       s            - (Dict) Current solution struct (target grid settings)
#       solution_old - (Dict) Previously saved solution struct to restart from
#       chemistry    - (Dict) Chemistry model struct for thermodynamic evaluations
#       disturbances - Disturbance settings passed to GENERATE_MESH
#
#   Outputs:
#       s            - (Dict) Updated solution struct with restarted fields
#
#   Notes:
#       - If Nchi or Neta differ between old and new solutions, remeshing is
#         triggered and all conserved fields are interpolated onto the new
#         grid (REMESH_INTERPOLATION: cubic interpolation along the old grid
#         lines, smooth on stretched and curved meshes).
#       - When a shock is present, only post-shock data is used for
#         interpolation to avoid smearing the discontinuity.
#       - A remesh moves the shock relative to the cells (it is re-centred
#         in its row). The interpolated state is not a steady state of the
#         new mesh: march further before linearizing about it.
#       - After interpolation, boundary conditions, chemistry, and
#         thermodynamic properties are recomputed on the new mesh.
#       - If no remeshing is needed, the old solution is copied directly.
#
# Part of: Hypersonics Stability MATLAB Solver - Initialization Module

function RESTART_SOLUTION(s::Dict{String,Any}, solution_old::Dict{String,Any},
                          chemistry::Dict{String,Any}, disturbances)

    ## Ensure solution_old is an independent copy to prevent aliasing issues
    ## (Julia Dicts are mutable reference types, unlike MATLAB structs which are value types)
    if s === solution_old
        solution_old = deepcopy(solution_old)
    end

    ## Check if remeshing is required
    if s["mesh"]["Nchi"] != solution_old["mesh"]["Nchi"] || s["mesh"]["Neta"] != solution_old["mesh"]["Neta"]
        s["remesh"] = true
    end

    if s["restart"] && s["remesh"]

        ## Initialize variables on new mesh
        s = VARIABLES_INITIALIZATION(s)
        s["linearize"] = false

        ## Precompute chi_wall for shock interpolation (same as GENERATE_MESH)
        chi_wall_temp = collect(LinRange(0, 1, s["mesh"]["Nchi"] + 1))
        if s["curvilinear_mapping"]["refinement_stagnation"]["state"]
            chi_wall_temp = REFINEMENT(chi_wall_temp,
                s["curvilinear_mapping"]["refinement_stagnation"]["BL_thickness"],
                0, s["curvilinear_mapping"]["refinement_stagnation"]["intensity"], 1, 0)
        end
        s["mesh"]["chi_wall"] = chi_wall_temp

        ## Interpolate shock position to new mesh
        if s["shock"]["enabled"] == true
            s["shock"]["points_chi"] = (s["mesh"]["chi_wall"][1:end-1] .+ s["mesh"]["chi_wall"][2:end]) ./ 2
            shock_points_chi_temp = (s["mesh"]["chi_wall"][1:end-1] .+ s["mesh"]["chi_wall"][2:end]) ./ 2
            s["shock"]["points_x"] = ppval(solution_old["shock"]["spline_func_x"], shock_points_chi_temp)
            s["shock"]["points_y"] = ppval(solution_old["shock"]["spline_func_y"], shock_points_chi_temp)
            s["shock"]["spline_func_x"] = solution_old["shock"]["spline_func_x"]
            s["shock"]["spline_func_y"] = solution_old["shock"]["spline_func_y"]
            (s["shock"]["points_chi"], s["shock"]["points_eta"]) = GO_TO_ELEMENT_SPACE(
                s["shock"]["points_x"], s["shock"]["points_y"], s)
        end

        ## Carry forward time and generate new mesh
        s["time_integration"]["t"] = solution_old["time_integration"]["t"]
        s = GENERATE_MESH(s, disturbances)

        ## Interpolate fields onto the new mesh (cubic, in the old grid index space)
        s = REMESH_INTERPOLATION(s, solution_old)

        ## Update freestream and shock boundary conditions
        s = SET_FREESTREAM_PROPERTIES(s, chemistry)

        if s["shock"]["enabled"] == true
            s = UPDATE_FLOW_CELLS(s, chemistry)
            s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry)
            s = UPDATE_SOUND_SPEED(s, chemistry)
            s = UPDATE_SHOCK_BC(s, chemistry)
            s = SET_FREESTREAM_PROPERTIES(s, chemistry)
            s = UPDATE_FIELD_UPSTREAM(s)
        end

        ## Extend upstream flow and finalize boundary conditions
        s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry)
        s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry)
        s = APPLY_BOUNDARY_CONDITIONS(s, chemistry)
        s = UPDATE_SOUND_SPEED(s, chemistry)

    elseif s["restart"] && !s["remesh"]
        ## No remesh: copy old solution directly
        s = deepcopy(solution_old)
    end

    return s
end
