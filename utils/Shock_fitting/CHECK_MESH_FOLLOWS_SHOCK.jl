# CHECK_MESH_FOLLOWS_SHOCK - Verify that mesh rows are aligned with the shock.
#
#   s = CHECK_MESH_FOLLOWS_SHOCK(s, chemistry)
#
#   Checks whether the shocked-cell index is uniform across all streamwise
#   stations. If any station has a different shocked-cell index than the
#   first, the mesh is considered misaligned with the shock and the
#   s is restarted on a new mesh via RESTART_SOLUTION.
#
#   Inputs:
#       s  (Dict) - Solution structure containing at minimum:
#                            s["mesh"]["Nchi"]             - Number of streamwise cells
#                            s["shock"]["cell_indices"]    - (Nx x 1) indices of shocked cells
#       chemistry (Dict) - Chemistry/thermodynamic model structure passed
#                            to RESTART_SOLUTION if remeshing is required.
#
#   Outputs:
#       s  (Dict) - Unchanged if mesh is aligned; otherwise the
#                            restarted s on a corrected mesh.
#
#   Notes:
#       - A uniform shocked_cell_indices column is required for proper
#         stability analysis on a shock-aligned grid.
#       - When remeshing is triggered, disturbances are disabled (set to false).
#
#   See also: RESTART_SOLUTION
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

function CHECK_MESH_FOLLOWS_SHOCK(s::Dict{String, Any}, chemistry::Dict{String, Any})
    shock = s["shock"]::Dict{String, Any}

    if shock["enabled"]
        mesh = s["mesh"]::Dict{String, Any}
        cell_indices = shock["cell_indices"]
        Nchi = mesh["Nchi"]::Int
        @inbounds ref_idx = cell_indices[1, 1]

        @inbounds for i in 1:Nchi
            if cell_indices[i, 1] != ref_idx
                println("\n Mesh not aligned with shock. Remesh for proper stability analysis. \n")
                solution_old = s
                disturbances = false
                s = RESTART_SOLUTION(s, solution_old, chemistry, disturbances)
                return s
            end
        end
        println("\n Mesh aligned, not remeshing. \n")
    end

    return s
end
