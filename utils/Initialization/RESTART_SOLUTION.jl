using NaturalNeighbours

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
#         grid using linear interpolation.
#       - When a shock is present, only data from inner (sub-shock) cells
#         is used for interpolation to avoid smearing the discontinuity.
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

        ## Interpolate fields without shock (full-domain Delaunay linear interpolation)
        if !s["shock"]["enabled"]
            x_old = vec(solution_old["mesh"]["x_Ext"])
            y_old = vec(solution_old["mesh"]["y_Ext"])

            x_new = vec(s["mesh"]["x_Ext"])
            y_new = vec(s["mesh"]["y_Ext"])

            itp_rho = NaturalNeighbours.interpolate(x_old, y_old, vec(solution_old["var"]["rho"]))
            s["var"]["rho"] = reshape(itp_rho(x_new, y_new; method=Triangle()), size(s["mesh"]["x_Ext"]))

            itp_rho_u = NaturalNeighbours.interpolate(x_old, y_old, vec(solution_old["var"]["rho_u"]))
            s["var"]["rho_u"] = reshape(itp_rho_u(x_new, y_new; method=Triangle()), size(s["mesh"]["x_Ext"]))

            itp_rho_v = NaturalNeighbours.interpolate(x_old, y_old, vec(solution_old["var"]["rho_v"]))
            s["var"]["rho_v"] = reshape(itp_rho_v(x_new, y_new; method=Triangle()), size(s["mesh"]["x_Ext"]))

            itp_rho_E = NaturalNeighbours.interpolate(x_old, y_old, vec(solution_old["var"]["rho_E"]))
            s["var"]["rho_E"] = reshape(itp_rho_E(x_new, y_new; method=Triangle()), size(s["mesh"]["x_Ext"]))

            if !s["chemistry"]["chemical_equilibrium"]
                itp_gamma = NaturalNeighbours.interpolate(x_old, y_old, vec(solution_old["var"]["gamma_star"]))
                s["var"]["gamma_star"] = reshape(itp_gamma(x_new, y_new; method=Triangle()), size(s["mesh"]["x_Ext"]))

                itp_cv = NaturalNeighbours.interpolate(x_old, y_old, vec(solution_old["var"]["cv_star"]))
                s["var"]["cv_star"] = reshape(itp_cv(x_new, y_new; method=Triangle()), size(s["mesh"]["x_Ext"]))
            end
        end

        ## Interpolate fields with shock (sub-shock cells only)
        # Store only data from inner shock cells to avoid building an
        # interpolant across the discontinuity.
        if s["shock"]["enabled"]
            x_temp = Float64[]
            y_temp = Float64[]
            temp_rho = Float64[]
            temp_rho_u = Float64[]
            temp_rho_v = Float64[]
            temp_rho_E = Float64[]
            temp_gamma_star = Float64[]
            temp_cv_star = Float64[]

            for i in 1:(solution_old["mesh"]["Nchi"] + 2)
                if i > 1 && i < solution_old["mesh"]["Nchi"] + 2
                    jmax = solution_old["shock"]["cell_indices"][i - 1, 1] + 2
                elseif i == 1
                    jmax = solution_old["shock"]["cell_indices"][1, 1] + 2
                else
                    jmax = solution_old["shock"]["cell_indices"][solution_old["mesh"]["Nchi"], 1] + 2
                end

                for j in 1:Int(jmax)
                    push!(x_temp, solution_old["mesh"]["x_Ext"][i, j])
                    push!(y_temp, solution_old["mesh"]["y_Ext"][i, j])
                    push!(temp_rho, solution_old["var"]["rho"][i, j])
                    push!(temp_rho_u, solution_old["var"]["rho_u"][i, j])
                    push!(temp_rho_v, solution_old["var"]["rho_v"][i, j])
                    push!(temp_rho_E, solution_old["var"]["rho_E"][i, j])
                    if !s["chemistry"]["chemical_equilibrium"]
                        push!(temp_gamma_star, solution_old["var"]["gamma_star"][i, j])
                        push!(temp_cv_star, solution_old["var"]["cv_star"][i, j])
                    end
                end
            end

            ## Fit and evaluate interpolants on new mesh (Delaunay linear)
            x_new = vec(s["mesh"]["x_Ext"])
            y_new = vec(s["mesh"]["y_Ext"])

            itp_rho = NaturalNeighbours.interpolate(x_temp, y_temp, temp_rho)
            s["var"]["rho"] = reshape(itp_rho(x_new, y_new; method=Triangle()), size(s["mesh"]["x_Ext"]))

            itp_rho_u = NaturalNeighbours.interpolate(x_temp, y_temp, temp_rho_u)
            s["var"]["rho_u"] = reshape(itp_rho_u(x_new, y_new; method=Triangle()), size(s["mesh"]["x_Ext"]))

            itp_rho_v = NaturalNeighbours.interpolate(x_temp, y_temp, temp_rho_v)
            s["var"]["rho_v"] = reshape(itp_rho_v(x_new, y_new; method=Triangle()), size(s["mesh"]["x_Ext"]))

            itp_rho_E = NaturalNeighbours.interpolate(x_temp, y_temp, temp_rho_E)
            s["var"]["rho_E"] = reshape(itp_rho_E(x_new, y_new; method=Triangle()), size(s["mesh"]["x_Ext"]))

            if !s["chemistry"]["chemical_equilibrium"]
                itp_gamma = NaturalNeighbours.interpolate(x_temp, y_temp, temp_gamma_star)
                s["var"]["gamma_star"] = reshape(itp_gamma(x_new, y_new; method=Triangle()), size(s["mesh"]["x_Ext"]))

                itp_cv = NaturalNeighbours.interpolate(x_temp, y_temp, temp_cv_star)
                s["var"]["cv_star"] = reshape(itp_cv(x_new, y_new; method=Triangle()), size(s["mesh"]["x_Ext"]))
            end
        end

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
