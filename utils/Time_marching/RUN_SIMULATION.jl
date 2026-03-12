using Printf

# RUN_SIMULATION  Execute the main time-marching loop for the flow solver.
#
#   s = RUN_SIMULATION(s, solution_base, chemistry)
#
#   Runs the unsteady flow simulation for s["time_integration"]["N_iter"] time steps.
#   At each iteration the flow is advected, diagnostic information is
#   printed, running plots are updated, and adaptive mesh refinement is
#   triggered when the shock moves too close to the grid boundary or
#   drifts from its last refined position. After the main loop, boundary
#   conditions and thermodynamic properties are finalised, and shock
#   evolution diagnostics are plotted.
#
#   Inputs:
#       s              (Dict{String,Any}) - Flow s state (conservative
#                                           variables, grid, parameters, etc.).
#       solution_base  (Dict{String,Any}) - Reference / base-flow s passed
#                                           to the plotting routine PLOT_RUNNING.
#       chemistry      (Dict{String,Any}) - Chemistry / thermodynamic model data.
#
#   Outputs:
#       s              (Dict{String,Any}) - Flow s at the final time step,
#                                           fully updated with BCs and chemistry.
#
#   Notes:
#       - Stability analysis is temporarily disabled during the run and
#         restored to its original flag value at the end.
#       - Adaptive remeshing (RESTART_SOLUTION) is invoked only when
#         s["shock"]["enabled"] and s["remesh"] are both true and one of
#         three triggering conditions is met.
#       - Running plots support single or multiple variables via the
#         s["running_plot"]["variable"] field (string or cell array).
#
# Part of: Hypersonics Stability Julia Solver - Time Marching Module

function RUN_SIMULATION(s::Dict{String,Any}, solution_base::Dict{String,Any}, chemistry::Dict{String,Any})

    ## Print simulation banner
    println()
    println("------------------------")
    println("     Run simulation     ")
    println("------------------------")
    println()

    ## Initialise graphics handles
    if s["running_plot"]["enabled"]
        if isa(s["running_plot"]["variable"], Vector)
            num_variables_to_plot = length(s["running_plot"]["variable"])
        else
            num_variables_to_plot = 1
        end

        plot_handles = Dict{String,Any}(
            "fig"                  => nothing,
            "ax"                   => Vector{Any}(nothing, num_variables_to_plot),
            "plot_object"          => Vector{Any}(nothing, num_variables_to_plot),
            "colorbar"             => Vector{Any}(nothing, num_variables_to_plot),
            "shock_line"           => Vector{Any}(nothing, num_variables_to_plot),
            "caxis_limits_per_var" => Dict{Int, Any}()
        )
    else
        plot_handles = Dict{String,Any}()
    end

    ## Initialise iteration and mesh-refinement tracking
    s["time_integration"]["iter"] = 0
    refine_interval               = 100
    refine_movement_cells         = 3
    s["index_refined"]     = s["shock"]["cell_indices"][1, 1]
    s["last_refined_iter"] = -1

    ## Declare solution_old in outer scope
    solution_old = Dict{String,Any}()

    ## Main time-marching loop
    for i in 1:s["time_integration"]["N_iter"]
        # Advect flow one time step

        s            = ADVECT_FLOW(s, chemistry)

        # Print iteration diagnostics every 10 steps
        if mod(s["time_integration"]["iter"], 10) == 0
            @printf("Iteration = %8i, Timestep = %.5f, Time = %.5f",
                    s["time_integration"]["iter"], s["time_integration"]["dt"], s["time_integration"]["t"])
            if haskey(s["time_integration"], "time_integrator") && s["time_integration"]["time_integrator"] == "Implicit_Euler"
                if haskey(s, "count_implicit_iterations")
                    @printf(", Implicit iterations = %i\n", s["count_implicit_iterations"])
                else
                    @printf("\n")
                end
            else
                @printf("\n")
            end
        end

        # Update running plots
        if s["running_plot"]["enabled"]
            plot_handles = PLOT_RUNNING(s, solution_base, plot_handles, chemistry)
        end

        # Adaptive mesh refinement when shock drifts or approaches boundary
        if s["shock"]["enabled"] && s["remesh"]
            if (maximum(s["shock"]["cell_indices"][:, 1]) > s["mesh"]["Neta"] - 3) ||
               (maximum(abs.(s["shock"]["cell_indices"][:, 1] .- s["index_refined"])) >= refine_movement_cells) ||
               (mod(s["last_refined_iter"], refine_interval) == 0)

                s["last_refined_iter"] = i
                s["restart"] = true
                s["remesh"]  = true
                println()
                print("------- Adapt mesh -------")
                disturbances = false
                solution_old = deepcopy(s)
                s = RESTART_SOLUTION(s, solution_old, chemistry, disturbances)
                s["index_refined"] = s["shock"]["cell_indices"][1, 1]
                print("-- Adapted  succesfully ------")
                println()
            end
        end
    end

    ## Finalise s
    print("Final time = ")
    @printf("%.3f \n", s["time_integration"]["t"])

    if s["shock"]["enabled"]
        s = UPDATE_FLOW_CELLS(s, chemistry)
        s = UPDATE_SHOCK_BC(s, chemistry)
    end
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry)
    s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry)
    s = UPDATE_SOUND_SPEED(s, chemistry)
    s = APPLY_BOUNDARY_CONDITIONS(s, chemistry)

    return s
end
