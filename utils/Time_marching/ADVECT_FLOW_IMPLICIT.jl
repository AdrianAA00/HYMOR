# ADVECT_FLOW_IMPLICIT  Advance the flow s by one implicit time step.
#
#   s = ADVECT_FLOW_IMPLICIT(s, chemistry)
#
#   Performs a single implicit advection step using iterative relaxation.
#   The routine repeatedly evaluates the PDE right-hand side and applies
#   implicit time integration until the residual drops below the specified
#   tolerance or the maximum number of iterations is reached. After
#   convergence the boundary conditions are updated, a minimum-density
#   limiter is enforced, and the CFL-based time step is recomputed.
#
#   Inputs:
#       s          (Dict{String,Any}) - Flow s state containing conservative
#                                       variables, grid data, time step, tolerance,
#                                       max_iter, and relaxation parameters.
#       chemistry  (Dict{String,Any}) - Chemistry / thermodynamic model data used by
#                                       the PDE right-hand side and boundary updates.
#
#   Outputs:
#       s          (Dict{String,Any}) - Updated flow s advanced by one time
#                                       step (s["time_integration"]["t"] incremented
#                                       by s["time_integration"]["dt"]).
#
#   Notes:
#       - Uses adaptive relaxation (ADAPT_RELAXATION) to accelerate or
#         stabilise the implicit iteration.
#       - Calls PDE, INTEGRATION_IMPLICIT, UPDATE_THERMODYNAMIC_PROPERTIES,
#         UPDATE_SHOCK_BC, APPLY_BOUNDARY_CONDITIONS, UPDATE_FIELD_UPSTREAM,
#         MIN_RHO, and CFL_TIMESTEP as sub-routines.
#
# Part of: Hypersonics Stability Julia Solver - Time Marching Module

function ADVECT_FLOW_IMPLICIT(s::Dict{String,Any}, chemistry::Dict{String,Any})


    ## Initialise implicit iteration
    residual          = 100.0
    residual_previous = 0.0
    
    s = UPDATE_FLOW_CELLS(s, chemistry)
    s = CFL_TIMESTEP(s) # Recompute CFL-based time step for the next iteration
    solution_temp     = deepcopy(s)
    count_implicit_iterations = 0

    ## Implicit iteration loop
    while residual > s["time_integration"]["tolerance"] && count_implicit_iterations < s["time_integration"]["max_iter_implicit"]
        s = ADAPT_RELAXATION(residual, residual_previous, s)

        # Compute PDE right-hand side
        solution_temp = PDE(solution_temp, chemistry)

        # Implicit time integration and residual evaluation
        solution_temp, residual = INTEGRATION_IMPLICIT(s, solution_temp)
        count_implicit_iterations = count_implicit_iterations + 1
        residual_previous = residual
    end
    s = deepcopy(solution_temp)

    ## Update boundary conditions and thermodynamic properties
    s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry)
    if s["shock"]["enabled"]
        s = UPDATE_SHOCK_BC(s, chemistry)
    end
    s = APPLY_BOUNDARY_CONDITIONS(s, chemistry)

    ## Enforce minimum density limiter
    s = MIN_RHO(s)

    ## Advance time and iteration counters

    s["time_integration"]["t"]    = s["time_integration"]["t"] + s["time_integration"]["dt"]
    s["time_integration"]["iter"] = s["time_integration"]["iter"] + 1
    s["count_implicit_iterations"] = count_implicit_iterations

    return s
end


# ADAPT_RELAXATION  Adjust the variable relaxation factor based on residual trend.
#
#   s = ADAPT_RELAXATION(residual, residual_previous, s)
#
#   Uses a sigmoid mapping of the residual change to adapt the relaxation
#   factor between 0.5 and 0.99 for the implicit iteration.

function ADAPT_RELAXATION(residual::Float64, residual_previous::Float64, s::Dict{String,Any})

    C    = 0.5
    a    = C * (residual - residual_previous)
    temp = 1.7 * s["time_integration"]["relax_factor"] / (1 + exp(-a))
    m    = min(0.99, temp)
    s["relax_factor_variable"] = max(0.5, m)

    return s
end
