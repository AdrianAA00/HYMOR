# INITIALIZATION  Master initialization routine for the flow solver.
#
#   s = INITIALIZATION(s, solution_save, chemistry)
#
#   Orchestrates the complete initialization of the flow solver. Depending
#   on the configuration flags within the s structure, this function
#   either starts a fresh simulation via START_SOLUTION or restarts from a
#   previously saved state via RESTART_SOLUTION. It also stores the initial
#   shock geometry and disables linearization at the end of initialization.
#
#   Inputs:
#       s             - (Dict) Main solver data structure containing all
#                       simulation parameters, mesh data, and flow field
#                       variables. Must include the flags:
#                         s["restart"] - (Bool) Whether to restart
#       solution_save - (Dict or nothing) Previously saved s structure
#                       for restart. Pass nothing for a fresh start.
#       chemistry     - (Dict) Chemistry model containing thermodynamic
#                       evaluation functions (e.g., eval_e, eval_a, eval_mu).
#
#   Outputs:
#       s             - (Dict) Fully initialized solver data structure
#                       with all fields populated and ready for time
#                       integration.
#
#   Notes:
#       - When s["restart"] is true, a valid solution_save must be
#         provided; otherwise an error is raised.
#       - The field s["linearize"] is always set to false upon exit.
#
#   See also: START_SOLUTION, RESTART_SOLUTION
#
#   Part of: Hypersonics Stability MATLAB Solver - Initialization Module

function INITIALIZATION(s::Dict{String,Any}, solution_save::Union{Dict{String,Any},Nothing}=nothing, chemistry::Dict{String,Any}=Dict{String,Any}())

    ## Verify input data
    s = VERIFY_INPUT_DATA(s)

    ## Validate restart inputs
    if isnothing(solution_save)
        solution_save = nothing
    end

    if s["restart"] && (isnothing(solution_save))
        error("Restart option activated, but there is no solution_save structure to restart from.")
    end

    ## Initialize or restart the s
    if s["restart"] == false
        s = START_SOLUTION(s, chemistry)
    else
        disturbances = false
        s = RESTART_SOLUTION(s, solution_save, chemistry, disturbances)
    end

    ## Disable linearization
    s["linearize"] = false

    return s
end
