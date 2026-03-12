using JLD2

# RESTART_FROM_FILE  Load a saved simulation state from a JLD2 file.
#
#   (s, solution_save, chemistry) = RESTART_FROM_FILE(s, solution_save, chemistry; Re=nothing, rho=nothing, T=nothing)
#
#   When restart-from-file mode is enabled, this function loads the solver
#   state from the JLD2 file specified in s["filename_restart"]. After
#   loading, any freestream parameters (Reynolds number, density,
#   temperature) provided as optional keyword arguments will override the
#   values stored in the file, enabling parameter sweeps from a common
#   base state.
#
#   Inputs:
#       s             - (Dict) Solver data structure with fields:
#                         s["restart"]           - (Bool) Global restart flag
#                         s["restart_from_file"] - (Bool) Load from file flag
#                         s["filename_restart"]  - (String) Path to the JLD2 file
#       solution_save - (Dict) Previous s state (overwritten on load).
#       chemistry     - (Dict) Chemistry model (overwritten on load).
#       Re            - (Float64 or nothing) Reynolds number override.
#       rho           - (Float64 or nothing) Freestream density override [kg/m^3].
#       T             - (Float64 or nothing) Freestream temperature override [K].
#
#   Outputs:
#       s             - (Dict) Loaded and optionally modified solver state.
#       solution_save - (Dict) Loaded previous s state.
#       chemistry     - (Dict) Loaded chemistry model.
#
#   Notes:
#       - If s["restart"] or s["restart_from_file"] is false,
#         the function returns immediately without modification.
#       - The JLD2 file is expected to contain the variables: s,
#         solution_save, and chemistry.
#
#   Part of: Hypersonics Stability MATLAB Solver - Initialization Module

function RESTART_FROM_FILE(s::Dict{String,Any}, solution_save::Union{Dict{String,Any},Nothing},
                           chemistry::Dict{String,Any};
                           Re::Union{Float64,Nothing}=nothing,
                           rho::Union{Float64,Nothing}=nothing,
                           T::Union{Float64,Nothing}=nothing)

    ## Load state from file if restart-from-file is active
    if s["restart"] && s["restart_from_file"]
        data = JLD2.load(s["filename_restart"])
        s = data["s"]
        solution_save = data["solution_save"]
        chemistry = data["chemistry"]
        println("Restarting from file: " * s["filename_restart"])
    else
        return s, solution_save, chemistry
    end

    ## Override freestream parameters if provided
    println(" ")
    println("-----------------------")
    println("Forced input parameters")

    if !isnothing(Re)
        s["freestream"]["Re"] = Re
        println("Re = $(s["freestream"]["Re"])")
    end
    if !isnothing(rho)
        s["freestream"]["rho"] = rho
        println("rho = $(s["freestream"]["rho"])")
    end
    if !isnothing(T)
        s["freestream"]["T"] = T
        println("T = $(s["freestream"]["T"])")
    end

    println("-----------------------")
    println(" ")

    return s, solution_save, chemistry
end
