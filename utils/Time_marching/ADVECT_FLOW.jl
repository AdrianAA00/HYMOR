# ADVECT_FLOW  Advance the flow s by one time step.
#
#   s = ADVECT_FLOW(s, chemistry) dispatches to the
#   appropriate time-integration routine based on the scheme stored in
#   s["time_integration"]["time_integrator"] and returns the updated flow state.
#
#   Inputs:
#       s          (Dict{String,Any}) - Solution structure containing all flow field
#                                       variables, grid data, and solver parameters.
#                                       Must include the field "time_integrator" with
#                                       one of the supported values listed below.
#       chemistry  (Dict{String,Any}) - Chemistry/thermodynamic model parameters used
#                                       by the downstream integration routines.
#
#   Outputs:
#       s          (Dict{String,Any}) - Updated s structure after one full
#                                       time step has been completed.
#
#   Supported time integrators:
#       "Explicit_RK4"   - Classical 4th-order Runge-Kutta (explicit).
#       "Implicit_Euler" - First-order backward Euler (implicit) with
#                          adaptive relaxation.
#
#   Notes:
#       This function serves as a dispatcher only; all integration work
#       is performed inside ADVECT_FLOW_EXPLICIT or ADVECT_FLOW_IMPLICIT.
#
# Part of: Hypersonics Stability Julia Solver - Time Marching Module

function ADVECT_FLOW(s::Dict{String,Any}, chemistry::Dict{String,Any})

    if s["time_integration"]["time_integrator"] == "Explicit_RK4"
        s = ADVECT_FLOW_EXPLICIT(s, chemistry)
    elseif s["time_integration"]["time_integrator"] == "Implicit_Euler"
        s = ADVECT_FLOW_IMPLICIT(s, chemistry)
    end

    return s
end
