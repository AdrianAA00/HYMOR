# ADVECT_FLOW_EXPLICIT  Advance the flow one time step using explicit RK4.
#
#   s = ADVECT_FLOW_EXPLICIT(s, chemistry) performs a single
#   time step of the classical four-stage Runge-Kutta scheme on the
#   conservative flow variables, then applies boundary conditions and
#   updates the CFL-based time step.
#
#   Inputs:
#       s          (Dict{String,Any}) - Solution structure with flow field arrays
#                                       (rho, rho_u, rho_v, rho_E, etc.), grid data,
#                                       time step dt, and solver parameters.
#       chemistry  (Dict{String,Any}) - Chemistry/thermodynamic model parameters
#                                       required by PDE and boundary-condition routines.
#
#   Outputs:
#       s          (Dict{String,Any}) - Updated s after one RK4 time step.
#                                       Fields t, iter, and dt are incremented /
#                                       recomputed before return.
#
#   Algorithm (classical RK4):
#       k1 = PDE(y_n)
#       k2 = PDE(y_n + dt/2 * k1)
#       k3 = PDE(y_n + dt/2 * k2)
#       k4 = PDE(y_n + dt   * k3)
#       y_{n+1} = y_n + (dt/6)(k1 + 2*k2 + 2*k3 + k4)
#
#   The weighted accumulation is carried out through successive calls to
#   INTEGRATION, which updates the interior cells of the conservative
#   variables.
#
#   See also: ADVECT_FLOW, ADVECT_FLOW_IMPLICIT, PDE, INTEGRATION,
#             CFL_TIMESTEP
#
# Part of: Hypersonics Stability Julia Solver - Time Marching Module

function ADVECT_FLOW_EXPLICIT(s::Dict{String,Any}, chemistry::Dict{String,Any})

    ## Prepare initial state
    s = UPDATE_FLOW_CELLS(s, chemistry)
    s = CFL_TIMESTEP(s)
    dt = s["time_integration"]["dt"]
    solution_0 = deepcopy(s)
    
    ## Stage 1 (k1)
    temp = PDE(solution_0, chemistry)
    s    = INTEGRATION(s, temp, dt / 6)
    temp = INTEGRATION(solution_0, temp, dt / 2)

    ## Stage 2 (k2)
    temp = PDE(temp, chemistry)
    s    = INTEGRATION(s, temp, dt / 3)
    temp = INTEGRATION(solution_0, temp, dt / 2)

    ## Stage 3 (k3)
    temp = PDE(temp, chemistry)
    s    = INTEGRATION(s, temp, dt / 3)
    temp = INTEGRATION(solution_0, temp, dt)

    ## Stage 4 (k4)
    temp = PDE(temp, chemistry)
    s    = INTEGRATION(s, temp, dt / 6)
    
    ## Post-step: boundary conditions and corrections
    s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry)
    if s["shock"]["enabled"]
        s = UPDATE_SHOCK_BC(s, chemistry)
    end
    s = APPLY_BOUNDARY_CONDITIONS(s, chemistry)
    if s["shock"]["enabled"]
        s = UPDATE_FIELD_UPSTREAM(s)
    end

    # Enforce minimum density to avoid negative values
    s = MIN_RHO(s)

    ## Advance time counters and recompute time step
    s["time_integration"]["t"]    = s["time_integration"]["t"] + s["time_integration"]["dt"]
    s["time_integration"]["iter"] = s["time_integration"]["iter"] + 1

    return s
end
