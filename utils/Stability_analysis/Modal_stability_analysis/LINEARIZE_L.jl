using Printf

"""
    LINEARIZE_L(s, chemistry)

Construct the linearized Jacobian matrix of the flow operator.

Builds the sparse Jacobian matrix L that linearizes the nonlinear dynamics
operator about the current base flow. The linearization is performed via
numerical finite differences using LINEARIZE_EQUATIONS_NO_DISCONTINUITY,
and the result is validated against the nonlinear operator using
CHECK_LINEARIZATION_L.

# Arguments
- `s`: Solution Dict{String,Any} containing the base flow, grid parameters,
  and analysis flags (stability_analysis, transient_growth_analysis).
- `chemistry`: Chemistry model Dict{String,Any} for thermodynamic evaluations.

# Returns
- `s`: Updated solution Dict (linearize flag toggled during computation,
  returned with linearize = false).
- `L`: Sparse Jacobian matrix. Returns 0 if neither stability analysis
  nor transient growth analysis is requested.

# Notes
- The mesh is verified to follow the shock before linearization via
  CHECK_MESH_FOLLOWS_SHOCK.
- Chemistry equilibrium state is updated before constructing L.
- The s["linearize"] flag is set to true during construction
  and restored to false upon return.

Part of: Hypersonics Stability Julia Solver - Stability Analysis / Linear Stability Module
"""
function LINEARIZE_L(s::Dict{String,Any}, chemistry::Dict{String,Any})
    s["linearize"] = true

    ## Linearization setup
    println("------------------------")
    println("      Linearize L       ")
    println("------------------------")

    ## Prepare base flow
    if s["shock"]["enabled"]
        s = CHECK_MESH_FOLLOWS_SHOCK(s, chemistry)
    end
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry)

    ## Construct Jacobian via finite differences
    time_start = time()
    L = LINEARIZE_EQUATIONS_NO_DISCONTINUITY(s, chemistry)
    elapsed = time() - time_start
    @printf("Elapsed time: %.6f seconds\n", elapsed)

    ## Validate linearization
    CHECK_LINEARIZATION_L(L, s, chemistry)

    s["linearize"] = false
    return (s, L)
end
