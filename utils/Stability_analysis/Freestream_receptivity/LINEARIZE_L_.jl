using SparseArrays
using LinearAlgebra
using Printf

"""
    LINEARIZE_L_(L, s, chemistry, w_infty)

Construct the extended system matrix L_ including freestream coupling.

Builds the augmented system matrix L_ that couples the downstream
Jacobian A with the upstream boundary Jacobian B and the freestream
frequency oscillation terms. The resulting block matrix has the form:

    L_ = [L,             B,          B,          ..., B
          0_(N*p x n),   i*w_1*I_p,  0,          ..., 0
          0_(N*p x n),   0,          i*w_2*I_p,  ..., 0
          ...
          0_(N*p x n),   0,          0,          ..., i*w_N*I_p]

# Arguments
- `L`: Downstream Jacobian matrix (sparse, m x n).
- `s`: Solution Dict{String,Any} with flow fields, mesh, and analysis flags.
- `chemistry`: Chemistry model Dict{String,Any}.
- `w_infty`: Vector of freestream disturbance frequencies.

# Returns
- `s`: Updated solution Dict (linearize flag reset to false).
- `L_`: Extended system matrix (sparse); returns 0 if analysis is disabled.

# Notes
- The upstream boundary operator B is computed via finite-difference
  linearization (LINEARIZE_B) and validated (CHECK_LINEARIZATION_B).

Part of: Hypersonics Stability Julia Solver - Transient Growth Freestream Module
"""
function LINEARIZE_L_(L, s::Dict{String,Any}, chemistry::Dict{String,Any}, w_infty)
    s["linearize"] = true

    ## Display header
    println()
    println("------------------------")
    println("     Linearize L_       ")
    println("------------------------")

    ## Prepare base flow
    if s["shock"]["enabled"]
        s = CHECK_MESH_FOLLOWS_SHOCK(s, chemistry)
    end
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry)

    ## Linearize upstream boundary operator B
    time_start = time()
    B = LINEARIZE_B(s, s["stability_analysis"]["perturbation_magnitude"], chemistry)
    elapsed = time() - time_start
    @printf("Elapsed time: %.6f seconds\n", elapsed)
    CHECK_LINEARIZATION_B(B, s, chemistry)

    ## Assemble extended system matrix L_
    (m, n) = size(L)
    p_size = size(B, 2)
    N = length(w_infty)

    # Frequency diagonal: kron(diag(-i*w), I_p)
    freq_diag = kron(spdiagm(0 => -1im .* w_infty), sparse(1.0I, p_size, p_size))

    L_ = [L    repeat(B, 1, N);
          spzeros(N * p_size, n)  freq_diag]

    s["linearize"] = false
    return (s, L_)
end
