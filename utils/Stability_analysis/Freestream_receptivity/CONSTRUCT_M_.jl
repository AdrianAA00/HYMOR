using SparseArrays
using LinearAlgebra

"""
    CONSTRUCT_M_(s, norms, w_infty)

Build the extended energy-weight matrix M_ for the coupled system.

Constructs the block-diagonal energy-weight matrix for the extended
state vector that includes both downstream and freestream degrees of
freedom. The downstream block uses the standard M matrix; the
freestream frequency blocks are zeroed out (energy is measured only
in the downstream domain).

Structure:
    M_ = blkdiag(M, 0, 0, ..., 0)
    where there are N = length(w_infty) zero blocks of size 4*Nx.

# Arguments
- `s`: Solution Dict{String,Any} with flow fields and mesh data.
- `norms`: 4-element vector [pressure, u-momentum, v-momentum, entropy].
- `w_infty`: Vector of freestream disturbance frequencies.

# Returns
- `M_`: Sparse block-diagonal energy-weight matrix for the extended system.

Part of: Hypersonics Stability Julia Solver - Transient Growth Freestream Module
"""
function CONSTRUCT_M_(s::Dict{String,Any}, norms, w_infty)
    ## Build downstream energy-weight matrix
    M = CONSTRUCT_M(s, norms)

    ## Append zero blocks for freestream frequency DoFs
    N = length(w_infty)
    Nchi = s["mesh"]["Nchi"]
    M_ = blockdiag(M, spzeros(4 * Nchi * N, 4 * Nchi * N))
    return M_
end
