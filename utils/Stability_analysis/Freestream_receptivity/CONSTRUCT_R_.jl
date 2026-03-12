using SparseArrays
using LinearAlgebra

"""
    CONSTRUCT_R_(s, w_infty)

Build the extended variable-transformation matrix R_ for the coupled system.

Constructs the block-diagonal transformation matrix for the extended
state vector that includes both downstream and freestream degrees of
freedom. The downstream block uses the standard R matrix; each
freestream frequency block uses R_infty (uniform freestream transform).

Structure:
    R_ = blkdiag(R, R_infty, R_infty, ..., R_infty)
    with N = length(w_infty) copies of R_infty.

# Arguments
- `s`: Solution Dict{String,Any} with flow fields and mesh data.
- `w_infty`: Vector of freestream disturbance frequencies.

# Returns
- `R_`: Sparse block-diagonal transformation matrix for the extended system.

Part of: Hypersonics Stability Julia Solver - Transient Growth Freestream Module
"""
function CONSTRUCT_R_(s::Dict{String,Any}, w_infty)
    ## Build downstream and freestream transformation matrices
    R = CONSTRUCT_R(s)
    R_infty = CONSTRUCT_R_INFTY(s)

    ## Assemble block-diagonal matrix
    N = length(w_infty)
    R_blocks = [R_infty for _ in 1:N]
    R_ = blockdiag(R, R_blocks...)
    return R_
end


"""
    CONSTRUCT_R_INFTY(s)

Build the variable-transformation matrix for uniform freestream.

Constructs the sparse transformation matrix from conservative to
primitive/entropy variables evaluated at the uniform freestream state.
This is the freestream analogue of CONSTRUCT_R.

# Arguments
- `s`: Solution Dict{String,Any} with upstream flow state.

# Returns
- `R_infty`: Sparse transformation matrix of size (4*Nx x 4*Nx).
"""
function CONSTRUCT_R_INFTY(s::Dict{String,Any})
    ## Extract uniform freestream quantities
    Nx = s["mesh"]["Nchi"]
    gamma_star = s["freestream"]["gamma_star"] .* ones(Nx, 1)
    p_val = s["freestream"]["p_0"] .* ones(Nx, 1)
    rho = s["freestream"]["rho_0"] .* ones(Nx, 1)
    u = s["freestream"]["rho_u_0"] / s["freestream"]["rho_0"] .* ones(Nx, 1)
    v = s["freestream"]["rho_v_0"] / s["freestream"]["rho_0"] .* ones(Nx, 1)

    ## Define index offsets
    idx_rho   = 0
    idx_rho_u = Nx
    idx_rho_v = 2 * Nx
    idx_rho_E = 3 * Nx
    idx_p = 0
    idx_u = Nx
    idx_v = 2 * Nx
    idx_S = 3 * Nx
    indices = reshape(1:Nx, Nx, 1)

    ## Initialize storage for sparse triplets
    R_vals = zeros(4, 4, Nx)
    R_n = zeros(Int, 4, 4, Nx)
    R_m = zeros(Int, 4, 4, Nx)

    ## Row 1: Pressure from conservative variables
    R_vals[1, 1, :] .= vec((gamma_star .- 1) .* (u.^2 .+ v.^2) ./ 2)
    R_n[1, 1, :] .= vec(idx_rho .+ indices)
    R_m[1, 1, :] .= vec(idx_p .+ indices)

    R_vals[1, 2, :] .= vec(-u .* (gamma_star .- 1))
    R_n[1, 2, :] .= vec(idx_rho_u .+ indices)
    R_m[1, 2, :] .= vec(idx_p .+ indices)

    R_vals[1, 3, :] .= vec(-v .* (gamma_star .- 1))
    R_n[1, 3, :] .= vec(idx_rho_v .+ indices)
    R_m[1, 3, :] .= vec(idx_p .+ indices)

    R_vals[1, 4, :] .= vec(gamma_star .- 1)
    R_n[1, 4, :] .= vec(idx_rho_E .+ indices)
    R_m[1, 4, :] .= vec(idx_p .+ indices)

    ## Row 2: X-velocity from conservative variables
    R_vals[2, 1, :] .= vec(-u ./ rho)
    R_n[2, 1, :] .= vec(idx_rho .+ indices)
    R_m[2, 1, :] .= vec(idx_u .+ indices)

    R_vals[2, 2, :] .= vec(1.0 ./ rho)
    R_n[2, 2, :] .= vec(idx_rho_u .+ indices)
    R_m[2, 2, :] .= vec(idx_u .+ indices)

    ## Row 3: Y-velocity from conservative variables
    R_vals[3, 1, :] .= vec(-v ./ rho)
    R_n[3, 1, :] .= vec(idx_rho .+ indices)
    R_m[3, 1, :] .= vec(idx_v .+ indices)

    R_vals[3, 3, :] .= vec(1.0 ./ rho)
    R_n[3, 3, :] .= vec(idx_rho_v .+ indices)
    R_m[3, 3, :] .= vec(idx_v .+ indices)

    ## Row 4: Entropy from conservative variables
    R_vals[4, 1, :] .= vec((u.^2 .+ v.^2) ./ (2 .* p_val) .- gamma_star ./ (gamma_star .- 1) ./ rho)
    R_n[4, 1, :] .= vec(idx_rho .+ indices)
    R_m[4, 1, :] .= vec(idx_S .+ indices)

    R_vals[4, 2, :] .= vec(-u ./ p_val)
    R_n[4, 2, :] .= vec(idx_rho_u .+ indices)
    R_m[4, 2, :] .= vec(idx_S .+ indices)

    R_vals[4, 3, :] .= vec(-v ./ p_val)
    R_n[4, 3, :] .= vec(idx_rho_v .+ indices)
    R_m[4, 3, :] .= vec(idx_S .+ indices)

    R_vals[4, 4, :] .= vec(1.0 ./ p_val)
    R_n[4, 4, :] .= vec(idx_rho_E .+ indices)
    R_m[4, 4, :] .= vec(idx_S .+ indices)

    ## Assemble sparse matrix from triplets
    values_flat = vec(R_vals)
    vim = vec(R_m)
    vin = vec(R_n)

    # Keep only non-zero index entries
    mask = (vim .!= 0) .& (vin .!= 0)
    values_flat = values_flat[mask]
    vim = vim[mask]
    vin = vin[mask]

    R_infty = sparse(vim, vin, values_flat, 4 * Nx, 4 * Nx)
    return R_infty
end
