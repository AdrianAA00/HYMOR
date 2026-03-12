using SparseArrays
using LinearAlgebra

"""
    CONSTRUCT_M_INFTY_(s, norms, T, w_infty, scaling_non_temporal)

Build the extended initial-condition energy-weight matrix.

Constructs the block-diagonal energy-weight matrix for the denominator
of the transient growth gain in the freestream-coupled system. The
downstream block uses M_hat (initial-condition norm); each freestream
frequency block uses M_infty (inflow energy norm), with a factor of
1/2 applied for non-zero frequencies (since the time-averaged energy
of sin^2 or cos^2 is 1/2).

Structure:
    M_infty_ = blkdiag(M_hat, M_infty_1, M_infty_2, ..., M_infty_N)
    where M_infty_k = M_infty for w=0, or M_infty/2 for w != 0.

# Arguments
- `s`: Solution Dict{String,Any} with flow fields and mesh data.
- `norms`: 4-element vector [pressure, u-momentum, v-momentum, entropy].
- `T`: Time horizon for temporal scaling.
- `w_infty`: Vector of freestream disturbance frequencies.
- `scaling_non_temporal`: If true, use reference time T_ref instead of T for scaling.

# Returns
- `M_infty_`: Sparse block-diagonal matrix for the extended system.

Part of: Hypersonics Stability Julia Solver - Transient Growth Freestream Module
"""
function CONSTRUCT_M_INFTY_(s::Dict{String,Any}, norms, T, w_infty, scaling_non_temporal::Bool)
    ## Build downstream and freestream weight matrices
    M_hat = CONSTRUCT_M_HAT(s, norms)
    M_infty = CONSTRUCT_M_INFTY(s, norms, T, scaling_non_temporal)

    ## Create frequency-dependent blocks
    N = length(w_infty)
    M_blocks = Vector{typeof(M_infty)}(undef, N)
    for k in 1:N
        if w_infty[k] == 0
            M_blocks[k] = M_infty
        else
            M_blocks[k] = M_infty / 2
        end
    end

    ## Assemble block-diagonal matrix
    M_infty_ = blockdiag(M_hat, M_blocks...)
    return M_infty_
end


"""
    CONSTRUCT_M_INFTY(s, norms, T, scaling_non_temporal)

Build the freestream inflow energy-weight matrix.

Constructs a sparse diagonal matrix that weights the freestream
upstream state by the Chu energy norm scaled by the shock arc length
and a time or reference-time factor.

# Arguments
- `s`: Solution Dict{String,Any}.
- `norms`: 4-element vector.
- `T`: Time horizon.
- `scaling_non_temporal`: If true, use advection time T_ref for scaling.

# Returns
- `M_infty`: Sparse diagonal matrix of size (4*Nx x 4*Nx).
"""
function CONSTRUCT_M_INFTY(s::Dict{String,Any}, norms, T, scaling_non_temporal::Bool)
    ## Extract freestream quantities
    Nx = s["mesh"]["Nchi"]
    Ny = s["mesh"]["Neta"]
    gamma_star = s["freestream"]["gamma_star"] .* ones(Nx, 1)
    p_val = s["freestream"]["p_0"] .* ones(Nx, 1)
    rho = s["freestream"]["rho_0"] .* ones(Nx, 1)
    u = s["freestream"]["rho_u_0"] / s["freestream"]["rho_0"] .* ones(Nx, 1)
    v = s["freestream"]["rho_v_0"] / s["freestream"]["rho_0"] .* ones(Nx, 1)
    a = s["freestream"]["a_0"] .* ones(Nx, 1)
    volume = s["mesh"]["volume"]
    flow_cells = s["shock"]["flow_cells"]
    indices = reshape(1:Nx, Nx, 1)

    ## Compute shock arc lengths and post-shock volume
    shock_arc_length = zeros(Nx, 1)
    total_arc_length = 0.0
    total_post_shock_volume = 0.0
    volume_advect = zeros(Nx, 1)

    for i in 1:Nx
        shock_arc_length[i] = s["mesh"]["bt_area"][i, s["shock"]["cell_indices"][i, 1]] * s["mesh"]["bt_y_normal"][i, s["shock"]["cell_indices"][i, 1]]
        total_arc_length += shock_arc_length[i]
    end
    for i in 1:Nx
        for j in 1:Ny
            total_post_shock_volume += volume[i, j] * flow_cells[i, j]
        end
        volume_advect[i, 1] = total_post_shock_volume
    end

    ## Determine scaling factor
    if scaling_non_temporal
        T_ref = GET_T_REF(s)
        scaling = T_ref
    else
        scaling = T
    end

    ## Initialize storage for sparse triplets
    M_vals = zeros(4, Nx)
    M_m = zeros(Int, 4, Nx)
    M_n = zeros(Int, 4, Nx)

    ## Pressure norm contribution
    if norms[1] == true || norms[1] == 1
        M_vals[1, :] .= vec(rho .* a.^2 ./ (2 .* (gamma_star .* p_val).^2) .* shock_arc_length .* scaling)
        M_n[1, :] .= vec(0 * Nx .+ indices)
        M_m[1, :] .= vec(0 * Nx .+ indices)
    end

    ## X-momentum norm contribution
    if norms[2] == true || norms[2] == 1
        M_vals[2, :] .= vec(rho ./ 2 .* shock_arc_length .* scaling)
        M_n[2, :] .= vec(1 * Nx .+ indices)
        M_m[2, :] .= vec(1 * Nx .+ indices)
    end

    ## Y-momentum norm contribution
    if norms[3] == true || norms[3] == 1
        M_vals[3, :] .= vec(rho ./ 2 .* shock_arc_length .* scaling)
        M_n[3, :] .= vec(2 * Nx .+ indices)
        M_m[3, :] .= vec(2 * Nx .+ indices)
    end

    ## Entropy norm contribution
    if norms[4] == true || norms[4] == 1
        M_vals[4, :] .= vec((gamma_star .- 1) .* p_val ./ (2 .* gamma_star) .* shock_arc_length .* scaling)
        M_n[4, :] .= vec(3 * Nx .+ indices)
        M_m[4, :] .= vec(3 * Nx .+ indices)
    end

    ## Assemble sparse matrix from triplets
    values_flat = vec(M_vals)
    vim = vec(M_m)
    vin = vec(M_n)

    # Keep only non-zero index entries
    mask = (vim .!= 0) .& (vin .!= 0)
    values_flat = values_flat[mask]
    vim = vim[mask]
    vin = vin[mask]

    M_infty = sparse(vim, vin, values_flat, 4 * Nx, 4 * Nx)
    return M_infty
end
