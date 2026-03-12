using SparseArrays
using LinearAlgebra

"""
    CONSTRUCT_M(s, norms)

Build the diagonal energy-weight matrix M for downstream disturbances.

Constructs a sparse diagonal matrix M that defines the Chu energy norm
for the downstream flow domain. Each of the four conserved-variable
contributions (pressure, x-momentum, y-momentum, entropy) can be
toggled independently via the norms vector. A spatial mask may be
applied to exclude the stagnation region or select/exclude the
boundary layer.

# Arguments
- `s`: Solution Dict{String,Any} containing flow fields and mesh data.
- `norms`: 4-element vector [pressure, u-momentum, v-momentum, entropy]
  selecting which energy components to include.

# Returns
- `M`: Sparse diagonal energy-weight matrix of size (4*Nx*Ny) or
  (4*Nx*Ny+Nx) when shock perturbation DoFs are included.

Part of: Hypersonics Stability Julia Solver - Transient Growth Downstream Module
"""
function CONSTRUCT_M(s::Dict{String,Any}, norms)
    ## Extract flow quantities (interior cells only)
    Nx = s["mesh"]["Nchi"]
    Ny = s["mesh"]["Neta"]
    gamma_star = s["var"]["gamma_star"][2:end-1, 2:end-1]
    p_val = s["var"]["p"][2:end-1, 2:end-1]
    u = s["var"]["rho_u"][2:end-1, 2:end-1] ./ s["var"]["rho"][2:end-1, 2:end-1]
    v_vel = s["var"]["rho_v"][2:end-1, 2:end-1] ./ s["var"]["rho"][2:end-1, 2:end-1]
    rho = s["var"]["rho"][2:end-1, 2:end-1]
    a = s["var"]["a"][2:end-1, 2:end-1]
    flow_cells = s["shock"]["flow_cells"]
    volume = s["mesh"]["volume"] .* flow_cells
    indices = reshape(1:Nx*Ny, Nx, Ny)

    ## Initialize storage for sparse triplets
    M_vals = zeros(4, Nx, Ny)
    M_m = zeros(Int, 4, Nx, Ny)
    M_n = zeros(Int, 4, Nx, Ny)

    ## Build spatial mask
    mask = ones(Nx, Ny)

    ## Pressure norm contribution
    if norms[1] == true || norms[1] == 1
        M_vals[1, :, :] .= rho .* a.^2 ./ (2 .* (gamma_star .* p_val).^2) .* volume .* mask
        M_n[1, :, :] .= 0 * Nx * Ny .+ indices
        M_m[1, :, :] .= 0 * Nx * Ny .+ indices
    end

    ## X-momentum norm contribution
    if norms[2] == true || norms[2] == 1
        M_vals[2, :, :] .= rho ./ 2 .* volume .* mask
        M_n[2, :, :] .= 1 * Nx * Ny .+ indices
        M_m[2, :, :] .= 1 * Nx * Ny .+ indices
    end

    ## Y-momentum norm contribution
    if norms[3] == true || norms[3] == 1
        M_vals[3, :, :] .= rho ./ 2 .* volume .* mask
        M_n[3, :, :] .= 2 * Nx * Ny .+ indices
        M_m[3, :, :] .= 2 * Nx * Ny .+ indices
    end

    ## Entropy norm contribution
    if norms[4] == true || norms[4] == 1
        M_vals[4, :, :] .= (gamma_star .- 1) .* p_val ./ (2 .* gamma_star) .* volume .* mask
        M_n[4, :, :] .= 3 * Nx * Ny .+ indices
        M_m[4, :, :] .= 3 * Nx * Ny .+ indices
    end

    ## Assemble sparse matrix from triplets
    values_flat = vec(M_vals)
    vim = vec(M_m)
    vin = vec(M_n)

    # Keep only non-zero index entries
    mask_idx = (vim .!= 0) .& (vin .!= 0)
    values_flat = values_flat[mask_idx]
    vim = vim[mask_idx]
    vin = vin[mask_idx]

    if s["stability_analysis"]["perturb_shock"]
        M = sparse(vim, vin, values_flat, 4 * Nx * Ny + Nx, 4 * Nx * Ny + Nx)
    else
        M = sparse(vim, vin, values_flat, 4 * Nx * Ny, 4 * Nx * Ny)
    end
    return M
end
