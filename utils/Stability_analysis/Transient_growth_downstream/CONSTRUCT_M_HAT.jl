using SparseArrays
using LinearAlgebra

"""
    CONSTRUCT_M_HAT(s, norms)

Build the initial-condition energy-weight matrix M_hat.

Constructs a sparse diagonal matrix M_hat that defines the energy norm
for the initial condition (denominator of the gain). Unlike CONSTRUCT_M,
this matrix does not apply a spatial mask and includes shock-position
degrees of freedom when enabled.

# Arguments
- `s`: Solution Dict{String,Any} containing flow fields and mesh data.
- `norms`: 4-element vector [pressure, u-momentum, v-momentum, entropy]
  selecting which energy components to include.

# Returns
- `M_hat`: Sparse diagonal energy-weight matrix of size (4*Nx*Ny) or
  (4*Nx*Ny+Nx) when shock perturbation DoFs are included.

Part of: Hypersonics Stability Julia Solver - Transient Growth Downstream Module
"""
function CONSTRUCT_M_HAT(s::Dict{String,Any}, norms)
    ## Extract flow quantities (interior cells only)
    Nx = s["mesh"]["Nchi"]
    Ny = s["mesh"]["Neta"]
    gamma_star = s["var"]["gamma_star"][2:end-1, 2:end-1]
    p_val = s["var"]["p"][2:end-1, 2:end-1]
    u = s["var"]["rho_u"][2:end-1, 2:end-1] ./ s["var"]["rho"][2:end-1, 2:end-1]
    v_vel = s["var"]["rho_v"][2:end-1, 2:end-1] ./ s["var"]["rho"][2:end-1, 2:end-1]
    rho = s["var"]["rho"][2:end-1, 2:end-1]
    a = s["var"]["a"][2:end-1, 2:end-1]
    volume = s["mesh"]["volume"]
    indices = reshape(1:Nx*Ny, Nx, Ny)
    indices_r = reshape(1:Nx, Nx, 1)

    ## Initialize storage for sparse triplets
    M_vals = zeros(4, Nx, Ny)
    M_m = zeros(Int, 4, Nx, Ny)
    M_n = zeros(Int, 4, Nx, Ny)

    ## Pressure norm contribution
    if norms[1] == true || norms[1] == 1
        M_vals[1, :, :] .= rho .* a.^2 ./ (2 .* (gamma_star .* p_val).^2) .* volume
        M_n[1, :, :] .= 0 * Nx * Ny .+ indices
        M_m[1, :, :] .= 0 * Nx * Ny .+ indices
    end

    ## X-momentum norm contribution
    if norms[2] == true || norms[2] == 1
        M_vals[2, :, :] .= rho ./ 2 .* volume
        M_n[2, :, :] .= 1 * Nx * Ny .+ indices
        M_m[2, :, :] .= 1 * Nx * Ny .+ indices
    end

    ## Y-momentum norm contribution
    if norms[3] == true || norms[3] == 1
        M_vals[3, :, :] .= rho ./ 2 .* volume
        M_n[3, :, :] .= 2 * Nx * Ny .+ indices
        M_m[3, :, :] .= 2 * Nx * Ny .+ indices
    end

    ## Entropy norm contribution
    if norms[4] == true || norms[4] == 1
        M_vals[4, :, :] .= (gamma_star .- 1) .* p_val ./ (2 .* gamma_star) .* volume
        M_n[4, :, :] .= 3 * Nx * Ny .+ indices
        M_m[4, :, :] .= 3 * Nx * Ny .+ indices
    end

    ## Flatten triplets
    values_flat = vec(M_vals)
    vim = vec(M_m)
    vin = vec(M_n)

    ## Shock-position degrees of freedom
    if s["stability_analysis"]["perturb_shock"]
        M_r = zeros(Nx, 1)
        for i in 1:Nx
            M_r[i] = s["mesh"]["bt_area"][i, s["shock"]["cell_indices"][i, 1]] / 2
        end

        M_r_n = 4 * Nx * Ny .+ vec(indices_r)
        M_r_m = 4 * Nx * Ny .+ vec(indices_r)

        values_flat = vcat(values_flat, vec(M_r))
        vim = vcat(vim, M_r_m)
        vin = vcat(vin, M_r_n)
    end

    ## Assemble sparse matrix from triplets
    # Keep only non-zero index entries
    mask_idx = (vim .!= 0) .& (vin .!= 0)
    values_flat = values_flat[mask_idx]
    vim = vim[mask_idx]
    vin = vin[mask_idx]

    if s["stability_analysis"]["perturb_shock"]
        M_hat = sparse(vim, vin, values_flat, 4 * Nx * Ny + Nx, 4 * Nx * Ny + Nx)
    else
        M_hat = sparse(vim, vin, values_flat, 4 * Nx * Ny, 4 * Nx * Ny)
    end
    return M_hat
end
