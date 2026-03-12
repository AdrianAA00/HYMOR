using SparseArrays
using LinearAlgebra

"""
    CONSTRUCT_R(s)

Build the variable-transformation matrix R.

Constructs a sparse matrix R that transforms the conservative state
vector [rho, rho*u, rho*v, rho*E]' into the primitive/entropy
variables [p, u, v, S]' used in the Chu energy norm. When shock
perturbation is enabled, an identity block is appended for the
shock-position degrees of freedom.

# Arguments
- `s`: Solution Dict{String,Any} containing flow fields and mesh data.

# Returns
- `R`: Sparse transformation matrix of size (4*Nx*Ny) or
  (4*Nx*Ny+Nx) when shock perturbation DoFs are included.

Part of: Hypersonics Stability Julia Solver - Transient Growth Downstream Module
"""
function CONSTRUCT_R(s::Dict{String,Any})
    ## Extract flow quantities (interior cells only)
    Nx = s["mesh"]["Nchi"]
    Ny = s["mesh"]["Neta"]
    gamma_star = s["var"]["gamma_star"][2:end-1, 2:end-1]
    p_val = s["var"]["p"][2:end-1, 2:end-1]
    u = s["var"]["rho_u"][2:end-1, 2:end-1] ./ s["var"]["rho"][2:end-1, 2:end-1]
    v_vel = s["var"]["rho_v"][2:end-1, 2:end-1] ./ s["var"]["rho"][2:end-1, 2:end-1]
    rho = s["var"]["rho"][2:end-1, 2:end-1]

    ## Define index offsets for conservative and primitive variables
    idx_rho   = 0
    idx_rho_u = Nx * Ny
    idx_rho_v = 2 * Nx * Ny
    idx_rho_E = 3 * Nx * Ny
    idx_p = 0
    idx_u = Nx * Ny
    idx_v = 2 * Nx * Ny
    idx_S = 3 * Nx * Ny
    idx_r = 4 * Nx * Ny
    indices   = reshape(1:Nx*Ny, Nx, Ny)
    indices_r = reshape(1:Nx, Nx, 1)

    ## Initialize storage for sparse triplets
    R_vals = zeros(4, 4, Nx, Ny)
    R_n = zeros(Int, 4, 4, Nx, Ny)
    R_m = zeros(Int, 4, 4, Nx, Ny)

    ## Row 1: Pressure from conservative variables
    R_vals[1, 1, :, :] .= (gamma_star .- 1) .* (u.^2 .+ v_vel.^2) ./ 2
    R_n[1, 1, :, :] .= idx_rho .+ indices
    R_m[1, 1, :, :] .= idx_p .+ indices

    R_vals[1, 2, :, :] .= -u .* (gamma_star .- 1)
    R_n[1, 2, :, :] .= idx_rho_u .+ indices
    R_m[1, 2, :, :] .= idx_p .+ indices

    R_vals[1, 3, :, :] .= -v_vel .* (gamma_star .- 1)
    R_n[1, 3, :, :] .= idx_rho_v .+ indices
    R_m[1, 3, :, :] .= idx_p .+ indices

    R_vals[1, 4, :, :] .= gamma_star .- 1
    R_n[1, 4, :, :] .= idx_rho_E .+ indices
    R_m[1, 4, :, :] .= idx_p .+ indices

    ## Row 2: X-velocity from conservative variables
    R_vals[2, 1, :, :] .= -u ./ rho
    R_n[2, 1, :, :] .= idx_rho .+ indices
    R_m[2, 1, :, :] .= idx_u .+ indices

    R_vals[2, 2, :, :] .= 1.0 ./ rho
    R_n[2, 2, :, :] .= idx_rho_u .+ indices
    R_m[2, 2, :, :] .= idx_u .+ indices

    ## Row 3: Y-velocity from conservative variables
    R_vals[3, 1, :, :] .= -v_vel ./ rho
    R_n[3, 1, :, :] .= idx_rho .+ indices
    R_m[3, 1, :, :] .= idx_v .+ indices

    R_vals[3, 3, :, :] .= 1.0 ./ rho
    R_n[3, 3, :, :] .= idx_rho_v .+ indices
    R_m[3, 3, :, :] .= idx_v .+ indices

    ## Row 4: Entropy from conservative variables
    R_vals[4, 1, :, :] .= (u.^2 .+ v_vel.^2) ./ (2 .* p_val) .- gamma_star ./ (gamma_star .- 1) ./ rho
    R_n[4, 1, :, :] .= idx_rho .+ indices
    R_m[4, 1, :, :] .= idx_S .+ indices

    R_vals[4, 2, :, :] .= -u ./ p_val
    R_n[4, 2, :, :] .= idx_rho_u .+ indices
    R_m[4, 2, :, :] .= idx_S .+ indices

    R_vals[4, 3, :, :] .= -v_vel ./ p_val
    R_n[4, 3, :, :] .= idx_rho_v .+ indices
    R_m[4, 3, :, :] .= idx_S .+ indices

    R_vals[4, 4, :, :] .= 1.0 ./ p_val
    R_n[4, 4, :, :] .= idx_rho_E .+ indices
    R_m[4, 4, :, :] .= idx_S .+ indices

    ## Flatten triplets
    values_flat = vec(R_vals)
    vim = vec(R_m)
    vin = vec(R_n)

    ## Shock-position degrees of freedom (identity block)
    if s["stability_analysis"]["perturb_shock"]
        R_r = ones(Nx, 1)
        R_r_n = idx_r .+ vec(indices_r)
        R_r_m = idx_r .+ vec(indices_r)

        values_flat = vcat(values_flat, vec(R_r))
        vim = vcat(vim, R_r_m)
        vin = vcat(vin, R_r_n)
    end

    ## Assemble sparse matrix from triplets
    # Keep only non-zero index entries
    mask_idx = (vim .!= 0) .& (vin .!= 0)
    values_flat = values_flat[mask_idx]
    vim = vim[mask_idx]
    vin = vin[mask_idx]

    if s["stability_analysis"]["perturb_shock"]
        R = sparse(vim, vin, values_flat, 4 * Nx * Ny + Nx, 4 * Nx * Ny + Nx)
    else
        R = sparse(vim, vin, values_flat, 4 * Nx * Ny, 4 * Nx * Ny)
    end
    return R
end
