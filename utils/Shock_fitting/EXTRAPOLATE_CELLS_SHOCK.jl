# EXTRAPOLATE_CELLS_SHOCK - Populate auxiliary ghost cells beyond the shock.
#
#   s = EXTRAPOLATE_CELLS_SHOCK(s)
#
#   Fills the two auxiliary (ghost) cell layers immediately upstream of
#   each shocked cell using extrapolation from the downstream flow field.
#   The extrapolation order is controlled by s["shock"]["interpolate"]:
#     "1st" - zeroth-order (constant) copy of the shocked-cell value
#     "2nd" - first-order (linear) extrapolation using two downstream points
#     "3rd" - second-order (quadratic) extrapolation using three downstream points
#
#   This ensures that flux computations near the shock do not encounter
#   discontinuous jumps in the reconstructed stencil.
#
#   Inputs:
#       s (Dict) - Solution structure containing at minimum:
#                    s["mesh"]["Nchi"], s["mesh"]["Neta"]   - Grid dimensions
#                    s["shock"]["cell_indices"]             - (Nx x 1) column index of shocked cell
#                    s["shock"]["interpolate"]              - Extrapolation order: "1st", "2nd", or "3rd"
#                    s["var"]["rho"], s["var"]["rho_u"], s["var"]["rho_v"],
#                    s["var"]["rho_E"], s["var"]["p"]       - (Nx+2 x Ny+2) conservative/pressure fields
#                    s["chemistry"]["chemical_equilibrium"] - (Bool) equilibrium chemistry flag
#                    s["chemistry"]["is_chemistry_enabled"] - (Bool) chemistry enabled flag
#                    s["var"]["gamma_star"], s["var"]["cv_star"] - Thermodynamic property fields (if chemistry)
#
#   Outputs:
#       s (Dict) - Solution with auxiliary ghost cells filled by extrapolation.
#
#   Notes:
#       - Ghost cell indices are offset by +1 due to the presence of boundary
#         ghost cells in the s arrays.
#       - Uses CartesianIndex / linear indexing for vectorized indexing across all streamwise stations.
#       - KEY OPTIMIZATION: Instead of copying entire NxM arrays (copy(s["var"]["rho"])),
#         we extract only the values at the specific CartesianIndex positions needed.
#         This replaces 5-7 full-array deep copies with small vectors of length Nchi.
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

function EXTRAPOLATE_CELLS_SHOCK(s::Dict{String, Any})
    ## Extract dict refs
    var  = s["var"]::Dict{String, Any}
    chem = s["chemistry"]::Dict{String, Any}
    shock = s["shock"]::Dict{String, Any}
    Nchi = s["mesh"]["Nchi"]::Int

    rho_arr       = var["rho"]::Matrix{Float64}
    rho_u_arr     = var["rho_u"]::Matrix{Float64}
    rho_v_arr     = var["rho_v"]::Matrix{Float64}
    rho_E_arr     = var["rho_E"]::Matrix{Float64}
    p_arr         = var["p"]::Matrix{Float64}
    chem_eq       = chem["chemical_equilibrium"]::Bool
    is_chem       = chem["is_chemistry_enabled"]::Bool
    interp_order  = shock["interpolate"]::String

    ## Compute CartesianIndex arrays for shocked and neighboring cells
    i_values_s = collect(1:Nchi) .+ 1                        # Row indices (offset for ghost cells)
    j_values_s = shock["cell_indices"][:, 1]                 # Column indices of shocked cells

    idx_2  = CartesianIndex.(i_values_s, j_values_s .+ 3)  # 2nd auxiliary cell
    idx_1  = CartesianIndex.(i_values_s, j_values_s .+ 2)  # 1st auxiliary cell
    idx_0  = CartesianIndex.(i_values_s, j_values_s .+ 1)  # Shocked cell
    idx_m1 = CartesianIndex.(i_values_s, j_values_s)        # 1 cell downstream
    idx_m2 = CartesianIndex.(i_values_s, j_values_s .- 1)   # 2 cells downstream

    ## Extract values at needed positions (avoids full-array deep copies)
    rho_at_0  = rho_arr[idx_0]
    rhou_at_0 = rho_u_arr[idx_0]
    rhov_at_0 = rho_v_arr[idx_0]
    rhoE_at_0 = rho_E_arr[idx_0]
    p_at_0    = p_arr[idx_0]

    need_neq = !chem_eq && is_chem
    if need_neq
        gamma_star_arr = var["gamma_star"]::Matrix{Float64}
        cv_star_arr    = var["cv_star"]::Matrix{Float64}
        gs_at_0  = gamma_star_arr[idx_0]
        cvs_at_0 = cv_star_arr[idx_0]
    end

    ## Apply extrapolation based on selected order
    if interp_order == "1st"
        # -- Zeroth-order (constant) extrapolation --
        rho_arr[idx_1]   .= rho_at_0
        rho_u_arr[idx_1] .= rhou_at_0
        rho_v_arr[idx_1] .= rhov_at_0
        rho_E_arr[idx_1] .= rhoE_at_0
        p_arr[idx_1]     .= p_at_0
        if need_neq
            gamma_star_arr[idx_1] .= gs_at_0
            cv_star_arr[idx_1]    .= cvs_at_0
        end

        rho_arr[idx_2]   .= rho_at_0
        rho_u_arr[idx_2] .= rhou_at_0
        rho_v_arr[idx_2] .= rhov_at_0
        rho_E_arr[idx_2] .= rhoE_at_0
        p_arr[idx_2]     .= p_at_0

    elseif interp_order == "2nd"
        # -- First-order (linear) extrapolation --
        rho_at_m1  = rho_arr[idx_m1]
        rhou_at_m1 = rho_u_arr[idx_m1]
        rhov_at_m1 = rho_v_arr[idx_m1]
        rhoE_at_m1 = rho_E_arr[idx_m1]
        p_at_m1    = p_arr[idx_m1]

        rho_arr[idx_1]   .= @. 2 * rho_at_0   - rho_at_m1
        rho_u_arr[idx_1] .= @. 2 * rhou_at_0  - rhou_at_m1
        rho_v_arr[idx_1] .= @. 2 * rhov_at_0  - rhov_at_m1
        rho_E_arr[idx_1] .= @. 2 * rhoE_at_0  - rhoE_at_m1
        p_arr[idx_1]     .= @. 2 * p_at_0     - p_at_m1
        if need_neq
            gs_at_m1  = gamma_star_arr[idx_m1]
            cvs_at_m1 = cv_star_arr[idx_m1]
            gamma_star_arr[idx_1] .= @. 2 * gs_at_0  - gs_at_m1
            cv_star_arr[idx_1]    .= @. 2 * cvs_at_0 - cvs_at_m1
        end

        # For idx_2, use the now-written idx_1 values
        rho_at_1  = rho_arr[idx_1]
        rhou_at_1 = rho_u_arr[idx_1]
        rhov_at_1 = rho_v_arr[idx_1]
        rhoE_at_1 = rho_E_arr[idx_1]

        rho_arr[idx_2]   .= @. 2 * rho_at_1  - rho_at_0
        rho_u_arr[idx_2] .= @. 2 * rhou_at_1 - rhou_at_0
        rho_v_arr[idx_2] .= @. 2 * rhov_at_1 - rhov_at_0
        rho_E_arr[idx_2] .= @. 2 * rhoE_at_1 - rhoE_at_0
        if need_neq
            gs_at_1  = gamma_star_arr[idx_1]
            cvs_at_1 = cv_star_arr[idx_1]
            gamma_star_arr[idx_2] .= @. 2 * gs_at_1  - gs_at_0
            cv_star_arr[idx_2]    .= @. 2 * cvs_at_1 - cvs_at_0
        end

    elseif interp_order == "3rd"
        # -- Second-order (quadratic) extrapolation --
        rho_at_m1  = rho_arr[idx_m1]
        rhou_at_m1 = rho_u_arr[idx_m1]
        rhov_at_m1 = rho_v_arr[idx_m1]
        rhoE_at_m1 = rho_E_arr[idx_m1]
        p_at_m1    = p_arr[idx_m1]
        rho_at_m2  = rho_arr[idx_m2]
        rhou_at_m2 = rho_u_arr[idx_m2]
        rhov_at_m2 = rho_v_arr[idx_m2]
        rhoE_at_m2 = rho_E_arr[idx_m2]
        p_at_m2    = p_arr[idx_m2]

        rho_arr[idx_1]   .= @. 3 * rho_at_0   - 3 * rho_at_m1   + rho_at_m2
        rho_u_arr[idx_1] .= @. 3 * rhou_at_0  - 3 * rhou_at_m1  + rhou_at_m2
        rho_v_arr[idx_1] .= @. 3 * rhov_at_0  - 3 * rhov_at_m1  + rhov_at_m2
        rho_E_arr[idx_1] .= @. 3 * rhoE_at_0  - 3 * rhoE_at_m1  + rhoE_at_m2
        p_arr[idx_1]     .= @. 3 * p_at_0     - 3 * p_at_m1     + p_at_m2
        if need_neq
            gamma_star_arr[idx_1] .= gs_at_0
            cv_star_arr[idx_1]    .= cvs_at_0
        end

        # For idx_2, use the now-written idx_1 values
        rho_at_1  = rho_arr[idx_1]
        rhou_at_1 = rho_u_arr[idx_1]
        rhov_at_1 = rho_v_arr[idx_1]
        rhoE_at_1 = rho_E_arr[idx_1]

        rho_arr[idx_2]   .= @. 3 * rho_at_1  - 3 * rho_at_0  + rho_at_m1
        rho_u_arr[idx_2] .= @. 3 * rhou_at_1 - 3 * rhou_at_0 + rhou_at_m1
        rho_v_arr[idx_2] .= @. 3 * rhov_at_1 - 3 * rhov_at_0 + rhov_at_m1
        rho_E_arr[idx_2] .= @. 3 * rhoE_at_1 - 3 * rhoE_at_0 + rhoE_at_m1
        if need_neq
            gamma_star_arr[idx_2] .= gamma_star_arr[idx_1]
            cv_star_arr[idx_2]    .= cv_star_arr[idx_1]
        end
    end

    return s
end
