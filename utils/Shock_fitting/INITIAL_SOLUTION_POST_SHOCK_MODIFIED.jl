# INITIAL_SOLUTION_POST_SHOCK_MODIFIED - Apply a blunt-body profile correction to the post-shock field.
#
#   s = INITIAL_SOLUTION_POST_SHOCK_MODIFIED(s, chemistry)
#
#   Modifies the initial post-shock s by applying a parabolic
#   density/energy amplification and an exponentially decaying normal-
#   velocity reduction from the shock toward the body. This produces a
#   more physical initial guess for blunt-body flows that accelerates
#   convergence of the flow solver.
#
#   Inputs:
#       s  (Dict) - Solution structure (already initialized by
#                            INITIAL_SOLUTION_POST_SHOCK) with fields:
#                            s["mesh"]["Nchi"], s["mesh"]["Neta"]  - Grid dimensions
#                            s["shock"]["cell_indices"]     - (Nx x 1) shocked-cell indices
#                            s["var"]["rho"], s["var"]["rho_u"], s["var"]["rho_v"], s["var"]["rho_E"]
#                            s["mesh"]["bt_x_normal"], s["mesh"]["bt_y_normal"] - Shock-normal unit vectors
#       chemistry (Dict) - Chemistry model (unused in this function but
#                            kept for interface consistency).
#
#   Outputs:
#       s  (Dict) - Modified s with adjusted density, momentum,
#                            and total energy fields.
#
#   Notes:
#       - Density and energy are amplified by a factor of
#         (1 + 0.5 * (1 - (j/(N-1))^2)), which peaks at the shock and
#         vanishes at the body.
#       - Normal velocity is reduced exponentially with decay index = 5.
#       - Kinetic energy is corrected so that total energy remains consistent
#         with the modified momentum field.
#
#   See also: INITIAL_SOLUTION_POST_SHOCK
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

function INITIAL_SOLUTION_POST_SHOCK_MODIFIED(s::Dict{String, Any}, chemistry::Dict{String, Any})
    ## Extract dict refs
    var   = s["var"]::Dict{String, Any}
    mesh  = s["mesh"]::Dict{String, Any}
    shock = s["shock"]::Dict{String, Any}
    Nchi  = mesh["Nchi"]::Int

    rho_arr   = var["rho"]::Matrix{Float64}
    rho_u_arr = var["rho_u"]::Matrix{Float64}
    rho_v_arr = var["rho_v"]::Matrix{Float64}
    rho_E_arr = var["rho_E"]::Matrix{Float64}
    bt_x_normal = mesh["bt_x_normal"]::Matrix{Float64}
    bt_y_normal = mesh["bt_y_normal"]::Matrix{Float64}
    cell_indices = shock["cell_indices"]

    ## Modification loop with inline kinetic energy correction.
    ## Eliminates two full-array temporaries (k0, k1) by computing the
    ## kinetic energy difference per modified entry within the loop.
    decay_index = 5

    for i in 1:Nchi
        N_end = cell_indices[i, 1]
        inv_Nm1 = 1.0 / (N_end - 1)
        bx_i = bt_x_normal[i, 1]
        by_i = bt_y_normal[i, 1]
        @inbounds for j in 1:N_end
            ii = i + 1
            jj = j + 1
            # Save old kinetic energy before modification
            k_old = (rho_u_arr[ii, jj] * rho_u_arr[ii, jj] + rho_v_arr[ii, jj] * rho_v_arr[ii, jj]) / rho_arr[ii, jj] * 0.5

            normal = rho_u_arr[ii, jj] * bx_i + rho_v_arr[ii, jj] * by_i

            fac = 1.0 + 0.5 * (1.0 - (j * inv_Nm1)^2)
            rho_arr[ii, jj]   *= fac
            rho_E_arr[ii, jj] *= fac

            decay = exp(-decay_index * j * inv_Nm1)
            rho_u_arr[ii, jj] -= normal * bx_i * decay
            rho_v_arr[ii, jj] -= normal * by_i * decay

            # Correct total energy for modified kinetic energy
            k_new = (rho_u_arr[ii, jj] * rho_u_arr[ii, jj] + rho_v_arr[ii, jj] * rho_v_arr[ii, jj]) / rho_arr[ii, jj] * 0.5
            rho_E_arr[ii, jj] += k_new - k_old
        end
    end

    return s
end
