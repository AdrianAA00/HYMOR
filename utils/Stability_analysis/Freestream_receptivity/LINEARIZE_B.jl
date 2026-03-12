using SparseArrays

"""
    LINEARIZE_B(s, perturbation, chemistry)

Compute the linearized upstream boundary operator B via finite differences.

Constructs the Jacobian matrix B that maps upstream (freestream)
perturbations in [rho, rho*u, rho*v, rho*E] to flux perturbations in
the downstream domain. The linearization uses column-wise finite
differences with a 3-color stencil to account for the tridiagonal
coupling along the shock front.

# Arguments
- `s`: Solution Dict{String,Any} with base flow, mesh, and parameters.
- `perturbation`: Finite-difference step size for Jacobian approximation.
- `chemistry`: Chemistry model Dict{String,Any}.

# Returns
- `B`: Sparse Jacobian matrix mapping upstream state to downstream fluxes.
  Size is (4*Nx*Ny x 4*Nx) or (4*Nx*Ny+Nx x 4*Nx) when shock
  perturbation DoFs are included.

# Notes
- Uses a 3-color finite-difference scheme to reduce the number of
  nonlinear evaluations from O(Nx) to O(1) per variable.
- Contains local helper functions GET_BLOCK and REMOVE_BASE_OUTPUT_B.

Part of: Hypersonics Stability Julia Solver - Transient Growth Freestream Module
"""
function LINEARIZE_B(s::Dict{String,Any}, perturbation, chemistry::Dict{String,Any})
    ## Work on a local copy to avoid mutating the caller's dict.
    ## In MATLAB structs are value types (pass-by-value), so the caller is
    ## unaffected by local modifications. Julia Dicts are mutable references,
    ## so without this deepcopy the fields added below (rho_0_upstream_p,
    ## freestream.*_p, disturbance, etc.) would leak back to the caller and
    ## persist across grid-size changes, causing dimension mismatches.
    s = deepcopy(s)

    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]

    ## Initialize sparse triplet storage
    values = Float64[]
    vector_indices_m = Int[]
    vector_indices_n = Int[]

    ## Compute base flow
    s["freestream"]["disturbance"] = Dict{String,Any}("amplitude" => [0, 0, 0, 0], "k_x" => 1.0, "k_y" => 1.0)
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry)

    s["freestream"]["rho_0_p"] = s["freestream"]["rho_0"] .* ones(Nchi, 1)
    s["freestream"]["rho_u_0_p"] = s["freestream"]["rho_u_0"] .* ones(Nchi, 1)
    s["freestream"]["rho_v_0_p"] = s["freestream"]["rho_v_0"] .* ones(Nchi, 1)
    s["freestream"]["rho_E_0_p"] = s["freestream"]["rho_E_0"] .* ones(Nchi, 1)
    s["rho_0_upstream_p"] = true

    # Use deepcopy(s) to prevent NON_LINEAR_DYNAMICS_NO_DISCONTINUITY from
    # mutating s in-place (Julia Dicts are mutable references, unlike MATLAB
    # structs which use copy-on-write). The original s must remain unchanged
    # so that all subsequent perturbations start from the correct base state.
    s_new = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(deepcopy(s), chemistry)

    ## Linearize rho column block
    count_values = 0
    count_n = 0

    for iter in 1:3
        s_perturbed = deepcopy(s)
        pert_val = zeros(Nchi, 1)
        pert_val[iter:3:end] .= 1
        s_perturbed["freestream"]["rho_0_p"] = s["freestream"]["rho_0_p"] .+ perturbation .* pert_val

        # West cell
        if iter == 1
            pert_ind_left = collect(iter+2:3:Nchi-1)
            pert_ind = pert_ind_left .+ 1
        else
            pert_ind_left = collect(iter-1:3:Nchi-1)
            pert_ind = pert_ind_left .+ 1
        end
        (values, vector_indices_n, vector_indices_m, count_values, count_n) =
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_left, values, vector_indices_n, vector_indices_m, chemistry)

        # Mid cell
        pert_ind_mid = collect(iter:3:Nchi)
        pert_ind = pert_ind_mid
        (values, vector_indices_n, vector_indices_m, count_values, count_n) =
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_mid, values, vector_indices_n, vector_indices_m, chemistry)

        # East cell
        pert_ind_right = collect(iter+1:3:Nchi)
        pert_ind = pert_ind_right .- 1
        (values, vector_indices_n, vector_indices_m, count_values, count_n) =
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_right, values, vector_indices_n, vector_indices_m, chemistry)
    end

    ## Linearize rho_u column block
    count_n += Nchi

    for iter in 1:3
        s_perturbed = deepcopy(s)
        pert_val = zeros(Nchi, 1)
        pert_val[iter:3:end] .= 1
        s_perturbed["freestream"]["rho_u_0_p"] = s["freestream"]["rho_u_0_p"] .+ perturbation .* pert_val

        # West cell
        if iter == 1
            pert_ind_left = collect(iter+2:3:Nchi-1)
            pert_ind = pert_ind_left .+ 1
        else
            pert_ind_left = collect(iter-1:3:Nchi-1)
            pert_ind = pert_ind_left .+ 1
        end
        (values, vector_indices_n, vector_indices_m, count_values, count_n) =
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_left, values, vector_indices_n, vector_indices_m, chemistry)

        # Mid cell
        pert_ind_mid = collect(iter:3:Nchi)
        pert_ind = pert_ind_mid
        (values, vector_indices_n, vector_indices_m, count_values, count_n) =
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_mid, values, vector_indices_n, vector_indices_m, chemistry)

        # East cell
        pert_ind_right = collect(iter+1:3:Nchi)
        pert_ind = pert_ind_right .- 1
        (values, vector_indices_n, vector_indices_m, count_values, count_n) =
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_right, values, vector_indices_n, vector_indices_m, chemistry)
    end

    ## Linearize rho_v column block
    count_n += Nchi

    for iter in 1:3
        s_perturbed = deepcopy(s)
        pert_val = zeros(Nchi, 1)
        pert_val[iter:3:end] .= 1
        s_perturbed["freestream"]["rho_v_0_p"] = s["freestream"]["rho_v_0_p"] .+ perturbation .* pert_val

        # West cell
        if iter == 1
            pert_ind_left = collect(iter+2:3:Nchi-1)
            pert_ind = pert_ind_left .+ 1
        else
            pert_ind_left = collect(iter-1:3:Nchi-1)
            pert_ind = pert_ind_left .+ 1
        end
        (values, vector_indices_n, vector_indices_m, count_values, count_n) =
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_left, values, vector_indices_n, vector_indices_m, chemistry)

        # Mid cell
        pert_ind_mid = collect(iter:3:Nchi)
        pert_ind = pert_ind_mid
        (values, vector_indices_n, vector_indices_m, count_values, count_n) =
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_mid, values, vector_indices_n, vector_indices_m, chemistry)

        # East cell
        pert_ind_right = collect(iter+1:3:Nchi)
        pert_ind = pert_ind_right .- 1
        (values, vector_indices_n, vector_indices_m, count_values, count_n) =
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_right, values, vector_indices_n, vector_indices_m, chemistry)
    end

    ## Linearize rho_E column block
    count_n += Nchi

    for iter in 1:3
        s_perturbed = deepcopy(s)
        pert_val = zeros(Nchi, 1)
        pert_val[iter:3:end] .= 1
        s_perturbed["freestream"]["rho_E_0_p"] = s["freestream"]["rho_E_0_p"] .+ perturbation .* pert_val

        # West cell
        if iter == 1
            pert_ind_left = collect(iter+2:3:Nchi-1)
            pert_ind = pert_ind_left .+ 1
        else
            pert_ind_left = collect(iter-1:3:Nchi-1)
            pert_ind = pert_ind_left .+ 1
        end
        (values, vector_indices_n, vector_indices_m, count_values, count_n) =
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_left, values, vector_indices_n, vector_indices_m, chemistry)

        # Mid cell
        pert_ind_mid = collect(iter:3:Nchi)
        pert_ind = pert_ind_mid
        (values, vector_indices_n, vector_indices_m, count_values, count_n) =
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_mid, values, vector_indices_n, vector_indices_m, chemistry)

        # East cell
        pert_ind_right = collect(iter+1:3:Nchi)
        pert_ind = pert_ind_right .- 1
        (values, vector_indices_n, vector_indices_m, count_values, count_n) =
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_right, values, vector_indices_n, vector_indices_m, chemistry)
    end

    ## Assemble sparse matrix
    if !s["stability_analysis"]["perturb_shock"]
        B = sparse(vector_indices_m, vector_indices_n, values, 4 * Nchi * Neta, 4 * Nchi)
    else
        B = sparse(vector_indices_m, vector_indices_n, values, 4 * Nchi * Neta + Nchi, 4 * Nchi)
    end
    return B
end


"""
    GET_BLOCK(s_perturbed, s, count_values, count_n, perturbation, pert_n, pert_m, values, vector_indices_n, vector_indices_m, chemistry)

Compute one finite-difference block of the B Jacobian.
"""
function GET_BLOCK(s_perturbed::Dict{String,Any}, s::Dict{String,Any}, count_values::Int, count_n::Int,
                   perturbation, pert_n, pert_m, values, vector_indices_n, vector_indices_m, chemistry::Dict{String,Any})
    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]
    count_m = 0
    # Use deepcopy(s_perturbed) to prevent NON_LINEAR_DYNAMICS_NO_DISCONTINUITY
    # from mutating the caller's dict in-place (Julia Dicts are mutable
    # references, unlike MATLAB structs which use copy-on-write). The same
    # s_perturbed is passed to GET_BLOCK multiple times per iter loop (West,
    # Mid, East cells), so each call must start from the unmodified state.
    s_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(deepcopy(s_perturbed), chemistry)
    s_perturbed = REMOVE_BASE_OUTPUT_B(s_perturbed, s, perturbation)
    size_ind = length(pert_m)

    # Compute linear indices: sub2ind equivalent for (pert_m, cell_indices[pert_m]-1)
    ind = [CartesianIndex(pert_m[k], s["shock"]["cell_indices"][pert_m[k]] - 1) for k in 1:size_ind]

    ## rho flux block
    for k in 1:size_ind
        push!(values, s_perturbed["flux"]["rho"][ind[k]])
        push!(vector_indices_m, (s["shock"]["cell_indices"][pert_m[k]] - 2) * Nchi + pert_m[k] + count_m)
        push!(vector_indices_n, pert_n[k] + count_n)
    end
    count_values += size_ind
    count_m += Nchi * Neta

    ## rho_u flux block
    for k in 1:size_ind
        push!(values, s_perturbed["flux"]["rho_u"][ind[k]])
        push!(vector_indices_m, (s["shock"]["cell_indices"][pert_m[k]] - 2) * Nchi + pert_m[k] + count_m)
        push!(vector_indices_n, pert_n[k] + count_n)
    end
    count_values += size_ind
    count_m += Nchi * Neta

    ## rho_v flux block
    for k in 1:size_ind
        push!(values, s_perturbed["flux"]["rho_v"][ind[k]])
        push!(vector_indices_m, (s["shock"]["cell_indices"][pert_m[k]] - 2) * Nchi + pert_m[k] + count_m)
        push!(vector_indices_n, pert_n[k] + count_n)
    end
    count_values += size_ind
    count_m += Nchi * Neta

    ## rho_E flux block
    for k in 1:size_ind
        push!(values, s_perturbed["flux"]["rho_E"][ind[k]])
        push!(vector_indices_m, (s["shock"]["cell_indices"][pert_m[k]] - 2) * Nchi + pert_m[k] + count_m)
        push!(vector_indices_n, pert_n[k] + count_n)
    end
    count_values += size_ind
    count_m += Nchi * Neta

    ## Shock velocity block (if enabled)
    if s["stability_analysis"]["perturb_shock"]
        for k in 1:size_ind
            push!(values, s_perturbed["shock"]["relative_increase_velocity"][pert_m[k], 1])
            push!(vector_indices_m, pert_m[k] + count_m)
            push!(vector_indices_n, pert_n[k] + count_n)
        end
        count_values += size_ind
        count_m += Nchi
    end

    return (values, vector_indices_n, vector_indices_m, count_values, count_n)
end


"""
    REMOVE_BASE_OUTPUT_B(solution_perturbed, s, perturbation_magnitude)

Extract the linearized flux perturbation from the nonlinear evaluation.
Computes (f(x + dx) - f(x)) / dx for each flux component.
"""
function REMOVE_BASE_OUTPUT_B(solution_perturbed::Dict{String,Any}, s::Dict{String,Any}, perturbation_magnitude)
    solution_perturbation_only = deepcopy(s)
    solution_perturbation_only["flux"]["rho"] = (solution_perturbed["flux"]["rho"] .- s["flux"]["rho"]) ./ perturbation_magnitude
    solution_perturbation_only["flux"]["rho_u"] = (solution_perturbed["flux"]["rho_u"] .- s["flux"]["rho_u"]) ./ perturbation_magnitude
    solution_perturbation_only["flux"]["rho_v"] = (solution_perturbed["flux"]["rho_v"] .- s["flux"]["rho_v"]) ./ perturbation_magnitude
    solution_perturbation_only["flux"]["rho_E"] = (solution_perturbed["flux"]["rho_E"] .- s["flux"]["rho_E"]) ./ perturbation_magnitude
    if s["shock"]["enabled"]
        solution_perturbation_only["shock"]["relative_increase_velocity"] = (solution_perturbed["shock"]["relative_increase_velocity"] .- s["shock"]["relative_increase_velocity"]) ./ perturbation_magnitude
    end
    return solution_perturbation_only
end
