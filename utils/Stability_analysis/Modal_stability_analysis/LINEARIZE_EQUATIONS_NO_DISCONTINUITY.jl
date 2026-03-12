using SparseArrays
using LinearAlgebra

"""
    LINEARIZE_EQUATIONS_NO_DISCONTINUITY(s, chemistry)

Build the sparse Jacobian matrix via numerical finite differences of the
nonlinear dynamics operator.

Constructs the linearized operator A = df/dq using column-wise finite
differences. The matrix maps perturbations in conservative variables (and
optionally shock position) to flux perturbations:

    [flux_rho; flux_rho_u; flux_rho_v; flux_rho_E; dx_shock] =
        A * [rho; rho_u; rho_v; rho_E; x_shock]

The Jacobian has block structure:
    A = [A11, A12;   where A11: flow -> flow,  A12: shock -> flow
         A21, A22]         A21: flow -> shock,  A22: shock -> shock

# Arguments
- `s`: Solution Dict{String,Any} with base flow, grid, and solver parameters.
- `chemistry`: Chemistry model Dict{String,Any} for thermodynamic evaluations.

# Returns
- `A`: Sparse Jacobian matrix of size:
      (4*Nx*Ny) x (4*Nx*Ny)              without shock perturbation
      (4*Nx*Ny+Nx) x (4*Nx*Ny+Nx)        with shock perturbation

# Notes
- The stencil has a 3x3 domain of influence (9 neighbors per cell).
- Nine disjoint perturbation patterns are swept to cover all cells
  without overlap in the stencil.
- Supports periodic boundary conditions for channel flows.
- Shock perturbation blocks (A12, A21, A22) are only assembled when
  s["stability_analysis"]["perturb_shock"] is true.
- Uses a spline_param switch to select between per-point and
  all-at-once shock perturbation strategies (A21 block).

Part of: Hypersonics Stability Julia Solver - Stability Analysis / Linear Stability Module
"""
function LINEARIZE_EQUATIONS_NO_DISCONTINUITY(s::Dict{String,Any}, chemistry::Dict{String,Any})
    ## Evaluate base flow operator f(x)
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry)
    # Use deepcopy(s) to prevent NON_LINEAR_DYNAMICS_NO_DISCONTINUITY from
    # mutating s in-place (Julia Dicts are mutable references, unlike MATLAB
    # structs which use copy-on-write). The original s must remain unchanged
    # so that all subsequent perturbations start from the correct base state.
    solution_new = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(deepcopy(s), chemistry)

    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]

    ## Initialize sparse matrix triplet storage
    vector_indices_m = Int[]
    vector_indices_n = Int[]
    values = Float64[]

    ## A11: Flow field perturbations -> flux responses
    #  Perturb rho, rho_u, rho_v, rho_E independently and collect
    #  flux_rho, flux_rho_u, flux_rho_v, flux_rho_E responses.

    # --- rho perturbation ---
    for index in 1:9
        (matrix_perturbation, matrix_indices_2D) = GET_PERTURBATION_MATRIX(index, s)
        solution_perturbed = deepcopy(s)
        if s["shock"]["enabled"]
            solution_perturbed["var"]["rho"][2:end-1, 2:end-1] .= s["var"]["rho"][2:end-1, 2:end-1] .+ s["shock"]["flow_cells"] .* matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        else
            solution_perturbed["var"]["rho"][2:end-1, 2:end-1] .= s["var"]["rho"][2:end-1, 2:end-1] .+ matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        end
        solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)
        solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s["stability_analysis"]["perturbation_magnitude"])
        (vim_temp, vin_temp, val_temp) = GET_OUTPUT_INDEXING_A11(solution_perturbed, matrix_indices_2D)
        append!(vector_indices_m, vim_temp)
        append!(vector_indices_n, vin_temp)
        append!(values, val_temp)
    end

    # --- rho_u perturbation ---
    for index in 1:9
        (matrix_perturbation, matrix_indices_2D) = GET_PERTURBATION_MATRIX(index, s)
        solution_perturbed = deepcopy(s)
        if s["shock"]["enabled"]
            solution_perturbed["var"]["rho_u"][2:end-1, 2:end-1] .= s["var"]["rho_u"][2:end-1, 2:end-1] .+ s["shock"]["flow_cells"] .* matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        else
            solution_perturbed["var"]["rho_u"][2:end-1, 2:end-1] .= s["var"]["rho_u"][2:end-1, 2:end-1] .+ matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        end
        solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)
        solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s["stability_analysis"]["perturbation_magnitude"])
        (vim_temp, vin_temp, val_temp) = GET_OUTPUT_INDEXING_A11(solution_perturbed, matrix_indices_2D)
        append!(vector_indices_m, vim_temp)
        append!(vector_indices_n, vin_temp .+ Nchi * Neta)
        append!(values, val_temp)
    end

    # --- rho_v perturbation ---
    for index in 1:9
        (matrix_perturbation, matrix_indices_2D) = GET_PERTURBATION_MATRIX(index, s)
        solution_perturbed = deepcopy(s)
        if s["shock"]["enabled"]
            solution_perturbed["var"]["rho_v"][2:end-1, 2:end-1] .= s["var"]["rho_v"][2:end-1, 2:end-1] .+ s["shock"]["flow_cells"] .* matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        else
            solution_perturbed["var"]["rho_v"][2:end-1, 2:end-1] .= s["var"]["rho_v"][2:end-1, 2:end-1] .+ matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        end
        solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)
        solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s["stability_analysis"]["perturbation_magnitude"])
        (vim_temp, vin_temp, val_temp) = GET_OUTPUT_INDEXING_A11(solution_perturbed, matrix_indices_2D)
        append!(vector_indices_m, vim_temp)
        append!(vector_indices_n, vin_temp .+ 2 * Nchi * Neta)
        append!(values, val_temp)
    end

    # --- rho_E perturbation ---
    for index in 1:9
        (matrix_perturbation, matrix_indices_2D) = GET_PERTURBATION_MATRIX(index, s)
        solution_perturbed = deepcopy(s)
        if s["shock"]["enabled"] && s["shock"]["feedback"]
            solution_perturbed["var"]["rho_E"][2:end-1, 2:end-1] .= s["var"]["rho_E"][2:end-1, 2:end-1] .+ s["shock"]["flow_cells"] .* matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        elseif s["shock"]["enabled"] && !s["shock"]["feedback"]
            solution_perturbed["var"]["rho_E"][2:end-1, 2:end-1] .= s["var"]["rho_E"][2:end-1, 2:end-1] .+ s["shock"]["flow_cells"] .* matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        else
            solution_perturbed["var"]["rho_E"][2:end-1, 2:end-1] .= s["var"]["rho_E"][2:end-1, 2:end-1] .+ matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        end
        solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)
        solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s["stability_analysis"]["perturbation_magnitude"])
        (vim_temp, vin_temp, val_temp) = GET_OUTPUT_INDEXING_A11(solution_perturbed, matrix_indices_2D)
        append!(vector_indices_m, vim_temp)
        append!(vector_indices_n, vin_temp .+ 3 * Nchi * Neta)
        append!(values, val_temp)
    end

    # Remove zero-index entries from boundary cells
    mask = [all(vector_indices_m[k] .!= 0) for k in eachindex(vector_indices_m)]
    vector_indices_n = vector_indices_n[mask]
    values = values[mask]
    vector_indices_m = vector_indices_m[mask]

    ## A21: Flow field perturbations -> shock speed response
    if s["stability_analysis"]["perturb_shock"]
        if s["shock"]["spline_param"] != 1
            # Per-point perturbation strategy (spline smoothing active)

            # --- rho -> shock speed ---
            for i in 1:Nchi
                for j in 1:2
                    solution_perturbed = deepcopy(s)
                    solution_perturbed["var"]["rho"][i+1, s["shock"]["cell_indices"][i, 1]-j+1] += s["stability_analysis"]["perturbation_magnitude"]
                    solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)
                    solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s["stability_analysis"]["perturbation_magnitude"])
                    (vim_temp, vin_temp, val_temp) = GET_OUTPUT_INDEXING_A21(solution_perturbed, i, s["shock"]["cell_indices"][i, 1]-j)
                    append!(vector_indices_m, vim_temp)
                    append!(vector_indices_n, vin_temp)
                    append!(values, val_temp)
                end
            end

            # --- rho_u -> shock speed ---
            for i in 1:Nchi
                for j in 1:2
                    solution_perturbed = deepcopy(s)
                    solution_perturbed["var"]["rho_u"][i+1, s["shock"]["cell_indices"][i, 1]-j+1] += s["stability_analysis"]["perturbation_magnitude"]
                    solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)
                    solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s["stability_analysis"]["perturbation_magnitude"])
                    (vim_temp, vin_temp, val_temp) = GET_OUTPUT_INDEXING_A21(solution_perturbed, i, s["shock"]["cell_indices"][i, 1]-j)
                    append!(vector_indices_m, vim_temp)
                    append!(vector_indices_n, vin_temp .+ Nchi * Neta)
                    append!(values, val_temp)
                end
            end

            # --- rho_v -> shock speed ---
            for i in 1:Nchi
                for j in 1:2
                    solution_perturbed = deepcopy(s)
                    solution_perturbed["var"]["rho_v"][i+1, s["shock"]["cell_indices"][i, 1]-j+1] += s["stability_analysis"]["perturbation_magnitude"]
                    solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)
                    solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s["stability_analysis"]["perturbation_magnitude"])
                    (vim_temp, vin_temp, val_temp) = GET_OUTPUT_INDEXING_A21(solution_perturbed, i, s["shock"]["cell_indices"][i, 1]-j)
                    append!(vector_indices_m, vim_temp)
                    append!(vector_indices_n, vin_temp .+ 2 * Nchi * Neta)
                    append!(values, val_temp)
                end
            end

            # --- rho_E -> shock speed ---
            for i in 1:Nchi
                for j in 1:2
                    solution_perturbed = deepcopy(s)
                    solution_perturbed["var"]["rho_E"][i+1, s["shock"]["cell_indices"][i, 1]-j+1] += s["stability_analysis"]["perturbation_magnitude"]
                    solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)
                    solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s["stability_analysis"]["perturbation_magnitude"])
                    (vim_temp, vin_temp, val_temp) = GET_OUTPUT_INDEXING_A21(solution_perturbed, i, s["shock"]["cell_indices"][i, 1]-j)
                    append!(vector_indices_m, vim_temp)
                    append!(vector_indices_n, vin_temp .+ 3 * Nchi * Neta)
                    append!(values, val_temp)
                end
            end

        else
            # All-at-once perturbation strategy (no spline smoothing)

            # --- rho -> shock speed ---
            for j in 1:2
                solution_perturbed = deepcopy(s)
                for i in 1:Nchi
                    solution_perturbed["var"]["rho"][i+1, s["shock"]["cell_indices"][i, 1]-j+1] += s["stability_analysis"]["perturbation_magnitude"]
                end
                solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)
                solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s["stability_analysis"]["perturbation_magnitude"])
                (vim_temp, vin_temp, val_temp) = GET_OUTPUT_INDEXING_A21_ALL(solution_perturbed, j)
                append!(vector_indices_m, vim_temp)
                append!(vector_indices_n, vin_temp)
                append!(values, val_temp)
            end

            # --- rho_u -> shock speed ---
            for j in 1:2
                solution_perturbed = deepcopy(s)
                for i in 1:Nchi
                    solution_perturbed["var"]["rho_u"][i+1, s["shock"]["cell_indices"][i, 1]-j+1] += s["stability_analysis"]["perturbation_magnitude"]
                end
                solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)
                solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s["stability_analysis"]["perturbation_magnitude"])
                (vim_temp, vin_temp, val_temp) = GET_OUTPUT_INDEXING_A21_ALL(solution_perturbed, j)
                append!(vector_indices_m, vim_temp)
                append!(vector_indices_n, vin_temp .+ Nchi * Neta)
                append!(values, val_temp)
            end

            # --- rho_v -> shock speed ---
            for j in 1:2
                solution_perturbed = deepcopy(s)
                for i in 1:Nchi
                    solution_perturbed["var"]["rho_v"][i+1, s["shock"]["cell_indices"][i, 1]-j+1] += s["stability_analysis"]["perturbation_magnitude"]
                end
                solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)
                solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s["stability_analysis"]["perturbation_magnitude"])
                (vim_temp, vin_temp, val_temp) = GET_OUTPUT_INDEXING_A21_ALL(solution_perturbed, j)
                append!(vector_indices_m, vim_temp)
                append!(vector_indices_n, vin_temp .+ 2 * Nchi * Neta)
                append!(values, val_temp)
            end

            # --- rho_E -> shock speed ---
            for j in 1:2
                solution_perturbed = deepcopy(s)
                for i in 1:Nchi
                    solution_perturbed["var"]["rho_E"][i+1, s["shock"]["cell_indices"][i, 1]-j+1] += s["stability_analysis"]["perturbation_magnitude"]
                end
                solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)
                solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s["stability_analysis"]["perturbation_magnitude"])
                (vim_temp, vin_temp, val_temp) = GET_OUTPUT_INDEXING_A21_ALL(solution_perturbed, j)
                append!(vector_indices_m, vim_temp)
                append!(vector_indices_n, vin_temp .+ 3 * Nchi * Neta)
                append!(values, val_temp)
            end
        end
    end

    ## A12 and A22: Shock position perturbation -> flux and shock speed responses
    if s["stability_analysis"]["perturb_shock"]
        for shock_point in 1:Nchi
            solution_perturbed = deepcopy(s)
            solution_perturbed = PERTURB_SHOCK(solution_perturbed, s["stability_analysis"]["perturbation_magnitude"], shock_point)
            solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)
            solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s["stability_analysis"]["perturbation_magnitude"])

            (vim_temp, vin_temp, val_temp) = GET_OUTPUT_INDEXING_A12_A22(solution_perturbed, shock_point)
            append!(vector_indices_m, vim_temp)
            append!(vector_indices_n, vin_temp .+ 4 * Nchi * Neta)
            append!(values, val_temp)
        end
    end

    ## Assemble sparse matrix
    if !s["stability_analysis"]["perturb_shock"]
        A = sparse(vector_indices_m, vector_indices_n, values,
            4 * Nchi * Neta,
            4 * Nchi * Neta)
    else
        A = sparse(vector_indices_m, vector_indices_n, values,
            4 * Nchi * Neta + Nchi,
            4 * Nchi * Neta + Nchi)
    end
    return A
end


## ========================================================================
#  A11 Block: Stencil-Based Output Indexing
#  ========================================================================
"""
    GET_OUTPUT_INDEXING_A11(s, matrix_indices_2D)

Map perturbation-response pairs into sparse matrix triplets for the A11
block (flow -> flow).

Uses a 3x3 stencil (9-point):
    x x x
    x o x
    x x x

Handles periodic boundary conditions.
"""
function GET_OUTPUT_INDEXING_A11(s::Dict{String,Any}, matrix_indices_2D::Matrix{Int})
    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]
    N = size(matrix_indices_2D, 1)

    # Pre-allocate: 4 flux outputs x 9 stencil cells x N perturbed cells
    vector_indices_n = zeros(Int, 4 * 9 * N)
    vector_indices_m = zeros(Int, 4 * 9 * N)
    values_arr = zeros(4 * 9 * N)

    # --- flux_rho ---
    for ele in 1:N
        i = matrix_indices_2D[ele, 1]
        j = matrix_indices_2D[ele, 2]
        index = 9 * (ele - 1)

        # Center cell
        vector_indices_n[index+1] = i + (j-1) * Nchi
        vector_indices_m[index+1] = i + (j-1) * Nchi
        values_arr[index+1] = s["flux"]["rho"][i, j]

        if i > 1  # East cell
            vector_indices_n[index+2] = i + (j-1) * Nchi
            vector_indices_m[index+2] = (i-1) + (j-1) * Nchi
            values_arr[index+2] = s["flux"]["rho"][i-1, j]
        end
        if i < Nchi  # West cell
            vector_indices_n[index+3] = i + (j-1) * Nchi
            vector_indices_m[index+3] = (i+1) + (j-1) * Nchi
            values_arr[index+3] = s["flux"]["rho"][i+1, j]
        end
        if j > 1  # South cell
            vector_indices_n[index+4] = i + (j-1) * Nchi
            vector_indices_m[index+4] = i + (j-2) * Nchi
            values_arr[index+4] = s["flux"]["rho"][i, j-1]
        end
        if j < Neta  # North cell
            vector_indices_n[index+5] = i + (j-1) * Nchi
            vector_indices_m[index+5] = i + (j) * Nchi
            values_arr[index+5] = s["flux"]["rho"][i, j+1]
        end
        if i > 1 && j < Neta  # North-east cell
            vector_indices_n[index+6] = i + (j-1) * Nchi
            vector_indices_m[index+6] = (i-1) + (j) * Nchi
            values_arr[index+6] = s["flux"]["rho"][i-1, j+1]
        end
        if i > 1 && j > 1  # South-east cell
            vector_indices_n[index+7] = i + (j-1) * Nchi
            vector_indices_m[index+7] = (i-1) + (j-2) * Nchi
            values_arr[index+7] = s["flux"]["rho"][i-1, j-1]
        end
        if i < Nchi && j > 1  # South-west cell
            vector_indices_n[index+8] = i + (j-1) * Nchi
            vector_indices_m[index+8] = (i+1) + (j-2) * Nchi
            values_arr[index+8] = s["flux"]["rho"][i+1, j-1]
        end
        if i < Nchi && j < Neta  # North-west cell
            vector_indices_n[index+9] = i + (j-1) * Nchi
            vector_indices_m[index+9] = (i+1) + (j) * Nchi
            values_arr[index+9] = s["flux"]["rho"][i+1, j+1]
        end

        # Periodic boundary corrections
        if s["boundary_conditions"]["boundary_chi0"]["name"] == "periodic" && s["boundary_conditions"]["boundary_chi1"]["name"] == "periodic"
            if i == 1
                vector_indices_n[index+2] = i + (j-1) * Nchi
                vector_indices_m[index+2] = (i-1) + (j-1) * Nchi + Nchi
                values_arr[index+2] = s["flux"]["rho"][Nchi, j]
                if j < Neta
                    vector_indices_n[index+6] = i + (j-1) * Nchi
                    vector_indices_m[index+6] = (i-1) + (j) * Nchi + Nchi
                    values_arr[index+6] = s["flux"]["rho"][Nchi, j+1]
                end
                if j > 1
                    vector_indices_n[index+7] = i + (j-1) * Nchi
                    vector_indices_m[index+7] = (i-1) + (j-2) * Nchi + Nchi
                    values_arr[index+7] = s["flux"]["rho"][Nchi, j-1]
                end
            end
            if i == Nchi
                vector_indices_n[index+3] = i + (j-1) * Nchi
                vector_indices_m[index+3] = (i+1) + (j-1) * Nchi - Nchi
                values_arr[index+3] = s["flux"]["rho"][1, j]
                if j > 1
                    vector_indices_n[index+8] = i + (j-1) * Nchi
                    vector_indices_m[index+8] = (i+1) + (j-2) * Nchi - Nchi
                    values_arr[index+8] = s["flux"]["rho"][1, j-1]
                end
                if j < Neta
                    vector_indices_n[index+9] = i + (j-1) * Nchi
                    vector_indices_m[index+9] = (i+1) + (j) * Nchi - Nchi
                    values_arr[index+9] = s["flux"]["rho"][1, j+1]
                end
            end
        end
        if s["boundary_conditions"]["boundary_eta0"]["name"] == "periodic" && s["boundary_conditions"]["boundary_eta1"]["name"] == "periodic"
            if j == 1
                vector_indices_n[index+4] = i + (j-1) * Nchi
                vector_indices_m[index+4] = i + (j-2) * Nchi + Nchi * Neta
                values_arr[index+4] = s["flux"]["rho"][i, Neta]
                if i > 1
                    vector_indices_n[index+7] = i + (j-1) * Nchi
                    vector_indices_m[index+7] = (i-1) + (j-2) * Nchi + Nchi * Neta
                    values_arr[index+7] = s["flux"]["rho"][i-1, Neta]
                end
                if i < Nchi
                    vector_indices_n[index+8] = i + (j-1) * Nchi
                    vector_indices_m[index+8] = (i+1) + (j-2) * Nchi + Nchi * Neta
                    values_arr[index+8] = s["flux"]["rho"][i+1, Neta]
                end
            end
            if j == Neta
                vector_indices_n[index+5] = i + (j-1) * Nchi
                vector_indices_m[index+5] = i + (j) * Nchi - Nchi * Neta
                values_arr[index+5] = s["flux"]["rho"][i, 1]
                if i > 1
                    vector_indices_n[index+6] = i + (j-1) * Nchi
                    vector_indices_m[index+6] = (i-1) + (j) * Nchi - Nchi * Neta
                    values_arr[index+6] = s["flux"]["rho"][i-1, 1]
                end
                if i < Nchi
                    vector_indices_n[index+9] = i + (j-1) * Nchi
                    vector_indices_m[index+9] = (i+1) + (j) * Nchi - Nchi * Neta
                    values_arr[index+9] = s["flux"]["rho"][i+1, 1]
                end
            end
        end
    end

    # --- flux_rho_u ---
    for ele in 1:N
        i = matrix_indices_2D[ele, 1]
        j = matrix_indices_2D[ele, 2]
        index = 9 * (ele - 1) + 9 * N

        vector_indices_n[index+1] = i + (j-1) * Nchi
        vector_indices_m[index+1] = i + (j-1) * Nchi + Nchi * Neta
        values_arr[index+1] = s["flux"]["rho_u"][i, j]

        if i > 1
            vector_indices_n[index+2] = i + (j-1) * Nchi
            vector_indices_m[index+2] = (i-1) + (j-1) * Nchi + Nchi * Neta
            values_arr[index+2] = s["flux"]["rho_u"][i-1, j]
        end
        if i < Nchi
            vector_indices_n[index+3] = i + (j-1) * Nchi
            vector_indices_m[index+3] = (i+1) + (j-1) * Nchi + Nchi * Neta
            values_arr[index+3] = s["flux"]["rho_u"][i+1, j]
        end
        if j > 1
            vector_indices_n[index+4] = i + (j-1) * Nchi
            vector_indices_m[index+4] = i + (j-2) * Nchi + Nchi * Neta
            values_arr[index+4] = s["flux"]["rho_u"][i, j-1]
        end
        if j < Neta
            vector_indices_n[index+5] = i + (j-1) * Nchi
            vector_indices_m[index+5] = i + (j) * Nchi + Nchi * Neta
            values_arr[index+5] = s["flux"]["rho_u"][i, j+1]
        end
        if i > 1 && j < Neta
            vector_indices_n[index+6] = i + (j-1) * Nchi
            vector_indices_m[index+6] = (i-1) + (j) * Nchi + Nchi * Neta
            values_arr[index+6] = s["flux"]["rho_u"][i-1, j+1]
        end
        if i > 1 && j > 1
            vector_indices_n[index+7] = i + (j-1) * Nchi
            vector_indices_m[index+7] = (i-1) + (j-2) * Nchi + Nchi * Neta
            values_arr[index+7] = s["flux"]["rho_u"][i-1, j-1]
        end
        if i < Nchi && j > 1
            vector_indices_n[index+8] = i + (j-1) * Nchi
            vector_indices_m[index+8] = (i+1) + (j-2) * Nchi + Nchi * Neta
            values_arr[index+8] = s["flux"]["rho_u"][i+1, j-1]
        end
        if i < Nchi && j < Neta
            vector_indices_n[index+9] = i + (j-1) * Nchi
            vector_indices_m[index+9] = (i+1) + (j) * Nchi + Nchi * Neta
            values_arr[index+9] = s["flux"]["rho_u"][i+1, j+1]
        end

        # Periodic boundary corrections
        if s["boundary_conditions"]["boundary_chi0"]["name"] == "periodic" && s["boundary_conditions"]["boundary_chi1"]["name"] == "periodic"
            if i == 1
                vector_indices_n[index+2] = i + (j-1) * Nchi
                vector_indices_m[index+2] = (i-1) + (j-1) * Nchi + Nchi + Nchi * Neta
                values_arr[index+2] = s["flux"]["rho_u"][Nchi, j]
                if j < Neta
                    vector_indices_n[index+6] = i + (j-1) * Nchi
                    vector_indices_m[index+6] = (i-1) + (j) * Nchi + Nchi + Nchi * Neta
                    values_arr[index+6] = s["flux"]["rho_u"][Nchi, j+1]
                end
                if j > 1
                    vector_indices_n[index+7] = i + (j-1) * Nchi
                    vector_indices_m[index+7] = (i-1) + (j-2) * Nchi + Nchi + Nchi * Neta
                    values_arr[index+7] = s["flux"]["rho_u"][Nchi, j-1]
                end
            end
            if i == Nchi
                vector_indices_n[index+3] = i + (j-1) * Nchi
                vector_indices_m[index+3] = (i+1) + (j-1) * Nchi - Nchi + Nchi * Neta
                values_arr[index+3] = s["flux"]["rho_u"][1, j]
                if j > 1
                    vector_indices_n[index+8] = i + (j-1) * Nchi
                    vector_indices_m[index+8] = (i+1) + (j-2) * Nchi - Nchi + Nchi * Neta
                    values_arr[index+8] = s["flux"]["rho_u"][1, j-1]
                end
                if j < Neta
                    vector_indices_n[index+9] = i + (j-1) * Nchi
                    vector_indices_m[index+9] = (i+1) + (j) * Nchi - Nchi + Nchi * Neta
                    values_arr[index+9] = s["flux"]["rho_u"][1, j+1]
                end
            end
        end
        if s["boundary_conditions"]["boundary_eta0"]["name"] == "periodic" && s["boundary_conditions"]["boundary_eta1"]["name"] == "periodic"
            if j == 1
                vector_indices_n[index+4] = i + (j-1) * Nchi
                vector_indices_m[index+4] = i + (j-2) * Nchi + Nchi * Neta + Nchi * Neta
                values_arr[index+4] = s["flux"]["rho_u"][i, Neta]
                if i > 1
                    vector_indices_n[index+7] = i + (j-1) * Nchi
                    vector_indices_m[index+7] = (i-1) + (j-2) * Nchi + Nchi * Neta + Nchi * Neta
                    values_arr[index+7] = s["flux"]["rho_u"][i-1, Neta]
                end
                if i < Nchi
                    vector_indices_n[index+8] = i + (j-1) * Nchi
                    vector_indices_m[index+8] = (i+1) + (j-2) * Nchi + Nchi * Neta + Nchi * Neta
                    values_arr[index+8] = s["flux"]["rho_u"][i+1, Neta]
                end
            end
            if j == Neta
                vector_indices_n[index+5] = i + (j-1) * Nchi
                vector_indices_m[index+5] = i + (j) * Nchi - Nchi * Neta + Nchi * Neta
                values_arr[index+5] = s["flux"]["rho_u"][i, 1]
                if i > 1
                    vector_indices_n[index+6] = i + (j-1) * Nchi
                    vector_indices_m[index+6] = (i-1) + (j) * Nchi - Nchi * Neta + Nchi * Neta
                    values_arr[index+6] = s["flux"]["rho_u"][i-1, 1]
                end
                if i < Nchi
                    vector_indices_n[index+9] = i + (j-1) * Nchi
                    vector_indices_m[index+9] = (i+1) + (j) * Nchi - Nchi * Neta + Nchi * Neta
                    values_arr[index+9] = s["flux"]["rho_u"][i+1, 1]
                end
            end
        end
    end

    # --- flux_rho_v ---
    for ele in 1:N
        i = matrix_indices_2D[ele, 1]
        j = matrix_indices_2D[ele, 2]
        index = 9 * (ele - 1) + 2 * 9 * N

        vector_indices_n[index+1] = i + (j-1) * Nchi
        vector_indices_m[index+1] = i + (j-1) * Nchi + 2 * Nchi * Neta
        values_arr[index+1] = s["flux"]["rho_v"][i, j]

        if i > 1
            vector_indices_n[index+2] = i + (j-1) * Nchi
            vector_indices_m[index+2] = (i-1) + (j-1) * Nchi + 2 * Nchi * Neta
            values_arr[index+2] = s["flux"]["rho_v"][i-1, j]
        end
        if i < Nchi
            vector_indices_n[index+3] = i + (j-1) * Nchi
            vector_indices_m[index+3] = (i+1) + (j-1) * Nchi + 2 * Nchi * Neta
            values_arr[index+3] = s["flux"]["rho_v"][i+1, j]
        end
        if j > 1
            vector_indices_n[index+4] = i + (j-1) * Nchi
            vector_indices_m[index+4] = i + (j-2) * Nchi + 2 * Nchi * Neta
            values_arr[index+4] = s["flux"]["rho_v"][i, j-1]
        end
        if j < Neta
            vector_indices_n[index+5] = i + (j-1) * Nchi
            vector_indices_m[index+5] = i + (j) * Nchi + 2 * Nchi * Neta
            values_arr[index+5] = s["flux"]["rho_v"][i, j+1]
        end
        if i > 1 && j < Neta
            vector_indices_n[index+6] = i + (j-1) * Nchi
            vector_indices_m[index+6] = (i-1) + (j) * Nchi + 2 * Nchi * Neta
            values_arr[index+6] = s["flux"]["rho_v"][i-1, j+1]
        end
        if i > 1 && j > 1
            vector_indices_n[index+7] = i + (j-1) * Nchi
            vector_indices_m[index+7] = (i-1) + (j-2) * Nchi + 2 * Nchi * Neta
            values_arr[index+7] = s["flux"]["rho_v"][i-1, j-1]
        end
        if i < Nchi && j > 1
            vector_indices_n[index+8] = i + (j-1) * Nchi
            vector_indices_m[index+8] = (i+1) + (j-2) * Nchi + 2 * Nchi * Neta
            values_arr[index+8] = s["flux"]["rho_v"][i+1, j-1]
        end
        if i < Nchi && j < Neta
            vector_indices_n[index+9] = i + (j-1) * Nchi
            vector_indices_m[index+9] = (i+1) + (j) * Nchi + 2 * Nchi * Neta
            values_arr[index+9] = s["flux"]["rho_v"][i+1, j+1]
        end

        # Periodic boundary corrections
        if s["boundary_conditions"]["boundary_chi0"]["name"] == "periodic" && s["boundary_conditions"]["boundary_chi1"]["name"] == "periodic"
            if i == 1
                vector_indices_n[index+2] = i + (j-1) * Nchi
                vector_indices_m[index+2] = (i-1) + (j-1) * Nchi + Nchi + 2 * Nchi * Neta
                values_arr[index+2] = s["flux"]["rho_v"][Nchi, j]
                if j < Neta
                    vector_indices_n[index+6] = i + (j-1) * Nchi
                    vector_indices_m[index+6] = (i-1) + (j) * Nchi + Nchi + 2 * Nchi * Neta
                    values_arr[index+6] = s["flux"]["rho_v"][Nchi, j+1]
                end
                if j > 1
                    vector_indices_n[index+7] = i + (j-1) * Nchi
                    vector_indices_m[index+7] = (i-1) + (j-2) * Nchi + Nchi + 2 * Nchi * Neta
                    values_arr[index+7] = s["flux"]["rho_v"][Nchi, j-1]
                end
            end
            if i == Nchi
                vector_indices_n[index+3] = i + (j-1) * Nchi
                vector_indices_m[index+3] = (i+1) + (j-1) * Nchi - Nchi + 2 * Nchi * Neta
                values_arr[index+3] = s["flux"]["rho_v"][1, j]
                if j > 1
                    vector_indices_n[index+8] = i + (j-1) * Nchi
                    vector_indices_m[index+8] = (i+1) + (j-2) * Nchi - Nchi + 2 * Nchi * Neta
                    values_arr[index+8] = s["flux"]["rho_v"][1, j-1]
                end
                if j < Neta
                    vector_indices_n[index+9] = i + (j-1) * Nchi
                    vector_indices_m[index+9] = (i+1) + (j) * Nchi - Nchi + 2 * Nchi * Neta
                    values_arr[index+9] = s["flux"]["rho_v"][1, j+1]
                end
            end
        end
        if s["boundary_conditions"]["boundary_eta0"]["name"] == "periodic" && s["boundary_conditions"]["boundary_eta1"]["name"] == "periodic"
            if j == 1
                vector_indices_n[index+4] = i + (j-1) * Nchi
                vector_indices_m[index+4] = i + (j-2) * Nchi + Nchi * Neta + 2 * Nchi * Neta
                values_arr[index+4] = s["flux"]["rho_v"][i, Neta]
                if i > 1
                    vector_indices_n[index+7] = i + (j-1) * Nchi
                    vector_indices_m[index+7] = (i-1) + (j-2) * Nchi + Nchi * Neta + 2 * Nchi * Neta
                    values_arr[index+7] = s["flux"]["rho_v"][i-1, Neta]
                end
                if i < Nchi
                    vector_indices_n[index+8] = i + (j-1) * Nchi
                    vector_indices_m[index+8] = (i+1) + (j-2) * Nchi + Nchi * Neta + 2 * Nchi * Neta
                    values_arr[index+8] = s["flux"]["rho_v"][i+1, Neta]
                end
            end
            if j == Neta
                vector_indices_n[index+5] = i + (j-1) * Nchi
                vector_indices_m[index+5] = i + (j) * Nchi - Nchi * Neta + 2 * Nchi * Neta
                values_arr[index+5] = s["flux"]["rho_v"][i, 1]
                if i > 1
                    vector_indices_n[index+6] = i + (j-1) * Nchi
                    vector_indices_m[index+6] = (i-1) + (j) * Nchi - Nchi * Neta + 2 * Nchi * Neta
                    values_arr[index+6] = s["flux"]["rho_v"][i-1, 1]
                end
                if i < Nchi
                    vector_indices_n[index+9] = i + (j-1) * Nchi
                    vector_indices_m[index+9] = (i+1) + (j) * Nchi - Nchi * Neta + 2 * Nchi * Neta
                    values_arr[index+9] = s["flux"]["rho_v"][i+1, 1]
                end
            end
        end
    end

    # --- flux_rho_E ---
    for ele in 1:N
        i = matrix_indices_2D[ele, 1]
        j = matrix_indices_2D[ele, 2]
        index = 9 * (ele - 1) + 3 * 9 * N

        vector_indices_n[index+1] = i + (j-1) * Nchi
        vector_indices_m[index+1] = i + (j-1) * Nchi + 3 * Nchi * Neta
        values_arr[index+1] = s["flux"]["rho_E"][i, j]

        if i > 1
            vector_indices_n[index+2] = i + (j-1) * Nchi
            vector_indices_m[index+2] = (i-1) + (j-1) * Nchi + 3 * Nchi * Neta
            values_arr[index+2] = s["flux"]["rho_E"][i-1, j]
        end
        if i < Nchi
            vector_indices_n[index+3] = i + (j-1) * Nchi
            vector_indices_m[index+3] = (i+1) + (j-1) * Nchi + 3 * Nchi * Neta
            values_arr[index+3] = s["flux"]["rho_E"][i+1, j]
        end
        if j > 1
            vector_indices_n[index+4] = i + (j-1) * Nchi
            vector_indices_m[index+4] = i + (j-2) * Nchi + 3 * Nchi * Neta
            values_arr[index+4] = s["flux"]["rho_E"][i, j-1]
        end
        if j < Neta
            vector_indices_n[index+5] = i + (j-1) * Nchi
            vector_indices_m[index+5] = i + (j) * Nchi + 3 * Nchi * Neta
            values_arr[index+5] = s["flux"]["rho_E"][i, j+1]
        end
        if i > 1 && j < Neta
            vector_indices_n[index+6] = i + (j-1) * Nchi
            vector_indices_m[index+6] = (i-1) + (j) * Nchi + 3 * Nchi * Neta
            values_arr[index+6] = s["flux"]["rho_E"][i-1, j+1]
        end
        if i > 1 && j > 1
            vector_indices_n[index+7] = i + (j-1) * Nchi
            vector_indices_m[index+7] = (i-1) + (j-2) * Nchi + 3 * Nchi * Neta
            values_arr[index+7] = s["flux"]["rho_E"][i-1, j-1]
        end
        if i < Nchi && j > 1
            vector_indices_n[index+8] = i + (j-1) * Nchi
            vector_indices_m[index+8] = (i+1) + (j-2) * Nchi + 3 * Nchi * Neta
            values_arr[index+8] = s["flux"]["rho_E"][i+1, j-1]
        end
        if i < Nchi && j < Neta
            vector_indices_n[index+9] = i + (j-1) * Nchi
            vector_indices_m[index+9] = (i+1) + (j) * Nchi + 3 * Nchi * Neta
            values_arr[index+9] = s["flux"]["rho_E"][i+1, j+1]
        end

        # Periodic boundary corrections
        if s["boundary_conditions"]["boundary_chi0"]["name"] == "periodic" && s["boundary_conditions"]["boundary_chi1"]["name"] == "periodic"
            if i == 1
                vector_indices_n[index+2] = i + (j-1) * Nchi
                vector_indices_m[index+2] = (i-1) + (j-1) * Nchi + Nchi + 3 * Nchi * Neta
                values_arr[index+2] = s["flux"]["rho_E"][Nchi, j]
                if j < Neta
                    vector_indices_n[index+6] = i + (j-1) * Nchi
                    vector_indices_m[index+6] = (i-1) + (j) * Nchi + Nchi + 3 * Nchi * Neta
                    values_arr[index+6] = s["flux"]["rho_E"][Nchi, j+1]
                end
                if j > 1
                    vector_indices_n[index+7] = i + (j-1) * Nchi
                    vector_indices_m[index+7] = (i-1) + (j-2) * Nchi + Nchi + 3 * Nchi * Neta
                    values_arr[index+7] = s["flux"]["rho_E"][Nchi, j-1]
                end
            end
            if i == Nchi
                vector_indices_n[index+3] = i + (j-1) * Nchi
                vector_indices_m[index+3] = (i+1) + (j-1) * Nchi - Nchi + 3 * Nchi * Neta
                values_arr[index+3] = s["flux"]["rho_E"][1, j]
                if j > 1
                    vector_indices_n[index+8] = i + (j-1) * Nchi
                    vector_indices_m[index+8] = (i+1) + (j-2) * Nchi - Nchi + 3 * Nchi * Neta
                    values_arr[index+8] = s["flux"]["rho_E"][1, j-1]
                end
                if j < Neta
                    vector_indices_n[index+9] = i + (j-1) * Nchi
                    vector_indices_m[index+9] = (i+1) + (j) * Nchi - Nchi + 3 * Nchi * Neta
                    values_arr[index+9] = s["flux"]["rho_E"][1, j+1]
                end
            end
        end
        if s["boundary_conditions"]["boundary_eta0"]["name"] == "periodic" && s["boundary_conditions"]["boundary_eta1"]["name"] == "periodic"
            if j == 1
                vector_indices_n[index+4] = i + (j-1) * Nchi
                vector_indices_m[index+4] = i + (j-2) * Nchi + Nchi * Neta + 3 * Nchi * Neta
                values_arr[index+4] = s["flux"]["rho_E"][i, Neta]
                if i > 1
                    vector_indices_n[index+7] = i + (j-1) * Nchi
                    vector_indices_m[index+7] = (i-1) + (j-2) * Nchi + Nchi * Neta + 3 * Nchi * Neta
                    values_arr[index+7] = s["flux"]["rho_E"][i-1, Neta]
                end
                if i < Nchi
                    vector_indices_n[index+8] = i + (j-1) * Nchi
                    vector_indices_m[index+8] = (i+1) + (j-2) * Nchi + Nchi * Neta + 3 * Nchi * Neta
                    values_arr[index+8] = s["flux"]["rho_E"][i+1, Neta]
                end
            end
            if j == Neta
                vector_indices_n[index+5] = i + (j-1) * Nchi
                vector_indices_m[index+5] = i + (j) * Nchi - Nchi * Neta + 3 * Nchi * Neta
                values_arr[index+5] = s["flux"]["rho_E"][i, 1]
                if i > 1
                    vector_indices_n[index+6] = i + (j-1) * Nchi
                    vector_indices_m[index+6] = (i-1) + (j) * Nchi - Nchi * Neta + 3 * Nchi * Neta
                    values_arr[index+6] = s["flux"]["rho_E"][i-1, 1]
                end
                if i < Nchi
                    vector_indices_n[index+9] = i + (j-1) * Nchi
                    vector_indices_m[index+9] = (i+1) + (j) * Nchi - Nchi * Neta + 3 * Nchi * Neta
                    values_arr[index+9] = s["flux"]["rho_E"][i+1, 1]
                end
            end
        end
    end

    return (vector_indices_m, vector_indices_n, values_arr)
end


## ========================================================================
#  A21 Block: Per-Point Shock Speed Indexing
#  ========================================================================
"""
    GET_OUTPUT_INDEXING_A21(s, i, j)

Map a single flow-field perturbation at (i,j) to shock speed responses
for all shock points (per-point strategy).
"""
function GET_OUTPUT_INDEXING_A21(s::Dict{String,Any}, i::Int, j::Int)
    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]

    if s["stability_analysis"]["perturb_shock"]
        vector_indices_n = zeros(Int, Nchi)
        vector_indices_m = zeros(Int, Nchi)
        values_arr = zeros(Nchi)

        for index in 1:Nchi
            vector_indices_n[index] = i + (j-1) * Nchi
            vector_indices_m[index] = index + 4 * Nchi * Neta
            values_arr[index] = s["shock"]["relative_increase_velocity"][index, 1]
        end
        return (vector_indices_m, vector_indices_n, values_arr)
    end
    return (Int[], Int[], Float64[])
end

## ========================================================================
#  A21 Block: All-at-Once Shock Speed Indexing
#  ========================================================================
"""
    GET_OUTPUT_INDEXING_A21_ALL(s, j)

Map simultaneous flow-field perturbations at all shocked cells to shock
speed responses (no-spline strategy).
"""
function GET_OUTPUT_INDEXING_A21_ALL(s::Dict{String,Any}, j::Int)
    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]

    if s["stability_analysis"]["perturb_shock"]
        vector_indices_n = zeros(Int, Nchi)
        vector_indices_m = zeros(Int, Nchi)
        values_arr = zeros(Nchi)

        for index in 1:Nchi
            vector_indices_n[index] = index + (s["shock"]["cell_indices"][index, 1] - j - 1) * Nchi
            vector_indices_m[index] = index + 4 * Nchi * Neta
            values_arr[index] = s["shock"]["relative_increase_velocity"][index, 1]
        end
        return (vector_indices_m, vector_indices_n, values_arr)
    end
    return (Int[], Int[], Float64[])
end


## ========================================================================
#  A12 & A22 Blocks: Shock Position Perturbation Indexing
#  ========================================================================
"""
    GET_OUTPUT_INDEXING_A12_A22(s, shock_point)

Map a shock position perturbation to flux responses (A12) and shock speed
responses (A22).

Dimensions per shock_point:
    4*Nx entries for A12 (rho, rho_u, rho_v, rho_E fluxes)
    Nx entries for A22 (shock speed)
"""
function GET_OUTPUT_INDEXING_A12_A22(s::Dict{String,Any}, shock_point::Int)
    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]

    if s["stability_analysis"]["perturb_shock"]
        vector_indices_n = zeros(Int, 5 * Nchi)
        vector_indices_m = zeros(Int, 5 * Nchi)
        values_arr = zeros(5 * Nchi)

        for i in 1:Nchi
            j = s["shock"]["cell_indices"][i, 1] - 1

            # rho flux
            index = i
            vector_indices_n[index] = shock_point
            vector_indices_m[index] = i + (j-1) * Nchi
            values_arr[index] = s["flux"]["rho"][i, j]

            # rho_u flux
            index = i + Nchi
            vector_indices_n[index] = shock_point
            vector_indices_m[index] = i + (j-1) * Nchi + Nchi * Neta
            values_arr[index] = s["flux"]["rho_u"][i, j]

            # rho_v flux
            index = i + 2 * Nchi
            vector_indices_n[index] = shock_point
            vector_indices_m[index] = i + (j-1) * Nchi + 2 * Nchi * Neta
            values_arr[index] = s["flux"]["rho_v"][i, j]

            # rho_E flux
            index = i + 3 * Nchi
            vector_indices_n[index] = shock_point
            vector_indices_m[index] = i + (j-1) * Nchi + 3 * Nchi * Neta
            values_arr[index] = s["flux"]["rho_E"][i, j]

            # Shock speed
            index = i + 4 * Nchi
            vector_indices_n[index] = shock_point
            vector_indices_m[index] = i + 4 * Nchi * Neta
            values_arr[index] = s["shock"]["relative_increase_velocity"][i, 1]
        end
        return (vector_indices_m, vector_indices_n, values_arr)
    end
    return (Int[], Int[], Float64[])
end

## ========================================================================
#  Shock Perturbation Helper
#  ========================================================================
"""
    PERTURB_SHOCK(s, perturbation_magnitude, shock_point)

Displace a single shock point in the shock-normal direction.
"""
function PERTURB_SHOCK(s::Dict{String,Any}, perturbation_magnitude::Real, shock_point::Int)
    ang = s["shock"]["beta"][shock_point, 1] - pi/2 - atan(s["freestream"]["rho_v_0"], s["freestream"]["rho_u_0"])

    s["shock"]["points_x"][shock_point, 1] += perturbation_magnitude * cos(ang)
    s["shock"]["points_y"][shock_point, 1] += perturbation_magnitude * sin(ang)
    return s
end

## ========================================================================
#  Base Flow Subtraction
#  ========================================================================
"""
    REMOVE_BASE_OUTPUT(solution_perturbed, s, perturbation_magnitude)

Compute the finite-difference Jacobian column:
(f(x+Dx) - f(x)) / Dx for all flux outputs.
"""
function REMOVE_BASE_OUTPUT(solution_perturbed::Dict{String,Any}, s::Dict{String,Any}, perturbation_magnitude::Real)
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


## ========================================================================
#  Perturbation Pattern Generator
#  ========================================================================
"""
    GET_PERTURBATION_MATRIX(index, s)

Generate a disjoint perturbation pattern for the given index (1..9).
Each pattern activates every 3rd cell in both i and j directions,
producing 9 non-overlapping patterns that together cover the full domain
without stencil interference.
"""
function GET_PERTURBATION_MATRIX(index::Int, s::Dict{String,Any})
    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]

    matrix_perturbation = zeros(Nchi, Neta)
    vector_indices = Matrix{Int}(undef, 0, 2)
    jump = 3

    if mod(index, 3) == 0
        start_i = 3
    else
        start_i = mod(index, 3)
    end

    start_j = ceil(Int, index / 3)

    for j in start_j:jump:Neta
        for i in start_i:jump:Nchi
            vector_indices = vcat(vector_indices, [i j])
            matrix_perturbation[i, j] = 1
        end
    end

    return (matrix_perturbation, vector_indices)
end
