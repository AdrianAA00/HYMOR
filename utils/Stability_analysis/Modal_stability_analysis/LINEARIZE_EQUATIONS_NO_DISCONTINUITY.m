function A = LINEARIZE_EQUATIONS_NO_DISCONTINUITY(s, chemistry)
% LINEARIZE_EQUATIONS_NO_DISCONTINUITY  Build the sparse Jacobian matrix via
%   numerical finite differences of the nonlinear dynamics operator.
%
%   A = LINEARIZE_EQUATIONS_NO_DISCONTINUITY(s, chemistry) constructs
%   the linearized operator A = df/dq using column-wise finite differences.
%   The matrix maps perturbations in conservative variables (and optionally
%   shock position) to flux perturbations:
%
%       [flux_rho; flux_rho_u; flux_rho_v; flux_rho_E; dx_shock] =
%           A * [rho; rho_u; rho_v; rho_E; x_shock]
%
%   The Jacobian has block structure:
%       A = [A11, A12;   where A11: flow -> flow,  A12: shock -> flow
%            A21, A22]         A21: flow -> shock,  A22: shock -> shock
%
%   Inputs:
%       s  - Solution structure with base flow, grid, and solver
%                   parameters.
%       chemistry - Chemistry model structure for thermodynamic evaluations.
%
%   Outputs:
%       A - Sparse Jacobian matrix of size:
%             (4*Nx*Ny) x (4*Nx*Ny)              without shock perturbation
%             (4*Nx*Ny+Nx) x (4*Nx*Ny+Nx)        with shock perturbation
%
%   Notes:
%       - The stencil has a 3x3 domain of influence (9 neighbors per cell).
%       - Nine disjoint perturbation patterns are swept to cover all cells
%         without overlap in the stencil.
%       - Supports periodic boundary conditions for channel flows.
%       - Shock perturbation blocks (A12, A21, A22) are only assembled when
%         s.stability_analysis.perturb_shock is true.
%       - Uses a spline_param switch to select between per-point and
%         all-at-once shock perturbation strategies (A21 block).
%
% Part of: Hypersonics Stability MATLAB Solver - Stability Analysis / Linear Stability Module

    %% Evaluate base flow operator f(x)
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry);
    [solution_new] = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s, chemistry);

    %% Initialize sparse matrix triplet storage
    vector_indices_m = [];
    vector_indices_n = [];
    values = [];

    %% A11: Flow field perturbations -> flux responses
    %  Perturb rho, rho_u, rho_v, rho_E independently and collect
    %  flux_rho, flux_rho_u, flux_rho_v, flux_rho_E responses.

    % --- rho perturbation ---
    for index = 1:9
        [matrix_perturbation, matrix_indices_2D] = GET_PERTURBATION_MATRIX(index, s);
        solution_perturbed = s;
        if s.shock.enabled
            solution_perturbed.var.rho(2:end-1, 2:end-1) = s.var.rho(2:end-1, 2:end-1) + s.shock.flow_cells .* matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        else
            solution_perturbed.var.rho(2:end-1, 2:end-1) = s.var.rho(2:end-1, 2:end-1) + matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        end
        solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
        solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s.stability_analysis.perturbation_magnitude);
        [vector_indices_m_temp, vector_indices_n_temp, values_temp] = GET_OUTPUT_INDEXING_A11(solution_perturbed, matrix_indices_2D);
        vector_indices_m = cat(1, vector_indices_m, vector_indices_m_temp);
        vector_indices_n = cat(1, vector_indices_n, vector_indices_n_temp);
        values = cat(1, values, values_temp);
    end

    % --- rho_u perturbation ---
    for index = 1:9
        [matrix_perturbation, matrix_indices_2D] = GET_PERTURBATION_MATRIX(index, s);
        solution_perturbed = s;
        if s.shock.enabled
            solution_perturbed.var.rho_u(2:end-1, 2:end-1) = s.var.rho_u(2:end-1, 2:end-1) + s.shock.flow_cells .* matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        else
            solution_perturbed.var.rho_u(2:end-1, 2:end-1) = s.var.rho_u(2:end-1, 2:end-1) + matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        end
        solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
        solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s.stability_analysis.perturbation_magnitude);
        [vector_indices_m_temp, vector_indices_n_temp, values_temp] = GET_OUTPUT_INDEXING_A11(solution_perturbed, matrix_indices_2D);
        vector_indices_m = cat(1, vector_indices_m, vector_indices_m_temp);
        vector_indices_n = cat(1, vector_indices_n, vector_indices_n_temp + s.mesh.Nchi * s.mesh.Neta);
        values = cat(1, values, values_temp);
    end

    % --- rho_v perturbation ---
    for index = 1:9
        [matrix_perturbation, matrix_indices_2D] = GET_PERTURBATION_MATRIX(index, s);
        solution_perturbed = s;
        if s.shock.enabled
            solution_perturbed.var.rho_v(2:end-1, 2:end-1) = s.var.rho_v(2:end-1, 2:end-1) + s.shock.flow_cells .* matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        else
            solution_perturbed.var.rho_v(2:end-1, 2:end-1) = s.var.rho_v(2:end-1, 2:end-1) + matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        end
        solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
        solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s.stability_analysis.perturbation_magnitude);
        [vector_indices_m_temp, vector_indices_n_temp, values_temp] = GET_OUTPUT_INDEXING_A11(solution_perturbed, matrix_indices_2D);
        vector_indices_m = cat(1, vector_indices_m, vector_indices_m_temp);
        vector_indices_n = cat(1, vector_indices_n, vector_indices_n_temp + 2 * s.mesh.Nchi * s.mesh.Neta);
        values = cat(1, values, values_temp);
    end

    % --- rho_E perturbation ---
    for index = 1:9
        [matrix_perturbation, matrix_indices_2D] = GET_PERTURBATION_MATRIX(index, s);
        solution_perturbed = s;
        if s.shock.enabled && s.shock.feedback
            solution_perturbed.var.rho_E(2:end-1, 2:end-1) = s.var.rho_E(2:end-1, 2:end-1) + s.shock.flow_cells .* matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        elseif s.shock.enabled && ~s.shock.feedback
            solution_perturbed.var.rho_E(2:end-1, 2:end-1) = s.var.rho_E(2:end-1, 2:end-1) + s.shock.flow_cells .* matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        else
            solution_perturbed.var.rho_E(2:end-1, 2:end-1) = s.var.rho_E(2:end-1, 2:end-1) + matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        end
        solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
        solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s.stability_analysis.perturbation_magnitude);
        [vector_indices_m_temp, vector_indices_n_temp, values_temp] = GET_OUTPUT_INDEXING_A11(solution_perturbed, matrix_indices_2D);
        vector_indices_m = cat(1, vector_indices_m, vector_indices_m_temp);
        vector_indices_n = cat(1, vector_indices_n, vector_indices_n_temp + 3 * s.mesh.Nchi * s.mesh.Neta);
        values = cat(1, values, values_temp);
    end

    % Remove zero-index entries from boundary cells
    vector_indices_n = vector_indices_n(all(vector_indices_m, 2), :);
    values = values(all(vector_indices_m, 2), :);
    vector_indices_m = vector_indices_m(all(vector_indices_m, 2), :);

    %% A21: Flow field perturbations -> shock speed response
    if s.stability_analysis.perturb_shock
        if s.shock.spline_param ~= 1
            % Per-point perturbation strategy (spline smoothing active)

            % --- rho -> shock speed ---
            for i = 1:s.mesh.Nchi
                for j = 1:2
                    solution_perturbed = s;
                    solution_perturbed.var.rho(i+1, s.shock.cell_indices(i, 1)-j+1) = solution_perturbed.var.rho(i+1, s.shock.cell_indices(i, 1)-j+1) + s.stability_analysis.perturbation_magnitude;
                    solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
                    solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s.stability_analysis.perturbation_magnitude);
                    [vector_indices_m_temp, vector_indices_n_temp, values_temp] = GET_OUTPUT_INDEXING_A21(solution_perturbed, i, s.shock.cell_indices(i, 1)-j);
                    vector_indices_m = cat(1, vector_indices_m, vector_indices_m_temp);
                    vector_indices_n = cat(1, vector_indices_n, vector_indices_n_temp);
                    values = cat(1, values, values_temp);
                end
            end

            % --- rho_u -> shock speed ---
            for i = 1:s.mesh.Nchi
                for j = 1:2
                    solution_perturbed = s;
                    solution_perturbed.var.rho_u(i+1, s.shock.cell_indices(i, 1)-j+1) = solution_perturbed.var.rho_u(i+1, s.shock.cell_indices(i, 1)-j+1) + s.stability_analysis.perturbation_magnitude;
                    solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
                    solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s.stability_analysis.perturbation_magnitude);
                    [vector_indices_m_temp, vector_indices_n_temp, values_temp] = GET_OUTPUT_INDEXING_A21(solution_perturbed, i, s.shock.cell_indices(i, 1)-j);
                    vector_indices_m = cat(1, vector_indices_m, vector_indices_m_temp);
                    vector_indices_n = cat(1, vector_indices_n, vector_indices_n_temp + s.mesh.Nchi * s.mesh.Neta);
                    values = cat(1, values, values_temp);
                end
            end

            % --- rho_v -> shock speed ---
            for i = 1:s.mesh.Nchi
                for j = 1:2
                    solution_perturbed = s;
                    solution_perturbed.var.rho_v(i+1, s.shock.cell_indices(i, 1)-j+1) = solution_perturbed.var.rho_v(i+1, s.shock.cell_indices(i, 1)-j+1) + s.stability_analysis.perturbation_magnitude;
                    solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
                    solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s.stability_analysis.perturbation_magnitude);
                    [vector_indices_m_temp, vector_indices_n_temp, values_temp] = GET_OUTPUT_INDEXING_A21(solution_perturbed, i, s.shock.cell_indices(i, 1)-j);
                    vector_indices_m = cat(1, vector_indices_m, vector_indices_m_temp);
                    vector_indices_n = cat(1, vector_indices_n, vector_indices_n_temp + 2 * s.mesh.Nchi * s.mesh.Neta);
                    values = cat(1, values, values_temp);
                end
            end

            % --- rho_E -> shock speed ---
            for i = 1:s.mesh.Nchi
                for j = 1:2
                    solution_perturbed = s;
                    solution_perturbed.var.rho_E(i+1, s.shock.cell_indices(i, 1)-j+1) = solution_perturbed.var.rho_E(i+1, s.shock.cell_indices(i, 1)-j+1) + s.stability_analysis.perturbation_magnitude;
                    solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
                    solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s.stability_analysis.perturbation_magnitude);
                    [vector_indices_m_temp, vector_indices_n_temp, values_temp] = GET_OUTPUT_INDEXING_A21(solution_perturbed, i, s.shock.cell_indices(i, 1)-j);
                    vector_indices_m = cat(1, vector_indices_m, vector_indices_m_temp);
                    vector_indices_n = cat(1, vector_indices_n, vector_indices_n_temp + 3 * s.mesh.Nchi * s.mesh.Neta);
                    values = cat(1, values, values_temp);
                end
            end

        else
            % All-at-once perturbation strategy (no spline smoothing)

            % --- rho -> shock speed ---
            for j = 1:2
                solution_perturbed = s;
                for i = 1:s.mesh.Nchi
                    solution_perturbed.var.rho(i+1, s.shock.cell_indices(i, 1)-j+1) = solution_perturbed.var.rho(i+1, s.shock.cell_indices(i, 1)-j+1) + s.stability_analysis.perturbation_magnitude;
                end
                solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
                solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s.stability_analysis.perturbation_magnitude);
                [vector_indices_m_temp, vector_indices_n_temp, values_temp] = GET_OUTPUT_INDEXING_A21_ALL(solution_perturbed, j);
                vector_indices_m = cat(1, vector_indices_m, vector_indices_m_temp);
                vector_indices_n = cat(1, vector_indices_n, vector_indices_n_temp);
                values = cat(1, values, values_temp);
            end

            % --- rho_u -> shock speed ---
            for j = 1:2
                solution_perturbed = s;
                for i = 1:s.mesh.Nchi
                    solution_perturbed.var.rho_u(i+1, s.shock.cell_indices(i, 1)-j+1) = solution_perturbed.var.rho_u(i+1, s.shock.cell_indices(i, 1)-j+1) + s.stability_analysis.perturbation_magnitude;
                end
                solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
                solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s.stability_analysis.perturbation_magnitude);
                [vector_indices_m_temp, vector_indices_n_temp, values_temp] = GET_OUTPUT_INDEXING_A21_ALL(solution_perturbed, j);
                vector_indices_m = cat(1, vector_indices_m, vector_indices_m_temp);
                vector_indices_n = cat(1, vector_indices_n, vector_indices_n_temp + s.mesh.Nchi * s.mesh.Neta);
                values = cat(1, values, values_temp);
            end

            % --- rho_v -> shock speed ---
            for j = 1:2
                solution_perturbed = s;
                for i = 1:s.mesh.Nchi
                    solution_perturbed.var.rho_v(i+1, s.shock.cell_indices(i, 1)-j+1) = solution_perturbed.var.rho_v(i+1, s.shock.cell_indices(i, 1)-j+1) + s.stability_analysis.perturbation_magnitude;
                end
                solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
                solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s.stability_analysis.perturbation_magnitude);
                [vector_indices_m_temp, vector_indices_n_temp, values_temp] = GET_OUTPUT_INDEXING_A21_ALL(solution_perturbed, j);
                vector_indices_m = cat(1, vector_indices_m, vector_indices_m_temp);
                vector_indices_n = cat(1, vector_indices_n, vector_indices_n_temp + 2 * s.mesh.Nchi * s.mesh.Neta);
                values = cat(1, values, values_temp);
            end

            % --- rho_E -> shock speed ---
            for j = 1:2
                solution_perturbed = s;
                for i = 1:s.mesh.Nchi
                    solution_perturbed.var.rho_E(i+1, s.shock.cell_indices(i, 1)-j+1) = solution_perturbed.var.rho_E(i+1, s.shock.cell_indices(i, 1)-j+1) + s.stability_analysis.perturbation_magnitude;
                end
                solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
                solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s.stability_analysis.perturbation_magnitude);
                [vector_indices_m_temp, vector_indices_n_temp, values_temp] = GET_OUTPUT_INDEXING_A21_ALL(solution_perturbed, j);
                vector_indices_m = cat(1, vector_indices_m, vector_indices_m_temp);
                vector_indices_n = cat(1, vector_indices_n, vector_indices_n_temp + 3 * s.mesh.Nchi * s.mesh.Neta);
                values = cat(1, values, values_temp);
            end
        end
    end

    %% A12 and A22: Shock position perturbation -> flux and shock speed responses
    if s.stability_analysis.perturb_shock
        for shock_point = 1:s.mesh.Nchi
            solution_perturbed = s;
            solution_perturbed = PERTURB_SHOCK(solution_perturbed, s.stability_analysis.perturbation_magnitude, shock_point);
            solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
            solution_perturbed = REMOVE_BASE_OUTPUT(solution_perturbed, solution_new, s.stability_analysis.perturbation_magnitude);

            [vector_indices_m_temp, vector_indices_n_temp, values_temp] = GET_OUTPUT_INDEXING_A12_A22(solution_perturbed, shock_point);
            vector_indices_m = cat(1, vector_indices_m, vector_indices_m_temp);
            vector_indices_n = cat(1, vector_indices_n, vector_indices_n_temp + 4 * s.mesh.Nchi * s.mesh.Neta);
            values = cat(1, values, values_temp);
        end
    end

    %% Assemble sparse matrix
    if ~s.stability_analysis.perturb_shock
        A = sparse(vector_indices_m, vector_indices_n, values, ...
            4 * s.mesh.Nchi * s.mesh.Neta, ...
            4 * s.mesh.Nchi * s.mesh.Neta, ...
            4 * 4 * 9 * s.mesh.Nchi * s.mesh.Neta);
    else
        A = sparse(vector_indices_m, vector_indices_n, values, ...
            4 * s.mesh.Nchi * s.mesh.Neta + s.mesh.Nchi, ...
            4 * s.mesh.Nchi * s.mesh.Neta + s.mesh.Nchi, ...
            4 * 4 * 9 * s.mesh.Nchi * s.mesh.Neta + 13 * s.mesh.Nchi * s.mesh.Nchi);
    end
end


%% ========================================================================
%  A11 Block: Stencil-Based Output Indexing
%  ========================================================================
function [vector_indices_m, vector_indices_n, values] = GET_OUTPUT_INDEXING_A11(s, matrix_indices_2D)
% GET_OUTPUT_INDEXING_A11  Map perturbation-response pairs into sparse
%   matrix triplets for the A11 block (flow -> flow).
%
%   Uses a 3x3 stencil (9-point):
%       x x x
%       x o x
%       x x x
%
%   Handles periodic boundary conditions

    N = size(matrix_indices_2D, 1);

    % Pre-allocate: 4 flux outputs x 9 stencil cells x N perturbed cells
    vector_indices_n = zeros(4 * 9 * N, 1);
    vector_indices_m = zeros(4 * 9 * N, 1);
    values = zeros(4 * 9 * N, 1);

    % --- flux_rho ---
    for ele = 1:N
        i = matrix_indices_2D(ele, 1);
        j = matrix_indices_2D(ele, 2);
        index = 9 * (ele - 1);

        % Center cell
        vector_indices_n(index+1) = i + (j-1) * s.mesh.Nchi;
        vector_indices_m(index+1) = i + (j-1) * s.mesh.Nchi;
        values(index+1) = s.flux.rho(i, j);

        if i > 1  % East cell
            vector_indices_n(index+2) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+2) = (i-1) + (j-1) * s.mesh.Nchi;
            values(index+2) = s.flux.rho(i-1, j);
        end
        if i < s.mesh.Nchi  % West cell
            vector_indices_n(index+3) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+3) = (i+1) + (j-1) * s.mesh.Nchi;
            values(index+3) = s.flux.rho(i+1, j);
        end
        if j > 1  % South cell
            vector_indices_n(index+4) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+4) = i + (j-2) * s.mesh.Nchi;
            values(index+4) = s.flux.rho(i, j-1);
        end
        if j < s.mesh.Neta  % North cell
            vector_indices_n(index+5) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+5) = i + (j) * s.mesh.Nchi;
            values(index+5) = s.flux.rho(i, j+1);
        end
        if i > 1 && j < s.mesh.Neta  % North-east cell
            vector_indices_n(index+6) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+6) = (i-1) + (j) * s.mesh.Nchi;
            values(index+6) = s.flux.rho(i-1, j+1);
        end
        if i > 1 && j > 1  % South-east cell
            vector_indices_n(index+7) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+7) = (i-1) + (j-2) * s.mesh.Nchi;
            values(index+7) = s.flux.rho(i-1, j-1);
        end
        if i < s.mesh.Nchi && j > 1  % South-west cell
            vector_indices_n(index+8) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+8) = (i+1) + (j-2) * s.mesh.Nchi;
            values(index+8) = s.flux.rho(i+1, j-1);
        end
        if i < s.mesh.Nchi && j < s.mesh.Neta  % North-west cell
            vector_indices_n(index+9) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+9) = (i+1) + (j) * s.mesh.Nchi;
            values(index+9) = s.flux.rho(i+1, j+1);
        end

        % Periodic boundary corrections
        if s.boundary_conditions.boundary_chi0.name == "periodic" && s.boundary_conditions.boundary_chi1.name == "periodic"
            if i == 1
                vector_indices_n(index+2) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+2) = (i-1) + (j-1) * s.mesh.Nchi + s.mesh.Nchi;
                values(index+2) = s.flux.rho(s.mesh.Nchi, j);
                if j < s.mesh.Neta
                    vector_indices_n(index+6) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+6) = (i-1) + (j) * s.mesh.Nchi + s.mesh.Nchi;
                    values(index+6) = s.flux.rho(s.mesh.Nchi, j+1);
                end
                if j > 1
                    vector_indices_n(index+7) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+7) = (i-1) + (j-2) * s.mesh.Nchi + s.mesh.Nchi;
                    values(index+7) = s.flux.rho(s.mesh.Nchi, j-1);
                end
            end
            if i == s.mesh.Nchi
                vector_indices_n(index+3) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+3) = (i+1) + (j-1) * s.mesh.Nchi - s.mesh.Nchi;
                values(index+3) = s.flux.rho(1, j);
                if j > 1
                    vector_indices_n(index+8) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+8) = (i+1) + (j-2) * s.mesh.Nchi - s.mesh.Nchi;
                    values(index+8) = s.flux.rho(1, j-1);
                end
                if j < s.mesh.Neta
                    vector_indices_n(index+9) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+9) = (i+1) + (j) * s.mesh.Nchi - s.mesh.Nchi;
                    values(index+9) = s.flux.rho(1, j+1);
                end
            end
        end
        if s.boundary_conditions.boundary_eta0.name == "periodic" && s.boundary_conditions.boundary_eta1.name == "periodic"
            if j == 1
                vector_indices_n(index+4) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+4) = i + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
                values(index+4) = s.flux.rho(i, s.mesh.Neta);
                if i > 1
                    vector_indices_n(index+7) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+7) = (i-1) + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
                    values(index+7) = s.flux.rho(i-1, s.mesh.Neta);
                end
                if i < s.mesh.Nchi
                    vector_indices_n(index+8) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+8) = (i+1) + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
                    values(index+8) = s.flux.rho(i+1, s.mesh.Neta);
                end
            end
            if j == s.mesh.Neta
                vector_indices_n(index+5) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+5) = i + (j) * s.mesh.Nchi - s.mesh.Nchi * s.mesh.Neta;
                values(index+5) = s.flux.rho(i, 1);
                if i > 1
                    vector_indices_n(index+6) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+6) = (i-1) + (j) * s.mesh.Nchi - s.mesh.Nchi * s.mesh.Neta;
                    values(index+6) = s.flux.rho(i-1, 1);
                end
                if i < s.mesh.Nchi
                    vector_indices_n(index+9) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+9) = (i+1) + (j) * s.mesh.Nchi - s.mesh.Nchi * s.mesh.Neta;
                    values(index+9) = s.flux.rho(i+1, 1);
                end
            end
        end
    end

    % --- flux_rho_u ---
    for ele = 1:N
        i = matrix_indices_2D(ele, 1);
        j = matrix_indices_2D(ele, 2);
        index = 9 * (ele - 1) + 9 * N;

        vector_indices_n(index+1) = i + (j-1) * s.mesh.Nchi;
        vector_indices_m(index+1) = i + (j-1) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
        values(index+1) = s.flux.rho_u(i, j);

        if i > 1
            vector_indices_n(index+2) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+2) = (i-1) + (j-1) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
            values(index+2) = s.flux.rho_u(i-1, j);
        end
        if i < s.mesh.Nchi
            vector_indices_n(index+3) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+3) = (i+1) + (j-1) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
            values(index+3) = s.flux.rho_u(i+1, j);
        end
        if j > 1
            vector_indices_n(index+4) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+4) = i + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
            values(index+4) = s.flux.rho_u(i, j-1);
        end
        if j < s.mesh.Neta
            vector_indices_n(index+5) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+5) = i + (j) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
            values(index+5) = s.flux.rho_u(i, j+1);
        end
        if i > 1 && j < s.mesh.Neta
            vector_indices_n(index+6) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+6) = (i-1) + (j) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
            values(index+6) = s.flux.rho_u(i-1, j+1);
        end
        if i > 1 && j > 1
            vector_indices_n(index+7) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+7) = (i-1) + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
            values(index+7) = s.flux.rho_u(i-1, j-1);
        end
        if i < s.mesh.Nchi && j > 1
            vector_indices_n(index+8) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+8) = (i+1) + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
            values(index+8) = s.flux.rho_u(i+1, j-1);
        end
        if i < s.mesh.Nchi && j < s.mesh.Neta
            vector_indices_n(index+9) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+9) = (i+1) + (j) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
            values(index+9) = s.flux.rho_u(i+1, j+1);
        end

        % Periodic boundary corrections
        if s.boundary_conditions.boundary_chi0.name == "periodic" && s.boundary_conditions.boundary_chi1.name == "periodic"
            if i == 1
                vector_indices_n(index+2) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+2) = (i-1) + (j-1) * s.mesh.Nchi + s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
                values(index+2) = s.flux.rho_u(s.mesh.Nchi, j);
                if j < s.mesh.Neta
                    vector_indices_n(index+6) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+6) = (i-1) + (j) * s.mesh.Nchi + s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
                    values(index+6) = s.flux.rho_u(s.mesh.Nchi, j+1);
                end
                if j > 1
                    vector_indices_n(index+7) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+7) = (i-1) + (j-2) * s.mesh.Nchi + s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
                    values(index+7) = s.flux.rho_u(s.mesh.Nchi, j-1);
                end
            end
            if i == s.mesh.Nchi
                vector_indices_n(index+3) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+3) = (i+1) + (j-1) * s.mesh.Nchi - s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
                values(index+3) = s.flux.rho_u(1, j);
                if j > 1
                    vector_indices_n(index+8) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+8) = (i+1) + (j-2) * s.mesh.Nchi - s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
                    values(index+8) = s.flux.rho_u(1, j-1);
                end
                if j < s.mesh.Neta
                    vector_indices_n(index+9) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+9) = (i+1) + (j) * s.mesh.Nchi - s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
                    values(index+9) = s.flux.rho_u(1, j+1);
                end
            end
        end
        if s.boundary_conditions.boundary_eta0.name == "periodic" && s.boundary_conditions.boundary_eta1.name == "periodic"
            if j == 1
                vector_indices_n(index+4) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+4) = i + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta + s.mesh.Nchi * s.mesh.Neta;
                values(index+4) = s.flux.rho_u(i, s.mesh.Neta);
                if i > 1
                    vector_indices_n(index+7) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+7) = (i-1) + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta + s.mesh.Nchi * s.mesh.Neta;
                    values(index+7) = s.flux.rho_u(i-1, s.mesh.Neta);
                end
                if i < s.mesh.Nchi
                    vector_indices_n(index+8) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+8) = (i+1) + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta + s.mesh.Nchi * s.mesh.Neta;
                    values(index+8) = s.flux.rho_u(i+1, s.mesh.Neta);
                end
            end
            if j == s.mesh.Neta
                vector_indices_n(index+5) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+5) = i + (j) * s.mesh.Nchi - s.mesh.Nchi * s.mesh.Neta + s.mesh.Nchi * s.mesh.Neta;
                values(index+5) = s.flux.rho_u(i, 1);
                if i > 1
                    vector_indices_n(index+6) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+6) = (i-1) + (j) * s.mesh.Nchi - s.mesh.Nchi * s.mesh.Neta + s.mesh.Nchi * s.mesh.Neta;
                    values(index+6) = s.flux.rho_u(i-1, 1);
                end
                if i < s.mesh.Nchi
                    vector_indices_n(index+9) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+9) = (i+1) + (j) * s.mesh.Nchi - s.mesh.Nchi * s.mesh.Neta + s.mesh.Nchi * s.mesh.Neta;
                    values(index+9) = s.flux.rho_u(i+1, 1);
                end
            end
        end
    end

    % --- flux_rho_v ---
    for ele = 1:N
        i = matrix_indices_2D(ele, 1);
        j = matrix_indices_2D(ele, 2);
        index = 9 * (ele - 1) + 2 * 9 * N;

        vector_indices_n(index+1) = i + (j-1) * s.mesh.Nchi;
        vector_indices_m(index+1) = i + (j-1) * s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
        values(index+1) = s.flux.rho_v(i, j);

        if i > 1
            vector_indices_n(index+2) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+2) = (i-1) + (j-1) * s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
            values(index+2) = s.flux.rho_v(i-1, j);
        end
        if i < s.mesh.Nchi
            vector_indices_n(index+3) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+3) = (i+1) + (j-1) * s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
            values(index+3) = s.flux.rho_v(i+1, j);
        end
        if j > 1
            vector_indices_n(index+4) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+4) = i + (j-2) * s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
            values(index+4) = s.flux.rho_v(i, j-1);
        end
        if j < s.mesh.Neta
            vector_indices_n(index+5) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+5) = i + (j) * s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
            values(index+5) = s.flux.rho_v(i, j+1);
        end
        if i > 1 && j < s.mesh.Neta
            vector_indices_n(index+6) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+6) = (i-1) + (j) * s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
            values(index+6) = s.flux.rho_v(i-1, j+1);
        end
        if i > 1 && j > 1
            vector_indices_n(index+7) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+7) = (i-1) + (j-2) * s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
            values(index+7) = s.flux.rho_v(i-1, j-1);
        end
        if i < s.mesh.Nchi && j > 1
            vector_indices_n(index+8) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+8) = (i+1) + (j-2) * s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
            values(index+8) = s.flux.rho_v(i+1, j-1);
        end
        if i < s.mesh.Nchi && j < s.mesh.Neta
            vector_indices_n(index+9) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+9) = (i+1) + (j) * s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
            values(index+9) = s.flux.rho_v(i+1, j+1);
        end

        % Periodic boundary corrections
        if s.boundary_conditions.boundary_chi0.name == "periodic" && s.boundary_conditions.boundary_chi1.name == "periodic"
            if i == 1
                vector_indices_n(index+2) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+2) = (i-1) + (j-1) * s.mesh.Nchi + s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
                values(index+2) = s.flux.rho_v(s.mesh.Nchi, j);
                if j < s.mesh.Neta
                    vector_indices_n(index+6) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+6) = (i-1) + (j) * s.mesh.Nchi + s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+6) = s.flux.rho_v(s.mesh.Nchi, j+1);
                end
                if j > 1
                    vector_indices_n(index+7) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+7) = (i-1) + (j-2) * s.mesh.Nchi + s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+7) = s.flux.rho_v(s.mesh.Nchi, j-1);
                end
            end
            if i == s.mesh.Nchi
                vector_indices_n(index+3) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+3) = (i+1) + (j-1) * s.mesh.Nchi - s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
                values(index+3) = s.flux.rho_v(1, j);
                if j > 1
                    vector_indices_n(index+8) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+8) = (i+1) + (j-2) * s.mesh.Nchi - s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+8) = s.flux.rho_v(1, j-1);
                end
                if j < s.mesh.Neta
                    vector_indices_n(index+9) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+9) = (i+1) + (j) * s.mesh.Nchi - s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+9) = s.flux.rho_v(1, j+1);
                end
            end
        end
        if s.boundary_conditions.boundary_eta0.name == "periodic" && s.boundary_conditions.boundary_eta1.name == "periodic"
            if j == 1
                vector_indices_n(index+4) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+4) = i + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta + 2 * s.mesh.Nchi * s.mesh.Neta;
                values(index+4) = s.flux.rho_v(i, s.mesh.Neta);
                if i > 1
                    vector_indices_n(index+7) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+7) = (i-1) + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta + 2 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+7) = s.flux.rho_v(i-1, s.mesh.Neta);
                end
                if i < s.mesh.Nchi
                    vector_indices_n(index+8) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+8) = (i+1) + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta + 2 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+8) = s.flux.rho_v(i+1, s.mesh.Neta);
                end
            end
            if j == s.mesh.Neta
                vector_indices_n(index+5) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+5) = i + (j) * s.mesh.Nchi - s.mesh.Nchi * s.mesh.Neta + 2 * s.mesh.Nchi * s.mesh.Neta;
                values(index+5) = s.flux.rho_v(i, 1);
                if i > 1
                    vector_indices_n(index+6) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+6) = (i-1) + (j) * s.mesh.Nchi - s.mesh.Nchi * s.mesh.Neta + 2 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+6) = s.flux.rho_v(i-1, 1);
                end
                if i < s.mesh.Nchi
                    vector_indices_n(index+9) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+9) = (i+1) + (j) * s.mesh.Nchi - s.mesh.Nchi * s.mesh.Neta + 2 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+9) = s.flux.rho_v(i+1, 1);
                end
            end
        end
    end

    % --- flux_rho_E ---
    for ele = 1:N
        i = matrix_indices_2D(ele, 1);
        j = matrix_indices_2D(ele, 2);
        index = 9 * (ele - 1) + 3 * 9 * N;

        vector_indices_n(index+1) = i + (j-1) * s.mesh.Nchi;
        vector_indices_m(index+1) = i + (j-1) * s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
        values(index+1) = s.flux.rho_E(i, j);

        if i > 1
            vector_indices_n(index+2) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+2) = (i-1) + (j-1) * s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
            values(index+2) = s.flux.rho_E(i-1, j);
        end
        if i < s.mesh.Nchi
            vector_indices_n(index+3) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+3) = (i+1) + (j-1) * s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
            values(index+3) = s.flux.rho_E(i+1, j);
        end
        if j > 1
            vector_indices_n(index+4) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+4) = i + (j-2) * s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
            values(index+4) = s.flux.rho_E(i, j-1);
        end
        if j < s.mesh.Neta
            vector_indices_n(index+5) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+5) = i + (j) * s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
            values(index+5) = s.flux.rho_E(i, j+1);
        end
        if i > 1 && j < s.mesh.Neta
            vector_indices_n(index+6) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+6) = (i-1) + (j) * s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
            values(index+6) = s.flux.rho_E(i-1, j+1);
        end
        if i > 1 && j > 1
            vector_indices_n(index+7) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+7) = (i-1) + (j-2) * s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
            values(index+7) = s.flux.rho_E(i-1, j-1);
        end
        if i < s.mesh.Nchi && j > 1
            vector_indices_n(index+8) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+8) = (i+1) + (j-2) * s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
            values(index+8) = s.flux.rho_E(i+1, j-1);
        end
        if i < s.mesh.Nchi && j < s.mesh.Neta
            vector_indices_n(index+9) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index+9) = (i+1) + (j) * s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
            values(index+9) = s.flux.rho_E(i+1, j+1);
        end

        % Periodic boundary corrections
        if s.boundary_conditions.boundary_chi0.name == "periodic" && s.boundary_conditions.boundary_chi1.name == "periodic"
            if i == 1
                vector_indices_n(index+2) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+2) = (i-1) + (j-1) * s.mesh.Nchi + s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
                values(index+2) = s.flux.rho_E(s.mesh.Nchi, j);
                if j < s.mesh.Neta
                    vector_indices_n(index+6) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+6) = (i-1) + (j) * s.mesh.Nchi + s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+6) = s.flux.rho_E(s.mesh.Nchi, j+1);
                end
                if j > 1
                    vector_indices_n(index+7) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+7) = (i-1) + (j-2) * s.mesh.Nchi + s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+7) = s.flux.rho_E(s.mesh.Nchi, j-1);
                end
            end
            if i == s.mesh.Nchi
                vector_indices_n(index+3) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+3) = (i+1) + (j-1) * s.mesh.Nchi - s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
                values(index+3) = s.flux.rho_E(1, j);
                if j > 1
                    vector_indices_n(index+8) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+8) = (i+1) + (j-2) * s.mesh.Nchi - s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+8) = s.flux.rho_E(1, j-1);
                end
                if j < s.mesh.Neta
                    vector_indices_n(index+9) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+9) = (i+1) + (j) * s.mesh.Nchi - s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+9) = s.flux.rho_E(1, j+1);
                end
            end
        end
        if s.boundary_conditions.boundary_eta0.name == "periodic" && s.boundary_conditions.boundary_eta1.name == "periodic"
            if j == 1
                vector_indices_n(index+4) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+4) = i + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta + 3 * s.mesh.Nchi * s.mesh.Neta;
                values(index+4) = s.flux.rho_E(i, s.mesh.Neta);
                if i > 1
                    vector_indices_n(index+7) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+7) = (i-1) + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta + 3 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+7) = s.flux.rho_E(i-1, s.mesh.Neta);
                end
                if i < s.mesh.Nchi
                    vector_indices_n(index+8) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+8) = (i+1) + (j-2) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta + 3 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+8) = s.flux.rho_E(i+1, s.mesh.Neta);
                end
            end
            if j == s.mesh.Neta
                vector_indices_n(index+5) = i + (j-1) * s.mesh.Nchi;
                vector_indices_m(index+5) = i + (j) * s.mesh.Nchi - s.mesh.Nchi * s.mesh.Neta + 3 * s.mesh.Nchi * s.mesh.Neta;
                values(index+5) = s.flux.rho_E(i, 1);
                if i > 1
                    vector_indices_n(index+6) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+6) = (i-1) + (j) * s.mesh.Nchi - s.mesh.Nchi * s.mesh.Neta + 3 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+6) = s.flux.rho_E(i-1, 1);
                end
                if i < s.mesh.Nchi
                    vector_indices_n(index+9) = i + (j-1) * s.mesh.Nchi;
                    vector_indices_m(index+9) = (i+1) + (j) * s.mesh.Nchi - s.mesh.Nchi * s.mesh.Neta + 3 * s.mesh.Nchi * s.mesh.Neta;
                    values(index+9) = s.flux.rho_E(i+1, 1);
                end
            end
        end
    end
end


%% ========================================================================
%  A21 Block: Per-Point Shock Speed Indexing
%  ========================================================================
function [vector_indices_m, vector_indices_n, values] = GET_OUTPUT_INDEXING_A21(s, i, j)
% GET_OUTPUT_INDEXING_A21  Map a single flow-field perturbation at (i,j) to
%   shock speed responses for all shock points (per-point strategy).

    if s.stability_analysis.perturb_shock
        vector_indices_n = zeros(s.mesh.Nchi, 1);
        vector_indices_m = zeros(s.mesh.Nchi, 1);
        values = zeros(s.mesh.Nchi, 1);

        for index = 1:s.mesh.Nchi
            vector_indices_n(index) = i + (j-1) * s.mesh.Nchi;
            vector_indices_m(index) = index + 4 * s.mesh.Nchi * s.mesh.Neta;
            values(index) = s.shock.relative_increase_velocity(index, 1);
        end
    end
end

%% ========================================================================
%  A21 Block: All-at-Once Shock Speed Indexing
%  ========================================================================
function [vector_indices_m, vector_indices_n, values] = GET_OUTPUT_INDEXING_A21_ALL(s, j)
% GET_OUTPUT_INDEXING_A21_ALL  Map simultaneous flow-field perturbations at
%   all shocked cells to shock speed responses (no-spline strategy).

    if s.stability_analysis.perturb_shock
        vector_indices_n = zeros(s.mesh.Nchi, 1);
        vector_indices_m = zeros(s.mesh.Nchi, 1);
        values = zeros(s.mesh.Nchi, 1);

        for index = 1:s.mesh.Nchi
            vector_indices_n(index) = index + (s.shock.cell_indices(index, 1)-j-1) * s.mesh.Nchi;
            vector_indices_m(index) = index + 4 * s.mesh.Nchi * s.mesh.Neta;
            values(index) = s.shock.relative_increase_velocity(index, 1);
        end
    end
end


%% ========================================================================
%  A12 & A22 Blocks: Shock Position Perturbation Indexing
%  ========================================================================
function [vector_indices_m, vector_indices_n, values] = GET_OUTPUT_INDEXING_A12_A22(s, shock_point)
% GET_OUTPUT_INDEXING_A12_A22  Map a shock position perturbation to flux
%   responses (A12) and shock speed responses (A22).
%
%   Dimensions per shock_point:
%       4*Nx entries for A12 (rho, rho_u, rho_v, rho_E fluxes)
%       Nx entries for A22 (shock speed)

    if s.stability_analysis.perturb_shock
        vector_indices_n = zeros(5 * s.mesh.Nchi, 1);
        vector_indices_m = zeros(5 * s.mesh.Nchi, 1);
        values = zeros(5 * s.mesh.Nchi, 1);

        for i = 1:s.mesh.Nchi
            j = s.shock.cell_indices(i, 1) - 1;

            % rho flux
            index = i;
            vector_indices_n(index) = shock_point;
            vector_indices_m(index) = i + (j-1) * s.mesh.Nchi;
            values(index) = s.flux.rho(i, j);

            % rho_u flux
            index = i + s.mesh.Nchi;
            vector_indices_n(index) = shock_point;
            vector_indices_m(index) = i + (j-1) * s.mesh.Nchi + s.mesh.Nchi * s.mesh.Neta;
            values(index) = s.flux.rho_u(i, j);

            % rho_v flux
            index = i + 2 * s.mesh.Nchi;
            vector_indices_n(index) = shock_point;
            vector_indices_m(index) = i + (j-1) * s.mesh.Nchi + 2 * s.mesh.Nchi * s.mesh.Neta;
            values(index) = s.flux.rho_v(i, j);

            % rho_E flux
            index = i + 3 * s.mesh.Nchi;
            vector_indices_n(index) = shock_point;
            vector_indices_m(index) = i + (j-1) * s.mesh.Nchi + 3 * s.mesh.Nchi * s.mesh.Neta;
            values(index) = s.flux.rho_E(i, j);

            % Shock speed
            index = i + 4 * s.mesh.Nchi;
            vector_indices_n(index) = shock_point;
            vector_indices_m(index) = i + 4 * s.mesh.Nchi * s.mesh.Neta;
            values(index) = s.shock.relative_increase_velocity(i, 1);
        end
    end
end

%% ========================================================================
%  Shock Perturbation Helper
%  ========================================================================
function s = PERTURB_SHOCK(s, perturbation_magnitude, shock_point)
% PERTURB_SHOCK  Displace a single shock point in the shock-normal direction.

    ang = s.shock.beta(shock_point, 1) - pi/2 - atan2(s.freestream.rho_v_0, s.freestream.rho_u_0);

    s.shock.points_x(shock_point, 1) = s.shock.points_x(shock_point, 1) + perturbation_magnitude * cos(ang);
    s.shock.points_y(shock_point, 1) = s.shock.points_y(shock_point, 1) + perturbation_magnitude * sin(ang);
end

%% ========================================================================
%  Base Flow Subtraction
%  ========================================================================
function solution_perturbation_only = REMOVE_BASE_OUTPUT(solution_perturbed, s, perturbation_magnitude)
% REMOVE_BASE_OUTPUT  Compute the finite-difference Jacobian column:
%   (f(x+Dx) - f(x)) / Dx for all flux outputs.

    solution_perturbation_only = s;
    solution_perturbation_only.flux.rho = (solution_perturbed.flux.rho - s.flux.rho) / perturbation_magnitude;
    solution_perturbation_only.flux.rho_u = (solution_perturbed.flux.rho_u - s.flux.rho_u) / perturbation_magnitude;
    solution_perturbation_only.flux.rho_v = (solution_perturbed.flux.rho_v - s.flux.rho_v) / perturbation_magnitude;
    solution_perturbation_only.flux.rho_E = (solution_perturbed.flux.rho_E - s.flux.rho_E) / perturbation_magnitude;
    if s.shock.enabled
        solution_perturbation_only.shock.relative_increase_velocity = (solution_perturbed.shock.relative_increase_velocity - s.shock.relative_increase_velocity) / perturbation_magnitude;
    end
end


%% ========================================================================
%  Perturbation Pattern Generator
%  ========================================================================
function [matrix_perturbation, vector_indices] = GET_PERTURBATION_MATRIX(index, s)
% GET_PERTURBATION_MATRIX  Generate a disjoint perturbation pattern for the
%   given index (1..9). Each pattern activates every 3rd cell in both i and
%   j directions, producing 9 non-overlapping patterns that together cover
%   the full domain without stencil interference.

    matrix_perturbation = zeros(s.mesh.Nchi, s.mesh.Neta);
    count = 1;
    jump = 3;

    if (mod(index, 3) == 0)
        start_i = 3;
    else
        start_i = mod(index, 3);
    end

    start_j = ceil(index / 3);

    for j = start_j:jump:s.mesh.Neta
        for i = start_i:jump:s.mesh.Nchi
            vector_indices(count, 1) = i;
            vector_indices(count, 2) = j;
            count = count + 1;
            matrix_perturbation(i, j) = 1;
        end
    end
end
