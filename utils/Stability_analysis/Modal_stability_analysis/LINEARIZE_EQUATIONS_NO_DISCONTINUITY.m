function A = LINEARIZE_EQUATIONS_NO_DISCONTINUITY(s, chemistry)
% LINEARIZE_EQUATIONS_NO_DISCONTINUITY  Build the sparse Jacobian matrix via
%   numerical finite differences of the nonlinear dynamics operator.
%
%   A = LINEARIZE_EQUATIONS_NO_DISCONTINUITY(s, chemistry) constructs
%   the linearized operator A = df/dq using column-wise finite differences.
%   The matrix maps perturbations in conservative variables (and optionally
%   shock position) to flux perturbations:
%
%       [flux_rho; flux_rho_u; flux_rho_v; flux_rho_E; (flux_rho_w); dx_shock] =
%           A * [rho; rho_u; rho_v; rho_E; (rho_w); x_shock]
%
%   The azimuthal momentum rho_w is present only for an azimuthal mode
%   (s.azimuthal, see SETUP_AZIMUTHAL_MODE); PERTURBATION_FIELDS lists the
%   variables in order.
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
%             (nv*Nx*Ny) x (nv*Nx*Ny)            without shock perturbation
%             (nv*Nx*Ny+Nx) x (nv*Nx*Ny+Nx)      with shock perturbation
%           where nv = numel(PERTURBATION_FIELDS(s)) (4, or 5 with rho_w).
%
%   Notes:
%       - The stencil has a 3x3 domain of influence (9 neighbors per cell).
%       - Nine disjoint perturbation patterns are swept to cover all cells
%         without overlap in the stencil.
%       - Supports periodic boundary conditions (chi and/or eta pairs).
%       - Shock perturbation blocks (A12, A21, A22) are only assembled when
%         s.stability_analysis.perturb_shock is true.
%       - Uses a spline_param switch to select between per-point and
%         all-at-once shock perturbation strategies (A21 block).
%       - rho_w does not enter the Rankine-Hugoniot normal jump, so its
%         A21 columns are zero and are not evaluated.
%       - The shock DoFs are upstream normal displacements; their rows are
%         the shock speed relative_increase_velocity (positive upstream).
%
% Part of: Hypersonics Stability MATLAB Solver - Stability Analysis / Linear Stability Module

    %% Evaluate base flow operator f(x)
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry);
    s_base = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s, chemistry);

    fields = PERTURBATION_FIELDS(s);
    nv = numel(fields);
    Nchi = s.mesh.Nchi;
    Neta = s.mesh.Neta;
    N = Nchi * Neta;
    n_dof = nv * N + s.stability_analysis.perturb_shock * Nchi;
    epsilon = s.stability_analysis.perturbation_magnitude;
    bc = s.boundary_conditions;
    periodic_chi = bc.boundary_chi0.name == "periodic" && bc.boundary_chi1.name == "periodic";
    periodic_eta = bc.boundary_eta0.name == "periodic" && bc.boundary_eta1.name == "periodic";

    triplets = cell(0, 3);

    %% A11: Flow field perturbations -> flux responses
    for k = 1:nv
        for colour = 1:9
            [matrix_perturbation, cells] = GET_PERTURBATION_MATRIX(colour, s);
            if s.shock.enabled
                matrix_perturbation = s.shock.flow_cells .* matrix_perturbation;
            end
            s_pert = s;
            s_pert.var.(fields(k))(2:end-1, 2:end-1) = s.var.(fields(k))(2:end-1, 2:end-1) ...
                + matrix_perturbation * epsilon;
            response = FLUX_RESPONSE(s_pert, s_base, fields, epsilon, chemistry);
            triplets(end+1, :) = STENCIL_ENTRIES(response, cells, k, Nchi, Neta, ...
                                                 periodic_chi, periodic_eta); %#ok<AGROW>
        end
    end

    if s.stability_analysis.perturb_shock
        cell_j = s.shock.cell_indices(:, 1);
        shock_rows = nv * N + (1:Nchi)';

        %% A21: Flow cells next to the shock -> shock speed
        for k = 1:4
            for jj = 1:2
                cols = (k - 1) * N + (1:Nchi)' + (cell_j - jj - 1) * Nchi;
                idx = sub2ind(size(s.var.rho), (2:Nchi+1)', cell_j - jj + 1);
                if s.shock.spline_param ~= 1
                    % Spline-coupled shock points: one cell per evaluation
                    for i = 1:Nchi
                        s_pert = s;
                        s_pert.var.(fields(k))(idx(i)) = s_pert.var.(fields(k))(idx(i)) + epsilon;
                        response = FLUX_RESPONSE(s_pert, s_base, fields, epsilon, chemistry);
                        triplets(end+1, :) = {shock_rows, repmat(cols(i), Nchi, 1), ...
                                              response.shock_speed}; %#ok<AGROW>
                    end
                else
                    % Uncoupled shock points: perturb the whole row at once
                    s_pert = s;
                    s_pert.var.(fields(k))(idx) = s_pert.var.(fields(k))(idx) + epsilon;
                    response = FLUX_RESPONSE(s_pert, s_base, fields, epsilon, chemistry);
                    triplets(end+1, :) = {shock_rows, cols, response.shock_speed}; %#ok<AGROW>
                end
            end
        end

        %% A12 / A22: Shock point perturbations -> post-shock cells and shock speed
        cells_below = (1:Nchi)' + (cell_j - 2) * Nchi;
        for shock_point = 1:Nchi
            s_pert = PERTURB_SHOCK(s, epsilon, shock_point);
            response = FLUX_RESPONSE(s_pert, s_base, fields, epsilon, chemistry);
            col = repmat(nv * N + shock_point, Nchi, 1);
            for o = 1:nv
                triplets(end+1, :) = {(o - 1) * N + cells_below, col, ...
                                      response.flux{o}(cells_below)}; %#ok<AGROW>
            end
            triplets(end+1, :) = {shock_rows, col, response.shock_speed}; %#ok<AGROW>
        end
    end

    A = sparse(vertcat(triplets{:, 1}), vertcat(triplets{:, 2}), vertcat(triplets{:, 3}), ...
               n_dof, n_dof);
end

%% ========================================================================
%  Stencil Triplets of One Colour Sweep
%  ========================================================================
function entries = STENCIL_ENTRIES(response, cells, k, Nchi, Neta, periodic_chi, periodic_eta)
% STENCIL_ENTRIES  Sparse triplets {rows, cols, values} coupling every
%   perturbed cell (a column of variable k) to the 3x3 neighbourhood of
%   every flux field. Neighbours outside the grid are dropped unless the
%   direction is periodic, in which case they wrap around.

    N = Nchi * Neta;
    nv = numel(response.flux);
    [di, dj] = ndgrid(-1:1, -1:1);
    ti = cells(:, 1) + di(:)';        % (n_cells x 9) neighbour indices
    tj = cells(:, 2) + dj(:)';
    if periodic_chi
        ti = mod(ti - 1, Nchi) + 1;
    end
    if periodic_eta
        tj = mod(tj - 1, Neta) + 1;
    end
    valid = ti >= 1 & ti <= Nchi & tj >= 1 & tj <= Neta;
    target = ti(valid) + (tj(valid) - 1) * Nchi;
    col = repmat((k - 1) * N + cells(:, 1) + (cells(:, 2) - 1) * Nchi, 1, 9);
    col = col(valid);

    rows = cell(nv, 1);
    vals = cell(nv, 1);
    for o = 1:nv
        rows{o} = (o - 1) * N + target;
        vals{o} = response.flux{o}(target);
    end
    entries = {vertcat(rows{:}), repmat(col, nv, 1), vertcat(vals{:})};
end

%% ========================================================================
%  Shock Perturbation Helper
%  ========================================================================
function s = PERTURB_SHOCK(s, perturbation_magnitude, shock_point)
% PERTURB_SHOCK  Displace a single shock point upstream, normal to the shock.
%   The shock DoF is the upstream normal displacement, whose rate is the
%   shock speed relative_increase_velocity (positive upstream, MOVE_SHOCK).

    ang = s.shock.beta(shock_point, 1) + pi/2 - atan2(s.freestream.rho_v_0, s.freestream.rho_u_0);

    s.shock.points_x(shock_point, 1) = s.shock.points_x(shock_point, 1) + perturbation_magnitude * cos(ang);
    s.shock.points_y(shock_point, 1) = s.shock.points_y(shock_point, 1) + perturbation_magnitude * sin(ang);
end

%% ========================================================================
%  Perturbation Pattern Generator
%  ========================================================================
function [matrix_perturbation, vector_indices] = GET_PERTURBATION_MATRIX(index, s)
% GET_PERTURBATION_MATRIX  Generate a disjoint perturbation pattern for the
%   given index (1..9). Each pattern activates every 3rd cell in both i and
%   j directions, producing 9 non-overlapping patterns that together cover
%   the full domain without stencil interference.

    start_i = mod(index - 1, 3) + 1;
    start_j = ceil(index / 3);
    [ii, jj] = ndgrid(start_i:3:s.mesh.Nchi, start_j:3:s.mesh.Neta);
    vector_indices = [ii(:), jj(:)];
    matrix_perturbation = zeros(s.mesh.Nchi, s.mesh.Neta);
    matrix_perturbation(sub2ind(size(matrix_perturbation), ii(:), jj(:))) = 1;
end
