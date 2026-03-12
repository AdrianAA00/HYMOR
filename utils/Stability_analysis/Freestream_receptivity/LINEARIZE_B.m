function B = LINEARIZE_B(s, perturbation, chemistry)
% LINEARIZE_B  Compute the linearized upstream boundary operator B via finite differences.
%
%   B = LINEARIZE_B(s, perturbation, chemistry)
%
%   Constructs the Jacobian matrix B that maps upstream (freestream)
%   perturbations in [rho, rho*u, rho*v, rho*E] to flux perturbations in
%   the downstream domain. The linearization uses column-wise finite
%   differences with a 3-color stencil to account for the tridiagonal
%   coupling along the shock front.
%
%   Inputs:
%       s            - Solution structure with base flow, mesh, and parameters
%       perturbation - Finite-difference step size for Jacobian approximation
%       chemistry    - Chemistry model structure
%
%   Outputs:
%       B - Sparse Jacobian matrix mapping upstream state to downstream fluxes.
%           Size is (4*Nx*Ny x 4*Nx) or (4*Nx*Ny+Nx x 4*Nx) when shock
%           perturbation DoFs are included.
%
%   Notes:
%       - Uses a 3-color finite-difference scheme to reduce the number of
%         nonlinear evaluations from O(Nx) to O(1) per variable.
%       - Contains local helper functions GET_BLOCK and REMOVE_BASE_OUTPUT.
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Freestream Module

    %% Initialize sparse triplet storage
    values = 0;
    vector_indices_m = 0;
    vector_indices_n = 0;

    %% Compute base flow
    s.freestream.disturbance.amplitude = [0, 0, 0, 0];
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry);

    s.freestream.rho_0_p = s.freestream.rho_0 * ones(s.mesh.Nchi, 1);
    s.freestream.rho_u_0_p = s.freestream.rho_u_0 * ones(s.mesh.Nchi, 1);
    s.freestream.rho_v_0_p = s.freestream.rho_v_0 * ones(s.mesh.Nchi, 1);
    s.freestream.rho_E_0_p = s.freestream.rho_E_0 * ones(s.mesh.Nchi, 1);
    s.rho_0_upstream_p = true;

    s_new = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s, chemistry);

    %% Linearize rho column block
    count_values = 0;
    count_n = 0;

    for iter = 1:3
        s_perturbed = s;
        pert_val = zeros(s.mesh.Nchi, 1);
        pert_val(iter:3:end) = 1;
        s_perturbed.freestream.rho_0_p = s.freestream.rho_0_p + perturbation * pert_val;

        % West cell
        if iter == 1
            pert_ind_left = (iter+2:3:s.mesh.Nchi-1)';
            pert_ind = pert_ind_left + 1;
        else
            pert_ind_left = (iter-1:3:s.mesh.Nchi-1)';
            pert_ind = pert_ind_left + 1;
        end
        [values, vector_indices_n, vector_indices_m, count_values, count_n] = ...
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_left, values, vector_indices_n, vector_indices_m, chemistry);

        % Mid cell
        pert_ind_mid = (iter:3:s.mesh.Nchi)';
        pert_ind = pert_ind_mid;
        [values, vector_indices_n, vector_indices_m, count_values, count_n] = ...
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_mid, values, vector_indices_n, vector_indices_m, chemistry);

        % East cell
        pert_ind_right = (iter+1:3:s.mesh.Nchi)';
        pert_ind = pert_ind_right - 1;
        [values, vector_indices_n, vector_indices_m, count_values, count_n] = ...
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_right, values, vector_indices_n, vector_indices_m, chemistry);
    end

    %% Linearize rho_u column block
    count_n = count_n + s.mesh.Nchi;

    for iter = 1:3
        s_perturbed = s;
        pert_val = zeros(s.mesh.Nchi, 1);
        pert_val(iter:3:end) = 1;
        s_perturbed.freestream.rho_u_0_p = s.freestream.rho_u_0_p + perturbation * pert_val;

        % West cell
        if iter == 1
            pert_ind_left = (iter+2:3:s.mesh.Nchi-1)';
            pert_ind = pert_ind_left + 1;
        else
            pert_ind_left = (iter-1:3:s.mesh.Nchi-1)';
            pert_ind = pert_ind_left + 1;
        end
        [values, vector_indices_n, vector_indices_m, count_values, count_n] = ...
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_left, values, vector_indices_n, vector_indices_m, chemistry);

        % Mid cell
        pert_ind_mid = (iter:3:s.mesh.Nchi)';
        pert_ind = pert_ind_mid;
        [values, vector_indices_n, vector_indices_m, count_values, count_n] = ...
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_mid, values, vector_indices_n, vector_indices_m, chemistry);

        % East cell
        pert_ind_right = (iter+1:3:s.mesh.Nchi)';
        pert_ind = pert_ind_right - 1;
        [values, vector_indices_n, vector_indices_m, count_values, count_n] = ...
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_right, values, vector_indices_n, vector_indices_m, chemistry);
    end

    %% Linearize rho_v column block
    count_n = count_n + s.mesh.Nchi;

    for iter = 1:3
        s_perturbed = s;
        pert_val = zeros(s.mesh.Nchi, 1);
        pert_val(iter:3:end) = 1;
        s_perturbed.freestream.rho_v_0_p = s.freestream.rho_v_0_p + perturbation * pert_val;

        % West cell
        if iter == 1
            pert_ind_left = (iter+2:3:s.mesh.Nchi-1)';
            pert_ind = pert_ind_left + 1;
        else
            pert_ind_left = (iter-1:3:s.mesh.Nchi-1)';
            pert_ind = pert_ind_left + 1;
        end
        [values, vector_indices_n, vector_indices_m, count_values, count_n] = ...
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_left, values, vector_indices_n, vector_indices_m, chemistry);

        % Mid cell
        pert_ind_mid = (iter:3:s.mesh.Nchi)';
        pert_ind = pert_ind_mid;
        [values, vector_indices_n, vector_indices_m, count_values, count_n] = ...
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_mid, values, vector_indices_n, vector_indices_m, chemistry);

        % East cell
        pert_ind_right = (iter+1:3:s.mesh.Nchi)';
        pert_ind = pert_ind_right - 1;
        [values, vector_indices_n, vector_indices_m, count_values, count_n] = ...
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_right, values, vector_indices_n, vector_indices_m, chemistry);
    end

    %% Linearize rho_E column block
    count_n = count_n + s.mesh.Nchi;

    for iter = 1:3
        s_perturbed = s;
        pert_val = zeros(s.mesh.Nchi, 1);
        pert_val(iter:3:end) = 1;
        s_perturbed.freestream.rho_E_0_p = s.freestream.rho_E_0_p + perturbation * pert_val;

        % West cell
        if iter == 1
            pert_ind_left = (iter+2:3:s.mesh.Nchi-1)';
            pert_ind = pert_ind_left + 1;
        else
            pert_ind_left = (iter-1:3:s.mesh.Nchi-1)';
            pert_ind = pert_ind_left + 1;
        end
        [values, vector_indices_n, vector_indices_m, count_values, count_n] = ...
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_left, values, vector_indices_n, vector_indices_m, chemistry);

        % Mid cell
        pert_ind_mid = (iter:3:s.mesh.Nchi)';
        pert_ind = pert_ind_mid;
        [values, vector_indices_n, vector_indices_m, count_values, count_n] = ...
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_mid, values, vector_indices_n, vector_indices_m, chemistry);

        % East cell
        pert_ind_right = (iter+1:3:s.mesh.Nchi)';
        pert_ind = pert_ind_right - 1;
        [values, vector_indices_n, vector_indices_m, count_values, count_n] = ...
            GET_BLOCK(s_perturbed, s_new, count_values, count_n, perturbation, pert_ind, pert_ind_right, values, vector_indices_n, vector_indices_m, chemistry);
    end

    %% Assemble sparse matrix
    if ~s.stability_analysis.perturb_shock
        B = sparse(vector_indices_m, vector_indices_n, values, 4*s.mesh.Nchi*s.mesh.Neta, 4*s.mesh.Nchi, size(vector_indices_m, 2));
    else
        B = sparse(vector_indices_m, vector_indices_n, values, 4*s.mesh.Nchi*s.mesh.Neta + s.mesh.Nchi, size(vector_indices_m, 2));
    end
end


function [values, vector_indices_n, vector_indices_m, count_values, count_n] = GET_BLOCK(s_perturbed, s, count_values, count_n, perturbation, pert_n, pert_m, values, vector_indices_n, vector_indices_m, chemistry)
% GET_BLOCK  Compute one finite-difference block of the B Jacobian.
%
%   Evaluates the perturbed nonlinear dynamics and extracts the Jacobian
%   entries for the specified perturbation and response indices.
%
%   Inputs:
%       s_perturbed  - Perturbed s structure
%       s            - Base s structure (with pre-computed fluxes)
%       count_values - Running counter for sparse triplet storage
%       count_n      - Column offset for the current variable block
%       perturbation - Finite-difference step size
%       pert_n       - Column indices (perturbation source cells)
%       pert_m       - Row indices (response cells)
%       values, vector_indices_n, vector_indices_m - Sparse triplet arrays
%       chemistry    - Chemistry model structure
%
%   Outputs:
%       Updated triplet arrays and counters

    count_m = 0;
    s_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s_perturbed, chemistry);
    s_perturbed = REMOVE_BASE_OUTPUT(s_perturbed, s, perturbation);
    size_ind = size(pert_m, 1);
    ind = sub2ind(size(s_perturbed.flux.rho), pert_m, s.shock.cell_indices(pert_m) - 1);

    %% rho flux block
    values(1+count_values:size_ind+count_values) = s_perturbed.flux.rho(ind);
    vector_indices_m(1+count_values:size_ind+count_values) = (s.shock.cell_indices(pert_m) - 2) * s.mesh.Nchi + pert_m + count_m;
    vector_indices_n(1+count_values:size_ind+count_values) = pert_n + count_n;
    count_values = count_values + size_ind;
    count_m = count_m + s.mesh.Nchi * s.mesh.Neta;

    %% rho_u flux block
    values(1+count_values:size_ind+count_values) = s_perturbed.flux.rho_u(ind);
    vector_indices_m(1+count_values:size_ind+count_values) = (s.shock.cell_indices(pert_m) - 2) * s.mesh.Nchi + pert_m + count_m;
    vector_indices_n(1+count_values:size_ind+count_values) = pert_n + count_n;
    count_values = count_values + size_ind;
    count_m = count_m + s.mesh.Nchi * s.mesh.Neta;

    %% rho_v flux block
    values(1+count_values:size_ind+count_values) = s_perturbed.flux.rho_v(ind);
    vector_indices_m(1+count_values:size_ind+count_values) = (s.shock.cell_indices(pert_m) - 2) * s.mesh.Nchi + pert_m + count_m;
    vector_indices_n(1+count_values:size_ind+count_values) = pert_n + count_n;
    count_values = count_values + size_ind;
    count_m = count_m + s.mesh.Nchi * s.mesh.Neta;

    %% rho_E flux block
    values(1+count_values:size_ind+count_values) = s_perturbed.flux.rho_E(ind);
    vector_indices_m(1+count_values:size_ind+count_values) = (s.shock.cell_indices(pert_m) - 2) * s.mesh.Nchi + pert_m + count_m;
    vector_indices_n(1+count_values:size_ind+count_values) = pert_n + count_n;
    count_values = count_values + size_ind;
    count_m = count_m + s.mesh.Nchi * s.mesh.Neta;

    %% Shock velocity block (if enabled)
    if s.stability_analysis.perturb_shock
        values(1+count_values:size_ind+count_values) = s_perturbed.shock.relative_increase_velocity(pert_m, 1);
        vector_indices_m(1+count_values:size_ind+count_values) = pert_m + count_m;
        vector_indices_n(1+count_values:size_ind+count_values) = pert_n + count_n;
        count_values = count_values + size_ind;
        count_m = count_m + s.mesh.Nchi;
    end
end


function solution_perturbation_only = REMOVE_BASE_OUTPUT(solution_perturbed, s, perturbation_magnitude)
% REMOVE_BASE_OUTPUT  Extract the linearized flux perturbation from the nonlinear evaluation.
%
%   solution_perturbation_only = REMOVE_BASE_OUTPUT(solution_perturbed, s, perturbation_magnitude)
%
%   Computes (f(x + dx) - f(x)) / dx for each flux component.
%
%   Inputs:
%       solution_perturbed    - Solution after perturbed nonlinear evaluation
%       s              - Solution after base nonlinear evaluation
%       perturbation_magnitude - Finite-difference step size
%
%   Outputs:
%       solution_perturbation_only - Structure with linearized flux fields

    solution_perturbation_only = s;
    solution_perturbation_only.flux.rho = (solution_perturbed.flux.rho - s.flux.rho) / perturbation_magnitude;
    solution_perturbation_only.flux.rho_u = (solution_perturbed.flux.rho_u - s.flux.rho_u) / perturbation_magnitude;
    solution_perturbation_only.flux.rho_v = (solution_perturbed.flux.rho_v - s.flux.rho_v) / perturbation_magnitude;
    solution_perturbation_only.flux.rho_E = (solution_perturbed.flux.rho_E - s.flux.rho_E) / perturbation_magnitude;
    if s.shock.enabled
        solution_perturbation_only.shock.relative_increase_velocity = (solution_perturbed.shock.relative_increase_velocity - s.shock.relative_increase_velocity) / perturbation_magnitude;
    end
end
