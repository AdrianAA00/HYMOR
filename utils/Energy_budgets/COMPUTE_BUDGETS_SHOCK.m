function budgets_shock = COMPUTE_BUDGETS_SHOCK(V, w_infty, s, chemistry)
% COMPUTE_BUDGETS_SHOCK  Compute energy budget contributions at the shock interface.
%
%   budgets_shock = COMPUTE_BUDGETS_SHOCK(V, w_infty, s, chemistry)
%
%   Evaluates the perturbation energy fluxes through the shock surface by
%   applying the Rankine-Hugoniot jump conditions to the perturbed upstream
%   state, computing the resulting post-shock perturbation, and integrating
%   the Chu energy flux across the shock. The upstream perturbation may
%   consist of multiple frequency modes that are summed before the jump
%   calculation.
%
%   Inputs:
%       V         - (double column vector) Global perturbation state vector.
%                   The first 4*Nx*Ny entries are interior perturbations; the
%                   remaining entries contain upstream (freestream) modal
%                   perturbations arranged as N_f blocks of 4*Nx entries each.
%       w_infty   - (double vector) Freestream modal frequency/weight vector
%                   of length N_f (number of frequency modes).
%       s  - (struct) Solution structure with base-flow fields, grid
%                   data (Nx, Ny, x, x_Ext, shocked_cell_indices, bt_area,
%                   bt_y_normal, shock_beta), and upstream reference states
%                   (rho_0_upstream, rho_u_0_upstream, etc.).
%       chemistry - (struct) Chemistry model passed to
%                   UPDATE_SHOCK_JUMP_PROPERTIES for equation-of-state
%                   evaluation across the shock.
%
%   Outputs:
%       budgets_shock - (struct) Shock-surface integrated energy budget:
%           .adv_acoustic  - Advective acoustic energy flux through shock.
%           .adv_kinetic   - Advective kinetic energy flux through shock.
%           .work_kinetic  - Pressure work (p' * u'_n) at the shock.
%           .adv_entropic  - Advective entropic energy flux through shock.
%
%   Notes:
%       - The upstream perturbation is obtained by summing all N_f frequency
%         modes and adding them to the base upstream state before calling
%         UPDATE_SHOCK_JUMP_PROPERTIES to compute perturbed post-shock
%         conditions via the Rankine-Hugoniot relations.
%       - Post-shock energy norms are computed with GET_ENERGY_NORM (Chu norm).
%       - The velocity perturbation normal to the shock is projected using
%         the local shock angle (shock_beta).
%       - Requires helper functions: UPDATE_SHOCK_JUMP_PROPERTIES,
%         GET_ENERGY_NORM, GET_CHU_COMPONENTS.
%
%   References:
%       Chu, B.-T. (1965), "On the energy transfer to small disturbances in
%       fluid flow (Part I)", Acta Mechanica, 1(3), 215-234.
%
% Part of: Hypersonics Stability MATLAB Solver - Energy Budgets Module

    %% Unpack grid dimensions
    N_f = length(w_infty);
    Nx  = s.mesh.Nchi;
    Ny  = s.mesh.Neta;
    base_offset = 4 * Nx * Ny;

    %% Extract upstream perturbations for each frequency mode (vectorized)
    rho_infty_all   = zeros(Nx * N_f, 1);
    rho_u_infty_all = zeros(Nx * N_f, 1);
    rho_v_infty_all = zeros(Nx * N_f, 1);
    rho_E_infty_all = zeros(Nx * N_f, 1);

    for i = 0:N_f-1
        row_range = i*Nx + 1 : i*Nx + Nx;
        col_base  = base_offset + 4*i*Nx;
        rho_infty_all(row_range, 1)   = V(col_base + 1         : col_base + Nx,   1);
        rho_u_infty_all(row_range, 1) = V(col_base + Nx + 1    : col_base + 2*Nx, 1);
        rho_v_infty_all(row_range, 1) = V(col_base + 2*Nx + 1  : col_base + 3*Nx, 1);
        rho_E_infty_all(row_range, 1) = V(col_base + 3*Nx + 1  : col_base + 4*Nx, 1);
    end

    %% Sum perturbations across all frequency modes
    indices = (1:Nx)';

    pert_rho   = zeros(Nx, 1);
    pert_rho_u = zeros(Nx, 1);
    pert_rho_v = zeros(Nx, 1);
    pert_rho_E = zeros(Nx, 1);

    for i = 1:N_f
        idx = indices + (i - 1) * Nx;
        pert_rho   = pert_rho   + real(rho_infty_all(idx, 1));
        pert_rho_u = pert_rho_u + real(rho_u_infty_all(idx, 1));
        pert_rho_v = pert_rho_v + real(rho_v_infty_all(idx, 1));
        pert_rho_E = pert_rho_E + real(rho_E_infty_all(idx, 1));
    end

    %% Compute base and perturbed post-shock states via Rankine-Hugoniot
    s = UPDATE_SHOCK_JUMP_PROPERTIES(s, chemistry);
    base_properties = s.shock.properties;

    s.freestream.rho_0_p   = s.freestream.rho_0   + pert_rho;
    s.freestream.rho_u_0_p = s.freestream.rho_u_0 + pert_rho_u;
    s.freestream.rho_v_0_p = s.freestream.rho_v_0 + pert_rho_v;
    s.freestream.rho_E_0_p = s.freestream.rho_E_0 + pert_rho_E;

    s = UPDATE_SHOCK_JUMP_PROPERTIES(s, chemistry);
    pert_properties = s.shock.properties;

    % Post-shock perturbation = perturbed state minus base state
    pert_rho_2   = pert_properties.rho   - base_properties.rho;
    pert_rho_u_2 = pert_properties.rho_u - base_properties.rho_u;
    pert_rho_v_2 = pert_properties.rho_v - base_properties.rho_v;
    pert_rho_E_2 = pert_properties.rho_E - base_properties.rho_E;

    %% Compute energy flux at post-shock interface
    % Locate shocked cell within the Nx-by-Ny grid
    valid_ix = (1:Nx)';
    sc_idy   = s.shock.cell_indices(valid_ix, 1);
    lin_idx_shock     = sub2ind(size(s.mesh.x),     valid_ix, sc_idy - 1);
    lin_idx_shock_Ext = sub2ind(size(s.mesh.x_Ext), valid_ix + 1, sc_idy);

    % Effective shock surface area
    shock_area_effective = s.mesh.bt_area(lin_idx_shock) .* s.mesh.bt_y_normal(lin_idx_shock);

    % Post-shock base-flow quantities
    rho_2        = s.var.rho(lin_idx_shock_Ext);
    rho_u_2      = s.var.rho_u(lin_idx_shock_Ext);
    rho_v_2      = s.var.rho_v(lin_idx_shock_Ext);
    gamma_star_2 = s.var.gamma_star(lin_idx_shock_Ext);
    p_2          = s.var.p(lin_idx_shock_Ext);
    a_2          = s.var.a(lin_idx_shock_Ext);

    %% Evaluate Chu energy norm and components at post-shock
    E_2 = GET_ENERGY_NORM(rho_2, rho_u_2, rho_v_2, p_2, a_2, gamma_star_2, ...
                          pert_rho_2, pert_rho_u_2, pert_rho_v_2, pert_rho_E_2);

    [pert_p, pert_u, pert_v, ~] = GET_CHU_COMPONENTS( ...
        rho_2, rho_u_2, rho_v_2, p_2, a_2, gamma_star_2, ...
        pert_rho_2, pert_rho_u_2, pert_rho_v_2, pert_rho_E_2);

    % Shock-normal velocity perturbation
    u_pert_norm_shock = -pert_u .* sin(s.shock.beta(valid_ix, 1)) + ...
                        pert_v .* cos(s.shock.beta(valid_ix, 1));

    %% Integrate energy fluxes over shock surface
    budgets_shock.adv_acoustic = sum(E_2.acoustic .* abs(rho_v_2 ./ rho_2) .* shock_area_effective, "all");
    budgets_shock.adv_kinetic  = sum(E_2.kinetic  .* abs(rho_v_2 ./ rho_2) .* shock_area_effective, "all");
    budgets_shock.work_kinetic = -sum(pert_p .* u_pert_norm_shock .* shock_area_effective, "all");
    budgets_shock.adv_entropic = sum(E_2.entropic .* abs(rho_v_2 ./ rho_2) .* shock_area_effective, "all");
end
