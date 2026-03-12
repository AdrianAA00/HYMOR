function budgets_entropy = COMPUTE_BUDGETS_ENTROPY(V, s, chemistry, output_flow)
% COMPUTE_BUDGETS_ENTROPY  Compute entropy energy budget terms for a perturbation field.
%
%   budgets_entropy = COMPUTE_BUDGETS_ENTROPY(V, s, chemistry, output_flow)
%
%   Computes the individual terms in the entropy component of the Chu energy
%   budget equation. The entropy budget decomposes the time rate of change of
%   entropy-related perturbation energy into advection, production (momentum,
%   mass, and temperature), heat-flux transport, thermal dissipation, and
%   viscous source contributions. Both 2-D Cartesian and 3-D axisymmetric
%   formulations are supported.
%
%   Inputs:
%       V           - (double column vector) Global perturbation state vector
%                     containing [rho'; rho_u'; rho_v'; rho_E'] stacked over
%                     the Nx-by-Ny grid (length >= 4*Nx*Ny).
%       s    - (struct) Solution structure with base-flow fields
%                     (rho, rho_u, rho_v, rho_E, p, T, mu, k, a, gamma_star,
%                     cv_star), grid information (Nx, Ny, x, volume,
%                     flow_cells), freestream scaling factors, and the
%                     dimension flag ("3D-axisymmetric" or "2D").
%       chemistry   - (struct) Chemistry model providing eval_s(rho, e) for
%                     entropy evaluation when chemistry_state is active.
%       output_flow - (logical) If true, the full 2-D field of every budget
%                     term is stored in the output structure in addition to
%                     the volume-integrated scalars.
%
%   Outputs:
%       budgets_entropy - (struct) Structure containing volume-integrated
%                         budget terms (suffixed _sum) and, optionally, the
%                         spatially resolved 2-D fields:
%           .A_adv_sum       - Integrated advection of entropy energy.
%           .P_mom_sum       - Integrated momentum production.
%           .P_mass_sum      - Integrated mass (density) production.
%           .P_T_sum         - Integrated temperature production.
%           .Transport_sum   - Integrated heat-flux transport.
%           .Dissipation_sum - Integrated thermal dissipation.
%           .Source_sum      - Integrated viscous source (Phi').
%           .A_adv, .P_mom, .P_mass, .P_T, .Transport, .Dissipation,
%           .Source          - (Nx-by-Ny, only when output_flow is true)
%
%   Notes:
%       - The entropy perturbation is defined through the Chu decomposition
%         (see GET_CHU_COMPONENTS) and scaled by the gas constant R_0.
%       - The scaling coefficient C_scale = (gamma*-1)*p_0 / (gamma* * R_0^2)
%         ensures consistency with the Chu energy norm.
%       - Viscosity and thermal conductivity are non-dimensionalised using
%         freestream reference quantities.
%       - Requires helper functions: EXTEND_TO_GHOST_POINTS, DERIVATIVE_EXT,
%         GET_CHU_COMPONENTS.
%
%   References:
%       Chu, B.-T. (1965), "On the energy transfer to small disturbances in
%       fluid flow (Part I)", Acta Mechanica, 1(3), 215-234.
%
% Part of: Hypersonics Stability MATLAB Solver - Energy Budgets Module

    N_v = s.mesh.Nchi * s.mesh.Neta;

    %% Extract perturbations from state vector
    pert_rho   = real(reshape(V(0*N_v+1 : 1*N_v, 1), s.mesh.Nchi, s.mesh.Neta));
    pert_rho_u = real(reshape(V(1*N_v+1 : 2*N_v, 1), s.mesh.Nchi, s.mesh.Neta));
    pert_rho_v = real(reshape(V(2*N_v+1 : 3*N_v, 1), s.mesh.Nchi, s.mesh.Neta));
    pert_rho_E = real(reshape(V(3*N_v+1 : 4*N_v, 1), s.mesh.Nchi, s.mesh.Neta));

    %% Extract base-flow variables
    rho_0   = s.var.rho(2:end-1, 2:end-1);
    rho_u_0 = s.var.rho_u(2:end-1, 2:end-1);
    rho_v_0 = s.var.rho_v(2:end-1, 2:end-1);
    rho_E_0 = s.var.rho_E(2:end-1, 2:end-1);
    u_0 = rho_u_0 ./ rho_0;
    v_0 = s.var.rho_v(2:end-1, 2:end-1) ./ rho_0;
    e_0 = rho_E_0 ./ rho_0 - (u_0.^2 + v_0.^2) / 2;
    T_0 = s.var.T(2:end-1, 2:end-1);
    p_0 = s.var.p(2:end-1, 2:end-1);
    a_0 = s.var.a(2:end-1, 2:end-1);

    u_0_Ext = EXTEND_TO_GHOST_POINTS(s.var.rho_u(2:end-1, 2:end-1) ./ s.var.rho(2:end-1, 2:end-1), s);
    v_0_Ext = EXTEND_TO_GHOST_POINTS(s.var.rho_v(2:end-1, 2:end-1) ./ s.var.rho(2:end-1, 2:end-1), s);

    %% Non-dimensional scaling and material properties
    mu_0    = s.var.mu_star(2:end-1, 2:end-1) / s.freestream.Re;
    kappa_0 = s.var.k_star(2:end-1, 2:end-1) * s.freestream.gamma_star / s.freestream.Re / s.freestream.Pr;

    gamma_star_0 = s.var.gamma_star(2:end-1, 2:end-1);
    cv_star_0    = s.var.cv_star(2:end-1, 2:end-1);
    R_0 = (gamma_star_0 - 1) .* cv_star_0;

    %% Compute base-flow entropy
    if s.chemistry.is_chemistry_enabled
        S_0 = chemistry.eval_s(s.freestream.rho_factor * rho_0, ...
                               s.freestream.energy_factor * e_0) ./ s.freestream.cv;
    else
        S_0 = cv_star_0 .* log(T_0) + R_0 .* log(p_0);
    end
    S_0_ext = EXTEND_TO_GHOST_POINTS(S_0, s);

    volume     = s.mesh.volume;
    flow_cells = s.shock.flow_cells;

    %% Compute base-flow derivatives
    [du0_dx, du0_dy] = DERIVATIVE_EXT(u_0_Ext, s);
    [dv0_dx, dv0_dy] = DERIVATIVE_EXT(v_0_Ext, s);
    [dS0_dx, dS0_dy] = DERIVATIVE_EXT(S_0_ext, s);

    %% Transform perturbations to Chu variables
    [pert_p, pert_u, pert_v, pert_S_R] = GET_CHU_COMPONENTS( ...
        rho_0, rho_u_0, rho_v_0, p_0, a_0, gamma_star_0, ...
        pert_rho, pert_rho_u, pert_rho_v, pert_rho_E);
    pert_S = pert_S_R .* R_0;

    pert_u_ext = EXTEND_TO_GHOST_POINTS(pert_u, s);
    pert_v_ext = EXTEND_TO_GHOST_POINTS(pert_v, s);
    pert_p_ext = EXTEND_TO_GHOST_POINTS(pert_p, s);

    [dpert_u_dx, dpert_u_dy] = DERIVATIVE_EXT(pert_u_ext, s);
    [dpert_v_dx, dpert_v_dy] = DERIVATIVE_EXT(pert_v_ext, s);

    x = s.mesh.x;  % Radial coordinate at cell centers

    %% Entropy scaling coefficients
    % C_scale = (gamma*_0 - 1) * p_0 / (gamma*_0 * R_0^2)
    C_scale      = ((gamma_star_0 - 1) .* p_0) ./ (gamma_star_0 .* R_0.^2);
    C_scale_therm = C_scale ./ (rho_0 .* T_0);

    %% Entropy and temperature perturbation derivatives
    pert_S_ext = EXTEND_TO_GHOST_POINTS(pert_S, s);
    [dpert_S_dx, dpert_S_dy] = DERIVATIVE_EXT(pert_S_ext, s);

    % Temperature perturbation: T' = T_0 * (p'/p_0 - rho'/rho_0)
    pert_T = T_0 .* (pert_p ./ p_0 - pert_rho ./ rho_0);
    pert_T_ext = EXTEND_TO_GHOST_POINTS(pert_T, s);
    [dpert_T_dx, dpert_T_dy] = DERIVATIVE_EXT(pert_T_ext, s);

    %% Compute budget terms
    if s.PDE_dimension == "3D-axisymmetric"
        % ========== ADVECTION TERM ==========
        % A^s = -C_scale * u_{j,0} * d/dx_j(s'^2/2)
        A_adv = -C_scale .* u_0 .* pert_S .* dpert_S_dx ...
                - C_scale .* v_0 .* pert_S .* dpert_S_dy;

        % ========== PRODUCTION TERMS ==========
        % Momentum production: -C_scale * u'_j * s' * dS_0/dx_j
        P_mom = -C_scale .* pert_u .* pert_S .* dS0_dx ...
                - C_scale .* pert_v .* pert_S .* dS0_dy;

        % Mass production: -C_scale * (rho'/rho_0) * u_{j,0} * s' * dS_0/dx_j
        P_mass = -C_scale .* (pert_rho ./ rho_0) .* u_0 .* pert_S .* dS0_dx ...
                 - C_scale .* (pert_rho ./ rho_0) .* v_0 .* pert_S .* dS0_dy;

        % Temperature production: -C_scale * (T'/T_0) * u_{j,0} * s' * dS_0/dx_j
        P_T = -C_scale .* (pert_T ./ T_0) .* u_0 .* pert_S .* dS0_dx ...
              - C_scale .* (pert_T ./ T_0) .* v_0 .* pert_S .* dS0_dy;

        % ========== HEAT FLUX PERTURBATION ==========
        % q'_j = -kappa_0 * dT'/dx_j
        q_x_prime = -kappa_0 .* dpert_T_dx;
        q_y_prime = -kappa_0 .* dpert_T_dy;

        % ========== VISCOUS DISSIPATION FUNCTION ==========
        % Base-flow stress tensor (axisymmetric)
        div_base = du0_dx + u_0 ./ x + dv0_dy;
        tau_xx_0          = mu_0 .* (2 * du0_dx    - 2/3 * div_base);
        tau_yy_0          = mu_0 .* (2 * dv0_dy    - 2/3 * div_base);
        tau_theta_theta_0 = mu_0 .* (2 * u_0 ./ x  - 2/3 * div_base);
        tau_xy_0          = mu_0 .* (du0_dy + dv0_dx);

        % Perturbation stress tensor (axisymmetric)
        div_pert = dpert_u_dx + pert_u ./ x + dpert_v_dy;
        pert_tau_xx          = mu_0 .* (2 * dpert_u_dx    - 2/3 * div_pert);
        pert_tau_yy          = mu_0 .* (2 * dpert_v_dy    - 2/3 * div_pert);
        pert_tau_theta_theta = mu_0 .* (2 * pert_u ./ x   - 2/3 * div_pert);
        pert_tau_xy          = mu_0 .* (dpert_u_dy + dpert_v_dx);

        % Phi' = tau_{ij,0} * du'_i/dx_j + tau'_{ij} * du_{i,0}/dx_j
        Phi_prime = tau_xx_0          .* dpert_u_dx           + ...
                    tau_theta_theta_0 .* (pert_u ./ x)        + ...
                    tau_yy_0          .* dpert_v_dy           + ...
                    tau_xy_0          .* (dpert_u_dy + dpert_v_dx) + ...
                    pert_tau_xx          .* du0_dx             + ...
                    pert_tau_theta_theta .* (u_0 ./ x)        + ...
                    pert_tau_yy          .* dv0_dy             + ...
                    pert_tau_xy          .* (du0_dy + dv0_dx);

        % ========== TRANSPORT TERM ==========
        % T^s = -C_scale_therm * div(s' * q'_j) with axisymmetric geometry
        Flux_entropy_x = x .* pert_S .* q_x_prime;
        Flux_entropy_y = pert_S .* q_y_prime;

        Flux_entropy_x_Ext = EXTEND_TO_GHOST_POINTS(Flux_entropy_x, s);
        Flux_entropy_y_Ext = EXTEND_TO_GHOST_POINTS(Flux_entropy_y, s);
        [dFlux_entropy_x_dx, ~] = DERIVATIVE_EXT(Flux_entropy_x_Ext, s);
        [~, dFlux_entropy_y_dy] = DERIVATIVE_EXT(Flux_entropy_y_Ext, s);

        Transport = -C_scale_therm .* ((1 ./ x) .* dFlux_entropy_x_dx + dFlux_entropy_y_dy);

        % ========== DISSIPATION TERM ==========
        % D^s = C_scale_therm * (ds'/dx_j * q'_j)
        Dissipation = C_scale_therm .* (dpert_S_dx .* q_x_prime + ...
                                        dpert_S_dy .* q_y_prime);

        % ========== SOURCE TERM ==========
        % S^s = C_scale_therm * (s' * Phi')
        Source = C_scale_therm .* (pert_S .* Phi_prime);

    else  % Cartesian 2D
        % ========== ADVECTION TERM ==========
        A_adv = -C_scale .* u_0 .* pert_S .* dpert_S_dx ...
                - C_scale .* v_0 .* pert_S .* dpert_S_dy;

        % ========== PRODUCTION TERMS ==========
        P_mom = -C_scale .* pert_u .* pert_S .* dS0_dx ...
                - C_scale .* pert_v .* pert_S .* dS0_dy;

        P_mass = -C_scale .* (pert_rho ./ rho_0) .* u_0 .* pert_S .* dS0_dx ...
                 - C_scale .* (pert_rho ./ rho_0) .* v_0 .* pert_S .* dS0_dy;

        P_T = -C_scale .* (pert_T ./ T_0) .* u_0 .* pert_S .* dS0_dx ...
              - C_scale .* (pert_T ./ T_0) .* v_0 .* pert_S .* dS0_dy;

        % ========== HEAT FLUX PERTURBATION ==========
        q_x_prime = -kappa_0 .* dpert_T_dx;
        q_y_prime = -kappa_0 .* dpert_T_dy;

        % ========== VISCOUS DISSIPATION FUNCTION ==========
        % Base-flow stress tensor (Cartesian)
        div_base = du0_dx + dv0_dy;
        tau_xx_0 = mu_0 .* (2 * du0_dx - 2/3 * div_base);
        tau_yy_0 = mu_0 .* (2 * dv0_dy - 2/3 * div_base);
        tau_xy_0 = mu_0 .* (du0_dy + dv0_dx);

        % Perturbation stress tensor (Cartesian)
        div_pert = dpert_u_dx + dpert_v_dy;
        pert_tau_xx = mu_0 .* (2 * dpert_u_dx - 2/3 * div_pert);
        pert_tau_yy = mu_0 .* (2 * dpert_v_dy - 2/3 * div_pert);
        pert_tau_xy = mu_0 .* (dpert_u_dy + dpert_v_dx);

        Phi_prime = tau_xx_0 .* dpert_u_dx                    + ...
                    tau_yy_0 .* dpert_v_dy                    + ...
                    tau_xy_0 .* (dpert_u_dy + dpert_v_dx)     + ...
                    pert_tau_xx .* du0_dx                     + ...
                    pert_tau_yy .* dv0_dy                     + ...
                    pert_tau_xy .* (du0_dy + dv0_dx);

        % ========== TRANSPORT TERM ==========
        Flux_entropy_x = pert_S .* q_x_prime;
        Flux_entropy_y = pert_S .* q_y_prime;

        Flux_entropy_x_Ext = EXTEND_TO_GHOST_POINTS(Flux_entropy_x, s);
        Flux_entropy_y_Ext = EXTEND_TO_GHOST_POINTS(Flux_entropy_y, s);
        [dFlux_entropy_x_dx, ~] = DERIVATIVE_EXT(Flux_entropy_x_Ext, s);
        [~, dFlux_entropy_y_dy] = DERIVATIVE_EXT(Flux_entropy_y_Ext, s);

        Transport = -C_scale_therm .* (dFlux_entropy_x_dx + dFlux_entropy_y_dy);

        % ========== DISSIPATION TERM ==========
        Dissipation = C_scale_therm .* (dpert_S_dx .* q_x_prime + ...
                                        dpert_S_dy .* q_y_prime);

        % ========== SOURCE TERM ==========
        Source = C_scale_therm .* (pert_S .* Phi_prime);
    end

    %% Compute volume-integrated budgets
    budgets_entropy.A_adv_sum       = sum(A_adv       .* volume .* flow_cells, "all");
    budgets_entropy.P_mom_sum       = sum(P_mom       .* volume .* flow_cells, "all");
    budgets_entropy.P_mass_sum      = sum(P_mass      .* volume .* flow_cells, "all");
    budgets_entropy.P_T_sum         = sum(P_T         .* volume .* flow_cells, "all");
    budgets_entropy.Transport_sum   = sum(Transport   .* volume .* flow_cells, "all");
    budgets_entropy.Dissipation_sum = sum(Dissipation .* volume .* flow_cells, "all");
    budgets_entropy.Source_sum      = sum(Source      .* volume .* flow_cells, "all");

    %% Optionally store full spatial fields
    if output_flow
        budgets_entropy.A_adv       = A_adv;
        budgets_entropy.P_mom       = P_mom;
        budgets_entropy.P_mass      = P_mass;
        budgets_entropy.P_T         = P_T;
        budgets_entropy.Transport   = Transport;
        budgets_entropy.Dissipation = Dissipation;
        budgets_entropy.Source      = Source;
    end
end
