function s = CFL_TIMESTEP(s)
% CFL_TIMESTEP  Compute an adaptive time step based on the CFL condition.
%
%   s = CFL_TIMESTEP(s)
%
%   Evaluates the combined convective + viscous spectral radius per cell on
%   the structured grid and returns the most restrictive stable time step.
%   The final dt is the minimum of:
%     1) Combined convective + viscous CFL limit (per-cell Blazek form),
%     2) Flux-magnitude limit on each conserved variable,
%     3) User-specified maximum time step (s.time_integration.max_dt).
%
%   Inputs:
%       s - struct : Solution structure containing conservative
%                           variables (s.var.*), grid metrics (s.mesh.lr_area,
%                           s.mesh.bt_area, s.mesh.lr_x_normal, s.mesh.lr_y_normal,
%                           s.mesh.bt_x_normal, s.mesh.bt_y_normal, s.mesh.volume),
%                           sound speed (s.var.a), local gamma, viscosity,
%                           Prandtl number, CFL number (s.time_integration.CFL),
%                           and freestream reference quantities.
%
%   Outputs:
%       s - struct : Same struct with s.time_integration.dt updated to the
%                           new stable time step.
%
%   Notes:
%       - Convective spectral radius uses the proper face-normal velocity
%         |u*n_x + v*n_y| (not the conservative bound |u||n_x|+|v||n_y|).
%       - Viscous spectral radius follows Blazek's compressible-NS form:
%         lambda_v = max(4/3, gamma/Pr) / (rho * Re_flow) * (A_xi^2 + A_eta^2) / V,
%         which includes the 1/rho factor (kinematic viscosity), the 4/3
%         deviatoric-stress prefactor, and the gamma/Pr thermal-diffusion
%         eigenvalue (whichever is larger).
%       - Convective and viscous contributions are summed per cell before
%         taking the global minimum dt (the proper combined explicit-CFL
%         constraint), rather than min(dt_conv, dt_visc).
%       - For shock-fitted simulations (s.shock.enabled == true) the limits
%         are masked by s.shock.flow_cells so that only active flow cells
%         contribute.
%
% Part of: Hypersonics Stability MATLAB Solver - Time Marching Module

    %% Sound speed estimate
    c = max(abs(s.var.a), [], "all");

    %% Convective spectral radius (per cell)
    % Left-right face-centered velocities (mass-weighted average of neighbours)
    rho_sum_lr = s.var.rho(1:end-1,2:end-1) + s.var.rho(2:end,2:end-1);
    u_face_lr  = (s.var.rho_u(1:end-1,2:end-1) + s.var.rho_u(2:end,2:end-1)) ./ rho_sum_lr;
    v_face_lr  = (s.var.rho_v(1:end-1,2:end-1) + s.var.rho_v(2:end,2:end-1)) ./ rho_sum_lr;
    Vn_lr       = abs(u_face_lr .* s.mesh.lr_x_normal + v_face_lr .* s.mesh.lr_y_normal);
    lam_face_lr = (Vn_lr + c) .* s.mesh.lr_area;

    % Bottom-top face-centered velocities
    rho_sum_bt = s.var.rho(2:end-1,1:end-1) + s.var.rho(2:end-1,2:end);
    u_face_bt  = (s.var.rho_u(2:end-1,1:end-1) + s.var.rho_u(2:end-1,2:end)) ./ rho_sum_bt;
    v_face_bt  = (s.var.rho_v(2:end-1,1:end-1) + s.var.rho_v(2:end-1,2:end)) ./ rho_sum_bt;
    Vn_bt       = abs(u_face_bt .* s.mesh.bt_x_normal + v_face_bt .* s.mesh.bt_y_normal);
    lam_face_bt = (Vn_bt + c) .* s.mesh.bt_area;

    % Collapse opposing faces to per-cell xi/eta spectral radii (max of the two)
    lam_xi  = max(lam_face_lr(1:end-1,:), lam_face_lr(2:end,:));
    lam_eta = max(lam_face_bt(:,1:end-1), lam_face_bt(:,2:end));
    lam_c   = lam_xi + lam_eta;

    %% Viscous spectral radius (per cell, Blazek-style)
    rho_int   = s.var.rho(2:end-1,2:end-1);
    gamma_int = s.var.gamma_star(2:end-1,2:end-1);
    Re_flow   = s.var.Re_flow(2:end-1,2:end-1);
    Pr_flow   = s.var.Pr_flow(2:end-1,2:end-1);

    A_xi_cell  = 0.5 * (s.mesh.lr_area(1:end-1,:) + s.mesh.lr_area(2:end,:));
    A_eta_cell = 0.5 * (s.mesh.bt_area(:,1:end-1) + s.mesh.bt_area(:,2:end));

    nu_eff = max(4/3, gamma_int ./ Pr_flow) ./ (rho_int .* Re_flow);
    lam_v  = 2 * nu_eff .* (A_xi_cell.^2 + A_eta_cell.^2) ./ s.mesh.volume;

    %% Combined convective + viscous time step (per-cell sum, global min)
    lam_over_V = (lam_c + lam_v) .* s.shock.flow_cells ./ s.mesh.volume;
    dt1 = s.time_integration.CFL / max(lam_over_V, [], "all");

    %% Flux-magnitude time step
    if s.shock.enabled == true
        max_val_u   = max(((abs(s.flux.rho_v) + abs(s.flux.rho_u)) ./ ...
                           s.var.rho(2:end-1,2:end-1) .* s.shock.flow_cells), [], "all");
        max_val_rho = max(abs(s.flux.rho ./ s.var.rho(2:end-1,2:end-1) .* s.shock.flow_cells), [], "all");
        max_val_E   = max(abs(s.flux.rho_E ./ s.var.rho_E(2:end-1,2:end-1) .* s.shock.flow_cells), [], "all");
    else
        max_val_u   = max(((abs(s.flux.rho_v) + abs(s.flux.rho_u)) ./ ...
                           s.var.rho(2:end-1,2:end-1)), [], "all");
        max_val_rho = max(abs(s.flux.rho ./ s.var.rho(2:end-1,2:end-1)), [], "all");
        max_val_E   = max(abs(s.flux.rho_E ./ s.var.rho_E(2:end-1,2:end-1)), [], "all");
    end
    dt2 = s.time_integration.CFL / max([max_val_u, max_val_rho, max_val_E], [], "all");

    %% Final time step: minimum of all constraints
    s.time_integration.dt = min([dt1, dt2, s.time_integration.max_dt]);

end
