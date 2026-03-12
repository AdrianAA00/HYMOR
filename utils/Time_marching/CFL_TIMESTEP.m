function s = CFL_TIMESTEP(s)
% CFL_TIMESTEP  Compute an adaptive time step based on the CFL condition.
%
%   s = CFL_TIMESTEP(s)
%
%   Evaluates convective and viscous CFL constraints on the structured
%   grid to determine a stable time step. The final dt is the minimum of:
%     1) Convective CFL limit (left-right and south-north face sweeps),
%     2) Flux-magnitude limit on each conserved variable,
%     3) Viscous diffusion limit (Reynolds and Prandtl-scaled),
%     4) User-specified maximum time step (s.time_integration.max_dt).
%
%   Inputs:
%       s - struct : Solution structure containing conservative
%                           variables (s.var.*), grid metrics (s.mesh.lr_area,
%                           s.mesh.bt_area, s.mesh.lr_x_normal, s.mesh.lr_y_normal,
%                           s.mesh.bt_x_normal, s.mesh.bt_y_normal, s.mesh.volume,
%                           s.mesh.x_Ext, s.mesh.y_Ext), sound speed (s.var.a),
%                           viscosity, Prandtl number, CFL number
%                           (s.time_integration.CFL), and freestream reference
%                           quantities.
%
%   Outputs:
%       s - struct : Same struct with s.time_integration.dt updated to the
%                           new stable time step.
%
%   Notes:
%       - For shock-fitted simulations (s.shock.enabled == true) the flux-
%         magnitude limit is masked by s.shock.flow_cells so that only
%         active flow cells contribute.
%       - The viscous limit uses the minimum centroid distance between
%         south-north and left-right neighbours.
%
% Part of: Hypersonics Stability MATLAB Solver - Time Marching Module

    %% Sound speed estimate
    c = max(abs(s.var.a), [], "all");

    %% Left-Right convective CFL
    u_lr = abs((s.var.rho_u(1:end-1,2:end-1) + s.var.rho_u(2:end,2:end-1)) ./ ...
               (s.var.rho(2:end,2:end-1) + s.var.rho(1:end-1,2:end-1)) .* s.mesh.lr_x_normal);
    v_lr = abs((s.var.rho_v(1:end-1,2:end-1) + s.var.rho_v(2:end,2:end-1)) ./ ...
               (s.var.rho(2:end,2:end-1) + s.var.rho(1:end-1,2:end-1)) .* s.mesh.lr_y_normal);
    U_lr     = u_lr + v_lr;
    U_lr_vol = max(U_lr(1:end-1,:), U_lr(2:end,:));
    A_lr     = max(s.mesh.lr_area(1:end-1,:), s.mesh.lr_area(2:end,:));
    U_lr_A   = (U_lr_vol + c) .* A_lr;
    U_lr_A_eig = max(U_lr_A .* s.shock.flow_cells ./ s.mesh.volume, [], "all");

    %% South-North convective CFL
    u_bt = abs((s.var.rho_u(2:end-1,1:end-1) + s.var.rho_u(2:end-1,2:end)) ./ ...
               (s.var.rho(2:end-1,2:end) + s.var.rho(2:end-1,1:end-1)) .* s.mesh.bt_x_normal);
    v_bt = abs((s.var.rho_v(2:end-1,1:end-1) + s.var.rho_v(2:end-1,2:end)) ./ ...
               (s.var.rho(2:end-1,2:end) + s.var.rho(2:end-1,1:end-1)) .* s.mesh.bt_y_normal);
    U_bt     = u_bt + v_bt;
    U_bt_vol = max(U_bt(:,1:end-1), U_bt(:,2:end));
    A_bt     = max(s.mesh.bt_area(:,1:end-1), s.mesh.bt_area(:,2:end));
    U_bt_A   = (U_bt_vol + c) .* A_bt;
    U_bt_A_eig = max(U_bt_A .* s.shock.flow_cells ./ s.mesh.volume, [], "all");

    %% Convective time step
    dt1 = s.time_integration.CFL / (U_lr_A_eig + U_bt_A_eig);

    %% Viscous time step
    centroids_bt_distance = sqrt((s.mesh.y_Ext(2:end-1,2:end) - s.mesh.y_Ext(2:end-1,1:end-1)).^2 + ...
                                 (s.mesh.x_Ext(2:end-1,2:end) - s.mesh.x_Ext(2:end-1,1:end-1)).^2);
    centroids_lr_distance = sqrt((s.mesh.y_Ext(2:end,2:end-1) - s.mesh.y_Ext(1:end-1,2:end-1)).^2 + ...
                                 (s.mesh.x_Ext(2:end,2:end-1) - s.mesh.x_Ext(1:end-1,2:end-1)).^2);

    Re_flow = s.var.Re_flow(2:end-1,2:end-1);
    Pr_flow = s.var.Pr_flow(2:end-1,2:end-1);
    d = min(centroids_bt_distance(:,2:end), centroids_lr_distance(2:end,:));

    one_Re    = max(abs(1 ./ Re_flow) ./ d.^2 .* s.shock.flow_cells, [], "all") * 4;
    one_Re_Pr = max(abs(1 ./ (Re_flow .* Pr_flow) ./ d.^2 .* s.shock.flow_cells), [], "all") * 4;
    dt1_visc  = min(s.time_integration.CFL / one_Re, s.time_integration.CFL / one_Re_Pr);

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
    dt3 = min(dt1, dt2);
    dt4 = min(dt3, dt1_visc);
    s.time_integration.dt = min(dt4, s.time_integration.max_dt);

end
