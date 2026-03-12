function s = VISCOUS_FLUX(s)
% VISCOUS_FLUX  Compute 2D viscous flux contributions for momentum and energy equations.
%
%   s = VISCOUS_FLUX(s)
%
%   Evaluates viscous stress tensor components and heat flux on a general
%   non-orthogonal quadrilateral mesh using non-orthogonal correction.
%   Velocity and temperature gradients are decomposed into normal and
%   tangential components at each cell face. The viscous stress work and
%   Fourier heat conduction are added to the energy flux.
%
%   Inputs:
%       s - Solution struct containing:
%                  .var.rho, .var.rho_u, .var.rho_v, .var.rho_E - Conserved variables
%                  .mesh.x_Ext, .mesh.y_Ext       - Extended cell centroid coordinates
%                  .mesh.x_corner, .mesh.y_corner - Cell corner coordinates
%                  .mesh.bt_x_normal, .mesh.bt_y_normal - Bottom-top face unit normals
%                  .mesh.lr_x_normal, .mesh.lr_y_normal - Left-right face unit normals
%                  .mesh.bt_area, .mesh.lr_area   - Face areas
%                  .bt_length, .lr_length         - Face lengths
%                  .mesh.volume                   - Cell volumes
%                  .var.mu_star                   - Dynamic viscosity field
%                  .var.k_star                    - Thermal conductivity field
%                  .var.T                         - Temperature field
%                  .freestream.*                  - Non-dimensionalization parameters
%
%   Outputs:
%       s - Updated s struct with viscous contributions added
%                  to .flux.rho_u, .flux.rho_v, and .flux.rho_E
%
%   Notes:
%       - Uses non-orthogonal correction for gradient reconstruction.
%       - The viscous tensor uses the deviatoric form with the 2/3
%         divergence term.
%       - Heat flux is computed via Fourier's law with corrected normal
%         temperature gradients.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% Pre-calculate velocity fields at cell centers
    u = s.var.rho_u ./ s.var.rho;
    v = s.var.rho_v ./ s.var.rho;

    %% Velocity values at cell corners (4-point average)
    u_cell_corners = s.var.rho_u ./ s.var.rho;
    v_cell_corners = s.var.rho_v ./ s.var.rho;

    u_c = (u_cell_corners(1:end-1, 1:end-1) + u_cell_corners(2:end, 1:end-1) + ...
           u_cell_corners(1:end-1, 2:end)   + u_cell_corners(2:end, 2:end)) / 4;
    v_c = (v_cell_corners(1:end-1, 1:end-1) + v_cell_corners(2:end, 1:end-1) + ...
           v_cell_corners(1:end-1, 2:end)   + v_cell_corners(2:end, 2:end)) / 4;

    %% Distances between cell centroids
    centroids_bt_distance = sqrt((s.mesh.y_Ext(2:end-1, 2:end) - s.mesh.y_Ext(2:end-1, 1:end-1)).^2 + ...
                                 (s.mesh.x_Ext(2:end-1, 2:end) - s.mesh.x_Ext(2:end-1, 1:end-1)).^2);
    centroids_lr_distance = sqrt((s.mesh.y_Ext(2:end, 2:end-1) - s.mesh.y_Ext(1:end-1, 2:end-1)).^2 + ...
                                 (s.mesh.x_Ext(2:end, 2:end-1) - s.mesh.x_Ext(1:end-1, 2:end-1)).^2);

    %% Edge lengths of cell faces
    bt_distance = sqrt((s.mesh.x_corner(2:end, :) - s.mesh.x_corner(1:end-1, :)).^2 + ...
                       (s.mesh.y_corner(2:end, :) - s.mesh.y_corner(1:end-1, :)).^2);
    lr_distance = sqrt((s.mesh.x_corner(:, 2:end) - s.mesh.x_corner(:, 1:end-1)).^2 + ...
                       (s.mesh.y_corner(:, 2:end) - s.mesh.y_corner(:, 1:end-1)).^2);

    %% Unit vector components for non-orthogonal correction
    c1_x = -(s.mesh.x_Ext(2:end, 2:end-1) - s.mesh.x_Ext(1:end-1, 2:end-1)) ./ centroids_lr_distance;
    c1_y = -(s.mesh.y_Ext(2:end, 2:end-1) - s.mesh.y_Ext(1:end-1, 2:end-1)) ./ centroids_lr_distance;
    c2_x =  (s.mesh.x_Ext(2:end-1, 2:end) - s.mesh.x_Ext(2:end-1, 1:end-1)) ./ centroids_bt_distance;
    c2_y =  (s.mesh.y_Ext(2:end-1, 2:end) - s.mesh.y_Ext(2:end-1, 1:end-1)) ./ centroids_bt_distance;

    %% Projection matrices for reference system transformation
    Ct2_Fn2 = zeros(size(s.mesh.bt_x_normal));
    Ct2_Ft2 = ones(size(s.mesh.bt_x_normal));
    Cn2_Fn2 = c2_x .* s.mesh.bt_x_normal + c2_y .* s.mesh.bt_y_normal;
    Cn2_Ft2 = c2_x .* s.mesh.bt_y_normal - c2_y .* s.mesh.bt_x_normal;

    Cn1_Fn1 = -c1_x .* s.mesh.lr_x_normal - c1_y .* s.mesh.lr_y_normal;
    Cn1_Ft1 =  c1_x .* s.mesh.lr_y_normal - c1_y .* s.mesh.lr_x_normal;
    Ct1_Fn1 = zeros(size(s.mesh.lr_x_normal));
    Ct1_Ft1 = ones(size(s.mesh.lr_y_normal));

    %% Velocity derivatives along cell-center lines
    du_dc2 =  (u(2:end-1, 2:end) - u(2:end-1, 1:end-1)) ./ centroids_bt_distance;
    dv_dc2 =  (v(2:end-1, 2:end) - v(2:end-1, 1:end-1)) ./ centroids_bt_distance;
    du_dc1 = -(u(2:end, 2:end-1) - u(1:end-1, 2:end-1)) ./ centroids_lr_distance;
    dv_dc1 = -(v(2:end, 2:end-1) - v(1:end-1, 2:end-1)) ./ centroids_lr_distance;

    %% Tangential velocity derivatives at faces
    du_dt2 = -(u_c(2:end, :) - u_c(1:end-1, :)) ./ bt_distance;
    dv_dt2 = -(v_c(2:end, :) - v_c(1:end-1, :)) ./ bt_distance;
    du_dt1 =  (u_c(:, 2:end) - u_c(:, 1:end-1)) ./ lr_distance;
    dv_dt1 =  (v_c(:, 2:end) - v_c(:, 1:end-1)) ./ lr_distance;

    %% Normal derivatives at face 2 (bottom-top faces)
    divide = Ct2_Fn2 .* Cn2_Ft2 - Cn2_Fn2 .* Ct2_Ft2;
    du_dn2 = (Cn2_Ft2 .* du_dt2 - du_dc2 .* Ct2_Ft2) ./ divide;
    dv_dn2 = (Cn2_Ft2 .* dv_dt2 - dv_dc2 .* Ct2_Ft2) ./ divide;

    %% Transform to local coordinate system (n2, t2) at face 2
    dut_dt2 = du_dt2 .* s.mesh.bt_y_normal - dv_dt2 .* s.mesh.bt_x_normal;
    dun_dt2 = du_dt2 .* s.mesh.bt_x_normal + dv_dt2 .* s.mesh.bt_y_normal;
    dut_dn2 = du_dn2 .* s.mesh.bt_y_normal - dv_dn2 .* s.mesh.bt_x_normal;
    dun_dn2 = du_dn2 .* s.mesh.bt_x_normal + dv_dn2 .* s.mesh.bt_y_normal;

    %% Normal derivatives at face 1 (left-right faces)
    divide = Cn1_Fn1 .* Ct1_Ft1 - Ct1_Fn1 .* Cn1_Ft1;
    du_dn1 = (Ct1_Ft1 .* du_dc1 - du_dt1 .* Cn1_Ft1) ./ divide;
    dv_dn1 = (Ct1_Ft1 .* dv_dc1 - dv_dt1 .* Cn1_Ft1) ./ divide;

    %% Transform to local coordinate system (n1, t1) at face 1
    dun_dn1 = -du_dn1 .* s.mesh.lr_x_normal - dv_dn1 .* s.mesh.lr_y_normal;
    dut_dn1 =  du_dn1 .* s.mesh.lr_y_normal - dv_dn1 .* s.mesh.lr_x_normal;
    dun_dt1 = -du_dt1 .* s.mesh.lr_x_normal - dv_dt1 .* s.mesh.lr_y_normal;
    dut_dt1 =  du_dt1 .* s.mesh.lr_y_normal - dv_dt1 .* s.mesh.lr_x_normal;

    %% Face-averaged viscosity (non-dimensionalized)
    mu_1 = (s.var.mu_star(2:end, 2:end-1) + s.var.mu_star(1:end-1, 2:end-1)) / 2;
    mu_2 = (s.var.mu_star(2:end-1, 2:end) + s.var.mu_star(2:end-1, 1:end-1)) / 2;

    %% Viscous stress tensor (deviatoric form)
    Re_infty = s.freestream.Re;
    div_1 = dut_dt1 + dun_dn1;
    tau_11 = mu_1/Re_infty .* (2 * dun_dn1 - 2/3 * div_1);
    tau_12 = mu_1/Re_infty .* (dun_dt1 + dut_dn1);
    div_2 = dut_dt2 + dun_dn2;
    tau_22 = mu_2/Re_infty .* (2 * dun_dn2 - 2/3 * div_2);
    tau_21 = mu_2/Re_infty .* (dut_dn2 + dun_dt2);

    %% Viscous stress in global (x-y) coordinates
    tau_1x = -(tau_11 .* s.mesh.lr_x_normal - tau_12 .* s.mesh.lr_y_normal);
    tau_1y = -(tau_11 .* s.mesh.lr_y_normal + tau_12 .* s.mesh.lr_x_normal);
    tau_2x =  (tau_22 .* s.mesh.bt_x_normal + tau_21 .* s.mesh.bt_y_normal);
    tau_2y =  (tau_22 .* s.mesh.bt_y_normal - tau_21 .* s.mesh.bt_x_normal);

    %% Viscous momentum flux (x-direction)
    rho_u_temp = (tau_1x(1:end-1, :) .* s.mesh.lr_area(1:end-1, :) - tau_1x(2:end, :) .* s.mesh.lr_area(2:end, :)) + ...
                 (tau_2x(:, 2:end)   .* s.mesh.bt_area(:, 2:end)   - tau_2x(:, 1:end-1) .* s.mesh.bt_area(:, 1:end-1));
    s.flux.rho_u = s.flux.rho_u + rho_u_temp ./ s.mesh.volume;

    %% Viscous momentum flux (y-direction)
    rho_v_temp = (tau_1y(1:end-1, :) .* s.mesh.lr_area(1:end-1, :) - tau_1y(2:end, :) .* s.mesh.lr_area(2:end, :)) + ...
                 (tau_2y(:, 2:end)   .* s.mesh.bt_area(:, 2:end)   - tau_2y(:, 1:end-1) .* s.mesh.bt_area(:, 1:end-1));
    s.flux.rho_v = s.flux.rho_v + rho_v_temp ./ s.mesh.volume;

    %% Viscous work (stress power) contribution to energy flux
    inv_rho_1 = 1 ./ (s.var.rho(1:end-1, 2:end-1) + s.var.rho(2:end, 2:end-1));
    inv_rho_2 = 1 ./ (s.var.rho(2:end-1, 1:end-1) + s.var.rho(2:end-1, 2:end));

    u_1 = (s.var.rho_u(1:end-1, 2:end-1) + s.var.rho_u(2:end, 2:end-1)) .* inv_rho_1;
    u_2 = (s.var.rho_u(2:end-1, 1:end-1) + s.var.rho_u(2:end-1, 2:end)) .* inv_rho_2;
    v_1 = (s.var.rho_v(1:end-1, 2:end-1) + s.var.rho_v(2:end, 2:end-1)) .* inv_rho_1;
    v_2 = (s.var.rho_v(2:end-1, 1:end-1) + s.var.rho_v(2:end-1, 2:end)) .* inv_rho_2;

    w1 = u_1 .* tau_1x + v_1 .* tau_1y;
    w2 = u_2 .* tau_2x + v_2 .* tau_2y;

    rho_E_temp = (w1(1:end-1, :) .* s.mesh.lr_area(1:end-1, :) - w1(2:end, :) .* s.mesh.lr_area(2:end, :)) + ...
                 (w2(:, 2:end)   .* s.mesh.bt_area(:, 2:end)   - w2(:, 1:end-1) .* s.mesh.bt_area(:, 1:end-1));
    s.flux.rho_E = s.flux.rho_E + rho_E_temp ./ s.mesh.volume;

    %% Heat flux with non-orthogonal correction (Fourier's law)
    k_bt = (s.var.k_star(2:end-1, 2:end) + s.var.k_star(2:end-1, 1:end-1)) / 2 * s.freestream.gamma_star / s.freestream.Re / s.freestream.Pr;
    k_lr = (s.var.k_star(2:end, 2:end-1) + s.var.k_star(1:end-1, 2:end-1)) / 2 * s.freestream.gamma_star / s.freestream.Re / s.freestream.Pr;

    %% Temperature at cell corners (4-point average)
    T_cell_corners = s.var.T;

    T_c = (T_cell_corners(1:end-1, 1:end-1) + T_cell_corners(2:end, 1:end-1) + ...
           T_cell_corners(1:end-1, 2:end)   + T_cell_corners(2:end, 2:end)) / 4;

    %% Temperature derivatives along cell-center lines
    dT_dc2 =  (s.var.T(2:end-1, 2:end) - s.var.T(2:end-1, 1:end-1)) ./ centroids_bt_distance;
    dT_dc1 = -(s.var.T(2:end, 2:end-1) - s.var.T(1:end-1, 2:end-1)) ./ centroids_lr_distance;

    %% Tangential temperature derivatives at faces
    dT_dt2 = -(T_c(2:end, :) - T_c(1:end-1, :)) ./ bt_distance;
    dT_dt1 =  (T_c(:, 2:end) - T_c(:, 1:end-1)) ./ lr_distance;

    %% Corrected normal temperature derivative at face 2
    divide2 = Ct2_Fn2 .* Cn2_Ft2 - Cn2_Fn2 .* Ct2_Ft2;
    dT_dn2 = (Cn2_Ft2 .* dT_dt2 - dT_dc2 .* Ct2_Ft2) ./ divide2;

    %% Corrected normal temperature derivative at face 1
    divide1 = Cn1_Fn1 .* Ct1_Ft1 - Ct1_Fn1 .* Cn1_Ft1;
    dT_dn1 = (Ct1_Ft1 .* dT_dc1 - dT_dt1 .* Cn1_Ft1) ./ divide1;

    %% Assemble heat flux and add to energy equation
    q_n2 = k_bt .* dT_dn2;
    q_n1 = k_lr .* dT_dn1;

    rho_E_temp = (q_n2(:, 2:end) .* s.mesh.bt_area(:, 2:end) - q_n2(:, 1:end-1) .* s.mesh.bt_area(:, 1:end-1)) ...
               - (q_n1(2:end, :) .* s.mesh.lr_area(2:end, :) - q_n1(1:end-1, :) .* s.mesh.lr_area(1:end-1, :));
    s.flux.rho_E = s.flux.rho_E + rho_E_temp ./ s.mesh.volume;
end
