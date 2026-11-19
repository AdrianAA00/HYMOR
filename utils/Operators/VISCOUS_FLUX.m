function s = VISCOUS_FLUX(s)
% VISCOUS_FLUX  Compute viscous flux contributions for momentum and energy equations.
%
%   s = VISCOUS_FLUX(s)
%
%   Evaluates the viscous stress tensor and the heat flux on every face of
%   a general (non-orthogonal) quadrilateral mesh and adds their face
%   integrals to the momentum and energy fluxes. Velocity and temperature
%   gradients come from FACE_GRADIENTS (exact for linear fields). The
%   viscous stress work and Fourier heat conduction are added to the
%   energy flux.
%
%   For PDE_dimension = "3D-axisymmetric" (y = r, v the radial velocity)
%   the face areas and volumes of the mesh carry the radius, the divergence
%   gains v/r and the hoop stress tau_theta_theta = mu (2 v/r - 2/3 div)
%   enters the radial momentum as the source -tau_theta_theta / r.
%
%   Inputs:
%       s - Solution struct containing:
%                  .var.rho, .var.rho_u, .var.rho_v - Conserved variables
%                  .var.mu_star, .var.k_star       - Viscosity and conductivity
%                  .var.T                          - Temperature field
%                  .mesh.*                         - Centroids, corners, face
%                                                    normals, areas, volumes
%                  .PDE_dimension                  - "2D" or "3D-axisymmetric"
%                  .freestream.*                   - Non-dimensionalization parameters
%
%   Outputs:
%       s - Updated s struct with viscous contributions added
%                  to .flux.rho_u, .flux.rho_v, and .flux.rho_E
%
%   Notes:
%       - The viscous tensor uses the deviatoric form with the 2/3
%         divergence term.
%       - v/r on a face is the mean of v/r in the two adjacent cells, so it
%         stays finite on faces lying on the symmetry axis; the hoop stress
%         uses v/r of the cell itself.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    mesh = s.mesh;
    Re_infty = s.freestream.Re;
    axisymmetric = s.PDE_dimension == "3D-axisymmetric";

    %% Face gradients of velocity and temperature
    u = s.var.rho_u ./ s.var.rho;
    v = s.var.rho_v ./ s.var.rho;
    [gu, gv, gT] = FACE_GRADIENTS(mesh, u, v, s.var.T);

    %% Face-averaged viscosity (non-dimensionalized)
    mu_lr = FACE_LR(s.var.mu_star) / Re_infty;
    mu_bt = FACE_BT(s.var.mu_star) / Re_infty;

    %% Velocity divergence on the faces
    div_lr = gu.x_lr + gv.y_lr;
    div_bt = gu.x_bt + gv.y_bt;
    if axisymmetric
        v_over_r = v ./ mesh.y_Ext;      % cell centroids, ghosts included, never on r = 0
        div_lr = div_lr + FACE_LR(v_over_r);
        div_bt = div_bt + FACE_BT(v_over_r);
    end

    %% Viscous traction tau . n on every face (normal towards increasing index)
    [tx_lr, ty_lr] = TRACTION(mu_lr, gu.x_lr, gu.y_lr, gv.x_lr, gv.y_lr, div_lr, mesh.lr_x_normal, mesh.lr_y_normal);
    [tx_bt, ty_bt] = TRACTION(mu_bt, gu.x_bt, gu.y_bt, gv.x_bt, gv.y_bt, div_bt, mesh.bt_x_normal, mesh.bt_y_normal);

    %% Viscous momentum flux
    s.flux.rho_u = s.flux.rho_u + FACE_DIVERGENCE(tx_lr .* mesh.lr_area, tx_bt .* mesh.bt_area, mesh);
    s.flux.rho_v = s.flux.rho_v + FACE_DIVERGENCE(ty_lr .* mesh.lr_area, ty_bt .* mesh.bt_area, mesh);

    %% Hoop stress (axisymmetric): tau_tt = mu (2 v/r - 2/3 div) in the cell
    if axisymmetric
        tau_theta_theta = HOOP_STRESS(s.var.mu_star(2:end-1, 2:end-1) / Re_infty, gu, gv, ...
                                      v_over_r(2:end-1, 2:end-1));
        s.flux.rho_v = s.flux.rho_v - tau_theta_theta ./ mesh.y;
    end

    %% Viscous work (stress power) contribution to energy flux
    rho_lr = FACE_LR(s.var.rho);
    rho_bt = FACE_BT(s.var.rho);
    u_lr = FACE_LR(s.var.rho_u) ./ rho_lr;
    u_bt = FACE_BT(s.var.rho_u) ./ rho_bt;
    v_lr = FACE_LR(s.var.rho_v) ./ rho_lr;
    v_bt = FACE_BT(s.var.rho_v) ./ rho_bt;
    s.flux.rho_E = s.flux.rho_E + FACE_DIVERGENCE((u_lr .* tx_lr + v_lr .* ty_lr) .* mesh.lr_area, ...
                                                  (u_bt .* tx_bt + v_bt .* ty_bt) .* mesh.bt_area, mesh);

    %% Heat flux (Fourier's law): k grad(T) . n
    k_scale = s.freestream.gamma_star / s.freestream.Re / s.freestream.Pr;
    k_lr = FACE_LR(s.var.k_star) * k_scale;
    k_bt = FACE_BT(s.var.k_star) * k_scale;
    q_lr = k_lr .* (gT.x_lr .* mesh.lr_x_normal + gT.y_lr .* mesh.lr_y_normal) .* mesh.lr_area;
    q_bt = k_bt .* (gT.x_bt .* mesh.bt_x_normal + gT.y_bt .* mesh.bt_y_normal) .* mesh.bt_area;
    s.flux.rho_E = s.flux.rho_E + FACE_DIVERGENCE(q_lr, q_bt, mesh);
end

%% ========================================================================
%  Viscous Traction
%  ========================================================================
function [tx, ty] = TRACTION(mu, ux, uy, vx, vy, div, nx, ny)
% TRACTION  Components of tau . n for tau = mu (grad u + grad u^T - 2/3 div I).
    txx = mu .* (2 * ux - 2/3 * div);
    tyy = mu .* (2 * vy - 2/3 * div);
    txy = mu .* (uy + vx);
    tx = txx .* nx + txy .* ny;
    ty = txy .* nx + tyy .* ny;
end
