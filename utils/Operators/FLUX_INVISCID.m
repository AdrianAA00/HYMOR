function s = FLUX_INVISCID(s)
% FLUX_INVISCID  Inviscid mass, momentum and total-energy fluxes.
%
%   s = FLUX_INVISCID(s)
%
%   Central finite-volume fluxes in divergence form. With averages of the
%   two cells adjacent to each face, avg(.), and the face velocity
%   v_n = (avg(rho_u) n_x + avg(rho_v) n_y) / avg(rho):
%       mass:      avg(rho_u) n_x + avg(rho_v) n_y
%       momentum:  avg(rho_u_i) v_n + avg(p) n_i
%       energy:    (avg(rho_E) + avg(p)) v_n
%   The face averages are formed once and shared by the four equations.
%   For PDE_dimension = "3D-axisymmetric" the face areas and volumes carry
%   the radius and the front/back faces of the wedge add the hoop pressure
%   to the radial momentum (the azimuthal velocity is zero).
%
%   Inputs:
%       s - Solution struct with .var (rho, rho_u, rho_v, rho_E, p on the
%           extended grid), .mesh (normals, areas, volumes) and
%           .PDE_dimension.
%
%   Outputs:
%       s - Updated s struct with .flux.rho, .flux.rho_u, .flux.rho_v and
%           .flux.rho_E set to the inviscid contributions.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    mesh = s.mesh;

    %% Face averages and face velocity
    rho_u_lr = FACE_LR(s.var.rho_u);   rho_u_bt = FACE_BT(s.var.rho_u);
    rho_v_lr = FACE_LR(s.var.rho_v);   rho_v_bt = FACE_BT(s.var.rho_v);
    p_lr = FACE_LR(s.var.p);           p_bt = FACE_BT(s.var.p);
    mass_lr = rho_u_lr .* mesh.lr_x_normal + rho_v_lr .* mesh.lr_y_normal;
    mass_bt = rho_u_bt .* mesh.bt_x_normal + rho_v_bt .* mesh.bt_y_normal;
    vn_lr = mass_lr ./ FACE_LR(s.var.rho);
    vn_bt = mass_bt ./ FACE_BT(s.var.rho);

    %% Net outflow per unit volume, sign flipped into the right-hand side
    s.flux.rho = -FACE_DIVERGENCE(mass_lr .* mesh.lr_area, mass_bt .* mesh.bt_area, mesh);
    s.flux.rho_u = -FACE_DIVERGENCE((rho_u_lr .* vn_lr + p_lr .* mesh.lr_x_normal) .* mesh.lr_area, ...
                                    (rho_u_bt .* vn_bt + p_bt .* mesh.bt_x_normal) .* mesh.bt_area, mesh);
    s.flux.rho_v = -FACE_DIVERGENCE((rho_v_lr .* vn_lr + p_lr .* mesh.lr_y_normal) .* mesh.lr_area, ...
                                    (rho_v_bt .* vn_bt + p_bt .* mesh.bt_y_normal) .* mesh.bt_area, mesh);
    s.flux.rho_E = -FACE_DIVERGENCE((FACE_LR(s.var.rho_E) + p_lr) .* vn_lr .* mesh.lr_area, ...
                                    (FACE_BT(s.var.rho_E) + p_bt) .* vn_bt .* mesh.bt_area, mesh);

    %% Hoop pressure on the front/back faces (axisymmetric)
    if s.PDE_dimension == "3D-axisymmetric"
        p_cell = s.var.p(2:end-1, 2:end-1);
        s.flux.rho_v = s.flux.rho_v + p_cell * (mesh.fb_y_normal_front - mesh.fb_y_normal_back) ...
            .* mesh.fb_area .* sign(mesh.y) ./ mesh.volume;
    end
end
