function s = FLUX_AZIMUTHAL_MODE(s)
% FLUX_AZIMUTHAL_MODE  Add the azimuthal-mode terms of the linearized
%   axisymmetric Navier-Stokes equations to the fluxes.
%
%   s = FLUX_AZIMUTHAL_MODE(s)
%
%   For a swirl-free axisymmetric base flow (u0, v0, 0) and a perturbation
%   with azimuthal wavenumber m,
%       (rho', rho_u', rho_v', rho_E')  =  hat(r, x) * cos(m*theta),
%       rho_w'                          =  hat(r, x) * sin(m*theta),
%   the linear operator is real. The cos-amplitude rows (mass, x/r momentum,
%   energy) receive the theta-derivatives of the azimuthal fluxes; the
%   sin-amplitude row (azimuthal momentum, s.flux.rho_w) receives the full
%   azimuthal momentum equation:
%       d(rho w)/dt = -(1/r^2) d_r(r^2 F_r) - d_x F_x  +  (m/r) (p' - tau'_tt)
%   with F = rho w u - (tau_xt, tau_rt) written in angular-momentum form.
%
%   Every term is written as a nonlinear expression of the current state
%   that vanishes on the base flow and is exact to first order in the
%   perturbation (theta-derivatives: d_theta(cos-type q) -> -m (q - q0),
%   d_theta(w) -> m w). The finite-difference linearization therefore picks
%   up exactly the linear azimuthal operator. The thermodynamic closure
%   (p, T, mu, k) comes from the regular state update, so every gas model
%   is supported.
%
%   Inputs:
%       s - Solution structure in linearization mode, with s.var.rho_w and
%           s.azimuthal (base fields incl. the hoop stress, m) from SETUP_AZIMUTHAL_MODE;
%           boundary conditions already applied.
%
%   Outputs:
%       s - Same structure with s.flux.rho_w set and the azimuthal terms
%           added to s.flux.rho, .rho_u, .rho_v and .rho_E.
%
%   Notes:
%       - x is the axial and y the radial coordinate (y >= 0).
%       - For m = 0 only the decoupled swirl equation is evaluated: the
%         cos-amplitude rows are left untouched (bitwise).
%       - Face gradients, face means and the hoop stress use the helpers of
%         VISCOUS_FLUX (FACE_GRADIENTS, FACE_LR/FACE_BT, HOOP_STRESS).
%       - Quotients by r that are singular on the axis for odd m appear
%         only in regular combinations, (v + m w)/r and (w + m v)/r, or are
%         multiplied by the face radius (angular-momentum form).
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    mesh = s.mesh;
    base = s.azimuthal.base;
    m = s.azimuthal.m;
    Re = s.freestream.Re;
    r = mesh.y;                                                  % cell radius
    r_lr = (mesh.y_corner(:, 1:end-1) + mesh.y_corner(:, 2:end)) / 2;   % face radii
    r_bt = (mesh.y_corner(1:end-1, :) + mesh.y_corner(2:end, :)) / 2;
    nx_lr = mesh.lr_x_normal;  ny_lr = mesh.lr_y_normal;
    nx_bt = mesh.bt_x_normal;  ny_bt = mesh.bt_y_normal;

    %% Primitive fields on the extended grid (ghost cells included)
    rho = s.var.rho;
    rho_w = s.var.rho_w;
    u = s.var.rho_u ./ rho;
    v = s.var.rho_v ./ rho;
    w = rho_w ./ rho;
    mu = s.var.mu_star / Re;
    du = u - base.u;
    dv = v - base.v;

    %% Face values and in-plane gradients
    [gu, gv, gw] = FACE_GRADIENTS(mesh, u, v, w);
    mu_lr = FACE_LR(mu);  mu_bt = FACE_BT(mu);
    un_lr = (FACE_LR(s.var.rho_u) .* nx_lr + FACE_LR(s.var.rho_v) .* ny_lr) ./ FACE_LR(rho);
    un_bt = (FACE_BT(s.var.rho_u) .* nx_bt + FACE_BT(s.var.rho_v) .* ny_bt) ./ FACE_BT(rho);

    %% Azimuthal momentum (sin amplitude), angular-momentum form
    % r_f * (tau . n)_theta with r tau_xt = mu (r dw/dx - m du),
    %                            r tau_rt = mu (r dw/dr - w - m dv)
    G_lr = r_lr .* (mu_lr .* (gw.x_lr .* nx_lr + gw.y_lr .* ny_lr) - FACE_LR(rho_w) .* un_lr) ...
         - mu_lr .* (FACE_LR(w) .* ny_lr + m * (FACE_LR(du) .* nx_lr + FACE_LR(dv) .* ny_lr));
    G_bt = r_bt .* (mu_bt .* (gw.x_bt .* nx_bt + gw.y_bt .* ny_bt) - FACE_BT(rho_w) .* un_bt) ...
         - mu_bt .* (FACE_BT(w) .* ny_bt + m * (FACE_BT(du) .* nx_bt + FACE_BT(dv) .* ny_bt));
    s.flux.rho_w = FACE_DIVERGENCE(G_lr .* mesh.lr_area, G_bt .* mesh.bt_area, mesh) ./ r;

    if m == 0
        return
    end

    I = {2:size(rho, 1) - 1, 2:size(rho, 2) - 1};
    tau_tt = HOOP_STRESS(mu(I{:}), gu, gv, (v(I{:}) + m * w(I{:})) ./ r);
    dp = s.var.p(I{:}) - base.p;
    s.flux.rho_w = s.flux.rho_w + m ./ r .* (dp - (tau_tt - base.tau_tt));

    %% Cos-amplitude rows
    m_r = m ./ r;
    rho_w_c = rho_w(I{:});
    u_c = u(I{:});
    v_c = v(I{:});
    w_c = w(I{:});
    mu_c = mu(I{:});
    H_c = (s.var.rho_E(I{:}) + s.var.p(I{:})) ./ rho(I{:});
    tau_xt = mu_c .* (CELL_AVERAGE(gw.x_lr, gw.x_bt) - m_r .* du(I{:}));     % sin amplitudes
    tau_rt = mu_c .* (CELL_AVERAGE(gw.y_lr, gw.y_bt) - (w_c + m * dv(I{:})) ./ r);

    % Isotropic part of the stress from (1/r) d_theta(w) = D: -2/3 mu D on the
    % in-plane faces; its share 2 mu D - 2/3 mu D of tau_tt (hoop source)
    % completes the v/r part of VISCOUS_FLUX
    D = m * w ./ mesh.y_Ext;
    sigma_lr = -2/3 * mu_lr .* FACE_LR(D);
    sigma_bt = -2/3 * mu_bt .* FACE_BT(D);

    s.flux.rho = s.flux.rho - m_r .* rho_w_c;
    s.flux.rho_u = s.flux.rho_u + m_r .* (tau_xt - rho_w_c .* u_c) ...
        + FACE_DIVERGENCE(sigma_lr .* nx_lr .* mesh.lr_area, sigma_bt .* nx_bt .* mesh.bt_area, mesh);
    s.flux.rho_v = s.flux.rho_v + m_r .* (tau_rt - rho_w_c .* v_c) - 4/3 * mu_c .* D(I{:}) ./ r ...
        + FACE_DIVERGENCE(sigma_lr .* ny_lr .* mesh.lr_area, sigma_bt .* ny_bt .* mesh.bt_area, mesh);

    k_c = s.var.k_star(I{:}) * s.freestream.gamma_star / Re / s.freestream.Pr;
    dT = s.var.T(I{:}) - base.T;
    work_lr = sigma_lr .* (FACE_LR(u) .* nx_lr + FACE_LR(v) .* ny_lr) .* mesh.lr_area;
    work_bt = sigma_bt .* (FACE_BT(u) .* nx_bt + FACE_BT(v) .* ny_bt) .* mesh.bt_area;
    s.flux.rho_E = s.flux.rho_E ...
        + m_r .* (u_c .* tau_xt + v_c .* tau_rt + tau_tt .* w_c - rho_w_c .* H_c) ...
        - m_r.^2 .* k_c .* dT + FACE_DIVERGENCE(work_lr, work_bt, mesh);
end
