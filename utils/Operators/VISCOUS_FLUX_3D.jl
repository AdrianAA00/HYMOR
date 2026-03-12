function VISCOUS_FLUX_3D(s::Dict{String, Any})
# VISCOUS_FLUX_3D - Compute 3D axisymmetric viscous flux contributions for momentum and energy.
#
# Extends VISCOUS_FLUX to axisymmetric (3D) geometries by including the
# hoop stress (tau_theta_theta) and the geometric u/r source term in the
# divergence. Velocity and temperature gradients are computed on a general
# non-orthogonal quadrilateral mesh with non-orthogonal correction.
# The radial (x) direction carries the axisymmetric correction terms.
#
# Inputs:
#     s - Solution dictionary containing:
#                .var.rho, .var.rho_u, .var.rho_v, .var.rho_E - Conserved variables
#                .mesh.x_Ext, .mesh.y_Ext       - Extended cell centroid coordinates
#                .mesh.x_corner, .mesh.y_corner - Cell corner coordinates
#                .mesh.bt_x_normal, .mesh.bt_y_normal - Bottom-top face unit normals
#                .mesh.lr_x_normal, .mesh.lr_y_normal - Left-right face unit normals
#                .mesh.bt_area, .mesh.lr_area   - Face areas
#                .mesh.volume                   - Cell volumes
#                .var.mu_star                   - Dynamic viscosity field
#                .var.k_star                    - Thermal conductivity field
#                .var.T                         - Temperature field
#                .freestream.*                  - Non-dimensionalization parameters
#
# Outputs:
#     s - Updated dictionary with axisymmetric viscous
#                contributions added to .flux.rho_u, .flux.rho_v, .flux.rho_E
#
# Notes:
#     - The divergence is modified: div = du_t/dt + du_n/dn + u/r
#       where r is the radial coordinate (x-direction).
#     - Hoop stress tau_theta_theta = mu*(2*u/r - 2/3*div) is added as
#       a source term to the radial momentum equation.
#     - Heat flux uses non-orthogonal correction identical to VISCOUS_FLUX.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract arrays to avoid repeated dictionary lookups (with type assertions)
    rho_arr     = s["var"]["rho"]::Matrix{Float64}
    rho_u_arr   = s["var"]["rho_u"]::Matrix{Float64}
    rho_v_arr   = s["var"]["rho_v"]::Matrix{Float64}
    mu_star     = s["var"]["mu_star"]::Matrix{Float64}
    k_star_arr  = s["var"]["k_star"]::Matrix{Float64}
    T_arr       = s["var"]["T"]::Matrix{Float64}
    x_Ext       = s["mesh"]["x_Ext"]::Matrix{Float64}
    y_Ext       = s["mesh"]["y_Ext"]::Matrix{Float64}
    x_corner    = s["mesh"]["x_corner"]::Matrix{Float64}
    y_corner    = s["mesh"]["y_corner"]::Matrix{Float64}
    lr_xn       = s["mesh"]["lr_x_normal"]::Matrix{Float64}
    lr_yn       = s["mesh"]["lr_y_normal"]::Matrix{Float64}
    lr_a        = s["mesh"]["lr_area"]::Matrix{Float64}
    bt_xn       = s["mesh"]["bt_x_normal"]::Matrix{Float64}
    bt_yn       = s["mesh"]["bt_y_normal"]::Matrix{Float64}
    bt_a        = s["mesh"]["bt_area"]::Matrix{Float64}
    vol         = s["mesh"]["volume"]::Matrix{Float64}
    flux_rho_u  = s["flux"]["rho_u"]::Matrix{Float64}
    flux_rho_v  = s["flux"]["rho_v"]::Matrix{Float64}
    flux_rho_E  = s["flux"]["rho_E"]::Matrix{Float64}
    Re_infty    = s["freestream"]["Re"]::Float64
    gamma_star_fs = s["freestream"]["gamma_star"]::Float64
    Pr_fs       = s["freestream"]["Pr"]::Float64

    ## Pre-calculate velocity fields at cell centers
    u = @. rho_u_arr / rho_arr
    v = @. rho_v_arr / rho_arr

    ## Velocity values at cell corners (4-point average)
    @views u_c = @. (u[1:end-1, 1:end-1] + u[2:end, 1:end-1] +
                     u[1:end-1, 2:end]   + u[2:end, 2:end]) / 4
    @views v_c = @. (v[1:end-1, 1:end-1] + v[2:end, 1:end-1] +
                     v[1:end-1, 2:end]   + v[2:end, 2:end]) / 4

    ## Distances between cell centroids
    @views centroids_bt_distance = @. sqrt((y_Ext[2:end-1, 2:end] - y_Ext[2:end-1, 1:end-1])^2 +
                                           (x_Ext[2:end-1, 2:end] - x_Ext[2:end-1, 1:end-1])^2)
    @views centroids_lr_distance = @. sqrt((y_Ext[2:end, 2:end-1] - y_Ext[1:end-1, 2:end-1])^2 +
                                           (x_Ext[2:end, 2:end-1] - x_Ext[1:end-1, 2:end-1])^2)

    ## Edge lengths of cell faces
    @views bt_distance = @. sqrt((x_corner[2:end, :] - x_corner[1:end-1, :])^2 +
                                 (y_corner[2:end, :] - y_corner[1:end-1, :])^2)
    @views lr_distance = @. sqrt((x_corner[:, 2:end] - x_corner[:, 1:end-1])^2 +
                                 (y_corner[:, 2:end] - y_corner[:, 1:end-1])^2)

    ## Unit vector components for non-orthogonal correction
    @views c1_x = @. -(x_Ext[2:end, 2:end-1] - x_Ext[1:end-1, 2:end-1]) / centroids_lr_distance
    @views c1_y = @. -(y_Ext[2:end, 2:end-1] - y_Ext[1:end-1, 2:end-1]) / centroids_lr_distance
    @views c2_x = @.  (x_Ext[2:end-1, 2:end] - x_Ext[2:end-1, 1:end-1]) / centroids_bt_distance
    @views c2_y = @.  (y_Ext[2:end-1, 2:end] - y_Ext[2:end-1, 1:end-1]) / centroids_bt_distance

    ## Projection matrices for reference system transformation
    Cn2_Fn2 = @. c2_x * bt_xn + c2_y * bt_yn
    Cn2_Ft2 = @. c2_x * bt_yn - c2_y * bt_xn

    Cn1_Fn1 = @. -c1_x * lr_xn - c1_y * lr_yn
    Cn1_Ft1 = @.  c1_x * lr_yn - c1_y * lr_xn

    ## Velocity derivatives along cell-center lines
    @views du_dc2 = @.  (u[2:end-1, 2:end] - u[2:end-1, 1:end-1]) / centroids_bt_distance
    @views dv_dc2 = @.  (v[2:end-1, 2:end] - v[2:end-1, 1:end-1]) / centroids_bt_distance
    @views du_dc1 = @. -(u[2:end, 2:end-1] - u[1:end-1, 2:end-1]) / centroids_lr_distance
    @views dv_dc1 = @. -(v[2:end, 2:end-1] - v[1:end-1, 2:end-1]) / centroids_lr_distance

    ## Tangential velocity derivatives at faces
    @views du_dt2 = @. -(u_c[2:end, :] - u_c[1:end-1, :]) / bt_distance
    @views dv_dt2 = @. -(v_c[2:end, :] - v_c[1:end-1, :]) / bt_distance
    @views du_dt1 = @.  (u_c[:, 2:end] - u_c[:, 1:end-1]) / lr_distance
    @views dv_dt1 = @.  (v_c[:, 2:end] - v_c[:, 1:end-1]) / lr_distance

    ## Normal derivatives at face 2 (bottom-top faces)
    # Original: divide = Ct2_Fn2 .* Cn2_Ft2 .- Cn2_Fn2 .* Ct2_Ft2
    # With Ct2_Fn2=0 and Ct2_Ft2=1: divide = -Cn2_Fn2
    divide = @. -Cn2_Fn2
    du_dn2 = @. (Cn2_Ft2 * du_dt2 - du_dc2) / divide
    dv_dn2 = @. (Cn2_Ft2 * dv_dt2 - dv_dc2) / divide

    ## Transform to local coordinate system (n2, t2) at face 2
    dut_dt2 = @. du_dt2 * bt_yn - dv_dt2 * bt_xn
    dun_dt2 = @. du_dt2 * bt_xn + dv_dt2 * bt_yn
    dut_dn2 = @. du_dn2 * bt_yn - dv_dn2 * bt_xn
    dun_dn2 = @. du_dn2 * bt_xn + dv_dn2 * bt_yn

    ## Normal derivatives at face 1 (left-right faces)
    # Original: divide = Cn1_Fn1 .* Ct1_Ft1 .- Ct1_Fn1 .* Cn1_Ft1
    # With Ct1_Fn1=0 and Ct1_Ft1=1: divide = Cn1_Fn1
    divide1 = @. Cn1_Fn1
    du_dn1 = @. (du_dc1 - du_dt1 * Cn1_Ft1) / divide1
    dv_dn1 = @. (dv_dc1 - dv_dt1 * Cn1_Ft1) / divide1

    ## Transform to local coordinate system (n1, t1) at face 1
    dun_dn1 = @. -du_dn1 * lr_xn - dv_dn1 * lr_yn
    dut_dn1 = @.  du_dn1 * lr_yn - dv_dn1 * lr_xn
    dun_dt1 = @. -du_dt1 * lr_xn - dv_dt1 * lr_yn
    dut_dt1 = @.  du_dt1 * lr_yn - dv_dt1 * lr_xn

    ## Face-averaged viscosity (non-dimensionalized)
    @views mu_1 = @. (mu_star[2:end, 2:end-1] + mu_star[1:end-1, 2:end-1]) / 2
    @views mu_2 = @. (mu_star[2:end-1, 2:end] + mu_star[2:end-1, 1:end-1]) / 2

    ## Axisymmetric geometry: radial coordinates and velocities at faces
    @views x_face_1 = @. (x_Ext[2:end, 2:end-1] + x_Ext[1:end-1, 2:end-1]) / 2
    @views x_face_2 = @. (x_Ext[2:end-1, 2:end] + x_Ext[2:end-1, 1:end-1]) / 2

    @views u_face_1 = @. (u[2:end, 2:end-1] + u[1:end-1, 2:end-1]) / 2
    @views u_face_2 = @. (u[2:end-1, 2:end] + u[2:end-1, 1:end-1]) / 2

    ## Viscous stress tensor with axisymmetric divergence correction (u/r term)
    div_1 = @. dut_dt1 + dun_dn1 + u_face_1 / x_face_1
    div_2 = @. dut_dt2 + dun_dn2 + u_face_2 / x_face_2
    tau_11 = @. mu_1 / Re_infty * (2 * dun_dn1 - (2/3) * div_1)
    tau_12 = @. mu_1 / Re_infty * (dun_dt1 + dut_dn1)
    tau_22 = @. mu_2 / Re_infty * (2 * dun_dn2 - (2/3) * div_2)
    tau_21 = @. mu_2 / Re_infty * (dut_dn2 + dun_dt2)

    ## Hoop stress for axisymmetric correction
    tau_theta_theta_1 = @. mu_1 / Re_infty * (2 * u_face_1 / x_face_1 - (2/3) * div_1)
    tau_theta_theta_2 = @. mu_2 / Re_infty * (2 * u_face_2 / x_face_2 - (2/3) * div_2)

    # Average hoop stress to cell centers
    @views tau_theta_theta_cell = @. (tau_theta_theta_1[1:end-1, :] + tau_theta_theta_1[2:end, :] +
                                      tau_theta_theta_2[:, 1:end-1] + tau_theta_theta_2[:, 2:end]) / 4

    # Radial coordinate at cell centers
    @views x_cell = x_Ext[2:end-1, 2:end-1]

    ## Viscous stress in global (x-y) coordinates
    tau_1x = @. -(tau_11 * lr_xn - tau_12 * lr_yn)
    tau_1y = @. -(tau_11 * lr_yn + tau_12 * lr_xn)
    tau_2x = @.  (tau_22 * bt_xn + tau_21 * bt_yn)
    tau_2y = @.  (tau_22 * bt_yn - tau_21 * bt_xn)

    ## Viscous momentum flux (x-direction) with hoop stress source term
    @views rho_u_temp = @. (tau_1x[1:end-1, :] * lr_a[1:end-1, :] - tau_1x[2:end, :] * lr_a[2:end, :]) +
                           (tau_2x[:, 2:end]   * bt_a[:, 2:end]   - tau_2x[:, 1:end-1] * bt_a[:, 1:end-1])
    @. flux_rho_u += rho_u_temp / vol

    # Add hoop stress source term to radial momentum
    @. flux_rho_u -= tau_theta_theta_cell / x_cell

    ## Viscous momentum flux (y-direction)
    @views rho_v_temp = @. (tau_1y[1:end-1, :] * lr_a[1:end-1, :] - tau_1y[2:end, :] * lr_a[2:end, :]) +
                           (tau_2y[:, 2:end]   * bt_a[:, 2:end]   - tau_2y[:, 1:end-1] * bt_a[:, 1:end-1])
    @. flux_rho_v += rho_v_temp / vol

    ## Viscous work (stress power) contribution to energy flux
    @views inv_rho_1 = @. 1 / (rho_arr[1:end-1, 2:end-1] + rho_arr[2:end, 2:end-1])
    @views inv_rho_2 = @. 1 / (rho_arr[2:end-1, 1:end-1] + rho_arr[2:end-1, 2:end])

    @views u_1 = @. (rho_u_arr[1:end-1, 2:end-1] + rho_u_arr[2:end, 2:end-1]) * inv_rho_1
    @views u_2 = @. (rho_u_arr[2:end-1, 1:end-1] + rho_u_arr[2:end-1, 2:end]) * inv_rho_2
    @views v_1 = @. (rho_v_arr[1:end-1, 2:end-1] + rho_v_arr[2:end, 2:end-1]) * inv_rho_1
    @views v_2 = @. (rho_v_arr[2:end-1, 1:end-1] + rho_v_arr[2:end-1, 2:end]) * inv_rho_2

    w1 = @. u_1 * tau_1x + v_1 * tau_1y
    w2 = @. u_2 * tau_2x + v_2 * tau_2y

    @views rho_E_temp = @. (w1[1:end-1, :] * lr_a[1:end-1, :] - w1[2:end, :] * lr_a[2:end, :]) +
                           (w2[:, 2:end]   * bt_a[:, 2:end]   - w2[:, 1:end-1] * bt_a[:, 1:end-1])
    @. flux_rho_E += rho_E_temp / vol

    ## Heat flux with non-orthogonal correction (Fourier's law)
    @views k_bt = @. (k_star_arr[2:end-1, 2:end] + k_star_arr[2:end-1, 1:end-1]) / 2 * gamma_star_fs / Re_infty / Pr_fs
    @views k_lr = @. (k_star_arr[2:end, 2:end-1] + k_star_arr[1:end-1, 2:end-1]) / 2 * gamma_star_fs / Re_infty / Pr_fs

    ## Temperature at cell corners (4-point average)
    @views T_c = @. (T_arr[1:end-1, 1:end-1] + T_arr[2:end, 1:end-1] +
                     T_arr[1:end-1, 2:end]   + T_arr[2:end, 2:end]) / 4

    ## Temperature derivatives along cell-center lines
    @views dT_dc2 = @.  (T_arr[2:end-1, 2:end] - T_arr[2:end-1, 1:end-1]) / centroids_bt_distance
    @views dT_dc1 = @. -(T_arr[2:end, 2:end-1] - T_arr[1:end-1, 2:end-1]) / centroids_lr_distance

    ## Tangential temperature derivatives at faces
    @views dT_dt2 = @. -(T_c[2:end, :] - T_c[1:end-1, :]) / bt_distance
    @views dT_dt1 = @.  (T_c[:, 2:end] - T_c[:, 1:end-1]) / lr_distance

    ## Corrected normal temperature derivative at face 2
    dT_dn2 = @. (Cn2_Ft2 * dT_dt2 - dT_dc2) / divide

    ## Corrected normal temperature derivative at face 1
    dT_dn1 = @. (dT_dc1 - dT_dt1 * Cn1_Ft1) / divide1

    ## Assemble heat flux and add to energy equation
    q_n2 = @. k_bt * dT_dn2
    q_n1 = @. k_lr * dT_dn1

    @views rho_E_temp2 = @. (q_n2[:, 2:end] * bt_a[:, 2:end] - q_n2[:, 1:end-1] * bt_a[:, 1:end-1]) -
                            (q_n1[2:end, :] * lr_a[2:end, :] - q_n1[1:end-1, :] * lr_a[1:end-1, :])
    @. flux_rho_E += rho_E_temp2 / vol

    return s
end
