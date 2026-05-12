# CFL_TIMESTEP  Compute an adaptive time step based on the CFL condition.
#
#   s = CFL_TIMESTEP(s)
#
#   Evaluates the combined convective + viscous spectral radius per cell on
#   the structured grid and returns the most restrictive stable time step.
#   The final dt is the minimum of:
#     1) Combined convective + viscous CFL limit (per-cell Blazek form),
#     2) Flux-magnitude limit on each conserved variable,
#     3) User-specified maximum time step (s["time_integration"]["max_dt"]).
#
#   Inputs:
#       s (Dict{String,Any}) - Solution structure containing conservative
#                              variables (s["var"].*), grid metrics (s["mesh"].lr_area,
#                              s["mesh"].bt_area, etc.), sound speed (s["var"]["a"]),
#                              local gamma, viscosity, Prandtl number,
#                              CFL number (s["time_integration"]["CFL"]),
#                              and freestream reference quantities.
#
#   Outputs:
#       s (Dict{String,Any}) - Same dict with s["time_integration"]["dt"] updated
#                              to the new stable time step.
#
#   Notes:
#       - Convective spectral radius uses the proper face-normal velocity
#         |u*n_x + v*n_y| (not the conservative bound |u||n_x|+|v||n_y|).
#       - Viscous spectral radius follows Blazek's compressible-NS form:
#         lambda_v = max(4/3, gamma/Pr) / (rho * Re_flow) * (A_xi^2 + A_eta^2) / V,
#         which includes the 1/rho factor (kinematic viscosity), the 4/3
#         deviatoric-stress prefactor, and the gamma/Pr thermal-diffusion
#         eigenvalue (whichever is larger).
#       - Convective and viscous contributions are summed per cell before
#         taking the global minimum dt (the proper combined explicit-CFL
#         constraint), rather than min(dt_conv, dt_visc).
#       - For shock-fitted simulations (s["shock"]["enabled"] == true) the
#         limits are masked by s["shock"]["flow_cells"] so that only active
#         flow cells contribute.
#
# Part of: Hypersonics Stability Julia Solver - Time Marching Module

function CFL_TIMESTEP(s::Dict{String,Any})

    ## Extract arrays with type assertions for compiler optimization
    var  = s["var"]::Dict{String,Any}
    mesh = s["mesh"]::Dict{String,Any}
    ti   = s["time_integration"]::Dict{String,Any}
    shock = s["shock"]::Dict{String,Any}
    flux  = s["flux"]::Dict{String,Any}

    rho_arr     = var["rho"]::Matrix{Float64}
    rho_u_arr   = var["rho_u"]::Matrix{Float64}
    rho_v_arr   = var["rho_v"]::Matrix{Float64}
    rho_E_arr   = var["rho_E"]::Matrix{Float64}
    a_arr       = var["a"]::Matrix{Float64}
    gamma_arr   = var["gamma_star"]::Matrix{Float64}
    Re_flow_arr = var["Re_flow"]::Matrix{Float64}
    Pr_flow_arr = var["Pr_flow"]::Matrix{Float64}

    lr_xn  = mesh["lr_x_normal"]::Matrix{Float64}
    lr_yn  = mesh["lr_y_normal"]::Matrix{Float64}
    lr_a   = mesh["lr_area"]::Matrix{Float64}
    bt_xn  = mesh["bt_x_normal"]::Matrix{Float64}
    bt_yn  = mesh["bt_y_normal"]::Matrix{Float64}
    bt_a   = mesh["bt_area"]::Matrix{Float64}
    vol    = mesh["volume"]::Matrix{Float64}

    flow_cells = shock["flow_cells"]::Matrix{Float64}
    shock_enabled = shock["enabled"]::Bool
    CFL_num    = ti["CFL"]::Float64
    max_dt     = ti["max_dt"]::Float64

    flux_rho   = flux["rho"]::Matrix{Float64}
    flux_rho_u = flux["rho_u"]::Matrix{Float64}
    flux_rho_v = flux["rho_v"]::Matrix{Float64}
    flux_rho_E = flux["rho_E"]::Matrix{Float64}

    @views begin

    ## Sound speed estimate
    c = maximum(abs, a_arr)

    ## Convective spectral radius (per cell)
    # Left-right face-centered velocities (mass-weighted average of neighbours)
    rho_sum_lr = @. rho_arr[1:end-1, 2:end-1] + rho_arr[2:end, 2:end-1]
    u_face_lr  = @. (rho_u_arr[1:end-1, 2:end-1] + rho_u_arr[2:end, 2:end-1]) / rho_sum_lr
    v_face_lr  = @. (rho_v_arr[1:end-1, 2:end-1] + rho_v_arr[2:end, 2:end-1]) / rho_sum_lr
    Vn_lr       = @. abs(u_face_lr * lr_xn + v_face_lr * lr_yn)
    lam_face_lr = @. (Vn_lr + c) * lr_a

    # Bottom-top face-centered velocities
    rho_sum_bt = @. rho_arr[2:end-1, 1:end-1] + rho_arr[2:end-1, 2:end]
    u_face_bt  = @. (rho_u_arr[2:end-1, 1:end-1] + rho_u_arr[2:end-1, 2:end]) / rho_sum_bt
    v_face_bt  = @. (rho_v_arr[2:end-1, 1:end-1] + rho_v_arr[2:end-1, 2:end]) / rho_sum_bt
    Vn_bt       = @. abs(u_face_bt * bt_xn + v_face_bt * bt_yn)
    lam_face_bt = @. (Vn_bt + c) * bt_a

    # Collapse opposing faces to per-cell xi/eta spectral radii (max of the two)
    lam_xi  = @. max(lam_face_lr[1:end-1, :], lam_face_lr[2:end, :])
    lam_eta = @. max(lam_face_bt[:, 1:end-1], lam_face_bt[:, 2:end])
    lam_c   = @. lam_xi + lam_eta

    ## Viscous spectral radius (per cell, Blazek-style)
    rho_int   = rho_arr[2:end-1, 2:end-1]
    gamma_int = gamma_arr[2:end-1, 2:end-1]
    Re_flow   = Re_flow_arr[2:end-1, 2:end-1]
    Pr_flow   = Pr_flow_arr[2:end-1, 2:end-1]

    A_xi_cell  = @. 0.5 * (lr_a[1:end-1, :] + lr_a[2:end, :])
    A_eta_cell = @. 0.5 * (bt_a[:, 1:end-1] + bt_a[:, 2:end])

    nu_eff = @. max(4.0/3.0, gamma_int / Pr_flow) / (rho_int * Re_flow)
    lam_v  = @. 2 * nu_eff * (A_xi_cell^2 + A_eta_cell^2) / vol

    ## Combined convective + viscous time step (per-cell sum, global min)
    lam_over_V = @. (lam_c + lam_v) * flow_cells / vol
    dt1 = CFL_num / maximum(lam_over_V)

    ## Flux-magnitude time step
    rho_int_flux   = rho_arr[2:end-1, 2:end-1]
    rho_E_int_flux = rho_E_arr[2:end-1, 2:end-1]
    if shock_enabled
        max_val_u   = maximum(@. (abs(flux_rho_v) + abs(flux_rho_u)) / rho_int_flux * flow_cells)
        max_val_rho = maximum(@. abs(flux_rho / rho_int_flux * flow_cells))
        max_val_E   = maximum(@. abs(flux_rho_E / rho_E_int_flux * flow_cells))
    else
        max_val_u   = maximum(@. (abs(flux_rho_v) + abs(flux_rho_u)) / rho_int_flux)
        max_val_rho = maximum(@. abs(flux_rho / rho_int_flux))
        max_val_E   = maximum(@. abs(flux_rho_E / rho_E_int_flux))
    end
    dt2 = CFL_num / max(max_val_u, max_val_rho, max_val_E)

    ## Final time step: minimum of all constraints
    ti["dt"] = min(dt1, dt2, max_dt)

    end # @views

    return s
end
