# CFL_TIMESTEP  Compute an adaptive time step based on the CFL condition.
#
#   s = CFL_TIMESTEP(s)
#
#   Evaluates convective and viscous CFL constraints on the structured
#   grid to determine a stable time step. The final dt is the minimum of:
#     1) Convective CFL limit (left-right and south-north face sweeps),
#     2) Flux-magnitude limit on each conserved variable,
#     3) Viscous diffusion limit (Reynolds and Prandtl-scaled),
#     4) User-specified maximum time step (s["time_integration"]["max_dt"]).
#
#   Inputs:
#       s (Dict{String,Any}) - Solution structure containing conservative
#                              variables (s["var"].*), grid metrics (s["mesh"].lr_area,
#                              s["mesh"].bt_area, etc.), sound speed (s["var"]["a"]),
#                              viscosity, Prandtl number, CFL number
#                              (s["time_integration"]["CFL"]), and freestream reference
#                              quantities.
#
#   Outputs:
#       s (Dict{String,Any}) - Same dict with s["time_integration"]["dt"] updated to the
#                              new stable time step.
#
#   Notes:
#       - For shock-fitted simulations (s["shock"]["enabled"] == true) the flux-
#         magnitude limit is masked by s["shock"]["flow_cells"] so that only
#         active flow cells contribute.
#       - The viscous limit uses the minimum centroid distance between
#         south-north and left-right neighbours.
#
# Part of: Hypersonics Stability Julia Solver - Time Marching Module

function CFL_TIMESTEP(s::Dict{String,Any})

    ## Extract arrays with type assertions for compiler optimization
    var  = s["var"]::Dict{String,Any}
    mesh = s["mesh"]::Dict{String,Any}
    ti   = s["time_integration"]::Dict{String,Any}
    shock = s["shock"]::Dict{String,Any}
    flux  = s["flux"]::Dict{String,Any}

    rho_arr   = var["rho"]::Matrix{Float64}
    rho_u_arr = var["rho_u"]::Matrix{Float64}
    rho_v_arr = var["rho_v"]::Matrix{Float64}
    rho_E_arr = var["rho_E"]::Matrix{Float64}
    a_arr     = var["a"]::Matrix{Float64}
    Re_flow_arr = var["Re_flow"]::Matrix{Float64}
    Pr_flow_arr = var["Pr_flow"]::Matrix{Float64}

    lr_xn  = mesh["lr_x_normal"]::Matrix{Float64}
    lr_yn  = mesh["lr_y_normal"]::Matrix{Float64}
    lr_a   = mesh["lr_area"]::Matrix{Float64}
    bt_xn  = mesh["bt_x_normal"]::Matrix{Float64}
    bt_yn  = mesh["bt_y_normal"]::Matrix{Float64}
    bt_a   = mesh["bt_area"]::Matrix{Float64}
    vol    = mesh["volume"]::Matrix{Float64}
    x_Ext  = mesh["x_Ext"]::Matrix{Float64}
    y_Ext  = mesh["y_Ext"]::Matrix{Float64}

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

    ## Left-Right convective CFL
    u_lr = @. abs((rho_u_arr[1:end-1, 2:end-1] + rho_u_arr[2:end, 2:end-1]) /
                  (rho_arr[2:end, 2:end-1] + rho_arr[1:end-1, 2:end-1]) * lr_xn)
    v_lr = @. abs((rho_v_arr[1:end-1, 2:end-1] + rho_v_arr[2:end, 2:end-1]) /
                  (rho_arr[2:end, 2:end-1] + rho_arr[1:end-1, 2:end-1]) * lr_yn)
    U_lr     = @. u_lr + v_lr
    U_lr_vol = @. max(U_lr[1:end-1, :], U_lr[2:end, :])
    A_lr     = @. max(lr_a[1:end-1, :], lr_a[2:end, :])
    U_lr_A   = @. (U_lr_vol + c) * A_lr
    U_lr_A_eig = maximum(@. U_lr_A * flow_cells / vol)

    ## South-North convective CFL
    u_bt = @. abs((rho_u_arr[2:end-1, 1:end-1] + rho_u_arr[2:end-1, 2:end]) /
                  (rho_arr[2:end-1, 2:end] + rho_arr[2:end-1, 1:end-1]) * bt_xn)
    v_bt = @. abs((rho_v_arr[2:end-1, 1:end-1] + rho_v_arr[2:end-1, 2:end]) /
                  (rho_arr[2:end-1, 2:end] + rho_arr[2:end-1, 1:end-1]) * bt_yn)
    U_bt     = @. u_bt + v_bt
    U_bt_vol = @. max(U_bt[:, 1:end-1], U_bt[:, 2:end])
    A_bt     = @. max(bt_a[:, 1:end-1], bt_a[:, 2:end])
    U_bt_A   = @. (U_bt_vol + c) * A_bt
    U_bt_A_eig = maximum(@. U_bt_A * flow_cells / vol)

    ## Convective time step
    dt1 = CFL_num / (U_lr_A_eig + U_bt_A_eig)

    ## Viscous time step
    centroids_bt_distance = @. sqrt((y_Ext[2:end-1, 2:end] - y_Ext[2:end-1, 1:end-1])^2 +
                                    (x_Ext[2:end-1, 2:end] - x_Ext[2:end-1, 1:end-1])^2)
    centroids_lr_distance = @. sqrt((y_Ext[2:end, 2:end-1] - y_Ext[1:end-1, 2:end-1])^2 +
                                    (x_Ext[2:end, 2:end-1] - x_Ext[1:end-1, 2:end-1])^2)

    Re_flow = Re_flow_arr[2:end-1, 2:end-1]
    Pr_flow = Pr_flow_arr[2:end-1, 2:end-1]
    d = @. min(centroids_bt_distance[:, 2:end], centroids_lr_distance[2:end, :])

    one_Re    = maximum(@. (1.0 / Re_flow) / d^2 * flow_cells) * 4
    one_Re_Pr = maximum(@. (1.0 / (Re_flow * Pr_flow)) / d^2 * flow_cells) * 4
    dt1_visc  = min(CFL_num / one_Re, CFL_num / one_Re_Pr)

    ## Flux-magnitude time step
    rho_int   = rho_arr[2:end-1, 2:end-1]
    rho_E_int = rho_E_arr[2:end-1, 2:end-1]
    if shock_enabled
        max_val_u   = maximum(@. (abs(flux_rho_v) + abs(flux_rho_u)) / rho_int * flow_cells)
        max_val_rho = maximum(@. abs(flux_rho / rho_int * flow_cells))
        max_val_E   = maximum(@. abs(flux_rho_E / rho_E_int * flow_cells))
    else
        max_val_u   = maximum(@. (abs(flux_rho_v) + abs(flux_rho_u)) / rho_int)
        max_val_rho = maximum(@. abs(flux_rho / rho_int))
        max_val_E   = maximum(@. abs(flux_rho_E / rho_E_int))
    end
    dt2 = CFL_num / max(max_val_u, max_val_rho, max_val_E)

    ## Final time step: minimum of all constraints
    ti["dt"] = min(dt1, dt2, dt1_visc, max_dt)

    end # @views

    return s
end
