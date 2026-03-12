# COMPUTE_UPSTREAM_CONDITIONS - Evaluate freestream primitive variables at shock points.
#
#   (p_infty, u_inf, v_inf, rho_inf, e_inf) = COMPUTE_UPSTREAM_CONDITIONS(s)
#
#   Retrieves the upstream (pre-shock) conservative variables from
#   UPDATE_SHOCK_UPSTREAM, converts them to primitive form (density,
#   velocity components, internal energy), and computes the upstream
#   pressure using the equation of state.
#
#   Inputs:
#       s (Dict) - Solution structure containing upstream flow parameters,
#                    freestream conditions, perturbation amplitudes, and
#                    shock-point coordinates (see UPDATE_SHOCK_UPSTREAM).
#
#   Outputs:
#       p_infty (Nx x 1) - Upstream pressure at each shock point.
#       u_inf   (Nx x 1) - Upstream x-velocity at each shock point.
#       v_inf   (Nx x 1) - Upstream y-velocity at each shock point.
#       rho_inf (Nx x 1) - Upstream density at each shock point.
#       e_inf   (Nx x 1) - Upstream specific internal energy at each
#                                  shock point.
#
#   Notes:
#       - Pressure is computed from the calorically perfect gas relation:
#         p = rho * e * (gamma_star - 1).
#
#   See also: UPDATE_SHOCK_UPSTREAM (local subfunction)
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

function COMPUTE_UPSTREAM_CONDITIONS(s::Dict{String, Any})
    rho, rho_u, rho_v, rho_E = UPDATE_SHOCK_UPSTREAM(s)

    fs = s["freestream"]::Dict{String, Any}
    gamma_star = fs["gamma_star"]::Float64
    gm1 = gamma_star - 1.0

    rho_inf = rho
    @inbounds begin
        v_inf = @. rho_v / rho_inf
        u_inf = @. rho_u / rho_inf
        u_mag = @. sqrt(v_inf^2 + u_inf^2)
        e_inf = @. (rho_E - rho_inf * u_mag^2 * 0.5) / rho_inf
        p_infty = @. e_inf * rho_inf * gm1
    end

    return p_infty, u_inf, v_inf, rho_inf, e_inf
end

# UPDATE_SHOCK_UPSTREAM - Build upstream conservative state with optional perturbations.
#
#   Computes upstream density, momentum, and total energy at the shock
#   points, either from a prescribed perturbed field (rho_0_upstream_p)
#   or by superimposing a travelling-wave perturbation on the uniform
#   freestream state. Enforces axisymmetric boundary conditions when
#   the dimension flag is set to "3D-axisymmetric".

function UPDATE_SHOCK_UPSTREAM(s::Dict{String, Any})
    ## Extract dict refs
    fs    = s["freestream"]::Dict{String, Any}
    shock = s["shock"]::Dict{String, Any}
    ti    = s["time_integration"]::Dict{String, Any}

    rho_0   = fs["rho_0"]::Float64
    rho_u_0 = fs["rho_u_0"]::Float64
    rho_v_0 = fs["rho_v_0"]::Float64
    rho_E_0 = fs["rho_E_0"]::Float64

    disturbance = fs["disturbance"]::Dict{String, Any}
    k_y       = disturbance["k_y"]::Float64
    k_x       = disturbance["k_x"]::Float64
    amplitude = disturbance["amplitude"]

    pts_x = shock["points_x"]
    pts_y = shock["points_y"]
    t     = ti["t"]::Float64

    u_y = rho_v_0 / rho_0
    u_x = rho_u_0 / rho_0
    two_pi_ky = 2.0 * pi * k_y
    two_pi_kx = 2.0 * pi * k_x
    u_y_t = u_y * t
    u_x_t = u_x * t
    @inbounds perturbation = @. cos(two_pi_ky * (pts_y - u_y_t)) *
                                sin(two_pi_kx * (pts_x - u_x_t))

    if haskey(s, "rho_0_upstream_p")
        rho   = fs["rho_0_p"]
        rho_u = fs["rho_u_0_p"]
        rho_v = fs["rho_v_0_p"]
        rho_E = fs["rho_E_0_p"]
    else
        amp1 = amplitude[1]
        amp2 = amplitude[2]
        amp3 = amplitude[3]
        amp4 = amplitude[4]
        @inbounds begin
            rho   = @. rho_0   + perturbation * amp1
            rho_u = @. rho_u_0 + perturbation * amp2
            rho_v = @. rho_v_0 + perturbation * amp3
            rho_E = @. rho_E_0 + perturbation * amp4
        end
    end

    ## Apply symmetry boundary conditions for axisymmetric flows
    pde_dim = s["PDE_dimension"]::String
    if pde_dim == "3D-axisymmetric"
        @inbounds rho_v[end, 1] = 0
    end

    return rho, rho_u, rho_v, rho_E
end
