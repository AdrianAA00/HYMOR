# COMPUTE_BETA - Compute the local shock wave angle at each sampled point.
#
#   s = COMPUTE_BETA(s)
#
#   Evaluates the first derivative of the shock spline in both x and y
#   directions, then combines with the upstream flow direction to obtain
#   the local oblique-shock angle beta at every shock sample point.
#
#   Inputs:
#       s (Dict) - Solution structure containing at minimum:
#                    s["shock"]["spline_func_x"] - Spline for shock x-coordinates
#                    s["shock"]["spline_func_y"] - Spline for shock y-coordinates
#                    s["mesh"]["chi"]            - (Nx x Ny) parametric coordinate of shock points
#                    (plus fields required by COMPUTE_UPSTREAM_CONDITIONS)
#
#   Outputs:
#       s (Dict) - Updated s with:
#                    s["shock"]["beta"] - (Nx x 1) local shock wave angle [rad]
#
#   Notes:
#       - The shock angle beta is defined as the angle between the upstream
#         velocity vector and the shock-tangent vector.
#       - Uses _spline_derivative which applies BSplineKit's Derivative(1)
#         operator to obtain an exact analytical derivative spline (see CSAPS.jl).
#
#   See also: COMPUTE_UPSTREAM_CONDITIONS, _spline_derivative
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

function COMPUTE_BETA(s::Dict{String, Any})
    ## Extract dict refs
    shock = s["shock"]::Dict{String, Any}
    mesh  = s["mesh"]::Dict{String, Any}
    chi   = mesh["chi"]::Matrix{Float64}

    ## Compute shock tangent derivatives
    @views deriv_x = _spline_derivative(shock["spline_func_x"], chi[:, 1])
    @views deriv_y = _spline_derivative(shock["spline_func_y"], chi[:, 1])

    ## Compute shock angle (independent of chi traversal direction)
    _, u_inf, v_inf, _, _ = COMPUTE_UPSTREAM_CONDITIONS(s)

    # Cross and dot products between shock tangent and upstream velocity
    @inbounds begin
        cross_tv = @. deriv_x * v_inf - deriv_y * u_inf
        dot_tv   = @. deriv_x * u_inf + deriv_y * v_inf

        # Acute angle between shock tangent line and upstream velocity
        shock["beta"] = @. atan(abs(cross_tv), abs(dot_tv))
    end

    return s
end
