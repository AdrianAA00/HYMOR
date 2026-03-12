# MOVE_SHOCK - Advance shock point positions by the computed shock speed.
#
#   s = MOVE_SHOCK(s, dsolution_dt, w)
#
#   Updates the shock point coordinates by adding the product of the
#   shock speed, the Runge-Kutta weight, and the user-defined shock
#   relaxation factor. This is called during the time-integration
#   loop to evolve the shock location.
#
#   Inputs:
#       s     (Dict) - Solution structure containing:
#                               s["shock"]["points_x"]/["points_y"]    - (Nx x 1) current shock coordinates
#                               s["shock"]["relaxation"]    - (scalar) under-relaxation factor
#       dsolution_dt (Dict) - Time-derivative structure containing:
#                               dsolution_dt["shock"]["speed_x"]/["speed_y"]     - (Nx x 1) shock velocity components
#       w            (Float64) - Runge-Kutta stage weight.
#
#   Outputs:
#       s     (Dict) - Updated s with new shock coordinates.
#
#   Notes:
#       - The relaxation factor allows under-relaxation of the shock
#         motion for stability (shock_relaxation < 1 damps oscillations).
#
#   See also: UPDATE_SHOCK_JUMP_PROPERTIES
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

function MOVE_SHOCK(s::Dict{String, Any}, dsolution_dt::Dict{String, Any}, w::Float64)
    ## Extract dict refs
    shock     = s["shock"]::Dict{String, Any}
    ds_shock  = dsolution_dt["shock"]::Dict{String, Any}
    relax     = shock["relaxation"]::Float64
    w_relax   = w * relax

    pts_x  = shock["points_x"]
    pts_y  = shock["points_y"]
    spd_x  = ds_shock["speed_x"]
    spd_y  = ds_shock["speed_y"]

    @inbounds @. pts_x = pts_x + spd_x * w_relax
    @inbounds @. pts_y = pts_y + spd_y * w_relax

    return s
end
