# LEAST_SQUARES_SHOCK_POINTS - Fit and resample shock points using spline interpolation.
#
#   s = LEAST_SQUARES_SHOCK_POINTS(s)
#
#   Fits a cubic smoothing spline through the current shock point
#   coordinates using BSplineKit.jl, re-evaluates the spline at the grid
#   parametric locations (chi), and converts the new Cartesian coordinates
#   back to the element-space representation (eta, chi).
#
#   The smoothing parameter uses the MATLAB csaps convention:
#     min  p Σ(yi - f(xi))² + (1-p) ∫(f''(x))² dx
#   Internally converted to BSplineKit's λ = (1 - p) / p.
#
#   Inputs:
#       s (Dict) - Solution structure containing at minimum:
#                    s["shock"]["spline_param"]   - Smoothing parameter,
#                                                   uses MATLAB csaps convention: p ∈ [0, 1]
#                    s["shock"]["points_chi"]     - (Nx x 1) parametric coordinates of shock points
#                    s["shock"]["points_x"]/["points_y"]     - (Nx x 1) Cartesian shock coordinates
#                    s["mesh"]["chi"]             - (Nx x Ny) grid parametric coordinates
#
#   Outputs:
#       s (Dict) - Updated s with:
#                    s["shock"]["spline_func_x"]/["spline_func_y"] - Stored BSplineKit spline objects
#                    s["shock"]["points_x"]/["points_y"]           - Resampled Cartesian coordinates
#                    s["shock"]["points_eta"]                      - Radial coordinates in element space
#                    s["shock"]["points_chi"]                      - Updated parametric coordinates
#
#   See also: GO_TO_ELEMENT_SPACE
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

using BSplineKit
using StatsAPI

function LEAST_SQUARES_SHOCK_POINTS(s::Dict{String, Any})
    ## Extract dict refs
    shock = s["shock"]::Dict{String, Any}
    mesh  = s["mesh"]::Dict{String, Any}

    ## Fit cubic smoothing spline to shock points using BSplineKit
    # Convert MATLAB csaps parameter p to BSplineKit λ.
    #   MATLAB csaps minimizes:  p Σ(yi-f(xi))² + (1-p) ∫(f''(x))² dx
    #   BSplineKit minimizes:    Σ(yi-f(xi))² + λ ∫(f''(x))² dx
    #   Dividing MATLAB's objective by p ⟹ λ = (1 - p) / p
    p = shock["spline_param"]
    λ = (1 - p) / p

    chi = Float64.(vec(shock["points_chi"]))

    # 1. Fit the base cubic smoothing splines
    spl_x_base = StatsAPI.fit(BSplineOrder(4), chi, Float64.(vec(shock["points_x"])), λ)
    spl_y_base = StatsAPI.fit(BSplineOrder(4), chi, Float64.(vec(shock["points_y"])), λ)

    # 2. Wrap them to automatically extrapolate outside the 'chi' bounds.
    # Smooth() continues the cubic polynomial. Use Linear() if Smooth() curves too wildly.
    spl_x = BSplineKit.extrapolate(spl_x_base, BSplineKit.Linear())
    spl_y = BSplineKit.extrapolate(spl_y_base, BSplineKit.Linear())

    # 3. Store the extrapolated versions in your dictionary
    shock["spline_func_x"] = spl_x
    shock["spline_func_y"] = spl_y

    ## Resample shock points at grid parametric locations
    chi_grid = mesh["chi"]::Matrix{Float64}
    @views shock["points_x"] = spl_x.(chi_grid[:, 1])
    @views shock["points_y"] = spl_y.(chi_grid[:, 1])

    ## Convert back to element space
    shock["points_chi"], shock["points_eta"] = GO_TO_ELEMENT_SPACE(shock["points_x"], shock["points_y"], s)

    return s
end
