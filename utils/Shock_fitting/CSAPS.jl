# CSAPS.jl - BSplineKit derivative helpers for shock spline differentiation
#
#   Provides _spline_derivative, a helper that computes first derivatives
#   for BSplineKit spline objects using the Derivative operator, so
#   COMPUTE_BETA can transparently differentiate the shock spline.
#
#   BSplineKit treats derivatives as first-class operators: applying
#   Derivative(1) * spline produces a new spline representing the exact
#   analytical derivative of the original curve.
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

using BSplineKit

"""
    _spline_derivative(spl, xi)

Evaluate the first derivative of a BSplineKit spline at points `xi`.

Uses the BSplineKit `Derivative` operator to construct the exact analytical
derivative spline, then evaluates it at the requested points.
"""
function _spline_derivative(spl, xi::AbstractVector)
    spl_deriv = Derivative(1) * spl
    return spl_deriv.(xi)
end

function _spline_derivative(spl, xi::Real)
    spl_deriv = Derivative(1) * spl
    return spl_deriv(xi)
end
