function REFINEMENT(y, boundary_layer_thickness, wall_coordinate, refinement_intensity, max_val, min_val)
# REFINEMENT - Apply hyperbolic-tangent grid clustering near a wall coordinate.
#
# Transforms a uniform coordinate distribution into a non-uniform one
# clustered around a specified wall coordinate using an antisymmetric
# tanh-based stretching function.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Apply tanh-based stretching
    y_scaled = @. y - tanh(refinement_intensity * boundary_layer_thickness * (y - wall_coordinate)) / boundary_layer_thickness
    y_max = max_val - tanh(refinement_intensity * boundary_layer_thickness * (max_val - wall_coordinate)) / boundary_layer_thickness
    y_min = min_val - tanh(refinement_intensity * boundary_layer_thickness * (min_val - wall_coordinate)) / boundary_layer_thickness

    ## Rescale to desired output range
    y_range = y_max - y_min
    @. y_scaled = (y_scaled - y_min) / y_range * (max_val - min_val) + min_val

    return y_scaled
end
