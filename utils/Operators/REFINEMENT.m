function y_scaled = REFINEMENT(y, boundary_layer_thickness, wall_coordinate, refinement_intensity, max_val, min_val)
% REFINEMENT  Apply hyperbolic-tangent grid clustering near a wall coordinate.
%
%   y_scaled = REFINEMENT(y, boundary_layer_thickness, wall_coordinate,
%                         refinement_intensity, max_val, min_val)
%
%   Transforms a uniform coordinate distribution into a non-uniform one
%   clustered around a specified wall coordinate using an antisymmetric
%   tanh-based stretching function. The output is rescaled to fit within
%   the specified [min_val, max_val] range.
%
%   Inputs:
%       y                        - Original coordinate array (uniform or otherwise)
%       boundary_layer_thickness - Controls the width of the clustering region
%                                  (higher values produce tighter clustering)
%       wall_coordinate          - Coordinate around which clustering is centered
%       refinement_intensity     - Stretching strength, 0 = no refinement,
%                                  1 = maximum refinement
%       max_val                  - Upper bound of the output coordinate range
%       min_val                  - Lower bound of the output coordinate range
%
%   Outputs:
%       y_scaled - Refined coordinate array rescaled to [min_val, max_val]
%
%   Notes:
%       - The clustering is antisymmetric about wall_coordinate.
%       - When refinement_intensity = 0, the output is a linearly rescaled
%         version of the input_file.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% Apply tanh-based stretching
    y_scaled = y - tanh(refinement_intensity * boundary_layer_thickness * (y - wall_coordinate)) / boundary_layer_thickness;
    y_max = max_val - tanh(refinement_intensity * boundary_layer_thickness * (max_val - wall_coordinate)) / boundary_layer_thickness;
    y_min = min_val - tanh(refinement_intensity * boundary_layer_thickness * (min_val - wall_coordinate)) / boundary_layer_thickness;

    %% Rescale to desired output range
    y_range = y_max - y_min;
    y_scaled = (y_scaled - y_min) ./ y_range .* (max_val - min_val) + min_val;
end
