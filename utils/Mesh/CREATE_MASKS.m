function [mask_no_stagnation, mask_no_BL, mask_only_BL] = CREATE_MASKS(s)
% CREATE_MASKS - Generate spatial masks for stagnation and boundary-layer regions.
%
%   [mask_no_stagnation, mask_no_BL, mask_only_BL] = CREATE_MASKS(s)
%
%   Creates smooth sigmoid-based masks that separate the computational
%   domain into stagnation region, boundary-layer region, and outer flow
%   region. These masks are used for selective application of numerical
%   treatments (e.g., artificial dissipation, limiting) in different
%   parts of the flow field.
%
%   Inputs:
%       s  (struct) - Solution structure containing:
%           .curvilinear_mapping.boundary_type - Geometry type ("MSL", "blunt_cone", etc.)
%           .curvilinear_mapping.L             - Reference body length
%           .curvilinear_mapping.R             - Nose radius
%           .curvilinear_mapping.theta         - Half-cone angle [rad]
%           .mesh.chi           - (Nchi x Neta) element-space chi coordinates
%           .mesh.Nchi, .mesh.Neta - Number of cells in each direction
%           .var.T              - (Nchi+2 x Neta+2) temperature field with ghosts
%           .shock.flow_cells   - (Nchi x Neta) logical mask of active flow cells
%
%   Outputs:
%       mask_no_stagnation  (Nx x 1)  - Sigmoid mask that is ~0 near the
%                                        stagnation point and ~1 elsewhere
%       mask_no_BL          (Nx x Ny) - Sigmoid mask that is ~0 inside the
%                                        thermal boundary layer and ~1 outside
%       mask_only_BL        (Nx x Ny) - Complement of mask_no_BL (1 - mask_no_BL)
%
%   Notes:
%       - For "MSL" and "blunt_cone" geometries, the stagnation mask uses
%         the surface coordinate chi and a logistic sigmoid function.
%       - The boundary-layer mask is based on the thermal field: cells
%         with temperature below 80% of the mean are considered outside
%         the boundary layer.
%       - For all other geometry types, trivial masks are returned (no
%         stagnation exclusion, no BL exclusion).
%
% Part of: Hypersonics Stability MATLAB Solver - Mesh Module

    if s.curvilinear_mapping.boundary_type == "MSL" || s.curvilinear_mapping.boundary_type == "blunt_cone"
        %% Geometry parameters
        s_Straight = ((s.curvilinear_mapping.L - s.curvilinear_mapping.R) + s.curvilinear_mapping.R * cos(s.curvilinear_mapping.theta)) / sin(s.curvilinear_mapping.theta);
        s_Curve = s.curvilinear_mapping.theta * s.curvilinear_mapping.R;
        s_tot = s_Curve + s_Straight;

        %% Stagnation mask (sigmoid along surface coordinate)
        amount = 0.5;
        length = 0.01;
        coordinate = (1 - 2 * s.mesh.chi - amount * s_Curve / s_tot) / length;
        mask_no_stagnation = exp(coordinate) ./ (1 + exp(coordinate));

        %% Boundary-layer mask (sigmoid based on thermal field)
        length = 0.02;
        T_field = s.var.T(2:end-1, 2:end-1) .* s.shock.flow_cells;
        mean_T = mean(T_field(T_field ~= 0));
        coordinate = (T_field - mean_T * 0.8) ./ mean_T / length;
        mask_no_BL = exp(coordinate) ./ (1 + exp(coordinate)) .* s.shock.flow_cells ...
                     + ones(s.mesh.Nchi, s.mesh.Neta) - s.shock.flow_cells;

        mask_only_BL = ones(s.mesh.Nchi, s.mesh.Neta) - mask_no_BL;
    else
        %% Trivial masks for non-blunt geometries
        mask_no_stagnation = ones(s.mesh.Nchi, s.mesh.Neta);
        mask_only_BL = zeros(s.mesh.Nchi, s.mesh.Neta);
        mask_no_BL = ones(s.mesh.Nchi, s.mesh.Neta);
    end

end
