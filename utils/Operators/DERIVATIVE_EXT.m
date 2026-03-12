function [du_dx, du_dy] = DERIVATIVE_EXT(u, s)
% DERIVATIVE_EXT - Compute spatial derivatives on the extended grid.
%
%   Computes the partial derivatives du/dx and du/dy of a scalar field u
%   using central finite differences between cell centroids on the extended
%   grid (x_Ext, y_Ext). The derivatives are first computed in the local
%   cell-aligned curvilinear coordinate system (c1, c2) and then
%   transformed to the Cartesian (x, y) system via a coordinate Jacobian.
%
%   When a shock is present, the derivative values near the shock front are
%   replaced with second-order extrapolations from cells farther away to
%   avoid differencing across the discontinuity.
%
% Syntax:
%   [du_dx, du_dy] = DERIVATIVE_EXT(u, s)
%
% Inputs:
%   u        - 2D array of the scalar field (including ghost cells).
%   s - Structure containing at minimum:
%                .mesh.x_Ext, .mesh.y_Ext - Extended-grid cell centroid coordinates.
%                .shock.enabled          - Logical flag; true if a shock is present.
%                .mesh.Nchi              - Number of interior cells in x.
%                .shock.cell_indices     - (Nx x 1) array of shock-cell j-indices.
%
% Outputs:
%   du_dx - 2D array of du/dx on interior cells.
%   du_dy - 2D array of du/dy on interior cells.
%
% Notes:
%   - Uses second-order central differencing on cell centroid distances.
%   - Near-shock extrapolation uses a linear (2nd-order) scheme to maintain
%     accuracy while avoiding the discontinuity jump.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% Compute centroid-to-centroid distances on extended grid
    centroids_bt_distance = sqrt((s.mesh.y_Ext(2:end-1,3:end) - s.mesh.y_Ext(2:end-1,1:end-2)).^2 + ...
                                 (s.mesh.x_Ext(2:end-1,3:end) - s.mesh.x_Ext(2:end-1,1:end-2)).^2);
    centroids_lr_distance = sqrt((s.mesh.y_Ext(3:end,2:end-1) - s.mesh.y_Ext(1:end-2,2:end-1)).^2 + ...
                                 (s.mesh.x_Ext(3:end,2:end-1) - s.mesh.x_Ext(1:end-2,2:end-1)).^2);

    %% Compute derivatives in cell-aligned coordinate system
    du_dc1 = -(u(3:end,2:end-1) - u(1:end-2,2:end-1)) ./ centroids_lr_distance;
    du_dc2 =  (u(2:end-1,3:end) - u(2:end-1,1:end-2)) ./ centroids_bt_distance;

    %% Project cell coordinate system onto Cartesian (x, y)
    c1_x = -(s.mesh.x_Ext(3:end,2:end-1) - s.mesh.x_Ext(1:end-2,2:end-1)) ./ centroids_lr_distance;
    c1_y = -(s.mesh.y_Ext(3:end,2:end-1) - s.mesh.y_Ext(1:end-2,2:end-1)) ./ centroids_lr_distance;
    c2_x =  (s.mesh.x_Ext(2:end-1,3:end) - s.mesh.x_Ext(2:end-1,1:end-2)) ./ centroids_bt_distance;
    c2_y =  (s.mesh.y_Ext(2:end-1,3:end) - s.mesh.y_Ext(2:end-1,1:end-2)) ./ centroids_bt_distance;

    %% Transform derivatives to Cartesian system
    divide = c1_x .* c2_y - c2_x .* c1_y;
    du_dx = (c2_y .* du_dc1 - c1_y .* du_dc2) ./ divide;
    du_dy = (c1_x .* du_dc2 - c2_x .* du_dc1) ./ divide;

    %% Extrapolate derivatives near shocked cells
    if s.shock.enabled
        % Compute linear indices for shocked cells and their neighbors
        i_values_s = (1:s.mesh.Nchi)';
        j_values_s = s.shock.cell_indices(:,1);
        idx_3 = sub2ind(size(du_dy), i_values_s, j_values_s + 3);
        idx_2 = sub2ind(size(du_dy), i_values_s, j_values_s + 2);
        idx_1 = sub2ind(size(du_dy), i_values_s, j_values_s + 1);
        idx_0 = sub2ind(size(du_dy), i_values_s, j_values_s);

        % Extrapolate du_dx near shock (2nd-order linear extrapolation)
        temp_x = du_dx;
        du_dx(idx_1) = 2 * temp_x(idx_2) - temp_x(idx_3);
        du_dx(idx_0) = 2 * du_dx(idx_1) - du_dx(idx_2);

        % Extrapolate du_dy near shock (2nd-order linear extrapolation)
        temp_y = du_dy;
        du_dy(idx_1) = 2 * temp_y(idx_2) - temp_y(idx_3);
        du_dy(idx_0) = 2 * du_dy(idx_1) - du_dy(idx_2);
    end
end
