function [du_dx, du_dy] = DERIVATIVE(u, s)
% DERIVATIVE - Compute spatial derivatives on a structured 2D grid.
%
%   Computes the partial derivatives du/dx and du/dy of a scalar field u
%   using central finite differences between cell centroids. The derivatives
%   are first computed in the local cell-aligned curvilinear coordinate
%   system (c1, c2) and then transformed to the Cartesian (x, y) system
%   via a coordinate Jacobian.
%
% Syntax:
%   [du_dx, du_dy] = DERIVATIVE(u, s)
%
% Inputs:
%   u        - 2D array of the scalar field (including ghost cells).
%   s - Structure containing at minimum:
%                .mesh.x, .mesh.y        - Cell centroid coordinates.
%                .shock.enabled          - Logical flag; true if a shock is present.
%                .mesh.Nchi              - Number of interior cells in x.
%                .shock.cell_indices     - (Nx+2 x 1) array of shock-cell j-indices.
%
% Outputs:
%   du_dx - 2D array of du/dx on interior cells (size: Nx-2 x Ny-2).
%   du_dy - 2D array of du/dy on interior cells (size: Nx-2 x Ny-2).
%
% Notes:
%   - Uses second-order central differencing on cell centroid distances.
%   - When a shock is present, derivatives in and beyond shocked cells are
%     zeroed out to avoid differencing across the discontinuity.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% Compute centroid-to-centroid distances
    centroids_bt_distance = sqrt((s.mesh.y(2:end-1,3:end) - s.mesh.y(2:end-1,1:end-2)).^2 + ...
                                 (s.mesh.x(2:end-1,3:end) - s.mesh.x(2:end-1,1:end-2)).^2);
    centroids_lr_distance = sqrt((s.mesh.y(3:end,2:end-1) - s.mesh.y(1:end-2,2:end-1)).^2 + ...
                                 (s.mesh.x(3:end,2:end-1) - s.mesh.x(1:end-2,2:end-1)).^2);

    %% Compute derivatives in cell-aligned coordinate system
    du_dc1 = -(u(3:end,2:end-1) - u(1:end-2,2:end-1)) ./ centroids_lr_distance;
    du_dc2 =  (u(2:end-1,3:end) - u(2:end-1,1:end-2)) ./ centroids_bt_distance;

    %% Project cell coordinate system onto Cartesian (x, y)
    c1_x = -(s.mesh.x(3:end,2:end-1) - s.mesh.x(1:end-2,2:end-1)) ./ centroids_lr_distance;
    c1_y = -(s.mesh.y(3:end,2:end-1) - s.mesh.y(1:end-2,2:end-1)) ./ centroids_lr_distance;
    c2_x =  (s.mesh.x(2:end-1,3:end) - s.mesh.x(2:end-1,1:end-2)) ./ centroids_bt_distance;
    c2_y =  (s.mesh.y(2:end-1,3:end) - s.mesh.y(2:end-1,1:end-2)) ./ centroids_bt_distance;

    %% Transform derivatives to Cartesian system
    divide = c1_x .* c2_y - c2_x .* c1_y;
    du_dx = (c2_y .* du_dc1 - c1_y .* du_dc2) ./ divide;
    du_dy = (c1_x .* du_dc2 - c2_x .* du_dc1) ./ divide;

    %% Zero out derivatives at and beyond shocked cells
    if s.shock.enabled
        for i = 1:s.mesh.Nchi-2
            du_dx(i, s.shock.cell_indices(i+1,1)-1:end) = 0;
            du_dy(i, s.shock.cell_indices(i+1,1)-1:end) = 0;
        end
    end
end
