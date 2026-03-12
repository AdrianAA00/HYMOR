function s = GET_EXTENDED_CELL_CENTERS(s)
% GET_EXTENDED_CELL_CENTERS  Compute cell centers and build extended grid.
%
%   s = GET_EXTENDED_CELL_CENTERS(s)
%
%   Computes cell centers as the average of the four surrounding corner
%   coordinates, converts them to element space (chi, eta), and builds
%   an extended coordinate array (x_Ext, y_Ext) that includes one layer
%   of ghost cells on each boundary. Ghost cell positions are extrapolated
%   by linear reflection from the boundary face centers.
%
%   Inputs:
%       s - Struct containing at minimum:
%           .mesh.Nchi, .mesh.Neta             - Number of cells in x and y
%           .mesh.x_corner, .mesh.y_corner     - (Nx+1 x Ny+1) corner coordinates
%           .curvilinear_mapping.boundary_type  - String for GO_TO_ELEMENT_SPACE
%
%   Outputs:
%       s - Updated struct with added/modified fields:
%           .mesh.x, .mesh.y         - (Nx x Ny) cell center coordinates
%           .mesh.chi, .mesh.eta     - (Nx x Ny) element-space coordinates
%           .mesh.x_Ext, .mesh.y_Ext - (Nx+2 x Ny+2) extended cell centers
%                                      including ghost cells
%
%   Notes:
%       - Interior cell centers occupy indices (2:end-1, 2:end-1) in the
%         extended arrays.
%       - Boundary face centers (midpoints of adjacent corners) are placed
%         at the extended array edges, then reflected outward to create
%         ghost cell positions for second-order boundary treatment.
%       - Element-space conversion is performed by GO_TO_ELEMENT_SPACE.
%
% Part of: Hypersonics Stability MATLAB Solver - Mesh Module

    %% Compute interior cell centers
    for i = 1:s.mesh.Nchi
        for j = 1:s.mesh.Neta
            s.mesh.x(i,j) = (s.mesh.x_corner(i,j) + s.mesh.x_corner(i+1,j) + s.mesh.x_corner(i,j+1) + s.mesh.x_corner(i+1,j+1)) / 4;
            s.mesh.y(i,j) = (s.mesh.y_corner(i,j) + s.mesh.y_corner(i+1,j) + s.mesh.y_corner(i,j+1) + s.mesh.y_corner(i+1,j+1)) / 4;
        end
    end

    %% Convert to element space
    [s.mesh.chi, s.mesh.eta] = GO_TO_ELEMENT_SPACE(s.mesh.x, s.mesh.y, s);

    %% Populate interior of extended arrays
    s.mesh.x_Ext(2:end-1, 2:end-1) = s.mesh.x;
    s.mesh.y_Ext(2:end-1, 2:end-1) = s.mesh.y;

    %% Boundary face centers (midpoints of adjacent corners)
    s.mesh.x_Ext(2:end-1, 1)   = (s.mesh.x_corner(1:end-1, 1)   + s.mesh.x_corner(2:end, 1))   / 2; % Bottom boundary
    s.mesh.y_Ext(2:end-1, 1)   = (s.mesh.y_corner(1:end-1, 1)   + s.mesh.y_corner(2:end, 1))   / 2;
    s.mesh.x_Ext(2:end-1, end) = (s.mesh.x_corner(1:end-1, end) + s.mesh.x_corner(2:end, end)) / 2; % Top boundary
    s.mesh.y_Ext(2:end-1, end) = (s.mesh.y_corner(1:end-1, end) + s.mesh.y_corner(2:end, end)) / 2;
    s.mesh.x_Ext(1, 2:end-1)   = (s.mesh.x_corner(1, 1:end-1)  + s.mesh.x_corner(1, 2:end))   / 2; % Right boundary
    s.mesh.y_Ext(1, 2:end-1)   = (s.mesh.y_corner(1, 1:end-1)  + s.mesh.y_corner(1, 2:end))   / 2;
    s.mesh.x_Ext(end, 2:end-1) = (s.mesh.x_corner(end, 1:end-1) + s.mesh.x_corner(end, 2:end)) / 2; % Left boundary
    s.mesh.y_Ext(end, 2:end-1) = (s.mesh.y_corner(end, 1:end-1) + s.mesh.y_corner(end, 2:end)) / 2;

    %% Corner points of extended array
    s.mesh.x_Ext(1, 1)     = s.mesh.x_corner(1, 1);       % Bottom right corner
    s.mesh.y_Ext(1, 1)     = s.mesh.y_corner(1, 1);
    s.mesh.x_Ext(1, end)   = s.mesh.x_corner(1, end);     % Top right corner
    s.mesh.y_Ext(1, end)   = s.mesh.y_corner(1, end);
    s.mesh.x_Ext(end, 1)   = s.mesh.x_corner(end, 1);     % Bottom left corner
    s.mesh.y_Ext(end, 1)   = s.mesh.y_corner(end, 1);
    s.mesh.x_Ext(end, end) = s.mesh.x_corner(end, end);   % Top left corner
    s.mesh.y_Ext(end, end) = s.mesh.y_corner(end, end);

    %% Ghost cells (linear extrapolation from boundary face centers)
    s.mesh.x_Ext(:, 1)   = 2*s.mesh.x_Ext(:, 1)   - s.mesh.x_Ext(:, 2);
    s.mesh.x_Ext(:, end) = 2*s.mesh.x_Ext(:, end) - s.mesh.x_Ext(:, end-1);
    s.mesh.x_Ext(1, :)   = 2*s.mesh.x_Ext(1, :)   - s.mesh.x_Ext(2, :);
    s.mesh.x_Ext(end, :) = 2*s.mesh.x_Ext(end, :) - s.mesh.x_Ext(end-1, :);
    s.mesh.y_Ext(:, 1)   = 2*s.mesh.y_Ext(:, 1)   - s.mesh.y_Ext(:, 2);
    s.mesh.y_Ext(:, end) = 2*s.mesh.y_Ext(:, end) - s.mesh.y_Ext(:, end-1);
    s.mesh.y_Ext(1, :)   = 2*s.mesh.y_Ext(1, :)   - s.mesh.y_Ext(2, :);
    s.mesh.y_Ext(end, :) = 2*s.mesh.y_Ext(end, :) - s.mesh.y_Ext(end-1, :);
end
