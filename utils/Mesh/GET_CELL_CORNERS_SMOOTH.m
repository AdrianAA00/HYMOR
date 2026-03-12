function s = GET_CELL_CORNERS_SMOOTH(s)
% GET_CELL_CORNERS_SMOOTH  Smooth mesh corners via weighted averaging.
%
%   s = GET_CELL_CORNERS_SMOOTH(s)
%
%   Recomputes interior cell corner positions as a weighted average of the
%   four surrounding extended cell centers. An exponential weight function
%   preserves the original wall boundary coordinates (j = 1 row) while
%   increasingly smoothing corners farther from the wall.
%
%   Inputs:
%       s - Struct containing at minimum:
%           .mesh.Nchi, .mesh.Neta       - Number of cells in x and y
%           .mesh.x_corner, .mesh.y_corner - (Nx+1 x Ny+1) corner coordinates
%           .mesh.x_Ext, .mesh.y_Ext      - (Nx+2 x Ny+2) extended cell centers
%
%   Outputs:
%       s - Updated struct with smoothed .mesh.x_corner and .mesh.y_corner
%           (wall row j=1 is unchanged).
%
%   Notes:
%       - The weight function is w = 1 - exp(-10*j/Ny), which is nearly
%         zero for j close to the wall and approaches 1 away from it.
%       - This function is typically called iteratively (e.g., 10 times)
%         in combination with GET_EXTENDED_CELL_CENTERS to progressively
%         smooth meshes with curvature discontinuities (MSL, blunt_cone).
%
% Part of: Hypersonics Stability MATLAB Solver - Mesh Module

    %% Smooth interior corners via weighted averaging
    for i = 1:s.mesh.Nchi+1
        for j = 2:s.mesh.Neta+1
            weight = 1 - exp(-10*j/s.mesh.Neta);
            s.mesh.x_corner(i,j) = (s.mesh.x_Ext(i,j) + s.mesh.x_Ext(i+1,j) + s.mesh.x_Ext(i,j+1) + s.mesh.x_Ext(i+1,j+1))/4 * weight ...
                                   + s.mesh.x_corner(i,j) * (1-weight);
            s.mesh.y_corner(i,j) = (s.mesh.y_Ext(i,j) + s.mesh.y_Ext(i+1,j) + s.mesh.y_Ext(i,j+1) + s.mesh.y_Ext(i+1,j+1))/4 * weight ...
                                   + s.mesh.y_corner(i,j) * (1-weight);
        end
    end
end
