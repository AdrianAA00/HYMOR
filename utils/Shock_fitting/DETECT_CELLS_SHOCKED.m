function s = DETECT_CELLS_SHOCKED(s)
% DETECT_CELLS_SHOCKED - Identify grid cells closest to the shock and build cell masks.
%
%   s = DETECT_CELLS_SHOCKED(s)
%
%   For each streamwise station, computes the Euclidean distance from the
%   shock point to every grid node in the normal direction, selects the
%   closest cell as the shocked cell, and constructs three logical masks:
%     - shocked_cells : marks exactly the shocked cell at each station
%     - flow_cells    : marks all cells downstream of the shock (body side)
%     - flow_cells_E  : marks flow cells plus the shocked cell itself
%
%   Inputs:
%       s (struct) - Solution structure containing at minimum:
%                    .mesh.Nchi, .mesh.Neta  - Grid dimensions
%                    .shock.points_x/y      - (Nx x 1) shock point coordinates
%                    .mesh.x, .mesh.y       - (Nx x Ny) grid coordinates
%                    .shock.cells            - (Nx x Ny) preallocated mask
%                    .shock.flow_cells       - (Nx x Ny) preallocated mask
%                    .shock.flow_cells_E     - (Nx x Ny) preallocated mask
%
%   Outputs:
%       s (struct) - Updated s with:
%                    .shock.cell_indices     - (Nx x 1) column index of shocked cell
%                    .shock.cells            - (Nx x Ny) binary mask (1 at shock)
%                    .shock.flow_cells       - (Nx x Ny) binary mask (1 for body-side cells)
%                    .shock.flow_cells_E     - (Nx x Ny) binary mask (flow_cells + shocked cell)
%
%   Notes:
%       - Distance is computed in a fully vectorized manner using broadcasting.
%
% Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

    %% Compute distances from shock points to all grid nodes
    distance = sqrt((s.shock.points_x .* ones(1, s.mesh.Neta) - s.mesh.x).^2 ...
                  + (s.shock.points_y .* ones(1, s.mesh.Neta) - s.mesh.y).^2);
    [~, s.shock.cell_indices] = min(distance, [], 2);

    %% Build cell masks
    s.shock.cells = s.shock.cells * 0;
    s.shock.flow_cells    = s.shock.flow_cells * 0;
    s.shock.flow_cells_E  = s.shock.flow_cells_E * 0;

    for i = 1:s.mesh.Nchi
        s.shock.cells(i, s.shock.cell_indices(i, 1)) = 1;
        s.shock.flow_cells(i, 1:s.shock.cell_indices(i, 1) - 1) = 1;
        s.shock.flow_cells_E(i, 1:s.shock.cell_indices(i, 1)) = 1;
    end
end
