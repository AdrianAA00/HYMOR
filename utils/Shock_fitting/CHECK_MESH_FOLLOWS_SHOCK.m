function s = CHECK_MESH_FOLLOWS_SHOCK(s, chemistry)
% CHECK_MESH_FOLLOWS_SHOCK - Verify that mesh rows are aligned with the shock.
%
%   s = CHECK_MESH_FOLLOWS_SHOCK(s, chemistry)
%
%   Checks whether the shocked-cell index is uniform across all streamwise
%   stations. If any station has a different shocked-cell index than the
%   first, the mesh is considered misaligned with the shock and the
%   s is restarted on a new mesh via RESTART_SOLUTION.
%
%   Inputs:
%       s  (struct) - Solution structure containing at minimum:
%                            .mesh.Nchi             - Number of streamwise cells
%                            .shock.cell_indices    - (Nx x 1) indices of shocked cells
%       chemistry (struct) - Chemistry/thermodynamic model structure passed
%                            to RESTART_SOLUTION if remeshing is required.
%
%   Outputs:
%       s  (struct) - Unchanged if mesh is aligned; otherwise the
%                            restarted s on a corrected mesh.
%
%   Notes:
%       - A uniform shocked_cell_indices column is required for proper
%         stability analysis on a shock-aligned grid.
%       - When remeshing is triggered, disturbances are disabled (set to false).
%
%   See also: RESTART_SOLUTION
%
% Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module
    if s.shock.enabled
        for i = 1:s.mesh.Nchi
            if s.shock.cell_indices(i, 1) ~= s.shock.cell_indices(1, 1)
                fprintf("\n Mesh not aligned with shock. Remesh for proper stability analysis. \n")
                solution_old = s;
                disturbances = false;
                s = RESTART_SOLUTION(s, solution_old, chemistry, disturbances);
                return
            end
        end
        fprintf("\n Mesh aligned, not remeshing. \n")
    end

end
