function s = UPDATE_SHOCK_BC(s, chemistry)
% UPDATE_SHOCK_BC - Apply shock boundary conditions via Rankine-Hugoniot jump relations.
%
%   s = UPDATE_SHOCK_BC(s, chemistry)
%
%   When shock fitting is active, this function performs the full shock
%   boundary-condition update sequence:
%     1. Fit/resample the shock points (LEAST_SQUARES_SHOCK_POINTS).
%     2. Compute the local shock angle beta (COMPUTE_BETA).
%     3. Enforce symmetry at the last station (beta = pi/2).
%     4. Update shocked-cell values from Rankine-Hugoniot relations
%        (UPDATE_SHOCK_JUMP_PROPERTIES).
%     5. Extrapolate into ghost cells (EXTRAPOLATE_CELLS_SHOCK).
%
%   Inputs:
%       s  (struct) - Solution structure with .shock.enabled flag and flow fields.
%       chemistry (struct) - Chemistry model structure.
%
%   Outputs:
%       s  (struct) - Updated s with shock boundary conditions applied.
%
%   See also: LEAST_SQUARES_SHOCK_POINTS, COMPUTE_BETA,
%             UPDATE_SHOCK_JUMP_PROPERTIES, EXTRAPOLATE_CELLS_SHOCK
%
% Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

    if s.shock.enabled
        s = LEAST_SQUARES_SHOCK_POINTS(s);
        s = COMPUTE_BETA(s);
        if strcmp(s.boundary_conditions.boundary_chi0.name, 'symmetry')
            s.shock.beta(1, 1) = pi / 2;  % Symmetry condition
        end
        s = UPDATE_SHOCK_JUMP_PROPERTIES(s, chemistry);
        s = EXTRAPOLATE_CELLS_SHOCK(s);
    end
end
