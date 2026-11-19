function s = UPDATE_SHOCK_BC(s, chemistry)
% UPDATE_SHOCK_BC - Apply shock boundary conditions via Rankine-Hugoniot jump relations.
%
%   s = UPDATE_SHOCK_BC(s, chemistry)
%
%   When shock fitting is active, this function performs the full shock
%   boundary-condition update sequence:
%     1. Fit/resample the shock points (LEAST_SQUARES_SHOCK_POINTS).
%     2. Compute the local shock angle beta (COMPUTE_BETA).
%     3. Enforce symmetry at the last station (beta = pi/2); an odd
%        azimuthal mode keeps its perturbation of the spline angle.
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
            if isfield(s, 'azimuthal') && s.azimuthal.odd
                % Odd azimuthal modes tilt the shock on the axis: keep the
                % perturbation of the spline angle about the symmetric shock
                s.shock.beta(1, 1) = pi / 2 + s.shock.beta(1, 1) ...
                    - s.azimuthal.base.shock_beta_spline_1;
            else
                s.shock.beta(1, 1) = pi / 2;  % Symmetry condition
            end
        end
        s = UPDATE_SHOCK_JUMP_PROPERTIES(s, chemistry);
        s = EXTRAPOLATE_CELLS_SHOCK(s);
    end
end
