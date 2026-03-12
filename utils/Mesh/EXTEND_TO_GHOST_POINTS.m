function v_ext = EXTEND_TO_GHOST_POINTS(v, s)
% EXTEND_TO_GHOST_POINTS - Extrapolate a field to include ghost-cell values.
%
%   v_ext = EXTEND_TO_GHOST_POINTS(v, s)
%
%   Takes an (Nx x Ny) interior field and produces an (Nx+2 x Ny+2)
%   extended field by adding one layer of ghost cells on each side.
%   Ghost values are set by linear extrapolation from the two nearest
%   interior values. If a shock is present, the ghost cells immediately
%   upstream of the shock are overwritten with values extrapolated from
%   the post-shock side to prevent spurious flux computation across the
%   discontinuity.
%
%   Inputs:
%       v         (Nx x Ny double) - Interior field values
%       s  (struct)         - Solution structure containing:
%           .mesh.Nchi, .mesh.Neta  - Number of interior cells
%           .shock.enabled          - (logical) Shock presence flag
%           .shock.cell_indices     - (Nx x 1+) column indices of shocked cells
%
%   Outputs:
%       v_ext  ((Nx+2) x (Ny+2) double) - Extended field with ghost cells
%
%   Notes:
%       - Ghost cells use second-order linear extrapolation:
%             v_ghost = 2*v_boundary - v_interior
%       - Corner ghost points are extrapolated from the already-computed
%         ghost rows/columns.
%       - For shocked meshes, the ghost cell at the shock interface is
%         replaced by extrapolation from the post-shock interior to
%         avoid mixing pre- and post-shock states.
%
% Part of: Hypersonics Stability MATLAB Solver - Mesh Module

    Nx = s.mesh.Nchi;
    Ny = s.mesh.Neta;

    %% Initialize extended matrix (Nx+2) x (Ny+2)
    v_ext = zeros(Nx + 2, Ny + 2);

    %% Copy interior values
    v_ext(2:end-1, 2:end-1) = v;

    %% Extrapolate ghost edges (linear extrapolation)
    v_ext(2:end-1, 1)   = 2 * v(:, 1)   - v(:, 2);         % Left column
    v_ext(2:end-1, end) = 2 * v(:, end) - v(:, end-1);      % Right column
    v_ext(1, 2:end-1)   = 2 * v(1, :)   - v(2, :);          % Top row
    v_ext(end, 2:end-1) = 2 * v(end, :) - v(end-1, :);      % Bottom row

    %% Extrapolate corner ghost points
    v_ext(1, 1)     = 2 * v_ext(1, 2)     - v_ext(1, 3);          % Top-left
    v_ext(1, end)   = 2 * v_ext(1, end-1) - v_ext(1, end-2);      % Top-right
    v_ext(end, 1)   = 2 * v_ext(end, 2)   - v_ext(end, 3);        % Bottom-left
    v_ext(end, end) = 2 * v_ext(end, end-1) - v_ext(end, end-2);  % Bottom-right

    %% Shock ghost-point override
    if s.shock.enabled
        temp = v_ext;

        % Compute linear indices for shocked cells and their neighbours
        i_values_s = (1:Nx)' + 1;                          % Row indices (+1 for ghost offset)
        j_values_s = s.shock.cell_indices(:, 1);   % Column indices of shocked cells

        idx_0  = sub2ind(size(v_ext), i_values_s, j_values_s + 1);  % Shocked cell
        idx_m1 = sub2ind(size(v_ext), i_values_s, j_values_s);      % One cell post-shock
        idx_m2 = sub2ind(size(v_ext), i_values_s, j_values_s - 1);  % Two cells post-shock

        % Extrapolate from post-shock side into the shocked cell
        v_ext(idx_0) = 2 * temp(idx_m1) - temp(idx_m2);
    end
end
