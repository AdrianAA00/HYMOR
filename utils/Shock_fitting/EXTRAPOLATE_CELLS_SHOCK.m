function [s] = EXTRAPOLATE_CELLS_SHOCK(s)
% EXTRAPOLATE_CELLS_SHOCK - Populate auxiliary ghost cells beyond the shock.
%
%   s = EXTRAPOLATE_CELLS_SHOCK(s)
%
%   Fills the two auxiliary (ghost) cell layers immediately upstream of
%   each shocked cell using extrapolation from the downstream flow field.
%   The extrapolation order is controlled by s.shock.interpolate:
%     "1st" - zeroth-order (constant) copy of the shocked-cell value
%     "2nd" - first-order (linear) extrapolation using two downstream points
%     "3rd" - second-order (quadratic) extrapolation using three downstream points
%
%   This ensures that flux computations near the shock do not encounter
%   discontinuous jumps in the reconstructed stencil.
%
%   Inputs:
%       s (struct) - Solution structure containing at minimum:
%                    .mesh.Nchi, .mesh.Neta   - Grid dimensions
%                    .shock.cell_indices      - (Nx x 1) column index of shocked cell
%                    .shock.interpolate       - Extrapolation order: "1st", "2nd", or "3rd"
%                    .var.rho, .var.rho_u, .var.rho_v,
%                    .var.rho_E, .var.p       - (Nx+2 x Ny+2) conservative/pressure fields
%                    .var.rho_w               - azimuthal momentum (azimuthal modes only)
%                    .chemistry.chemical_equilibrium - (logical) equilibrium chemistry flag
%                    .chemistry.is_chemistry_enabled - (logical) chemistry enabled flag
%                    .var.gamma_star, .var.cv_star   - Thermodynamic property fields (if chemistry)
%
%   Outputs:
%       s (struct) - Solution with auxiliary ghost cells filled by extrapolation.
%
%   Notes:
%       - Ghost cell indices are offset by +1 due to the presence of boundary
%         ghost cells in the s arrays.
%       - Uses sub2ind for vectorized indexing across all streamwise stations.
%
% Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

    %% Store original field values for extrapolation source
    temp_p = s.var.p;
    if ~s.chemistry.chemical_equilibrium && s.chemistry.is_chemistry_enabled
        temp_gamma_star = s.var.gamma_star;
        temp_cv_star    = s.var.cv_star;
    end

    %% Compute linear indices for shocked and neighboring cells
    i_values_s = (1:s.mesh.Nchi)' + 1;              % Row indices (offset for ghost cells)
    j_values_s = s.shock.cell_indices(:, 1);        % Column indices of shocked cells

    % OPTIMIZATION: Manual linear indexing
    Nx = size(s.var.rho, 1);
    idx = struct( ...
        'p2',  i_values_s + (j_values_s + 2) .* Nx, ...  % 2nd auxiliary cell
        'p1',  i_values_s + (j_values_s + 1) .* Nx, ...  % 1st auxiliary cell
        'c',   i_values_s + (j_values_s) .* Nx, ...      % Shocked cell
        'm1',  i_values_s + (j_values_s - 1) .* Nx, ...  % 1 cell downstream
        'm2',  i_values_s + (j_values_s - 2) .* Nx);     % 2 cells downstream

    %% Conservative variables
    for f = PERTURBATION_FIELDS(s)
        s.var.(f) = EXTRAPOLATE_FIELD(s.var.(f), s.shock.interpolate, idx);
    end

    %% Pressure and thermodynamic properties
    if s.shock.interpolate == "1st"
        s.var.p(idx.p1) = temp_p(idx.c);
        s.var.p(idx.p2) = temp_p(idx.c);
        if ~s.chemistry.chemical_equilibrium && s.chemistry.is_chemistry_enabled
            s.var.gamma_star(idx.p1) = temp_gamma_star(idx.c);
            s.var.cv_star(idx.p1)    = temp_cv_star(idx.c);
        end

    elseif s.shock.interpolate == "2nd"
        s.var.p(idx.p1) = 2 * temp_p(idx.c) - temp_p(idx.m1);
        if ~s.chemistry.chemical_equilibrium && s.chemistry.is_chemistry_enabled
            s.var.gamma_star(idx.p1) = 2 * temp_gamma_star(idx.c) - temp_gamma_star(idx.m1);
            s.var.cv_star(idx.p1)    = 2 * temp_cv_star(idx.c)    - temp_cv_star(idx.m1);
            s.var.gamma_star(idx.p2) = 2 * temp_gamma_star(idx.p1) - temp_gamma_star(idx.c);
            s.var.cv_star(idx.p2)    = 2 * temp_cv_star(idx.p1)    - temp_cv_star(idx.c);
        end

    elseif s.shock.interpolate == "3rd"
        s.var.p(idx.p1) = 3 * temp_p(idx.c) - 3 * temp_p(idx.m1) + temp_p(idx.m2);
        if ~s.chemistry.chemical_equilibrium && s.chemistry.is_chemistry_enabled
            s.var.gamma_star(idx.p1) = temp_gamma_star(idx.c);
            s.var.cv_star(idx.p1)    = temp_cv_star(idx.c);
            s.var.gamma_star(idx.p2) = temp_gamma_star(idx.p1);
            s.var.cv_star(idx.p2)    = temp_cv_star(idx.p1);
        end
    end
end

%% ========================================================================
%  Extrapolation of One Field
%  ========================================================================
function f = EXTRAPOLATE_FIELD(f, order, idx)
% EXTRAPOLATE_FIELD  Fill the two auxiliary cells of one field from the
%   values before the update (constant, linear or quadratic in the cell
%   index, as selected by order).
    t = f;
    switch order
        case "1st"
            f(idx.p1) = t(idx.c);
            f(idx.p2) = t(idx.c);
        case "2nd"
            f(idx.p1) = 2 * t(idx.c) - t(idx.m1);
            f(idx.p2) = 2 * t(idx.p1) - t(idx.c);
        case "3rd"
            f(idx.p1) = 3 * t(idx.c) - 3 * t(idx.m1) + t(idx.m2);
            f(idx.p2) = 3 * t(idx.p1) - 3 * t(idx.c) + t(idx.m1);
    end
end
