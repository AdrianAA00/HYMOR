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
    temp_rho   = s.var.rho;
    temp_rho_u = s.var.rho_u;
    temp_rho_v = s.var.rho_v;
    temp_rho_E = s.var.rho_E;
    temp_p     = s.var.p;
    if ~s.chemistry.chemical_equilibrium && s.chemistry.is_chemistry_enabled
        temp_gamma_star = s.var.gamma_star;
        temp_cv_star    = s.var.cv_star;
    end

    %% Compute linear indices for shocked and neighboring cells
    i_values_s = (1:s.mesh.Nchi)' + 1;              % Row indices (offset for ghost cells)
    j_values_s = s.shock.cell_indices(:, 1);        % Column indices of shocked cells

    % OPTIMIZATION: Manual linear indexing 
    Nx = size(s.var.rho, 1);
    idx_2  = i_values_s + (j_values_s + 2) .* Nx;  % 2nd auxiliary cell
    idx_1  = i_values_s + (j_values_s + 1) .* Nx;  % 1st auxiliary cell
    idx_0  = i_values_s + (j_values_s) .* Nx;      % Shocked cell
    idx_m1 = i_values_s + (j_values_s - 1) .* Nx;  % 1 cell downstream
    idx_m2 = i_values_s + (j_values_s - 2) .* Nx;  % 2 cells downstream

    %% Apply extrapolation based on selected order
    if s.shock.interpolate == "1st"
        % -- Zeroth-order (constant) extrapolation --
        s.var.rho(idx_1)   = temp_rho(idx_0);
        s.var.rho_u(idx_1) = temp_rho_u(idx_0);
        s.var.rho_v(idx_1) = temp_rho_v(idx_0);
        s.var.rho_E(idx_1) = temp_rho_E(idx_0);
        s.var.p(idx_1)     = temp_p(idx_0);
        if ~s.chemistry.chemical_equilibrium && s.chemistry.is_chemistry_enabled
            s.var.gamma_star(idx_1) = temp_gamma_star(idx_0);
            s.var.cv_star(idx_1)    = temp_cv_star(idx_0);
        end

        s.var.rho(idx_2)   = temp_rho(idx_0);
        s.var.rho_u(idx_2) = temp_rho_u(idx_0);
        s.var.rho_v(idx_2) = temp_rho_v(idx_0);
        s.var.rho_E(idx_2) = temp_rho_E(idx_0);
        s.var.p(idx_2)     = temp_p(idx_0);

    elseif s.shock.interpolate == "2nd"
        % -- First-order (linear) extrapolation --
        s.var.rho(idx_1)   = 2 * temp_rho(idx_0)   - temp_rho(idx_m1);
        s.var.rho_u(idx_1) = 2 * temp_rho_u(idx_0) - temp_rho_u(idx_m1);
        s.var.rho_v(idx_1) = 2 * temp_rho_v(idx_0) - temp_rho_v(idx_m1);
        s.var.rho_E(idx_1) = 2 * temp_rho_E(idx_0) - temp_rho_E(idx_m1);
        s.var.p(idx_1)     = 2 * temp_p(idx_0)     - temp_p(idx_m1);
        if ~s.chemistry.chemical_equilibrium && s.chemistry.is_chemistry_enabled
            s.var.gamma_star(idx_1) = 2 * temp_gamma_star(idx_0) - temp_gamma_star(idx_m1);
            s.var.cv_star(idx_1)    = 2 * temp_cv_star(idx_0)    - temp_cv_star(idx_m1);
        end

        s.var.rho(idx_2)   = 2 * temp_rho(idx_1)   - temp_rho(idx_0);
        s.var.rho_u(idx_2) = 2 * temp_rho_u(idx_1) - temp_rho_u(idx_0);
        s.var.rho_v(idx_2) = 2 * temp_rho_v(idx_1) - temp_rho_v(idx_0);
        s.var.rho_E(idx_2) = 2 * temp_rho_E(idx_1) - temp_rho_E(idx_0);
        if ~s.chemistry.chemical_equilibrium && s.chemistry.is_chemistry_enabled
            s.var.gamma_star(idx_2) = 2 * temp_gamma_star(idx_1) - temp_gamma_star(idx_0);
            s.var.cv_star(idx_2)    = 2 * temp_cv_star(idx_1)    - temp_cv_star(idx_0);
        end

    elseif s.shock.interpolate == "3rd"
        % -- Second-order (quadratic) extrapolation --
        s.var.rho(idx_1)   = 3 * temp_rho(idx_0)   - 3 * temp_rho(idx_m1)   + temp_rho(idx_m2);
        s.var.rho_u(idx_1) = 3 * temp_rho_u(idx_0) - 3 * temp_rho_u(idx_m1) + temp_rho_u(idx_m2);
        s.var.rho_v(idx_1) = 3 * temp_rho_v(idx_0) - 3 * temp_rho_v(idx_m1) + temp_rho_v(idx_m2);
        s.var.rho_E(idx_1) = 3 * temp_rho_E(idx_0) - 3 * temp_rho_E(idx_m1) + temp_rho_E(idx_m2);
        s.var.p(idx_1)     = 3 * temp_p(idx_0)     - 3 * temp_p(idx_m1)     + temp_p(idx_m2);
        if ~s.chemistry.chemical_equilibrium && s.chemistry.is_chemistry_enabled
            s.var.gamma_star(idx_1) = temp_gamma_star(idx_0);
            s.var.cv_star(idx_1)    = temp_cv_star(idx_0);
        end

        s.var.rho(idx_2)   = 3 * temp_rho(idx_1)   - 3 * temp_rho(idx_0)   + temp_rho(idx_m1);
        s.var.rho_u(idx_2) = 3 * temp_rho_u(idx_1) - 3 * temp_rho_u(idx_0) + temp_rho_u(idx_m1);
        s.var.rho_v(idx_2) = 3 * temp_rho_v(idx_1) - 3 * temp_rho_v(idx_0) + temp_rho_v(idx_m1);
        s.var.rho_E(idx_2) = 3 * temp_rho_E(idx_1) - 3 * temp_rho_E(idx_0) + temp_rho_E(idx_m1);
        if ~s.chemistry.chemical_equilibrium && s.chemistry.is_chemistry_enabled
            s.var.gamma_star(idx_2) = temp_gamma_star(idx_1);
            s.var.cv_star(idx_2)    = temp_cv_star(idx_1);
        end
    end
end
