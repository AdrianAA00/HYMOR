function [s] = REMOVE_AUX_SHOCK_CELLS(s)
% REMOVE_AUX_SHOCK_CELLS - Reset auxiliary upstream cells to freestream values.
%
%   s = REMOVE_AUX_SHOCK_CELLS(s)
%
%   Replaces all field values in cells upstream of the shocked cell
%   (including the shocked cell itself) with uniform freestream
%   conditions. This cleans up auxiliary/ghost cells that were
%   populated during the shock-fitting extrapolation step.
%
%   Inputs:
%       s (struct) - Solution structure containing at minimum:
%                           .shock.enabled         - (logical) shock-fitting flag
%                           .mesh.Nchi             - Number of streamwise cells
%                           .shock.cell_indices    - (Nx x 1) shocked-cell indices
%                           .freestream.rho_0, .freestream.rho_u_0, .freestream.rho_v_0,
%                           .freestream.rho_E_0, .freestream.p_0
%                                                  - Scalar freestream conservative values
%                           .freestream.gamma_star
%                           .freestream.T          - Freestream thermodynamic properties
%
%   Outputs:
%       s (struct) - Solution with upstream cells set to freestream values.
%
% Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

    if s.shock.enabled
        for i = 1:s.mesh.Nchi
            sc_i = s.shock.cell_indices(i, 1);
            s.var.rho(i + 1, sc_i + 1:end)        = s.freestream.rho_0;
            s.var.rho_u(i + 1, sc_i + 1:end)      = s.freestream.rho_u_0;
            s.var.rho_v(i + 1, sc_i + 1:end)      = s.freestream.rho_v_0;
            s.var.rho_E(i + 1, sc_i + 1:end)      = s.freestream.rho_E_0;
            s.var.p(i + 1, sc_i + 1:end)          = s.freestream.p_0;
            s.var.gamma_star(i + 1, sc_i + 1:end) = s.freestream.gamma_star;
            s.var.cv_star(i + 1, sc_i + 1:end)    = 1; % Non-dimensional cv_star is 1 for freestream
            s.var.T(i + 1, sc_i + 1:end)          = s.freestream.T;
        end
    end
end
