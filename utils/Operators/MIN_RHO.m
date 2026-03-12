function [s] = MIN_RHO(s)
% MIN_RHO  Enforce minimum density threshold on the s field.
%
%   [s] = MIN_RHO(s)
%
%   Clamps the density field to a prescribed minimum value (s.numerical_dissipation.rho_min)
%   within the interior cells. If any interior density values fall below the
%   minimum threshold, they are replaced and a warning message is printed.
%
%   Inputs:
%       s - Solution struct containing:
%                  .var.rho - Density field (2D array with ghost cells)
%                  .numerical_dissipation.rho_min - Minimum allowable density value (scalar)
%
%   Outputs:
%       s - Updated s struct with density limited to rho_min
%
%   Notes:
%       - Only interior cells (2:end-1, 2:end-1) are modified; ghost cells
%         are left unchanged.
%       - A diagnostic message is printed when the limiter activates.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% Apply minimum density threshold
    temp = max(s.var.rho, s.numerical_dissipation.rho_min);

    if (~isequal(s.var.rho(2:end-1, 2:end-1), temp(2:end-1, 2:end-1)))
        s.var.rho(2:end-1, 2:end-1) = temp(2:end-1, 2:end-1);
        fprintf('\n')
        fprintf(" Rho limited activated ")
        fprintf('\n')
    end
end
