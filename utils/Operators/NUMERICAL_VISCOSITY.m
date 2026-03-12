function s = NUMERICAL_VISCOSITY(s)
% NUMERICAL_VISCOSITY  Add artificial numerical viscosity to all conserved variable fluxes.
%
%   s = NUMERICAL_VISCOSITY(s)
%
%   Applies a Laplacian-type artificial dissipation to the density,
%   momentum, and total energy flux fields. The dissipation is computed
%   using a discrete Laplacian stencil weighted by face lengths and
%   cell volumes, scaled by field-specific numerical viscosity coefficients.
%
%   Inputs:
%       s - Solution struct containing:
%                  .var.rho, .var.rho_u, .var.rho_v, .var.rho_E - Conserved variable fields
%                  .bt_length, .lr_length       - Bottom-top and left-right face lengths
%                  .mesh.volume                 - Cell volumes
%                  .numerical_dissipation.mu_rho   - Numerical viscosity coeff. for density
%                  .numerical_dissipation.mu_rho_u - Numerical viscosity coeff. for x-momentum
%                  .numerical_dissipation.mu_rho_v - Numerical viscosity coeff. for y-momentum
%                  .numerical_dissipation.mu_rho_E - Numerical viscosity coeff. for energy
%
%   Outputs:
%       s - Updated s struct with numerical dissipation added
%                  to .flux.rho, .flux.rho_u, .flux.rho_v, .flux.rho_E
%
%   Notes:
%       - Uses a 5-point Laplacian stencil on the structured grid.
%       - Each conserved variable has its own viscosity coefficient.
%       - Ghost cell values participate in the stencil but only interior
%         fluxes are modified.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% Density numerical viscosity
    rho_temp = -s.var.rho(2:end-1, 2:end-1) .* (s.bt_length(:, 1:end-1) + s.bt_length(:, 2:end) + s.lr_length(1:end-1, :) + s.lr_length(2:end, :));
    rho_temp = rho_temp + s.var.rho(1:end-2, 2:end-1) .* s.lr_length(1:end-1, :);
    rho_temp = rho_temp + s.var.rho(2:end-1, 1:end-2) .* s.bt_length(:, 1:end-1);
    rho_temp = rho_temp + s.var.rho(3:end, 2:end-1)   .* s.lr_length(2:end, :);
    rho_temp = rho_temp + s.var.rho(2:end-1, 3:end)   .* s.bt_length(:, 2:end);
    s.flux.rho = s.flux.rho + s.numerical_dissipation.mu_rho * rho_temp ./ s.mesh.volume;

    %% x-momentum numerical viscosity
    rho_u_temp = -s.var.rho_u(2:end-1, 2:end-1) .* (s.bt_length(:, 1:end-1) + s.bt_length(:, 2:end) + s.lr_length(1:end-1, :) + s.lr_length(2:end, :));
    rho_u_temp = rho_u_temp + s.var.rho_u(1:end-2, 2:end-1) .* s.lr_length(1:end-1, :);
    rho_u_temp = rho_u_temp + s.var.rho_u(2:end-1, 1:end-2) .* s.bt_length(:, 1:end-1);
    rho_u_temp = rho_u_temp + s.var.rho_u(3:end, 2:end-1)   .* s.lr_length(2:end, :);
    rho_u_temp = rho_u_temp + s.var.rho_u(2:end-1, 3:end)   .* s.bt_length(:, 2:end);
    s.flux.rho_u = s.flux.rho_u + s.numerical_dissipation.mu_rho_u * rho_u_temp ./ s.mesh.volume;

    %% y-momentum numerical viscosity
    rho_v_temp = -s.var.rho_v(2:end-1, 2:end-1) .* (s.bt_length(:, 1:end-1) + s.bt_length(:, 2:end) + s.lr_length(1:end-1, :) + s.lr_length(2:end, :));
    rho_v_temp = rho_v_temp + s.var.rho_v(1:end-2, 2:end-1) .* s.lr_length(1:end-1, :);
    rho_v_temp = rho_v_temp + s.var.rho_v(2:end-1, 1:end-2) .* s.bt_length(:, 1:end-1);
    rho_v_temp = rho_v_temp + s.var.rho_v(3:end, 2:end-1)   .* s.lr_length(2:end, :);
    rho_v_temp = rho_v_temp + s.var.rho_v(2:end-1, 3:end)   .* s.bt_length(:, 2:end);
    s.flux.rho_v = s.flux.rho_v + s.numerical_dissipation.mu_rho_v * rho_v_temp ./ s.mesh.volume;

    %% Total energy numerical viscosity
    rho_E_temp = -s.var.rho_E(2:end-1, 2:end-1) .* (s.bt_length(:, 1:end-1) + s.bt_length(:, 2:end) + s.lr_length(1:end-1, :) + s.lr_length(2:end, :));
    rho_E_temp = rho_E_temp + s.var.rho_E(1:end-2, 2:end-1) .* s.lr_length(1:end-1, :);
    rho_E_temp = rho_E_temp + s.var.rho_E(2:end-1, 1:end-2) .* s.bt_length(:, 1:end-1);
    rho_E_temp = rho_E_temp + s.var.rho_E(3:end, 2:end-1)   .* s.lr_length(2:end, :);
    rho_E_temp = rho_E_temp + s.var.rho_E(2:end-1, 3:end)   .* s.bt_length(:, 2:end);
    s.flux.rho_E = s.flux.rho_E + s.numerical_dissipation.mu_rho_E * rho_E_temp ./ s.mesh.volume;
end
