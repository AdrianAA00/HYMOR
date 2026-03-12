function s = FILTERING(s)
% FILTERING - Apply numerical dissipation filtering to conserved variables.
%
%   Blends the computed flux with a 5-point Laplacian smoothing stencil to
%   provide numerical dissipation. Each conserved variable (rho, rho*u,
%   rho*v, rho*E) is filtered independently using its own numerical
%   viscosity coefficient (mu_numerical_*). The blending formula is:
%
%       flux_filtered = flux * (1 - mu) + mu * (stencil_avg / 8)
%
%   where the stencil average uses the center cell weighted by 4 and each
%   of the four cardinal neighbours weighted by 1.
%
% Syntax:
%   s = FILTERING(s)
%
% Inputs:
%   s - Structure containing at minimum:
%                .var.rho, .var.rho_u, .var.rho_v, .var.rho_E - Conserved variable fields.
%                .flux.rho, .flux.rho_u,
%                .flux.rho_v, .flux.rho_E      - Pre-computed flux fields.
%                .numerical_dissipation.mu_rho  - Filtering coefficient for density.
%                .numerical_dissipation.mu_rho_u - Filtering coefficient for x-momentum.
%                .numerical_dissipation.mu_rho_v - Filtering coefficient for y-momentum.
%                .numerical_dissipation.mu_rho_E - Filtering coefficient for energy.
%
% Outputs:
%   s - Updated structure with filtered flux fields.
%
% Notes:
%   - The stencil (4*center + N + S + E + W) / 8 corresponds to a standard
%     5-point averaging operator.
%   - Setting mu_numerical_* = 0 disables filtering for that variable.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% Filter density flux
    rho_temp = 4 * s.var.rho(2:end-1,2:end-1);
    rho_temp = rho_temp + s.var.rho(1:end-2,2:end-1);
    rho_temp = rho_temp + s.var.rho(2:end-1,1:end-2);
    rho_temp = rho_temp + s.var.rho(3:end,2:end-1);
    rho_temp = rho_temp + s.var.rho(2:end-1,3:end);
    s.flux.rho = s.flux.rho * (1 - s.numerical_dissipation.mu_rho) + ...
                        s.numerical_dissipation.mu_rho * rho_temp / 8;

    %% Filter x-momentum flux
    rho_u_temp = 4 * s.var.rho_u(2:end-1,2:end-1);
    rho_u_temp = rho_u_temp + s.var.rho_u(1:end-2,2:end-1);
    rho_u_temp = rho_u_temp + s.var.rho_u(2:end-1,1:end-2);
    rho_u_temp = rho_u_temp + s.var.rho_u(3:end,2:end-1);
    rho_u_temp = rho_u_temp + s.var.rho_u(2:end-1,3:end);
    s.flux.rho_u = s.flux.rho_u * (1 - s.numerical_dissipation.mu_rho_u) + ...
                          s.numerical_dissipation.mu_rho_u * rho_u_temp / 8;

    %% Filter y-momentum flux
    rho_v_temp = 4 * s.var.rho_v(2:end-1,2:end-1);
    rho_v_temp = rho_v_temp + s.var.rho_v(1:end-2,2:end-1);
    rho_v_temp = rho_v_temp + s.var.rho_v(2:end-1,1:end-2);
    rho_v_temp = rho_v_temp + s.var.rho_v(3:end,2:end-1);
    rho_v_temp = rho_v_temp + s.var.rho_v(2:end-1,3:end);
    s.flux.rho_v = s.flux.rho_v * (1 - s.numerical_dissipation.mu_rho_v) + ...
                          s.numerical_dissipation.mu_rho_v * rho_v_temp / 8;

    %% Filter total energy flux
    rho_E_temp = 4 * s.var.rho_E(2:end-1,2:end-1);
    rho_E_temp = rho_E_temp + s.var.rho_E(1:end-2,2:end-1);
    rho_E_temp = rho_E_temp + s.var.rho_E(2:end-1,1:end-2);
    rho_E_temp = rho_E_temp + s.var.rho_E(3:end,2:end-1);
    rho_E_temp = rho_E_temp + s.var.rho_E(2:end-1,3:end);
    s.flux.rho_E = s.flux.rho_E * (1 - s.numerical_dissipation.mu_rho_E) + ...
                          s.numerical_dissipation.mu_rho_E * rho_E_temp / 8;
end
