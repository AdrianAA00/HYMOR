function [dE_dt_inflow, dE_dt_V_inflow] = ENERGY_INFLOW(V, s, w_infty)
% ENERGY_INFLOW  Compute inflow energy flux and its per-unit-volume density.
%
%   [dE_dt_inflow, dE_dt_V_inflow] = ENERGY_INFLOW(V, s, w_infty)
%
%   Evaluates the total perturbation energy flux entering the domain through
%   the freestream boundary. The flux operator D is constructed from the
%   freestream coupling matrices (R_ and M_infty_) and applied to the global
%   perturbation vector V. Individual acoustic, kinetic, and entropic
%   contributions are computed and reported. The result is also normalised
%   by the total post-shock volume to yield an energy density metric.
%
%   Inputs:
%       V        - (double column vector) Global perturbation state vector.
%       s - (struct) Solution structure with grid information (Nx, Ny,
%                  volume, flow_cells) and freestream parameters needed by
%                  CONSTRUCT_R_ and CONSTRUCT_M_INFTY_.
%       w_infty  - (double vector) Freestream modal frequency/weight vector.
%
%   Outputs:
%       dE_dt_inflow   - (double scalar) Total perturbation energy flux
%                        through the inflow boundary (V' * D * V).
%       dE_dt_V_inflow - (double scalar) Energy flux normalised by the total
%                        post-shock volume (energy density rate).
%
%   Notes:
%       - The operator D = R_' * M_infty_ * R_ is assembled four times with
%         different norm selectors to isolate the total, acoustic (norms [1,0,0,0]),
%         kinetic (norms [0,1,1,0]), and entropic (norms [0,0,0,1]) contributions.
%       - Diagnostic messages are printed to the console showing the fractional
%         contribution of each energy component.
%       - Requires helper functions: CONSTRUCT_R_, CONSTRUCT_M_INFTY_.
%
% Part of: Hypersonics Stability MATLAB Solver - Energy Budgets Module

    %% Construct inflow flux operators
    R_ = CONSTRUCT_R_(s, w_infty);
    T = 1;  % Flux-of-energy mode
    scaling_non_temporal = false;

    % Total energy flux operator
    norms_total = [1, 1, 1, 1];
    M_infty_total = CONSTRUCT_M_INFTY_(s, norms_total, T, w_infty, scaling_non_temporal);
    D_ = R_' * M_infty_total * R_;

    % Acoustic-only operator
    norms_p = [1, 0, 0, 0];
    M_infty_p = CONSTRUCT_M_INFTY_(s, norms_p, T, w_infty, scaling_non_temporal);
    D_p = R_' * M_infty_p * R_;

    % Kinetic-only operator
    norms_u = [0, 1, 1, 0];
    M_infty_u = CONSTRUCT_M_INFTY_(s, norms_u, T, w_infty, scaling_non_temporal);
    D_u = R_' * M_infty_u * R_;

    % Entropic-only operator
    norms_s = [0, 0, 0, 1];
    M_infty_s = CONSTRUCT_M_INFTY_(s, norms_s, T, w_infty, scaling_non_temporal);
    D_s = R_' * M_infty_s * R_;

    %% Evaluate inflow energy fluxes
    dE_dt_inflow   = V' * D_  * V;
    dE_dt_inflow_p = V' * D_p * V;
    dE_dt_inflow_s = V' * D_s * V;
    dE_dt_inflow_u = V' * D_u * V;

    disp("Freestream acoustic: " + abs(dE_dt_inflow_p) / abs(dE_dt_inflow))
    disp("Freestream entropic: " + abs(dE_dt_inflow_s) / abs(dE_dt_inflow))
    disp("Freestream kinetic: "  + abs(dE_dt_inflow_u) / abs(dE_dt_inflow))

    %% Compute total post-shock volume (vectorized)
    volume     = s.mesh.volume;
    flow_cells = s.shock.flow_cells;
    total_post_shock_volume = sum(volume .* flow_cells, "all");

    % Energy density: inflow flux per unit volume
    dE_dt_V_inflow = dE_dt_inflow / total_post_shock_volume;
end
