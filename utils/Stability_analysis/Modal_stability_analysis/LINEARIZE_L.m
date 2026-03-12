function [s, L] = LINEARIZE_L(s, chemistry)
% LINEARIZE_L  Construct the linearized Jacobian matrix of the flow operator.
%
%   [s, L] = LINEARIZE_L(s, chemistry) builds the sparse
%   Jacobian matrix L that linearizes the nonlinear dynamics operator about
%   the current base flow. The linearization is performed via numerical
%   finite differences using LINEARIZE_EQUATIONS_NO_DISCONTINUITY, and the
%   result is validated against the nonlinear operator using
%   CHECK_LINEARIZATION_L.
%
%   Inputs:
%       s  - Solution structure containing the base flow, grid
%                   parameters, and analysis flags (stability_analysis,
%                   transient_growth_analysis).
%       chemistry - Chemistry model structure for thermodynamic evaluations.
%
%   Outputs:
%       s - Updated s structure (linearize flag toggled during
%                  computation, returned with linearize = false).
%       L        - Sparse Jacobian matrix. Returns 0 if neither stability
%                  analysis nor transient growth analysis is requested.
%
%   Notes:
%       - The mesh is verified to follow the shock before linearization via
%         CHECK_MESH_FOLLOWS_SHOCK.
%       - Chemistry equilibrium state is updated before constructing L.
%       - The s.linearize flag is set to true during construction
%         and restored to false upon return.
%
% Part of: Hypersonics Stability MATLAB Solver - Stability Analysis / Linear Stability Module

    s.linearize = true;

    %% Linearization setup
    disp("------------------------")
    disp("      Linearize L       ")
    disp("------------------------")

    %% Prepare base flow
    if s.shock.enabled
        s = CHECK_MESH_FOLLOWS_SHOCK(s, chemistry);
    end
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry);

    %% Construct Jacobian via finite differences
    tic
    L = LINEARIZE_EQUATIONS_NO_DISCONTINUITY(s, chemistry);
    toc

    %% Validate linearization
    CHECK_LINEARIZATION_L(L, s, chemistry);

    s.linearize = false;
end
