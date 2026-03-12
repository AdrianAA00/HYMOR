function [s, L_] = LINEARIZE_L_(L, s, chemistry, w_infty)
% LINEARIZE_L_  Construct the extended system matrix L_ including freestream coupling.
%
%   [s, L_] = LINEARIZE_L_(L, s, chemistry, w_infty)
%
%   Builds the augmented system matrix L_ that couples the downstream
%   Jacobian A with the upstream boundary Jacobian B and the freestream
%   frequency oscillation terms. The resulting block matrix has the form:
%
%       L_ = [L,             B,          B,          ..., B
%             0_(N*p x n),   i*w_1*I_p,  0,          ..., 0
%             0_(N*p x n),   0,          i*w_2*I_p,  ..., 0
%             ...
%             0_(N*p x n),   0,          0,          ..., i*w_N*I_p]
%
%   Inputs:
%       L         - Downstream Jacobian matrix (sparse, m x n)
%       s  - Solution structure with flow fields, mesh, and analysis flags
%       chemistry - Chemistry model structure
%       w_infty   - Vector of freestream disturbance frequencies
%
%   Outputs:
%       s - Updated s structure (linearize flag reset to false)
%       L_       - Extended system matrix (sparse); returns 0 if analysis is disabled
%
%   Notes:
%       - The upstream boundary operator B is computed via finite-difference
%         linearization (LINEARIZE_B) and validated (CHECK_LINEARIZATION_B).
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Freestream Module

    s.linearize = true;

    %% Display header
    fprintf("\n")
    disp("------------------------")
    disp("     Linearize L_       ")
    disp("------------------------")

    %% Prepare base flow
    if s.shock.enabled
        s = CHECK_MESH_FOLLOWS_SHOCK(s, chemistry);
    end
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry);

    %% Linearize upstream boundary operator B
    tic
    B = LINEARIZE_B(s, s.stability_analysis.perturbation_magnitude, chemistry);
    toc
    CHECK_LINEARIZATION_B(B, s, chemistry);

    %% Assemble extended system matrix L_
    [m, n] = size(L);
    [~, p] = size(B);
    N = length(w_infty);

    % Frequency diagonal: kron(diag(-i*w), I_p)
    freq_diag = kron(spdiags(-1i*w_infty, 0, N, N), speye(p));

    L_ = [L, repmat(B, 1, N);
            sparse(N*p, n), freq_diag];

    s.linearize = false;
end
