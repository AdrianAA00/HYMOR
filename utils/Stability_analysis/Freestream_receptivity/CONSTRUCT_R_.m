function R_ = CONSTRUCT_R_(s, w_infty)
% CONSTRUCT_R_  Build the extended variable-transformation matrix R_ for the coupled system.
%
%   R_ = CONSTRUCT_R_(s, w_infty)
%
%   Constructs the block-diagonal transformation matrix for the extended
%   state vector that includes both downstream and freestream degrees of
%   freedom. The downstream block uses the standard R matrix; each
%   freestream frequency block uses R_infty (uniform freestream transform).
%
%   Structure:
%       R_ = blkdiag(R, R_infty, R_infty, ..., R_infty)
%       with N = length(w_infty) copies of R_infty.
%
%   Inputs:
%       s - Solution structure with flow fields and mesh data
%       w_infty  - Vector of freestream disturbance frequencies
%
%   Outputs:
%       R_ - Sparse block-diagonal transformation matrix for the extended system
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Freestream Module

    %% Build downstream and freestream transformation matrices
    R = CONSTRUCT_R(s);
    R_infty = CONSTRUCT_R_INFTY(s);

    %% Assemble block-diagonal matrix
    N = length(w_infty);
    R_blocks = repmat({R_infty}, 1, N);
    R_ = blkdiag(R, R_blocks{:});
end


function R_infty = CONSTRUCT_R_INFTY(s)
% CONSTRUCT_R_INFTY  Build the variable-transformation matrix for uniform freestream.
%
%   R_infty = CONSTRUCT_R_INFTY(s)
%
%   Constructs the sparse transformation matrix from conservative to
%   primitive/entropy variables evaluated at the uniform freestream state.
%   This is the freestream analogue of CONSTRUCT_R.
%
%   Inputs:
%       s - Solution structure with upstream flow state
%
%   Outputs:
%       R_infty - Sparse transformation matrix of size (4*Nx x 4*Nx)

    %% Extract uniform freestream quantities
    Nx = s.mesh.Nchi;
    gamma_star = s.freestream.gamma_star * ones(Nx, 1);
    p = s.freestream.p_0 * ones(Nx, 1);
    rho = s.freestream.rho_0 * ones(Nx, 1);
    u = s.freestream.rho_u_0 / s.freestream.rho_0 * ones(Nx, 1);
    v = s.freestream.rho_v_0 / s.freestream.rho_0 * ones(Nx, 1);

    %% Define index offsets
    idx_rho   = 0;
    idx_rho_u = Nx;
    idx_rho_v = 2 * Nx;
    idx_rho_E = 3 * Nx;
    idx_p = 0;
    idx_u = Nx;
    idx_v = 2 * Nx;
    idx_S = 3 * Nx;
    indices = reshape(1:Nx, Nx, 1);

    %% Initialize storage for sparse triplets
    R   = zeros(4, 4, Nx);
    R_n = zeros(4, 4, Nx);
    R_m = zeros(4, 4, Nx);

    %% Row 1: Pressure from conservative variables
    R(1, 1, :) = (gamma_star - 1) .* (u.^2 + v.^2) ./ 2;
    R_n(1, 1, :) = idx_rho + indices;
    R_m(1, 1, :) = idx_p + indices;

    R(1, 2, :) = -u .* (gamma_star - 1);
    R_n(1, 2, :) = idx_rho_u + indices;
    R_m(1, 2, :) = idx_p + indices;

    R(1, 3, :) = -v .* (gamma_star - 1);
    R_n(1, 3, :) = idx_rho_v + indices;
    R_m(1, 3, :) = idx_p + indices;

    R(1, 4, :) = gamma_star - 1;
    R_n(1, 4, :) = idx_rho_E + indices;
    R_m(1, 4, :) = idx_p + indices;

    %% Row 2: X-velocity from conservative variables
    R(2, 1, :) = -u ./ rho;
    R_n(2, 1, :) = idx_rho + indices;
    R_m(2, 1, :) = idx_u + indices;

    R(2, 2, :) = 1 ./ rho;
    R_n(2, 2, :) = idx_rho_u + indices;
    R_m(2, 2, :) = idx_u + indices;

    %% Row 3: Y-velocity from conservative variables
    R(3, 1, :) = -v ./ rho;
    R_n(3, 1, :) = idx_rho + indices;
    R_m(3, 1, :) = idx_v + indices;

    R(3, 3, :) = 1 ./ rho;
    R_n(3, 3, :) = idx_rho_v + indices;
    R_m(3, 3, :) = idx_v + indices;

    %% Row 4: Entropy from conservative variables
    R(4, 1, :) = (u.^2 + v.^2) ./ (2*p) - gamma_star ./ (gamma_star - 1) ./ rho;
    R_n(4, 1, :) = idx_rho + indices;
    R_m(4, 1, :) = idx_S + indices;

    R(4, 2, :) = -u ./ p;
    R_n(4, 2, :) = idx_rho_u + indices;
    R_m(4, 2, :) = idx_S + indices;

    R(4, 3, :) = -v ./ p;
    R_n(4, 3, :) = idx_rho_v + indices;
    R_m(4, 3, :) = idx_S + indices;

    R(4, 4, :) = 1 ./ p;
    R_n(4, 4, :) = idx_rho_E + indices;
    R_m(4, 4, :) = idx_S + indices;

    %% Assemble sparse matrix from triplets
    values = R(:);
    vector_indices_m = R_m(:);
    vector_indices_n = R_n(:);

    % Keep only non-zero index entries
    mask = (vector_indices_m ~= 0) & (vector_indices_n ~= 0);
    values = values(mask);
    vector_indices_m = vector_indices_m(mask);
    vector_indices_n = vector_indices_n(mask);

    R_infty = sparse(vector_indices_m, vector_indices_n, values, 4*Nx, 4*Nx, size(vector_indices_m, 1));
end
