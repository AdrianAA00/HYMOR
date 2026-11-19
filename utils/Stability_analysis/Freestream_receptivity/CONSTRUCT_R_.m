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
%       R_infty - Sparse transformation matrix of size (nv*Nx x nv*Nx)

    Nx = s.mesh.Nchi;
    rho = s.freestream.rho_0 * ones(Nx, 1);
    R_infty = CHU_TRANSFORMATION(rho, s.freestream.rho_u_0 / s.freestream.rho_0 * ones(Nx, 1), ...
                                 s.freestream.rho_v_0 / s.freestream.rho_0 * ones(Nx, 1), ...
                                 s.freestream.p_0 * ones(Nx, 1), s.freestream.gamma_star * ones(Nx, 1), ...
                                 numel(PERTURBATION_FIELDS(s)));
end
