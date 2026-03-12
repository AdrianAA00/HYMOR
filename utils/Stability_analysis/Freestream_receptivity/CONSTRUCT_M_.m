function M_ = CONSTRUCT_M_(s, norms, w_infty)
% CONSTRUCT_M_  Build the extended energy-weight matrix M_ for the coupled system.
%
%   M_ = CONSTRUCT_M_(s, norms, w_infty)
%
%   Constructs the block-diagonal energy-weight matrix for the extended
%   state vector that includes both downstream and freestream degrees of
%   freedom. The downstream block uses the standard M matrix; the
%   freestream frequency blocks are zeroed out (energy is measured only
%   in the downstream domain).
%
%   Structure:
%       M_ = blkdiag(M, 0, 0, ..., 0)
%       where there are N = length(w_infty) zero blocks of size 4*Nx.
%
%   Inputs:
%       s - Solution structure with flow fields and mesh data
%       norms    - 4-element logical vector [pressure, u-momentum, v-momentum, entropy]
%       w_infty  - Vector of freestream disturbance frequencies
%
%   Outputs:
%       M_ - Sparse block-diagonal energy-weight matrix for the extended system
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Freestream Module

    %% Build downstream energy-weight matrix
    M = CONSTRUCT_M(s, norms);

    %% Append zero blocks for freestream frequency DoFs
    N = length(w_infty);
    M_ = blkdiag(M, sparse(4 * s.mesh.Nchi * N, 4 * s.mesh.Nchi * N));
end
