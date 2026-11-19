function M = CONSTRUCT_M(s, norms)
% CONSTRUCT_M  Build the diagonal energy-weight matrix M for downstream disturbances.
%
%   M = CONSTRUCT_M(s, norms)
%
%   Constructs a sparse diagonal matrix M that defines the Chu energy norm
%   for the downstream flow domain. Each contribution (pressure, x-momentum,
%   y-momentum, entropy and, for an azimuthal mode, azimuthal momentum) can
%   be toggled independently via the norms vector (CHU_NORM_WEIGHTS). A spatial mask may be
%   applied to exclude the stagnation region or select/exclude the
%   boundary layer.
%
%   Inputs:
%       s - Solution structure containing flow fields and mesh data:
%                    .mesh.Nchi, .mesh.Neta          - Grid dimensions
%                    .gamma_star       - Ratio of specific heats (with ghost cells)
%                    .p                - Pressure field (with ghost cells)
%                    .rho, .rho_u, .rho_v - Conservative variables (with ghost cells)
%                    .a                - Speed of sound (with ghost cells)
%                    .volume           - Cell volumes
%                    .flow_cells       - Active-cell mask
%                    .stability_analysis.perturb_shock - Flag for shock DoF
%       norms    - Logical vector [pressure, u-momentum, v-momentum, entropy, (w-momentum)]
%                  selecting which energy components to include
%
%   Outputs:
%       M - Sparse diagonal energy-weight matrix of size (nv*Nx*Ny) or
%           (nv*Nx*Ny+Nx) when shock perturbation DoFs are included
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Downstream Module

    nv = numel(PERTURBATION_FIELDS(s));
    I = {2:s.mesh.Nchi+1, 2:s.mesh.Neta+1};
    volume = s.mesh.volume .* s.shock.flow_cells;
    W = CHU_NORM_WEIGHTS(s.var.rho(I{:}), s.var.a(I{:}), s.var.gamma_star(I{:}), s.var.p(I{:}), norms, nv);
    values = W .* volume(:);
    n_shock = s.stability_analysis.perturb_shock * s.mesh.Nchi;
    M = spdiags([values(:); zeros(n_shock, 1)], 0, numel(values) + n_shock, numel(values) + n_shock);
end
