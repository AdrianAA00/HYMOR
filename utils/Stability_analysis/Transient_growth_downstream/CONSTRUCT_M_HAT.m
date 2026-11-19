function M_hat = CONSTRUCT_M_HAT(s, norms)
% CONSTRUCT_M_HAT  Build the initial-condition energy-weight matrix M_hat.
%
%   M_hat = CONSTRUCT_M_HAT(s, norms)
%
%   Constructs a sparse diagonal matrix M_hat that defines the energy norm
%   for the initial condition (denominator of the gain). Unlike CONSTRUCT_M,
%   this matrix does not apply a spatial mask and includes shock-position
%   degrees of freedom when enabled.
%
%   Inputs:
%       s - Solution structure containing flow fields and mesh data:
%                    .mesh.Nchi, .mesh.Neta          - Grid dimensions
%                    .gamma_star       - Ratio of specific heats (with ghost cells)
%                    .p                - Pressure field (with ghost cells)
%                    .rho, .rho_u, .rho_v - Conservative variables (with ghost cells)
%                    .a                - Speed of sound (with ghost cells)
%                    .volume           - Cell volumes
%                    .stability_analysis.perturb_shock - Flag for shock DoF
%                    .bt_area          - Cell boundary areas
%                    .shocked_cell_indices - Indices of shocked cells
%       norms    - Logical vector [pressure, u-momentum, v-momentum, entropy, (w-momentum)]
%                  selecting which energy components to include (CHU_NORM_WEIGHTS)
%
%   Outputs:
%       M_hat - Sparse diagonal energy-weight matrix of size (nv*Nx*Ny) or
%               (nv*Nx*Ny+Nx) when shock perturbation DoFs are included
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Downstream Module

    nv = numel(PERTURBATION_FIELDS(s));
    I = {2:s.mesh.Nchi+1, 2:s.mesh.Neta+1};
    W = CHU_NORM_WEIGHTS(s.var.rho(I{:}), s.var.a(I{:}), s.var.gamma_star(I{:}), s.var.p(I{:}), norms, nv);
    values = W .* s.mesh.volume(:);
    values = values(:);
    if s.stability_analysis.perturb_shock
        % Shock DoFs: half the area of the face crossed by the shock
        shock_faces = sub2ind(size(s.mesh.bt_area), (1:s.mesh.Nchi)', s.shock.cell_indices(:, 1));
        values = [values; s.mesh.bt_area(shock_faces) / 2];
    end
    M_hat = spdiags(values, 0, numel(values), numel(values));
end
