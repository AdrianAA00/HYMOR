function R = CONSTRUCT_R(s)
% CONSTRUCT_R  Build the variable-transformation matrix R.
%
%   R = CONSTRUCT_R(s)
%
%   Constructs a sparse matrix R that transforms the conservative state
%   vector [rho, rho*u, rho*v, rho*E, (rho*w)]' into the primitive/entropy
%   variables [p, u, v, S, (w)]' used in the Chu energy norm
%   (CHU_TRANSFORMATION). When shock
%   perturbation is enabled, an identity block is appended for the
%   shock-position degrees of freedom.
%
%   Inputs:
%       s - Solution structure containing flow fields and mesh data:
%                    .mesh.Nchi, .mesh.Neta          - Grid dimensions
%                    .gamma_star       - Ratio of specific heats (with ghost cells)
%                    .p                - Pressure field (with ghost cells)
%                    .rho, .rho_u, .rho_v - Conservative variables (with ghost cells)
%                    .stability_analysis.perturb_shock - Flag for shock DoF
%
%   Outputs:
%       R - Sparse transformation matrix of size (nv*Nx*Ny) or
%           (nv*Nx*Ny+Nx) when shock perturbation DoFs are included
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Downstream Module

    nv = numel(PERTURBATION_FIELDS(s));
    I = {2:s.mesh.Nchi+1, 2:s.mesh.Neta+1};
    rho = s.var.rho(I{:});
    R = CHU_TRANSFORMATION(rho, s.var.rho_u(I{:}) ./ rho, s.var.rho_v(I{:}) ./ rho, ...
                           s.var.p(I{:}), s.var.gamma_star(I{:}), nv);
    if s.stability_analysis.perturb_shock
        R = blkdiag(R, speye(s.mesh.Nchi));
    end
end
