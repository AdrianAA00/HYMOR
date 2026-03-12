function R = CONSTRUCT_R(s)
% CONSTRUCT_R  Build the variable-transformation matrix R.
%
%   R = CONSTRUCT_R(s)
%
%   Constructs a sparse matrix R that transforms the conservative state
%   vector [rho, rho*u, rho*v, rho*E]' into the primitive/entropy
%   variables [p, u, v, S]' used in the Chu energy norm. When shock
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
%       R - Sparse transformation matrix of size (4*Nx*Ny) or
%           (4*Nx*Ny+Nx) when shock perturbation DoFs are included
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Downstream Module

    %% Extract flow quantities (interior cells only)
    Nx = s.mesh.Nchi;
    Ny = s.mesh.Neta;
    gamma_star = s.var.gamma_star(2:end-1, 2:end-1);
    p = s.var.p(2:end-1, 2:end-1);
    u = s.var.rho_u(2:end-1, 2:end-1) ./ s.var.rho(2:end-1, 2:end-1);
    v = s.var.rho_v(2:end-1, 2:end-1) ./ s.var.rho(2:end-1, 2:end-1);
    rho = s.var.rho(2:end-1, 2:end-1);

    %% Define index offsets for conservative and primitive variables
    idx_rho   = 0;
    idx_rho_u = Nx * Ny;
    idx_rho_v = 2 * Nx * Ny;
    idx_rho_E = 3 * Nx * Ny;
    idx_p = 0;
    idx_u = Nx * Ny;
    idx_v = 2 * Nx * Ny;
    idx_S = 3 * Nx * Ny;
    idx_r = 4 * Nx * Ny;
    indices   = reshape(1:Nx*Ny, Nx, Ny);
    indices_r = reshape(1:Nx, Nx, 1);

    %% Initialize storage for sparse triplets
    R   = zeros(4, 4, Nx, Ny);
    R_n = zeros(4, 4, Nx, Ny);
    R_m = zeros(4, 4, Nx, Ny);

    %% Row 1: Pressure from conservative variables
    R(1, 1, :, :) = (gamma_star - 1) .* (u.^2 + v.^2) ./ 2;
    R_n(1, 1, :, :) = idx_rho + indices;
    R_m(1, 1, :, :) = idx_p + indices;

    R(1, 2, :, :) = -u .* (gamma_star - 1);
    R_n(1, 2, :, :) = idx_rho_u + indices;
    R_m(1, 2, :, :) = idx_p + indices;

    R(1, 3, :, :) = -v .* (gamma_star - 1);
    R_n(1, 3, :, :) = idx_rho_v + indices;
    R_m(1, 3, :, :) = idx_p + indices;

    R(1, 4, :, :) = gamma_star - 1;
    R_n(1, 4, :, :) = idx_rho_E + indices;
    R_m(1, 4, :, :) = idx_p + indices;

    %% Row 2: X-velocity from conservative variables
    R(2, 1, :, :) = -u ./ rho;
    R_n(2, 1, :, :) = idx_rho + indices;
    R_m(2, 1, :, :) = idx_u + indices;

    R(2, 2, :, :) = 1 ./ rho;
    R_n(2, 2, :, :) = idx_rho_u + indices;
    R_m(2, 2, :, :) = idx_u + indices;

    %% Row 3: Y-velocity from conservative variables
    R(3, 1, :, :) = -v ./ rho;
    R_n(3, 1, :, :) = idx_rho + indices;
    R_m(3, 1, :, :) = idx_v + indices;

    R(3, 3, :, :) = 1 ./ rho;
    R_n(3, 3, :, :) = idx_rho_v + indices;
    R_m(3, 3, :, :) = idx_v + indices;

    %% Row 4: Entropy from conservative variables
    R(4, 1, :, :) = (u.^2 + v.^2) ./ (2*p) - gamma_star ./ (gamma_star - 1) ./ rho;
    R_n(4, 1, :, :) = idx_rho + indices;
    R_m(4, 1, :, :) = idx_S + indices;

    R(4, 2, :, :) = -u ./ p;
    R_n(4, 2, :, :) = idx_rho_u + indices;
    R_m(4, 2, :, :) = idx_S + indices;

    R(4, 3, :, :) = -v ./ p;
    R_n(4, 3, :, :) = idx_rho_v + indices;
    R_m(4, 3, :, :) = idx_S + indices;

    R(4, 4, :, :) = 1 ./ p;
    R_n(4, 4, :, :) = idx_rho_E + indices;
    R_m(4, 4, :, :) = idx_S + indices;

    %% Flatten triplets
    values = R(:);
    vector_indices_m = R_m(:);
    vector_indices_n = R_n(:);

    %% Shock-position degrees of freedom (identity block)
    if s.stability_analysis.perturb_shock
        R_r = ones(Nx, 1);
        R_r_n = idx_r + indices_r(:);
        R_r_m = idx_r + indices_r(:);

        values = [values; R_r];
        vector_indices_m = [vector_indices_m; R_r_m];
        vector_indices_n = [vector_indices_n; R_r_n];
    end

    %% Assemble sparse matrix from triplets
    % Keep only non-zero index entries
    mask = (vector_indices_m ~= 0) & (vector_indices_n ~= 0);
    values = values(mask);
    vector_indices_m = vector_indices_m(mask);
    vector_indices_n = vector_indices_n(mask);

    if s.stability_analysis.perturb_shock
        R = sparse(vector_indices_m, vector_indices_n, values, 4*Nx*Ny+Nx, 4*Nx*Ny+Nx, size(vector_indices_m, 1));
    else
        R = sparse(vector_indices_m, vector_indices_n, values, 4*Nx*Ny, 4*Nx*Ny, size(vector_indices_m, 1));
    end
end
