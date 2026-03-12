function M = CONSTRUCT_M(s, norms)
% CONSTRUCT_M  Build the diagonal energy-weight matrix M for downstream disturbances.
%
%   M = CONSTRUCT_M(s, norms)
%
%   Constructs a sparse diagonal matrix M that defines the Chu energy norm
%   for the downstream flow domain. Each of the four conserved-variable
%   contributions (pressure, x-momentum, y-momentum, entropy) can be
%   toggled independently via the norms vector. A spatial mask may be
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
%       norms    - 4-element logical vector [pressure, u-momentum, v-momentum, entropy]
%                  selecting which energy components to include
%
%   Outputs:
%       M - Sparse diagonal energy-weight matrix of size (4*Nx*Ny) or
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
    a = s.var.a(2:end-1, 2:end-1);
    flow_cells = s.shock.flow_cells;
    volume = s.mesh.volume .* flow_cells;
    indices = reshape(1:Nx*Ny, Nx, Ny);

    %% Initialize storage for sparse triplets
    M = zeros(4, Nx, Ny);
    M_m = zeros(4, Nx, Ny);
    M_n = zeros(4, Nx, Ny);

    %% Build spatial mask
    % s.no_stagnation = false;
    % s.no_BL         = false;
    % s.only_BL       = false;
    % [mask_no_stagnation, mask_no_BL, mask_only_BL] = CREATE_MASKS(s);
    mask = ones(Nx, Ny);

    %% Pressure norm contribution
    if norms(1) == true
        M(1, :, :) = rho .* a.^2 ./ (2 .* (gamma_star .* p).^2) .* volume .* mask;
        M_n(1, :, :) = 0*Nx*Ny + indices;
        M_m(1, :, :) = 0*Nx*Ny + indices;
    end

    %% X-momentum norm contribution
    if norms(2) == true
        M(2, :, :) = rho/2 .* volume .* mask;
        M_n(2, :, :) = 1*Nx*Ny + indices;
        M_m(2, :, :) = 1*Nx*Ny + indices;
    end

    %% Y-momentum norm contribution
    if norms(3) == true
        M(3, :, :) = rho/2 .* volume .* mask;
        M_n(3, :, :) = 2*Nx*Ny + indices;
        M_m(3, :, :) = 2*Nx*Ny + indices;
    end

    %% Entropy norm contribution
    if norms(4) == true
        M(4, :, :) = (gamma_star - 1) .* p ./ (2 * gamma_star) .* volume .* mask;
        M_n(4, :, :) = 3*Nx*Ny + indices;
        M_m(4, :, :) = 3*Nx*Ny + indices;
    end

    %% Assemble sparse matrix from triplets
    values = M(:);
    vector_indices_m = M_m(:);
    vector_indices_n = M_n(:);

    % Keep only non-zero index entries
    mask = (vector_indices_m ~= 0) & (vector_indices_n ~= 0);
    values = values(mask);
    vector_indices_m = vector_indices_m(mask);
    vector_indices_n = vector_indices_n(mask);

    if s.stability_analysis.perturb_shock
        M = sparse(vector_indices_m, vector_indices_n, values, 4*Nx*Ny+Nx, 4*Nx*Ny+Nx, size(vector_indices_m, 1));
    else
        M = sparse(vector_indices_m, vector_indices_n, values, 4*Nx*Ny, 4*Nx*Ny, size(vector_indices_m, 1));
    end
end
