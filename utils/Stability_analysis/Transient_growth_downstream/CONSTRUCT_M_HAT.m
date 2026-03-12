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
%       norms    - 4-element logical vector [pressure, u-momentum, v-momentum, entropy]
%                  selecting which energy components to include
%
%   Outputs:
%       M_hat - Sparse diagonal energy-weight matrix of size (4*Nx*Ny) or
%               (4*Nx*Ny+Nx) when shock perturbation DoFs are included
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
    volume = s.mesh.volume;
    indices = reshape(1:Nx*Ny, Nx, Ny);
    indices_r = reshape(1:Nx, Nx, 1);

    %% Initialize storage for sparse triplets
    M = zeros(4, Nx, Ny);
    M_m = zeros(4, Nx, Ny);
    M_n = zeros(4, Nx, Ny);

    %% Pressure norm contribution
    if norms(1) == true
        M(1, :, :) = rho .* a.^2 ./ (2 .* (gamma_star .* p).^2) .* volume;
        M_n(1, :, :) = 0*Nx*Ny + indices;
        M_m(1, :, :) = 0*Nx*Ny + indices;
    end

    %% X-momentum norm contribution
    if norms(2) == true
        M(2, :, :) = rho/2 .* volume;
        M_n(2, :, :) = 1*Nx*Ny + indices;
        M_m(2, :, :) = 1*Nx*Ny + indices;
    end

    %% Y-momentum norm contribution
    if norms(3) == true
        M(3, :, :) = rho/2 .* volume;
        M_n(3, :, :) = 2*Nx*Ny + indices;
        M_m(3, :, :) = 2*Nx*Ny + indices;
    end

    %% Entropy norm contribution
    if norms(4) == true
        M(4, :, :) = (gamma_star - 1) .* p ./ (2 * gamma_star) .* volume;
        M_n(4, :, :) = 3*Nx*Ny + indices;
        M_m(4, :, :) = 3*Nx*Ny + indices;
    end

    %% Flatten triplets
    values = M(:);
    vector_indices_m = M_m(:);
    vector_indices_n = M_n(:);

    %% Shock-position degrees of freedom
    if s.stability_analysis.perturb_shock
        M_r = zeros(Nx, 1);
        for i = 1:Nx
            M_r(i) = s.mesh.bt_area(i, s.shock.cell_indices(i, 1)) / 2;
        end

        M_r_n = 4*Nx*Ny + indices_r(:);
        M_r_m = 4*Nx*Ny + indices_r(:);

        values = [values; M_r];
        vector_indices_m = [vector_indices_m; M_r_m];
        vector_indices_n = [vector_indices_n; M_r_n];
    end

    %% Assemble sparse matrix from triplets
    % Keep only non-zero index entries
    mask = (vector_indices_m ~= 0) & (vector_indices_n ~= 0);
    values = values(mask);
    vector_indices_m = vector_indices_m(mask);
    vector_indices_n = vector_indices_n(mask);

    if s.stability_analysis.perturb_shock
        M_hat = sparse(vector_indices_m, vector_indices_n, values, 4*Nx*Ny+Nx, 4*Nx*Ny+Nx, size(vector_indices_m, 1));
    else
        M_hat = sparse(vector_indices_m, vector_indices_n, values, 4*Nx*Ny, 4*Nx*Ny, size(vector_indices_m, 1));
    end
end
