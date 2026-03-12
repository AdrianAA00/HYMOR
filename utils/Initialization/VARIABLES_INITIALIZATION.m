function s = VARIABLES_INITIALIZATION(s)
% VARIABLES_INITIALIZATION  Pre-allocate all field arrays in the s struct.
%   Zeros out coordinate arrays, mesh metric arrays, conserved-variable
%   fields (with ghost cells), chemistry fields, and shock-tracking arrays
%   based on the grid dimensions Nx and Ny stored in the s struct.
%
%   s = VARIABLES_INITIALIZATION(s)
%
%   Inputs:
%       s - Solution struct containing at minimum the grid
%                  dimensions s.mesh.Nchi and s.mesh.Neta
%
%   Outputs:
%       s - Solution struct with all field arrays initialized to
%                  zero at the appropriate sizes
%
%   Notes:
%       - Wall arrays have size (Nx+1, 1) corresponding to cell edges.
%       - Interior cell-center arrays have size (Nx, Ny).
%       - Extended arrays (with ghost cells) have size (Nx+2, Ny+2).
%       - Cell-corner arrays have size (Nx+1, Ny+1).
%       - Face-area and normal arrays are sized for left-right (Nx+1, Ny),
%         bottom-top (Nx, Ny+1), and front-back (Nx, Ny) faces.
%       - Shock arrays are sized (Nx, 1) for per-column shock data.
%
% Part of: Hypersonics Stability MATLAB Solver - Initialization Module

    %% Wall coordinates and normals
    s.mesh.x_wall        = zeros(s.mesh.Nchi + 1, 1);
    s.mesh.y_wall        = zeros(s.mesh.Nchi + 1, 1);
    s.mesh.eta_wall        = zeros(s.mesh.Nchi + 1, 1);
    s.mesh.chi_wall    = zeros(s.mesh.Nchi + 1, 1);
    s.mesh.x_wall_normal = zeros(s.mesh.Nchi + 1, 1);
    s.mesh.y_wall_normal = zeros(s.mesh.Nchi + 1, 1);

    %% Cell-center coordinates (interior)
    s.mesh.x     = zeros(s.mesh.Nchi, s.mesh.Neta);
    s.mesh.y     = zeros(s.mesh.Nchi, s.mesh.Neta);
    s.mesh.eta     = zeros(s.mesh.Nchi, s.mesh.Neta);
    s.mesh.chi = zeros(s.mesh.Nchi, s.mesh.Neta);

    %% Cell-center coordinates (extended with ghost cells)
    s.mesh.x_Ext = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);
    s.mesh.y_Ext = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);

    %% Cell-corner coordinates
    s.mesh.x_corner = zeros(s.mesh.Nchi + 1, s.mesh.Neta + 1);
    s.mesh.y_corner = zeros(s.mesh.Nchi + 1, s.mesh.Neta + 1);

    %% Cell-face areas
    s.mesh.lr_area = zeros(s.mesh.Nchi + 1, s.mesh.Neta);    % left-right faces
    s.mesh.bt_area = zeros(s.mesh.Nchi, s.mesh.Neta + 1);    % bottom-top faces
    s.mesh.fb_area = zeros(s.mesh.Nchi, s.mesh.Neta);        % front-back faces (3D axisymmetric)

    %% Cell-face normals
    s.mesh.lr_x_normal = zeros(s.mesh.Nchi + 1, s.mesh.Neta);  % left-right x-component
    s.mesh.lr_y_normal = zeros(s.mesh.Nchi + 1, s.mesh.Neta);  % left-right y-component
    s.mesh.bt_x_normal = zeros(s.mesh.Nchi, s.mesh.Neta + 1);  % bottom-top x-component
    s.mesh.bt_y_normal = zeros(s.mesh.Nchi, s.mesh.Neta + 1);  % bottom-top y-component

    %% Cell volumes
    s.mesh.volume = zeros(s.mesh.Nchi, s.mesh.Neta);

    %% Conserved variables (extended with ghost cells for BC)
    s.var.rho        = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);
    s.var.rho_u      = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);
    s.var.rho_v      = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);
    s.var.rho_E      = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);
    s.var.p          = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);

    %% Conserved-variable fluxes (interior cells only)
    s.flux.rho   = zeros(s.mesh.Nchi, s.mesh.Neta);
    s.flux.rho_u = zeros(s.mesh.Nchi, s.mesh.Neta);
    s.flux.rho_v = zeros(s.mesh.Nchi, s.mesh.Neta);
    s.flux.rho_E = zeros(s.mesh.Nchi, s.mesh.Neta);

    %% Chemistry field variables (extended with ghost cells)
    s.var.gamma_star      = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);
    s.var.gamma_star_eq   = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);
    s.flux.gamma_star = zeros(s.mesh.Nchi, s.mesh.Neta);
    s.var.cv_star         = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);
    s.var.cv_star_eq      = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);
    s.flux.cv_star    = zeros(s.mesh.Nchi, s.mesh.Neta);

    %% Transport property fields (extended with ghost cells)
    s.var.mu_star  = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);
    s.var.k_star   = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);
    s.var.Re_flow  = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);
    s.var.Pr_flow  = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 2);

    %% Shock-tracking arrays
    s.shock.cell_indices = zeros(s.mesh.Nchi, 1);
    s.shock.cells        = zeros(s.mesh.Nchi, s.mesh.Neta);
    s.shock.flow_cells     = ones(s.mesh.Nchi, s.mesh.Neta);
    s.shock.flow_cells_E   = ones(s.mesh.Nchi, s.mesh.Neta);
    s.shock.speed_x        = zeros(s.mesh.Nchi, 1);
    s.shock.speed_y        = zeros(s.mesh.Nchi, 1);
    s.shock.points_x       = zeros(s.mesh.Nchi, 1);
    s.shock.points_y       = zeros(s.mesh.Nchi, 1);
    s.shock.points_eta       = zeros(s.mesh.Nchi, 1);
    s.shock.points_chi   = zeros(s.mesh.Nchi, 1);
    s.shock.beta           = zeros(s.mesh.Nchi, 1);
    s.shock.M              = zeros(s.mesh.Nchi, 1);
end
