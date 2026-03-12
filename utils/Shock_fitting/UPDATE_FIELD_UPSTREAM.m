function s = UPDATE_FIELD_UPSTREAM(s)
% UPDATE_FIELD_UPSTREAM - Set upstream field cells to freestream values with perturbations.
%
%   s = UPDATE_FIELD_UPSTREAM(s)
%
%   Fills all grid cells upstream of the shock (beyond the two ghost-cell
%   layers used for extrapolation) with either a prescribed perturbed
%   upstream state or a uniform freestream plus travelling-wave
%   perturbation. Also sets gamma_star and cv_star to freestream values
%   for non-equilibrium chemistry cases.
%
%   Inputs:
%       s (struct) - Solution structure containing at minimum:
%                    .mesh.Nchi              - Number of streamwise cells
%                    .shock.cell_indices     - (Nx x 1) shocked-cell indices
%                    .freestream.rho_0, .freestream.rho_u_0, .freestream.rho_v_0,
%                    .freestream.rho_E_0       - Scalar uniform freestream values
%                    .freestream.disturbance.k_y, .freestream.disturbance.k_x
%                                            - Perturbation wavenumbers (tangential, normal)
%                    .freestream.disturbance.amplitude     - (4 x 1) perturbation amplitudes
%                    .mesh.x_Ext, .mesh.y_Ext - Extended grid coordinates
%                    .time_integration.t      - Current simulation time
%                    .chemistry.chemical_equilibrium - (logical) equilibrium chemistry flag
%                    .freestream.gamma_star   - Freestream thermodynamic properties
%                    (optionally) .freestream.rho_0_p, etc. - Prescribed perturbed fields
%
%   Outputs:
%       s (struct) - Solution with upstream cells populated.
%
% Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

    u_y = s.freestream.rho_v_0 / s.freestream.rho_0;
    u_x = s.freestream.rho_u_0 / s.freestream.rho_0;
    
    % Create a 2D Logical Mask
    % Find grid dimensions
    [Nx, Ny] = size(s.var.rho);
    [col_grid, ~] = meshgrid(1:Ny, 1:Nx);
    
    % Initialize cutoff columns to infinity so we only target the active Nchi rows
    cutoff_cols = inf(Nx, 1);
    cutoff_cols(2:s.mesh.Nchi+1) = s.shock.cell_indices(:, 1) + 2;
    
    % Create boolean mask: true ONLY for cells upstream of the shock
    upstream_mask = col_grid > cutoff_cols;

    if isfield(s, 'rho_0_upstream_p')
        %% Prescribed perturbed upstream state
        
        % Expand the 1D prescribed vectors to full 2D grids (row-dependent)
        rho_p_full  = zeros(Nx, Ny); rho_p_full(2:s.mesh.Nchi+1, :)  = repmat(s.freestream.rho_0_p, 1, Ny);
        rhou_p_full = zeros(Nx, Ny); rhou_p_full(2:s.mesh.Nchi+1, :) = repmat(s.freestream.rho_u_0_p, 1, Ny);
        rhov_p_full = zeros(Nx, Ny); rhov_p_full(2:s.mesh.Nchi+1, :) = repmat(s.freestream.rho_v_0_p, 1, Ny);
        rhoE_p_full = zeros(Nx, Ny); rhoE_p_full(2:s.mesh.Nchi+1, :) = repmat(s.freestream.rho_E_0_p, 1, Ny);
        
        % Apply values using the logical mask (Blazing fast in MATLAB)
        s.var.rho(upstream_mask)   = rho_p_full(upstream_mask);
        s.var.rho_u(upstream_mask) = rhou_p_full(upstream_mask);
        s.var.rho_v(upstream_mask) = rhov_p_full(upstream_mask);
        s.var.rho_E(upstream_mask) = rhoE_p_full(upstream_mask);
        
        if ~s.chemistry.chemical_equilibrium
            s.var.gamma_star(upstream_mask) = s.freestream.gamma_star;
            s.var.cv_star(upstream_mask)    = 1;
        end
        
    else
        %% Uniform freestream plus travelling-wave perturbation
        
        % OPTIMIZATION 2: Extract ONLY the coordinates that need perturbation
        x_up = s.mesh.x_Ext(upstream_mask);
        y_up = s.mesh.y_Ext(upstream_mask);
        
        % OPTIMIZATION 3: Compute expensive trig functions ONLY on the upstream vector
        perturbation_up = cos(2 * pi * s.freestream.disturbance.k_y * (x_up - u_x * s.time_integration.t)) ...
                       .* sin(2 * pi * s.freestream.disturbance.k_x * (y_up - u_y * s.time_integration.t));
                       
        % Apply to the flow field arrays
        s.var.rho(upstream_mask)   = s.freestream.rho_0   + perturbation_up * s.freestream.disturbance.amplitude(1);
        s.var.rho_u(upstream_mask) = s.freestream.rho_u_0 + perturbation_up * s.freestream.disturbance.amplitude(2);
        s.var.rho_v(upstream_mask) = s.freestream.rho_v_0 + perturbation_up * s.freestream.disturbance.amplitude(3);
        s.var.rho_E(upstream_mask) = s.freestream.rho_E_0 + perturbation_up * s.freestream.disturbance.amplitude(4);
        
        if ~s.chemistry.chemical_equilibrium
            s.var.gamma_star(upstream_mask) = s.freestream.gamma_star;
            s.var.cv_star(upstream_mask)    = 1;
        end
    end
end