function s = APPLY_BOUNDARY_CONDITIONS(s, chemistry)
% APPLY_BOUNDARY_CONDITIONS - Enforce boundary conditions on all domain edges
%
% Syntax:
%   s = APPLY_BOUNDARY_CONDITIONS(s, chemistry)
%
% Description:
%   Applies ghost-cell boundary conditions on the four edges of a 2D
%   structured grid (eta0, eta1, chi0, chi1). Each edge can be independently
%   assigned a boundary condition type through the fields s.boundary_conditions.boundary_eta0,
%   s.boundary_conditions.boundary_eta1, s.boundary_conditions.boundary_chi0, and s.boundary_conditions.boundary_chi1. Ghost cell values
%   for density, momentum, total energy, pressure, temperature, and
%   non-equilibrium chemistry variables are set according to the
%   selected condition.
%
% Inputs:
%   s         - (struct) Solution structure containing:
%                 - Conservative variables: rho, rho_u, rho_v, rho_E
%                 - Thermodynamic variables: p, T, a
%                 - Non-equilibrium: gamma_star, cv_star (if applicable)
%                 - Geometry: bt_x_normal, bt_y_normal, lr_x_normal, lr_y_normal
%                 - Boundary type strings: boundary_eta0, boundary_eta1,
%                   boundary_chi0, boundary_chi1
%                 - Freestream/upstream references: rho_0, rho_u_0, etc.
%   chemistry - (struct) Chemistry model data (unused directly but
%                 passed for interface consistency)
%
% Outputs:
%   s         - (struct) Updated s structure with ghost cell
%                 values populated for all boundary edges
%
% Supported boundary conditions (common to all edges):
%   inflow_subsonic, inflow_supersonic, outflow_supersonic,
%   outflow_supersonic_1st, outflow_subsonic, periodic, symmetry,
%   outflow_NRCBC, no_slip_adiabatic, no_slip_isothermal
%
% Side-specific boundary conditions:
%   eta1 edge: shock, lid_driven_cavity
%
% Notes:
%   - Ghost cells use linear extrapolation (2nd order) or mirroring
%     depending on the boundary type.
%   - Periodic boundaries copy values from the opposite interior cell.
%   - Non-equilibrium chemistry variables (gamma_star, cv_star) are
%     handled when s.chemistry.chemical_equilibrium is false.
%   - For eta faces, bt_normal is used; for chi faces, lr_normal is used.
%   - NRCBC uses LODI (Poinsot & Lele): L1 = 0 => dp/dx_n = rho*c*du_n/dx_n.
%     Checks supersonic/subsonic per cell: u_n >= a vs u_n < a.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Left boundary (chi0)
    switch s.boundary_conditions.boundary_chi0.name
        case 'inflow_subsonic'
            s.var.rho(1,:) = s.freestream.rho_0;
            s.var.rho_u(1,:) = s.freestream.rho_u_0;
            s.var.rho_v(1,:) = s.freestream.rho_v_0;
            s.var.rho_E(1,:) = 2 * s.var.rho_E(2,:) - s.var.rho_E(3,:);
            s.var.p(1,:) = 2 * s.var.p(2,:) - s.var.p(3,:);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(1,:) = s.var.gamma_star_eq(1,:);
                s.var.cv_star(1,:) = s.var.cv_star_eq(1,:);
            end

        case 'inflow_supersonic'
            s.var.rho(1,:) = s.freestream.rho_0;
            s.var.rho_u(1,:) = s.freestream.rho_u_0;
            s.var.rho_v(1,:) = s.freestream.rho_v_0;
            s.var.rho_E(1,:) = s.freestream.rho_E_0;
            s.var.p(1,:) = s.freestream.p_0;

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(1,:) = s.var.gamma_star_eq(1,:);
                s.var.cv_star(1,:) = s.var.cv_star_eq(1,:);
            end

        case 'no_slip_adiabatic'
            s.var.rho(1,:) = s.var.rho(2,:);
            s.var.rho_u(1,:) = -s.var.rho_u(2,:);
            s.var.rho_v(1,:) = -s.var.rho_v(2,:);
            s.var.rho_E(1,:) = 2 * s.var.rho_E(2,:) - s.var.rho_E(3,:);
            s.var.p(1,:) = 2 * s.var.p(2,:) - s.var.p(3,:);
            s.var.T(1,:) = s.var.T(2,:);  % No temperature gradient for adiabatic wall

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(1,:) = 2 * s.var.gamma_star(2,:) - s.var.gamma_star(3,:);
                s.var.cv_star(1,:) = 2 * s.var.cv_star(2,:) - s.var.cv_star(3,:);
            end

        case 'no_slip_isothermal'
            s.var.rho(1,:) = s.var.rho(2,:);
            s.var.rho_u(1,:) = -s.var.rho_u(2,:);
            s.var.rho_v(1,:) = -s.var.rho_v(2,:);
            s.var.rho_E(1,:) = 2 * s.var.rho_E(2,:) - s.var.rho_E(3,:);
            s.var.p(1,:) = 2 * s.var.p(2,:) - s.var.p(3,:);
            T_w = s.boundary_conditions.boundary_chi0.Tw * s.freestream.cv / s.freestream.energy_factor; % Non-dimensional wall temperature
            s.var.T(1,:) = 2 * T_w - s.var.T(2,:);  % Set wall temperature

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(1,:) = s.var.gamma_star(2,:);
                s.var.cv_star(1,:) = s.var.cv_star(2,:);
            end

        case 'outflow_subsonic'
            s.var.rho(1,:) = 2 * s.var.rho(2,:) - s.var.rho(3,:);
            s.var.rho_E(1,:) = 2 * s.var.rho_E(2,:) - s.var.rho_E(3,:);
            s.var.rho_u(1,:) = 2 * s.var.rho_u(2,:) - s.var.rho_u(3,:);
            s.var.rho_v(1,:) = 2 * s.var.rho_v(2,:) - s.var.rho_v(3,:);
            s.var.p(1,:) = s.var.p(2,:); % zero pressure gradient

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(1,:) = s.var.gamma_star(2,:);
                s.var.cv_star(1,:) = s.var.cv_star(2,:);
            end

        case 'outflow_supersonic'
            s.var.rho(1,:) = 2 * s.var.rho(2,:) - s.var.rho(3,:);
            s.var.rho_E(1,:) = 2 * s.var.rho_E(2,:) - s.var.rho_E(3,:);
            s.var.rho_u(1,:) = 2 * s.var.rho_u(2,:) - s.var.rho_u(3,:);
            s.var.rho_v(1,:) = 2 * s.var.rho_v(2,:) - s.var.rho_v(3,:);
            s.var.p(1,:) = 2 * s.var.p(2,:) - s.var.p(3,:);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(1,:) = 2 * s.var.gamma_star(2,:) - s.var.gamma_star(3,:);
                s.var.cv_star(1,:) = 2 * s.var.cv_star(2,:) - s.var.cv_star(3,:);
            end

        case 'outflow_supersonic_1st'
            s.var.rho_u(1,:) = s.var.rho_u(2,:);
            s.var.rho_v(1,:) = s.var.rho_v(2,:);
            s.var.rho(1,:) = s.var.rho(2,:);
            s.var.rho_E(1,:) = s.var.rho_E(2,:);
            s.var.p(1,:) = s.var.p(2,:);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(1,:) = s.var.gamma_star(2,:);
                s.var.cv_star(1,:) = s.var.cv_star(2,:);
            end

        case 'outflow_NRCBC'
            % LODI boundary conditions (Poinsot & Lele): L1 = 0
            % Extrapolate all variables (supersonic default)
            s.var.rho(1,:) = 2 * s.var.rho(2,:) - s.var.rho(3,:);
            s.var.rho_u(1,:) = 2 * s.var.rho_u(2,:) - s.var.rho_u(3,:);
            s.var.rho_v(1,:) = 2 * s.var.rho_v(2,:) - s.var.rho_v(3,:);
            s.var.rho_E(1,:) = 2 * s.var.rho_E(2,:) - s.var.rho_E(3,:);
            s.var.p(1,:) = 2 * s.var.p(2,:) - s.var.p(3,:);

            % Sound speed and density at first interior cell
            a = s.var.a(2,2:end-1);
            rho_int = s.var.rho(2,2:end-1);

            % Normal velocity at ghost cell, positive exiting domain
            u_n_ghost = -(s.var.rho_u(1,2:end-1) .* s.mesh.lr_x_normal(1,:) + ...
                          s.var.rho_v(1,2:end-1) .* s.mesh.lr_y_normal(1,:)) ./ s.var.rho(1,2:end-1);

            % Normal velocity at first interior cell
            u_n_int = -(s.var.rho_u(2,2:end-1) .* s.mesh.lr_x_normal(1,:) + ...
                        s.var.rho_v(2,2:end-1) .* s.mesh.lr_y_normal(1,:)) ./ rho_int;

            % LODI correction for subsonic cells (L1 = 0 => dp/dx_n = rho*c*du_n/dx_n)
            subsonic = u_n_ghost < a;
            if any(subsonic)
                p_int = s.var.p(2,2:end-1);
                p_lodi = p_int + rho_int .* a .* (u_n_ghost - u_n_int);
                s.var.p(1,2:end-1) = subsonic .* p_lodi + (~subsonic) .* s.var.p(1,2:end-1);
                rho_e = s.var.p(1,2:end-1) ./ (s.var.gamma_star(1,2:end-1) - 1);
                rho_E_lodi = rho_e + 0.5 * (s.var.rho_u(1,2:end-1).^2 + s.var.rho_v(1,2:end-1).^2) ./ s.var.rho(1,2:end-1);
                s.var.rho_E(1,2:end-1) = subsonic .* rho_E_lodi + (~subsonic) .* s.var.rho_E(1,2:end-1);
                            s.p_int = p_int;
                            s.p_lodi = p_lodi;
            end

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(1,:) = s.var.gamma_star_eq(1,:);
                s.var.cv_star(1,:) = s.var.cv_star_eq(1,:);
            end

        case 'periodic'
            s.var.rho_u(1,:) = s.var.rho_u(end-1,:);
            s.var.rho_v(1,:) = s.var.rho_v(end-1,:);
            s.var.rho(1,:) = s.var.rho(end-1,:);
            s.var.rho_E(1,:) = s.var.rho_E(end-1,:);
            s.var.p(1,:) = s.var.p(end-1,:);
            s.var.T(1,:) = s.var.T(end-1,:);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(1,:) = s.var.gamma_star(end-1,:);
                s.var.cv_star(1,:) = s.var.cv_star(end-1,:);
            end

        case 'symmetry'
            s.var.rho(1,:) = s.var.rho(2,:);
            s.var.p(1,:) = s.var.p(2,:);

            % Remove wall-normal velocity component via reflection
            flux_vel = s.var.rho_u(2,2:end-1) .* s.mesh.lr_x_normal(1,:) + ...
                       s.var.rho_v(2,2:end-1) .* s.mesh.lr_y_normal(1,:);
            s.var.rho_u(1,2:end-1) = s.var.rho_u(2,2:end-1) - 2 * flux_vel .* s.mesh.lr_x_normal(1,:);
            s.var.rho_v(1,2:end-1) = s.var.rho_v(2,2:end-1) - 2 * flux_vel .* s.mesh.lr_y_normal(1,:);
            s.var.rho_u(1,1) = s.var.rho_u(1,2);
            s.var.rho_v(1,1) = s.var.rho_v(1,2);
            s.var.rho_u(1,end) = s.var.rho_u(1,end-1);
            s.var.rho_v(1,end) = s.var.rho_v(1,end-1);

            % Remove heat flux (adiabatic symmetry)
            s.var.rho_E(1,:) = s.var.rho_E(2,:);
            s.var.T(1,:) = s.var.T(2,:);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(1,:) = s.var.gamma_star(2,:);
                s.var.cv_star(1,:) = s.var.cv_star(2,:);
            end

        otherwise
            error('Non defined Boundary-condition, chi0')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Right boundary (chi1)
    switch s.boundary_conditions.boundary_chi1.name
        case 'inflow_subsonic'
            s.var.rho_u(end,:) = s.freestream.rho_u_0;
            s.var.rho_v(end,:) = s.freestream.rho_v_0;
            s.var.rho(end,:) = s.freestream.rho_0;
            s.var.rho_E(end,:) = 2 * s.var.rho_E(end-1,:) - s.var.rho_E(end-2,:);
            s.var.p(end,:) = 2 * s.var.p(end-1,:) - s.var.p(end-2,:);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(end,:) = s.var.gamma_star_eq(end,:);
                s.var.cv_star(end,:) = s.var.cv_star_eq(end,:);
            end

        case 'inflow_supersonic'
            s.var.rho_u(end,:) = s.freestream.rho_u_0;
            s.var.rho_v(end,:) = s.freestream.rho_v_0;
            s.var.rho(end,:) = s.freestream.rho_0;
            s.var.rho_E(end,:) = s.freestream.rho_E_0;
            s.var.p(end,:) = s.freestream.p_0;

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(end,:) = s.var.gamma_star_eq(end,:);
                s.var.cv_star(end,:) = s.var.cv_star_eq(end,:);
            end

        case 'no_slip_adiabatic'
            s.var.rho(end,:) = s.var.rho(end-1,:);
            s.var.rho_u(end,:) = -s.var.rho_u(end-1,:);
            s.var.rho_v(end,:) = -s.var.rho_v(end-1,:);
            s.var.rho_E(end,:) = 2 * s.var.rho_E(end-1,:) - s.var.rho_E(end-2,:);
            s.var.p(end,:) = 2 * s.var.p(end-1,:) - s.var.p(end-2,:);
            s.var.T(end,:) = s.var.T(end-1,:);  % No temperature gradient for adiabatic wall

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(end,:) = 2 * s.var.gamma_star(end-1,:) - s.var.gamma_star(end-2,:);
                s.var.cv_star(end,:) = 2 * s.var.cv_star(end-1,:) - s.var.cv_star(end-2,:);
            end

        case 'no_slip_isothermal'
            s.var.rho(end,:) = s.var.rho(end-1,:);
            s.var.rho_u(end,:) = -s.var.rho_u(end-1,:);
            s.var.rho_v(end,:) = -s.var.rho_v(end-1,:);
            s.var.rho_E(end,:) = 2 * s.var.rho_E(end-1,:) - s.var.rho_E(end-2,:);
            s.var.p(end,:) = 2 * s.var.p(end-1,:) - s.var.p(end-2,:);
            T_w = s.boundary_conditions.boundary_chi0.Tw * s.freestream.cv / s.freestream.energy_factor; % Non-dimensional wall temperature
            s.var.T(end,:) = 2 * T_w - s.var.T(end-1,:);  % Set wall temperature

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(end,:) = s.var.gamma_star(end-1,:);
                s.var.cv_star(end,:) = s.var.cv_star(end-1,:);
            end

        case 'outflow_supersonic'
            s.var.rho(end,:) = 2 * s.var.rho(end-1,:) - s.var.rho(end-2,:);
            s.var.rho_E(end,:) = 2 * s.var.rho_E(end-1,:) - s.var.rho_E(end-2,:);
            s.var.rho_u(end,:) = 2 * s.var.rho_u(end-1,:) - s.var.rho_u(end-2,:);
            s.var.rho_v(end,:) = 2 * s.var.rho_v(end-1,:) - s.var.rho_v(end-2,:);
            s.var.p(end,:) = 2 * s.var.p(end-1,:) - s.var.p(end-2,:);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(end,:) = 2 * s.var.gamma_star(end-1,:) - s.var.gamma_star(end-2,:);
                s.var.cv_star(end,:) = 2 * s.var.cv_star(end-1,:) - s.var.cv_star(end-2,:);
            end

        case 'outflow_supersonic_1st'
            s.var.rho_u(end,:) = s.var.rho_u(end-1,:);
            s.var.rho_v(end,:) = s.var.rho_v(end-1,:);
            s.var.rho(end,:) = s.var.rho(end-1,:);
            s.var.rho_E(end,:) = s.var.rho_E(end-1,:);
            s.var.p(end,:) = s.var.p(end-1,:);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(end,:) = s.var.gamma_star(end-1,:);
                s.var.cv_star(end,:) = s.var.cv_star(end-1,:);
            end

        case 'outflow_subsonic'
            s.var.rho_u(end,:) = 2 * s.var.rho_u(end-1,:) - s.var.rho_u(end-2,:);
            s.var.rho_v(end,:) = 2 * s.var.rho_v(end-1,:) - s.var.rho_v(end-2,:);
            s.var.rho(end,:) = 2 * s.var.rho(end-1,:) - s.var.rho(end-2,:);
            s.var.rho_E(end,:) = 2 * s.var.rho_E(end-1,:) - s.var.rho_E(end-2,:);
            s.var.p(end,:) = s.var.p(end-1,:); % zero pressure gradient

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(end,:) = s.var.gamma_star(end-1,:);
                s.var.cv_star(end,:) = s.var.cv_star(end-1,:);
            end

        case 'outflow_NRCBC'
            % LODI boundary conditions (Poinsot & Lele): L1 = 0
            % Extrapolate all variables (supersonic default)
            s.var.rho(end,:) = 2 * s.var.rho(end-1,:) - s.var.rho(end-2,:);
            s.var.rho_u(end,:) = 2 * s.var.rho_u(end-1,:) - s.var.rho_u(end-2,:);
            s.var.rho_v(end,:) = 2 * s.var.rho_v(end-1,:) - s.var.rho_v(end-2,:);
            s.var.rho_E(end,:) = 2 * s.var.rho_E(end-1,:) - s.var.rho_E(end-2,:);
            s.var.p(end,:) = 2 * s.var.p(end-1,:) - s.var.p(end-2,:);

            % Sound speed and density at last interior cell
            a = s.var.a(end-1,2:end-1);
            rho_int = s.var.rho(end-1,2:end-1);

            % Normal velocity at ghost cell, positive exiting domain
            u_n_ghost = (s.var.rho_u(end,2:end-1) .* s.mesh.lr_x_normal(end,:) + ...
                        s.var.rho_v(end,2:end-1) .* s.mesh.lr_y_normal(end,:)) ./ s.var.rho(end,2:end-1);

            % Normal velocity at last interior cell
            u_n_int = (s.var.rho_u(end-1,2:end-1) .* s.mesh.lr_x_normal(end,:) + ...
                       s.var.rho_v(end-1,2:end-1) .* s.mesh.lr_y_normal(end,:)) ./ rho_int;

            % LODI correction for subsonic cells (L1 = 0 => dp/dx_n = rho*c*du_n/dx_n)
            subsonic = u_n_ghost < a;
            if any(subsonic)
                p_int = s.var.p(end-1,2:end-1);
                p_lodi = p_int + rho_int .* a .* (u_n_ghost - u_n_int);
                s.var.p(end,2:end-1) = subsonic .* p_lodi + (~subsonic) .* s.var.p(end,2:end-1);
                rho_e = s.var.p(end,2:end-1) ./ (s.var.gamma_star(end,2:end-1) - 1);
                rho_E_lodi = rho_e + 0.5 * (s.var.rho_u(end,2:end-1).^2 + s.var.rho_v(end,2:end-1).^2) ./ s.var.rho(end,2:end-1);
                s.var.rho_E(end,2:end-1) = subsonic .* rho_E_lodi + (~subsonic) .* s.var.rho_E(end,2:end-1);
            end

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(end,:) = s.var.gamma_star_eq(end,:);
                s.var.cv_star(end,:) = s.var.cv_star_eq(end,:);
            end

        case 'periodic'
            s.var.rho_u(end,:) = s.var.rho_u(2,:);
            s.var.rho_v(end,:) = s.var.rho_v(2,:);
            s.var.rho(end,:) = s.var.rho(2,:);
            s.var.rho_E(end,:) = s.var.rho_E(2,:);
            s.var.p(end,:) = s.var.p(2,:);
            s.var.T(end,:) = s.var.T(2,:);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(end,:) = s.var.gamma_star(2,:);
                s.var.cv_star(end,:) = s.var.cv_star(2,:);
            end

        case 'symmetry'
            s.var.rho(end,:) = s.var.rho(end-1,:);
            s.var.p(end,:) = s.var.p(end-1,:);

            % Remove wall-normal velocity component via reflection
            flux_vel = s.var.rho_u(end-1,2:end-1) .* s.mesh.lr_x_normal(end,:) + ...
                       s.var.rho_v(end-1,2:end-1) .* s.mesh.lr_y_normal(end,:);
            s.var.rho_u(end,2:end-1) = s.var.rho_u(end-1,2:end-1) - 2 * flux_vel .* s.mesh.lr_x_normal(end,:);
            s.var.rho_v(end,2:end-1) = s.var.rho_v(end-1,2:end-1) - 2 * flux_vel .* s.mesh.lr_y_normal(end,:);
            s.var.rho_u(end,1) = s.var.rho_u(end,2);
            s.var.rho_v(end,1) = s.var.rho_v(end,2);
            s.var.rho_u(end,end) = s.var.rho_u(end,end-1);
            s.var.rho_v(end,end) = s.var.rho_v(end,end-1);

            % Remove heat flux (adiabatic symmetry)
            s.var.rho_E(end,:) = s.var.rho_E(end-1,:);
            s.var.T(end,:) = s.var.T(end-1,:);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(end,:) = s.var.gamma_star(end-1,:);
                s.var.cv_star(end,:) = s.var.cv_star(end-1,:);
            end
        otherwise
            error('Non defined Boundary-condition, chi0')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Bottom boundary (eta0)
    switch s.boundary_conditions.boundary_eta0.name
        case 'symmetry'
            s.var.p(:,1) = s.var.p(:,2);
            s.var.rho(:,1) = s.var.rho(:,2);

            % Remove wall-normal velocity component via reflection
            flux_vel = s.var.rho_u(2:end-1,2) .* s.mesh.bt_x_normal(:,1) + ...
                       s.var.rho_v(2:end-1,2) .* s.mesh.bt_y_normal(:,1);
            s.var.rho_u(2:end-1,1) = s.var.rho_u(2:end-1,2) - 2 * flux_vel .* s.mesh.bt_x_normal(:,1);
            s.var.rho_v(2:end-1,1) = s.var.rho_v(2:end-1,2) - 2 * flux_vel .* s.mesh.bt_y_normal(:,1);
            s.var.rho_u(1,1) = s.var.rho_u(2,1);
            s.var.rho_v(1,1) = s.var.rho_v(2,1);
            s.var.rho_u(end,1) = s.var.rho_u(end-1,1);
            s.var.rho_v(end,1) = s.var.rho_v(end-1,1);

            % Remove heat flux (adiabatic symmetry)
            s.var.rho_E(:,1) = s.var.rho_E(:,2);
            s.var.T(:,1) = s.var.T(:,2);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,1) = s.var.gamma_star(:,2);
                s.var.cv_star(:,1) = s.var.cv_star(:,2);
            end

        case 'no_slip_adiabatic'
            s.var.rho(:,1) = s.var.rho(:,2);
            s.var.rho_u(:,1) = -s.var.rho_u(:,2);
            s.var.rho_v(:,1) = -s.var.rho_v(:,2);
            s.var.rho_E(:,1) = 2 * s.var.rho_E(:,2) - s.var.rho_E(:,3);
            s.var.p(:,1) = 2 * s.var.p(:,2) - s.var.p(:,3);
            s.var.T(:,1) = s.var.T(:,2);  % No temperature gradient for adiabatic wall

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,1) = 2 * s.var.gamma_star(:,2) - s.var.gamma_star(:,3);
                s.var.cv_star(:,1) = 2 * s.var.cv_star(:,2) - s.var.cv_star(:,3);
            end

        case 'no_slip_isothermal'
            s.var.rho(:,1) = s.var.rho(:,2);
            s.var.rho_u(:,1) = -s.var.rho_u(:,2);
            s.var.rho_v(:,1) = -s.var.rho_v(:,2);
            s.var.rho_E(:,1) = 2 * s.var.rho_E(:,2) - s.var.rho_E(:,3);
            s.var.p(:,1) = 2 * s.var.p(:,2) - s.var.p(:,3);
            T_w = s.boundary_conditions.boundary_eta0.Tw * s.freestream.cv / s.freestream.energy_factor; % Non-dimensional wall temperature
            s.var.T(:,1) = 2 * T_w - s.var.T(:,2);  % Set wall temperature

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,1) = s.var.gamma_star(:,2);
                s.var.cv_star(:,1) = s.var.cv_star(:,2);
            end

        case 'outflow_supersonic'
            s.var.rho(:,1) = 2 * s.var.rho(:,2) - s.var.rho(:,3);
            s.var.rho_E(:,1) = 2 * s.var.rho_E(:,2) - s.var.rho_E(:,3);
            s.var.rho_u(:,1) = 2 * s.var.rho_u(:,2) - s.var.rho_u(:,3);
            s.var.rho_v(:,1) = 2 * s.var.rho_v(:,2) - s.var.rho_v(:,3);
            s.var.p(:,1) = 2 * s.var.p(:,2) - s.var.p(:,3);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,1) = 2 * s.var.gamma_star(:,2) - s.var.gamma_star(:,3);
                s.var.cv_star(:,1) = 2 * s.var.cv_star(:,2) - s.var.cv_star(:,3);
            end

        case 'outflow_supersonic_1st'
            s.var.rho(:,1) = s.var.rho(:,2);
            s.var.rho_u(:,1) = s.var.rho_u(:,2);
            s.var.rho_v(:,1) = s.var.rho_v(:,2);
            s.var.rho_E(:,1) = s.var.rho_E(:,2);
            s.var.p(:,1) = s.var.p(:,2);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,1) = s.var.gamma_star(:,2);
                s.var.cv_star(:,1) = s.var.cv_star(:,2);
            end

        case 'outflow_subsonic'
            s.var.rho(:,1) = 2 * s.var.rho(:,2) - s.var.rho(:,3);
            s.var.rho_E(:,1) = 2 * s.var.rho_E(:,2) - s.var.rho_E(:,3);
            s.var.rho_u(:,1) = 2 * s.var.rho_u(:,2) - s.var.rho_u(:,3);
            s.var.rho_v(:,1) = 2 * s.var.rho_v(:,2) - s.var.rho_v(:,3);
            s.var.p(:,1) = s.var.p(:,2); % zero pressure gradient

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,1) = s.var.gamma_star(:,2);
                s.var.cv_star(:,1) = s.var.cv_star(:,2);
            end

        case 'outflow_NRCBC'
            % LODI boundary conditions (Poinsot & Lele): L1 = 0
            % Extrapolate all variables (supersonic default)
            s.var.rho(:,1) = 2 * s.var.rho(:,2) - s.var.rho(:,3);
            s.var.rho_u(:,1) = 2 * s.var.rho_u(:,2) - s.var.rho_u(:,3);
            s.var.rho_v(:,1) = 2 * s.var.rho_v(:,2) - s.var.rho_v(:,3);
            s.var.rho_E(:,1) = 2 * s.var.rho_E(:,2) - s.var.rho_E(:,3);
            s.var.p(:,1) = 2 * s.var.p(:,2) - s.var.p(:,3);

            % Sound speed and density at first interior cell
            a = s.var.a(2:end-1,2);
            rho_int = s.var.rho(2:end-1,2);

            % Normal velocity at ghost cell, positive exiting domain
            u_n_ghost = -(s.var.rho_u(2:end-1,1) .* s.mesh.bt_x_normal(:,1) + ...
                          s.var.rho_v(2:end-1,1) .* s.mesh.bt_y_normal(:,1)) ./ s.var.rho(2:end-1,1);

            % Normal velocity at first interior cell
            u_n_int = -(s.var.rho_u(2:end-1,2) .* s.mesh.bt_x_normal(:,1) + ...
                        s.var.rho_v(2:end-1,2) .* s.mesh.bt_y_normal(:,1)) ./ rho_int;

            % LODI correction for subsonic cells (L1 = 0 => dp/dx_n = rho*c*du_n/dx_n)
            subsonic = u_n_ghost < a;
            if any(subsonic)
                p_int = s.var.p(2:end-1,2);
                p_lodi = p_int + rho_int .* a .* (u_n_ghost - u_n_int);
                % Apply LODI pressure only to subsonic cells, keep extrapolated for supersonic
                s.var.p(2:end-1,1) = subsonic .* p_lodi + (~subsonic) .* s.var.p(2:end-1,1);
                % Recompute total energy from corrected pressure
                rho_e = s.var.p(2:end-1,1) ./ (s.var.gamma_star(2:end-1,1) - 1);
                rho_E_lodi = rho_e + 0.5 * (s.var.rho_u(2:end-1,1).^2 + s.var.rho_v(2:end-1,1).^2) ./ s.var.rho(2:end-1,1);
                s.var.rho_E(2:end-1,1) = subsonic .* rho_E_lodi + (~subsonic) .* s.var.rho_E(2:end-1,1);
            end

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,1) = s.var.gamma_star_eq(:,1);
                s.var.cv_star(:,1) = s.var.cv_star_eq(:,1);
            end

        case 'inflow_supersonic'
            s.var.rho(:,1) = s.freestream.rho_0;
            s.var.rho_u(:,1) = s.freestream.rho_u_0;
            s.var.rho_v(:,1) = s.freestream.rho_v_0;
            s.var.rho_E(:,1) = s.freestream.rho_E_0;
            s.var.p(:,1) = s.freestream.p_0;

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,1) = s.var.gamma_star_eq(:,1);
                s.var.cv_star(:,1) = s.var.cv_star_eq(:,1);
            end

        case 'inflow_subsonic'
            s.var.rho(:,1) = s.freestream.rho_0;
            s.var.rho_u(:,1) = s.freestream.rho_u_0;
            s.var.rho_v(:,1) = s.freestream.rho_v_0;
            s.var.rho_E(:,1) = 2 * s.var.rho_E(:,2) - s.var.rho_E(:,3);
            s.var.p(:,1) = 2 * s.var.p(:,2) - s.var.p(:,3); % Release acoustic characteristic for subsonic inflow

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,1) = s.var.gamma_star_eq(:,1);
                s.var.cv_star(:,1) = s.var.cv_star_eq(:,1);
            end

        case 'periodic'
            s.var.rho(:,1) = s.var.rho(:,end-1);
            s.var.rho_u(:,1) = s.var.rho_u(:,end-1);
            s.var.rho_v(:,1) = s.var.rho_v(:,end-1);
            s.var.rho_E(:,1) = s.var.rho_E(:,end-1);
            s.var.p(:,1) = s.var.p(:,end-1);
            s.var.T(:,1) = s.var.T(:,end-1);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,1) = s.var.gamma_star(:,end-1);
                s.var.cv_star(:,1) = s.var.cv_star(:,end-1);
            end

        otherwise
            error('Non defined Boundary-condition, eta0')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Top boundary (eta1)
    switch s.boundary_conditions.boundary_eta1.name
        case 'symmetry'
            s.var.p(:,end) = s.var.p(:,end-1);
            s.var.rho(:,end) = s.var.rho(:,end-1);

            % Remove wall-normal velocity component
            flux_vel = s.var.rho_u(2:end-1,end-1) .* s.mesh.bt_x_normal(:,end) + ...
                       s.var.rho_v(2:end-1,end-1) .* s.mesh.bt_y_normal(:,end);
            s.var.rho_u(2:end-1,end) = s.var.rho_u(2:end-1,end-1) - 2 * flux_vel .* s.mesh.bt_x_normal(:,end);
            s.var.rho_v(2:end-1,end) = s.var.rho_v(2:end-1,end-1) - 2 * flux_vel .* s.mesh.bt_y_normal(:,end);
            s.var.rho_u(1,end) = s.var.rho_u(2,end);
            s.var.rho_v(1,end) = s.var.rho_v(2,end);
            s.var.rho_u(end,end) = s.var.rho_u(end-1,end);
            s.var.rho_v(end,end) = s.var.rho_v(end-1,end);

            % Remove heat flux (adiabatic symmetry)
            s.var.rho_E(:,end) = s.var.rho_E(:,end-1);
            s.var.T(:,end) = s.var.T(:,end-1);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,end) = s.var.gamma_star(:,end-1);
                s.var.cv_star(:,end) = s.var.cv_star(:,end-1);
            end

        case 'no_slip_adiabatic'
            s.var.rho(:,end) = s.var.rho(:,end-1);
            s.var.rho_u(:,end) = -s.var.rho_u(:,end-1);
            s.var.rho_v(:,end) = -s.var.rho_v(:,end-1);
            s.var.rho_E(:,end) = 2 * s.var.rho_E(:,end-1) - s.var.rho_E(:,end-2);
            s.var.p(:,end) = 2 * s.var.p(:,end-1) - s.var.p(:,end-2);
            s.var.T(:,end) = s.var.T(:,end-1); % No temperature gradient for adiabatic wall

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,end) = s.var.gamma_star(:,end-1);
                s.var.cv_star(:,end) = s.var.cv_star(:,end-1);
            end

        case 'no_slip_isothermal'
            s.var.rho(:,end) = s.var.rho(:,end-1);
            s.var.rho_u(:,end) = -s.var.rho_u(:,end-1);
            s.var.rho_v(:,end) = -s.var.rho_v(:,end-1);
            s.var.rho_E(:,end) = 2 * s.var.rho_E(:,end-1) - s.var.rho_E(:,end-2);
            s.var.p(:,end) = 2 * s.var.p(:,end-1) - s.var.p(:,end-2);
            T_w = s.boundary_conditions.boundary_eta1.Tw * s.freestream.cv / s.freestream.energy_factor; % Non-dimensional wall temperature
            s.var.T(:,end) = 2 * T_w - s.var.T(:,end-1); % Set wall temperature

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,end) = s.var.gamma_star(:,end-1);
                s.var.cv_star(:,end) = s.var.cv_star(:,end-1);
            end

        case 'outflow_supersonic'
            s.var.rho(:,end) = 2 * s.var.rho(:,end-1) - s.var.rho(:,end-2);
            s.var.rho_E(:,end) = 2 * s.var.rho_E(:,end-1) - s.var.rho_E(:,end-2);
            s.var.rho_u(:,end) = 2 * s.var.rho_u(:,end-1) - s.var.rho_u(:,end-2);
            s.var.rho_v(:,end) = 2 * s.var.rho_v(:,end-1) - s.var.rho_v(:,end-2);
            s.var.p(:,end) = 2 * s.var.p(:,end-1) - s.var.p(:,end-2);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,end) = 2 * s.var.gamma_star(:,end-1) - s.var.gamma_star(:,end-2);
                s.var.cv_star(:,end) = 2 * s.var.cv_star(:,end-1) - s.var.cv_star(:,end-2);
            end

        case 'outflow_supersonic_1st'
            s.var.rho(:,end) = s.var.rho(:,end-1);
            s.var.rho_u(:,end) = s.var.rho_u(:,end-1);
            s.var.rho_v(:,end) = s.var.rho_v(:,end-1);
            s.var.rho_E(:,end) = s.var.rho_E(:,end-1);
            s.var.p(:,end) = s.var.p(:,end-1);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,end) = s.var.gamma_star(:,end-1);
                s.var.cv_star(:,end) = s.var.cv_star(:,end-1);
            end

        case 'outflow_subsonic'
            s.var.rho(:,end) = 2 * s.var.rho(:,end-1) - s.var.rho(:,end-2);
            s.var.rho_E(:,end) = 2 * s.var.rho_E(:,end-1) - s.var.rho_E(:,end-2);
            s.var.rho_u(:,end) = 2 * s.var.rho_u(:,end-1) - s.var.rho_u(:,end-2);
            s.var.rho_v(:,end) = 2 * s.var.rho_v(:,end-1) - s.var.rho_v(:,end-2);
            s.var.p(:,end) = s.var.p(:,end-1); % zero pressure gradient

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,end) = s.var.gamma_star(:,end-1);
                s.var.cv_star(:,end) = s.var.cv_star(:,end-1);
            end

        case 'outflow_NRCBC'
            % LODI boundary conditions (Poinsot & Lele): L1 = 0
            % Extrapolate all variables (supersonic default)
            s.var.rho(:,end) = 2 * s.var.rho(:,end-1) - s.var.rho(:,end-2);
            s.var.rho_u(:,end) = 2 * s.var.rho_u(:,end-1) - s.var.rho_u(:,end-2);
            s.var.rho_v(:,end) = 2 * s.var.rho_v(:,end-1) - s.var.rho_v(:,end-2);
            s.var.rho_E(:,end) = 2 * s.var.rho_E(:,end-1) - s.var.rho_E(:,end-2);
            s.var.p(:,end) = 2 * s.var.p(:,end-1) - s.var.p(:,end-2);

            % Sound speed and density at last interior cell
            a = s.var.a(2:end-1,end-1);
            rho_int = s.var.rho(2:end-1,end-1);

            % Normal velocity at ghost cell, positive exiting domain (bt_normal points +eta = outward at eta1)
            u_n_ghost = (s.var.rho_u(2:end-1,end) .* s.mesh.bt_x_normal(:,end) + ...
                         s.var.rho_v(2:end-1,end) .* s.mesh.bt_y_normal(:,end)) ./ s.var.rho(2:end-1,end);

            % Normal velocity at last interior cell
            u_n_int = (s.var.rho_u(2:end-1,end-1) .* s.mesh.bt_x_normal(:,end) + ...
                       s.var.rho_v(2:end-1,end-1) .* s.mesh.bt_y_normal(:,end)) ./ rho_int;

            % LODI correction for subsonic cells (L1 = 0 => dp/dx_n = rho*c*du_n/dx_n)
            subsonic = u_n_ghost < a;
            if any(subsonic)
                p_int = s.var.p(2:end-1,end-1);
                p_lodi = p_int + rho_int .* a .* (u_n_ghost - u_n_int);
                s.var.p(2:end-1,end) = subsonic .* p_lodi + (~subsonic) .* s.var.p(2:end-1,end);
                rho_e = s.var.p(2:end-1,end) ./ (s.var.gamma_star(2:end-1,end) - 1);
                rho_E_lodi = rho_e + 0.5 * (s.var.rho_u(2:end-1,end).^2 + s.var.rho_v(2:end-1,end).^2) ./ s.var.rho(2:end-1,end);
                s.var.rho_E(2:end-1,end) = subsonic .* rho_E_lodi + (~subsonic) .* s.var.rho_E(2:end-1,end);
            end

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,end) = s.var.gamma_star_eq(:,end);
                s.var.cv_star(:,end) = s.var.cv_star_eq(:,end);
            end

        case 'inflow_supersonic'
            s.var.rho_u(:,end) = s.freestream.rho_u_0;
            s.var.rho_v(:,end) = s.freestream.rho_v_0;
            s.var.rho(:,end) = s.freestream.rho_0;
            s.var.rho_E(:,end) = s.freestream.rho_E_0;
            s.var.p(:,end) = s.freestream.p_0;

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,end) = s.var.gamma_star_eq(:,end);
                s.var.cv_star(:,end) = s.var.cv_star_eq(:,end);
            end

        case 'inflow_subsonic'
            s.var.rho_u(:,end) = s.freestream.rho_u_0;
            s.var.rho_v(:,end) = s.freestream.rho_v_0;
            s.var.rho(:,end) = s.freestream.rho_0;
            s.var.rho_E(:,end) = 2 * s.var.rho_E(:,end-1) - s.var.rho_E(:,end-2);
            s.var.p(:,end) = 2 * s.var.p(:,end-1) - s.var.p(:,end-2);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,end) = s.var.gamma_star_eq(:,end);
                s.var.cv_star(:,end) = s.var.cv_star_eq(:,end);
            end

        case 'periodic'
            s.var.rho(:,end) = s.var.rho(:,2);
            s.var.rho_u(:,end) = s.var.rho_u(:,2);
            s.var.rho_v(:,end) = s.var.rho_v(:,2);
            s.var.rho_E(:,end) = s.var.rho_E(:,2);
            s.var.p(:,end) = s.var.p(:,2);
            s.var.T(:,end) = s.var.T(:,2);

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,end) = s.var.gamma_star(:,2);
                s.var.cv_star(:,end) = s.var.cv_star(:,2);
            end

        case 'shock' % Only for eta1
            % Shock update handled within UPDATE_SHOCK_BC

        case 'lid_driven_cavity' % Only for eta1
            % Zero pressure gradient at lid
            s.var.p(:,end) = s.var.p(:,end-1);
            e2 = (s.var.rho_E(:,end-1) - (s.var.rho_u(:,end-1).^2 + s.var.rho_v(:,end-1).^2) ./ ...
                (s.var.rho(:,end-1)) / 2) ./ s.var.rho(:,end-1);
            s.var.rho(:,end) = s.var.rho(:,end-1) .* e2 / s.e_ref;

            % Impose lid velocity profile (quartic bump)
            u_tan = 16 * (s.mesh.x_Ext(:,end).^2) .* (1 - s.mesh.x_Ext(:,end)).^2;
            u_ghost = 2 * u_tan - s.var.rho_u(:,end-1) ./ s.var.rho(:,end-1);
            s.var.rho_v(:,end) = -s.var.rho_v(:,end-1);
            s.var.rho_u(:,end) = u_ghost .* s.var.rho(:,end);

            % Isothermal wall condition
            eghost = 2 * s.e_ref - e2;
            s.var.rho_E(:,end) = eghost .* s.var.rho(:,end) + ...
                (s.var.rho_u(:,end).^2 + s.var.rho_v(:,end).^2) ./ (s.var.rho(:,end)) / 2;

            if ~s.chemistry.chemical_equilibrium
                s.var.gamma_star(:,end) = 2 * s.var.gamma_star(:,end-1) - s.var.gamma_star(:,end-2);
                s.var.cv_star(:,end) = 2 * s.var.cv_star(:,end-1) - s.var.cv_star(:,end-2);
            end

        otherwise
            error('Non defined Boundary-condition, eta1')
    end
end
