function s = FLUX_TOTAL_ENERGY_3D(s)
% FLUX_TOTAL_ENERGY_3D - Compute the inviscid total energy flux for the 3D axisymmetric Euler equations.
%
%   Evaluates the net total energy flux through the east, west, south, and
%   north faces of each interior cell using midpoint interpolation and
%   midpoint quadrature, both second-order accurate on rectangular cells.
%   This is the axisymmetric (3D) variant. The front/back faces do not
%   contribute because the azimuthal velocity is zero. Results are
%   normalised by cell volume.
%
% Syntax:
%   s = FLUX_TOTAL_ENERGY_3D(s)
%
% Inputs:
%   s - Structure containing at minimum:
%                .var.rho                     - Density field (including ghosts).
%                .var.rho_u, .var.rho_v      - Momentum components (including ghosts).
%                .var.rho_E                  - Total energy field (including ghosts).
%                .var.p                       - Pressure field (including ghosts).
%                .mesh.lr_x_normal, .mesh.lr_y_normal - Left/right face unit normals.
%                .mesh.bt_x_normal, .mesh.bt_y_normal - Bottom/top face unit normals.
%                .mesh.lr_area, .mesh.bt_area - Face areas.
%                .mesh.volume                - Cell volumes.
%
% Outputs:
%   s - Updated structure with:
%                .flux.rho_E - Cell-averaged total energy flux (interior cells).
%
% Notes:
%   - The energy flux is (rho_E + p) * V_n, i.e. the enthalpy flux.
%   - Front/back faces do not contribute since azimuthal velocity is zero.
%   - Sign convention: outgoing flux is negative, incoming flux is positive.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% East and west faces
    % Midpoint interpolation (2nd order on rectangles)
    temp_rho_E_face = (s.var.rho_E(1:end-1,2:end-1) + s.var.rho_E(2:end,2:end-1)) / 2;
    temp_P_face = (s.var.p(1:end-1,2:end-1) + s.var.p(2:end,2:end-1)) / 2;
    temp_vel_face = (s.var.rho_u(1:end-1,2:end-1) + s.var.rho_u(2:end,2:end-1)) ./ ...
                    (s.var.rho(2:end,2:end-1) + s.var.rho(1:end-1,2:end-1)) .* s.mesh.lr_x_normal;
    temp_vel_face = temp_vel_face + (s.var.rho_v(1:end-1,2:end-1) + s.var.rho_v(2:end,2:end-1)) ./ ...
                    (s.var.rho(2:end,2:end-1) + s.var.rho(1:end-1,2:end-1)) .* s.mesh.lr_y_normal;

    % Midpoint quadrature (2nd order on rectangles)
    s.flux.rho_E = -((temp_rho_E_face(2:end,:) + temp_P_face(2:end,:)) .* temp_vel_face(2:end,:)) .* s.mesh.lr_area(2:end,:);            % Outgoing fluxes
    s.flux.rho_E = s.flux.rho_E + ((temp_rho_E_face(1:end-1,:) + temp_P_face(1:end-1,:)) .* temp_vel_face(1:end-1,:)) .* s.mesh.lr_area(1:end-1,:);  % Incoming fluxes

    %% South and north faces
    % Midpoint interpolation (2nd order on rectangles)
    temp_rho_E_face2 = (s.var.rho_E(2:end-1,1:end-1) + s.var.rho_E(2:end-1,2:end)) / 2;
    temp_P_face2 = (s.var.p(2:end-1,1:end-1) + s.var.p(2:end-1,2:end)) / 2;
    temp_vel_face2 = (s.var.rho_u(2:end-1,1:end-1) + s.var.rho_u(2:end-1,2:end)) ./ ...
                     (s.var.rho(2:end-1,2:end) + s.var.rho(2:end-1,1:end-1)) .* s.mesh.bt_x_normal;
    temp_vel_face2 = temp_vel_face2 + (s.var.rho_v(2:end-1,1:end-1) + s.var.rho_v(2:end-1,2:end)) ./ ...
                     (s.var.rho(2:end-1,2:end) + s.var.rho(2:end-1,1:end-1)) .* s.mesh.bt_y_normal;

    % Midpoint quadrature (2nd order on rectangles)
    s.flux.rho_E = s.flux.rho_E - ((temp_rho_E_face2(:,2:end) + temp_P_face2(:,2:end)) .* temp_vel_face2(:,2:end)) .* s.mesh.bt_area(:,2:end);            % Outgoing fluxes
    s.flux.rho_E = s.flux.rho_E + ((temp_rho_E_face2(:,1:end-1) + temp_P_face2(:,1:end-1)) .* temp_vel_face2(:,1:end-1)) .* s.mesh.bt_area(:,1:end-1);  % Incoming fluxes

    %% Front and back faces (3D axisymmetric)
    % No contribution: azimuthal velocity is zero

    %% Normalise by cell volume
    s.flux.rho_E = s.flux.rho_E ./ s.mesh.volume;
end
