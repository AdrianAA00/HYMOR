function s = FLUX_MOMENTUM(s)
% FLUX_MOMENTUM - Compute inviscid momentum fluxes for the 2D Euler equations.
%
%   Evaluates the net x- and y-momentum fluxes through all four cell faces
%   (east, west, south, north) using midpoint interpolation and midpoint
%   quadrature, both second-order accurate on rectangular cells. The
%   momentum flux includes both the convective term (rho*u_i * V_n) and the
%   pressure contribution (p * n_i). Results are normalised by cell volume.
%
% Syntax:
%   s = FLUX_MOMENTUM(s)
%
% Inputs:
%   s - Structure containing at minimum:
%                .var.rho                     - Density field (including ghosts).
%                .var.rho_u, .var.rho_v      - Momentum components (including ghosts).
%                .var.p                       - Pressure field (including ghosts).
%                .mesh.lr_x_normal, .mesh.lr_y_normal - Left/right face unit normals.
%                .mesh.bt_x_normal, .mesh.bt_y_normal - Bottom/top face unit normals.
%                .mesh.lr_area, .mesh.bt_area - Face areas.
%                .mesh.volume                - Cell volumes.
%
% Outputs:
%   s - Updated structure with:
%                .flux.rho_u - Cell-averaged x-momentum flux (interior cells).
%                .flux.rho_v - Cell-averaged y-momentum flux (interior cells).
%
% Notes:
%   - Face velocity is computed as (rho_u / rho) projected onto the face
%     normal, using the sum of neighbour momentum and density to avoid an
%     extra division.
%   - Sign convention: outgoing flux is negative, incoming flux is positive.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% X-momentum (u-component)

    % East and west faces - midpoint interpolation (2nd order on rectangles)
    temp_rho_u_face = (s.var.rho_u(1:end-1,2:end-1) + s.var.rho_u(2:end,2:end-1)) / 2;
    temp_P_face = (s.var.p(1:end-1,2:end-1) + s.var.p(2:end,2:end-1)) .* s.mesh.lr_x_normal / 2;
    temp_vel_face = (s.var.rho_u(1:end-1,2:end-1) + s.var.rho_u(2:end,2:end-1)) ./ ...
                    (s.var.rho(2:end,2:end-1) + s.var.rho(1:end-1,2:end-1)) .* s.mesh.lr_x_normal;
    temp_vel_face = temp_vel_face + (s.var.rho_v(1:end-1,2:end-1) + s.var.rho_v(2:end,2:end-1)) ./ ...
                    (s.var.rho(2:end,2:end-1) + s.var.rho(1:end-1,2:end-1)) .* s.mesh.lr_y_normal;

    % East and west faces - midpoint quadrature (2nd order on rectangles)
    s.flux.rho_u = -(temp_rho_u_face(2:end,:) .* temp_vel_face(2:end,:) + temp_P_face(2:end,:)) .* s.mesh.lr_area(2:end,:);            % Outgoing fluxes
    s.flux.rho_u = s.flux.rho_u + (temp_rho_u_face(1:end-1,:) .* temp_vel_face(1:end-1,:) + temp_P_face(1:end-1,:)) .* s.mesh.lr_area(1:end-1,:);  % Incoming fluxes

    % South and north faces - midpoint interpolation (2nd order on rectangles)
    temp_rho_u_face2 = (s.var.rho_u(2:end-1,1:end-1) + s.var.rho_u(2:end-1,2:end)) / 2;
    temp_P_face2 = (s.var.p(2:end-1,1:end-1) + s.var.p(2:end-1,2:end)) .* s.mesh.bt_x_normal / 2;
    temp_vel_face2 = (s.var.rho_u(2:end-1,1:end-1) + s.var.rho_u(2:end-1,2:end)) ./ ...
                     (s.var.rho(2:end-1,2:end) + s.var.rho(2:end-1,1:end-1)) .* s.mesh.bt_x_normal;
    temp_vel_face2 = temp_vel_face2 + (s.var.rho_v(2:end-1,1:end-1) + s.var.rho_v(2:end-1,2:end)) ./ ...
                     (s.var.rho(2:end-1,2:end) + s.var.rho(2:end-1,1:end-1)) .* s.mesh.bt_y_normal;

    % South and north faces - midpoint quadrature (2nd order on rectangles)
    s.flux.rho_u = s.flux.rho_u - (temp_rho_u_face2(:,2:end) .* temp_vel_face2(:,2:end) + temp_P_face2(:,2:end)) .* s.mesh.bt_area(:,2:end);            % Outgoing fluxes
    s.flux.rho_u = s.flux.rho_u + (temp_rho_u_face2(:,1:end-1) .* temp_vel_face2(:,1:end-1) + temp_P_face2(:,1:end-1)) .* s.mesh.bt_area(:,1:end-1);  % Incoming fluxes

    % Normalise by cell volume
    s.flux.rho_u = s.flux.rho_u ./ s.mesh.volume;

    %% Y-momentum (v-component)

    % East and west faces - midpoint interpolation (2nd order on rectangles)
    temp_rho_v_face = (s.var.rho_v(1:end-1,2:end-1) + s.var.rho_v(2:end,2:end-1)) / 2;
    temp_P_face = (s.var.p(1:end-1,2:end-1) + s.var.p(2:end,2:end-1)) .* s.mesh.lr_y_normal / 2;

    % East and west faces - midpoint quadrature (2nd order on rectangles)
    s.flux.rho_v = -(temp_rho_v_face(2:end,:) .* temp_vel_face(2:end,:) + temp_P_face(2:end,:)) .* s.mesh.lr_area(2:end,:);            % Outgoing fluxes
    s.flux.rho_v = s.flux.rho_v + (temp_rho_v_face(1:end-1,:) .* temp_vel_face(1:end-1,:) + temp_P_face(1:end-1,:)) .* s.mesh.lr_area(1:end-1,:);  % Incoming fluxes

    % South and north faces - midpoint interpolation (2nd order on rectangles)
    temp_rho_v_face2 = (s.var.rho_v(2:end-1,1:end-1) + s.var.rho_v(2:end-1,2:end)) / 2;
    temp_P_face2 = (s.var.p(2:end-1,1:end-1) + s.var.p(2:end-1,2:end)) .* s.mesh.bt_y_normal / 2;

    % South and north faces - midpoint quadrature (2nd order on rectangles)
    s.flux.rho_v = s.flux.rho_v - (temp_rho_v_face2(:,2:end) .* temp_vel_face2(:,2:end) + temp_P_face2(:,2:end)) .* s.mesh.bt_area(:,2:end);            % Outgoing fluxes
    s.flux.rho_v = s.flux.rho_v + (temp_rho_v_face2(:,1:end-1) .* temp_vel_face2(:,1:end-1) + temp_P_face2(:,1:end-1)) .* s.mesh.bt_area(:,1:end-1);  % Incoming fluxes

    % Normalise by cell volume
    s.flux.rho_v = s.flux.rho_v ./ s.mesh.volume;
end
