function s = FLUX_MOMENTUM_3D(s)
% FLUX_MOMENTUM_3D - Compute inviscid momentum fluxes for the 3D axisymmetric Euler equations.
%
%   Evaluates the net x- and y-momentum fluxes through all cell faces
%   (east, west, south, north, and front/back for the x-component) using
%   midpoint interpolation and midpoint quadrature, both second-order
%   accurate on rectangular cells. The momentum flux includes the convective
%   term (rho*u_i * V_n) and the pressure contribution (p * n_i). The
%   x-momentum equation additionally receives a front/back pressure term
%   arising from the axisymmetric geometry. Results are normalised by cell
%   volume.
%
% Syntax:
%   s = FLUX_MOMENTUM_3D(s)
%
% Inputs:
%   s - Structure containing at minimum:
%                .var.rho                     - Density field (including ghosts).
%                .var.rho_u, .var.rho_v      - Momentum components (including ghosts).
%                .var.p                       - Pressure field (including ghosts).
%                .mesh.lr_x_normal, .mesh.lr_y_normal - Left/right face unit normals.
%                .mesh.bt_x_normal, .mesh.bt_y_normal - Bottom/top face unit normals.
%                .mesh.lr_area, .mesh.bt_area - Face areas.
%                .mesh.fb_y_normal_front,
%                .mesh.fb_y_normal_back      - Front/back face y-normals (3D).
%                .mesh.fb_area               - Front/back face areas (3D).
%                .mesh.y                     - Cell centroid y-coordinates.
%                .mesh.volume                - Cell volumes.
%
% Outputs:
%   s - Updated structure with:
%                .flux.rho_u - Cell-averaged x-momentum flux (interior cells).
%                .flux.rho_v - Cell-averaged y-momentum flux (interior cells).
%
% Notes:
%   - The front/back face pressure term only contributes to x-momentum.
%   - The sign(s.mesh.x) factor accounts for the direction of the
%     axisymmetric pressure contribution.
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
    
    % Front and back faces (3D axisymmetric pressure contribution)
    temp_P_face_front = s.var.p(2:end-1,2:end-1) * s.mesh.fb_y_normal_front;
    temp_P_face_back = s.var.p(2:end-1,2:end-1) * s.mesh.fb_y_normal_back;
    s.flux.rho_v = s.flux.rho_v + (temp_P_face_front - temp_P_face_back) .* s.mesh.fb_area .* sign(s.mesh.y);

    % Normalise by cell volume
    s.flux.rho_v = s.flux.rho_v ./ s.mesh.volume;
end
