function s = FLUX_MASS(s)
% FLUX_MASS - Compute the inviscid mass flux for the 2D continuity equation.
%
%   Evaluates the net mass flux through all four faces (east, west, south,
%   north) of each interior cell using midpoint interpolation and midpoint
%   quadrature, both second-order accurate on rectangular cells. The
%   resulting flux is normalised by cell volume to yield the cell-averaged
%   rate of change of density.
%
% Syntax:
%   s = FLUX_MASS(s)
%
% Inputs:
%   s - Structure containing at minimum:
%                .var.rho_u, .var.rho_v      - Momentum components (including ghosts).
%                .mesh.lr_x_normal, .mesh.lr_y_normal - Left/right face unit normals.
%                .mesh.bt_x_normal, .mesh.bt_y_normal - Bottom/top face unit normals.
%                .mesh.lr_area, .mesh.bt_area - Face areas.
%                .mesh.volume                - Cell volumes.
%
% Outputs:
%   s - Updated structure with:
%                .flux.rho - Cell-averaged mass flux (interior cells).
%
% Notes:
%   - Midpoint interpolation: face value = average of adjacent cell values.
%   - Midpoint quadrature: face flux = face value * face area.
%   - Sign convention: outgoing flux is negative, incoming flux is positive.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% East and west faces
    % Midpoint interpolation (2nd order on rectangles)
    temp_rho_u_face = (s.var.rho_u(1:end-1,2:end-1) + s.var.rho_u(2:end,2:end-1)) .* s.mesh.lr_x_normal / 2;
    temp_rho_u_face = temp_rho_u_face + (s.var.rho_v(1:end-1,2:end-1) + s.var.rho_v(2:end,2:end-1)) .* s.mesh.lr_y_normal / 2;

    % Midpoint quadrature (2nd order on rectangles)
    s.flux.rho = -temp_rho_u_face(2:end,:) .* s.mesh.lr_area(2:end,:);            % Outgoing fluxes
    s.flux.rho = s.flux.rho + temp_rho_u_face(1:end-1,:) .* s.mesh.lr_area(1:end-1,:);  % Incoming fluxes

    %% South and north faces
    % Midpoint interpolation (2nd order on rectangles)
    temp_rho_u_face2 = (s.var.rho_u(2:end-1,1:end-1) + s.var.rho_u(2:end-1,2:end)) .* s.mesh.bt_x_normal / 2;
    temp_rho_u_face2 = temp_rho_u_face2 + (s.var.rho_v(2:end-1,1:end-1) + s.var.rho_v(2:end-1,2:end)) .* s.mesh.bt_y_normal / 2;

    % Midpoint quadrature (2nd order on rectangles)
    s.flux.rho = s.flux.rho - temp_rho_u_face2(:,2:end) .* s.mesh.bt_area(:,2:end);            % Outgoing fluxes
    s.flux.rho = s.flux.rho + temp_rho_u_face2(:,1:end-1) .* s.mesh.bt_area(:,1:end-1);  % Incoming fluxes

    %% Normalise by cell volume
    s.flux.rho = s.flux.rho ./ s.mesh.volume;
end
