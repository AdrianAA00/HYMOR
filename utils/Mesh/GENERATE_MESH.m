function s = GENERATE_MESH(s, disturbances)
% GENERATE_MESH  Build structured mesh from uniform element space [0,1]x[0,1].
%
%   s = GENERATE_MESH(s, disturbances)
%
%   Constructs a body-fitted structured grid by creating a uniform mesh in
%   element space [0,1]x[0,1] and mapping it to physical space via
%   GO_TO_PHYSICAL_SPACE. All geometry-specific logic is contained within
%   the coordinate transform functions. This function handles stagnation
%   refinement, shock alignment, mesh smoothing, and computes all cell-face
%   areas, normals, and volumes needed by the finite-volume solver.
%
%   Inputs:
%       s            - Struct containing mesh and curvilinear_mapping params
%       disturbances - Logical flag; when true, mesh rows upstream of the
%                      shock are made parallel to y = 0.
%
%   Outputs:
%       s - Updated struct with mesh coordinates, face areas/normals,
%           volumes, and geometry factors.
%
%   Notes:
%       - Element space is always [0,1]x[0,1] regardless of geometry.
%       - Stagnation refinement clusters chi near chi=1 (stagnation point).
%       - For shock restarts, the mesh is extruded from the wall along
%         shock-aligned normals instead of using the transform directly.
%       - Wall refinement is applied only during shock restarts.
%       - Mesh smoothing is controlled by s.curvilinear_mapping.smooth_mesh.
%
% Part of: Hypersonics Stability MATLAB Solver - Mesh Module

    %% Create uniform element space [0,1]
    chi_1d = linspace(0, 1, s.mesh.Nchi+1)';  % (Nchi+1 x 1)
    eta_1d = linspace(0, 1, s.mesh.Neta+1);   % (1 x Neta+1)

    %% Apply stagnation refinement to chi (clusters near chi=1)
    if s.curvilinear_mapping.refinement_stagnation.state
        chi_1d = REFINEMENT(chi_1d, ...
            s.curvilinear_mapping.refinement_stagnation.BL_thickness, ...
            0, ...
            s.curvilinear_mapping.refinement_stagnation.intensity, ...
            1, 0);
    end

    %% Store wall element-space coordinates
    s.mesh.chi_wall = chi_1d;
    s.mesh.eta_wall = zeros(s.mesh.Nchi+1, 1);

    if s.restart && s.shock.enabled
        %% Restart with shock: extrude from wall along shock-aligned normals

        % Get wall positions from transform
        [s.mesh.x_wall, s.mesh.y_wall] = GO_TO_PHYSICAL_SPACE( ...
            chi_1d, zeros(s.mesh.Nchi+1, 1), s);

        % Get shock positions from spline at wall chi values
        shocked_points_x = ppval(s.shock.spline_func_x, chi_1d);
        shocked_points_y = ppval(s.shock.spline_func_y, chi_1d);

        % Compute shock-wall distances and reoriented normals
        dx_shock = shocked_points_x - s.mesh.x_wall;
        dy_shock = shocked_points_y - s.mesh.y_wall;
        shocked_r = sqrt(dx_shock.^2 + dy_shock.^2);
        s.mesh.x_wall_normal = dx_shock ./ shocked_r;
        s.mesh.y_wall_normal = dy_shock ./ shocked_r;

        % Compute wall-normal distances with optional refinement
        index = round((s.mesh.Neta+1) / s.shock.remesh_shock_distance);
        wall_distance = zeros(s.mesh.Nchi+1, s.mesh.Neta+1);

        for i = 1:s.mesh.Nchi+1
            wall_distance_line = (0:s.mesh.Neta) * shocked_r(i) / (index - 1/2);

            if s.curvilinear_mapping.refinement_wall.state
                wall_distance_line = REFINEMENT(wall_distance_line, ...
                    s.curvilinear_mapping.refinement_wall.BL_thickness, ...
                    0, ...
                    s.curvilinear_mapping.refinement_wall.intensity, ...
                    max(wall_distance_line), 0);

                % Rescale so shock passes through cell center at index
                shocked_cell_center = (wall_distance_line(index) ...
                    + wall_distance_line(index+1)) / 2;
                factor = shocked_r(i) / shocked_cell_center;
                wall_distance_line = factor * wall_distance_line;
            end

            wall_distance(i,:) = wall_distance_line;
        end

        % Compute corners from wall + normal extrusion
        s.mesh.x_corner = s.mesh.x_wall + wall_distance .* s.mesh.x_wall_normal;
        s.mesh.y_corner = s.mesh.y_wall + wall_distance .* s.mesh.y_wall_normal;

    else
        %% Non-restart: map element space to physical space via transform
        chi_2D = repmat(chi_1d, 1, s.mesh.Neta+1);
        eta_2D = repmat(eta_1d, s.mesh.Nchi+1, 1);
        [s.mesh.x_corner, s.mesh.y_corner] = GO_TO_PHYSICAL_SPACE(chi_2D, eta_2D, s);

        % Extract wall coordinates
        s.mesh.x_wall = s.mesh.x_corner(:, 1);
        s.mesh.y_wall = s.mesh.y_corner(:, 1);

        % Compute wall normals from transform direction (eta=0 to eta=delta)
        dx = s.mesh.x_corner(:, 2) - s.mesh.x_corner(:, 1);
        dy = s.mesh.y_corner(:, 2) - s.mesh.y_corner(:, 1);
        norm_len = sqrt(dx.^2 + dy.^2);
        s.mesh.x_wall_normal = dx ./ norm_len;
        s.mesh.y_wall_normal = dy ./ norm_len;
    end

    %% Make cells upstream from shock parallel to y = 0
    if s.shock.enabled
        index = round((s.mesh.Neta+1) / s.shock.remesh_shock_distance);
        if disturbances
            s.mesh.y_corner(:, index:end) = repmat(s.mesh.y_corner(:,index), ...
                1, size(s.mesh.y_corner,2)-index+1);
            s.mesh.x_corner(:, index:end) = s.mesh.x_corner(end, index:end) ...
                + s.mesh.x_corner(:,index) - s.mesh.x_corner(end,index);
        end
    end

    %% Compute cell centers and extended grid
    s = GET_EXTENDED_CELL_CENTERS(s);

    s.mesh.arc_3D_Dtheta_Dr = 0.001; % Arc angle for 3D axisymmetric wedge

    %% Mesh smoothing for curvature discontinuities
    if isfield(s.curvilinear_mapping, 'smooth_mesh') ...
            && s.curvilinear_mapping.smooth_mesh && ~s.restart
        for iter = 1:10
            s = GET_CELL_CORNERS_SMOOTH(s);
            s = GET_EXTENDED_CELL_CENTERS(s);
        end
    end

    %% Compute left-right face areas and normals
    if s.PDE_dimension == "2D"
        for i = 1:s.mesh.Nchi+1
            for j = 1:s.mesh.Neta
                s.mesh.lr_area(i,j) = sqrt((s.mesh.x_corner(i,j+1) - s.mesh.x_corner(i,j))^2 ...
                                           + (s.mesh.y_corner(i,j+1) - s.mesh.y_corner(i,j))^2);
                s.mesh.lr_x_normal(i,j) = -(s.mesh.y_corner(i,j+1) - s.mesh.y_corner(i,j)) / s.mesh.lr_area(i,j);
                s.mesh.lr_y_normal(i,j) =  (s.mesh.x_corner(i,j+1) - s.mesh.x_corner(i,j)) / s.mesh.lr_area(i,j);
                if i < s.mesh.Nchi+1
                    sgn = sign(s.mesh.lr_x_normal(i,j) * (s.mesh.x_corner(i+1,j) - s.mesh.x_corner(i,j)) ...
                             + s.mesh.lr_y_normal(i,j) * (s.mesh.y_corner(i+1,j) - s.mesh.y_corner(i,j)));
                else
                    sgn = sign(s.mesh.lr_x_normal(i,j) * (s.mesh.x_corner(i,j) - s.mesh.x_corner(i-1,j)) ...
                             + s.mesh.lr_y_normal(i,j) * (s.mesh.y_corner(i,j) - s.mesh.y_corner(i-1,j)));
                end
                s.mesh.lr_x_normal(i,j) =  s.mesh.lr_x_normal(i,j) * sgn;
                s.mesh.lr_y_normal(i,j) =  s.mesh.lr_y_normal(i,j) * sgn;
            end
        end
    elseif s.PDE_dimension == "3D-axisymmetric"
        for i = 1:s.mesh.Nchi+1
            for j = 1:s.mesh.Neta
                depth = abs(s.mesh.y_corner(i,j+1) + s.mesh.y_corner(i,j)) / 2 .* s.mesh.arc_3D_Dtheta_Dr;
                edge_len = sqrt((s.mesh.x_corner(i,j+1) - s.mesh.x_corner(i,j))^2 ...
                              + (s.mesh.y_corner(i,j+1) - s.mesh.y_corner(i,j))^2);
                s.mesh.lr_area(i,j) = edge_len * depth;
                s.mesh.lr_x_normal(i,j) = -(s.mesh.y_corner(i,j+1) - s.mesh.y_corner(i,j)) / edge_len;
                s.mesh.lr_y_normal(i,j) =  (s.mesh.x_corner(i,j+1) - s.mesh.x_corner(i,j)) / edge_len;
                if i < s.mesh.Nchi+1
                    sgn = sign(s.mesh.lr_x_normal(i,j) * (s.mesh.x_corner(i+1,j) - s.mesh.x_corner(i,j)) ...
                             + s.mesh.lr_y_normal(i,j) * (s.mesh.y_corner(i+1,j) - s.mesh.y_corner(i,j)));
                else
                    sgn = sign(s.mesh.lr_x_normal(i,j) * (s.mesh.x_corner(i,j) - s.mesh.x_corner(i-1,j)) ...
                             + s.mesh.lr_y_normal(i,j) * (s.mesh.y_corner(i,j) - s.mesh.y_corner(i-1,j)));
                end
                s.mesh.lr_x_normal(i,j) =  s.mesh.lr_x_normal(i,j) * sgn;
                s.mesh.lr_y_normal(i,j) =  s.mesh.lr_y_normal(i,j) * sgn;
            end
        end
    else
        error("Dimension selection not recognized")
    end

    %% Compute bottom-top face areas and normals
    if s.PDE_dimension == "2D"
        for i = 1:s.mesh.Nchi
            for j = 1:s.mesh.Neta+1
                s.mesh.bt_area(i,j) = sqrt((s.mesh.x_corner(i+1,j) - s.mesh.x_corner(i,j))^2 ...
                                           + (s.mesh.y_corner(i+1,j) - s.mesh.y_corner(i,j))^2);
                s.mesh.bt_x_normal(i,j) =  (s.mesh.y_corner(i+1,j) - s.mesh.y_corner(i,j)) / s.mesh.bt_area(i,j);
                s.mesh.bt_y_normal(i,j) = -(s.mesh.x_corner(i+1,j) - s.mesh.x_corner(i,j)) / s.mesh.bt_area(i,j);
                if j < s.mesh.Neta+1
                    sgn = sign(s.mesh.bt_x_normal(i,j) * (s.mesh.x_corner(i,j+1) - s.mesh.x_corner(i,j)) ...
                             + s.mesh.bt_y_normal(i,j) * (s.mesh.y_corner(i,j+1) - s.mesh.y_corner(i,j)));
                else
                    sgn = sign(s.mesh.bt_x_normal(i,j) * (s.mesh.x_corner(i,j) - s.mesh.x_corner(i,j-1)) ...
                             + s.mesh.bt_y_normal(i,j) * (s.mesh.y_corner(i,j) - s.mesh.y_corner(i,j-1)));
                end
                s.mesh.bt_x_normal(i,j) =  s.mesh.bt_x_normal(i,j) * sgn;
                s.mesh.bt_y_normal(i,j) =  s.mesh.bt_y_normal(i,j) * sgn;
            end
        end
    elseif s.PDE_dimension == "3D-axisymmetric"
        for i = 1:s.mesh.Nchi
            for j = 1:s.mesh.Neta+1
                depth = abs(s.mesh.y_corner(i+1,j) + s.mesh.y_corner(i,j)) / 2 .* s.mesh.arc_3D_Dtheta_Dr;
                edge_len = sqrt((s.mesh.x_corner(i+1,j) - s.mesh.x_corner(i,j))^2 ...
                              + (s.mesh.y_corner(i+1,j) - s.mesh.y_corner(i,j))^2);
                s.mesh.bt_area(i,j) = edge_len * depth;
                s.mesh.bt_x_normal(i,j) =  (s.mesh.y_corner(i+1,j) - s.mesh.y_corner(i,j)) / edge_len;
                s.mesh.bt_y_normal(i,j) = -(s.mesh.x_corner(i+1,j) - s.mesh.x_corner(i,j)) / edge_len;
                if j < s.mesh.Neta+1
                    sgn = sign(s.mesh.bt_x_normal(i,j) * (s.mesh.x_corner(i,j+1) - s.mesh.x_corner(i,j)) ...
                             + s.mesh.bt_y_normal(i,j) * (s.mesh.y_corner(i,j+1) - s.mesh.y_corner(i,j)));
                else
                    sgn = sign(s.mesh.bt_x_normal(i,j) * (s.mesh.x_corner(i,j) - s.mesh.x_corner(i,j-1)) ...
                             + s.mesh.bt_y_normal(i,j) * (s.mesh.y_corner(i,j) - s.mesh.y_corner(i,j-1)));
                end
                s.mesh.bt_x_normal(i,j) =  s.mesh.bt_x_normal(i,j) * sgn;
                s.mesh.bt_y_normal(i,j) =  s.mesh.bt_y_normal(i,j) * sgn;
            end
        end
    else
        error("Dimension selection not recognized")
    end

    %% Compute front-back face areas (3D-axisymmetric only)
    if s.PDE_dimension == "3D-axisymmetric"
        for i = 1:s.mesh.Nchi
            for j = 1:s.mesh.Neta
                % Front-back area equals the 2D cell area (Shoelace formula)
                s.mesh.fb_area(i,j) = -s.mesh.x_corner(i,j)     * s.mesh.y_corner(i+1,j)   + s.mesh.x_corner(i+1,j)   * s.mesh.y_corner(i,j);
                s.mesh.fb_area(i,j) = s.mesh.fb_area(i,j) - s.mesh.x_corner(i+1,j)   * s.mesh.y_corner(i+1,j+1) + s.mesh.x_corner(i+1,j+1) * s.mesh.y_corner(i+1,j);
                s.mesh.fb_area(i,j) = s.mesh.fb_area(i,j) - s.mesh.x_corner(i+1,j+1) * s.mesh.y_corner(i,j+1)   + s.mesh.x_corner(i,j+1)   * s.mesh.y_corner(i+1,j+1);
                s.mesh.fb_area(i,j) = s.mesh.fb_area(i,j) - s.mesh.x_corner(i,j+1)   * s.mesh.y_corner(i,j)     + s.mesh.x_corner(i,j)     * s.mesh.y_corner(i,j+1);
                s.mesh.fb_area(i,j) = abs(s.mesh.fb_area(i,j) / 2);
                theta_angle_fb_faces = atan(s.mesh.arc_3D_Dtheta_Dr / 2);
                s.mesh.fb_x_normal_front =  0;
                s.mesh.fb_y_normal_front =  sin(theta_angle_fb_faces);
                s.mesh.fb_x_normal_back  =  0;
                s.mesh.fb_y_normal_back  = -sin(theta_angle_fb_faces);
            end
        end
    end

    %% Compute cell volumes
    if s.PDE_dimension == "2D"
        for i = 1:s.mesh.Nchi
            for j = 1:s.mesh.Neta
                % Shoelace formula for quadrilateral area
                s.mesh.volume(i,j) = -s.mesh.x_corner(i,j)     * s.mesh.y_corner(i+1,j)   + s.mesh.x_corner(i+1,j)   * s.mesh.y_corner(i,j);
                s.mesh.volume(i,j) = s.mesh.volume(i,j) - s.mesh.x_corner(i+1,j)   * s.mesh.y_corner(i+1,j+1) + s.mesh.x_corner(i+1,j+1) * s.mesh.y_corner(i+1,j);
                s.mesh.volume(i,j) = s.mesh.volume(i,j) - s.mesh.x_corner(i+1,j+1) * s.mesh.y_corner(i,j+1)   + s.mesh.x_corner(i,j+1)   * s.mesh.y_corner(i+1,j+1);
                s.mesh.volume(i,j) = s.mesh.volume(i,j) - s.mesh.x_corner(i,j+1)   * s.mesh.y_corner(i,j)     + s.mesh.x_corner(i,j)     * s.mesh.y_corner(i,j+1);
                s.mesh.volume(i,j) = abs(s.mesh.volume(i,j) / 2);
            end
        end
    elseif s.PDE_dimension == "3D-axisymmetric"
        for i = 1:s.mesh.Nchi
            for j = 1:s.mesh.Neta
                % Compute corners of wedge volume element
                c1 = [s.mesh.x_corner(i,j),     s.mesh.y_corner(i,j),      s.mesh.y_corner(i,j)/2     .* s.mesh.arc_3D_Dtheta_Dr];
                c2 = [s.mesh.x_corner(i+1,j),   s.mesh.y_corner(i+1,j),    s.mesh.y_corner(i+1,j)/2   .* s.mesh.arc_3D_Dtheta_Dr];
                c3 = [s.mesh.x_corner(i+1,j+1), s.mesh.y_corner(i+1,j+1),  s.mesh.y_corner(i+1,j+1)/2 .* s.mesh.arc_3D_Dtheta_Dr];
                c4 = [s.mesh.x_corner(i,j+1),   s.mesh.y_corner(i,j+1),    s.mesh.y_corner(i,j+1)/2   .* s.mesh.arc_3D_Dtheta_Dr];
                c5 = [s.mesh.x_corner(i,j),     s.mesh.y_corner(i,j),     -s.mesh.y_corner(i,j)/2     .* s.mesh.arc_3D_Dtheta_Dr];
                c6 = [s.mesh.x_corner(i+1,j),   s.mesh.y_corner(i+1,j),   -s.mesh.y_corner(i+1,j)/2   .* s.mesh.arc_3D_Dtheta_Dr];
                c7 = [s.mesh.x_corner(i+1,j+1), s.mesh.y_corner(i+1,j+1), -s.mesh.y_corner(i+1,j+1)/2 .* s.mesh.arc_3D_Dtheta_Dr];
                c8 = [s.mesh.x_corner(i,j+1),   s.mesh.y_corner(i,j+1),   -s.mesh.y_corner(i,j+1)/2   .* s.mesh.arc_3D_Dtheta_Dr];

                % Decompose hexahedron into six tetrahedra and sum volumes
                tet1 = tetrahedronVolume(c1, c2, c4, c6);
                tet2 = tetrahedronVolume(c2, c6, c4, c3);
                tet3 = tetrahedronVolume(c6, c4, c3, c7);
                tet4 = tetrahedronVolume(c1, c4, c8, c6);
                tet5 = tetrahedronVolume(c8, c4, c6, c7);
                tet6 = tetrahedronVolume(c5, c1, c6, c8);
                s.mesh.volume(i,j) = tet1 + tet2 + tet3 + tet4 + tet5 + tet6;
            end
        end
    else
        error("Dimension selection not recognized")
    end

    %% Compute geometry factors for CFL estimation
    max_1 = max(s.mesh.bt_area(:,2:end) ./ s.mesh.volume, [], "all");
    max_2 = max(s.mesh.lr_area(2:end,:) ./ s.mesh.volume, [], "all");
    centroids_bt_distance = sqrt((s.mesh.y_Ext(2:end-1, 2:end) - s.mesh.y_Ext(2:end-1, 1:end-1)).^2 ...
                                + (s.mesh.x_Ext(2:end-1, 2:end) - s.mesh.x_Ext(2:end-1, 1:end-1)).^2);
    centroids_lr_distance = sqrt((s.mesh.y_Ext(2:end, 2:end-1) - s.mesh.y_Ext(1:end-1, 2:end-1)).^2 ...
                                + (s.mesh.x_Ext(2:end, 2:end-1) - s.mesh.x_Ext(1:end-1, 2:end-1)).^2);
    max_3 = max(s.mesh.bt_area(:,2:end) ./ s.mesh.volume ./ centroids_bt_distance(:,2:end), [], "all");
    max_4 = max(s.mesh.lr_area(2:end,:) ./ s.mesh.volume ./ centroids_lr_distance(2:end,:), [], "all");
    s.mesh.geometry_factor = max(max_1, max_2);
    s.mesh.volume_factor = max(max_3, max_4);
end


function volume = tetrahedronVolume(p1, p2, p3, p4)
% tetrahedronVolume  Volume of a tetrahedron from four vertex coordinates.
%
%   volume = tetrahedronVolume(p1, p2, p3, p4)
%
%   Computes the volume using the scalar triple product:
%       V = |v1 . (v2 x v3)| / 6
%   where v1, v2, v3 are edge vectors from p1 to the other vertices.

    v1 = p2 - p1;
    v2 = p3 - p1;
    v3 = p4 - p1;

    volume = abs(dot(v1, cross(v2, v3))) / 6;
end
