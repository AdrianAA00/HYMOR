function PLOT_INFLECTION_POINTS(s)
% PLOT_INFLECTION_POINTS  Compute and overlay generalized inflection-point
%   paths (extrema of rho * dU_t/dy_n) on the current figure.
%
%   PLOT_INFLECTION_POINTS(s)
%
%   Inputs:
%       s - Solution structure containing grid coordinates
%                  (x_Ext, y_Ext), conserved variables (rho, rho_u, rho_v),
%                  shock geometry (shock_beta, shocked_cell_indices,
%                  x_wall_normal, y_wall_normal), and mesh sizes (Nx, Ny).
%
%   Outputs:
%       (none) - Plots maxima (pink solid) and minima (brown dashed)
%                paths of rho * dU_t / dy_n onto the active figure.
%
%   Notes:
%       The generalized inflection criterion rho * dU_t/dy_n is evaluated
%       at every streamwise station. Local extrema are detected with a
%       threshold epsilon, tracked across stations by nearest-neighbor
%       Euclidean matching, and smoothed with a moving-average filter.
%       Only paths longer than 40 points are displayed.
%
% Part of: Hypersonics Stability MATLAB Solver - Postprocessing Module

    do_plot = true;
    rho_du = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 1);
    grid_d_1 = zeros(s.mesh.Nchi + 2, s.mesh.Neta + 1);
    signs = ones(size(s.mesh.x_Ext(1,:)));
    signs(1) = -1;

    %% Compute rho * dU_t / dy_n at each streamwise station
    for r = 2:s.mesh.Nchi - 1
        shock_index = min(r - 1, s.mesh.Nchi);
        x_data = signs .* sqrt((s.mesh.x_Ext(r,:) - s.mesh.x_Ext(r,1)/2 - s.mesh.x_Ext(r,2)/2).^2 + ...
                               (s.mesh.y_Ext(r,:) - s.mesh.y_Ext(r,1)/2 - s.mesh.y_Ext(r,2)/2).^2);
        cell_x_wall_normal = (s.mesh.x_wall_normal(1:end-1,1) + s.mesh.x_wall_normal(2:end,1)) / 2;
        cell_y_wall_normal = (s.mesh.y_wall_normal(1:end-1,1) + s.mesh.y_wall_normal(2:end,1)) / 2;
        cell_x_wall_tan = cell_y_wall_normal;
        cell_y_wall_tan = -cell_x_wall_normal;
        u_tan_shock = s.var.rho_u(r,:) ./ s.var.rho(r,:) .* cos(s.shock.beta(shock_index,1)) + ...
                      s.var.rho_v(r,:) ./ s.var.rho(r,:) .* sin(s.shock.beta(shock_index,1));
        u_tan_wall = s.var.rho_u(r,:) ./ s.var.rho(r,:) .* cell_x_wall_tan(shock_index,1) + ...
                     s.var.rho_v(r,:) ./ s.var.rho(r,:) .* cell_y_wall_tan(shock_index,1);
        du_tan_shock = (u_tan_wall(1,2:end) - u_tan_wall(1,1:end-1)) ./ (x_data(1,2:end) - x_data(1,1:end-1));
        du_tan_shock(1, s.shock.cell_indices(shock_index,1):end) = 0;
        rho_d = (s.var.rho(r,2:end) + s.var.rho(r,1:end-1)) / 2;
        rho_du(r,:) = rho_d .* du_tan_shock;
        grid_d_1(r,:) = (x_data(1,2:end) + x_data(1,1:end-1)) / 2;
    end

    %% Find and plot extrema on top of density plot
    epsilon = 0.0002;
    [max_paths, min_paths] = TRACK_EXTREMA(s, rho_du, epsilon, do_plot);
end

%% ========================================================================
function [max_paths, min_paths] = TRACK_EXTREMA(s, rho_du, epsilon, do_plot)
% TRACK_EXTREMA  Find and track local extrema of rho * dU/dy across
%   streamwise stations, then optionally overlay them on the current figure.
%
%   Inputs:
%       s - Solution structure with grid and flow data.
%       rho_du   - Array of rho * dU/dy values (Nx+2 x Ny+1).
%       epsilon  - Threshold above/below neighbors for extremum detection.
%       do_plot  - Logical flag; if true, plot paths on current figure.
%
%   Outputs:
%       max_paths - Cell array of smoothed maxima paths {i} = [x, y].
%       min_paths - Cell array of smoothed minima paths {i} = [x, y].

    if nargin < 4
        do_plot = false;
    end

    %% Detect extrema at each streamwise station
    max_points = cell(s.mesh.Nchi, 1);
    min_points = cell(s.mesh.Nchi, 1);

    for r = 2:s.mesh.Nchi - 1
        shock_index = min(r - 1, s.mesh.Nchi);
        max_j = round(s.shock.cell_indices(shock_index) * 3 / 4);

        % Local maxima
        max_idx = [];
        for j = 2:min(length(rho_du(r,:)) - 1, max_j)
            if rho_du(r,j) > rho_du(r,j-1) + epsilon && rho_du(r,j) > rho_du(r,j+1) + epsilon
                max_idx = [max_idx; j];
            end
        end
        if ~isempty(max_idx)
            max_points{r} = [s.mesh.x_Ext(r, max_idx)' + s.curvilinear_mapping.L, s.mesh.y_Ext(r, max_idx)'];
        end

        % Local minima
        min_idx = [];
        for j = 2:min(length(rho_du(r,:)) - 1, max_j)
            if rho_du(r,j) < rho_du(r,j-1) - epsilon && rho_du(r,j) < rho_du(r,j+1) - epsilon
                min_idx = [min_idx; j];
            end
        end
        if ~isempty(min_idx)
            min_points{r} = [s.mesh.x_Ext(r, min_idx)' + s.curvilinear_mapping.L, s.mesh.y_Ext(r, min_idx)'];
        end
    end

    %% Track and smooth paths
    max_paths = TRACK_POINTS(max_points, s.mesh.Nchi);
    min_paths = TRACK_POINTS(min_points, s.mesh.Nchi);

    max_paths_smooth = SMOOTH_PATHS(max_paths);
    min_paths_smooth = SMOOTH_PATHS(min_paths);

    %% Overlay on current figure
    if do_plot
        % Plot maxima paths (pink solid)
        for i = 1:length(max_paths_smooth)
            path = max_paths_smooth{i};
            if size(path, 1) > 40
                plot(path(:,1), path(:,2), '-', 'Color', [1, 0.41, 0.71], 'LineWidth', 1.5);
            end
        end

        % Plot minima paths (brown dashed)
        for i = 1:length(min_paths_smooth)
            path = min_paths_smooth{i};
            if size(path, 1) > 40
                plot(path(:,1), path(:,2), '--', 'Color', [0.55, 0.27, 0.07], 'LineWidth', 1.5);
            end
        end

        % Legend entries via invisible dummy lines
        h_dummy_max = plot(NaN, NaN, '-', 'Color', [1, 0.41, 0.71], 'LineWidth', 1.5);
        h_dummy_min = plot(NaN, NaN, '--', 'Color', [0.55, 0.27, 0.07], 'LineWidth', 1.5);
        legend([h_dummy_max, h_dummy_min], ...
            {'Maxima ($\rho dU_t/dy_n$)', 'Minima ($\rho dU_t/dy_n$)'}, ...
            'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 16);
    end

    % Return smoothed paths
    max_paths = max_paths_smooth;
    min_paths = min_paths_smooth;
end

%% ========================================================================
function paths = TRACK_POINTS(points, Nx)
% TRACK_POINTS  Connect extrema across streamwise stations into continuous
%   paths using nearest-neighbor Euclidean-distance matching.
%
%   Inputs:
%       points - Cell array (Nx x 1) of [x, y] extremum coordinates per
%                station.
%       Nx     - Number of streamwise stations.
%
%   Outputs:
%       paths  - Cell array of tracked paths; each element is an N x 2
%                matrix of [x, y] coordinates.

    paths = {};
    used = cell(Nx, 1);
    for r = 1:Nx
        if ~isempty(points{r})
            used{r} = false(size(points{r}, 1), 1);
        end
    end

    %% Adaptive maximum tracking distance
    avg_spacing = 0;
    count = 0;
    for r = 2:Nx - 2
        if ~isempty(points{r}) && ~isempty(points{r+1})
            for p1 = 1:size(points{r}, 1)
                for p2 = 1:size(points{r+1}, 1)
                    dx = points{r+1}(p2,1) - points{r}(p1,1);
                    dy = points{r+1}(p2,2) - points{r}(p1,2);
                    avg_spacing = avg_spacing + sqrt(dx^2 + dy^2);
                    count = count + 1;
                end
            end
        end
    end

    if count > 0
        avg_spacing = avg_spacing / count;
        max_distance = 3 * avg_spacing;
    else
        max_distance = 0.5;
    end

    fprintf('Adaptive max_distance set to: %.4f\n', max_distance);

    %% Forward tracking from each unassigned point
    for r = 2:Nx - 1
        if isempty(points{r})
            continue;
        end

        for p = 1:size(points{r}, 1)
            if used{r}(p)
                continue;
            end

            path = points{r}(p, :);
            used{r}(p) = true;
            current_point = points{r}(p, :);

            for r_next = r + 1:Nx - 1
                if isempty(points{r_next})
                    continue;
                end

                dx = points{r_next}(:,1) - current_point(1);
                dy = points{r_next}(:,2) - current_point(2);
                distances = sqrt(dx.^2 + dy.^2);

                [sorted_dist, sorted_idx] = sort(distances);
                found = false;
                for k = 1:length(sorted_idx)
                    if ~used{r_next}(sorted_idx(k)) && sorted_dist(k) < max_distance
                        path = [path; points{r_next}(sorted_idx(k), :)];
                        used{r_next}(sorted_idx(k)) = true;
                        current_point = points{r_next}(sorted_idx(k), :);
                        found = true;
                        break;
                    end
                end

                if ~found
                    break;
                end
            end

            if size(path, 1) >= 2
                paths{end+1} = path;
            end
        end
    end
end

%% ========================================================================
function smooth_paths = SMOOTH_PATHS(paths)
% SMOOTH_PATHS  Apply a moving-average filter to tracked extrema paths
%   to eliminate discrete cell-to-cell jumps.
%
%   Inputs:
%       paths - Cell array of raw paths with [x, y] coordinates.
%
%   Outputs:
%       smooth_paths - Cell array of smoothed paths (same structure).

    smooth_paths = cell(size(paths));

    for i = 1:length(paths)
        path = paths{i};
        n_points = size(path, 1);

        if n_points < 5
            smooth_paths{i} = path;
            continue;
        end

        % Adaptive moving-average window (odd, at most n_points/3)
        window_size = min(5, floor(n_points / 3));
        if mod(window_size, 2) == 0
            window_size = window_size + 1;
        end

        x_smooth = movmean(path(:,1), window_size);
        y_smooth = movmean(path(:,2), window_size);

        smooth_paths{i} = [x_smooth, y_smooth];
    end
end
