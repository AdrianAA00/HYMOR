using Plots, Printf, LinearAlgebra, LoopVectorization

# PLOT_INFLECTION_POINTS  Compute and overlay generalized inflection-point
#   paths (extrema of rho * dU_t/dy_n) on the current figure.
#
#   PLOT_INFLECTION_POINTS(s)
#
#   Inputs:
#       s - Solution Dict containing grid coordinates
#                  (x_Ext, y_Ext), conserved variables (rho, rho_u, rho_v),
#                  shock geometry (shock_beta, shocked_cell_indices,
#                  x_wall_normal, y_wall_normal), and mesh sizes (Nx, Ny).
#
#   Outputs:
#       (none) - Plots maxima (pink solid) and minima (brown dashed)
#                paths of rho * dU_t / dy_n onto the active figure.
#
#   Notes:
#       The generalized inflection criterion rho * dU_t/dy_n is evaluated
#       at every streamwise station. Local extrema are detected with a
#       threshold epsilon, tracked across stations by nearest-neighbor
#       Euclidean matching, and smoothed with a moving-average filter.
#       Only paths longer than 40 points are displayed.
#
# Part of: Hypersonics Stability MATLAB Solver - Postprocessing Module

function PLOT_INFLECTION_POINTS(s::Dict{String, Any}; ax=nothing)
    do_plot = true
    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]
    rho_du = zeros(Nchi + 2, Neta + 1)
    grid_d_1 = zeros(Nchi + 2, Neta + 1)
    signs = ones(size(s["mesh"]["x_Ext"], 2))
    signs[1] = -1.0

    ## Compute rho * dU_t / dy_n at each streamwise station
    for r in 2:(Nchi - 1)
        shock_index = min(r - 1, Nchi)
        x_data = signs .* sqrt.((s["mesh"]["x_Ext"][r, :] .- s["mesh"]["x_Ext"][r, 1] / 2 .- s["mesh"]["x_Ext"][r, 2] / 2).^2 .+
                                (s["mesh"]["y_Ext"][r, :] .- s["mesh"]["y_Ext"][r, 1] / 2 .- s["mesh"]["y_Ext"][r, 2] / 2).^2)
        cell_x_wall_normal = (s["mesh"]["x_wall_normal"][1:end-1, 1] .+ s["mesh"]["x_wall_normal"][2:end, 1]) ./ 2
        cell_y_wall_normal = (s["mesh"]["y_wall_normal"][1:end-1, 1] .+ s["mesh"]["y_wall_normal"][2:end, 1]) ./ 2
        cell_x_wall_tan = cell_y_wall_normal
        cell_y_wall_tan = -cell_x_wall_normal
        u_tan_shock = s["var"]["rho_u"][r, :] ./ s["var"]["rho"][r, :] .* cos(s["shock"]["beta"][shock_index, 1]) .+
                      s["var"]["rho_v"][r, :] ./ s["var"]["rho"][r, :] .* sin(s["shock"]["beta"][shock_index, 1])
        u_tan_wall = s["var"]["rho_u"][r, :] ./ s["var"]["rho"][r, :] .* cell_x_wall_tan[shock_index, 1] .+
                     s["var"]["rho_v"][r, :] ./ s["var"]["rho"][r, :] .* cell_y_wall_tan[shock_index, 1]
        du_tan_shock = (u_tan_wall[2:end] .- u_tan_wall[1:end-1]) ./ (x_data[2:end] .- x_data[1:end-1])
        du_tan_shock[s["shock"]["cell_indices"][shock_index, 1]:end] .= 0.0
        rho_d = (s["var"]["rho"][r, 2:end] .+ s["var"]["rho"][r, 1:end-1]) ./ 2
        rho_du[r, :] = rho_d .* du_tan_shock
        grid_d_1[r, :] = (x_data[2:end] .+ x_data[1:end-1]) ./ 2
    end

    ## Find and plot extrema on top of density plot
    epsilon = 0.0002
    (max_paths, min_paths) = TRACK_EXTREMA(s, rho_du, epsilon, do_plot; ax=ax)
    return nothing
end


# ========================================================================
# TRACK_EXTREMA  Find and track local extrema of rho * dU/dy across
#   streamwise stations, then optionally overlay them on the current figure.
#
#   Inputs:
#       s - Solution Dict with grid and flow data.
#       rho_du   - Array of rho * dU/dy values (Nx+2 x Ny+1).
#       epsilon  - Threshold above/below neighbors for extremum detection.
#       do_plot  - Bool flag; if true, plot paths on current figure.
#
#   Outputs:
#       max_paths - Vector of smoothed maxima paths, each element = [x y] matrix.
#       min_paths - Vector of smoothed minima paths, each element = [x y] matrix.

function TRACK_EXTREMA(s::Dict{String, Any}, rho_du::Matrix{Float64}, epsilon::Float64, do_plot::Bool=false; ax=nothing)
    Nchi = s["mesh"]["Nchi"]

    ## Detect extrema at each streamwise station
    max_points = Vector{Union{Nothing, Matrix{Float64}}}(nothing, Nchi)
    min_points = Vector{Union{Nothing, Matrix{Float64}}}(nothing, Nchi)

    for r in 2:(Nchi - 1)
        shock_index = min(r - 1, Nchi)
        max_j = round(Int, s["shock"]["cell_indices"][shock_index] * 3 / 4)

        # Local maxima
        max_idx = Int[]
        for j in 2:min(size(rho_du, 2) - 1, max_j)
            if rho_du[r, j] > rho_du[r, j-1] + epsilon && rho_du[r, j] > rho_du[r, j+1] + epsilon
                push!(max_idx, j)
            end
        end
        if !isempty(max_idx)
            max_points[r] = hcat(s["mesh"]["x_Ext"][r, max_idx] .+ s["curvilinear_mapping"]["L"],
                                 s["mesh"]["y_Ext"][r, max_idx])
        end

        # Local minima
        min_idx = Int[]
        for j in 2:min(size(rho_du, 2) - 1, max_j)
            if rho_du[r, j] < rho_du[r, j-1] - epsilon && rho_du[r, j] < rho_du[r, j+1] - epsilon
                push!(min_idx, j)
            end
        end
        if !isempty(min_idx)
            min_points[r] = hcat(s["mesh"]["x_Ext"][r, min_idx] .+ s["curvilinear_mapping"]["L"],
                                 s["mesh"]["y_Ext"][r, min_idx])
        end
    end

    ## Track and smooth paths
    max_paths = TRACK_POINTS(max_points, Nchi)
    min_paths = TRACK_POINTS(min_points, Nchi)

    max_paths_smooth = SMOOTH_PATHS(max_paths)
    min_paths_smooth = SMOOTH_PATHS(min_paths)

    ## Overlay on current PyPlot figure
    if do_plot
        if ax === nothing
            ax = PyPlot.gca()
        end

        # Plot maxima paths (pink solid)
        for i in 1:length(max_paths_smooth)
            path = max_paths_smooth[i]
            if size(path, 1) > 40
                ax.plot(path[:, 1], path[:, 2], linestyle="-",
                        color=(1.0, 0.41, 0.71), linewidth=1.5)
            end
        end

        # Plot minima paths (brown dashed)
        for i in 1:length(min_paths_smooth)
            path = min_paths_smooth[i]
            if size(path, 1) > 40
                ax.plot(path[:, 1], path[:, 2], linestyle="--",
                        color=(0.55, 0.27, 0.07), linewidth=1.5)
            end
        end

        # Legend entries via invisible dummy lines
        ax.plot([NaN], [NaN], linestyle="-", color=(1.0, 0.41, 0.71),
                linewidth=1.5, label=L"Maxima ($\rho dU_t/dy_n$)")
        ax.plot([NaN], [NaN], linestyle="--", color=(0.55, 0.27, 0.07),
                linewidth=1.5, label=L"Minima ($\rho dU_t/dy_n$)")
        ax.legend(loc="upper left")
    end

    # Return smoothed paths
    return (max_paths_smooth, min_paths_smooth)
end


# ========================================================================
# TRACK_POINTS  Connect extrema across streamwise stations into continuous
#   paths using nearest-neighbor Euclidean-distance matching.
#
#   Inputs:
#       points - Vector (Nx) of Union{Nothing, Matrix{Float64}} of [x, y]
#                extremum coordinates per station.
#       Nx     - Number of streamwise stations.
#
#   Outputs:
#       paths  - Vector of tracked paths; each element is an N x 2
#                matrix of [x, y] coordinates.

function TRACK_POINTS(points::Vector{Union{Nothing, Matrix{Float64}}}, Nx::Int)
    paths = Vector{Matrix{Float64}}()
    used = Vector{Union{Nothing, Vector{Bool}}}(nothing, Nx)
    for r in 1:Nx
        if points[r] !== nothing
            used[r] = falses(size(points[r], 1))
        end
    end

    ## Adaptive maximum tracking distance
    avg_spacing = 0.0
    count = 0
    for r in 2:(Nx - 2)
        if points[r] !== nothing && points[r+1] !== nothing
            for p1 in 1:size(points[r], 1)
                for p2 in 1:size(points[r+1], 1)
                    dx = points[r+1][p2, 1] - points[r][p1, 1]
                    dy = points[r+1][p2, 2] - points[r][p1, 2]
                    avg_spacing += sqrt(dx^2 + dy^2)
                    count += 1
                end
            end
        end
    end

    if count > 0
        avg_spacing = avg_spacing / count
        max_distance = 3 * avg_spacing
    else
        max_distance = 0.5
    end

    @printf("Adaptive max_distance set to: %.4f\n", max_distance)

    ## Forward tracking from each unassigned point
    for r in 2:(Nx - 1)
        if points[r] === nothing
            continue
        end

        for p in 1:size(points[r], 1)
            if used[r][p]
                continue
            end

            path = points[r][p:p, :]  # 1x2 matrix
            used[r][p] = true
            current_point = points[r][p, :]

            for r_next in (r + 1):(Nx - 1)
                if points[r_next] === nothing
                    continue
                end

                dx_vec = points[r_next][:, 1] .- current_point[1]
                dy_vec = points[r_next][:, 2] .- current_point[2]
                distances = sqrt.(dx_vec.^2 .+ dy_vec.^2)

                sorted_idx = sortperm(distances)
                sorted_dist = distances[sorted_idx]
                found = false
                for k in 1:length(sorted_idx)
                    if !used[r_next][sorted_idx[k]] && sorted_dist[k] < max_distance
                        path = vcat(path, points[r_next][sorted_idx[k]:sorted_idx[k], :])
                        used[r_next][sorted_idx[k]] = true
                        current_point = points[r_next][sorted_idx[k], :]
                        found = true
                        break
                    end
                end

                if !found
                    break
                end
            end

            if size(path, 1) >= 2
                push!(paths, path)
            end
        end
    end

    return paths
end


# ========================================================================
# SMOOTH_PATHS  Apply a moving-average filter to tracked extrema paths
#   to eliminate discrete cell-to-cell jumps.
#
#   Inputs:
#       paths - Vector of raw paths with [x, y] coordinates (each Nx2 matrix).
#
#   Outputs:
#       smooth_paths - Vector of smoothed paths (same structure).

function SMOOTH_PATHS(paths::Vector{Matrix{Float64}})
    smooth_paths = Vector{Matrix{Float64}}(undef, length(paths))

    for i in 1:length(paths)
        path = paths[i]
        n_points = size(path, 1)

        if n_points < 5
            smooth_paths[i] = path
            continue
        end

        # Adaptive moving-average window (odd, at most n_points/3)
        window_size = min(5, floor(Int, n_points / 3))
        if mod(window_size, 2) == 0
            window_size = window_size + 1
        end

        # movmean equivalent: simple moving average with centered window
        half_w = div(window_size, 2)
        x_smooth = zeros(n_points)
        y_smooth = zeros(n_points)
        for j in 1:n_points
            lo = max(1, j - half_w)
            hi = min(n_points, j + half_w)
            x_smooth[j] = sum(path[lo:hi, 1]) / (hi - lo + 1)
            y_smooth[j] = sum(path[lo:hi, 2]) / (hi - lo + 1)
        end

        smooth_paths[i] = hcat(x_smooth, y_smooth)
    end

    return smooth_paths
end
