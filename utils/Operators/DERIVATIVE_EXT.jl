function DERIVATIVE_EXT(u::AbstractMatrix, s::Dict{String, Any})
# DERIVATIVE_EXT - Compute spatial derivatives on the extended grid.
#
# Computes the partial derivatives du/dx and du/dy of a scalar field u
# using central finite differences between cell centroids on the extended
# grid (x_Ext, y_Ext). Near-shock extrapolation uses a linear (2nd-order)
# scheme to avoid differencing across the discontinuity.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract extended mesh coordinate arrays (ghost-padded)
    mesh_x = s["mesh"]["x_Ext"]::Matrix{Float64}
    mesh_y = s["mesh"]["y_Ext"]::Matrix{Float64}

    ## Compute centroid-to-centroid distances on extended grid
    @views centroids_bt_distance = @. sqrt((mesh_y[2:end-1, 3:end] - mesh_y[2:end-1, 1:end-2])^2 +
                                           (mesh_x[2:end-1, 3:end] - mesh_x[2:end-1, 1:end-2])^2)
    @views centroids_lr_distance = @. sqrt((mesh_y[3:end, 2:end-1] - mesh_y[1:end-2, 2:end-1])^2 +
                                           (mesh_x[3:end, 2:end-1] - mesh_x[1:end-2, 2:end-1])^2)

    ## Compute derivatives in cell-aligned coordinate system
    @views du_dc1 = @. -(u[3:end, 2:end-1] - u[1:end-2, 2:end-1]) / centroids_lr_distance
    @views du_dc2 = @.  (u[2:end-1, 3:end] - u[2:end-1, 1:end-2]) / centroids_bt_distance

    ## Project cell coordinate system onto Cartesian (x, y)
    @views c1_x = @. -(mesh_x[3:end, 2:end-1] - mesh_x[1:end-2, 2:end-1]) / centroids_lr_distance
    @views c1_y = @. -(mesh_y[3:end, 2:end-1] - mesh_y[1:end-2, 2:end-1]) / centroids_lr_distance
    @views c2_x = @.  (mesh_x[2:end-1, 3:end] - mesh_x[2:end-1, 1:end-2]) / centroids_bt_distance
    @views c2_y = @.  (mesh_y[2:end-1, 3:end] - mesh_y[2:end-1, 1:end-2]) / centroids_bt_distance

    ## Transform derivatives to Cartesian system
    divide = @. c1_x * c2_y - c2_x * c1_y
    du_dx = @. (c2_y * du_dc1 - c1_y * du_dc2) / divide
    du_dy = @. (c1_x * du_dc2 - c2_x * du_dc1) / divide

    ## Extrapolate derivatives near shocked cells
    if s["shock"]["enabled"]
        # Compute linear indices for shocked cells and their neighbors
        Nchi = s["mesh"]["Nchi"]::Int
        j_values_s = s["shock"]["cell_indices"]::Matrix{Int}

        idx_3 = [CartesianIndex(k, j_values_s[k, 1] + 3) for k in 1:Nchi]
        idx_2 = [CartesianIndex(k, j_values_s[k, 1] + 2) for k in 1:Nchi]
        idx_1 = [CartesianIndex(k, j_values_s[k, 1] + 1) for k in 1:Nchi]
        idx_0 = [CartesianIndex(k, j_values_s[k, 1])     for k in 1:Nchi]

        # Extrapolate du_dx near shock (2nd-order linear extrapolation)
        temp_x = copy(du_dx)
        du_dx[idx_1] .= 2 .* temp_x[idx_2] .- temp_x[idx_3]
        du_dx[idx_0] .= 2 .* du_dx[idx_1] .- du_dx[idx_2]

        # Extrapolate du_dy near shock (2nd-order linear extrapolation)
        temp_y = copy(du_dy)
        du_dy[idx_1] .= 2 .* temp_y[idx_2] .- temp_y[idx_3]
        du_dy[idx_0] .= 2 .* du_dy[idx_1] .- du_dy[idx_2]
    end

    return du_dx, du_dy
end
