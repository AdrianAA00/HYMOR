function DERIVATIVE(u::AbstractMatrix, s::Dict{String, Any})
# DERIVATIVE - Compute spatial derivatives on a structured 2D grid.
#
# Computes the partial derivatives du/dx and du/dy of a scalar field u
# using central finite differences between cell centroids. The derivatives
# are first computed in the local cell-aligned curvilinear coordinate
# system (c1, c2) and then transformed to the Cartesian (x, y) system
# via a coordinate Jacobian.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract arrays to avoid repeated dictionary lookups
    mesh_x = s["mesh"]["x"]::Matrix{Float64}
    mesh_y = s["mesh"]["y"]::Matrix{Float64}

    ## Compute centroid-to-centroid distances
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

    ## Zero out derivatives at and beyond shocked cells
    if s["shock"]["enabled"]
        for i in 1:s["mesh"]["Nchi"]-2
            @views du_dx[i, s["shock"]["cell_indices"][i+1, 1]-1:end] .= 0
            @views du_dy[i, s["shock"]["cell_indices"][i+1, 1]-1:end] .= 0
        end
    end

    return du_dx, du_dy
end
