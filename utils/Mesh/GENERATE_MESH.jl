# GENERATE_MESH  Build structured mesh from uniform element space [0,1]x[0,1].
#
#   s = GENERATE_MESH(s, disturbances)
#
#   Constructs a body-fitted structured grid by creating a uniform mesh in
#   element space [0,1]x[0,1] and mapping it to physical space via
#   GO_TO_PHYSICAL_SPACE. All geometry-specific logic is contained within
#   the coordinate transform functions. This function handles stagnation
#   refinement, shock alignment, mesh smoothing, and computes all cell-face
#   areas, normals, and volumes needed by the finite-volume solver.
#
#   Inputs:
#       s            - Dict containing mesh and curvilinear_mapping params
#       disturbances - Bool flag; when true, mesh rows upstream of the
#                      shock are made parallel to y = 0.
#
#   Outputs:
#       s - Updated Dict with mesh coordinates, face areas/normals,
#           volumes, and geometry factors.
#
#   Notes:
#       - Element space is always [0,1]x[0,1] regardless of geometry.
#       - Stagnation refinement clusters chi near chi=1 (stagnation point).
#       - For shock restarts, the mesh is extruded from the wall along
#         shock-aligned normals instead of using the transform directly.
#       - Wall refinement is applied only during shock restarts.
#       - Mesh smoothing is controlled by s["curvilinear_mapping"]["smooth_mesh"].
#
# Part of: Hypersonics Stability Julia Solver - Mesh Module

using LinearAlgebra

function GENERATE_MESH(s::Dict{String,Any}, disturbances::Bool)

    ## Create uniform element space [0,1]
    chi_1d = collect(range(0, 1, length=s["mesh"]["Nchi"]+1))   # (Nchi+1,)
    eta_1d = collect(range(0, 1, length=s["mesh"]["Neta"]+1))'  # (1 x Neta+1)

    ## Apply stagnation refinement to chi (clusters near chi=1)
    if s["curvilinear_mapping"]["refinement_stagnation"]["state"]
        chi_1d = REFINEMENT(chi_1d,
            s["curvilinear_mapping"]["refinement_stagnation"]["BL_thickness"],
            0,
            s["curvilinear_mapping"]["refinement_stagnation"]["intensity"],
            1, 0)
    end

    ## Store wall element-space coordinates
    s["mesh"]["chi_wall"] = chi_1d
    s["mesh"]["eta_wall"] = zeros(s["mesh"]["Nchi"]+1)

    if s["restart"] && s["shock"]["enabled"]
        ## Restart with shock: extrude from wall along shock-aligned normals

        # Get wall positions from transform
        (s["mesh"]["x_wall"], s["mesh"]["y_wall"]) = GO_TO_PHYSICAL_SPACE(
            chi_1d, zeros(s["mesh"]["Nchi"]+1), s)

        # Get shock positions from spline at wall chi values
        shocked_points_x = s["shock"]["spline_func_x"].(chi_1d)
        shocked_points_y = s["shock"]["spline_func_y"].(chi_1d)

        # Compute shock-wall distances and reoriented normals
        dx_shock = shocked_points_x .- s["mesh"]["x_wall"]
        dy_shock = shocked_points_y .- s["mesh"]["y_wall"]
        shocked_r = sqrt.(dx_shock.^2 .+ dy_shock.^2)
        s["mesh"]["x_wall_normal"] = dx_shock ./ shocked_r
        s["mesh"]["y_wall_normal"] = dy_shock ./ shocked_r

        # Compute wall-normal distances with optional refinement
        index = round(Int, (s["mesh"]["Neta"]+1) / s["shock"]["remesh_shock_distance"])
        wall_distance = zeros(s["mesh"]["Nchi"]+1, s["mesh"]["Neta"]+1)

        for i = 1:s["mesh"]["Nchi"]+1
            wall_distance_line = collect(0:s["mesh"]["Neta"]) .* shocked_r[i] / (index - 0.5)

            if s["curvilinear_mapping"]["refinement_wall"]["state"]
                wall_distance_line = REFINEMENT(wall_distance_line,
                    s["curvilinear_mapping"]["refinement_wall"]["BL_thickness"],
                    0,
                    s["curvilinear_mapping"]["refinement_wall"]["intensity"],
                    maximum(wall_distance_line), 0)

                # Rescale so shock passes through cell center at index
                shocked_cell_center = (wall_distance_line[index] +
                    wall_distance_line[index+1]) / 2
                factor = shocked_r[i] / shocked_cell_center
                wall_distance_line = factor .* wall_distance_line
            end

            wall_distance[i,:] = wall_distance_line
        end

        # Compute corners from wall + normal extrusion
        s["mesh"]["x_corner"] = s["mesh"]["x_wall"] .+ wall_distance .* s["mesh"]["x_wall_normal"]
        s["mesh"]["y_corner"] = s["mesh"]["y_wall"] .+ wall_distance .* s["mesh"]["y_wall_normal"]

    else
        ## Non-restart: map element space to physical space via transform
        chi_2D = repeat(chi_1d, 1, s["mesh"]["Neta"]+1)
        eta_2D = repeat(eta_1d, s["mesh"]["Nchi"]+1, 1)
        (s["mesh"]["x_corner"], s["mesh"]["y_corner"]) = GO_TO_PHYSICAL_SPACE(chi_2D, eta_2D, s)

        # Extract wall coordinates
        s["mesh"]["x_wall"] = s["mesh"]["x_corner"][:, 1]
        s["mesh"]["y_wall"] = s["mesh"]["y_corner"][:, 1]

        # Compute wall normals from transform direction (eta=0 to eta=delta)
        dx = s["mesh"]["x_corner"][:, 2] .- s["mesh"]["x_corner"][:, 1]
        dy = s["mesh"]["y_corner"][:, 2] .- s["mesh"]["y_corner"][:, 1]
        norm_len = sqrt.(dx.^2 .+ dy.^2)
        s["mesh"]["x_wall_normal"] = dx ./ norm_len
        s["mesh"]["y_wall_normal"] = dy ./ norm_len
    end

    ## Make cells upstream from shock parallel to y = 0
    if s["shock"]["enabled"]
        index = round(Int, (s["mesh"]["Neta"]+1) / s["shock"]["remesh_shock_distance"])
        if disturbances
            s["mesh"]["y_corner"][:, index:end] = repeat(s["mesh"]["y_corner"][:, index],
                1, size(s["mesh"]["y_corner"], 2) - index + 1)
            s["mesh"]["x_corner"][:, index:end] = s["mesh"]["x_corner"][end, index:end]' .+
                s["mesh"]["x_corner"][:, index] .- s["mesh"]["x_corner"][end, index]
        end
    end

    ## Compute cell centers and extended grid
    s = GET_EXTENDED_CELL_CENTERS(s)

    s["mesh"]["arc_3D_Dtheta_Dr"] = 0.001  # Arc angle for 3D axisymmetric wedge

    ## Mesh smoothing for curvature discontinuities
    if haskey(s["curvilinear_mapping"], "smooth_mesh") &&
            s["curvilinear_mapping"]["smooth_mesh"] && !s["restart"]
        for iter = 1:10
            s = GET_CELL_CORNERS_SMOOTH(s)
            s = GET_EXTENDED_CELL_CENTERS(s)
        end
    end

    ## Compute left-right face areas and normals
    Nchi = s["mesh"]["Nchi"]::Int
    Neta = s["mesh"]["Neta"]::Int

    # Extract array references to avoid repeated dict lookups in loops
    xc = s["mesh"]["x_corner"]::Matrix{Float64}
    yc = s["mesh"]["y_corner"]::Matrix{Float64}
    arc_3D = s["mesh"]["arc_3D_Dtheta_Dr"]::Float64

    # Pre-allocate lr arrays
    lr_area     = zeros(Nchi+1, Neta)
    lr_x_normal = zeros(Nchi+1, Neta)
    lr_y_normal = zeros(Nchi+1, Neta)

    if s["PDE_dimension"] == "2D"
        @inbounds for i = 1:Nchi+1
            for j = 1:Neta
                lr_area[i,j] = sqrt((xc[i,j+1] - xc[i,j])^2 + (yc[i,j+1] - yc[i,j])^2)
                lr_x_normal[i,j] = -(yc[i,j+1] - yc[i,j]) / lr_area[i,j]
                lr_y_normal[i,j] =  (xc[i,j+1] - xc[i,j]) / lr_area[i,j]
                if i < Nchi+1
                    sgn = sign(lr_x_normal[i,j] * (xc[i+1,j] - xc[i,j]) + lr_y_normal[i,j] * (yc[i+1,j] - yc[i,j]))
                else
                    sgn = sign(lr_x_normal[i,j] * (xc[i,j] - xc[i-1,j]) + lr_y_normal[i,j] * (yc[i,j] - yc[i-1,j]))
                end
                lr_x_normal[i,j] *= sgn
                lr_y_normal[i,j] *= sgn
            end
        end
    elseif s["PDE_dimension"] == "3D-axisymmetric"
        @inbounds for i = 1:Nchi+1
            for j = 1:Neta
                depth = abs(yc[i,j+1] + yc[i,j]) / 2 * arc_3D
                edge_len = sqrt((xc[i,j+1] - xc[i,j])^2 + (yc[i,j+1] - yc[i,j])^2)
                lr_area[i,j] = edge_len * depth
                lr_x_normal[i,j] = -(yc[i,j+1] - yc[i,j]) / edge_len
                lr_y_normal[i,j] =  (xc[i,j+1] - xc[i,j]) / edge_len
                if i < Nchi+1
                    sgn = sign(lr_x_normal[i,j] * (xc[i+1,j] - xc[i,j]) + lr_y_normal[i,j] * (yc[i+1,j] - yc[i,j]))
                else
                    sgn = sign(lr_x_normal[i,j] * (xc[i,j] - xc[i-1,j]) + lr_y_normal[i,j] * (yc[i,j] - yc[i-1,j]))
                end
                lr_x_normal[i,j] *= sgn
                lr_y_normal[i,j] *= sgn
            end
        end
    else
        error("Dimension selection not recognized")
    end
    s["mesh"]["lr_area"]     = lr_area
    s["mesh"]["lr_x_normal"] = lr_x_normal
    s["mesh"]["lr_y_normal"] = lr_y_normal

    ## Compute bottom-top face areas and normals
    # Pre-allocate bt arrays
    bt_area     = zeros(Nchi, Neta+1)
    bt_x_normal = zeros(Nchi, Neta+1)
    bt_y_normal = zeros(Nchi, Neta+1)

    if s["PDE_dimension"] == "2D"
        @inbounds for i = 1:Nchi
            for j = 1:Neta+1
                bt_area[i,j] = sqrt((xc[i+1,j] - xc[i,j])^2 + (yc[i+1,j] - yc[i,j])^2)
                bt_x_normal[i,j] =  (yc[i+1,j] - yc[i,j]) / bt_area[i,j]
                bt_y_normal[i,j] = -(xc[i+1,j] - xc[i,j]) / bt_area[i,j]
                if j < Neta+1
                    sgn = sign(bt_x_normal[i,j] * (xc[i,j+1] - xc[i,j]) + bt_y_normal[i,j] * (yc[i,j+1] - yc[i,j]))
                else
                    sgn = sign(bt_x_normal[i,j] * (xc[i,j] - xc[i,j-1]) + bt_y_normal[i,j] * (yc[i,j] - yc[i,j-1]))
                end
                bt_x_normal[i,j] *= sgn
                bt_y_normal[i,j] *= sgn
            end
        end
    elseif s["PDE_dimension"] == "3D-axisymmetric"
        @inbounds for i = 1:Nchi
            for j = 1:Neta+1
                depth = abs(yc[i+1,j] + yc[i,j]) / 2 * arc_3D
                edge_len = sqrt((xc[i+1,j] - xc[i,j])^2 + (yc[i+1,j] - yc[i,j])^2)
                bt_area[i,j] = edge_len * depth
                bt_x_normal[i,j] =  (yc[i+1,j] - yc[i,j]) / edge_len
                bt_y_normal[i,j] = -(xc[i+1,j] - xc[i,j]) / edge_len
                if j < Neta+1
                    sgn = sign(bt_x_normal[i,j] * (xc[i,j+1] - xc[i,j]) + bt_y_normal[i,j] * (yc[i,j+1] - yc[i,j]))
                else
                    sgn = sign(bt_x_normal[i,j] * (xc[i,j] - xc[i,j-1]) + bt_y_normal[i,j] * (yc[i,j] - yc[i,j-1]))
                end
                bt_x_normal[i,j] *= sgn
                bt_y_normal[i,j] *= sgn
            end
        end
    else
        error("Dimension selection not recognized")
    end
    s["mesh"]["bt_area"]     = bt_area
    s["mesh"]["bt_x_normal"] = bt_x_normal
    s["mesh"]["bt_y_normal"] = bt_y_normal

    ## Compute front-back face areas (3D-axisymmetric only)
    if s["PDE_dimension"] == "3D-axisymmetric"
        fb_area = zeros(Nchi, Neta)
        @inbounds for i = 1:Nchi
            for j = 1:Neta
                # Front-back area equals the 2D cell area (Shoelace formula)
                fb = -xc[i,j]     * yc[i+1,j]   + xc[i+1,j]   * yc[i,j]
                fb = fb - xc[i+1,j]   * yc[i+1,j+1] + xc[i+1,j+1] * yc[i+1,j]
                fb = fb - xc[i+1,j+1] * yc[i,j+1]   + xc[i,j+1]   * yc[i+1,j+1]
                fb = fb - xc[i,j+1]   * yc[i,j]     + xc[i,j]     * yc[i,j+1]
                fb_area[i,j] = abs(fb / 2)
            end
        end
        s["mesh"]["fb_area"] = fb_area
        theta_angle_fb_faces = atan(arc_3D / 2)
        s["mesh"]["fb_x_normal_front"] =  0
        s["mesh"]["fb_y_normal_front"] =  sin(theta_angle_fb_faces)
        s["mesh"]["fb_x_normal_back"]  =  0
        s["mesh"]["fb_y_normal_back"]  = -sin(theta_angle_fb_faces)
    end

    ## Compute cell volumes
    volume = zeros(Nchi, Neta)

    if s["PDE_dimension"] == "2D"
        @inbounds for i = 1:Nchi
            for j = 1:Neta
                # Shoelace formula for quadrilateral area
                vol = -xc[i,j]     * yc[i+1,j]   + xc[i+1,j]   * yc[i,j]
                vol = vol - xc[i+1,j]   * yc[i+1,j+1] + xc[i+1,j+1] * yc[i+1,j]
                vol = vol - xc[i+1,j+1] * yc[i,j+1]   + xc[i,j+1]   * yc[i+1,j+1]
                vol = vol - xc[i,j+1]   * yc[i,j]     + xc[i,j]     * yc[i,j+1]
                volume[i,j] = abs(vol / 2)
            end
        end
    elseif s["PDE_dimension"] == "3D-axisymmetric"
        @inbounds for i = 1:Nchi
            for j = 1:Neta
                # Compute corners of wedge volume element
                c1 = [xc[i,j],     yc[i,j],      yc[i,j]/2     * arc_3D]
                c2 = [xc[i+1,j],   yc[i+1,j],    yc[i+1,j]/2   * arc_3D]
                c3 = [xc[i+1,j+1], yc[i+1,j+1],  yc[i+1,j+1]/2 * arc_3D]
                c4 = [xc[i,j+1],   yc[i,j+1],    yc[i,j+1]/2   * arc_3D]
                c5 = [xc[i,j],     yc[i,j],     -yc[i,j]/2     * arc_3D]
                c6 = [xc[i+1,j],   yc[i+1,j],   -yc[i+1,j]/2   * arc_3D]
                c7 = [xc[i+1,j+1], yc[i+1,j+1], -yc[i+1,j+1]/2 * arc_3D]
                c8 = [xc[i,j+1],   yc[i,j+1],   -yc[i,j+1]/2   * arc_3D]

                # Decompose hexahedron into six tetrahedra and sum volumes
                tet1 = tetrahedronVolume(c1, c2, c4, c6)
                tet2 = tetrahedronVolume(c2, c6, c4, c3)
                tet3 = tetrahedronVolume(c6, c4, c3, c7)
                tet4 = tetrahedronVolume(c1, c4, c8, c6)
                tet5 = tetrahedronVolume(c8, c4, c6, c7)
                tet6 = tetrahedronVolume(c5, c1, c6, c8)
                volume[i,j] = tet1 + tet2 + tet3 + tet4 + tet5 + tet6
            end
        end
    else
        error("Dimension selection not recognized")
    end
    s["mesh"]["volume"] = volume

    ## Compute geometry factors for CFL estimation
    x_Ext = s["mesh"]["x_Ext"]::Matrix{Float64}
    y_Ext = s["mesh"]["y_Ext"]::Matrix{Float64}
    @views begin
    max_1 = maximum(bt_area[:, 2:end] ./ volume)
    max_2 = maximum(lr_area[2:end, :] ./ volume)
    centroids_bt_distance = sqrt.((y_Ext[2:end-1, 2:end] .- y_Ext[2:end-1, 1:end-1]).^2 .+
                                  (x_Ext[2:end-1, 2:end] .- x_Ext[2:end-1, 1:end-1]).^2)
    centroids_lr_distance = sqrt.((y_Ext[2:end, 2:end-1] .- y_Ext[1:end-1, 2:end-1]).^2 .+
                                  (x_Ext[2:end, 2:end-1] .- x_Ext[1:end-1, 2:end-1]).^2)
    max_3 = maximum(bt_area[:, 2:end] ./ volume ./ centroids_bt_distance[:, 2:end])
    max_4 = maximum(lr_area[2:end, :] ./ volume ./ centroids_lr_distance[2:end, :])
    end # @views
    s["mesh"]["geometry_factor"] = max(max_1, max_2)
    s["mesh"]["volume_factor"] = max(max_3, max_4)

    return s
end


# tetrahedronVolume  Volume of a tetrahedron from four vertex coordinates.
#
#   volume = tetrahedronVolume(p1, p2, p3, p4)
#
#   Computes the volume using the scalar triple product:
#       V = |v1 . (v2 x v3)| / 6
#   where v1, v2, v3 are edge vectors from p1 to the other vertices.

function tetrahedronVolume(p1::Vector{Float64}, p2::Vector{Float64},
                           p3::Vector{Float64}, p4::Vector{Float64})
    v1 = p2 .- p1
    v2 = p3 .- p1
    v3 = p4 .- p1

    volume = abs(dot(v1, cross(v2, v3))) / 6

    return volume
end
