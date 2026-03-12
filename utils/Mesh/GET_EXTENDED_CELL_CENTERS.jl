# GET_EXTENDED_CELL_CENTERS  Compute cell centers and build extended grid.
#
#   s = GET_EXTENDED_CELL_CENTERS(s)
#
#   Computes cell centers as the average of the four surrounding corner
#   coordinates, converts them to element space (chi, eta), and builds
#   an extended coordinate array (x_Ext, y_Ext) that includes one layer
#   of ghost cells on each boundary. Ghost cell positions are extrapolated
#   by linear reflection from the boundary face centers.
#
#   Inputs:
#       s - Dict containing at minimum:
#           s["mesh"]["Nchi"], s["mesh"]["Neta"]               - Number of cells in x and y
#           s["mesh"]["x_corner"], s["mesh"]["y_corner"]       - (Nx+1 x Ny+1) corner coordinates
#           s["curvilinear_mapping"]["boundary_type"]          - String for GO_TO_ELEMENT_SPACE
#
#   Outputs:
#       s - Updated Dict with added/modified fields:
#           s["mesh"]["x"], s["mesh"]["y"]           - (Nx x Ny) cell center coordinates
#           s["mesh"]["chi"], s["mesh"]["eta"]        - (Nx x Ny) element-space coordinates
#           s["mesh"]["x_Ext"], s["mesh"]["y_Ext"]    - (Nx+2 x Ny+2) extended cell centers
#                                                       including ghost cells
#
#   Notes:
#       - Interior cell centers occupy indices (2:end-1, 2:end-1) in the
#         extended arrays.
#       - Boundary face centers (midpoints of adjacent corners) are placed
#         at the extended array edges, then reflected outward to create
#         ghost cell positions for second-order boundary treatment.
#       - Element-space conversion is performed by GO_TO_ELEMENT_SPACE.
#
# Part of: Hypersonics Stability Julia Solver - Mesh Module

function GET_EXTENDED_CELL_CENTERS(s::Dict{String,Any})

    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]

    ## Compute interior cell centers
    x_corner = s["mesh"]["x_corner"]::Matrix{Float64}
    y_corner = s["mesh"]["y_corner"]::Matrix{Float64}
    x_center = zeros(Nchi, Neta)
    y_center = zeros(Nchi, Neta)
    @inbounds for i = 1:Nchi
        for j = 1:Neta
            x_center[i,j] = (x_corner[i,j] + x_corner[i+1,j] + x_corner[i,j+1] + x_corner[i+1,j+1]) / 4
            y_center[i,j] = (y_corner[i,j] + y_corner[i+1,j] + y_corner[i,j+1] + y_corner[i+1,j+1]) / 4
        end
    end
    s["mesh"]["x"] = x_center
    s["mesh"]["y"] = y_center

    ## Convert to element space
    (s["mesh"]["chi"], s["mesh"]["eta"]) = GO_TO_ELEMENT_SPACE(s["mesh"]["x"], s["mesh"]["y"], s)

    ## Populate interior of extended arrays
    x_Ext = zeros(Nchi+2, Neta+2)
    y_Ext = zeros(Nchi+2, Neta+2)

    @views begin
    x_Ext[2:end-1, 2:end-1] = s["mesh"]["x"]
    y_Ext[2:end-1, 2:end-1] = s["mesh"]["y"]

    ## Boundary face centers (midpoints of adjacent corners)
    @. x_Ext[2:end-1, 1]   = (x_corner[1:end-1, 1]   + x_corner[2:end, 1])   / 2  # Bottom boundary
    @. y_Ext[2:end-1, 1]   = (y_corner[1:end-1, 1]   + y_corner[2:end, 1])   / 2
    @. x_Ext[2:end-1, end] = (x_corner[1:end-1, end] + x_corner[2:end, end]) / 2  # Top boundary
    @. y_Ext[2:end-1, end] = (y_corner[1:end-1, end] + y_corner[2:end, end]) / 2
    @. x_Ext[1, 2:end-1]   = (x_corner[1, 1:end-1]  + x_corner[1, 2:end])   / 2  # Right boundary
    @. y_Ext[1, 2:end-1]   = (y_corner[1, 1:end-1]  + y_corner[1, 2:end])   / 2
    @. x_Ext[end, 2:end-1] = (x_corner[end, 1:end-1] + x_corner[end, 2:end]) / 2  # Left boundary
    @. y_Ext[end, 2:end-1] = (y_corner[end, 1:end-1] + y_corner[end, 2:end]) / 2

    ## Corner points of extended array
    x_Ext[1, 1]     = x_corner[1, 1]        # Bottom right corner
    y_Ext[1, 1]     = y_corner[1, 1]
    x_Ext[1, end]   = x_corner[1, end]      # Top right corner
    y_Ext[1, end]   = y_corner[1, end]
    x_Ext[end, 1]   = x_corner[end, 1]      # Bottom left corner
    y_Ext[end, 1]   = y_corner[end, 1]
    x_Ext[end, end] = x_corner[end, end]    # Top left corner
    y_Ext[end, end] = y_corner[end, end]

    ## Ghost cells (linear extrapolation from boundary face centers)
    @. x_Ext[:, 1]   = 2 * x_Ext[:, 1]   - x_Ext[:, 2]
    @. x_Ext[:, end] = 2 * x_Ext[:, end] - x_Ext[:, end-1]
    @. x_Ext[1, :]   = 2 * x_Ext[1, :]   - x_Ext[2, :]
    @. x_Ext[end, :] = 2 * x_Ext[end, :] - x_Ext[end-1, :]
    @. y_Ext[:, 1]   = 2 * y_Ext[:, 1]   - y_Ext[:, 2]
    @. y_Ext[:, end] = 2 * y_Ext[:, end] - y_Ext[:, end-1]
    @. y_Ext[1, :]   = 2 * y_Ext[1, :]   - y_Ext[2, :]
    @. y_Ext[end, :] = 2 * y_Ext[end, :] - y_Ext[end-1, :]
    end # @views

    s["mesh"]["x_Ext"] = x_Ext
    s["mesh"]["y_Ext"] = y_Ext

    return s
end
