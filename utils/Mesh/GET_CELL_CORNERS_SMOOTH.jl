# GET_CELL_CORNERS_SMOOTH  Smooth mesh corners via weighted averaging.
#
#   s = GET_CELL_CORNERS_SMOOTH(s)
#
#   Recomputes interior cell corner positions as a weighted average of the
#   four surrounding extended cell centers. An exponential weight function
#   preserves the original wall boundary coordinates (j = 1 row) while
#   increasingly smoothing corners farther from the wall.
#
#   Inputs:
#       s - Dict containing at minimum:
#           s["mesh"]["Nchi"], s["mesh"]["Neta"]         - Number of cells in x and y
#           s["mesh"]["x_corner"], s["mesh"]["y_corner"] - (Nx+1 x Ny+1) corner coordinates
#           s["mesh"]["x_Ext"], s["mesh"]["y_Ext"]       - (Nx+2 x Ny+2) extended cell centers
#
#   Outputs:
#       s - Updated Dict with smoothed s["mesh"]["x_corner"] and s["mesh"]["y_corner"]
#           (wall row j=1 is unchanged).
#
#   Notes:
#       - The weight function is w = 1 - exp(-10*j/Ny), which is nearly
#         zero for j close to the wall and approaches 1 away from it.
#       - This function is typically called iteratively (e.g., 10 times)
#         in combination with GET_EXTENDED_CELL_CENTERS to progressively
#         smooth meshes with curvature discontinuities (MSL, blunt_cone).
#
# Part of: Hypersonics Stability Julia Solver - Mesh Module

function GET_CELL_CORNERS_SMOOTH(s::Dict{String,Any})

    ## Smooth interior corners via weighted averaging
    Nchi = s["mesh"]["Nchi"]::Int
    Neta = s["mesh"]["Neta"]::Int
    x_corner = s["mesh"]["x_corner"]::Matrix{Float64}
    y_corner = s["mesh"]["y_corner"]::Matrix{Float64}
    x_Ext = s["mesh"]["x_Ext"]::Matrix{Float64}
    y_Ext = s["mesh"]["y_Ext"]::Matrix{Float64}

    @inbounds for i = 1:Nchi+1
        for j = 2:Neta+1
            weight = 1 - exp(-10*j/Neta)
            x_corner[i,j] = (x_Ext[i,j] + x_Ext[i+1,j] + x_Ext[i,j+1] + x_Ext[i+1,j+1])/4 * weight +
                              x_corner[i,j] * (1-weight)
            y_corner[i,j] = (y_Ext[i,j] + y_Ext[i+1,j] + y_Ext[i,j+1] + y_Ext[i+1,j+1])/4 * weight +
                              y_corner[i,j] * (1-weight)
        end
    end

    return s
end
