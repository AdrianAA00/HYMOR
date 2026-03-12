# DETECT_CELLS_SHOCKED - Identify grid cells closest to the shock and build cell masks.
#
#   s = DETECT_CELLS_SHOCKED(s)
#
#   For each streamwise station, computes the Euclidean distance from the
#   shock point to every grid node in the normal direction, selects the
#   closest cell as the shocked cell, and constructs three logical masks:
#     - shocked_cells : marks exactly the shocked cell at each station
#     - flow_cells    : marks all cells downstream of the shock (body side)
#     - flow_cells_E  : marks flow cells plus the shocked cell itself
#
#   Inputs:
#       s (Dict) - Solution structure containing at minimum:
#                    s["mesh"]["Nchi"], s["mesh"]["Neta"]  - Grid dimensions
#                    s["shock"]["points_x"]/["points_y"]  - (Nx x 1) shock point coordinates
#                    s["mesh"]["x"], s["mesh"]["y"]       - (Nx x Ny) grid coordinates
#                    s["shock"]["cells"]                  - (Nx x Ny) preallocated mask
#                    s["shock"]["flow_cells"]             - (Nx x Ny) preallocated mask
#                    s["shock"]["flow_cells_E"]           - (Nx x Ny) preallocated mask
#
#   Outputs:
#       s (Dict) - Updated s with:
#                    s["shock"]["cell_indices"]     - (Nx x 1) column index of shocked cell
#                    s["shock"]["cells"]            - (Nx x Ny) binary mask (1 at shock)
#                    s["shock"]["flow_cells"]       - (Nx x Ny) binary mask (1 for body-side cells)
#                    s["shock"]["flow_cells_E"]     - (Nx x Ny) binary mask (flow_cells + shocked cell)
#
#   Notes:
#       - Distance is computed using broadcasting, avoiding allocation of
#         full ones() matrices for broadcasting shock points.
#
# Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

function DETECT_CELLS_SHOCKED(s::Dict{String, Any})
    ## Extract dict refs
    shock = s["shock"]::Dict{String, Any}
    mesh  = s["mesh"]::Dict{String, Any}
    Nchi  = mesh["Nchi"]::Int
    Neta  = mesh["Neta"]::Int

    pts_x = shock["points_x"]
    pts_y = shock["points_y"]
    x     = mesh["x"]::Matrix{Float64}
    y     = mesh["y"]::Matrix{Float64}

    ## Find closest cell per row without allocating Nchi×Neta distance array.
    ## Uses squared distance (no sqrt needed for argmin).
    ## Replaces slow mapslices with direct loop.
    cell_indices = shock["cell_indices"]
    @inbounds for i in 1:Nchi
        px_i = Float64(pts_x[i])
        py_i = Float64(pts_y[i])
        min_dist_sq = Inf
        min_j = 1
        for j in 1:Neta
            dx = px_i - x[i, j]
            dy = py_i - y[i, j]
            d_sq = dx * dx + dy * dy
            if d_sq < min_dist_sq
                min_dist_sq = d_sq
                min_j = j
            end
        end
        cell_indices[i, 1] = min_j
    end

    ## Build cell masks - zero out using fill! instead of allocating new arrays
    cells       = shock["cells"]::Matrix{Float64}
    flow_cells  = shock["flow_cells"]::Matrix{Float64}
    flow_cells_E = shock["flow_cells_E"]::Matrix{Float64}

    fill!(cells, 0.0)
    fill!(flow_cells, 0.0)
    fill!(flow_cells_E, 0.0)

    @inbounds for i in 1:Nchi
        ci = cell_indices[i, 1]
        cells[i, ci] = 1.0
        for j in 1:ci - 1
            flow_cells[i, j] = 1.0
        end
        for j in 1:ci
            flow_cells_E[i, j] = 1.0
        end
    end

    return s
end
