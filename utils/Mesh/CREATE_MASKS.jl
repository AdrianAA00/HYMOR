# CREATE_MASKS - Generate spatial masks for stagnation and boundary-layer regions.
#
#   (mask_no_stagnation, mask_no_BL, mask_only_BL) = CREATE_MASKS(s)
#
#   Creates smooth sigmoid-based masks that separate the computational
#   domain into stagnation region, boundary-layer region, and outer flow
#   region. These masks are used for selective application of numerical
#   treatments (e.g., artificial dissipation, limiting) in different
#   parts of the flow field.
#
#   Inputs:
#       s  (Dict{String,Any}) - Solution structure containing:
#           s["curvilinear_mapping"]["boundary_type"] - Geometry type ("MSL", "blunt_cone", etc.)
#           s["curvilinear_mapping"]["L"]             - Reference body length
#           s["curvilinear_mapping"]["R"]             - Nose radius
#           s["curvilinear_mapping"]["theta"]         - Half-cone angle [rad]
#           s["mesh"]["chi"]           - (Nchi x Neta) element-space chi coordinates
#           s["mesh"]["Nchi"], s["mesh"]["Neta"] - Number of cells in each direction
#           s["var"]["T"]              - (Nchi+2 x Neta+2) temperature field with ghosts
#           s["shock"]["flow_cells"]   - (Nchi x Neta) logical mask of active flow cells
#
#   Outputs:
#       mask_no_stagnation  (Nx x 1)  - Sigmoid mask that is ~0 near the
#                                        stagnation point and ~1 elsewhere
#       mask_no_BL          (Nx x Ny) - Sigmoid mask that is ~0 inside the
#                                        thermal boundary layer and ~1 outside
#       mask_only_BL        (Nx x Ny) - Complement of mask_no_BL (1 - mask_no_BL)
#
#   Notes:
#       - For "MSL" and "blunt_cone" geometries, the stagnation mask uses
#         the surface coordinate chi and a logistic sigmoid function.
#       - The boundary-layer mask is based on the thermal field: cells
#         with temperature below 80% of the mean are considered outside
#         the boundary layer.
#       - For all other geometry types, trivial masks are returned (no
#         stagnation exclusion, no BL exclusion).
#
# Part of: Hypersonics Stability Julia Solver - Mesh Module

using Statistics

function CREATE_MASKS(s::Dict{String,Any})

    if s["curvilinear_mapping"]["boundary_type"] == "MSL" || s["curvilinear_mapping"]["boundary_type"] == "blunt_cone"
        ## Geometry parameters
        s_Straight = ((s["curvilinear_mapping"]["L"] - s["curvilinear_mapping"]["R"]) + s["curvilinear_mapping"]["R"] * cos(s["curvilinear_mapping"]["theta"])) / sin(s["curvilinear_mapping"]["theta"])
        s_Curve = s["curvilinear_mapping"]["theta"] * s["curvilinear_mapping"]["R"]
        s_tot = s_Curve + s_Straight

        ## Stagnation mask (sigmoid along surface coordinate)
        amount = 0.5
        len = 0.01
        chi = s["mesh"]["chi"]
        @. coordinate = (1 - 2 * chi - amount * s_Curve / s_tot) / len
        @. mask_no_stagnation = exp(coordinate) / (1 + exp(coordinate))

        ## Boundary-layer mask (sigmoid based on thermal field)
        len = 0.02
        flow_cells = s["shock"]["flow_cells"]
        @views T_field = s["var"]["T"][2:end-1, 2:end-1] .* flow_cells
        mean_T = mean(T_field[T_field .!= 0])
        @. coordinate = (T_field - mean_T * 0.8) / mean_T / len
        @. mask_no_BL = exp(coordinate) / (1 + exp(coordinate)) * flow_cells +
                     1 - flow_cells

        @. mask_only_BL = 1 - mask_no_BL
    else
        ## Trivial masks for non-blunt geometries
        Nchi = s["mesh"]["Nchi"]::Int
        Neta = s["mesh"]["Neta"]::Int
        mask_no_stagnation = ones(Nchi, Neta)
        mask_only_BL = zeros(Nchi, Neta)
        mask_no_BL = ones(Nchi, Neta)
    end

    return mask_no_stagnation, mask_no_BL, mask_only_BL
end
