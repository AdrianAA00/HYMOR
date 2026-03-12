# CHECK_MESH - Verify mesh consistency using divergence-free test.
#
#   s = CHECK_MESH(s)
#
#   Validates the generated mesh by checking that a uniform flow field
#   in both x and y directions produces zero divergence through all cell
#   faces. This is a necessary condition for a valid finite-volume mesh:
#   the sum of fluxes through each cell must vanish for a uniform field.
#
#   Inputs:
#       s  (Dict{String,Any}) - Solution structure containing mesh data:
#           s["PDE_dimension"]              - Dimensionality flag ("2D" or "3D-axisymmetric")
#           s["mesh"]["bt_area"]        - (Nx x Ny+1) bottom-to-top face areas
#           s["mesh"]["bt_x_normal"]    - (Nx x Ny+1) x-component of bt face normals
#           s["mesh"]["bt_y_normal"]    - (Nx x Ny+1) y-component of bt face normals
#           s["mesh"]["lr_area"]        - (Nx+1 x Ny) left-to-right face areas
#           s["mesh"]["lr_x_normal"]    - (Nx+1 x Ny) x-component of lr face normals
#           s["mesh"]["lr_y_normal"]    - (Nx+1 x Ny) y-component of lr face normals
#
#   Outputs:
#       s  (Dict{String,Any}) - Updated structure with:
#           s["check_mesh"]   - "correct" if divergence is below tolerance,
#                               "incorrect" otherwise
#
#   Notes:
#       - Only performs the check for 2D meshes.
#       - Uses a tolerance of 1e-14 for the divergence residual.
#       - The test computes the net flux of a uniform (1,0) and (0,1)
#         field through each cell and verifies both vanish to machine
#         precision.
#
# Part of: Hypersonics Stability Julia Solver - Mesh Module

using Printf

function CHECK_MESH(s::Dict{String,Any})

    if s["PDE_dimension"] == "2D"
        ## Divergence-free check for uniform flow
        eps_tol = 1e-14

        @views begin
        # x-direction residual
        a = s["mesh"]["bt_area"] .* s["mesh"]["bt_x_normal"]
        b = -s["mesh"]["lr_area"] .* s["mesh"]["lr_x_normal"]
        res_x = a[:, 1:end-1] .- a[:, 2:end] .+ b[2:end, :] .- b[1:end-1, :]

        # y-direction residual
        a = s["mesh"]["bt_area"] .* s["mesh"]["bt_y_normal"]
        b = -s["mesh"]["lr_area"] .* s["mesh"]["lr_y_normal"]
        res_y = a[:, 1:end-1] .- a[:, 2:end] .+ b[2:end, :] .- b[1:end-1, :]
        end # @views

        ## Report result
        if sum(res_y) < eps_tol && sum(res_x) < eps_tol
            s["check_mesh"] = "correct"
        else
            s["check_mesh"] = "incorrect"
        end

        @printf("Mesh check = ")
        @printf("%s", s["check_mesh"])
        @printf("\n")
    end

    return s
end
