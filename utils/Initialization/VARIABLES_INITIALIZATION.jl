# VARIABLES_INITIALIZATION  Pre-allocate all field arrays in the solution struct.
#   Zeros out coordinate arrays, mesh metric arrays, conserved-variable
#   fields (with ghost cells), chemistry fields, and shock-tracking arrays
#   based on the grid dimensions Nchi and Neta stored in the solution struct.
#
#   s = VARIABLES_INITIALIZATION(s)
#
#   Inputs:
#       s - (Dict) Solution struct containing at minimum the grid
#           dimensions s["mesh"]["Nchi"] and s["mesh"]["Neta"]
#
#   Outputs:
#       s - (Dict) Solution struct with all field arrays initialized to
#           zero at the appropriate sizes
#
#   Notes:
#       - Wall arrays have size (Nchi+1, 1) corresponding to cell edges.
#       - Interior cell-center arrays have size (Nchi, Neta).
#       - Extended arrays (with ghost cells) have size (Nchi+2, Neta+2).
#       - Cell-corner arrays have size (Nchi+1, Neta+1).
#       - Face-area and normal arrays are sized for left-right (Nchi+1, Neta),
#         bottom-top (Nchi, Neta+1), and front-back (Nchi, Neta) faces.
#       - Shock arrays are sized (Nchi, 1) for per-column shock data.
#
# Part of: Hypersonics Stability MATLAB Solver - Initialization Module

function VARIABLES_INITIALIZATION(s::Dict{String,Any})

    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]

    # Ensure nested dicts exist
    if !haskey(s, "mesh")
        s["mesh"] = Dict{String,Any}()
    end
    if !haskey(s, "var")
        s["var"] = Dict{String,Any}()
    end
    if !haskey(s, "flux")
        s["flux"] = Dict{String,Any}()
    end
    if !haskey(s, "shock")
        s["shock"] = Dict{String,Any}()
    end

    ## Wall coordinates and normals
    s["mesh"]["x_wall"]        = zeros(Nchi + 1, 1)
    s["mesh"]["y_wall"]        = zeros(Nchi + 1, 1)
    s["mesh"]["eta_wall"]      = zeros(Nchi + 1, 1)
    s["mesh"]["chi_wall"]      = zeros(Nchi + 1, 1)
    s["mesh"]["x_wall_normal"] = zeros(Nchi + 1, 1)
    s["mesh"]["y_wall_normal"] = zeros(Nchi + 1, 1)

    ## Cell-center coordinates (interior)
    s["mesh"]["x"]   = zeros(Nchi, Neta)
    s["mesh"]["y"]   = zeros(Nchi, Neta)
    s["mesh"]["eta"] = zeros(Nchi, Neta)
    s["mesh"]["chi"] = zeros(Nchi, Neta)

    ## Cell-center coordinates (extended with ghost cells)
    s["mesh"]["x_Ext"] = zeros(Nchi + 2, Neta + 2)
    s["mesh"]["y_Ext"] = zeros(Nchi + 2, Neta + 2)

    ## Cell-corner coordinates
    s["mesh"]["x_corner"] = zeros(Nchi + 1, Neta + 1)
    s["mesh"]["y_corner"] = zeros(Nchi + 1, Neta + 1)

    ## Cell-face areas
    s["mesh"]["lr_area"] = zeros(Nchi + 1, Neta)     # left-right faces
    s["mesh"]["bt_area"] = zeros(Nchi, Neta + 1)     # bottom-top faces
    s["mesh"]["fb_area"] = zeros(Nchi, Neta)          # front-back faces (3D axisymmetric)

    ## Cell-face normals
    s["mesh"]["lr_x_normal"] = zeros(Nchi + 1, Neta)  # left-right x-component
    s["mesh"]["lr_y_normal"] = zeros(Nchi + 1, Neta)  # left-right y-component
    s["mesh"]["bt_x_normal"] = zeros(Nchi, Neta + 1)  # bottom-top x-component
    s["mesh"]["bt_y_normal"] = zeros(Nchi, Neta + 1)  # bottom-top y-component

    ## Cell volumes
    s["mesh"]["volume"] = zeros(Nchi, Neta)

    ## Conserved variables (extended with ghost cells for BC)
    s["var"]["rho"]   = zeros(Nchi + 2, Neta + 2)
    s["var"]["rho_u"] = zeros(Nchi + 2, Neta + 2)
    s["var"]["rho_v"] = zeros(Nchi + 2, Neta + 2)
    s["var"]["rho_E"] = zeros(Nchi + 2, Neta + 2)
    s["var"]["p"]     = zeros(Nchi + 2, Neta + 2)

    ## Conserved-variable fluxes (interior cells only)
    s["flux"]["rho"]   = zeros(Nchi, Neta)
    s["flux"]["rho_u"] = zeros(Nchi, Neta)
    s["flux"]["rho_v"] = zeros(Nchi, Neta)
    s["flux"]["rho_E"] = zeros(Nchi, Neta)

    ## Chemistry field variables (extended with ghost cells)
    s["var"]["gamma_star"]    = zeros(Nchi + 2, Neta + 2)
    s["var"]["gamma_star_eq"] = zeros(Nchi + 2, Neta + 2)
    s["flux"]["gamma_star"]   = zeros(Nchi, Neta)
    s["var"]["cv_star"]       = zeros(Nchi + 2, Neta + 2)
    s["var"]["cv_star_eq"]    = zeros(Nchi + 2, Neta + 2)
    s["flux"]["cv_star"]      = zeros(Nchi, Neta)

    ## Transport property fields (extended with ghost cells)
    s["var"]["mu_star"]  = zeros(Nchi + 2, Neta + 2)
    s["var"]["k_star"]   = zeros(Nchi + 2, Neta + 2)
    s["var"]["Re_flow"]  = zeros(Nchi + 2, Neta + 2)
    s["var"]["Pr_flow"]  = zeros(Nchi + 2, Neta + 2)

    ## Shock-tracking arrays
    s["shock"]["properties"]   = Dict{String,Any}()
    s["shock"]["cell_indices"] = zeros(Int, Nchi, 1)
    s["shock"]["cells"]        = zeros(Nchi, Neta)
    s["shock"]["flow_cells"]   = ones(Nchi, Neta)
    s["shock"]["flow_cells_E"] = ones(Nchi, Neta)
    s["shock"]["speed_x"]      = zeros(Nchi, 1)
    s["shock"]["speed_y"]      = zeros(Nchi, 1)
    s["shock"]["points_x"]     = zeros(Nchi, 1)
    s["shock"]["points_y"]     = zeros(Nchi, 1)
    s["shock"]["points_eta"]   = zeros(Nchi, 1)
    s["shock"]["points_chi"]   = zeros(Nchi, 1)
    s["shock"]["beta"]         = zeros(Nchi, 1)
    s["shock"]["M"]            = zeros(Nchi, 1)

    return s
end
