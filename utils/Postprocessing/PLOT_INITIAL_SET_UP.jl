using Plots, Printf, LinearAlgebra, LoopVectorization

# PLOT_INITIAL_SET_UP  Visualize the computational mesh, shocked cells,
#   cell corners, and (optionally) the fitted shock spline.
#
#   PLOT_INITIAL_SET_UP(s)
#
#   Inputs:
#       s - Solution Dict containing:
#                    .mesh.x_Ext, .mesh.y_Ext        - Extended cell-centre coordinates.
#                    .mesh.x, .mesh.y                - Interior cell-centre coordinates.
#                    .mesh.x_corner, .mesh.y_corner  - Cell corner coordinates.
#                    .shock.cells            - Logical mask for shocked cells.
#                    .shock.enabled          - Logical flag for shock presence.
#                    .shock.points_chi       - Chi values at shock points.
#                    .shock.spline_func_x    - Piecewise polynomial (x vs chi).
#                    .shock.spline_func_y    - Piecewise polynomial (y vs chi).
#                    .shock.flow_cells       - (set to 1 if no shock present).
#
#   Outputs:
#       (none) - Produces figure with the mesh scatter plot, shocked
#                cells highlighted in red, grid lines in blue, and the
#                shock spline in black. Calls CHECK_MESH for validation.
#
#   Notes:
#       This function is intended to be called once after mesh generation
#       to verify grid quality before running the solver.
#
# Part of: Hypersonics Stability MATLAB Solver - Postprocessing Module

function PLOT_INITIAL_SET_UP(s::Dict{String, Any})
    SET_PLOT_DEFAULTS()

    ## Prepare shock spline for plotting
    if s["shock"]["enabled"] == true
        N_plot_shock_poly = 1000
        chi_p = range(s["shock"]["points_chi"][1, 1],
                      s["shock"]["points_chi"][end, 1],
                      length=N_plot_shock_poly)
        x_p = ppval(s["shock"]["spline_func_x"], collect(chi_p))
        y_p = ppval(s["shock"]["spline_func_y"], collect(chi_p))
    end

    ## Plot mesh, shocked cells, grid lines, and shock
    fig_size = auto_figure_size(s["mesh"]["x_Ext"], s["mesh"]["y_Ext"])
    fig, ax = PyPlot.subplots()
    fig.set_size_inches(fig_size[1] / 100.0, fig_size[2] / 100.0)

    ax.scatter(vec(s["mesh"]["x_Ext"]), vec(s["mesh"]["y_Ext"]),
               s=1, c="black", zorder=1)

    ax.scatter(vec(s["mesh"]["x"] .* s["shock"]["cells"]),
               vec(s["mesh"]["y"] .* s["shock"]["cells"]),
               s=1, c="red", zorder=2)

    # Plot grid lines along rows
    for i in 1:size(s["mesh"]["x_corner"], 1)
        ax.plot(s["mesh"]["x_corner"][i, :], s["mesh"]["y_corner"][i, :],
                color="blue", linewidth=0.5)
    end
    # Plot grid lines along columns
    for j in 1:size(s["mesh"]["x_corner"], 2)
        ax.plot(s["mesh"]["x_corner"][:, j], s["mesh"]["y_corner"][:, j],
                color="blue", linewidth=0.5)
    end

    if s["shock"]["enabled"] == true
        ax.plot(x_p, y_p, color="black", linewidth=2)
    end

    ax.set_aspect("equal")
    fig.tight_layout()
    save_or_display_figure(fig, name="initial_setup")

    ## Handle non-shock case and validate mesh
    if !s["shock"]["enabled"]
        s["shock"]["flow_cells"] = 1
    end

    CHECK_MESH(s)

    return nothing
end
