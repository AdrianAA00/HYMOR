function PLOT_INITIAL_SET_UP(s)
% PLOT_INITIAL_SET_UP  Visualize the computational mesh, shocked cells,
%   cell corners, and (optionally) the fitted shock spline.
%
%   PLOT_INITIAL_SET_UP(s)
%
%   Inputs:
%       s - Solution structure containing:
%                    .mesh.x_Ext, .mesh.y_Ext        - Extended cell-centre coordinates.
%                    .mesh.x, .mesh.y                - Interior cell-centre coordinates.
%                    .mesh.x_corner, .mesh.y_corner  - Cell corner coordinates.
%                    .shock.cells            - Logical mask for shocked cells.
%                    .shock.enabled          - Logical flag for shock presence.
%                    .shock.points_chi       - Chi values at shock points.
%                    .shock.spline_func_x    - Piecewise polynomial (x vs chi).
%                    .shock.spline_func_y    - Piecewise polynomial (y vs chi).
%                    .shock.flow_cells       - (set to 1 if no shock present).
%
%   Outputs:
%       (none) - Produces figure 100 with the mesh scatter plot, shocked
%                cells highlighted in red, grid lines in blue, and the
%                shock spline in black. Calls CHECK_MESH for validation.
%
%   Notes:
%       This function is intended to be called once after mesh generation
%       to verify grid quality before running the solver.
%
% Part of: Hypersonics Stability MATLAB Solver - Postprocessing Module

    %% Prepare shock spline for plotting
    if s.shock.enabled == true
        N_plot_shock_poly = 1000;
        chi_p = linspace(s.shock.points_chi(1,1), s.shock.points_chi(end,1), N_plot_shock_poly);
        x_p = ppval(s.shock.spline_func_x, chi_p);
        y_p = ppval(s.shock.spline_func_y, chi_p);
    end

    %% Plot mesh, shocked cells, grid lines, and shock
    figure(100)
    scatter(s.mesh.x_Ext, s.mesh.y_Ext, 10, "black")
    hold on
    scatter(s.mesh.x .* s.shock.cells, s.mesh.y .* s.shock.cells, 10, "red", "o", "filled")
    hold on
    plot(s.mesh.x_corner, s.mesh.y_corner, "blue")
    hold on
    plot(s.mesh.x_corner', s.mesh.y_corner', "blue")
    if s.shock.enabled == true
        hold on
        plot(x_p, y_p, "black", LineWidth=2)
    end
    axis equal
    hold off

    %% Handle non-shock case and validate mesh
    if ~s.shock.enabled
        s.shock.flow_cells = 1;
    end

    CHECK_MESH(s);
end
