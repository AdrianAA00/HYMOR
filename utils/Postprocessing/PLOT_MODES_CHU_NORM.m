function PLOT_MODES_CHU_NORM(freestream_disturbances, A, s, chemistry, V, T_plot, w_infty)
% PLOT_MODES_CHU_NORM  Visualise energy-budget fields.
%
%   PLOT_MODES_CHU_NORM(freestream_disturbances, A, s, chemistry, V,
%        T_plot, w_infty)
%
%   Generates publication-quality 2-D pseudocolor plots of perturbation
%   quantities (pressure, density, velocity components, vorticity,
%   divergence, entropy, energy flux) and, when applicable, kinetic and
%   entropic energy-budget maps.  Optionally evolves the perturbation state
%   forward in time on the GPU using a Taylor-series approximation to the
%   matrix exponential.
%
%   Inputs:
%       freestream_disturbances - (logical) true when upstream disturbance
%                                 modes are included in the state vector V.
%       A                       - (sparse matrix) linearised operator used
%                                 for GPU time integration.
%       s                - (struct) base-flow s containing
%                                 grid coordinates, thermodynamic fields,
%                                 shock geometry, and solver parameters.
%       chemistry               - (struct) chemistry model data (needed for
%                                 entropy budgets and remeshing).
%       V                       - (column vector) perturbation state vector
%                                 [rho'; rho_u'; rho_v'; rho_E'; ...]; may
%                                 include shock displacement and freestream
%                                 mode amplitudes.
%       T_plot                  - (scalar) physical time at which to
%                                 evaluate. T_plot = 0 plots the initial
%                                 perturbation; T_plot > 0 triggers GPU
%                                 time integration.
%       w_infty                 - (vector) angular frequencies of the
%                                 freestream disturbance modes (empty when
%                                 not used).
%
%   Outputs:
%       (none) - all output is graphical (MATLAB figure windows).
%
%   Notes:
%       - Auxiliary ghost / shock cells are stripped via
%         REMOVE_AUX_SHOCK_CELLS before any processing.
%       - When T_plot == 0 and freestream_disturbances is true the s
%         is remeshed with a wider wall-normal extent so the freestream
%         region is visible.
%       - GPU time integration uses LINEAR_INTEGRATION_GPU (local
%         subfunction) with a user-specified Taylor-series order.
%       - Scaling factors are derived from either the energy-flux inflow
%         rate or the maximum total energy rate of change at t = 0.
%       - All contour plots are created by CREATE_PUBLICATION_PLOT.
%       - Figure export is available through EXPORT_ALL_FIGURES.
%
% Part of: Hypersonics Stability MATLAB Solver - Postprocessing Module
    if s.shock.enabled
        [s] = REMOVE_AUX_SHOCK_CELLS(s);
    end

    N_v = s.mesh.Nchi * s.mesh.Neta;
    if freestream_disturbances
        N_f = length(w_infty);
    else
        N_f = 0;
    end

    %% Extract perturbation fields from eigenvector
    pert_rho   = real(reshape(V((0)*N_v+1:1*N_v, 1), s.mesh.Nchi, s.mesh.Neta));
    pert_rho_u = real(reshape(V((1)*N_v+1:2*N_v, 1), s.mesh.Nchi, s.mesh.Neta));
    pert_rho_v = real(reshape(V((2)*N_v+1:3*N_v, 1), s.mesh.Nchi, s.mesh.Neta));
    pert_rho_E = real(reshape(V((3)*N_v+1:4*N_v, 1), s.mesh.Nchi, s.mesh.Neta));

    %% Extract freestream disturbance amplitudes (if present)
    size_freestream_problem_1 = 4 * s.mesh.Nchi * s.mesh.Neta + 4 * s.mesh.Nchi * N_f + s.mesh.Nchi;
    size_freestream_problem_2 = 4 * s.mesh.Nchi * s.mesh.Neta + 4 * s.mesh.Nchi * N_f;

    if size(V, 1) == size_freestream_problem_1
        idx_start_infty = 4 * s.mesh.Nchi * s.mesh.Neta + s.mesh.Nchi;
        pert_infty_all = V(idx_start_infty+1:size_freestream_problem_1, 1);
    elseif size(V, 1) == size_freestream_problem_2
        idx_start_infty = 4 * s.mesh.Nchi * s.mesh.Neta;
        pert_infty_all = V(idx_start_infty+1:size_freestream_problem_2, 1);
    end

    %% Compute non-dimensional scaling factors
    if s.shock.enabled && freestream_disturbances && T_plot > 0
        solution_remesh = false;
        [pert_rho, pert_rho_u, pert_rho_v, pert_rho_E] = ADD_FREESTREAM_DISTURBANCES( ...
            pert_rho, pert_rho_u, pert_rho_v, pert_rho_E, ...
            pert_infty_all, s, w_infty, T_plot, solution_remesh);
        [~, dE_dt_V_inflow] = ENERGY_INFLOW(V, s, w_infty);
        scale_kinetic  = dE_dt_V_inflow;
        scale_entropic = dE_dt_V_inflow;

    elseif s.shock.enabled && freestream_disturbances && T_plot == 0
        % Remesh to visualise freestream region
        s.shock.remesh_shock_distance = 8;
        s.mesh.Neta = s.mesh.Neta * 5;
        solution_old = s;
        disturbances = true;
        [s] = RESTART_SOLUTION(s, solution_old, chemistry, disturbances);
        [s] = REMOVE_AUX_SHOCK_CELLS(s);
        solution_remesh = true;
        [pert_rho, pert_rho_u, pert_rho_v, pert_rho_E] = ADD_FREESTREAM_DISTURBANCES( ...
            pert_rho, pert_rho_u, pert_rho_v, pert_rho_E, ...
            pert_infty_all, s, w_infty, 0, solution_remesh);
    else
        % No freestream disturbances -- scale by max total energy rate
        output_flow = true;
        budgets_kinetic = COMPUTE_BUDGETS_KINETIC(V, s, output_flow);
        temp_kinetic = budgets_kinetic.A_adv + budgets_kinetic.P_mom ...
            + budgets_kinetic.P_mass + budgets_kinetic.Dilat_P ...
            + budgets_kinetic.Transport_tau + budgets_kinetic.Transport_p ...
            + budgets_kinetic.Dissipation;
        budgets_entropy = COMPUTE_BUDGETS_ENTROPY(V, s, chemistry, output_flow);
        temp_entropic = budgets_entropy.A_adv + budgets_entropy.P_mom ...
            + budgets_entropy.P_mass + budgets_entropy.P_T ...
            + budgets_entropy.Transport + budgets_entropy.Dissipation ...
            + budgets_entropy.Source;
        temp = temp_entropic + temp_kinetic;
        scale_kinetic  = max(abs(temp), [], "all");
        scale_entropic = max(abs(temp), [], "all");
    end

    %% Time-advance perturbation on GPU (optional)
    if T_plot > 0
        order = 5;

        s = CFL_TIMESTEP(s);
        t_temp = s.time_integration.dt / s.time_integration.CFL;
        n_t = round(T_plot / t_temp);
        dt  = T_plot / n_t;
        fprintf("\n")
        fprintf("Number timesteps = %f ", n_t)
        fprintf("\n")

        %% Detect GPU availability
        use_gpu = (gpuDeviceCount > 0);
        n_vec = size(A, 1);

        if use_gpu
            gpu = gpuDevice;
            disp(gpu.Name + " GPU selected.")
            gA       = gpuArray(A * dt);
            v_buf    = gpuArray(zeros(n_vec, 1));
            vout_buf = gpuArray(zeros(n_vec, 1));
            term_buf = gpuArray(zeros(n_vec, 1));
        else
            disp("No CUDA device found. Running on CPU.")
            gA       = A * dt;
            v_buf    = zeros(n_vec, 1);
            vout_buf = zeros(n_vec, 1);
            term_buf = zeros(n_vec, 1);
        end

        tic
        [V] = LINEAR_INTEGRATION_GPU(V, n_t, order, gA, use_gpu, v_buf, vout_buf, term_buf);
        toc

        %% GPU synchronization and memory cleanup
        if use_gpu
            wait(gpu);
            clear gA v_buf vout_buf term_buf;
            wait(gpu);
        else
            clear gA v_buf vout_buf term_buf;
        end

        pert_rho   = real(reshape(V((0)*N_v+1:1*N_v, 1), s.mesh.Nchi, s.mesh.Neta));
        pert_rho_u = real(reshape(V((1)*N_v+1:2*N_v, 1), s.mesh.Nchi, s.mesh.Neta));
        pert_rho_v = real(reshape(V((2)*N_v+1:3*N_v, 1), s.mesh.Nchi, s.mesh.Neta));
        pert_rho_E = real(reshape(V((3)*N_v+1:4*N_v, 1), s.mesh.Nchi, s.mesh.Neta));

        solution_remesh = false;
        if freestream_disturbances
            [pert_rho, pert_rho_u, pert_rho_v, pert_rho_E] = ADD_FREESTREAM_DISTURBANCES( ...
                pert_rho, pert_rho_u, pert_rho_v, pert_rho_E, ...
                pert_infty_all, s, w_infty, T_plot, solution_remesh);
        end
    end

    %% Derived perturbation quantities at time T
    p_0   = s.var.p(2:end-1, 2:end-1);
    rho_0 = s.var.rho(2:end-1, 2:end-1);
    u_0   = s.var.rho_u(2:end-1, 2:end-1) ./ s.var.rho(2:end-1, 2:end-1);
    v_0   = s.var.rho_v(2:end-1, 2:end-1) ./ s.var.rho(2:end-1, 2:end-1);
    pert_u = (pert_rho_u - u_0 .* pert_rho) ./ rho_0;
    pert_v = (pert_rho_v - v_0 .* pert_rho) ./ rho_0;

    cell_x_wall_normal = (s.mesh.x_wall_normal(1:end-1, 1) + s.mesh.x_wall_normal(2:end, 1)) / 2;
    cell_y_wall_normal = (s.mesh.y_wall_normal(1:end-1, 1) + s.mesh.y_wall_normal(2:end, 1)) / 2;
    cell_x_wall_tan    =  cell_y_wall_normal;
    cell_y_wall_tan    = -cell_x_wall_normal;

    pert_u_tan  = pert_u .* cell_x_wall_tan + pert_v .* cell_y_wall_tan;
    pert_u_norm = pert_u .* cell_x_wall_normal + pert_v .* cell_y_wall_normal;
    pert_u_mag  = (pert_u .* u_0 + pert_v .* v_0) ./ sqrt(u_0.^2 + v_0.^2);

    t1 = -u_0 .* pert_rho_u;
    t2 = -v_0 .* pert_rho_v;
    t3 = (u_0.^2 + v_0.^2) .* pert_rho / 2;
    pert_p = (s.var.gamma_star(2:end-1, 2:end-1) - 1) .* (pert_rho_E + t1 + t2 + t3);

    pert_u_Ext = EXTEND_TO_GHOST_POINTS(pert_u, s);
    pert_v_Ext = EXTEND_TO_GHOST_POINTS(pert_v, s);
    [dpert_u_dx, dpert_u_dy] = DERIVATIVE_EXT(pert_u_Ext, s);
    [dpert_v_dx, dpert_v_dy] = DERIVATIVE_EXT(pert_v_Ext, s);
    pert_div  = dpert_u_dx + dpert_v_dy;
    pert_vort = dpert_v_dx - dpert_u_dy;

    pert_entropy = pert_p ./ p_0 ./ (s.var.gamma_star(2:end-1, 2:end-1) - 1) ...
        - pert_rho ./ rho_0 .* s.var.gamma_star(2:end-1, 2:end-1) ...
        ./ (s.var.gamma_star(2:end-1, 2:end-1) - 1);

    %% Shock spline (possibly perturbed)
    N_plot_shock_poly = 1000;
    if s.shock.enabled == true
        solution2 = s;
        if s.stability_analysis.perturb_shock
            dr_shock = real(V((4)*N_v+1:4*N_v+s.mesh.Nchi, 1));
            ang = s.shock.beta;
            solution2.shock.points_x = s.shock.points_x ...
                - dr_shock .* sin(ang) / s.stability_analysis.perturbation_magnitude / 2000000;
            solution2.shock.points_y = s.shock.points_y ...
                + dr_shock .* cos(ang) / s.stability_analysis.perturbation_magnitude / 2000000;
            figure()
            plot(solution2.shock.points_x, solution2.shock.points_y)
            hold on
            solution2 = LEAST_SQUARES_SHOCK_POINTS(solution2);
            plot(solution2.shock.points_x, solution2.shock.points_y)
            legend('original', 'filtered')
            hold off
        elseif T_plot > 0
            valid_ix      = (1:s.mesh.Nchi)';
            sc_idy        = s.shock.cell_indices(valid_ix, 1);
            lin_idx_shock = sub2ind(size(s.mesh.x), valid_ix, sc_idy - 1);
            p_shock       = real(pert_p(lin_idx_shock));
            max_p_shock   = max(p_shock, [], "all");
            min_p_shock   = min(p_shock, [], "all");
            p_shock       = (p_shock - min_p_shock) / (max_p_shock - min_p_shock) - 1/2;
            displace      = -sin(pi * p_shock);
            solution2.shock.points_x = s.shock.points_x + displace / 500;
            solution2.shock.points_y = s.shock.points_y + displace / 500;
            solution2 = LEAST_SQUARES_SHOCK_POINTS(solution2);
        end

        chi_p = linspace(solution2.shock.points_chi(1, 1), ...
            solution2.shock.points_chi(end, 1), N_plot_shock_poly);
        x_p = ppval(solution2.shock.spline_func_x, chi_p);
        y_p = ppval(solution2.shock.spline_func_y, chi_p);
    end

    %% Configure plot parameters
    plot_params = struct();
    if T_plot == 0 && freestream_disturbances
        plot_params.fig_width = 10;
        plot_params.freestream_disturbance = true;
    else
        plot_params.fig_width = 5;
        plot_params.freestream_disturbance = false;
    end
    plot_params.fig_height      = 8;
    plot_params.font_size       = 20;
    plot_params.label_size      = 20;
    plot_params.cb_label_size   = 18;
    plot_params.line_width      = 1.5;
    plot_params.colormap_name   = 'turbo';
    plot_params.use_common_clim = false;

    %% Define energy-budget plot lists
    if T_plot > 0 || ~freestream_disturbances
        output_flow = true;
        budgets_kinetic = COMPUTE_BUDGETS_KINETIC(V, s, output_flow);

        kinetic_plots = {
            {budgets_kinetic.A_adv / scale_kinetic,          '$\displaystyle \mathcal{A}^k$'}
            {budgets_kinetic.P_mom / scale_kinetic,          '$\displaystyle \mathcal{P}_u^k$'}
            {budgets_kinetic.P_mass / scale_kinetic,         '$\displaystyle \mathcal{P}_\rho^k$'}
            {budgets_kinetic.Dilat_P / scale_kinetic,        '$\displaystyle \Pi_d^k$'}
            {budgets_kinetic.Transport_p / scale_kinetic,    '$\displaystyle \mathcal{T}^k_p$'}
            {budgets_kinetic.Transport_tau / scale_kinetic,  '$\displaystyle \mathcal{T}^k_\tau$'}
            {budgets_kinetic.Dissipation / scale_kinetic,    '$\displaystyle \mathcal{D}^k$'}
        };

        output_flow = true;
        budgets_entropy = COMPUTE_BUDGETS_ENTROPY(V, s, chemistry, output_flow);

        entropy_plots = {
            {budgets_entropy.A_adv / scale_entropic,         '$\displaystyle \mathcal{A}^s$'}
            {budgets_entropy.P_mom / scale_entropic,         '$\displaystyle \mathcal{P}_u^s$'}
            {budgets_entropy.P_mass / scale_entropic,        '$\displaystyle \mathcal{P}_\rho^s$'}
            {budgets_entropy.P_T / scale_entropic,           '$\displaystyle \mathcal{P}_T^s$'}
            {budgets_entropy.Transport / scale_entropic,     '$\displaystyle \mathcal{T}^s$'}
            {budgets_entropy.Dissipation / scale_entropic,   '$\displaystyle \mathcal{D}^s$'}
            {budgets_entropy.Source / scale_entropic,        '$\displaystyle \mathcal{S}^s$'}
        };
    end

    %% Build coordinate arrays and create all plots
    fprintf('Creating perturbation plots...\n');

    y_coord = s.mesh.y;
    x_coord = s.mesh.x;
    y_shock = y_p;
    x_shock = x_p;

    %% Render kinetic and entropic budget plots
    if T_plot > 0 || ~freestream_disturbances
        fprintf('Creating kinetic budget plots...\n');
        for i = 1:length(kinetic_plots)
            data      = kinetic_plots{i}{1};
            cb_label  = kinetic_plots{i}{2};
            data_plot = real(data);
            CREATE_PUBLICATION_PLOT(y_coord, x_coord, data_plot, cb_label, ...
                s.shock.enabled, y_shock, x_shock, plot_params, s);
        end

        fprintf('Creating entropy budget plots...\n');
        for i = 1:length(entropy_plots)
            data      = entropy_plots{i}{1};
            cb_label  = entropy_plots{i}{2};
            data_plot = real(data);
            CREATE_PUBLICATION_PLOT(y_coord, x_coord, data_plot, cb_label, ...
                s.shock.enabled, y_shock, x_shock, plot_params, s);
        end
    end

    fprintf('All plots created successfully!\n');

end


%% ========================================================================
%  LOCAL FUNCTION: LINEAR_INTEGRATION_GPU
%  ========================================================================
function [voutCPU] = LINEAR_INTEGRATION_GPU(x, n_t, order, Adt, use_gpu, v_buf, vout_buf, term_buf)
% LINEAR_INTEGRATION_GPU  Approximate exp(A*t)*x via truncated Taylor series.
%   Works on both GPU and CPU using pre-allocated buffers.
%
%   Inputs:
%       x        - (double array) initial state vector (CPU).
%       n_t      - (scalar) number of time steps.
%       order    - (scalar) Taylor-series truncation order.
%       Adt      - (gpuArray or sparse) system matrix A*dt.
%       use_gpu  - (logical) true for GPU path, false for CPU path.
%       v_buf    - Pre-allocated buffer for working vector.
%       vout_buf - Pre-allocated buffer for output accumulation.
%       term_buf - Pre-allocated buffer for Taylor term.
%
%   Outputs:
%       voutCPU - (double array) final state vector (CPU).

    if use_gpu
        v_buf(:) = x;
    else
        v_buf(:) = x;
    end
    vout_buf(:) = v_buf;

    for i = 1:n_t
        factorial_val = 1;
        for j = 1:order
            term_buf(:) = Adt * v_buf;
            v_buf(:) = term_buf;
            factorial_val = j * factorial_val;
            vout_buf(:) = vout_buf + v_buf / factorial_val;
        end
        v_buf(:) = vout_buf;
    end

    if use_gpu
        voutCPU = gather(vout_buf);
    else
        voutCPU = vout_buf;
    end
end


%% ========================================================================
%  LOCAL FUNCTION: CREATE_PUBLICATION_PLOT
%  ========================================================================
function CREATE_PUBLICATION_PLOT(x, y, data, cb_label, has_shock, x_shock, y_shock, params, s)
% CREATE_PUBLICATION_PLOT  Generate a single publication-quality pcolor
%   figure with a symmetric colorbar centred at zero.
%
%   Inputs:
%       x, y      - (2-D arrays) coordinate grids for pcolor.
%       data      - (2-D array) field data to plot.
%       cb_label  - (char) LaTeX string for the colorbar label.
%       has_shock - (logical) overlay shock line when true.
%       y_shock, x_shock - (vectors) shock-line coordinates.
%       params    - (struct) figure, font, and colormap parameters.
%       s  - (struct) s struct (boundary_type, L, etc.).
%
%   Outputs:
%       (none) - creates a new figure window.

    %% Figure sizing based on data aspect ratio
    height_range = max(y(:)) - min(y(:));
    width_range  = max(x(:)) - min(x(:));

    data_aspect = height_range / width_range / 1.1;

    plot_width  = params.fig_width;
    plot_height = plot_width * data_aspect;
    fig_height  = plot_height;

    %% Create figure window
    figure('Units', 'inches', ...
        'Position', [1, 1, plot_width, fig_height], ...
        'Color', 'white', ...
        'PaperPositionMode', 'auto', ...
        'PaperSize', [plot_width, fig_height]);

    %% Draw pcolor surface
    ax = pcolor(x, y, data);
    ax.FaceColor = 'interp';
    set(ax, 'EdgeColor', 'none');

    hold on
    if has_shock
        plot(y_shock, x_shock, 'k-', 'LineWidth', params.line_width)
    end

    colormap(params.colormap_name)

    %% Symmetric colour axis truncated to 1 significant figure
    max_abs_val = max(abs(data(:)));
    truncate_to_1sig = @(x) floor(x ./ 10.^floor(log10(abs(x)))) .* 10.^floor(log10(abs(x)));
    max_abs_val = truncate_to_1sig(max_abs_val);

    scaling_range = 1/2;
    caxis([-max_abs_val * scaling_range, max_abs_val * scaling_range]);

    %% Configure colorbar
    cb = colorbar;
    cb.Ticks = linspace(-max_abs_val * scaling_range, max_abs_val * scaling_range, 3);
    cb.Label.Interpreter       = 'latex';
    cb.Label.String            = cb_label;
    cb.Label.FontSize          = params.cb_label_size;
    cb.Label.Rotation          = 0;
    cb.Label.Position          = [6, 0];
    cb.Label.VerticalAlignment = 'top';
    cb.TickLabelInterpreter    = 'latex';
    cb.LineWidth               = 0.75;

    %% Axis ticks, labels, and limits
    xlabel('$x$', 'FontSize', params.label_size, 'Interpreter', 'latex')
    ylabel('$y$', 'FontSize', params.label_size, 'Interpreter', 'latex')

    set(gca, 'FontSize', params.font_size, ...
        'TickLabelInterpreter', 'latex', ...
        'LineWidth', 0.75, ...
        'Layer', 'top', ...
        'Box', 'on')

    daspect([1 1 1])

    hold off
end


%% ========================================================================
%  LOCAL FUNCTION: EXPORT_ALL_FIGURES
%  ========================================================================
function EXPORT_ALL_FIGURES(format, dpi)
% EXPORT_ALL_FIGURES  Export every open figure to a raster or vector file.
%
%   EXPORT_ALL_FIGURES(format)
%   EXPORT_ALL_FIGURES(format, dpi)
%
%   Inputs:
%       format - (char) 'png', 'pdf', 'eps', or 'tiff'.
%       dpi    - (scalar, optional) resolution for raster formats
%                (default 600).
%
%   Outputs:
%       (none) - writes figure_001, figure_002, ... files to disk.

    if nargin < 2
        dpi = 600;
    end

    fig_handles = findall(0, 'Type', 'figure');

    for i = 1:length(fig_handles)
        filename = sprintf('figure_%03d', i);
        figure(fig_handles(i));

        switch lower(format)
            case 'png'
                print(filename, '-dpng', sprintf('-r%d', dpi));
            case 'pdf'
                exportgraphics(gcf, [filename '.pdf'], 'ContentType', 'vector');
            case 'eps'
                print(filename, '-depsc', sprintf('-r%d', dpi));
            case 'tiff'
                print(filename, '-dtiff', sprintf('-r%d', dpi));
            otherwise
                warning('Unknown format: %s', format);
        end

        fprintf('Exported: %s.%s\n', filename, format);
    end
end


%% ========================================================================
%  LOCAL FUNCTION: SCALINGS
%  ========================================================================
function [scale_u, scale_vort, scale_div, scale_rho, scale_p, scale_s] = SCALINGS(s, pert_rho, pert_rho_u, pert_rho_v, pert_rho_E)
% SCALINGS  Compute maximum-absolute-value scaling factors for each
%   perturbation quantity, used to normalise contour plots.
%
%   Inputs:
%       s               - (struct) base-flow s.
%       pert_rho, pert_rho_u,
%       pert_rho_v, pert_rho_E - (2-D arrays) perturbation conservative
%                                 variables.
%
%   Outputs:
%       scale_u    - max |u'_mag|
%       scale_vort - max |vorticity'|
%       scale_div  - max |divergence'|
%       scale_rho  - max |rho'|
%       scale_p    - max |p'|
%       scale_s    - max |entropy'|

    p_0   = s.var.p(2:end-1, 2:end-1);
    rho_0 = s.var.rho(2:end-1, 2:end-1);
    u_0   = s.var.rho_u(2:end-1, 2:end-1) ./ s.var.rho(2:end-1, 2:end-1);
    v_0   = s.var.rho_v(2:end-1, 2:end-1) ./ s.var.rho(2:end-1, 2:end-1);

    pert_u     = (pert_rho_u - u_0 .* pert_rho) ./ rho_0;
    pert_v     = (pert_rho_v - v_0 .* pert_rho) ./ rho_0;
    pert_u_mag = (pert_u .* u_0 + pert_v .* v_0) ./ sqrt(u_0.^2 + v_0.^2);

    t1 = -u_0 .* pert_rho_u;
    t2 = -v_0 .* pert_rho_v;
    t3 = (u_0.^2 + v_0.^2) .* pert_rho / 2;
    pert_p = (s.var.gamma_star(2:end-1, 2:end-1) - 1) .* (pert_rho_E + t1 + t2 + t3);

    pert_u_Ext = EXTEND_TO_GHOST_POINTS(pert_u, s);
    pert_v_Ext = EXTEND_TO_GHOST_POINTS(pert_v, s);
    [dpert_u_dx, dpert_u_dy] = DERIVATIVE_EXT(pert_u_Ext, s);
    [dpert_v_dx, dpert_v_dy] = DERIVATIVE_EXT(pert_v_Ext, s);

    pert_div  = dpert_u_dx + dpert_v_dy;
    pert_vort = dpert_v_dx - dpert_u_dy;

    pert_entropy = pert_p ./ p_0 ./ (s.var.gamma_star(2:end-1, 2:end-1) - 1) ...
        - pert_rho ./ rho_0 .* s.var.gamma_star(2:end-1, 2:end-1) ...
        ./ (s.var.gamma_star(2:end-1, 2:end-1) - 1);

    scale_u    = max(abs(pert_u_mag),    [], "all");
    scale_vort = max(abs(pert_vort),     [], "all");
    scale_div  = max(abs(pert_div),      [], "all");
    scale_rho  = max(abs(pert_rho),      [], "all");
    scale_p    = max(abs(pert_p),        [], "all");
    scale_s    = max(abs(pert_entropy),  [], "all");
end


%% ========================================================================
%  LOCAL FUNCTION: ADD_FREESTREAM_DISTURBANCES
%  ========================================================================
function [pert_rho, pert_rho_u, pert_rho_v, pert_rho_E] = ADD_FREESTREAM_DISTURBANCES( ...
        pert_rho, pert_rho_u, pert_rho_v, pert_rho_E, ...
        pert_infty_all, s, w_infty, T_plot, solution_remesh)
% ADD_FREESTREAM_DISTURBANCES  Superpose freestream Fourier modes onto the
%   interior perturbation fields beyond the shock.
%
%   Inputs:
%       pert_rho, pert_rho_u,
%       pert_rho_v, pert_rho_E - (Nx x Ny) interior perturbation arrays.
%       pert_infty_all         - (vector) packed freestream amplitudes.
%       s               - (struct) grid and shock data.
%       w_infty                - (vector) freestream frequencies.
%       T_plot                 - (scalar) evaluation time.
%       solution_remesh        - (logical) if true, zero interior first.
%
%   Outputs:
%       pert_rho, pert_rho_u,
%       pert_rho_v, pert_rho_E - Updated perturbation arrays with
%                                 freestream contributions added.

    %% Unpack freestream disturbance coefficients
    N_f        = length(w_infty);
    Nx         = s.mesh.Nchi;
    Ny         = s.mesh.Neta;
    size_infty = s.mesh.Nchi * N_f;

    rho_infty_all   = zeros(size_infty, 1);
    rho_u_infty_all = zeros(size_infty, 1);
    rho_v_infty_all = zeros(size_infty, 1);
    rho_E_infty_all = zeros(size_infty, 1);

    for i = 0:N_f-1
        rho_infty_all(  i*Nx+1:i*Nx+Nx, 1) = pert_infty_all(4*i*Nx+1      : 4*i*Nx+Nx,   1);
        rho_u_infty_all(i*Nx+1:i*Nx+Nx, 1) = pert_infty_all(4*i*Nx+Nx+1   : 4*i*Nx+2*Nx, 1);
        rho_v_infty_all(i*Nx+1:i*Nx+Nx, 1) = pert_infty_all(4*i*Nx+2*Nx+1 : 4*i*Nx+3*Nx, 1);
        rho_E_infty_all(i*Nx+1:i*Nx+Nx, 1) = pert_infty_all(4*i*Nx+3*Nx+1 : 4*i*Nx+4*Nx, 1);
    end

    %% Zero freestream region and superpose Fourier modes
    indices   = (1:Nx)';
    shock_idx = s.shock.cell_indices(indices, 1);
    shock_x   = s.shock.points_x(indices, 1);

    if solution_remesh
        pert_rho   = zeros(Nx, Ny);
        pert_rho_u = zeros(Nx, Ny);
        pert_rho_v = zeros(Nx, Ny);
        pert_rho_E = zeros(Nx, Ny);
    else
        pert_rho(  indices, shock_idx:end) = 0;
        pert_rho_u(indices, shock_idx:end) = 0;
        pert_rho_v(indices, shock_idx:end) = 0;
        pert_rho_E(indices, shock_idx:end) = 0;
    end

    factor = 1;
    for i = 1:N_f
        y    = (s.mesh.x(indices, shock_idx:end) - shock_x);
        time = T_plot;

        pert_rho(  indices, shock_idx:end) = pert_rho(  indices, shock_idx:end) ...
            + rho_infty_all(  indices + (i-1)*Nx, 1) .* exp(1i*w_infty(i)*(y - time)) / factor;
        pert_rho_u(indices, shock_idx:end) = pert_rho_u(indices, shock_idx:end) ...
            + rho_u_infty_all(indices + (i-1)*Nx, 1) .* exp(1i*w_infty(i)*(y - time)) / factor;
        pert_rho_v(indices, shock_idx:end) = pert_rho_v(indices, shock_idx:end) ...
            + rho_v_infty_all(indices + (i-1)*Nx, 1) .* exp(1i*w_infty(i)*(y - time)) / factor;
        pert_rho_E(indices, shock_idx:end) = pert_rho_E(indices, shock_idx:end) ...
            + rho_E_infty_all(indices + (i-1)*Nx, 1) .* exp(1i*w_infty(i)*(y - time)) / factor;
    end
end

