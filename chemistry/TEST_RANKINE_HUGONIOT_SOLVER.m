function TEST_RANKINE_HUGONIOT_SOLVER(chemistry, s)
% TEST_RANKINE_HUGONIOT_SOLVER  Validate the Rankine-Hugoniot shock solver
%   against reference Python / SD Toolbox data.
%
%   TEST_RANKINE_HUGONIOT_SOLVER(chemistry, s) runs three test
%   cases:
%     1. Single moderate-shock solve (sanity check).
%     2. Mach-number sweep compared with SD Toolbox reference data.
%     3. Performance benchmark (100 simultaneous solves).
%
%   Inputs:
%       chemistry - Fitted chemistry structure from FIT_CHEMISTRY
%                   (must contain eval_gamma_star, eval_T, eval_e)
%       s  - Solution configuration structure with fields:
%                     .chemistry.is_chemistry_enabled  (logical)
%                     .chemistry.chemistry_type        (string, e.g. "Chemical-RTV")
%                     .chemistry.chemistry_composition (string, e.g. "Earth")
%
%   Outputs:
%       (none -- results are printed to the console and plotted)
%
%   Notes:
%       Reference shock-jump conditions are read from Python-generated
%       files via READ_SHOCK_JUMP_CONDITIONS.  If those files are absent,
%       default upstream conditions are used instead.
%
% Part of: Hypersonics Stability MATLAB Solver - Chemistry Module

    %% Guard: chemistry must be enabled
    if ~s.chemistry.is_chemistry_enabled
        fprintf('No chemistry enabled...\n');
        s.chemistry.is_chemistry_enabled  = false;
        s.chemistry.chemistry_type   = "None";
        s.chemistry.chemistry_composition = "None";
        return
    end

    if ~isfield(chemistry, 'eval_gamma_star')
        error('Chemistry structure must be fitted first. Run FIT_CHEMISTRY(chemistry)');
    end

    fprintf('Testing vectorized non-linear solver...\n\n');

    %% Read reference shock-jump data
    fprintf('=== Reading Python Shock Jump Conditions ===\n');
    planet = s.chemistry.chemistry_composition;
    chemistry_type = s.chemistry.chemistry_type;

    shock_jump_test = READ_SHOCK_JUMP_CONDITIONS(planet, chemistry_type);

    %% Test 1: single moderate shock
    fprintf('\n=== Test Case 1: Single Moderate Shock ===\n');
    w_1 = 1000;    % m/s
    rho_1 = 0.1;   % kg/m^3
    e_1 = 5e5;     % J/kg

    [rho_2, e_2, w_2] = SOLVE_RANKINE_HUGONIOT_CHEMISTRY(chemistry, w_1, rho_1, e_1);

    fprintf('Single case results: rho_2 = %.3e kg/m^3, e_2 = %.3e J/kg, w_2 = %.3e m/s\n\n', ...
        rho_2, e_2, w_2);

    %% Test 2: Mach-number sweep with Python data comparison
    fprintf('=== Test Case 2: Mach Number Sweep with Python Data Comparison ===\n');

    if ~isempty(fieldnames(shock_jump_test.upstream))
        rho_1_plot = shock_jump_test.upstream.rho1;
        e_1_plot   = shock_jump_test.upstream.e1;
        a_s        = shock_jump_test.upstream.a_s1;
        T_1_plot   = shock_jump_test.upstream.T1;

        fprintf('Using upstream conditions from Python files:\n');
        fprintf('  rho_1 = %.3e kg/m^3\n', rho_1_plot);
        fprintf('  e_1   = %.3e J/kg\n', e_1_plot);
        fprintf('  T_1   = %.1f K\n', T_1_plot);
        fprintf('  a_s   = %.1f m/s\n', a_s);
    else
        rho_1_plot = 0.1;   % kg/m^3
        e_1_plot   = 2e5;   % J/kg
        a_s        = 340;   % m/s
        T_1_plot   = 300;   % K

        fprintf('Using default upstream conditions (Python data not available):\n');
        fprintf('  rho_1 = %.3e kg/m^3\n', rho_1_plot);
        fprintf('  e_1   = %.3e J/kg\n', e_1_plot);
        fprintf('  T_1   = %.1f K\n', T_1_plot);
        fprintf('  a_s   = %.1f m/s\n', a_s);
    end

    PLOT_TEST_RANKINE_HUGONIOT(chemistry, rho_1_plot, e_1_plot, a_s, T_1_plot, shock_jump_test);

    %% Test 3: performance benchmark
    fprintf('\n=== Test Case 3: Performance Test ===\n');
    N_perf = 100;
    w_1_perf = linspace(500, 2000, N_perf);

    fprintf('Performance test with %d points...\n', N_perf);
    tic;
    [rho_2_perf, e_2_perf, w_2_perf] = SOLVE_RANKINE_HUGONIOT_CHEMISTRY( ...
        chemistry, w_1_perf, rho_1_plot, e_1_plot);
    elapsed_time = toc;

    fprintf('Performance results:\n');
    fprintf('  Total time: %.3f seconds\n', elapsed_time);
    fprintf('  Time per solve: %.1f ms\n', 1000 * elapsed_time / N_perf);

    fprintf('\n=== All Tests Completed Successfully ===\n');
end

%% ========================================================================
%  Plotting helper
%  ========================================================================
function PLOT_TEST_RANKINE_HUGONIOT(chemistry, rho_1, e_1, a_s, T_1, shock_jump_test)
% PLOT_TEST_RANKINE_HUGONIOT  Plot density and temperature ratios vs Mach
%   number with comparison to Python / SD Toolbox reference data.
%
%   PLOT_TEST_RANKINE_HUGONIOT(chemistry, rho_1, e_1, a_s, T_1, shock_jump_test)
%
%   Inputs:
%       chemistry       - Fitted chemistry structure from FIT_CHEMISTRY
%       rho_1           - Upstream density [kg/m^3] (scalar)
%       e_1             - Upstream specific energy [J/kg] (scalar)
%       a_s             - Upstream sound speed [m/s] (scalar)
%       T_1             - Upstream temperature [K] (scalar)
%       shock_jump_test - Structure with Python shock-jump data

    %% Input validation
    if ~isfield(chemistry, 'eval_gamma_star')
        error('Chemistry structure must be fitted first. Run FIT_CHEMISTRY(chemistry)');
    end

    if ~isscalar(rho_1) || ~isscalar(e_1) || ~isscalar(a_s) || ~isscalar(T_1)
        error('All inputs except shock_jump_test must be scalars for this plotting function');
    end

    if ~isfinite(T_1) || T_1 <= 0
        error('T_1 must be a positive finite number, got T_1 = %.2f', T_1);
    end

    if ~isfinite(rho_1) || rho_1 <= 0
        error('rho_1 must be a positive finite number, got rho_1 = %.2e', rho_1);
    end

    % Recompute e_1 from the interpolation for consistency
    e_1 = chemistry.eval_e(T_1, rho_1);

    %% Mach-number sweep
    M1_range = linspace(1.5, 30, 100);
    N_points = length(M1_range);
    w_1_range = M1_range * a_s;

    fprintf('Computing MATLAB Rankine-Hugoniot relations for M1 = %.1f to %.1f...\n', ...
        min(M1_range), max(M1_range));
    fprintf('Using %d Mach number points\n', N_points);
    fprintf('State 1 conditions:\n');
    fprintf('  rho_1 = %.3e kg/m^3\n', rho_1);
    fprintf('  e_1   = %.3e J/kg\n', e_1);
    fprintf('  T_1   = %.1f K\n', T_1);
    fprintf('  a_s   = %.1f m/s\n', a_s);

    fprintf('Solving vectorized Rankine-Hugoniot equations...\n');
    [rho_2_range, e_2_range, w_2_range] = SOLVE_RANKINE_HUGONIOT_CHEMISTRY( ...
        chemistry, w_1_range, rho_1, e_1);

    if any(isnan(rho_2_range)) || any(isinf(rho_2_range)) || ...
       any(isnan(e_2_range))   || any(isinf(e_2_range))
        warning('Solver produced NaN or Inf values. This may indicate convergence issues.');
    end

    %% Compute ratios
    density_ratio = rho_2_range ./ rho_1;

    T_2_range = chemistry.eval_T(rho_2_range, e_2_range);
    if all(isnan(T_2_range))
        error('All computed temperatures are NaN. Check chemistry interpolation data range.');
    end

    temperature_ratio = T_2_range ./ T_1;

    % Filter invalid values for plotting
    valid_idx = isfinite(density_ratio) & isfinite(temperature_ratio) & ...
                density_ratio > 0 & temperature_ratio > 0;

    if sum(valid_idx) == 0
        error('No valid solutions found. Check chemistry data and initial conditions.');
    end

    M1_plot = M1_range(valid_idx);
    density_ratio_plot = density_ratio(valid_idx);
    temperature_ratio_plot = temperature_ratio(valid_idx);

    fprintf('Using %d valid points out of %d total points for plotting\n', ...
        sum(valid_idx), length(valid_idx));

    %% Create figure
    figure('Position', [100, 100, 1000, 400]);

    % --- Subplot 1: density ratio ---
    subplot(1, 2, 1);
    hold on;
    plot(M1_plot, density_ratio_plot, 'b-', 'LineWidth', 2.5, ...
        'DisplayName', 'Present Solver');

    if isfield(shock_jump_test, 'chemistry_type') && isfield(shock_jump_test, 'data')
        chemistry_type = shock_jump_test.chemistry_type;
        field_name = strrep(chemistry_type, '-', '_');

        if isfield(shock_jump_test.data, field_name)
            data = shock_jump_test.data.(field_name);
            fprintf('Adding Python data to density plot for %s:\n', char(chemistry_type));

            python_density_ratio = data.rho2 / rho_1;
            scatter(data.M1(1:2:end), python_density_ratio(1:2:end), 30, 'r', 'o', ...
                'filled', 'DisplayName', sprintf('SD Toolbox %s', chemistry_type), ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.4);

            fprintf('  %s: %d points\n', char(chemistry_type), length(data.M1));
        end
    end

    grid on;
    xlabel('$M_1$', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('$\frac{\rho_2}{\rho_1}$', 'FontSize', 20, 'Interpreter', 'latex');
    xlim([1.5, 35]);

    density_max = max(density_ratio_plot);
    if isfinite(density_max) && density_max > 0
        ylim([0, density_max * 1.1]);
    else
        ylim([0, 10]);
    end
    legend('Location', 'southeast', 'FontSize', 10);

    % --- Subplot 2: temperature ratio ---
    subplot(1, 2, 2);
    hold on;
    plot(M1_plot, temperature_ratio_plot, 'b-', 'LineWidth', 2.5, ...
        'DisplayName', 'Present Solver');

    if isfield(shock_jump_test, 'chemistry_type') && isfield(shock_jump_test, 'data')
        chemistry_type = shock_jump_test.chemistry_type;
        field_name = strrep(chemistry_type, '-', '_');

        if isfield(shock_jump_test.data, field_name)
            data = shock_jump_test.data.(field_name);
            fprintf('Adding Python data to temperature plot for %s:\n', char(chemistry_type));

            python_temperature_ratio = data.T2 / T_1;
            scatter(data.M1(1:2:end), python_temperature_ratio(1:2:end), 30, 'r', 'o', ...
                'filled', 'DisplayName', sprintf('SD Toolbox %s', chemistry_type), ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.4);

            fprintf('  %s: %d points\n', char(chemistry_type), length(data.M1));
        end
    end

    grid on;
    xlabel('$M_1$', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('$\frac{T_2}{T_1}$', 'FontSize', 20, 'Interpreter', 'latex');
    xlim([1.5, 35]);

    temp_max = max(temperature_ratio_plot);
    if isfinite(temp_max) && temp_max > 0
        ylim([0, temp_max * 1.1]);
    else
        ylim([0, 20]);
    end
    legend('Location', 'southeast', 'FontSize', 10);

    %% Print summary statistics
    fprintf('\nMATLAB Results Summary:\n');
    fprintf('Mach number range: %.1f - %.1f\n', min(M1_plot), max(M1_plot));
    fprintf('Density ratio range: %.2f - %.2f\n', min(density_ratio_plot), max(density_ratio_plot));
    fprintf('Temperature ratio range: %.2f - %.2f\n', min(temperature_ratio_plot), max(temperature_ratio_plot));
    fprintf('Temperature range at state 2: %.1f - %.1f K\n', min(T_2_range), max(T_2_range));

    % Compare with Python reference data
    if isfield(shock_jump_test, 'chemistry_type') && isfield(shock_jump_test, 'data')
        chemistry_type = shock_jump_test.chemistry_type;
        field_name = strrep(chemistry_type, '-', '_');

        if isfield(shock_jump_test.data, field_name)
            data = shock_jump_test.data.(field_name);
            fprintf('\nPython Data Summary for %s:\n', char(chemistry_type));
            fprintf('  Mach range: %.2f - %.2f\n', min(data.M1), max(data.M1));
            fprintf('  Density ratio range: %.2f - %.2f\n', min(data.rho2) / rho_1, max(data.rho2) / rho_1);
            fprintf('  Temperature ratio range: %.2f - %.2f\n', min(data.T2) / T_1, max(data.T2) / T_1);
            fprintf('  Temperature range: %.1f - %.1f K\n', min(data.T2), max(data.T2));
        end
    end
end
