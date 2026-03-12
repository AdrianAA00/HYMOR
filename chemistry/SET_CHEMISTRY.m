function chemistry = SET_CHEMISTRY(s)
% SET_CHEMISTRY - Initialize and configure chemistry module for the solver
%
% Master setup function for the chemistry module. Loads gas property data
% from equilibrium model files, creates interpolation fits for all
% thermodynamic and transport properties, validates the fits, generates
% diagnostic plots, and optionally loads non-equilibrium relaxation models.
%
% Syntax:
%   chemistry = SET_CHEMISTRY(s)
%
% Inputs:
%   s - Solution structure containing configuration fields:
%       .chemistry.is_chemistry_enabled  - (logical) Enable/disable chemistry
%       .chemistry.chemistry_type        - (string) Chemistry model:
%                                          'Frozen-RTV', 'Chemical-RTV',
%                                          'Chemical-RTVE', 'Frozen-RTVE'
%       .chemistry.chemistry_composition - (string) Planet: 'Earth' or 'Mars'
%       .chemistry.chemical_equilibrium  - (logical) Equilibrium vs non-equilibrium
%       .solver_dir                      - (string) Root directory of the solver
%
% Outputs:
%   chemistry - Fully initialized chemistry structure containing:
%       .planet, .gas_model   - Planet and gas model identifiers
%       .rho, .e, .T, etc.    - Raw tabulated gas property data
%       .fit_T, .fit_gamma_star, .fit_cv_star, .fit_mu, .fit_k
%                              - griddedInterpolant objects for each property
%       .fit_a, .fit_s        - Sound speed and entropy fits (if available)
%       .eval_T, .eval_gamma_star, .eval_cv_star, .eval_mu, .eval_k
%                              - Function handles for property evaluation
%       .eval_e               - Inverse interpolant E(T, rho)
%       .eval_log10_<species> - Species concentration evaluators (if available)
%       .fit_info             - Metadata about the fitting process
%       .fit_validation       - Validation metrics (RMSE, R^2, MAE)
%       .frozen               - Frozen chemistry data (for non-equilibrium)
%       .neq                  - Non-equilibrium relaxation model (if applicable)
%
% Notes:
%   When chemistry is disabled (s.chemistry.is_chemistry_enabled = false),
%   returns a minimal structure with planet = "NONE" and gas_model = "NONE".
%   For non-equilibrium simulations (s.chemistry.chemical_equilibrium = false),
%   both equilibrium and frozen gas property tables are loaded and fitted
%   independently. The frozen model is selected as 'Frozen-RTV' or
%   'Frozen-RTVE' based on the current chemistry_type.
%   The function also calls diagnostic routines when detailed_output is true:
%     TEST_FIT_PERFORMANCE, PLOT_FIT_EVALUATION, PLOT_3D_CHEMISTRY_SURFACES
%
% Part of: Hypersonics Stability MATLAB Solver - Chemistry Module



    %% Early return if chemistry is disabled
    if ~s.chemistry.is_chemistry_enabled
        fprintf('No chemistry enabled...\n');
        s.chemistry.is_chemistry_enabled  = false;
        s.chemistry.chemistry_type   = "None";
        s.chemistry.chemical_equilibrium = false;
        s.chemistry.chemistry_composition = "None";
        chemistry.planet = "NONE";
        chemistry.gas_model = "NONE";
        return
    end

    %% Load equilibrium gas property data
    fprintf('Loading chemistry data...\n');
    name = s.solver_dir + "chemistry/equilibrium_models/gas_properties_" + s.chemistry.chemistry_type + "_" + s.chemistry.chemistry_composition + ".txt";

    chemistry = struct();
    chemistry = READ_GAS_PROPERTIES(chemistry, name);
    chemistry.detailed_output = false;
    fprintf('Initializing chemistry module with detailed output: %d\n', chemistry.detailed_output);

    %% Load frozen gas property data for non-equilibrium shock jumps
    if ~s.chemistry.chemical_equilibrium
        if strcmp(s.chemistry.chemistry_type, "Chemical-RTV") || strcmp(s.chemistry.chemistry_type, "Frozen-RTV")
            chemistry_frozen = "Frozen-RTV";
        elseif strcmp(s.chemistry.chemistry_type, "Chemical-RTVE") || strcmp(s.chemistry.chemistry_type, "Frozen-RTVE")
            chemistry_frozen = "Frozen-RTVE";
        end

        chemistry.frozen = struct();
        name_frozen = s.solver_dir + "chemistry/equilibrium_models/gas_properties_" + chemistry_frozen + "_" + s.chemistry.chemistry_composition + ".txt";
        chemistry.frozen = READ_GAS_PROPERTIES(chemistry.frozen, name_frozen);
    end

    %% Fit interpolation models
    fprintf('\nFitting chemistry models...\n');
    chemistry = FIT_CHEMISTRY(chemistry);
    if ~s.chemistry.chemical_equilibrium
        chemistry.frozen = FIT_CHEMISTRY(chemistry.frozen);
    end

    %% Run diagnostics
    if chemistry.detailed_output
        TEST_FIT_PERFORMANCE(chemistry);
        PLOT_FIT_EVALUATION(chemistry);
    
        fprintf('\nCreating 3D surface plots...\n');
        PLOT_3D_CHEMISTRY_SURFACES(chemistry);
        TEST_RANKINE_HUGONIOT_SOLVER(chemistry, s);
    end

    %% Load non-equilibrium relaxation model
    if ~s.chemistry.chemical_equilibrium
        if strcmp(s.chemistry.chemistry_type, "Frozen-RTV")
            error("Frozen-RTV: frozen chemistry cannot have non equilibrium model")
        end

        fprintf('\n');
        fprintf('Loading non-equilibrium model...\n');
        name = s.solver_dir + "chemistry/non_equilibrium_models/non_equi_model_" + s.chemistry.chemistry_composition + ".txt";
        chemistry = READ_NON_EQUILIBRIUM_MODEL(chemistry, name);
        fprintf('\n');
    end
end


%% ========================================================================
%  FIT_CHEMISTRY - Create griddedInterpolant models for all properties
%  ========================================================================
function chemistry = FIT_CHEMISTRY(chemistry)
% FIT_CHEMISTRY - Creates fast griddedInterpolant models with log transforms
%
% Builds interpolation models for T, gamma_star, cv_star, mu, k, a, s and
% log10 species concentrations as functions of (rho, e). Uses log10
% transformations for density and energy, with all variables scaled to
% [-1, 1] for numerical accuracy. Also creates an inverse fit E(T, rho)
% for the Rankine-Hugoniot solver.
%
% Inputs:
%   chemistry - Structure containing raw gas properties data with fields:
%               .rho, .e, .T, .gamma_star, .cv_star, .mu, .k
%               and optionally .a (sound speed), .s (entropy),
%               .has_sound_speed, .has_entropy, .has_molar_data
%
% Outputs:
%   chemistry - Enhanced structure with griddedInterpolant objects and
%               function handles for property evaluation
%
% Notes:
%   Uses griddata with cubic method to interpolate scattered data onto a
%   100x100 regular grid, then creates griddedInterpolant objects with
%   cubic interpolation and linear extrapolation.

    fprintf('Starting chemistry data fitting with log(rho) + scaled griddedInterpolant...\n');
    tic;

    %% Validate input data
    if ~isfield(chemistry, 'rho') || ~isfield(chemistry, 'e')
        error('Chemistry structure must contain rho and E fields');
    end

    required_fields = {'T', 'gamma_star', 'cv_star', 'mu', 'k'};
    for i = 1:length(required_fields)
        if ~isfield(chemistry, required_fields{i})
            error('Chemistry structure must contain %s field', required_fields{i});
        end
    end

    %% Extract data vectors
    rho = chemistry.rho(:);
    e = chemistry.e(:);
    T = chemistry.T(:);
    gamma_star = chemistry.gamma_star(:);
    cv_star = chemistry.cv_star(:);
    mu = chemistry.mu(:);
    k = chemistry.k(:);

    % Check for optional sound speed data
    if chemistry.has_sound_speed
        a = chemistry.a(:);
        has_a = true;
        fprintf('Sound speed data detected - will be included in fits.\n');
    else
        a = [];
        has_a = false;
        fprintf('Sound speed data not available - skipping sound speed fits.\n');
    end

    % Check for optional entropy data
    if chemistry.has_entropy
        s = chemistry.s(:);
        has_s = true;
        fprintf('Entropy data detected - will be included in fits.\n');
    else
        s = [];
        has_s = false;
        fprintf('Entropy data not available - skipping entropy fits.\n');
    end

    n_points = length(rho);
    fprintf('Fitting %d data points...\n', n_points);

    %% Apply logarithmic transformations
    fprintf('Applying logarithmic transformations: log10(rho), log10(e), log10(T)...\n');
    log_rho = log10(rho);
    log_e = log10(e);
    log_T = log10(T);

    %% Remove invalid data points
    if has_a && has_s
        valid_idx = ~(isnan(log_rho) | isnan(log_e) | isnan(log_T) | isnan(gamma_star) | ...
                      isnan(cv_star) | isnan(mu) | isnan(k) | isnan(a) | isnan(s) | ...
                      isinf(log_rho) | isinf(log_e) | isinf(log_T) | isinf(gamma_star) | ...
                      isinf(cv_star) | isinf(mu) | isinf(k) | isinf(a) | isinf(s));
    elseif has_a && ~has_s
        valid_idx = ~(isnan(log_rho) | isnan(log_e) | isnan(log_T) | isnan(gamma_star) | ...
                      isnan(cv_star) | isnan(mu) | isnan(k) | isnan(a) | ...
                      isinf(log_rho) | isinf(log_e) | isinf(log_T) | isinf(gamma_star) | ...
                      isinf(cv_star) | isinf(mu) | isinf(k) | isinf(a));
    elseif ~has_a && has_s
        valid_idx = ~(isnan(log_rho) | isnan(log_e) | isnan(log_T) | isnan(gamma_star) | ...
                      isnan(cv_star) | isnan(mu) | isnan(k) | isnan(s) | ...
                      isinf(log_rho) | isinf(log_e) | isinf(log_T) | isinf(gamma_star) | ...
                      isinf(cv_star) | isinf(mu) | isinf(k) | isinf(s));
    else
        valid_idx = ~(isnan(log_rho) | isnan(log_e) | isnan(log_T) | isnan(gamma_star) | ...
                      isnan(cv_star) | isnan(mu) | isnan(k) | ...
                      isinf(log_rho) | isinf(log_e) | isinf(log_T) | isinf(gamma_star) | ...
                      isinf(cv_star) | isinf(mu) | isinf(k));
    end

    if sum(valid_idx) < n_points
        fprintf('Warning: Removed %d invalid data points for basic properties\n', n_points - sum(valid_idx));
        rho = rho(valid_idx);
        log_rho = log_rho(valid_idx);
        e = e(valid_idx);
        log_e = log_e(valid_idx);
        T = T(valid_idx);
        log_T = log_T(valid_idx);
        gamma_star = gamma_star(valid_idx);
        cv_star = cv_star(valid_idx);
        mu = mu(valid_idx);
        k = k(valid_idx);
        if has_a
            a = a(valid_idx);
        end
        if has_s
            s = s(valid_idx);
        end
        n_points = length(rho);
    end

    %% Store data ranges for scaling
    chemistry.fit_info.rho_range = [min(rho), max(rho)];
    chemistry.fit_info.log_rho_range = [min(log_rho), max(log_rho)];
    chemistry.fit_info.e_range = [min(e), max(e)];
    chemistry.fit_info.log_e_range = [min(log_e), max(log_e)];
    chemistry.fit_info.T_range = [min(T), max(T)];
    chemistry.fit_info.log_T_range = [min(log_T), max(log_T)];
    chemistry.fit_info.gamma_star_range = [min(gamma_star), max(gamma_star)];
    chemistry.fit_info.cv_star_range = [min(cv_star), max(cv_star)];
    chemistry.fit_info.mu_range = [min(mu), max(mu)];
    chemistry.fit_info.k_range = [min(k), max(k)];
    if has_a
        chemistry.fit_info.a_range = [min(a), max(a)];
    end
    if has_s
        chemistry.fit_info.s_range = [min(s), max(s)];
    end
    chemistry.fit_info.n_points = n_points;

    %% Print data ranges
    fprintf('Original data ranges:\n');
    fprintf('  rho:        [%.3e, %.3e] kg/m^3\n', chemistry.fit_info.rho_range);
    fprintf('  log10(rho): [%.3f, %.3f]\n', chemistry.fit_info.log_rho_range);
    fprintf('  e:          [%.3e, %.3e] J/kg\n', chemistry.fit_info.e_range);
    fprintf('  log10(e):   [%.3f, %.3f]\n', chemistry.fit_info.log_e_range);
    fprintf('  T:          [%.1f, %.1f] K\n', chemistry.fit_info.T_range);
    fprintf('  log10(T):   [%.3f, %.3f]\n', chemistry.fit_info.log_T_range);
    fprintf('  gamma*:     [%.4f, %.4f]\n', chemistry.fit_info.gamma_star_range);
    fprintf('  cv*:        [%.1f, %.1f] J/kg/K\n', chemistry.fit_info.cv_star_range);
    fprintf('  mu:         [%.2e, %.2e] Pa*s\n', chemistry.fit_info.mu_range);
    fprintf('  k:          [%.4f, %.4f] W/m/K\n', chemistry.fit_info.k_range);
    if has_a
        fprintf('  a:          [%.1f, %.1f] m/s\n', chemistry.fit_info.a_range);
    end
    if has_s
        fprintf('  s:          [%.1f, %.1f] J/kg/K\n', chemistry.fit_info.s_range);
    end

    %% Scale all data to [-1, 1]
    if chemistry.detailed_output
        fprintf('\nScaling all data to [-1, 1] for improved interpolation...\n');
    end

    % Scale input variables
    log_rho_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_rho, chemistry.fit_info.log_rho_range);
    log_e_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_e, chemistry.fit_info.log_e_range);

    % Scale output variables
    T_scaled = SCALE_TO_MINUS_ONE_TO_ONE(T, chemistry.fit_info.T_range);
    gamma_star_scaled = SCALE_TO_MINUS_ONE_TO_ONE(gamma_star, chemistry.fit_info.gamma_star_range);
    cv_star_scaled = SCALE_TO_MINUS_ONE_TO_ONE(cv_star, chemistry.fit_info.cv_star_range);
    mu_scaled = SCALE_TO_MINUS_ONE_TO_ONE(mu, chemistry.fit_info.mu_range);
    k_scaled = SCALE_TO_MINUS_ONE_TO_ONE(k, chemistry.fit_info.k_range);
    if has_a
        a_scaled = SCALE_TO_MINUS_ONE_TO_ONE(a, chemistry.fit_info.a_range);
    end
    if has_s
        s_scaled = SCALE_TO_MINUS_ONE_TO_ONE(s, chemistry.fit_info.s_range);
    end

    %% Create regular grid for gridded interpolants
    n_grid_log_rho = 100;
    n_grid_log_e = 100;

    if chemistry.detailed_output
        fprintf('\nCreating regular grid (%dx%d = %d points)...\n', n_grid_log_rho, n_grid_log_e, n_grid_log_rho * n_grid_log_e);
    end

    log_rho_grid_vec = linspace(-1, 1, n_grid_log_rho);
    log_e_grid_vec = linspace(-1, 1, n_grid_log_e);

    % Use ndgrid (not meshgrid) for griddedInterpolant compatibility
    [log_e_grid, log_rho_grid] = ndgrid(log_e_grid_vec, log_rho_grid_vec);

    % Interpolation configuration
    griddata_method = 'cubic';
    interp_method = 'cubic';
    extrap_method = 'linear';

    if chemistry.detailed_output
        fprintf('Using griddata method: %s, griddedInterpolant method: %s\n', griddata_method, interp_method);
    end
    
    %% Interpolate scattered data onto regular grid
    fprintf('\nCreating griddedInterpolant objects for each property...\n');

    T_grid = griddata(log_rho_scaled, log_e_scaled, T_scaled, log_rho_grid, log_e_grid, griddata_method);
    fprintf('  T: done\n');

    gamma_star_grid = griddata(log_rho_scaled, log_e_scaled, gamma_star_scaled, log_rho_grid, log_e_grid, griddata_method);
    fprintf('  gamma_star: done\n');

    cv_star_grid = griddata(log_rho_scaled, log_e_scaled, cv_star_scaled, log_rho_grid, log_e_grid, griddata_method);
    fprintf('  cv_star: done\n');

    mu_grid = griddata(log_rho_scaled, log_e_scaled, mu_scaled, log_rho_grid, log_e_grid, griddata_method);
    fprintf('  mu: done\n');

    k_grid = griddata(log_rho_scaled, log_e_scaled, k_scaled, log_rho_grid, log_e_grid, griddata_method);
    fprintf('  k: done\n');

    if has_a
        a_grid = griddata(log_rho_scaled, log_e_scaled, a_scaled, log_rho_grid, log_e_grid, griddata_method);
        fprintf('  a (sound speed): done\n');
    end

    if has_s
        s_grid = griddata(log_rho_scaled, log_e_scaled, s_scaled, log_rho_grid, log_e_grid, griddata_method);
        fprintf('  s (entropy): done\n');
    end

    %% Create griddedInterpolant objects
    chemistry.fit_T = griddedInterpolant({log_e_grid_vec, log_rho_grid_vec}, T_grid, interp_method, extrap_method);
    chemistry.fit_gamma_star = griddedInterpolant({log_e_grid_vec, log_rho_grid_vec}, gamma_star_grid, interp_method, extrap_method);
    chemistry.fit_cv_star = griddedInterpolant({log_e_grid_vec, log_rho_grid_vec}, cv_star_grid, interp_method, extrap_method);
    chemistry.fit_mu = griddedInterpolant({log_e_grid_vec, log_rho_grid_vec}, mu_grid, interp_method, extrap_method);
    chemistry.fit_k = griddedInterpolant({log_e_grid_vec, log_rho_grid_vec}, k_grid, interp_method, extrap_method);
    if has_a
        chemistry.fit_a = griddedInterpolant({log_e_grid_vec, log_rho_grid_vec}, a_grid, interp_method, extrap_method);
    end
    if has_s
        chemistry.fit_s = griddedInterpolant({log_e_grid_vec, log_rho_grid_vec}, s_grid, interp_method, extrap_method);
    end

    %% Create inverse interpolant: e(T, rho)
    fprintf('  e: done\n');

    log_T_grid_vec = linspace(-1, 1, n_grid_log_e);
    log_rho_grid_vec_inv = linspace(-1, 1, n_grid_log_rho);

    log_T_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_T, chemistry.fit_info.log_T_range);

    [log_T_grid_inv, log_rho_grid_inv] = ndgrid(log_T_grid_vec, log_rho_grid_vec_inv);

    log_e_grid_inv = griddata(log_T_scaled, log_rho_scaled, log_e_scaled, log_T_grid_inv, log_rho_grid_inv, griddata_method);

    chemistry.fit_log_e_from_T_rho = griddedInterpolant({log_T_grid_vec, log_rho_grid_vec_inv}, log_e_grid_inv, interp_method, extrap_method);

    %% Create species interpolants
    if chemistry.has_molar_data
        if chemistry.detailed_output
            fprintf('\nSpecies data detected - creating interpolants for species concentrations...\n');
        end

        chemistry.species_fits = struct();
        chemistry.species_ranges = struct();

        for i = 1:length(chemistry.species_list)
            species = chemistry.species_list{i};
            log10_field = ['log10_' species];

            if isfield(chemistry, log10_field)
                log10_conc = chemistry.(log10_field)(valid_idx);

                % Additional filtering for species data
                species_valid_idx = ~(isnan(log10_conc) | isinf(log10_conc) | log10_conc < -50);

                if sum(species_valid_idx) > 10
                    % Store range for this species
                    range_name = ['range_log10_' species];
                    chemistry.species_ranges.(range_name) = [min(log10_conc(species_valid_idx)), max(log10_conc(species_valid_idx))];

                    % Scale species data to [-1, 1]
                    log10_conc_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log10_conc(species_valid_idx), chemistry.species_ranges.(range_name));

                    % Interpolate species data onto regular grid
                    log10_conc_grid = griddata(log_rho_scaled(species_valid_idx), log_e_scaled(species_valid_idx), log10_conc_scaled, log_rho_grid, log_e_grid, griddata_method);

                    % Create griddedInterpolant for this species
                    fit_name = ['fit_log10_' species];
                    chemistry.species_fits.(fit_name) = griddedInterpolant({log_e_grid_vec, log_rho_grid_vec}, log10_conc_grid, interp_method, extrap_method);
                    fprintf('  X_%s: done\n', species);

                    if chemistry.detailed_output
                        fprintf('  X_%s: %d valid points, range [%.2f, %.2f]\n', species, sum(species_valid_idx), chemistry.species_ranges.(range_name));
                    end
                else
                    fprintf('    %s: insufficient data (%d points), skipping\n', species, sum(species_valid_idx));
                end
            end
        end
    end

    %% Store grid information
    chemistry.fit_info.grid_size = [n_grid_log_e, n_grid_log_rho];
    chemistry.fit_info.log_rho_grid_vec = log_rho_grid_vec;
    chemistry.fit_info.log_e_grid_vec = log_e_grid_vec;

    %% Create evaluation function handles
    chemistry.eval_T = @(rho_new, e_new) EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(rho_new, e_new, chemistry.fit_T, chemistry.fit_info.log_rho_range, chemistry.fit_info.log_e_range, chemistry.fit_info.T_range);
    chemistry.eval_gamma_star = @(rho_new, e_new) EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(rho_new, e_new, chemistry.fit_gamma_star, chemistry.fit_info.log_rho_range, chemistry.fit_info.log_e_range, chemistry.fit_info.gamma_star_range);
    chemistry.eval_cv_star = @(rho_new, e_new) EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(rho_new, e_new, chemistry.fit_cv_star, chemistry.fit_info.log_rho_range, chemistry.fit_info.log_e_range, chemistry.fit_info.cv_star_range);
    chemistry.eval_mu = @(rho_new, e_new) EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(rho_new, e_new, chemistry.fit_mu, chemistry.fit_info.log_rho_range, chemistry.fit_info.log_e_range, chemistry.fit_info.mu_range);
    chemistry.eval_k = @(rho_new, e_new) EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(rho_new, e_new, chemistry.fit_k, chemistry.fit_info.log_rho_range, chemistry.fit_info.log_e_range, chemistry.fit_info.k_range);
    if has_a
        chemistry.eval_a = @(rho_new, e_new) EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(rho_new, e_new, chemistry.fit_a, chemistry.fit_info.log_rho_range, chemistry.fit_info.log_e_range, chemistry.fit_info.a_range);
    end
    if has_s
        chemistry.eval_s = @(rho_new, e_new) EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(rho_new, e_new, chemistry.fit_s, chemistry.fit_info.log_rho_range, chemistry.fit_info.log_e_range, chemistry.fit_info.s_range);
    end

    % Inverse evaluation function: E(T, rho)
    chemistry.eval_e = @(T_new, rho_new) EVAL_GRIDDED_INTERPOLANT_INVERSE_LOG_INPUTS(T_new, rho_new, chemistry.fit_log_e_from_T_rho, chemistry.fit_info.log_T_range, chemistry.fit_info.log_rho_range, chemistry.fit_info.log_e_range, chemistry.fit_info.e_range);

    %% Create species evaluation functions
    if chemistry.has_molar_data
        chemistry.species_eval_info = {};

        for i = 1:length(chemistry.species_list)
            species = chemistry.species_list{i};
            fit_name = ['fit_log10_' species];
            eval_name = ['eval_log10_' species];
            range_name = ['range_log10_' species];

            if isfield(chemistry.species_fits, fit_name) && isfield(chemistry.species_ranges, range_name)
                eval_info = struct();
                eval_info.species = species;
                eval_info.fit_object = chemistry.species_fits.(fit_name);
                eval_info.log_rho_range = chemistry.fit_info.log_rho_range;
                eval_info.log_e_range = chemistry.fit_info.log_e_range;
                eval_info.species_range = chemistry.species_ranges.(range_name);

                chemistry.species_eval_info{end+1} = eval_info;

                func_str = sprintf('chemistry.%s = @(rho_new, e_new) EVAL_SPECIES_GRIDDED_LOG_INPUTS(chemistry, ''%s'', rho_new, e_new);', eval_name, species);
                eval(func_str);
            end
        end

        chemistry.eval_all_species = @(rho_new, e_new) EVALUATE_ALL_SPECIES(chemistry, rho_new, e_new);
    end

    % Vectorized evaluation for all basic properties
    chemistry.eval_all = @(rho_new, e_new) EVALUATE_ALL_PROPERTIES(chemistry, rho_new, e_new, has_a);

    %% Validate fits on original data
    if chemistry.detailed_output
        fprintf('Validating fits on original data...\n');
        chemistry.fit_validation = VALIDATE_FITS(chemistry, rho, e, T, gamma_star, cv_star, mu, k, a, s);
    end

    %% Store fitting metadata
    chemistry.fit_info.griddata_method = griddata_method;
    chemistry.fit_info.interp_method = interp_method;
    chemistry.fit_info.extrapolation = extrap_method;
    chemistry.fit_info.fit_time = toc;
    chemistry.fit_info.date_fitted = datestr(now);
    chemistry.fit_info.has_sound_speed_fit = has_a;
    chemistry.fit_info.has_entropy_fit = has_s;
    chemistry.fit_info.interpolant_type = 'griddedInterpolant_with_log_rho_log_e_log_T';
    chemistry.fit_info.uses_log_rho = true;
    chemistry.fit_info.uses_log_e = true;
    chemistry.fit_info.uses_log_T = true;

    %% Print validation summary
    if chemistry.detailed_output
        fprintf('Fitting completed in %.2f seconds\n', chemistry.fit_info.fit_time);
        fprintf('Validation results for basic properties:\n');
        fields = fieldnames(chemistry.fit_validation);
        for i = 1:length(fields)
            if contains(fields{i}, 'rmse')
                fprintf('  %s: %.6e\n', fields{i}, chemistry.fit_validation.(fields{i}));
            elseif contains(fields{i}, 'r2')
                fprintf('  %s: %.6f\n', fields{i}, chemistry.fit_validation.(fields{i}));
            end
        end

        if chemistry.has_molar_data
            fprintf('Species fitting completed for: %s\n', ...
                strjoin(fieldnames(chemistry.species_fits), ', '));
        end

        fprintf('\nData approach: Scattered data interpolated onto %dx%d regular grid.\n', n_grid_log_e, n_grid_log_rho);
        fprintf('Data transformation: log10(rho), log10(e), log10(T) applied before [-1,1] scaling.\n');
        fprintf('Data scaling: All variables scaled to [-1,1] for improved numerical accuracy.\n');
        fprintf('Using griddedInterpolant with %s interpolation for fast, smooth evaluation.\n', interp_method);
        fprintf('Inverse fit E(T,rho) created in log-space for Rankine-Hugoniot solver integration.\n');
        if has_s
            fprintf('Entropy s(rho,E) interpolation included.\n');
        end
        fprintf('\nNote: log transformations improve interpolation over wide ranges.\n');
    end
end


%% ========================================================================
%  EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS - Evaluate with log transforms
%  ========================================================================
function result = EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(rho_new, e_new, interpolant, log_rho_range, log_e_range, values_range)
% EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS - Evaluate property using log-space interpolation
%
% Transforms query points to log10 space, scales to [-1, 1], evaluates
% the griddedInterpolant, and unscales the result to original units.
% Handles non-positive rho or e values by extrapolating from minimum
% log values instead of returning NaN.
%
% Inputs:
%   rho_new      - Query densities in original units (kg/m^3)
%   e_new        - Query energies in original units (J/kg)
%   interpolant  - griddedInterpolant object (expects scaled log inputs)
%   log_rho_range - [min, max] of log10(rho) for scaling
%   log_e_range   - [min, max] of log10(e) for scaling
%   values_range  - [min, max] of output variable for unscaling
%
% Outputs:
%   result - Interpolated values in original units

    % Handle non-positive values with extrapolation instead of NaN
    result = zeros(size(rho_new));
    valid_mask = (rho_new > 0) & (e_new > 0);

    if ~any(valid_mask(:))
        % All values are non-positive - use minimum log ranges
        log_rho_min = log_rho_range(1);
        log_e_min = log_e_range(1);
        log_rho_new_extrap = log_rho_min * ones(size(rho_new));
        log_e_new_extrap = log_e_min * ones(size(e_new));
        log_rho_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_rho_new_extrap, log_rho_range);
        log_e_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_e_new_extrap, log_e_range);
        result_scaled = interpolant(log_e_new_scaled, log_rho_new_scaled);
        result = SCALE_FROM_MINUS_ONE_TO_ONE(result_scaled, values_range);
        return;
    end

    % Process valid positive rho and e values
    if any(valid_mask(:))
        log_rho_new = log10(rho_new(valid_mask));
        log_e_new = log10(e_new(valid_mask));

        log_rho_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_rho_new, log_rho_range);
        log_e_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_e_new, log_e_range);

        result_scaled = interpolant(log_e_new_scaled, log_rho_new_scaled);
        result(valid_mask) = SCALE_FROM_MINUS_ONE_TO_ONE(result_scaled, values_range);
    end

    % For invalid values, extrapolate using minimum log values
    if any(~valid_mask(:))
        log_rho_min = log_rho_range(1);
        log_e_min = log_e_range(1);

        log_rho_extrap = log_rho_min * ones(sum(~valid_mask(:)), 1);
        log_e_extrap = log_e_min * ones(sum(~valid_mask(:)), 1);

        % Override with actual log values where possible
        invalid_rho = rho_new(~valid_mask);
        invalid_e = e_new(~valid_mask);
        valid_rho_only = invalid_rho > 0;
        valid_e_only = invalid_e > 0;

        if any(valid_rho_only)
            log_rho_extrap(valid_rho_only) = log10(invalid_rho(valid_rho_only));
        end
        if any(valid_e_only)
            log_e_extrap(valid_e_only) = log10(invalid_e(valid_e_only));
        end

        log_rho_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_rho_extrap, log_rho_range);
        log_e_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_e_extrap, log_e_range);

        result_scaled = interpolant(log_e_new_scaled, log_rho_new_scaled);
        result(~valid_mask) = SCALE_FROM_MINUS_ONE_TO_ONE(result_scaled, values_range);
    end
end


%% ========================================================================
%  EVAL_GRIDDED_INTERPOLANT_INVERSE_LOG_INPUTS - Inverse evaluation e(T, rho)
%  ========================================================================
function result = EVAL_GRIDDED_INTERPOLANT_INVERSE_LOG_INPUTS(T_new, rho_new, interpolant, log_T_range, log_rho_range, log_e_range, e_range)
% EVAL_GRIDDED_INTERPOLANT_INVERSE_LOG_INPUTS - Evaluate e(T, rho) in log-space
%
% Evaluates the inverse interpolant that maps (T, rho) -> e using
% log10 transformations on all variables. Returns energy in original
% units (not log-transformed).
%
% Inputs:
%   T_new         - Query temperatures in original units (K)
%   rho_new       - Query densities in original units (kg/m^3)
%   interpolant   - griddedInterpolant for log_e(log_T, log_rho)
%   log_T_range   - [min, max] of log10(T) for scaling
%   log_rho_range - [min, max] of log10(rho) for scaling
%   log_e_range   - [min, max] of log10(e) for scaling
%   e_range       - [min, max] of e for reference (unused, kept for API)
%
% Outputs:
%   result - Interpolated energy values in original units (J/kg)

    % Handle non-positive values with extrapolation
    result = zeros(size(rho_new));
    valid_mask = (T_new > 0) & (rho_new > 0);

    if ~any(valid_mask(:))
        % All values are non-positive - use minimum log ranges
        log_T_min = log_T_range(1);
        log_rho_min = log_rho_range(1);
        log_T_new_extrap = log_T_min * ones(size(T_new));
        log_rho_new_extrap = log_rho_min * ones(size(rho_new));
        log_T_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_T_new_extrap, log_T_range);
        log_rho_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_rho_new_extrap, log_rho_range);
        log_e_result_scaled = interpolant(log_T_new_scaled, log_rho_new_scaled);
        log_e_result = SCALE_FROM_MINUS_ONE_TO_ONE(log_e_result_scaled, log_e_range);
        result = 10.^log_e_result;
        return;
    end

    % Process valid positive T and rho values
    if any(valid_mask(:))
        log_T_new = log10(T_new(valid_mask));
        log_rho_new = log10(rho_new(valid_mask));

        log_T_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_T_new, log_T_range);
        log_rho_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_rho_new, log_rho_range);

        log_e_result_scaled = interpolant(log_T_new_scaled, log_rho_new_scaled);
        log_e_result = SCALE_FROM_MINUS_ONE_TO_ONE(log_e_result_scaled, log_e_range);
        result(valid_mask) = 10.^log_e_result;
    end

    % For invalid values, extrapolate using minimum log values
    if any(~valid_mask(:))
        log_T_min = log_T_range(1);
        log_rho_min = log_rho_range(1);

        log_T_extrap = log_T_min * ones(sum(~valid_mask(:)), 1);
        log_rho_extrap = log_rho_min * ones(sum(~valid_mask(:)), 1);

        % Override with actual log values where possible
        invalid_T = T_new(~valid_mask);
        invalid_rho = rho_new(~valid_mask);
        valid_T_only = invalid_T > 0;
        valid_rho_only = invalid_rho > 0;

        if any(valid_T_only)
            log_T_extrap(valid_T_only) = log10(invalid_T(valid_T_only));
        end
        if any(valid_rho_only)
            log_rho_extrap(valid_rho_only) = log10(invalid_rho(valid_rho_only));
        end

        log_T_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_T_extrap, log_T_range);
        log_rho_new_scaled = SCALE_TO_MINUS_ONE_TO_ONE(log_rho_extrap, log_rho_range);

        log_e_result_scaled = interpolant(log_T_new_scaled, log_rho_new_scaled);
        log_e_result = SCALE_FROM_MINUS_ONE_TO_ONE(log_e_result_scaled, log_e_range);
        result(~valid_mask) = 10.^log_e_result;
    end
end


%% ========================================================================
%  EVAL_SPECIES_GRIDDED_LOG_INPUTS - Evaluate species concentration
%  ========================================================================
function result = EVAL_SPECIES_GRIDDED_LOG_INPUTS(chemistry, species_name, rho_new, e_new)
% EVAL_SPECIES_GRIDDED_LOG_INPUTS - Evaluate log10 species concentration
%
% Looks up the species-specific interpolant and evaluates it using
% EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS. Returns NaN if species not found.
%
% Inputs:
%   chemistry    - Chemistry structure with species_eval_info
%   species_name - Species identifier string (e.g., 'O2', 'CO2')
%   rho_new      - Query densities (kg/m^3)
%   e_new        - Query energies (J/kg)
%
% Outputs:
%   result - log10 of species molar concentration

    for i = 1:length(chemistry.species_eval_info)
        if strcmp(chemistry.species_eval_info{i}.species, species_name)
            eval_info = chemistry.species_eval_info{i};
            result = EVAL_GRIDDED_INTERPOLANT_LOG_INPUTS(rho_new, e_new, eval_info.fit_object, ...
                eval_info.log_rho_range, eval_info.log_e_range, eval_info.species_range);
            return;
        end
    end
    % Species not found
    result = NaN(size(rho_new));
end


%% ========================================================================
%  SCALE_TO_MINUS_ONE_TO_ONE - Scale data to [-1, 1]
%  ========================================================================
function scaled_data = SCALE_TO_MINUS_ONE_TO_ONE(data, data_range)
% SCALE_TO_MINUS_ONE_TO_ONE - Linear scaling from original range to [-1, 1]
%
% Inputs:
%   data       - Data array to scale
%   data_range - [min, max] of the original data range
%
% Outputs:
%   scaled_data - Data scaled to [-1, 1]

    if data_range(2) == data_range(1)
        scaled_data = zeros(size(data));
    else
        scaled_data = 2 * (data - data_range(1)) / (data_range(2) - data_range(1)) - 1;
    end
end


%% ========================================================================
%  SCALE_FROM_MINUS_ONE_TO_ONE - Unscale data from [-1, 1]
%  ========================================================================
function original_data = SCALE_FROM_MINUS_ONE_TO_ONE(scaled_data, data_range)
% SCALE_FROM_MINUS_ONE_TO_ONE - Linear unscaling from [-1, 1] to original range
%
% Inputs:
%   scaled_data - Data in [-1, 1] range
%   data_range  - [min, max] of the target data range
%
% Outputs:
%   original_data - Data in original units

    if data_range(2) == data_range(1)
        original_data = data_range(1) * ones(size(scaled_data));
    else
        original_data = (scaled_data + 1) * (data_range(2) - data_range(1)) / 2 + data_range(1);
    end
end


%% ========================================================================
%  EVALUATE_ALL_PROPERTIES - Batch evaluate all basic properties
%  ========================================================================
function result = EVALUATE_ALL_PROPERTIES(chemistry, rho_new, e_new, has_a)
% EVALUATE_ALL_PROPERTIES - Evaluate T, gamma*, cv*, mu, k (and a) at once
%
% Inputs:
%   chemistry - Fitted chemistry structure
%   rho_new   - Query densities (kg/m^3)
%   e_new     - Query energies (J/kg)
%   has_a     - Logical flag for sound speed availability
%
% Outputs:
%   result - Structure with fields .T, .gamma_star, .cv_star, .mu, .k, .a

    result.T = chemistry.eval_T(rho_new, e_new);
    result.gamma_star = chemistry.eval_gamma_star(rho_new, e_new);
    result.cv_star = chemistry.eval_cv_star(rho_new, e_new);
    result.mu = chemistry.eval_mu(rho_new, e_new);
    result.k = chemistry.eval_k(rho_new, e_new);
    if has_a
        result.a = chemistry.eval_a(rho_new, e_new);
    end
end


%% ========================================================================
%  EVALUATE_ALL_SPECIES - Batch evaluate all species concentrations
%  ========================================================================
function result = EVALUATE_ALL_SPECIES(chemistry, rho_new, e_new)
% EVALUATE_ALL_SPECIES - Evaluate all species log10 concentrations at once
%
% Inputs:
%   chemistry - Fitted chemistry structure with species evaluators
%   rho_new   - Query densities (kg/m^3)
%   e_new     - Query energies (J/kg)
%
% Outputs:
%   result - Structure with one field per species containing log10 concentrations

    result = struct();

    for i = 1:length(chemistry.species_list)
        species = chemistry.species_list{i};
        eval_name = ['eval_log10_' species];

        if isfield(chemistry, eval_name)
            result.(species) = chemistry.(eval_name)(rho_new, e_new);
        else
            result.(species) = NaN;
        end
    end
end


%% ========================================================================
%  VALIDATE_FITS - Compute error metrics for all fitted properties
%  ========================================================================
function validation = VALIDATE_FITS(chemistry, rho, e, T_true, gamma_star_true, cv_star_true, mu_true, k_true, a_true, s_true)
% VALIDATE_FITS - Validate interpolation fits against original data
%
% Computes RMSE, MAE, and R^2 for each fitted property by comparing
% interpolated predictions against the original tabulated data. For
% entropy, performs both direct interpolation and thermodynamic
% integration validation.
%
% Inputs:
%   chemistry       - Fitted chemistry structure
%   rho, e          - Original density and energy data
%   T_true          - Original temperature data (K)
%   gamma_star_true - Original gamma* data
%   cv_star_true    - Original cv* data (J/kg/K)
%   mu_true         - Original viscosity data (Pa*s)
%   k_true          - Original thermal conductivity data (W/m/K)
%   a_true          - Original sound speed data (m/s), or empty
%   s_true          - Original entropy data (J/kg/K), or empty
%
% Outputs:
%   validation - Structure with RMSE, MAE, R^2 for each property and
%                predicted values in validation.predictions

    %% Predict using fits
    T_pred = chemistry.eval_T(rho, e);
    gamma_star_pred = chemistry.eval_gamma_star(rho, e);
    cv_star_pred = chemistry.eval_cv_star(rho, e);
    mu_pred = chemistry.eval_mu(rho, e);
    k_pred = chemistry.eval_k(rho, e);
    e_pred = chemistry.eval_e(T_true, rho);

    %% Create valid-data masks
    valid_T = ~isnan(T_pred);
    valid_gamma_star = ~isnan(gamma_star_pred);
    valid_cv_star = ~isnan(cv_star_pred);
    valid_mu = ~isnan(mu_pred);
    valid_k = ~isnan(k_pred);
    valid_e = ~isnan(e_pred);

    % Helper for safe total sum of squares
    safe_sst = @(true_data, valid_mask) sum((true_data(valid_mask) - mean(true_data(valid_mask))).^2);
    sst_threshold = 1e-12;

    %% Temperature validation
    validation.T_rmse = sqrt(mean((T_pred(valid_T) - T_true(valid_T)).^2));
    validation.T_mae = mean(abs(T_pred(valid_T) - T_true(valid_T)));
    sst_T = safe_sst(T_true, valid_T);
    if sst_T < sst_threshold
        validation.T_r2 = 1.0;
    else
        validation.T_r2 = 1 - sum((T_pred(valid_T) - T_true(valid_T)).^2) / sst_T;
    end

    %% gamma_star validation
    validation.gamma_star_rmse = sqrt(mean((gamma_star_pred(valid_gamma_star) - gamma_star_true(valid_gamma_star)).^2));
    validation.gamma_star_mae = mean(abs(gamma_star_pred(valid_gamma_star) - gamma_star_true(valid_gamma_star)));
    sst_gamma_star = safe_sst(gamma_star_true, valid_gamma_star);
    if sst_gamma_star < sst_threshold
        validation.gamma_star_r2 = 1.0;
    else
        validation.gamma_star_r2 = 1 - sum((gamma_star_pred(valid_gamma_star) - gamma_star_true(valid_gamma_star)).^2) / sst_gamma_star;
    end

    %% cv_star validation
    validation.cv_star_rmse = sqrt(mean((cv_star_pred(valid_cv_star) - cv_star_true(valid_cv_star)).^2));
    validation.cv_star_mae = mean(abs(cv_star_pred(valid_cv_star) - cv_star_true(valid_cv_star)));
    sst_cv_star = safe_sst(cv_star_true, valid_cv_star);
    if sst_cv_star < sst_threshold
        validation.cv_star_r2 = 1.0;
    else
        validation.cv_star_r2 = 1 - sum((cv_star_pred(valid_cv_star) - cv_star_true(valid_cv_star)).^2) / sst_cv_star;
    end

    %% Viscosity validation
    validation.mu_rmse = sqrt(mean((mu_pred(valid_mu) - mu_true(valid_mu)).^2));
    validation.mu_mae = mean(abs(mu_pred(valid_mu) - mu_true(valid_mu)));
    sst_mu = safe_sst(mu_true, valid_mu);
    if sst_mu < sst_threshold
        validation.mu_r2 = 1.0;
    else
        validation.mu_r2 = 1 - sum((mu_pred(valid_mu) - mu_true(valid_mu)).^2) / sst_mu;
    end

    %% Thermal conductivity validation
    validation.k_rmse = sqrt(mean((k_pred(valid_k) - k_true(valid_k)).^2));
    validation.k_mae = mean(abs(k_pred(valid_k) - k_true(valid_k)));
    sst_k = safe_sst(k_true, valid_k);
    if sst_k < sst_threshold
        validation.k_r2 = 1.0;
    else
        validation.k_r2 = 1 - sum((k_pred(valid_k) - k_true(valid_k)).^2) / sst_k;
    end

    %% Inverse fit E(T, rho) validation
    validation.e_rmse = sqrt(mean((e_pred(valid_e) - e(valid_e)).^2));
    validation.e_mae = mean(abs(e_pred(valid_e) - e(valid_e)));
    sst_e = safe_sst(e, valid_e);
    if sst_e < sst_threshold
        validation.e_r2 = 1.0;
    else
        validation.e_r2 = 1 - sum((e_pred(valid_e) - e(valid_e)).^2) / sst_e;
    end

    %% Sound speed validation (if available)
    if ~isempty(a_true)
        a_pred = chemistry.eval_a(rho, e);
        valid_a = ~isnan(a_pred);

        validation.a_rmse = sqrt(mean((a_pred(valid_a) - a_true(valid_a)).^2));
        validation.a_mae = mean(abs(a_pred(valid_a) - a_true(valid_a)));
        sst_a = safe_sst(a_true, valid_a);
        if sst_a < sst_threshold
            validation.a_r2 = 1.0;
        else
            validation.a_r2 = 1 - sum((a_pred(valid_a) - a_true(valid_a)).^2) / sst_a;
        end
        validation.predictions.a = a_pred;
    end

    %% Entropy validation (if available)
    if ~isempty(s_true)
        fprintf('\n=== Entropy Validation ===\n');

        % Method 1: Direct interpolation
        fprintf('Method 1: Direct interpolation s(rho,e)...\n');
        s_pred_interp = chemistry.eval_s(rho, e);
        valid_s_interp = ~isnan(s_pred_interp);

        validation.s_interp_rmse = sqrt(mean((s_pred_interp(valid_s_interp) - s_true(valid_s_interp)).^2));
        validation.s_interp_mae = mean(abs(s_pred_interp(valid_s_interp) - s_true(valid_s_interp)));
        sst_s_interp = safe_sst(s_true, valid_s_interp);
        if sst_s_interp < sst_threshold
            validation.s_interp_r2 = 1.0;
        else
            validation.s_interp_r2 = 1 - sum((s_pred_interp(valid_s_interp) - s_true(valid_s_interp)).^2) / sst_s_interp;
        end
        validation.predictions.s_interp = s_pred_interp;

        fprintf('  Interpolation R^2: %.6f\n', validation.s_interp_r2);
        fprintf('  Interpolation RMSE: %.4f J/kg/K\n', validation.s_interp_rmse);

        % Method 2: Thermodynamic integration
        fprintf('Method 2: Thermodynamic integration of ds = -(gamma*-1)cv*/(rho*e)drho + de/T...\n');
        tic;
        s_pred_thermo = INTEGRATE_ENTROPY_THERMODYNAMIC(chemistry, rho, e);
        t_integration = toc;
        fprintf('  Integration completed in %.3f seconds\n', t_integration);

        valid_s_thermo = ~isnan(s_pred_thermo) & ~isinf(s_pred_thermo);

        if sum(valid_s_thermo) > 0
            validation.s_thermo_rmse = sqrt(mean((s_pred_thermo(valid_s_thermo) - s_true(valid_s_thermo)).^2));
            validation.s_thermo_mae = mean(abs(s_pred_thermo(valid_s_thermo) - s_true(valid_s_thermo)));
            sst_s_thermo = safe_sst(s_true, valid_s_thermo);
            if sst_s_thermo < sst_threshold
                validation.s_thermo_r2 = 1.0;
            else
                validation.s_thermo_r2 = 1 - sum((s_pred_thermo(valid_s_thermo) - s_true(valid_s_thermo)).^2) / sst_s_thermo;
            end
            validation.predictions.s_thermo = s_pred_thermo;

            fprintf('  Thermodynamic R^2: %.6f\n', validation.s_thermo_r2);
            fprintf('  Thermodynamic RMSE: %.4f J/kg/K\n', validation.s_thermo_rmse);

            % Compare the two methods
            diff_methods = abs(s_pred_interp(valid_s_interp & valid_s_thermo) - ...
                              s_pred_thermo(valid_s_interp & valid_s_thermo));
            validation.s_method_diff_mean = mean(diff_methods);
            validation.s_method_diff_max = max(diff_methods);
            validation.s_method_diff_std = std(diff_methods);

            fprintf('  Difference between methods:\n');
            fprintf('    Mean: %.4f J/kg/K\n', validation.s_method_diff_mean);
            fprintf('    Max:  %.4f J/kg/K\n', validation.s_method_diff_max);
            fprintf('    Std:  %.4f J/kg/K\n', validation.s_method_diff_std);

            % Use better method as primary
            if validation.s_thermo_r2 > validation.s_interp_r2
                validation.s_rmse = validation.s_thermo_rmse;
                validation.s_mae = validation.s_thermo_mae;
                validation.s_r2 = validation.s_thermo_r2;
                validation.predictions.s = s_pred_thermo;
                fprintf('  -> Using thermodynamic integration as primary entropy validation\n');
            else
                validation.s_rmse = validation.s_interp_rmse;
                validation.s_mae = validation.s_interp_mae;
                validation.s_r2 = validation.s_interp_r2;
                validation.predictions.s = s_pred_interp;
                fprintf('  -> Using interpolation as primary entropy validation\n');
            end
        else
            fprintf('  Warning: Thermodynamic integration produced no valid values\n');
            validation.s_rmse = validation.s_interp_rmse;
            validation.s_mae = validation.s_interp_mae;
            validation.s_r2 = validation.s_interp_r2;
            validation.predictions.s = s_pred_interp;
        end

        fprintf('===========================\n');
    end

    %% Store predicted values
    validation.predictions.T = T_pred;
    validation.predictions.gamma_star = gamma_star_pred;
    validation.predictions.cv_star = cv_star_pred;
    validation.predictions.mu = mu_pred;
    validation.predictions.k = k_pred;
    validation.predictions.e = e_pred;
end


%% ========================================================================
%  INTEGRATE_ENTROPY_THERMODYNAMIC - Entropy via exact differential
%  ========================================================================
function s_integrated = INTEGRATE_ENTROPY_THERMODYNAMIC(chemistry, rho, e)
% INTEGRATE_ENTROPY_THERMODYNAMIC - Compute entropy by integrating exact differential
%
% Integrates: ds = -(gamma*-1)*cv* drho/rho + cv* de/e
%
% Uses a two-step path from a reference point (rho_ref, e_ref):
%   Step 1: (rho_ref, e_ref) -> (rho, e_ref)  [vary rho at constant e]
%   Step 2: (rho, e_ref) -> (rho, e)           [vary e at constant rho]
%
% Inputs:
%   chemistry - Fitted chemistry structure with eval functions
%   rho       - Density array (kg/m^3)
%   e         - Energy array (J/kg)
%
% Outputs:
%   s_integrated - Entropy from thermodynamic integration (J/kg/K)

    n_points = length(rho);
    s_integrated = zeros(n_points, 1);

    %% Set reference point
    [e_ref, idx_e_min] = min(e);
    rho_ref = rho(idx_e_min);
    s_ref = chemistry.s(idx_e_min);

    fprintf('    Reference point: rho_ref = %.4e kg/m^3, e_ref = %.4e J/kg, s_ref = %.2f J/kg/K\n', ...
            rho_ref, e_ref, s_ref);

    %% Create integration grid
    n_grid = 100;

    rho_grid = logspace(log10(rho_ref), log10(10), n_grid)';
    e_grid = logspace(log10(e_ref), log10(1e8), n_grid)';

    n_rho = length(rho_grid);
    n_e = length(e_grid);

    s_grid = zeros(n_rho, n_e);

    % Find reference point indices
    [~, i_ref] = min(abs(rho_grid - rho_ref));
    [~, j_ref] = min(abs(e_grid - e_ref));
    s_grid(i_ref, j_ref) = s_ref;

    fprintf('    Grid size: %d x %d = %d points\n', n_rho, n_e, n_rho * n_e);
    fprintf('    Integrating entropy field...\n');

    %% Step 1: Integrate along rows of constant e (varying rho)
    % ds = M(rho, e) drho, where M(rho, e) = -(gamma*-1)*cv*/rho
    for j = 1:n_e
        e_const = e_grid(j);

        for i = 1:n_rho
            if i == i_ref && j == j_ref
                continue;
            end

            if i < i_ref
                rho_path = flip(rho_grid(i:i_ref));
                sign_factor = -1;
            elseif i > i_ref
                rho_path = rho_grid(i_ref:i);
                sign_factor = 1;
            else
                if j == j_ref
                    s_grid(i, j) = s_ref;
                end
                continue;
            end

            % Evaluate properties along path at constant e
            gamma_star_path = chemistry.eval_gamma_star(rho_path, e_const * ones(size(rho_path)));
            cv_star_path = chemistry.eval_cv_star(rho_path, e_const * ones(size(rho_path)));

            M_path = -(gamma_star_path - 1) .* cv_star_path ./ rho_path;

            if length(rho_path) > 1
                ds_rho = sign_factor * trapz(rho_path, M_path);
            else
                ds_rho = 0;
            end

            if j == j_ref
                s_grid(i, j) = s_ref + ds_rho;
            end
        end
    end

    %% Step 2: Integrate along columns (varying e at constant rho)
    % ds = N(rho, e) de, where N(rho, e) = cv*/e
    for i = 1:n_rho
        rho_const = rho_grid(i);
        s_base = s_grid(i, j_ref);

        for j = 1:n_e
            if j == j_ref
                continue;
            end

            if j < j_ref
                e_path = flip(e_grid(j:j_ref));
                sign_factor = -1;
            else
                e_path = e_grid(j_ref:j);
                sign_factor = 1;
            end

            cv_star_path = chemistry.eval_cv_star(rho_const * ones(size(e_path)), e_path);
            N_path = cv_star_path ./ e_path;

            if length(e_path) > 1
                ds_e = sign_factor * trapz(e_path, N_path);
            else
                ds_e = 0;
            end

            s_grid(i, j) = s_base + ds_e;
        end
    end

    %% Step 3: Interpolate from grid to original data points
    fprintf('    Interpolating from grid to data points...\n');

    [E_GRID, RHO_GRID] = ndgrid(e_grid, rho_grid);
    s_interpolant = griddedInterpolant(E_GRID, RHO_GRID, s_grid', 'linear', 'linear');
    s_integrated = s_interpolant(e, rho);

    % Report statistics
    valid_integrated = ~isnan(s_integrated) & ~isinf(s_integrated);
    if sum(valid_integrated) > 0
        fprintf('    Valid integrated points: %d / %d (%.1f%%)\n', ...
                sum(valid_integrated), n_points, 100 * sum(valid_integrated) / n_points);
        fprintf('    Integrated entropy range: %.2f to %.2f J/kg/K\n', ...
                min(s_integrated(valid_integrated)), max(s_integrated(valid_integrated)));
    else
        fprintf('    Warning: No valid integrated values produced\n');
    end

    %% Compare with original interpolant
    fprintf('\n=== Comparing with Original Interpolant ===\n');

    s_original = chemistry.eval_s(rho, e);
    s_diff = s_integrated - s_original;
    s_rel_diff = 100 * s_diff ./ abs(s_original);

    valid = valid_integrated & isfinite(s_original);
    fprintf('Valid comparison points: %d / %d (%.1f%%)\n', ...
            sum(valid), n_points, 100 * sum(valid) / n_points);

    CREATE_COMPARISON_PLOTS(rho, e, s_original, s_integrated, s_diff, s_rel_diff, valid);
end


%% ========================================================================
%  CREATE_COMPARISON_PLOTS - 2D contour comparison of entropy fields
%  ========================================================================
function CREATE_COMPARISON_PLOTS(rho, e, s_original, s_integrated, s_diff, s_rel_diff, valid)
% CREATE_COMPARISON_PLOTS - Create 2D contour plots comparing entropy fields
%
% Generates side-by-side contour plots of the original interpolant,
% thermodynamic integration result, and their difference.
%
% Inputs:
%   rho, e         - Density and energy data points
%   s_original     - Entropy from direct interpolation
%   s_integrated   - Entropy from thermodynamic integration
%   s_diff         - Difference (integrated - original)
%   s_rel_diff     - Relative difference in percent
%   valid          - Logical mask for valid comparison points

    if length(unique(rho)) > 5 && length(unique(e)) > 5
        try
            rho_unique = unique(rho);
            e_unique = unique(e);

            % Check if data is approximately on a grid
            if length(rho_unique) * length(e_unique) >= 0.8 * length(rho)
                [RHO_GRID, E_GRID] = meshgrid(rho_unique, e_unique);

                F_orig = scatteredInterpolant(rho(valid), e(valid), s_original(valid), 'linear', 'none');
                F_integ = scatteredInterpolant(rho(valid), e(valid), s_integrated(valid), 'linear', 'none');
                F_diff = scatteredInterpolant(rho(valid), e(valid), s_diff(valid), 'linear', 'none');

                S_orig_grid = F_orig(RHO_GRID, E_GRID);
                S_integ_grid = F_integ(RHO_GRID, E_GRID);
                S_diff_grid = F_diff(RHO_GRID, E_GRID);

                figure('Position', [150, 150, 1400, 400], 'Name', 'Entropy Field Comparison');

                subplot(1, 3, 1);
                contourf(RHO_GRID, E_GRID, S_orig_grid, 20, 'LineColor', 'none');
                colorbar;
                xlabel('Density \rho (kg/m^3)');
                ylabel('Energy e (J/kg)');
                title('Original Interpolant');
                set(gca, 'FontSize', 10);

                subplot(1, 3, 2);
                contourf(RHO_GRID, E_GRID, S_integ_grid, 20, 'LineColor', 'none');
                colorbar;
                xlabel('Density \rho (kg/m^3)');
                ylabel('Energy e (J/kg)');
                title('Thermodynamic Integration');
                set(gca, 'FontSize', 10);

                subplot(1, 3, 3);
                contourf(RHO_GRID, E_GRID, S_diff_grid, 20, 'LineColor', 'none');
                colorbar;
                xlabel('Density \rho (kg/m^3)');
                ylabel('Energy e (J/kg)');
                title('Difference (Integrated - Original)');
                set(gca, 'FontSize', 10);
                caxis([-max(abs(S_diff_grid(:))), max(abs(S_diff_grid(:)))]);
                colormap(gca, 'redblue');

                sgtitle('2D Entropy Field Comparison', 'FontSize', 13, 'FontWeight', 'bold');
            end
        catch
            fprintf('    Note: Could not create contour plots (data not gridded)\n');
        end
    end
end


%% ========================================================================
%  TEST_FIT_PERFORMANCE - Benchmark interpolation speed
%  ========================================================================
function TEST_FIT_PERFORMANCE(chemistry)
% TEST_FIT_PERFORMANCE - Test evaluation speed of the fitted models
%
% Benchmarks vectorized and single-point evaluation times for all
% properties including the inverse E(T, rho) function.
%
% Inputs:
%   chemistry - Fitted chemistry structure
%
% Notes:
%   Tests 10000 random points for vectorized evaluation and 1000
%   repeated single-point evaluations.

    if ~isfield(chemistry, 'fit_T')
        error('Chemistry data must be fitted first. Run FIT_CHEMISTRY(chemistry)');
    end

    fprintf('\nTesting fit performance with scaled griddedInterpolant...\n');

    %% Generate test points
    n_test = 10000;
    rho_range = chemistry.fit_info.rho_range;
    e_range = chemistry.fit_info.e_range;
    T_range = chemistry.fit_info.T_range;

    rho_test = rho_range(1) + (rho_range(2) - rho_range(1)) * rand(n_test, 1);
    e_test = e_range(1) + (e_range(2) - e_range(1)) * rand(n_test, 1);
    T_test = T_range(1) + (T_range(2) - T_range(1)) * rand(n_test, 1);

    %% Time vectorized evaluations
    fprintf('Timing %d evaluations:\n', n_test);

    tic;
    T_eval = chemistry.eval_T(rho_test, e_test);
    t_T = toc;
    fprintf('  T evaluation: %.4f seconds (%.2e sec/point)\n', t_T, t_T / n_test);

    tic;
    gamma_star_test = chemistry.eval_gamma_star(rho_test, e_test);
    t_gamma_star = toc;
    fprintf('  gamma_star evaluation: %.4f seconds (%.2e sec/point)\n', t_gamma_star, t_gamma_star / n_test);

    if chemistry.fit_info.has_sound_speed_fit
        tic;
        a_test = chemistry.eval_a(rho_test, e_test);
        t_a = toc;
        fprintf('  Sound speed evaluation: %.4f seconds (%.2e sec/point)\n', t_a, t_a / n_test);
    end

    tic;
    all_props = chemistry.eval_all(rho_test, e_test);
    t_all = toc;
    fprintf('  All basic properties: %.4f seconds (%.2e sec/point)\n', t_all, t_all / n_test);

    tic;
    e_inverse = chemistry.eval_e(T_test, rho_test);
    t_e_inverse = toc;
    fprintf('  E(T,rho) inverse evaluation: %.4f seconds (%.2e sec/point)\n', t_e_inverse, t_e_inverse / n_test);

    if chemistry.has_molar_data
        tic;
        all_species = chemistry.eval_all_species(rho_test, e_test);
        t_species = toc;
        fprintf('  All species: %.4f seconds (%.2e sec/point)\n', t_species, t_species / n_test);
    end

    %% Time single-point evaluations
    rho_single = mean(rho_range);
    e_single = mean(e_range);
    T_single = mean(T_range);

    tic;
    for i = 1:1000
        T_eval_single = chemistry.eval_T(rho_single, e_single);
    end
    t_single = toc;
    fprintf('  Single point (1000 calls): %.4f seconds (%.2e sec/call)\n', t_single, t_single / 1000);

    tic;
    for i = 1:1000
        e_inverse_single = chemistry.eval_e(T_single, rho_single);
    end
    t_inverse_single = toc;
    fprintf('  Single inverse E(T,rho) (1000 calls): %.4f seconds (%.2e sec/call)\n', t_inverse_single, t_inverse_single / 1000);

    %% Print summary
    fprintf('\ngriddedInterpolant with [-1,1] scaling provides fast, smooth interpolation.\n');
    fprintf('Interpolant type: %s\n', chemistry.fit_info.interpolant_type);
    fprintf('Grid resolution: %dx%d points\n', chemistry.fit_info.grid_size(1), chemistry.fit_info.grid_size(2));
    fprintf('Griddata method: %s\n', chemistry.fit_info.griddata_method);
    fprintf('Interpolation method: %s\n', chemistry.fit_info.interp_method);
    fprintf('Inverse E(T,rho) function ready for Rankine-Hugoniot solver.\n');
    if chemistry.fit_info.has_sound_speed_fit
        fprintf('Sound speed interpolation included.\n');
    end
end


%% ========================================================================
%  PLOT_FIT_EVALUATION - Validation scatter plots
%  ========================================================================
function PLOT_FIT_EVALUATION(chemistry)
% PLOT_FIT_EVALUATION - Plot predicted vs. true values for all fitted properties
%
% Creates a multi-panel figure with scatter plots of predicted vs. true
% values for T, gamma*, cv*, mu, k, E(T,rho), and optionally sound speed
% and entropy. Each panel displays the R^2 value.
%
% Inputs:
%   chemistry - Fitted chemistry structure with fit_validation field

    if ~isfield(chemistry, 'fit_validation')
        error('Chemistry data must be fitted first. Run FIT_CHEMISTRY(chemistry)');
    end

    val = chemistry.fit_validation;
    has_a = isfield(val, 'a_rmse');
    has_s = isfield(val, 's_rmse');

    %% Determine subplot layout
    if has_a && has_s
        subplot_layout = [3, 3];
    elseif has_a || has_s
        subplot_layout = [3, 3];
    else
        subplot_layout = [2, 3];
    end

    %% Create validation plots
    figure('Name', 'Chemistry Fit Validation (scaled griddedInterpolant)', 'Position', [100, 100, 1200, 900]);

    % Temperature
    subplot(subplot_layout(1), subplot_layout(2), 1);
    plot(chemistry.T, val.predictions.T, 'b.', 'MarkerSize', 3);
    hold on;
    plot([min(chemistry.T), max(chemistry.T)], [min(chemistry.T), max(chemistry.T)], 'r-', 'LineWidth', 2);
    xlabel('True T (K)');
    ylabel('Predicted T (K)');
    title(['Temperature (R$^2$ = ' sprintf('%.4f', val.T_r2) ')'], 'Interpreter', 'latex');
    grid on;
    axis equal;

    % gamma_star
    subplot(subplot_layout(1), subplot_layout(2), 2);
    plot(chemistry.gamma_star, val.predictions.gamma_star, 'g.', 'MarkerSize', 3);
    hold on;
    plot([min(chemistry.gamma_star), max(chemistry.gamma_star)], [min(chemistry.gamma_star), max(chemistry.gamma_star)], 'r-', 'LineWidth', 2);
    xlabel('True $\gamma^*$', 'Interpreter', 'latex');
    ylabel('Predicted $\gamma^*$', 'Interpreter', 'latex');
    title(['Heat Capacity Ratio (R$^2$ = ' sprintf('%.4f', val.gamma_star_r2) ')'], 'Interpreter', 'latex');
    grid on;
    axis equal;

    % Specific Heat
    subplot(subplot_layout(1), subplot_layout(2), 3);
    plot(chemistry.cv_star, val.predictions.cv_star, 'm.', 'MarkerSize', 3);
    hold on;
    plot([min(chemistry.cv_star), max(chemistry.cv_star)], [min(chemistry.cv_star), max(chemistry.cv_star)], 'r-', 'LineWidth', 2);
    xlabel('True $c_v^*$ (J/kg/K)', 'Interpreter', 'latex');
    ylabel('Predicted $c_v^*$ (J/kg/K)', 'Interpreter', 'latex');
    title(['Specific Heat (R$^2$ = ' sprintf('%.4f', val.cv_star_r2) ')'], 'Interpreter', 'latex');
    grid on;
    axis equal;

    % Viscosity
    subplot(subplot_layout(1), subplot_layout(2), 4);
    plot(chemistry.mu, val.predictions.mu, 'c.', 'MarkerSize', 3);
    hold on;
    plot([min(chemistry.mu), max(chemistry.mu)], [min(chemistry.mu), max(chemistry.mu)], 'r-', 'LineWidth', 2);
    xlabel('True $\mu$ (Pa$\cdot$s)', 'Interpreter', 'latex');
    ylabel('Predicted $\mu$ (Pa$\cdot$s)', 'Interpreter', 'latex');
    title(['Viscosity (R$^2$ = ' sprintf('%.4f', val.mu_r2) ')'], 'Interpreter', 'latex');
    grid on;
    axis equal;

    % Thermal Conductivity
    subplot(subplot_layout(1), subplot_layout(2), 5);
    plot(chemistry.k, val.predictions.k, 'y.', 'MarkerSize', 3);
    hold on;
    plot([min(chemistry.k), max(chemistry.k)], [min(chemistry.k), max(chemistry.k)], 'r-', 'LineWidth', 2);
    xlabel('True k (W/m/K)');
    ylabel('Predicted k (W/m/K)');
    title(['Thermal Conductivity (R$^2$ = ' sprintf('%.4f', val.k_r2) ')'], 'Interpreter', 'latex');
    grid on;
    axis equal;

    % Energy inverse fit
    subplot(subplot_layout(1), subplot_layout(2), 6);
    plot(chemistry.e, val.predictions.e, 'k.', 'MarkerSize', 3);
    hold on;
    plot([min(chemistry.e), max(chemistry.e)], [min(chemistry.e), max(chemistry.e)], 'r-', 'LineWidth', 2);
    xlabel('True E (J/kg)');
    ylabel('Predicted E (J/kg)');
    title(['Energy $E(T,\rho)$ (R$^2$ = ' sprintf('%.4f', val.e_r2) ')'], 'Interpreter', 'latex');
    grid on;
    axis equal;

    % Sound speed (if available)
    plot_idx = 7;
    if has_a
        subplot(subplot_layout(1), subplot_layout(2), plot_idx);
        plot(chemistry.a, val.predictions.a, 'b.', 'MarkerSize', 3);
        hold on;
        plot([min(chemistry.a), max(chemistry.a)], [min(chemistry.a), max(chemistry.a)], 'r-', 'LineWidth', 2);
        xlabel('True a (m/s)');
        ylabel('Predicted a (m/s)');
        title(['Sound Speed (R$^2$ = ' sprintf('%.4f', val.a_r2) ')'], 'Interpreter', 'latex');
        grid on;
        axis equal;
        plot_idx = plot_idx + 1;
    end

    % Entropy (if available)
    if has_s
        subplot(subplot_layout(1), subplot_layout(2), plot_idx);
        plot(chemistry.s, val.predictions.s, 'r.', 'MarkerSize', 3);
        hold on;
        plot([min(chemistry.s), max(chemistry.s)], [min(chemistry.s), max(chemistry.s)], 'r-', 'LineWidth', 2);
        xlabel('True s (J/kg/K)');
        ylabel('Predicted s (J/kg/K)');
        title(['Entropy (R$^2$ = ' sprintf('%.4f', val.s_r2) ')'], 'Interpreter', 'latex');
        grid on;
        axis equal;
    end
end


%% ========================================================================
%  PLOT_3D_CHEMISTRY_SURFACES - 3D surface visualization
%  ========================================================================
function PLOT_3D_CHEMISTRY_SURFACES(chemistry)
% PLOT_3D_CHEMISTRY_SURFACES - Create 3D surface plots of all chemistry properties
%
% Generates surface plots showing gamma*, cv*, mu, k, a, s as functions
% of (T, rho) and (E, rho), with original scattered data points overlaid.
% Also creates standalone temperature T(E, rho) and energy E(T, rho) surfaces.
%
% Inputs:
%   chemistry - Fitted chemistry structure from FIT_CHEMISTRY

    if ~isfield(chemistry, 'fit_T')
        error('Chemistry data must be fitted first. Run FIT_CHEMISTRY(chemistry)');
    end

    has_a = chemistry.fit_info.has_sound_speed_fit;
    has_s = isfield(chemistry.fit_info, 'has_entropy_fit') && chemistry.fit_info.has_entropy_fit;

    fprintf('Creating 3D surface plots with original data points (scaled griddedInterpolant)...\n');

    %% Get data ranges and create evaluation grids
    rho_range = chemistry.fit_info.rho_range;
    e_range = chemistry.fit_info.e_range;
    T_range = chemistry.fit_info.T_range;

    n_points = 200;

    % (E, rho) grid - use ndgrid for griddedInterpolant compatibility
    e_grid = linspace(e_range(1), e_range(2), n_points);
    rho_grid = linspace(rho_range(1), rho_range(2), n_points);
    [e_mesh, rho_mesh] = ndgrid(e_grid, rho_grid);

    %% Evaluate properties on (E, rho) grid
    gamma_star_e_rho = reshape(chemistry.eval_gamma_star(rho_mesh, e_mesh), size(rho_mesh));
    cv_star_e_rho = reshape(chemistry.eval_cv_star(rho_mesh, e_mesh), size(rho_mesh));
    mu_e_rho = reshape(chemistry.eval_mu(rho_mesh, e_mesh), size(rho_mesh));
    k_e_rho = reshape(chemistry.eval_k(rho_mesh, e_mesh), size(rho_mesh));
    T_e_rho = reshape(chemistry.eval_T(rho_mesh, e_mesh), size(rho_mesh));
    if has_a
        a_e_rho = reshape(chemistry.eval_a(rho_mesh, e_mesh), size(rho_mesh));
    end
    if has_s
        s_e_rho = reshape(chemistry.eval_s(rho_mesh, e_mesh), size(rho_mesh));
    end

    %% Evaluate properties on (T, rho) grid
    fprintf('  Using built-in E(T,rho) interpolant for (T,rho) plots...\n');

    T_grid = linspace(T_range(1), T_range(2), n_points);
    [T_mesh, rho_mesh_T] = ndgrid(T_grid, rho_grid);

    e_from_T_rho = chemistry.eval_e(T_mesh, rho_mesh_T);

    gamma_star_T_rho = reshape(chemistry.eval_gamma_star(rho_mesh_T, e_from_T_rho), size(rho_mesh_T));
    cv_star_T_rho = reshape(chemistry.eval_cv_star(rho_mesh_T, e_from_T_rho), size(rho_mesh_T));
    mu_T_rho = reshape(chemistry.eval_mu(rho_mesh_T, e_from_T_rho), size(rho_mesh_T));
    k_T_rho = reshape(chemistry.eval_k(rho_mesh_T, e_from_T_rho), size(rho_mesh_T));
    if has_a
        a_T_rho = reshape(chemistry.eval_a(rho_mesh_T, e_from_T_rho), size(rho_mesh_T));
    end
    if has_s
        s_T_rho = reshape(chemistry.eval_s(rho_mesh_T, e_from_T_rho), size(rho_mesh_T));
    end

    %% Prepare original data points for overlay
    if has_a && has_s
        valid_idx = ~(isnan(chemistry.T) | isnan(chemistry.rho) | isnan(chemistry.e) | ...
                      isnan(chemistry.gamma_star) | isnan(chemistry.cv_star) | ...
                      isnan(chemistry.mu) | isnan(chemistry.k) | isnan(chemistry.a) | isnan(chemistry.s));
    elseif has_a
        valid_idx = ~(isnan(chemistry.T) | isnan(chemistry.rho) | isnan(chemistry.e) | ...
                      isnan(chemistry.gamma_star) | isnan(chemistry.cv_star) | ...
                      isnan(chemistry.mu) | isnan(chemistry.k) | isnan(chemistry.a));
    elseif has_s
        valid_idx = ~(isnan(chemistry.T) | isnan(chemistry.rho) | isnan(chemistry.e) | ...
                      isnan(chemistry.gamma_star) | isnan(chemistry.cv_star) | ...
                      isnan(chemistry.mu) | isnan(chemistry.k) | isnan(chemistry.s));
    else
        valid_idx = ~(isnan(chemistry.T) | isnan(chemistry.rho) | isnan(chemistry.e) | ...
                      isnan(chemistry.gamma_star) | isnan(chemistry.cv_star) | ...
                      isnan(chemistry.mu) | isnan(chemistry.k));
    end

    T_data = chemistry.T(valid_idx);
    rho_data = chemistry.rho(valid_idx);
    e_data = chemistry.e(valid_idx);
    gamma_star_data = chemistry.gamma_star(valid_idx);
    cv_star_data = chemistry.cv_star(valid_idx);
    mu_data = chemistry.mu(valid_idx);
    k_data = chemistry.k(valid_idx);
    if has_a
        a_data = chemistry.a(valid_idx);
    end
    if has_s
        s_data = chemistry.s(valid_idx);
    end

    fprintf('  Plotting %d valid data points\n', sum(valid_idx));

    % Determine subplot layout
    if has_a && has_s
        subplot_rows = 3;
    elseif has_a || has_s
        subplot_rows = 3;
    else
        subplot_rows = 2;
    end

    %% Figure 1: Properties vs (T, rho)
    figure('Name', 'Chemistry Properties vs (T, rho) - Scaled scatteredInterpolant', 'Position', [50, 50, 1000, 800]);

    % gamma_star vs (T, rho)
    subplot(subplot_rows, 2, 1);
    surf(T_mesh / 1000, rho_mesh_T, gamma_star_T_rho, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    scatter3(T_data / 1000, rho_data, gamma_star_data, 20, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    xlabel('T (kK)');
    ylabel('$\rho$ (kg/m$^3$)', 'Interpreter', 'latex');
    zlabel('$\gamma^*$', 'Interpreter', 'latex');
    title('Heat Capacity Ratio $\gamma^*(T,\rho)$', 'Interpreter', 'latex');
    colorbar; view(45, 30); grid on;
    legend('Interpolated Surface', 'Original Data', 'Location', 'best');

    % cv_star vs (T, rho)
    subplot(subplot_rows, 2, 2);
    surf(T_mesh / 1000, rho_mesh_T, cv_star_T_rho, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    scatter3(T_data / 1000, rho_data, cv_star_data, 20, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    xlabel('T (kK)');
    ylabel('$\rho$ (kg/m$^3$)', 'Interpreter', 'latex');
    zlabel('$c_v^*$ (J/kg/K)', 'Interpreter', 'latex');
    title('Specific Heat $c_v^*(T,\rho)$', 'Interpreter', 'latex');
    colorbar; view(45, 30); grid on;
    legend('Interpolated Surface', 'Original Data', 'Location', 'best');

    % mu vs (T, rho)
    subplot(subplot_rows, 2, 3);
    surf(T_mesh / 1000, rho_mesh_T, mu_T_rho * 1e6, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    scatter3(T_data / 1000, rho_data, mu_data * 1e6, 20, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    xlabel('T (kK)');
    ylabel('$\rho$ (kg/m$^3$)', 'Interpreter', 'latex');
    zlabel('$\mu$ ($\mu$Pa$\cdot$s)', 'Interpreter', 'latex');
    title('Viscosity $\mu(T,\rho)$', 'Interpreter', 'latex');
    colorbar; view(45, 30); grid on;
    legend('Interpolated Surface', 'Original Data', 'Location', 'best');

    % k vs (T, rho)
    subplot(subplot_rows, 2, 4);
    surf(T_mesh / 1000, rho_mesh_T, k_T_rho * 1000, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    scatter3(T_data / 1000, rho_data, k_data * 1000, 20, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    xlabel('T (kK)');
    ylabel('$\rho$ (kg/m$^3$)', 'Interpreter', 'latex');
    zlabel('k (mW/m/K)');
    title('Thermal Conductivity $k(T,\rho)$', 'Interpreter', 'latex');
    colorbar; view(45, 30); grid on;
    legend('Interpolated Surface', 'Original Data', 'Location', 'best');

    % Sound speed vs (T, rho)
    if has_a
        subplot(subplot_rows, 2, 5);
        surf(T_mesh / 1000, rho_mesh_T, a_T_rho, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        hold on;
        scatter3(T_data / 1000, rho_data, a_data, 20, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        xlabel('T (kK)');
        ylabel('$\rho$ (kg/m$^3$)', 'Interpreter', 'latex');
        zlabel('a (m/s)');
        title('Sound Speed $a(T,\rho)$', 'Interpreter', 'latex');
        colorbar; view(45, 30); grid on;
        legend('Interpolated Surface', 'Original Data', 'Location', 'best');
    end

    % Entropy vs (T, rho)
    if has_s
        if has_a
            subplot_idx = 6;
        else
            subplot_idx = 5;
        end
        subplot(subplot_rows, 2, subplot_idx);
        surf(T_mesh / 1000, rho_mesh_T, s_T_rho, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        hold on;
        scatter3(T_data / 1000, rho_data, s_data, 20, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        xlabel('T (kK)');
        ylabel('$\rho$ (kg/m$^3$)', 'Interpreter', 'latex');
        zlabel('s (J/kg/K)');
        title('Entropy $s(T,\rho)$', 'Interpreter', 'latex');
        colorbar; view(45, 30); grid on;
        legend('Interpolated Surface', 'Original Data', 'Location', 'best');
    end

    %% Figure 2: Properties vs (E, rho)
    figure('Name', 'Chemistry Properties vs (E, rho) - Scaled scatteredInterpolant', 'Position', [50, 50, 1000, 800]);

    % gamma_star vs (E, rho)
    subplot(subplot_rows, 2, 1);
    surf(e_mesh / 1e6, rho_mesh, gamma_star_e_rho, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    scatter3(e_data / 1e6, rho_data, gamma_star_data, 20, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    xlabel('E (MJ/kg)');
    ylabel('$\rho$ (kg/m$^3$)', 'Interpreter', 'latex');
    zlabel('$\gamma^*$', 'Interpreter', 'latex');
    title('Heat Capacity Ratio $\gamma^*(E,\rho)$', 'Interpreter', 'latex');
    colorbar; view(45, 30); grid on;
    legend('Interpolated Surface', 'Original Data', 'Location', 'best');

    % cv_star vs (E, rho)
    subplot(subplot_rows, 2, 2);
    surf(e_mesh / 1e6, rho_mesh, cv_star_e_rho, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    scatter3(e_data / 1e6, rho_data, cv_star_data, 20, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    xlabel('E (MJ/kg)');
    ylabel('$\rho$ (kg/m$^3$)', 'Interpreter', 'latex');
    zlabel('$c_v^*$ (J/kg/K)', 'Interpreter', 'latex');
    title('Specific Heat $c_v^*(E,\rho)$', 'Interpreter', 'latex');
    colorbar; view(45, 30); grid on;
    legend('Interpolated Surface', 'Original Data', 'Location', 'best');

    % mu vs (E, rho)
    subplot(subplot_rows, 2, 3);
    surf(e_mesh / 1e6, rho_mesh, mu_e_rho * 1e6, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    scatter3(e_data / 1e6, rho_data, mu_data * 1e6, 20, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    xlabel('E (MJ/kg)');
    ylabel('$\rho$ (kg/m$^3$)', 'Interpreter', 'latex');
    zlabel('$\mu$ ($\mu$Pa$\cdot$s)', 'Interpreter', 'latex');
    title('Viscosity $\mu(E,\rho)$', 'Interpreter', 'latex');
    colorbar; view(45, 30); grid on;
    legend('Interpolated Surface', 'Original Data', 'Location', 'best');

    % k vs (E, rho)
    subplot(subplot_rows, 2, 4);
    surf(e_mesh / 1e6, rho_mesh, k_e_rho * 1000, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    scatter3(e_data / 1e6, rho_data, k_data * 1000, 20, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    xlabel('E (MJ/kg)');
    ylabel('$\rho$ (kg/m$^3$)', 'Interpreter', 'latex');
    zlabel('k (mW/m/K)');
    title('Thermal Conductivity $k(E,\rho)$', 'Interpreter', 'latex');
    colorbar; view(45, 30); grid on;
    legend('Interpolated Surface', 'Original Data', 'Location', 'best');

    % Sound speed vs (E, rho)
    if has_a
        subplot(subplot_rows, 2, 5);
        surf(e_mesh / 1e6, rho_mesh, a_e_rho, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        hold on;
        scatter3(e_data / 1e6, rho_data, a_data, 20, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        xlabel('E (MJ/kg)');
        ylabel('$\rho$ (kg/m$^3$)', 'Interpreter', 'latex');
        zlabel('a (m/s)');
        title('Sound Speed $a(E,\rho)$', 'Interpreter', 'latex');
        colorbar; view(45, 30); grid on;
        legend('Interpolated Surface', 'Original Data', 'Location', 'best');
    end

    % Entropy vs (E, rho)
    if has_s
        if has_a
            subplot_idx = 6;
        else
            subplot_idx = 5;
        end
        subplot(subplot_rows, 2, subplot_idx);
        surf(e_mesh / 1e6, rho_mesh, s_e_rho, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        hold on;
        scatter3(e_data / 1e6, rho_data, s_data, 20, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        xlabel('E (MJ/kg)');
        ylabel('$\rho$ (kg/m$^3$)', 'Interpreter', 'latex');
        zlabel('s (J/kg/K)');
        title('Entropy $s(E,\rho)$', 'Interpreter', 'latex');
        colorbar; view(45, 30); grid on;
        legend('Interpolated Surface', 'Original Data', 'Location', 'best');
    end

    %% Figure 3: Temperature surface T(E, rho)
    figure('Name', 'Temperature Surface T(E,rho) - Scaled scatteredInterpolant', 'Position', [150, 150, 800, 600]);
    surf(e_mesh / 1e6, rho_mesh, T_e_rho / 1000, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    scatter3(e_data / 1e6, rho_data, T_data / 1000, 20, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    xlabel('E (MJ/kg)');
    ylabel('$\rho$ (kg/m$^3$)', 'Interpreter', 'latex');
    zlabel('T (kK)');
    title('Temperature $T(E,\rho)$', 'Interpreter', 'latex');
    colorbar; view(45, 30); grid on;
    legend('Interpolated Surface', 'Original Data', 'Location', 'best');

    %% Figure 4: Energy surface E(T, rho)
    figure('Name', 'Energy Surface E(T,rho) - Inverse Fit', 'Position', [200, 200, 800, 600]);
    surf(T_mesh / 1000, rho_mesh_T, e_from_T_rho / 1e6, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    scatter3(T_data / 1000, rho_data, e_data / 1e6, 20, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    xlabel('T (kK)');
    ylabel('$\rho$ (kg/m$^3$)', 'Interpreter', 'latex');
    zlabel('E (MJ/kg)');
    title('Energy $E(T,\rho)$ - Inverse Interpolant', 'Interpreter', 'latex');
    colorbar; view(45, 30); grid on;
    legend('Interpolated Surface', 'Original Data', 'Location', 'best');

    %% Print summary
    fprintf('3D surface plots created successfully with scaled scatteredInterpolant.\n');
    fprintf('Inverse E(T,rho) surface plot created to validate inverse interpolant.\n');
    if has_a
        fprintf('Sound speed surface plots included.\n');
    end
    if has_s
        fprintf('Entropy surface plots included.\n');
    end

    fprintf('\nData coverage statistics:\n');
    fprintf('  Total valid data points: %d\n', sum(valid_idx));
    fprintf('  Temperature range: %.1f - %.1f K\n', T_range);
    fprintf('  Energy range: %.2e - %.2e J/kg\n', e_range);
    fprintf('  Density range: %.2e - %.2e kg/m^3\n', rho_range);
    if has_a
        fprintf('  Sound speed range: %.1f - %.1f m/s\n', chemistry.fit_info.a_range);
    end
    if has_s
        fprintf('  Entropy range: %.2e - %.2e J/kg/K\n', chemistry.fit_info.s_range);
    end
end
