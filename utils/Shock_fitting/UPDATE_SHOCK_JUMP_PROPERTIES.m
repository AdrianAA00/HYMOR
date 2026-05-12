function s = UPDATE_SHOCK_JUMP_PROPERTIES(s, chemistry)
% UPDATE_SHOCK_JUMP_PROPERTIES  Update post-shock flow properties via Rankine-Hugoniot relations.
%
%   s = UPDATE_SHOCK_JUMP_PROPERTIES(s, chemistry) solves the Rankine-Hugoniot
%   jump conditions across the fitted shock using a Newton-Raphson iteration
%   on the Riemann invariant formulation. The routine interpolates the
%   downstream flow field to the shock location, iterates for the post-shock
%   pressure, and back-computes all conservative variables and shock speeds.
%
%   Inputs:
%       s         - Solution structure containing the flow field, grid, shock
%                   geometry, and solver parameters.
%       chemistry - Chemistry model structure providing thermodynamic property
%                   evaluators (gamma_star, cv_star) for equilibrium or frozen
%                   chemistry.
%
%   Outputs:
%       s - Updated s structure with:
%             .shock.properties  (rho, rho_u, rho_v, rho_E, p, gamma_star, cv_star)
%             .shock.speed_x, .shock.speed_y
%             .shock.relative_increase_velocity
%           and interpolated shocked-cell values in the flow field arrays.
%
%   Notes:
%       - Dispatches to UPDATE_SHOCK_JUMP_PROPERTIES_PRESSURE_3 (active method).
%       - Supports both Lagrangian and Eulerian shock-speed formulations.
%       - Handles calorically perfect gas, frozen chemistry, and equilibrium
%         chemistry via the chemistry argument.
%       - Ghost-cell extrapolation is performed via EXTRAPOLATE_CELLS_SHOCK.
%
% Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

    s = UPDATE_SHOCK_JUMP_PROPERTIES_PRESSURE_3(s, chemistry);
end

%% ========================================================================
%  Riemann-Invariant Rankine-Hugoniot Solver
%  ========================================================================
function s = UPDATE_SHOCK_JUMP_PROPERTIES_RIEMANN(s, chemistry)
% UPDATE_SHOCK_JUMP_PROPERTIES_RIEMANN  Solve R-H with Newton-Raphson on
%   the Riemann invariant (C- characteristic) to obtain post-shock pressure.

    %% Upstream conditions
    [p_infty, u_inf, v_inf, rho_inf, e_inf] = COMPUTE_UPSTREAM_CONDITIONS(s);
    s = UPDATE_FIELD_UPSTREAM(s);

    %% Index computation for shocked cells
    valid_ix = (1:s.mesh.Nchi)';
    sc_idy = s.shock.cell_indices(valid_ix, 1);

    lin_idx_0 = sub2ind(size(s.mesh.x), valid_ix, sc_idy);
    lin_idx_1 = sub2ind(size(s.mesh.x), valid_ix, sc_idy - 1);
    lin_idx_2 = sub2ind(size(s.mesh.x), valid_ix, sc_idy - 2);

    x_0 = s.mesh.x(lin_idx_0);
    y_0 = s.mesh.y(lin_idx_0);
    x_1 = s.mesh.x(lin_idx_1);
    y_1 = s.mesh.y(lin_idx_1);
    x_2 = s.mesh.x(lin_idx_2);
    y_2 = s.mesh.y(lin_idx_2);
    x_s = s.shock.points_x(valid_ix);
    y_s = s.shock.points_y(valid_ix);

    %% Distance computation for interpolation
    dist01 = sqrt((x_0 - x_1).^2 + (y_0 - y_1).^2);
    dist12 = sqrt((x_1 - x_2).^2 + (y_1 - y_2).^2);
    dist1s = ((x_s - x_1).*(x_1 - x_2) + (y_s - y_1).*(y_1 - y_2)) ./ dist12;
    dist0s = ((x_s - x_0).*(x_0 - x_1) + (y_s - y_0).*(y_0 - y_1)) ./ dist01;

    %% Flow field gradient computation for linear interpolation
    p_1 = s.var.p(sub2ind(size(s.var.p), valid_ix + 1, sc_idy + 1 - 1));
    p_2 = s.var.p(sub2ind(size(s.var.p), valid_ix + 1, sc_idy + 1 - 2));
    slope_p = (p_1 - p_2) ./ dist12;

    rho_1 = s.var.rho(sub2ind(size(s.var.rho), valid_ix + 1, sc_idy + 1 - 1));
    rho_2 = s.var.rho(sub2ind(size(s.var.rho), valid_ix + 1, sc_idy + 1 - 2));
    slope_rho = (rho_1 - rho_2) ./ dist12;

    rho_u_1 = s.var.rho_u(sub2ind(size(s.var.rho_u), valid_ix + 1, sc_idy + 1 - 1));
    rho_u_2 = s.var.rho_u(sub2ind(size(s.var.rho_u), valid_ix + 1, sc_idy + 1 - 2));
    slope_rho_u = (rho_u_1 - rho_u_2) ./ dist12;

    rho_v_1 = s.var.rho_v(sub2ind(size(s.var.rho_v), valid_ix + 1, sc_idy + 1 - 1));
    rho_v_2 = s.var.rho_v(sub2ind(size(s.var.rho_v), valid_ix + 1, sc_idy + 1 - 2));
    slope_rho_v = (rho_v_1 - rho_v_2) ./ dist12;

    rho_E_1 = s.var.rho_E(sub2ind(size(s.var.rho_E), valid_ix + 1, sc_idy + 1 - 1));
    rho_E_2 = s.var.rho_E(sub2ind(size(s.var.rho_E), valid_ix + 1, sc_idy + 1 - 2));
    slope_rho_E = (rho_E_1 - rho_E_2) ./ dist12;

    a_1 = s.var.a(sub2ind(size(s.var.a), valid_ix + 1, sc_idy + 1 - 1));
    a_2 = s.var.a(sub2ind(size(s.var.a), valid_ix + 1, sc_idy + 1 - 2));
    slope_a = (a_1 - a_2) ./ dist12;

    %% Internal energy gradient
    rho_velocity_squared = (s.var.rho_u.^2 + s.var.rho_v.^2) ./ s.var.rho;
    e = (s.var.rho_E - rho_velocity_squared / 2) ./ s.var.rho;
    e_1 = e(sub2ind(size(e), valid_ix + 1, sc_idy + 1 - 1));
    e_2 = e(sub2ind(size(e), valid_ix + 1, sc_idy + 1 - 2));
    slope_e = (e_1 - e_2) ./ dist12;

    %% Non-equilibrium chemistry gradients (gamma_star, cv_star)
    if ~s.chemistry.chemical_equilibrium && s.chemistry.is_chemistry_enabled
        gamma_star_1 = s.var.gamma_star(sub2ind(size(s.var.gamma_star), valid_ix + 1, sc_idy + 1 - 1));
        gamma_star_2 = s.var.gamma_star(sub2ind(size(s.var.gamma_star), valid_ix + 1, sc_idy + 1 - 2));
        slope_gamma_star = (gamma_star_1 - gamma_star_2) ./ dist12;
        cv_star_1 = s.var.cv_star(sub2ind(size(s.var.cv_star), valid_ix + 1, sc_idy + 1 - 1));
        cv_star_2 = s.var.cv_star(sub2ind(size(s.var.cv_star), valid_ix + 1, sc_idy + 1 - 2));
        slope_cv_star = (cv_star_1 - cv_star_2) ./ dist12;
    end

    %% Initial guess from flow field interpolation
    p_s = slope_p .* dist1s + p_1;
    e_s = slope_e .* dist1s + e_1;
    rho_s = slope_rho .* dist1s + rho_1;
    a_s = slope_a .* dist1s + a_1;
    rho_u_s = slope_rho_u .* dist1s + rho_u_1;
    rho_v_s = slope_rho_v .* dist1s + rho_v_1;
    u_s = rho_u_s ./ rho_s;
    v_s = rho_v_s ./ rho_s;
    ang_inf = s.shock.beta(valid_ix, 1) - atan2(v_inf, u_inf);
    ang_s = s.shock.beta(valid_ix, 1) - atan2(v_s, u_s);

    %% Thermodynamic property evaluation (gamma, cv)
    if s.chemistry.is_chemistry_enabled
        if s.chemistry.chemical_equilibrium
            gamma_s = chemistry.eval_gamma_star(rho_s * s.freestream.rho_factor, e_s * s.freestream.energy_factor);
            cv_s = chemistry.eval_cv_star(rho_s * s.freestream.rho_factor, e_s * s.freestream.energy_factor);
        else
            gamma_s = chemistry.frozen.eval_gamma_star(rho_s * s.freestream.rho_factor, e_s * s.freestream.energy_factor);
            cv_s = chemistry.frozen.eval_cv_star(rho_s * s.freestream.rho_factor, e_s * s.freestream.energy_factor);
        end
        gamma_inf = s.freestream.gamma_star;
    else
        cv_s = 1;
        gamma_s = s.freestream.gamma;
        gamma_inf = s.freestream.gamma;
    end

    %% Riemann invariant along C- characteristic
    normal_velocity_inf = sin(ang_inf(valid_ix, 1)) .* sqrt(u_inf.^2 + v_inf.^2);
    normal_velocity_s = sin(ang_s(valid_ix, 1)) .* sqrt(u_s.^2 + v_s.^2);
    J_minus = p_s - rho_s .* a_s .* normal_velocity_s;

    %% Newton-Raphson iteration for post-shock pressure
    max_iter = 1000;
    tol = 1e-12;
    epsilon = 1e-4 * min(abs(p_s), [], "all");

    n_points = length(p_s);
    converged = false(n_points, 1);

    for iter = 1:max_iter
        % Evaluate residual at current pressure
        [residual] = COMPUTE_RESIDUAL_RIEMANN(p_infty, normal_velocity_inf, ...
            rho_inf, gamma_inf, p_s, a_s, gamma_s, J_minus);

        % Check convergence
        not_converged = ~converged & (abs(residual) > tol);
        if ~any(not_converged)
            break;
        end

        % Numerical Jacobian via forward difference
        p_s_perturbed = p_s;
        p_s_perturbed(not_converged) = p_s(not_converged) + epsilon;

        [residual_perturbed] = COMPUTE_RESIDUAL_RIEMANN(p_infty, normal_velocity_inf, ...
            rho_inf, gamma_inf, p_s_perturbed, a_s, gamma_s, J_minus);

        dR_dP = (residual_perturbed - residual) / epsilon;

        % Newton update
        delta_p = -residual ./ dR_dP;
        p_s(not_converged) = p_s(not_converged) + delta_p(not_converged);

        % Update convergence flags
        converged = converged | (abs(residual) <= tol);
    end

    if iter == max_iter && any(~converged)
        warning('%d points did not converge after %d iterations', sum(~converged), max_iter);
        keyboard
    end

    %% Post-shock state from converged pressure
    y = p_s ./ p_infty;

    % Density ratio from Rankine-Hugoniot
    A = (gamma_s + 1) ./ (gamma_s - 1);
    B = (gamma_inf + 1) ./ (gamma_inf - 1);
    rho_s = rho_inf .* (y .* A + 1) ./ (y + B);

    % Internal energy
    e_s = p_s ./ ((gamma_s - 1) .* rho_s);

    %% Shock velocity computation
    normal_velocity_inf = sin(ang_inf(valid_ix, 1)) .* sqrt(u_inf.^2 + v_inf.^2);
    temp1 = 2 * y ./ (gamma_s - 1) - 2 ./ (gamma_inf - 1);
    temp2 = y .* A + 1;
    temp3 = p_infty .* (y - 1) ./ rho_inf;
    s.shock.relative_increase_velocity = -normal_velocity_inf + sqrt(temp3 .* temp2 ./ temp1);

    %% Lab-frame velocity decomposition
    tangential_velocity_s = cos(ang_inf(valid_ix, 1)) .* sqrt(u_inf.^2 + v_inf.^2);
    normal_velocity_s = -s.shock.relative_increase_velocity + rho_inf ./ rho_s .* (s.shock.relative_increase_velocity + normal_velocity_inf);

    %% Shock speed in grid coordinates
    if strcmp(s.shock.formulation, "Lagrangian")
        s.shock.speed_x = -s.shock.relative_increase_velocity .* sin(s.shock.beta);
        s.shock.speed_y = s.shock.relative_increase_velocity .* cos(s.shock.beta);
    elseif strcmp(s.shock.formulation, "Eulerian")
        s.shock.speed_x = -s.shock.relative_increase_velocity .* sin(s.shock.beta);
        s.shock.speed_y = s.shock.relative_increase_velocity .* cos(s.shock.beta);
        s.shock.speed_x = s.shock.speed_x - s.shock.speed_y .* tan(pi/2 - s.shock.beta);
        s.shock.speed_y = zeros(size(s.shock.relative_increase_velocity));
    end

    %% Conservative variable reconstruction at the shock
    rho_u_s = rho_s .* (normal_velocity_s .* sin(s.shock.beta(valid_ix, 1)) ...
        + tangential_velocity_s .* cos(s.shock.beta(valid_ix, 1)));
    rho_v_s = rho_s .* (-normal_velocity_s .* cos(s.shock.beta(valid_ix, 1)) ...
        + tangential_velocity_s .* sin(s.shock.beta(valid_ix, 1)));
    rho_E_s = e_s .* rho_s + (rho_u_s.^2 + rho_v_s.^2) ./ (2 * rho_s);

    %% Store shock properties
    s.shock.properties.rho = rho_s;
    s.shock.properties.rho_u = rho_u_s;
    s.shock.properties.rho_v = rho_v_s;
    s.shock.properties.rho_E = rho_E_s;
    s.shock.properties.p = p_s;
    s.shock.properties.gamma_star = gamma_s;
    s.shock.properties.cv_star = cv_s;

    %% Update shocked cells via linear interpolation
    idx_0 = sub2ind(size(s.var.rho), valid_ix + 1, sc_idy + 1);

    s.var.rho(idx_0) = -slope_rho .* dist0s + rho_s;
    s.var.rho_u(idx_0) = -slope_rho_u .* dist0s + rho_u_s;
    s.var.rho_v(idx_0) = -slope_rho_v .* dist0s + rho_v_s;
    s.var.rho_E(idx_0) = -slope_rho_E .* dist0s + rho_E_s;
    s.var.p(idx_0) = -slope_p .* dist0s + p_s;

    if ~s.chemistry.chemical_equilibrium && s.chemistry.is_chemistry_enabled
        s.var.gamma_star(idx_0) = -slope_gamma_star .* dist0s + gamma_s;
        s.var.cv_star(idx_0) = -slope_cv_star .* dist0s + cv_s;
    end

    %% Extrapolate to ghost cells
    s = EXTRAPOLATE_CELLS_SHOCK(s);
end

%% ========================================================================
%  Riemann Invariant Residual
%  ========================================================================
function [residual] = COMPUTE_RESIDUAL_RIEMANN(p_infty, normal_velocity_inf, ...
        rho_inf, gamma_inf, p_s, a_s, gamma_s, J_minus)
% COMPUTE_RESIDUAL_RIEMANN  Evaluate the Riemann-invariant residual for
%   the Newton-Raphson pressure iteration.
%
%   residual = P_s - rho_s * a_s * U_n,s - J^-

    % Density ratio from Rankine-Hugoniot
    y = p_s ./ p_infty;
    A = (gamma_s + 1) ./ (gamma_s - 1);
    B = (gamma_inf + 1) ./ (gamma_inf - 1);
    rho_s = rho_inf .* (y .* A + 1) ./ (y + B);
    temp1 = 2 * y ./ (gamma_s - 1) - 2 ./ (gamma_inf - 1);
    temp2 = y .* A + 1;
    temp3 = p_infty .* (y - 1) ./ rho_inf;
    relative_increase_velocity = -normal_velocity_inf + sqrt(temp3 .* temp2 ./ temp1);

    % Post-shock normal velocity in lab frame
    normal_velocity_s = -relative_increase_velocity + rho_inf ./ rho_s .* (relative_increase_velocity + normal_velocity_inf);

    % Residual: R = P_s - rho_s * a_s * U_n,s - J^-
    residual = p_s - rho_s .* a_s .* normal_velocity_s - J_minus;
end

%% ========================================================================
%  Pressure-Based Rankine-Hugoniot Solver
%  ========================================================================
function s = UPDATE_SHOCK_JUMP_PROPERTIES_PRESSURE_3(s, chemistry)
% UPDATE_SHOCK_JUMP_PROPERTIES_PRESSURE_3  Solve R-H based on direct
%   pressure-jump interpolation (alternative to Riemann invariant method).

    %% Upstream conditions
    [p_infty, u_inf, v_inf, rho_inf, e_inf] = COMPUTE_UPSTREAM_CONDITIONS(s);
    s = UPDATE_FIELD_UPSTREAM(s);

    %% Index computation for shocked cells
    valid_ix = (1:s.mesh.Nchi)';
    sc_idy = s.shock.cell_indices(valid_ix, 1);

    % Replace sub2ind with direct linear indexing
    Nx_mesh = size(s.mesh.x, 1);
    lin_idx_0 = valid_ix + (sc_idy - 1) .* Nx_mesh;
    lin_idx_1 = valid_ix + (sc_idy - 2) .* Nx_mesh;
    lin_idx_2 = valid_ix + (sc_idy - 3) .* Nx_mesh;

    x_0 = s.mesh.x(lin_idx_0);
    y_0 = s.mesh.y(lin_idx_0);
    x_1 = s.mesh.x(lin_idx_1);
    y_1 = s.mesh.y(lin_idx_1);
    x_2 = s.mesh.x(lin_idx_2);
    y_2 = s.mesh.y(lin_idx_2);
    x_s = s.shock.points_x(valid_ix);
    y_s = s.shock.points_y(valid_ix);

    %% Distance computation for interpolation
    dist01 = sqrt((x_0 - x_1).^2 + (y_0 - y_1).^2);
    dist12 = sqrt((x_1 - x_2).^2 + (y_1 - y_2).^2);
    dist1s = ((x_s - x_1).*(x_1 - x_2) + (y_s - y_1).*(y_1 - y_2)) ./ dist12;
    dist0s = ((x_s - x_0).*(x_0 - x_1) + (y_s - y_0).*(y_0 - y_1)) ./ dist01;

    %% Flow field gradient computation
    % Single computation for size, direct linear indexing
    Nx_var = size(s.var.p, 1);
    idx_v0 = (valid_ix + 1) + (sc_idy) .* Nx_var;       % sc_idy + 1 - 0
    idx_v1 = (valid_ix + 1) + (sc_idy - 1) .* Nx_var;  % sc_idy + 1 - 1
    idx_v2 = (valid_ix + 1) + (sc_idy - 2) .* Nx_var;  % sc_idy + 1 - 2

    p_0 = s.var.p(idx_v0);
    p_1 = s.var.p(idx_v1);
    p_2 = s.var.p(idx_v2);
    slope_p = (p_1 - p_2) ./ dist12;

    rho_1 = s.var.rho(idx_v1);
    rho_2 = s.var.rho(idx_v2);
    slope_rho = (rho_1 - rho_2) ./ dist12;

    rho_u_1 = s.var.rho_u(idx_v1);
    rho_u_2 = s.var.rho_u(idx_v2);
    slope_rho_u = (rho_u_1 - rho_u_2) ./ dist12;

    rho_v_1 = s.var.rho_v(idx_v1);
    rho_v_2 = s.var.rho_v(idx_v2);
    slope_rho_v = (rho_v_1 - rho_v_2) ./ dist12;

    rho_E_1 = s.var.rho_E(idx_v1);
    rho_E_2 = s.var.rho_E(idx_v2);
    slope_rho_E = (rho_E_1 - rho_E_2) ./ dist12;

    %% Internal energy gradient
    % Compute internal energy ONLY at the required 1D locations
    % instead of across the entire 2D flow field.
    rho_vsq_1 = (rho_u_1.^2 + rho_v_1.^2) ./ rho_1;
    e_1 = (rho_E_1 - rho_vsq_1 / 2) ./ rho_1;

    rho_vsq_2 = (rho_u_2.^2 + rho_v_2.^2) ./ rho_2;
    e_2 = (rho_E_2 - rho_vsq_2 / 2) ./ rho_2;

    slope_e = (e_1 - e_2) ./ dist12;

    %% Non-equilibrium chemistry gradients
    if ~s.chemistry.chemical_equilibrium && s.chemistry.is_chemistry_enabled
        gamma_star_1 = s.var.gamma_star(idx_v1);
        gamma_star_2 = s.var.gamma_star(idx_v2);
        slope_gamma_star = (gamma_star_1 - gamma_star_2) ./ dist12;
        
        cv_star_1 = s.var.cv_star(idx_v1);
        cv_star_2 = s.var.cv_star(idx_v2);
        slope_cv_star = (cv_star_1 - cv_star_2) ./ dist12;
    end

    %% Interpolated state at shock location (initial values for gamma evaluation)
    e_s = slope_e .* dist1s + e_1;
    rho_s = slope_rho .* dist1s + rho_1;

    %% Thermodynamic property evaluation (gamma, cv)
    if s.chemistry.is_chemistry_enabled
        if s.chemistry.chemical_equilibrium
            gamma_s = chemistry.eval_gamma_star(rho_s * s.freestream.rho_factor, e_s * s.freestream.energy_factor);
            cv_s = chemistry.eval_cv_star(rho_s * s.freestream.rho_factor, e_s * s.freestream.energy_factor);
        else
            gamma_s = chemistry.frozen.eval_gamma_star(rho_s * s.freestream.rho_factor, e_s * s.freestream.energy_factor);
            cv_s = chemistry.frozen.eval_cv_star(rho_s * s.freestream.rho_factor, e_s * s.freestream.energy_factor);
        end
        gamma_inf = s.freestream.gamma_star;
    else
        gamma_s = s.freestream.gamma;
        gamma_inf = s.freestream.gamma;
        cv_s = 1;
    end

    %% Pressure interpolation at shock
    p_s = slope_p .* dist1s + p_1;

    %% Pressure ratio and validity check
    y = p_s ./ p_infty;

    if any(y < 1)
        warning('ERROR: Mach number less than 1');
        warning('Add smoothing to shock and increase viscosity_scaling');
        keyboard
    end

    %% Density and energy from Rankine-Hugoniot
    A = (gamma_s + 1) ./ (gamma_s - 1);
    B = (gamma_inf + 1) ./ (gamma_inf - 1);
    rho_s = rho_inf .* (y .* A + 1) ./ (y + B);
    e_s = p_s ./ ((gamma_s - 1) .* rho_s);

    %% Shock velocity computation
    ang_inf = s.shock.beta(valid_ix, 1) - atan2(v_inf, u_inf);
    normal_velocity_inf = sin(ang_inf(valid_ix, 1)) .* sqrt(u_inf.^2 + v_inf.^2);
    temp1 = 2 * y ./ (gamma_s - 1) - 2 ./ (gamma_inf - 1);
    temp2 = y .* A + 1;
    temp3 = p_infty .* (y - 1) ./ rho_inf;
    s.shock.relative_increase_velocity = -normal_velocity_inf + sqrt(temp3 .* temp2 ./ temp1);

    %% Lab-frame velocity decomposition
    tangential_velocity_s = cos(ang_inf(valid_ix, 1)) .* sqrt(u_inf.^2 + v_inf.^2);
    normal_velocity_s = -s.shock.relative_increase_velocity + rho_inf ./ rho_s .* (s.shock.relative_increase_velocity + normal_velocity_inf);

    %% Shock speed in grid coordinates
    if strcmp(s.shock.formulation, "Lagrangian")
        s.shock.speed_x = -s.shock.relative_increase_velocity .* sin(s.shock.beta);
        s.shock.speed_y = s.shock.relative_increase_velocity .* cos(s.shock.beta);
    elseif strcmp(s.shock.formulation, "Eulerian")
        s.shock.speed_x = -s.shock.relative_increase_velocity .* sin(s.shock.beta);
        s.shock.speed_y = s.shock.relative_increase_velocity .* cos(s.shock.beta);
        s.shock.speed_x = s.shock.speed_x - s.shock.speed_y .* tan(pi/2 - s.shock.beta);
        s.shock.speed_y = zeros(size(s.shock.relative_increase_velocity));
    end

    %% Conservative variable reconstruction at the shock
    rho_u_s = rho_s .* (normal_velocity_s .* sin(s.shock.beta(valid_ix, 1)) ...
        + tangential_velocity_s .* cos(s.shock.beta(valid_ix, 1)));
    rho_v_s = rho_s .* (-normal_velocity_s .* cos(s.shock.beta(valid_ix, 1)) ...
        + tangential_velocity_s .* sin(s.shock.beta(valid_ix, 1)));
    rho_E_s = e_s .* rho_s + (rho_u_s.^2 + rho_v_s.^2) ./ (2 * rho_s);

    %% Store shock properties
    s.shock.properties.rho = rho_s;
    s.shock.properties.rho_u = rho_u_s;
    s.shock.properties.rho_v = rho_v_s;
    s.shock.properties.rho_E = rho_E_s;
    s.shock.properties.p = p_s;
    s.shock.properties.gamma_star = gamma_s;
    s.shock.properties.cv_star = cv_s;

    %% Update shocked cells via linear interpolation
    idx_update = (valid_ix + 1) + (sc_idy) .* Nx_var; % Same as idx_v0

    s.var.rho(idx_update) = -slope_rho .* dist0s + rho_s;
    s.var.rho_u(idx_update) = -slope_rho_u .* dist0s + rho_u_s;
    s.var.rho_v(idx_update) = -slope_rho_v .* dist0s + rho_v_s;
    s.var.rho_E(idx_update) = -slope_rho_E .* dist0s + rho_E_s;
    s.var.p(idx_update) = -slope_p .* dist0s + p_s;

    if ~s.chemistry.chemical_equilibrium && s.chemistry.is_chemistry_enabled
        s.var.gamma_star(idx_update) = -slope_gamma_star .* dist0s + gamma_s;
        s.var.cv_star(idx_update) = -slope_cv_star .* dist0s + cv_s;
    end

    %% Extrapolate to ghost cells
    s = EXTRAPOLATE_CELLS_SHOCK(s);
end

%% ========================================================================
%  Mach-Number-Based Pressure Solver (Legacy)
%  ========================================================================
function s = UPDATE_SHOCK_JUMP_PROPERTIES_PRESSURE(s)
% UPDATE_SHOCK_JUMP_PROPERTIES_PRESSURE  Legacy solver using Mach-number-
%   based pressure jump for shock fitting with or without shock feedback.

    if s.shock.feedback
        %% Precompute factors
        gamma = s.freestream.gamma;
        factor_M = (gamma + 1) ./ (2 * gamma);
        factor_add = (gamma - 1) ./ (2 * gamma);
        shock_diff_factor = s.mesh.Nchi^2 / s.Re_shock;

        [p_infty, u_inf, v_inf, rho_inf, e_inf] = COMPUTE_UPSTREAM_CONDITIONS(s);
        a_s_inf = sqrt(s.freestream.gamma * p_infty ./ rho_inf);
        u_mag = sqrt(v_inf.^2 + u_inf.^2);
        Mach_inf = u_mag ./ a_s_inf;

        %% Index computation for shocked cells
        valid_ix = (1:s.mesh.Nchi)';
        sc_idy = s.shock.cell_indices(valid_ix, 1);

        lin_idx_0 = sub2ind(size(s.mesh.x), valid_ix, sc_idy);
        lin_idx_1 = sub2ind(size(s.mesh.x), valid_ix, sc_idy - 1);
        lin_idx_2 = sub2ind(size(s.mesh.x), valid_ix, sc_idy - 2);

        x_0 = s.mesh.x(lin_idx_0);
        y_0 = s.mesh.y(lin_idx_0);
        x_1 = s.mesh.x(lin_idx_1);
        y_1 = s.mesh.y(lin_idx_1);
        x_2 = s.mesh.x(lin_idx_2);
        y_2 = s.mesh.y(lin_idx_2);
        x_s = s.shock.points_x(valid_ix);
        y_s = s.shock.points_y(valid_ix);

        %% Distance computation
        dist01 = sqrt((x_0 - x_1).^2 + (y_0 - y_1).^2);
        dist12 = sqrt((x_1 - x_2).^2 + (y_1 - y_2).^2);
        dist1s = ((x_s - x_1).*(x_1 - x_2) + (y_s - y_1).*(y_1 - y_2)) ./ dist12;
        dist0s = ((x_s - x_0).*(x_0 - x_1) + (y_s - y_0).*(y_0 - y_1)) ./ dist01;

        %% Pressure interpolation
        p_1 = s.var.p(sub2ind(size(s.var.p), valid_ix + 1, sc_idy + 1 - 1));
        p_2 = s.var.p(sub2ind(size(s.var.p), valid_ix + 1, sc_idy + 1 - 2));
        slope_p = (p_1 - p_2) ./ dist12;
        p_s = slope_p .* dist1s + p_1;

        %% Flow field gradients
        rho_1 = s.var.rho(sub2ind(size(s.var.rho), valid_ix + 1, sc_idy + 1 - 1));
        rho_2 = s.var.rho(sub2ind(size(s.var.rho), valid_ix + 1, sc_idy + 1 - 2));
        slope_rho = (rho_1 - rho_2) ./ dist12;

        rho_u_1 = s.var.rho_u(sub2ind(size(s.var.rho_u), valid_ix + 1, sc_idy + 1 - 1));
        rho_u_2 = s.var.rho_u(sub2ind(size(s.var.rho_u), valid_ix + 1, sc_idy + 1 - 2));
        slope_rho_u = (rho_u_1 - rho_u_2) ./ dist12;

        rho_v_1 = s.var.rho_v(sub2ind(size(s.var.rho_v), valid_ix + 1, sc_idy + 1 - 1));
        rho_v_2 = s.var.rho_v(sub2ind(size(s.var.rho_v), valid_ix + 1, sc_idy + 1 - 2));
        slope_rho_v = (rho_v_1 - rho_v_2) ./ dist12;

        rho_E_1 = s.var.rho_E(sub2ind(size(s.var.rho_E), valid_ix + 1, sc_idy + 1 - 1));
        rho_E_2 = s.var.rho_E(sub2ind(size(s.var.rho_E), valid_ix + 1, sc_idy + 1 - 2));
        slope_rho_E = (rho_E_1 - rho_E_2) ./ dist12;

        %% Shock Mach number and velocity
        s.shock.M = real(sqrt(p_s ./ p_infty .* factor_M + factor_add));
        if any(s.shock.M < 1)
            warning('ERROR: Mach number less than 1');
            warning('Add smoothing to shock and increase viscosity_scaling');
            keyboard
        end

        ang = s.shock.beta(valid_ix, 1) - atan2(v_inf, u_inf);
        s.shock.relative_increase_velocity = (s.shock.M - sin(ang) .* Mach_inf) .* a_s_inf;

        %% Shock speed in grid coordinates
        if strcmp(s.shock.formulation, "Lagrangian")
            s.shock.speed_x = -s.shock.relative_increase_velocity .* sin(s.shock.beta);
            s.shock.speed_y = s.shock.relative_increase_velocity .* cos(s.shock.beta);
        elseif strcmp(s.shock.formulation, "Eulerian")
            s.shock.speed_x = -s.shock.relative_increase_velocity .* sin(s.shock.beta);
            s.shock.speed_y = s.shock.relative_increase_velocity .* cos(s.shock.beta);
            s.shock.speed_x = s.shock.speed_x - s.shock.speed_y .* tan(pi/2 - s.shock.beta);
            s.shock.speed_y = zeros(size(s.shock.relative_increase_velocity));
        end

        %% Conservative variable reconstruction (supersonic cells only)
        valid_ix = find(s.shock.M(:, 1) > 1);

        if ~isempty(valid_ix)
            M = s.shock.M(valid_ix, 1);
            rho_s = (gamma + 1) .* M.^2 ./ ((gamma - 1) .* M.^2 + 2) .* rho_inf;

            tangential_velocity = cos(ang(valid_ix, 1)) .* sqrt(u_inf.^2 + v_inf.^2);
            normal_velocity = (M .* a_s_inf) .* (rho_inf ./ rho_s);

            shock_speed_x = s.shock.speed_x(valid_ix, 1);
            shock_speed_y = s.shock.speed_y(valid_ix, 1);
            rho_u_s = rho_s .* (normal_velocity .* sin(s.shock.beta(valid_ix, 1)) ...
                + tangential_velocity .* cos(s.shock.beta(valid_ix, 1)) + shock_speed_x);
            rho_v_s = rho_s .* (-normal_velocity .* cos(s.shock.beta(valid_ix, 1)) ...
                + tangential_velocity .* sin(s.shock.beta(valid_ix, 1)) + shock_speed_y);

            rho_E_s = p_s ./ (gamma - 1) + (rho_u_s.^2 + rho_v_s.^2) ./ (2 * rho_s);

            % Store properties
            s.shock.properties.rho = rho_s;
            s.shock.properties.rho_u = rho_u_s;
            s.shock.properties.rho_v = rho_v_s;
            s.shock.properties.rho_E = rho_E_s;
            s.shock.properties.p = p_s;
            s.shock.properties.gamma_star = gamma_s;
            s.shock.properties.cv_star = cv_s;

            % Update shocked cells via linear interpolation
            idx_0 = sub2ind(size(s.var.rho), valid_ix + 1, sc_idy + 1);

            s.var.rho(idx_0) = -slope_rho .* dist0s + rho_s;
            s.var.rho_u(idx_0) = -slope_rho_u .* dist0s + rho_u_s;
            s.var.rho_v(idx_0) = -slope_rho_v .* dist0s + rho_v_s;
            s.var.rho_E(idx_0) = -slope_rho_E .* dist0s + rho_E_s;
            s.var.p(idx_0) = -slope_p .* dist0s + p_s;
        end
    else
        %% Passive shock (no feedback) -- simple extrapolation
        valid_ix = (1:s.mesh.Nchi)';
        sc_idy = s.shock.cell_indices(valid_ix, 1);

        idx_0 = sub2ind(size(s.var.rho_E), valid_ix + 1, sc_idy + 1);
        idx_1 = sub2ind(size(s.var.rho_E), valid_ix + 1, sc_idy + 1 - 1);
        idx_2 = sub2ind(size(s.var.rho_E), valid_ix + 1, sc_idy + 1 - 2);

        s.var.rho_E(idx_0) = 2 * s.var.rho_E(idx_1) - s.var.rho_E(idx_2);
    end

    %% Extrapolate to ghost cells
    s = EXTRAPOLATE_CELLS_SHOCK(s);
end
