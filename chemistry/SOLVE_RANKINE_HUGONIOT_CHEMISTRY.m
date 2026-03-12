function [rho_2, e_2, w_2] = SOLVE_RANKINE_HUGONIOT_CHEMISTRY(chemistry, w_1, rho_1, e_1, varargin)
% SOLVE_RANKINE_HUGONIOT_CHEMISTRY  Vectorized Rankine-Hugoniot shock solver
%   with real-gas chemistry.
%
%   Solves N independent 2x2 nonlinear systems simultaneously using a
%   single fsolve call. The key optimisation is a single vectorised
%   chemistry evaluation per solver iteration instead of N separate calls.
%
%   [rho_2, e_2, w_2] = SOLVE_RANKINE_HUGONIOT_CHEMISTRY(chemistry, w_1, rho_1, e_1)
%   [rho_2, e_2, w_2] = SOLVE_RANKINE_HUGONIOT_CHEMISTRY(chemistry, w_1, rho_1, e_1, 'InitialGuess', [p0 q0])
%
%   Inputs:
%       chemistry  - Fitted chemistry structure from FIT_CHEMISTRY
%                    (must contain eval_gamma_star)
%       w_1        - Upstream velocity [m/s] (scalar or N-vector)
%       rho_1      - Upstream density  [kg/m^3] (scalar or N-vector)
%       e_1        - Upstream specific internal energy [J/kg] (scalar or N-vector)
%
%   Optional Name-Value Pairs:
%       'InitialGuess' - [p0, q0] initial density and energy ratios
%                        (default: [8, 1.2])
%
%   Outputs:
%       rho_2 - Post-shock density  [kg/m^3] (N-vector)
%       e_2   - Post-shock specific internal energy [J/kg] (N-vector)
%       w_2   - Post-shock velocity [m/s] (N-vector)
%
%   Notes:
%       The solver introduces normalised variables p = rho_2/rho_1 and
%       q = e_2/w_1^2 so that the conservation equations become
%       dimensionless.  A physics-based initial guess derived from the
%       ideal-gas normal-shock relations is used to accelerate convergence.
%
% Part of: Hypersonics Stability MATLAB Solver - Chemistry Module

    %% Parse optional inputs
    p = inputParser;
    addParameter(p, 'InitialGuess', [8, 1.2], @(x) isnumeric(x) && length(x) == 2);
    parse(p, varargin{:});

    initial_guess = p.Results.InitialGuess;

    %% Prepare input vectors
    w_1 = w_1(:);
    rho_1 = rho_1(:);
    e_1 = e_1(:);
    N = length(w_1);

    if isscalar(rho_1), rho_1 = repmat(rho_1, N, 1); end
    if isscalar(e_1), e_1 = repmat(e_1, N, 1); end

    %% Pre-compute upstream quantities (single vectorised chemistry call)
    gamma_1 = chemistry.eval_gamma_star(rho_1, e_1);

    e1_over_w1_sq = e_1 ./ (w_1.^2);
    const_1 = (gamma_1 - 1) .* e1_over_w1_sq + 1;
    const_2 = gamma_1 .* e1_over_w1_sq + 0.5;

    %% Build initial guess
    x0 = create_initial_guess(N, gamma_1, w_1, e_1);

    %% Configure and run fsolve
    options = optimoptions('fsolve', ...
        'Display', 'off', ...
        'Algorithm', 'trust-region-dogleg', ...
        'MaxIterations', 100, ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8);

    x_sol = fsolve(@equations_vectorized, x0, options);

    %% Extract and compute post-shock state
    vars_sol = reshape(x_sol, 2, N)';
    p_sol = vars_sol(:, 1);
    q_sol = vars_sol(:, 2);

    rho_2 = p_sol .* rho_1;
    e_2 = q_sol .* (w_1.^2);
    w_2 = w_1 ./ p_sol;

    %% Nested function: vectorised residual evaluation
    function residual = equations_vectorized(x)
        % Reshape x to [N x 2] matrix: [p1 q1; p2 q2; ... ; pN qN]
        vars = reshape(x, 2, N)';
        p_vals = vars(:, 1);
        q_vals = vars(:, 2);

        % Post-shock state for all systems
        rho_2_vals = p_vals .* rho_1;
        e_2_vals = q_vals .* (w_1.^2);

        % Single vectorised chemistry call -- major speed-up
        gamma_2 = chemistry.eval_gamma_star(rho_2_vals, e_2_vals);

        % Residuals for conservation equations
        eq1 = const_1 - ((gamma_2 - 1) .* p_vals .* q_vals + 1 ./ p_vals);
        eq2 = const_2 - (gamma_2 .* q_vals + 1 ./ (2 * p_vals.^2));

        % Interleave: [eq1_1, eq2_1, eq1_2, eq2_2, ..., eq1_N, eq2_N]
        residual = reshape([eq1'; eq2'], [], 1);
    end

end

%% ========================================================================
%  Helper: physics-based initial guess
%  ========================================================================
function x0 = create_initial_guess(N, gamma_1, w_1, e_1)
% CREATE_INITIAL_GUESS  Compute a physics-based starting point for fsolve.
%
%   Uses the ideal-gas normal-shock relations to estimate the density
%   ratio (p) and energy ratio (q) for improved convergence.

    % Approximate Mach number squared
    M1_sq_approx = w_1.^2 ./ (gamma_1 .* (gamma_1 - 1) .* e_1);

    % Density-ratio guess from normal-shock relation
    p_guess = (gamma_1 + 1) .* M1_sq_approx ./ ((gamma_1 - 1) .* M1_sq_approx + 2);

    % Energy-ratio guess
    q_guess = (1/2 + e_1 ./ w_1.^2) ./ gamma_1;

    if N == 1
        x0 = [p_guess; q_guess];
    else
        x0 = reshape([p_guess'; q_guess'], [], 1);
    end
end
