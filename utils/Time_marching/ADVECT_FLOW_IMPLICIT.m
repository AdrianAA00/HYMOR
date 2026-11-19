function s = ADVECT_FLOW_IMPLICIT(s, chemistry)
% ADVECT_FLOW_IMPLICIT  Advance the flow s by one implicit step (steady driver).
%
%   s = ADVECT_FLOW_IMPLICIT(s, chemistry)
%
%   Dispatches to one of two implicit solvers, selected by the field
%   s.time_integration.implicit_solver:
%
%     "Newton"      (default) - Pseudo-transient continuation: one damped
%                               Newton step per call on the backward-Euler /
%                               steady residual, using the sparse Jacobian L
%                               built by LINEARIZE_EQUATIONS_NO_DISCONTINUITY.
%                               The pseudo-time step (CFL) is ramped up as the
%                               residual decreases following the Switched
%                               Evolution Relaxation rule of Mulder & van Leer
%                               (1985), recovering Newton's method as CFL->inf
%                               (pseudo-transient continuation, Kelley & Keyes
%                               1998). This is the recommended way to obtain
%                               steady base flows.
%
%     "Relaxation"            - Legacy first-order backward Euler solved by a
%                               relaxed fixed-point (Picard) iteration with
%                               adaptive under-relaxation (ADAPT_RELAXATION +
%                               INTEGRATION_IMPLICIT). Kept for backward
%                               compatibility; cannot exploit large CFL.
%
%   Inputs:
%       s         - struct : Flow state (conservative variables, grid, shock,
%                            time-integration parameters).
%       chemistry - struct : Chemistry / thermodynamic model data.
%
%   Outputs:
%       s         - struct : Flow advanced by one (pseudo-)time step, with BCs,
%                            thermodynamic properties and the density limiter
%                            applied, and time/iteration counters incremented.
%
%   See also: ADVECT_FLOW_EXPLICIT, LINEARIZE_EQUATIONS_NO_DISCONTINUITY,
%             NON_LINEAR_DYNAMICS_NO_DISCONTINUITY, CFL_TIMESTEP.
%
% Part of: Hypersonics Stability MATLAB Solver - Time Marching Module

    solver = GET_OPT(s.time_integration, 'implicit_solver', "Newton");

    switch string(solver)
        case "Newton"
            s = ADVECT_FLOW_IMPLICIT_NEWTON(s, chemistry);
        case "Relaxation"
            s = ADVECT_FLOW_IMPLICIT_RELAXATION(s, chemistry);
        otherwise
            error("ADVECT_FLOW_IMPLICIT:unknownSolver", ...
                  "Unknown s.time_integration.implicit_solver = '%s' " + ...
                  "(use 'Newton' or 'Relaxation').", string(solver));
    end
end


% =========================================================================
%  Newton / pseudo-transient-continuation solver
% =========================================================================
function s = ADVECT_FLOW_IMPLICIT_NEWTON(s, chemistry)
% One damped-Newton (pseudo-transient continuation) step toward steady state.
%
%   Solves, for the increment dq of the interior conservative variables,
%
%       (I/dtau - L) dq = f(q^k),     q^{k+1} = q^k + alpha*dq
%
%   where f(q) = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY (the volume-normalised
%   flux residual, masked by the flow cells when shock fitting is active) and
%   L = df/dq is the sparse Jacobian assembled by
%   LINEARIZE_EQUATIONS_NO_DISCONTINUITY. The pseudo-time step dtau is the
%   explicit CFL-stable step scaled by the ramped CFL number; alpha is a
%   positivity safeguard (<=1). At the fixed point f(q)=0 (steady state) and,
%   as dtau->inf, the update reduces to Newton's method on f(q)=0.

    Nchi = s.mesh.Nchi;
    Neta = s.mesh.Neta;
    n    = Nchi * Neta;
    n4   = 4 * n;

    %% Solver / ramp options (defaults applied if absent from the input file)
    CFL_0     = s.time_integration.CFL;                                  % base CFL
    CFL_max   = GET_OPT(s.time_integration, 'CFL_max',            1e5);  % SER ceiling
    p_ramp    = GET_OPT(s.time_integration, 'CFL_ramp_exponent',  1.0);  % SER exponent
    growth    = GET_OPT(s.time_integration, 'CFL_growth',         2.0);  % max growth/step
    jac_every = GET_OPT(s.time_integration, 'jacobian_refresh',   1);    % rebuild L every k steps
    rho_frac  = GET_OPT(s.time_integration, 'newton_min_rho_fraction', 0.2); % positivity floor

    if ~s.chemistry.chemical_equilibrium
        warning("ADVECT_FLOW_IMPLICIT:nonEquilibrium", ...
            "Newton implicit solver linearises only (rho, rho_u, rho_v, rho_E); " + ...
            "non-equilibrium fields (gamma_star, cv_star) are held frozen. " + ...
            "Use Explicit_RK4 (or implicit_solver='Relaxation') for " + ...
            "non-equilibrium chemistry.");
    end

    %% 1) Flow cells + nonlinear steady residual f(q^k)
    s = UPDATE_FLOW_CELLS(s, chemistry);

    s.linearize = false;                                  % full physics in the residual
    s = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s, chemistry);   % sets s.flux, thermo, shock.speed

    if s.shock.enabled
        fc = s.shock.flow_cells > 0;                      % logical Nchi x Neta
    else
        fc = true(Nchi, Neta);
    end

    R_rho   = s.flux.rho   .* fc;
    R_rho_u = s.flux.rho_u .* fc;
    R_rho_v = s.flux.rho_v .* fc;
    R_rho_E = s.flux.rho_E .* fc;
    R_vec   = [R_rho(:); R_rho_u(:); R_rho_v(:); R_rho_E(:)];

    %% Per-variable relative residual (convergence monitor + plotting)
    rho_i  = s.var.rho(2:end-1,2:end-1);
    rhou_i = s.var.rho_u(2:end-1,2:end-1);
    rhov_i = s.var.rho_v(2:end-1,2:end-1);
    rhoE_i = s.var.rho_E(2:end-1,2:end-1);
    residual.rho   = rms(s.flux.rho(fc))   / max(rms(rho_i(fc)),  eps);
    residual.rho_u = rms(s.flux.rho_u(fc)) / max(rms(rhou_i(fc)), eps);
    residual.rho_v = rms(s.flux.rho_v(fc)) / max(rms(rhov_i(fc)), eps);
    residual.rho_E = rms(s.flux.rho_E(fc)) / max(rms(rhoE_i(fc)), eps);

    %% 2) Explicit-stable time step + SER ramp of the pseudo-time step
    s = CFL_TIMESTEP(s);                                  % dt at the base CFL
    dt_explicit = s.time_integration.dt;

    Rn = norm(R_vec) / sqrt(n4);                          % scalar residual norm
    if ~isfield(s.time_integration, 'residual_ref') || ...
            ~isfinite(s.time_integration.residual_ref) || s.time_integration.residual_ref <= 0
        s.time_integration.residual_ref = Rn;            % first step -> reference
        s.time_integration.CFL_current   = CFL_0;
    end
    R0       = s.time_integration.residual_ref;
    CFL_prev = s.time_integration.CFL_current;

    CFL_cur = CFL_0 * (R0 / max(Rn, realmin))^p_ramp;     % Switched Evolution Relaxation
    CFL_cur = min(CFL_cur, CFL_max);                      % ceiling
    CFL_cur = min(CFL_cur, growth * CFL_prev);            % limit growth per step
    CFL_cur = max(CFL_cur, CFL_0);                        % never below base
    s.time_integration.CFL_current = CFL_cur;

    dt_pstc = dt_explicit * (CFL_cur / CFL_0);            % ramped pseudo-time step

    %% 3) Sparse Jacobian L = df/dq (rebuilt every jac_every steps)
    L = NEWTON_JACOBIAN(s, chemistry, jac_every, s.time_integration.iter == 0);

    %% 4) Solve (I/dtau - L) dq = f(q^k)  with MATLAB's sparse LU (UMFPACK)
    M  = speye(n4) ./ dt_pstc - L;
    dq = M \ R_vec;

    if ~all(isfinite(dq))
        warning("ADVECT_FLOW_IMPLICIT:singularSolve", ...
            "Newton linear solve produced non-finite values; skipping update " + ...
            "this step (reduce CFL_max or CFL_growth).");
        dq = zeros(n4, 1);
    end

    %% 5) Positivity-safeguarded update of the interior conservative variables
    d_rho   = reshape(dq(0*n+1:1*n), Nchi, Neta) .* fc;
    d_rho_u = reshape(dq(1*n+1:2*n), Nchi, Neta) .* fc;
    d_rho_v = reshape(dq(2*n+1:3*n), Nchi, Neta) .* fc;
    d_rho_E = reshape(dq(3*n+1:4*n), Nchi, Neta) .* fc;

    alpha = 1;                                            % damped-Newton step length
    neg   = d_rho < 0;
    if any(neg(:))
        bound = (rho_frac - 1) * rho_i(neg) ./ d_rho(neg);   % keep rho >= rho_frac*rho_old
        alpha = max(0, min(1, min(bound)));
    end

    s.var.rho(2:end-1,2:end-1)   = rho_i  + alpha * d_rho;
    s.var.rho_u(2:end-1,2:end-1) = rhou_i + alpha * d_rho_u;
    s.var.rho_v(2:end-1,2:end-1) = rhov_i + alpha * d_rho_v;
    s.var.rho_E(2:end-1,2:end-1) = rhoE_i + alpha * d_rho_E;

    %% 6) Advance the shock explicitly on the (un-ramped) stable time scale
    if s.shock.enabled
        dt_shock = dt_explicit * s.shock.relaxation;
        s.shock.points_x = s.shock.points_x + s.shock.speed_x * dt_shock;
        s.shock.points_y = s.shock.points_y + s.shock.speed_y * dt_shock;
    end

    %% 7) Post-step: boundary conditions, thermodynamics, density limiter
    s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry);
    if s.shock.enabled
        s = UPDATE_SHOCK_BC(s, chemistry);
    end
    s = APPLY_BOUNDARY_CONDITIONS(s, chemistry);
    s = MIN_RHO(s);

    %% 8) Advance counters and store diagnostics
    s.time_integration.dt       = dt_pstc;
    s.time_integration.t        = s.time_integration.t + dt_pstc;
    s.time_integration.iter     = s.time_integration.iter + 1;
    s.time_integration.residual = residual;
    s.count_implicit_iterations = 1;     % one Newton step per call
end


function L = NEWTON_JACOBIAN(s, chemistry, jac_every, reset)
% NEWTON_JACOBIAN  Cached sparse flow Jacobian L = df/dq (A11 block only).
%
%   Rebuilds L (via LINEARIZE_EQUATIONS_NO_DISCONTINUITY, with the shock DOFs
%   disabled and chemistry frozen during the finite differences) every
%   jac_every calls, on a grid-size change, or when reset is true. Between
%   rebuilds the cached factorisable matrix is reused (modified Newton).

    persistent L_cache age sz

    n4 = 4 * s.mesh.Nchi * s.mesh.Neta;

    if reset
        L_cache = [];
        age     = [];
        sz      = [];
    end

    rebuild = isempty(L_cache) || isempty(sz) || sz ~= n4 || ...
              isempty(age) || mod(age, max(1, jac_every)) == 0;

    if rebuild
        s_lin = s;
        s_lin.linearize = true;                          % frozen-chemistry Jacobian
        if ~isfield(s_lin, 'stability_analysis') || ~isstruct(s_lin.stability_analysis)
            s_lin.stability_analysis = struct();
        end
        s_lin.stability_analysis.perturb_shock = false;  % flow-flow (A11) block only
        s_lin = CLEAR_AZIMUTHAL_MODE(s_lin);              % axisymmetric (4-variable) Jacobian
        if ~isfield(s_lin.stability_analysis, 'perturbation_magnitude')
            s_lin.stability_analysis.perturbation_magnitude = 1e-6;
        end
        L_cache = LINEARIZE_EQUATIONS_NO_DISCONTINUITY(s_lin, chemistry);
        sz      = n4;
        age     = 0;
    end

    age = age + 1;
    L   = L_cache;
end


function v = GET_OPT(strct, name, default)
% GET_OPT  Return strct.(name) if present, otherwise default.
    if isfield(strct, name)
        v = strct.(name);
    else
        v = default;
    end
end


% =========================================================================
%  Legacy relaxed fixed-point (Picard) backward-Euler solver
% =========================================================================
function s = ADVECT_FLOW_IMPLICIT_RELAXATION(s, chemistry)
% Backward Euler solved by an adaptively under-relaxed fixed-point iteration.
% (Original ADVECT_FLOW_IMPLICIT behaviour, retained for compatibility.)

    %% Initialise implicit iteration
    residual          = 100;
    residual_previous = 0;

    s          = UPDATE_FLOW_CELLS(s, chemistry);
    s = CFL_TIMESTEP(s); % Recompute CFL-based time step for the next iteration
    solution_temp     = s;
    count_implicit_iterations = 0;

    %% Implicit iteration loop
    while residual > s.time_integration.tolerance && count_implicit_iterations < s.time_integration.max_iter_implicit
        s = ADAPT_RELAXATION(residual, residual_previous, s);

        % Compute PDE right-hand side
        solution_temp = PDE(solution_temp, chemistry);

        % Implicit time integration and residual evaluation
        [solution_temp, residual] = INTEGRATION_IMPLICIT(s, solution_temp);
        count_implicit_iterations = count_implicit_iterations + 1;
        residual_previous = residual;
    end
    s      = solution_temp;

    %% Post-step: boundary conditions and corrections
    s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry);
    if s.shock.enabled
        s = UPDATE_SHOCK_BC(s, chemistry);
    end
    s = APPLY_BOUNDARY_CONDITIONS(s, chemistry);

    %% Enforce minimum density limiter
    s = MIN_RHO(s);

    %% Advance time and iteration counters
    s.time_integration.t    = s.time_integration.t + s.time_integration.dt;
    s.time_integration.iter = s.time_integration.iter + 1;
    s.count_implicit_iterations = count_implicit_iterations;
end


function s = ADAPT_RELAXATION(residual, residual_previous, s)
% ADAPT_RELAXATION  Adjust the variable relaxation factor based on residual trend.
%
%   s = ADAPT_RELAXATION(residual, residual_previous, s)
%
%   Uses a sigmoid mapping of the residual change to adapt the relaxation
%   factor between 0.5 and 0.99 for the implicit iteration.

    C    = 0.5;
    a    = C * (residual - residual_previous);
    temp = 1.7 * s.time_integration.relax_factor / (1 + exp(-a));
    m    = min(0.99, temp);
    s.relax_factor_variable = max(0.5, m);
end
