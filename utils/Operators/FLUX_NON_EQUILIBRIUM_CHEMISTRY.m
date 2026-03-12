function s = FLUX_NON_EQUILIBRIUM_CHEMISTRY(s, chemistry)
% FLUX_NON_EQUILIBRIUM_CHEMISTRY - Compute fluxes for non-equilibrium chemistry transport.
%
%   Handles the transport of effective thermodynamic variables (gamma_star,
%   cv_star) that track the non-equilibrium chemical state. In chemical
%   equilibrium mode the equilibrium values are assigned directly and zero
%   fluxes are returned. In non-equilibrium mode the function computes:
%     1. Advection of gamma_star and cv_star by the local velocity field.
%     2. Relaxation towards chemical equilibrium using an analytic
%        integration of the relaxation ODE to avoid stiffness.
%     3. Optional viscous-like diffusion to simulate species diffusion.
%
% Syntax:
%   s = FLUX_NON_EQUILIBRIUM_CHEMISTRY(s, chemistry)
%
% Inputs:
%   s  - Structure containing flow state data including:
%                 .chemistry.is_chemistry_enabled - Logical; true to enable chemistry.
%                 .chemistry.chemical_equilibrium - Logical; true for equilibrium mode.
%                 .var.gamma_star, .var.cv_star   - Current effective thermo variables.
%                 .var.gamma_star_eq, .var.cv_star_eq - Equilibrium reference values.
%                 .var.rho, .var.rho_u, .var.rho_v - Density and momentum fields.
%                 .time_integration.dt             - Current time step.
%                 .freestream.Re                   - Freestream Reynolds number.
%                 .mesh.x_Ext, .mesh.y_Ext        - Extended grid coordinates.
%   chemistry - Structure containing non-equilibrium model coefficients:
%                 .neq.gamma, .neq.cv   - Arrhenius fit coefficients.
%                 .neq.m_gamma, .neq.m_cv - Pressure exponents.
%
% Outputs:
%   s - Updated structure with:
%                .flux.gamma_star  - Rate of change of gamma_star (interior).
%                .flux.cv_star     - Rate of change of cv_star (interior).
%                .var.gamma_star   - Assigned to equilibrium if in eq. mode.
%                .var.cv_star      - Assigned to equilibrium if in eq. mode.
%
% Notes:
%   - The relaxation term uses an analytic exponential integrator:
%       flux += (eq - current) / dt * (1 - exp(-dt / tau))
%     to avoid numerical stiffness when tau << dt.
%   - The viscous diffusion term uses a standard second-order Laplacian
%     scaled by 1 / (2 * Re).
%
% See also: DERIVATIVE_EXT
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    if s.chemistry.is_chemistry_enabled
        if s.chemistry.chemical_equilibrium
            %% Chemical equilibrium mode
            s.var.gamma_star = s.var.gamma_star_eq;
            s.var.cv_star = s.var.cv_star_eq;
            s.flux.gamma_star = zeros(size(s.flux.rho));
            s.flux.cv_star = zeros(size(s.flux.rho));
        else
            %% Advection equations
            g = s.var.gamma_star;
            cv = s.var.cv_star;
            [dgamma_dx, dgamma_dy] = DERIVATIVE_EXT(g, s);
            [dcv_dx, dcv_dy] = DERIVATIVE_EXT(cv, s);
            u = s.var.rho_u(2:end-1,2:end-1) ./ s.var.rho(2:end-1,2:end-1);
            v = s.var.rho_v(2:end-1,2:end-1) ./ s.var.rho(2:end-1,2:end-1);
            s.flux.gamma_star = -u .* dgamma_dx - v .* dgamma_dy;
            s.flux.cv_star = -u .* dcv_dx - v .* dcv_dy;

            %% Relaxation towards chemical equilibrium
            s = UPDATE_TAU_CHEMISTRY(s, chemistry);

            % Analytic exponential integration to avoid stiffness
            s.flux.gamma_star = s.flux.gamma_star + ...
                (s.var.gamma_star_eq(2:end-1,2:end-1) - s.var.gamma_star(2:end-1,2:end-1)) ./ s.time_integration.dt .* ...
                (1 - exp(-s.time_integration.dt ./ s.tau_gamma_star(2:end-1,2:end-1)));
            s.flux.cv_star = s.flux.cv_star + ...
                (s.var.cv_star_eq(2:end-1,2:end-1) - s.var.cv_star(2:end-1,2:end-1)) ./ s.time_integration.dt .* ...
                (1 - exp(-s.time_integration.dt ./ s.tau_cv_star(2:end-1,2:end-1)));

            %% Viscous diffusion term (species diffusion analogue)
            centroids_bt_distance = sqrt((s.mesh.y_Ext(2:end-1,2:end) - s.mesh.y_Ext(2:end-1,1:end-1)).^2 + ...
                                         (s.mesh.x_Ext(2:end-1,2:end) - s.mesh.x_Ext(2:end-1,1:end-1)).^2);
            centroids_lr_distance = sqrt((s.mesh.y_Ext(2:end,2:end-1) - s.mesh.y_Ext(1:end-1,2:end-1)).^2 + ...
                                         (s.mesh.x_Ext(2:end,2:end-1) - s.mesh.x_Ext(1:end-1,2:end-1)).^2);

            dissip_gamma = (g(1:end-2,2:end-1) - 2 * g(2:end-1,2:end-1) + g(3:end,2:end-1)) ./ centroids_lr_distance(1:end-1,:).^2 + ...
                           (g(2:end-1,1:end-2) - 2 * g(2:end-1,2:end-1) + g(2:end-1,3:end)) ./ centroids_bt_distance(:,1:end-1).^2;
            dissip_cv = (cv(1:end-2,2:end-1) - 2 * cv(2:end-1,2:end-1) + cv(3:end,2:end-1)) ./ centroids_lr_distance(1:end-1,:).^2 + ...
                        (cv(2:end-1,1:end-2) - 2 * cv(2:end-1,2:end-1) + cv(2:end-1,3:end)) ./ centroids_bt_distance(:,1:end-1).^2;

            s.flux.gamma_star = s.flux.gamma_star + 1 / s.freestream.Re * dissip_gamma / 2;
            s.flux.cv_star = s.flux.cv_star + 1 / s.freestream.Re * dissip_cv / 2;
        end
    end
end

function s = UPDATE_TAU_CHEMISTRY(s, chemistry)
% UPDATE_TAU_CHEMISTRY - Update relaxation times for non-equilibrium chemistry.
%
%   Computes the characteristic relaxation time scales tau_gamma_star and
%   tau_cv_star from Arrhenius-type fits to Cantera equilibrium data. Two
%   regression models are supported:
%     Linear:    ln(tau * p^m) = ln(T) + a0 + a1 / T
%     Quadratic: ln(tau * p^m) = ln(T) + a0 + a1 / T + a2 / T^2
%
% Syntax:
%   s = UPDATE_TAU_CHEMISTRY(s, chemistry)
%
% Inputs:
%   s  - Structure containing T, p, freestream scaling factors, and
%               optionally .non_equilibrium_model ('linear' or 'quadratic').
%   chemistry - Structure with .neq sub-fields holding fit coefficients.
%
% Outputs:
%   s - Updated with .tau_gamma_star and .tau_cv_star (non-dimensional).
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    if s.chemistry.chemical_equilibrium
        %% Chemical equilibrium: infinite relaxation times
        s.tau_gamma_star = Inf;
        s.tau_cv_star = Inf;
    else
        %% Extract dimensional temperature and pressure
        T = s.var.T;
        p = s.var.p * s.freestream.velocity_factor^2 * s.freestream.rho_factor;

        % Determine model type (default to linear for backward compatibility)
        if isfield(s.chemistry, 'non_equilibrium_model')
            model_type = s.chemistry.non_equilibrium_model;
        else
            model_type = 'linear';
        end

        %% Compute relaxation times from Arrhenius fits
        if strcmpi(model_type, 'linear')
            % Linear model: ln(tau*p^m) = ln(T) + a0 + a1/T
            s.tau_gamma_star = exp(chemistry.neq.gamma.linear.a0 + ...
                                         chemistry.neq.gamma.linear.a1 ./ T + ...
                                         log(T)) ./ p.^chemistry.neq.m_gamma;

            s.tau_cv_star = exp(chemistry.neq.cv.linear.a0 + ...
                                      chemistry.neq.cv.linear.a1 ./ T + ...
                                      log(T)) ./ p.^chemistry.neq.m_cv;
        else
            % Quadratic model: ln(tau*p^m) = ln(T) + a0 + a1/T + a2/T^2
            s.tau_gamma_star = exp(chemistry.neq.gamma.quadratic.a0 + ...
                                         chemistry.neq.gamma.quadratic.a1 ./ T + ...
                                         chemistry.neq.gamma.quadratic.a2 ./ T^2 + ...
                                         log(T)) ./ p.^chemistry.neq.m_gamma;

            s.tau_cv_star = exp(chemistry.neq.cv.quadratic.a0 + ...
                                      chemistry.neq.cv.quadratic.a1 ./ T + ...
                                      chemistry.neq.cv.quadratic.a2 ./ T^2 + ...
                                      log(T)) ./ p.^chemistry.neq.m_cv;
        end

        %% Non-dimensionalise relaxation times
        s.tau_gamma_star = s.tau_gamma_star * s.freestream.U / s.freestream.L;
        s.tau_cv_star = s.tau_cv_star * s.freestream.U / s.freestream.L;
    end
end
