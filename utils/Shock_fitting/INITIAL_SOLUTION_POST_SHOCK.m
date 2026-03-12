function s = INITIAL_SOLUTION_POST_SHOCK(s, chemistry)
% INITIAL_SOLUTION_POST_SHOCK - Initialize post-shock flow field using Rankine-Hugoniot relations.
%
%   s = INITIAL_SOLUTION_POST_SHOCK(s, chemistry)
%
%   Computes the initial post-shock conservative variables (rho, rho_u,
%   rho_v, rho_E) at every streamwise station by applying oblique-shock
%   Rankine-Hugoniot jump conditions. Supports both calorically perfect
%   gas and finite-rate / equilibrium chemistry models.
%
%   Inputs:
%       s  (struct) - Solution structure containing at minimum:
%                            .mesh.Nchi, .mesh.Neta  - Grid dimensions
%                            .freestream.Mach        - Freestream Mach number
%                            .freestream.gamma       - Ratio of specific heats (perfect gas)
%                            .shock.beta             - (Nx x 1) local shock angles [rad]
%                            .shock.cell_indices     - (Nx x 1) shocked-cell column indices
%                            .chemistry.is_chemistry_enabled - (logical) chemistry enabled flag
%                            .chemistry.chemical_equilibrium - (logical) equilibrium chemistry flag
%                            .freestream             - Freestream reference quantities
%       chemistry (struct) - Chemistry model structure; used when
%                            chemistry_state is true. Contains either
%                            equilibrium or frozen sub-models with
%                            SOLVE_RANKINE_HUGONIOT_CHEMISTRY and
%                            eval_gamma_star / eval_cv_star methods.
%
%   Outputs:
%       s  (struct) - Updated s with post-shock fields populated
%                            from the body surface up to the shock location.
%
%   Notes:
%       - For perfect gas, the standard oblique-shock relations are used.
%         If M*sin(beta) <= 1 at a station, that station is treated as
%         unshocked and set to freestream conditions.
%       - For reacting flows, SOLVE_RANKINE_HUGONIOT_CHEMISTRY is called
%         to obtain the post-shock state, then results are non-dimensionalized.
%       - After the main loop, gamma_star and cv_star are evaluated for
%         frozen chemistry cases.
%
%   See also: SOLVE_RANKINE_HUGONIOT_CHEMISTRY, COMPUTE_BETA
%
% Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

    if ~s.chemistry.is_chemistry_enabled
        %% Perfect gas oblique-shock initialization
        M = s.freestream.Mach;

        for i = 1:s.mesh.Nchi
            beta = s.shock.beta(i, 1);

            if (M * sin(beta) > 1)
                ang = beta - atan2(s.freestream.rho_v_0, s.freestream.rho_u_0);
                temp = (s.freestream.gamma + 1) * (M * sin(beta))^2 ...
                     / ((s.freestream.gamma - 1) * (M * sin(beta))^2 + 2);
                s.var.rho(i + 1, 1:s.shock.cell_indices(i, 1) + 2) = temp;

                normal_velocity = sin(ang) / temp;
                tangential_velocity = cos(ang);
                s.var.rho_u(i + 1, 1:s.shock.cell_indices(i, 1) + 2) = ...
                    temp * (normal_velocity * sin(beta) + tangential_velocity * cos(beta));
                s.var.rho_v(i + 1, 1:s.shock.cell_indices(i, 1) + 2) = ...
                    temp * (-normal_velocity * cos(beta) + tangential_velocity * sin(beta));

                % Total enthalpy is preserved across the shock
                p_inf = 1 / s.freestream.gamma / M^2;
                k = (s.var.rho_u(i + 1, 1)^2 + s.var.rho_v(i + 1, 1)^2) / 2 / temp^2;
                p = p_inf * (2 * s.freestream.gamma * (M * sin(beta))^2 - (s.freestream.gamma - 1)) ...
                  / (s.freestream.gamma + 1);
                s.var.rho_E(i + 1, 1:s.shock.cell_indices(i, 1) + 2) = ...
                    p / (s.freestream.gamma - 1) + k * temp;
            else
                % No shock at this station: set to freestream
                s.var.rho(i + 1, 1:s.shock.cell_indices(i, 1) + 2) = 1;
                s.var.rho_u(i + 1, 1:s.shock.cell_indices(i, 1) + 2) = ...
                    s.var.rho_u(s.mesh.Nchi + 1, s.mesh.Neta + 1);
                s.var.rho_v(i + 1, 1:s.shock.cell_indices(i, 1) + 2) = ...
                    s.var.rho_v(s.mesh.Nchi + 1, s.mesh.Neta + 1);
                s.var.rho_E(i + 1, 1:s.shock.cell_indices(i, 1) + 2) = ...
                    (1/2 + 1 / (s.freestream.gamma - 1) / s.freestream.gamma / M^2);
            end
        end
    else
        %% Reacting-gas oblique-shock initialization
        for i = 1:s.mesh.Nchi
            beta = s.shock.beta(i, 1);
            ang = beta - atan2(s.freestream.rho_v_0, s.freestream.rho_u_0);
            w_1 = sin(ang) * s.freestream.U;
            e_1 = s.freestream.e;
            rho_1 = s.freestream.rho;

            if s.chemistry.chemical_equilibrium
                [rho_2, e_2, w_2] = SOLVE_RANKINE_HUGONIOT_CHEMISTRY(chemistry, w_1, rho_1, e_1);
            else
                [rho_2, e_2, w_2] = SOLVE_RANKINE_HUGONIOT_CHEMISTRY(chemistry.frozen, w_1, rho_1, e_1);
            end

            % Non-dimensionalize post-shock state
            rho_2_nd = rho_2 / s.freestream.rho_factor;
            e_2_nd   = e_2   / s.freestream.energy_factor;
            w_2_nd   = w_2   / s.freestream.velocity_factor;

            normal_velocity     = w_2_nd;
            tangential_velocity = cos(ang) * s.freestream.U / s.freestream.velocity_factor;

            s.var.rho(i + 1, 1:s.shock.cell_indices(i, 1) + 2) = rho_2_nd;
            s.var.rho_u(i + 1, 1:s.shock.cell_indices(i, 1) + 2) = ...
                rho_2_nd * (normal_velocity * sin(beta) + tangential_velocity * cos(beta));
            s.var.rho_v(i + 1, 1:s.shock.cell_indices(i, 1) + 2) = ...
                rho_2_nd * (-normal_velocity * cos(beta) + tangential_velocity * sin(beta));

            k = (normal_velocity^2 + tangential_velocity^2) / 2;
            s.var.rho_E(i + 1, 1:s.shock.cell_indices(i, 1) + 2) = ...
                rho_2_nd * e_2_nd + k * rho_2_nd;
        end
    end

    %% Evaluate thermodynamic properties for frozen chemistry
    if s.chemistry.is_chemistry_enabled
        if ~s.chemistry.chemical_equilibrium
            rho = s.var.rho;
            e = (s.var.rho_E - (s.var.rho_u.^2 + s.var.rho_v.^2) ./ rho / 2) ./ rho;
            s.var.gamma_star = chemistry.frozen.eval_gamma_star( ...
                s.freestream.rho_factor * rho, s.freestream.energy_factor * e);
            s.var.cv_star = chemistry.frozen.eval_cv_star( ...
                s.freestream.rho_factor * rho, s.freestream.energy_factor * e) ./ s.freestream.cv;
        end
    end
end
