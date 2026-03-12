function s = INITIAL_SOLUTION_POST_SHOCK_MODIFIED(s, chemistry)
% INITIAL_SOLUTION_POST_SHOCK_MODIFIED - Apply a blunt-body profile correction to the post-shock field.
%
%   s = INITIAL_SOLUTION_POST_SHOCK_MODIFIED(s, chemistry)
%
%   Modifies the initial post-shock s by applying a parabolic
%   density/energy amplification and an exponentially decaying normal-
%   velocity reduction from the shock toward the body. This produces a
%   more physical initial guess for blunt-body flows that accelerates
%   convergence of the flow solver.
%
%   Inputs:
%       s  (struct) - Solution structure (already initialized by
%                            INITIAL_SOLUTION_POST_SHOCK) with fields:
%                            .mesh.Nchi, .mesh.Neta  - Grid dimensions
%                            .shock.cell_indices     - (Nx x 1) shocked-cell indices
%                            .var.rho, .var.rho_u, .var.rho_v, .var.rho_E
%                            .mesh.bt_x_normal, .mesh.bt_y_normal - Shock-normal unit vectors
%       chemistry (struct) - Chemistry model (unused in this function but
%                            kept for interface consistency).
%
%   Outputs:
%       s  (struct) - Modified s with adjusted density, momentum,
%                            and total energy fields.
%
%   Notes:
%       - Density and energy are amplified by a factor of
%         (1 + 0.5 * (1 - (j/(N-1))^2)), which peaks at the shock and
%         vanishes at the body.
%       - Normal velocity is reduced exponentially with decay index = 5.
%       - Kinetic energy is corrected so that total energy remains consistent
%         with the modified momentum field.
%
%   See also: INITIAL_SOLUTION_POST_SHOCK
%
% Part of: Hypersonics Stability MATLAB Solver - Shock Fitting Module

    k0 = ((s.var.rho_u.^2 + s.var.rho_v.^2) ./ (s.var.rho) / 2);
    index = 5;

    for i = 1:s.mesh.Nchi
        N_end = s.shock.cell_indices(i, 1);
        for j = 1:N_end
            normal = s.var.rho_u(i + 1, j + 1) * s.mesh.bt_x_normal(i, 1) ...
                   + s.var.rho_v(i + 1, j + 1) * s.mesh.bt_y_normal(i, 1);

            s.var.rho(i + 1, j + 1) = ...
                s.var.rho(i + 1, j + 1) * (1 + 0.5 * (1 - (j / (N_end - 1))^2));
            s.var.rho_E(i + 1, j + 1) = ...
                s.var.rho_E(i + 1, j + 1) * (1 + 0.5 * (1 - (j / (N_end - 1))^2));
            s.var.rho_u(i + 1, j + 1) = ...
                s.var.rho_u(i + 1, j + 1) - normal * s.mesh.bt_x_normal(i, 1) * exp(-index * j / (N_end - 1));
            s.var.rho_v(i + 1, j + 1) = ...
                s.var.rho_v(i + 1, j + 1) - normal * s.mesh.bt_y_normal(i, 1) * exp(-index * j / (N_end - 1));
        end
    end

    %% Correct total energy for modified kinetic energy
    k1 = ((s.var.rho_u.^2 + s.var.rho_v.^2) ./ (s.var.rho) / 2);
    s.var.rho_E = s.var.rho_E + k1 - k0;
end
