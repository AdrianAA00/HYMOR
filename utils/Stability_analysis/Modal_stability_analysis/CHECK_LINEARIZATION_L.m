function CHECK_LINEARIZATION_L(A, s, chemistry)
% CHECK_LINEARIZATION_L  Verify the linearized Jacobian A against finite differences.
%
%   CHECK_LINEARIZATION_L(A, s, chemistry) validates the Jacobian
%   matrix A by comparing A*q (the linearized response) against a
%   finite-difference approximation obtained from the nonlinear operator:
%       dq = [f(x + Dx) - f(x)] / Dx
%   Random perturbations are applied to all conservative variables
%   (PERTURBATION_FIELDS: rho, rho_u, rho_v, rho_E and, for an azimuthal
%   mode, rho_w) and optionally to the shock position.
%
%   Inputs:
%       A         - Linearized Jacobian matrix (sparse).
%       s  - Solution structure containing the base flow, grid, and
%                   solver parameters.
%       chemistry - Chemistry model structure for thermodynamic evaluations.
%
%   Outputs:
%       (none) - Prints the relative linearization error to the console.
%
%   Notes:
%       - The perturbation magnitude is controlled by s.stability_analysis.perturbation_magnitude.
%       - When s.stability_analysis.perturb_shock is true, the shock
%         position is also perturbed along the upstream normal (the DoF of
%         LINEARIZE_EQUATIONS_NO_DISCONTINUITY).
%       - The nonlinear operator is evaluated via
%         NON_LINEAR_DYNAMICS_NO_DISCONTINUITY.
%       - A small value of the reported error (e.g., < 1e-3) indicates
%         that the linearization is consistent with the nonlinear operator.
%
% Part of: Hypersonics Stability MATLAB Solver - Stability Analysis / Linear Stability Module

    fields = PERTURBATION_FIELDS(s);
    nv = numel(fields);
    N = s.mesh.Nchi * s.mesh.Neta;
    epsilon = s.stability_analysis.perturbation_magnitude;
    q = zeros(nv * N + s.stability_analysis.perturb_shock * s.mesh.Nchi, 1);
    dq = zeros(size(q));

    %% Random perturbation of every variable (flow cells only behind a shock)
    solution_perturbed = s;
    for k = 1:nv
        perturbation = randn(s.mesh.Nchi, s.mesh.Neta) * epsilon;
        if s.shock.enabled
            perturbation = s.shock.flow_cells .* perturbation;
        end
        solution_perturbed.var.(fields(k))(2:end-1, 2:end-1) = ...
            s.var.(fields(k))(2:end-1, 2:end-1) + perturbation;
        q((k-1)*N+1:k*N) = perturbation(:);
    end
    if s.stability_analysis.perturb_shock
        vector_perturbation = randn(s.mesh.Nchi, 1);
        ang = s.shock.beta + pi/2 - atan2(s.freestream.rho_v_0, s.freestream.rho_u_0);   % upstream normal
        solution_perturbed.shock.points_x = solution_perturbed.shock.points_x + epsilon * cos(ang) .* vector_perturbation;
        solution_perturbed.shock.points_y = solution_perturbed.shock.points_y + epsilon * sin(ang) .* vector_perturbation;
        q(nv*N+1:end) = epsilon * vector_perturbation;
    end

    %% Nonlinear response
    s = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s, chemistry);
    solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
    for k = 1:nv
        response = solution_perturbed.flux.(fields(k)) - s.flux.(fields(k));
        if s.shock.enabled
            response = response .* s.shock.flow_cells;
        end
        dq((k-1)*N+1:k*N) = response(:);
    end
    if s.stability_analysis.perturb_shock
        dq(nv*N+1:end) = solution_perturbed.shock.relative_increase_velocity - s.shock.relative_increase_velocity;
    end

    %% Relative error
    dq_lin = A * q;
    error = norm((dq_lin - dq) ./ (abs(dq_lin) + 1e-8)) / sqrt(9 * nv * nv * N);
    fprintf("Error of linearization: %d . \n", error)
end
