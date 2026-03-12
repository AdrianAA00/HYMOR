function CHECK_LINEARIZATION_L(A, s, chemistry)
% CHECK_LINEARIZATION_L  Verify the linearized Jacobian A against finite differences.
%
%   CHECK_LINEARIZATION_L(A, s, chemistry) validates the Jacobian
%   matrix A by comparing A*q (the linearized response) against a
%   finite-difference approximation obtained from the nonlinear operator:
%       dq = [f(x + Dx) - f(x)] / Dx
%   Random perturbations are applied to all conservative variables
%   (rho, rho_u, rho_v, rho_E) and optionally to the shock position.
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
%         position is also perturbed in the normal direction.
%       - The nonlinear operator is evaluated via
%         NON_LINEAR_DYNAMICS_NO_DISCONTINUITY.
%       - A small value of the reported error (e.g., < 1e-3) indicates
%         that the linearization is consistent with the nonlinear operator.
%
% Part of: Hypersonics Stability MATLAB Solver - Stability Analysis / Linear Stability Module

    %% Initialize perturbation and response vectors
    if s.stability_analysis.perturb_shock
        q = zeros(4 * s.mesh.Nchi * s.mesh.Neta + s.mesh.Nchi, 1);
        dq = zeros(4 * s.mesh.Nchi * s.mesh.Neta + s.mesh.Nchi, 1);
    else
        q = zeros(4 * s.mesh.Nchi * s.mesh.Neta, 1);
        dq = zeros(4 * s.mesh.Nchi * s.mesh.Neta, 1);
    end

    solution_perturbed = s;

    %% Perturb rho
    matrix_perturbation = randn(s.mesh.Nchi, s.mesh.Neta);
    if s.shock.enabled
        solution_perturbed.var.rho(2:end-1, 2:end-1) = s.var.rho(2:end-1, 2:end-1) + s.shock.flow_cells .* matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        q(1:s.mesh.Nchi*s.mesh.Neta, 1) = reshape(s.shock.flow_cells .* matrix_perturbation .* s.stability_analysis.perturbation_magnitude, [s.mesh.Nchi*s.mesh.Neta, 1]);
    else
        solution_perturbed.var.rho(2:end-1, 2:end-1) = s.var.rho(2:end-1, 2:end-1) + matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        q(1:s.mesh.Nchi*s.mesh.Neta, 1) = reshape(matrix_perturbation * s.stability_analysis.perturbation_magnitude, [s.mesh.Nchi*s.mesh.Neta, 1]);
    end

    %% Perturb rho_u
    matrix_perturbation = randn(s.mesh.Nchi, s.mesh.Neta);
    if s.shock.enabled
        solution_perturbed.var.rho_u(2:end-1, 2:end-1) = s.var.rho_u(2:end-1, 2:end-1) + s.shock.flow_cells .* matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        q(s.mesh.Nchi*s.mesh.Neta+1:2*s.mesh.Nchi*s.mesh.Neta, 1) = reshape(matrix_perturbation * s.stability_analysis.perturbation_magnitude, [s.mesh.Nchi*s.mesh.Neta, 1]);
    else
        solution_perturbed.var.rho_u(2:end-1, 2:end-1) = s.var.rho_u(2:end-1, 2:end-1) + matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        q(s.mesh.Nchi*s.mesh.Neta+1:2*s.mesh.Nchi*s.mesh.Neta, 1) = reshape(matrix_perturbation * s.stability_analysis.perturbation_magnitude, [s.mesh.Nchi*s.mesh.Neta, 1]);
    end

    %% Perturb rho_v
    matrix_perturbation = randn(s.mesh.Nchi, s.mesh.Neta);
    if s.shock.enabled
        solution_perturbed.var.rho_v(2:end-1, 2:end-1) = s.var.rho_v(2:end-1, 2:end-1) + s.shock.flow_cells .* matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        q(2*s.mesh.Nchi*s.mesh.Neta+1:3*s.mesh.Nchi*s.mesh.Neta, 1) = reshape(matrix_perturbation * s.stability_analysis.perturbation_magnitude, [s.mesh.Nchi*s.mesh.Neta, 1]);
    else
        solution_perturbed.var.rho_v(2:end-1, 2:end-1) = s.var.rho_v(2:end-1, 2:end-1) + matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        q(2*s.mesh.Nchi*s.mesh.Neta+1:3*s.mesh.Nchi*s.mesh.Neta, 1) = reshape(matrix_perturbation * s.stability_analysis.perturbation_magnitude, [s.mesh.Nchi*s.mesh.Neta, 1]);
    end

    %% Perturb rho_E
    matrix_perturbation = randn(s.mesh.Nchi, s.mesh.Neta);
    if s.shock.enabled && s.shock.feedback
        solution_perturbed.var.rho_E(2:end-1, 2:end-1) = s.var.rho_E(2:end-1, 2:end-1) + s.shock.flow_cells .* matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        q(3*s.mesh.Nchi*s.mesh.Neta+1:4*s.mesh.Nchi*s.mesh.Neta, 1) = reshape(matrix_perturbation * s.stability_analysis.perturbation_magnitude, [s.mesh.Nchi*s.mesh.Neta, 1]);
    elseif s.shock.enabled && ~s.shock.feedback
        solution_perturbed.var.rho_E(2:end-1, 2:end-1) = s.var.rho_E(2:end-1, 2:end-1) + s.shock.flow_cells .* matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        q(3*s.mesh.Nchi*s.mesh.Neta+1:4*s.mesh.Nchi*s.mesh.Neta, 1) = reshape(matrix_perturbation * s.stability_analysis.perturbation_magnitude, [s.mesh.Nchi*s.mesh.Neta, 1]);
    else
        solution_perturbed.var.rho_E(2:end-1, 2:end-1) = s.var.rho_E(2:end-1, 2:end-1) + matrix_perturbation * s.stability_analysis.perturbation_magnitude;
        q(3*s.mesh.Nchi*s.mesh.Neta+1:4*s.mesh.Nchi*s.mesh.Neta, 1) = reshape(matrix_perturbation * s.stability_analysis.perturbation_magnitude, [s.mesh.Nchi*s.mesh.Neta, 1]);
    end

    %% Perturb shock position (if enabled)
    if s.stability_analysis.perturb_shock
        vector_perturbation = randn(s.mesh.Nchi, 1);
        ang = s.shock.beta - pi/2 - atan2(s.freestream.rho_v_0, s.freestream.rho_u_0);

        solution_perturbed.shock.points_x = solution_perturbed.shock.points_x + s.stability_analysis.perturbation_magnitude * cos(ang) .* vector_perturbation;
        solution_perturbed.shock.points_y = solution_perturbed.shock.points_y + s.stability_analysis.perturbation_magnitude * sin(ang) .* vector_perturbation;
        q(4*s.mesh.Nchi*s.mesh.Neta+1:4*s.mesh.Nchi*s.mesh.Neta + s.mesh.Nchi, 1) = s.stability_analysis.perturbation_magnitude * vector_perturbation;
    end

    %% Evaluate nonlinear operator on base and perturbed solutions
    s = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s, chemistry);
    solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);

    %% Construct finite-difference response vector dq
    if s.shock.enabled
        dq(0*s.mesh.Nchi*s.mesh.Neta+1:1*s.mesh.Nchi*s.mesh.Neta, 1) = reshape((solution_perturbed.flux.rho - s.flux.rho) .* s.shock.flow_cells, [s.mesh.Nchi*s.mesh.Neta, 1]);
        dq(1*s.mesh.Nchi*s.mesh.Neta+1:2*s.mesh.Nchi*s.mesh.Neta, 1) = reshape((solution_perturbed.flux.rho_u - s.flux.rho_u) .* s.shock.flow_cells, [s.mesh.Nchi*s.mesh.Neta, 1]);
        dq(2*s.mesh.Nchi*s.mesh.Neta+1:3*s.mesh.Nchi*s.mesh.Neta, 1) = reshape((solution_perturbed.flux.rho_v - s.flux.rho_v) .* s.shock.flow_cells, [s.mesh.Nchi*s.mesh.Neta, 1]);
        dq(3*s.mesh.Nchi*s.mesh.Neta+1:4*s.mesh.Nchi*s.mesh.Neta, 1) = reshape((solution_perturbed.flux.rho_E - s.flux.rho_E) .* s.shock.flow_cells, [s.mesh.Nchi*s.mesh.Neta, 1]);
    else
        dq(0*s.mesh.Nchi*s.mesh.Neta+1:1*s.mesh.Nchi*s.mesh.Neta, 1) = reshape((solution_perturbed.flux.rho - s.flux.rho), [s.mesh.Nchi*s.mesh.Neta, 1]);
        dq(1*s.mesh.Nchi*s.mesh.Neta+1:2*s.mesh.Nchi*s.mesh.Neta, 1) = reshape((solution_perturbed.flux.rho_u - s.flux.rho_u), [s.mesh.Nchi*s.mesh.Neta, 1]);
        dq(2*s.mesh.Nchi*s.mesh.Neta+1:3*s.mesh.Nchi*s.mesh.Neta, 1) = reshape((solution_perturbed.flux.rho_v - s.flux.rho_v), [s.mesh.Nchi*s.mesh.Neta, 1]);
        dq(3*s.mesh.Nchi*s.mesh.Neta+1:4*s.mesh.Nchi*s.mesh.Neta, 1) = reshape((solution_perturbed.flux.rho_E - s.flux.rho_E), [s.mesh.Nchi*s.mesh.Neta, 1]);        
    end
    
    if s.stability_analysis.perturb_shock
        dq(4*s.mesh.Nchi*s.mesh.Neta+1:4*s.mesh.Nchi*s.mesh.Neta + s.mesh.Nchi, 1) = (solution_perturbed.shock.relative_increase_velocity - s.shock.relative_increase_velocity);
    end

    %% Compare linearized response A*q against finite-difference dq
    dq_lin = A * q;

    temp = (dq_lin - dq) ./ (abs(dq_lin) + 1e-8);
    error = norm((dq_lin - dq) ./ (abs(dq_lin) + 1e-8)) / sqrt(9 * 4 * 4 * s.mesh.Nchi * s.mesh.Neta);
    fprintf("Error of linearization: %d . \n", error)

end
