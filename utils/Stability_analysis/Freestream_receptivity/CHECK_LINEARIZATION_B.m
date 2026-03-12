function CHECK_LINEARIZATION_B(B, s, chemistry)
% CHECK_LINEARIZATION_B  Verify the linearized upstream boundary operator B by finite differences.
%
%   CHECK_LINEARIZATION_B(B, s, chemistry)
%
%   Validates the Jacobian matrix B (upstream boundary linearization) by
%   comparing B*q against a finite-difference approximation computed from
%   the full nonlinear dynamics. Prints the relative error norm.
%
%   Inputs:
%       B         - Linearized upstream boundary operator (sparse matrix)
%       s  - Solution structure with base flow, mesh, and parameters:
%                     .mesh.Nchi, .mesh.Neta                 - Grid dimensions
%                     .freestream.rho_0, etc.     - Base upstream state
%                     .stability_analysis.perturbation_magnitude   - Finite-difference step size
%                     .stability_analysis.perturb_shock - Shock DoF flag
%       chemistry - Chemistry model structure
%
%   Outputs:
%       (none) - Prints the linearization error to the console
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Freestream Module

    %% Set uniform upstream base state
    s.freestream.rho_0_p = s.freestream.rho_0 * ones(s.mesh.Nchi, 1);
    s.freestream.rho_u_0_p = s.freestream.rho_u_0 * ones(s.mesh.Nchi, 1);
    s.freestream.rho_v_0_p = s.freestream.rho_v_0 * ones(s.mesh.Nchi, 1);
    s.freestream.rho_E_0_p = s.freestream.rho_E_0 * ones(s.mesh.Nchi, 1);
    s.rho_0_upstream_p = true;

    q = zeros(4*s.mesh.Nchi, 1);
    dq = zeros(4*s.mesh.Nchi, 1);

    %% Compute nonlinear finite-difference perturbation
    solution_perturbed = s;

    % rho perturbation
    p = randn(s.mesh.Nchi, 1);
    solution_perturbed.freestream.rho_0_p = s.freestream.rho_0_p + p .* s.stability_analysis.perturbation_magnitude;
    q(1:s.mesh.Nchi, 1) = p .* s.stability_analysis.perturbation_magnitude;

    % rho_u perturbation
    p = randn(s.mesh.Nchi, 1);
    solution_perturbed.freestream.rho_u_0_p = s.freestream.rho_u_0_p + p .* s.stability_analysis.perturbation_magnitude;
    q(1+s.mesh.Nchi:2*s.mesh.Nchi, 1) = p .* s.stability_analysis.perturbation_magnitude;

    % rho_v perturbation
    p = randn(s.mesh.Nchi, 1);
    solution_perturbed.freestream.rho_v_0_p = s.freestream.rho_v_0_p + p .* s.stability_analysis.perturbation_magnitude;
    q(1+2*s.mesh.Nchi:3*s.mesh.Nchi, 1) = p .* s.stability_analysis.perturbation_magnitude;

    % rho_E perturbation
    p = randn(s.mesh.Nchi, 1);
    solution_perturbed.freestream.rho_E_0_p = s.freestream.rho_E_0_p + p .* s.stability_analysis.perturbation_magnitude;
    q(1+3*s.mesh.Nchi:4*s.mesh.Nchi, 1) = p .* s.stability_analysis.perturbation_magnitude;

    %% Evaluate nonlinear dynamics for base and perturbed states
    s = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s, chemistry);
    solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);

    %% Build finite-difference Jacobian-vector product
    dq(0*s.mesh.Nchi*s.mesh.Neta+1:1*s.mesh.Nchi*s.mesh.Neta, 1) = reshape((solution_perturbed.flux.rho - s.flux.rho) .* s.shock.flow_cells, [s.mesh.Nchi*s.mesh.Neta, 1]);
    dq(1*s.mesh.Nchi*s.mesh.Neta+1:2*s.mesh.Nchi*s.mesh.Neta, 1) = reshape((solution_perturbed.flux.rho_u - s.flux.rho_u) .* s.shock.flow_cells, [s.mesh.Nchi*s.mesh.Neta, 1]);
    dq(2*s.mesh.Nchi*s.mesh.Neta+1:3*s.mesh.Nchi*s.mesh.Neta, 1) = reshape((solution_perturbed.flux.rho_v - s.flux.rho_v) .* s.shock.flow_cells, [s.mesh.Nchi*s.mesh.Neta, 1]);
    dq(3*s.mesh.Nchi*s.mesh.Neta+1:4*s.mesh.Nchi*s.mesh.Neta, 1) = reshape((solution_perturbed.flux.rho_E - s.flux.rho_E) .* s.shock.flow_cells, [s.mesh.Nchi*s.mesh.Neta, 1]);

    if s.stability_analysis.perturb_shock
        dq(4*s.mesh.Nchi*s.mesh.Neta+1:4*s.mesh.Nchi*s.mesh.Neta + s.mesh.Nchi, 1) = (solution_perturbed.shock.relative_increase_velocity - s.shock.relative_increase_velocity);
    end

    %% Compare linearized vs finite-difference result
    dq_lin = B * q;

    %% Compute and report error
    temp = (dq_lin - dq) ./ (abs(dq_lin) + 1e-8);
    [a, b] = max(temp);
    error = norm((dq_lin - dq) ./ (abs(dq_lin) + 1e-8)) / sqrt(3*4*4*s.mesh.Nchi);
    fprintf("Error of linearization: %d . \n", error)
end
