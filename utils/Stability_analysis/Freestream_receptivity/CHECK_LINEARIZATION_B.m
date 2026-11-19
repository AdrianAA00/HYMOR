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

    [s, upstream_fields] = PRESCRIBE_UPSTREAM_STATE(s);
    fields = PERTURBATION_FIELDS(s);
    nv = numel(fields);
    Nchi = s.mesh.Nchi;
    N = Nchi * s.mesh.Neta;
    epsilon = s.stability_analysis.perturbation_magnitude;

    %% Random perturbation of every upstream variable
    q = zeros(nv * Nchi, 1);
    solution_perturbed = s;
    for k = 1:nv
        perturbation = randn(Nchi, 1) * epsilon;
        solution_perturbed.freestream.(upstream_fields(k)) = s.freestream.(upstream_fields(k)) + perturbation;
        q((k-1)*Nchi+1:k*Nchi) = perturbation;
    end

    %% Nonlinear response
    s = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s, chemistry);
    solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry);
    dq = zeros(size(B, 1), 1);
    for k = 1:nv
        response = (solution_perturbed.flux.(fields(k)) - s.flux.(fields(k))) .* s.shock.flow_cells;
        dq((k-1)*N+1:k*N) = response(:);
    end
    if s.stability_analysis.perturb_shock
        dq(nv*N+1:end) = solution_perturbed.shock.relative_increase_velocity - s.shock.relative_increase_velocity;
    end

    %% Relative error
    dq_lin = B * q;
    error = norm((dq_lin - dq) ./ (abs(dq_lin) + 1e-8)) / sqrt(3 * nv * nv * Nchi);
    fprintf("Error of linearization: %d . \n", error)
end
