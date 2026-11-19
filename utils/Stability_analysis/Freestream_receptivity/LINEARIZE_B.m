function B = LINEARIZE_B(s, perturbation, chemistry)
% LINEARIZE_B  Compute the linearized upstream boundary operator B via finite differences.
%
%   B = LINEARIZE_B(s, perturbation, chemistry)
%
%   Constructs the Jacobian matrix B that maps upstream (freestream)
%   perturbations in [rho, rho*u, rho*v, rho*E, (rho*w)] (PRESCRIBE_UPSTREAM_STATE)
%   to flux perturbations in
%   the downstream domain. The linearization uses column-wise finite
%   differences with a 3-color stencil to account for the tridiagonal
%   coupling along the shock front.
%
%   Inputs:
%       s            - Solution structure with base flow, mesh, and parameters
%       perturbation - Finite-difference step size for Jacobian approximation
%       chemistry    - Chemistry model structure
%
%   Outputs:
%       B - Sparse Jacobian matrix mapping upstream state to downstream fluxes.
%           Size is (nv*Nx*Ny x nv*Nx) or (nv*Nx*Ny+Nx x nv*Nx) when shock
%           perturbation DoFs are included, nv = numel(PERTURBATION_FIELDS(s)).
%
%   Notes:
%       - Uses a 3-color finite-difference scheme to reduce the number of
%         nonlinear evaluations from O(Nx) to O(1) per variable.
%       - An upstream point couples to the post-shock cells of its own and
%         its two neighbouring shock points (the row below the shock).
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Freestream Module

    s.freestream.disturbance.amplitude = [0, 0, 0, 0];
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry);
    [s, upstream_fields] = PRESCRIBE_UPSTREAM_STATE(s);
    s_base = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s, chemistry);

    fields = PERTURBATION_FIELDS(s);
    nv = numel(fields);
    Nchi = s.mesh.Nchi;
    N = Nchi * s.mesh.Neta;
    n_rows = nv * N + s.stability_analysis.perturb_shock * Nchi;
    cells_below = (1:Nchi)' + (s.shock.cell_indices(:, 1) - 2) * Nchi;   % row below the shock

    triplets = cell(0, 3);
    for k = 1:nv
        for colour = 1:3
            points = (colour:3:Nchi)';
            s_pert = s;
            s_pert.freestream.(upstream_fields(k))(points) = ...
                s_pert.freestream.(upstream_fields(k))(points) + perturbation;
            response = FLUX_RESPONSE(s_pert, s_base, fields, perturbation, chemistry);

            % Responses at the shock points next to each perturbed point
            for offset = -1:1
                source = points(points + offset >= 1 & points + offset <= Nchi);
                target = source + offset;
                cols = (k - 1) * Nchi + source;
                for o = 1:nv
                    triplets(end+1, :) = {(o - 1) * N + cells_below(target), cols, ...
                                          response.flux{o}(cells_below(target))}; %#ok<AGROW>
                end
                if s.stability_analysis.perturb_shock
                    triplets(end+1, :) = {nv * N + target, cols, response.shock_speed(target)}; %#ok<AGROW>
                end
            end
        end
    end

    B = sparse(vertcat(triplets{:, 1}), vertcat(triplets{:, 2}), vertcat(triplets{:, 3}), ...
               n_rows, nv * Nchi);
end
