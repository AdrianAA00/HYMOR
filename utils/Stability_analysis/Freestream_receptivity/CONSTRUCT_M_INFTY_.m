function M_infty_ = CONSTRUCT_M_INFTY_(s, norms, T, w_infty, scaling_non_temporal)
% CONSTRUCT_M_INFTY_  Build the extended initial-condition energy-weight matrix.
%
%   M_infty_ = CONSTRUCT_M_INFTY_(s, norms, T, w_infty, scaling_non_temporal)
%
%   Constructs the block-diagonal energy-weight matrix for the denominator
%   of the transient growth gain in the freestream-coupled system. The
%   downstream block uses M_hat (initial-condition norm); each freestream
%   frequency block uses M_infty (inflow energy norm), with a factor of
%   1/2 applied for non-zero frequencies (since the time-averaged energy
%   of sin^2 or cos^2 is 1/2).
%
%   Structure:
%       M_infty_ = blkdiag(M_hat, M_infty_1, M_infty_2, ..., M_infty_N)
%       where M_infty_k = M_infty for w=0, or M_infty/2 for w != 0.
%
%   Inputs:
%       s             - Solution structure with flow fields and mesh data
%       norms                - Logical vector [pressure, u-momentum, v-momentum, entropy, (w-momentum)]
%       T                    - Time horizon for temporal scaling
%       w_infty              - Vector of freestream disturbance frequencies
%       scaling_non_temporal - If true, use reference time T_ref instead of T for scaling
%
%   Outputs:
%       M_infty_ - Sparse block-diagonal matrix for the extended system
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Freestream Module

    %% Build downstream and freestream weight matrices
    M_hat = CONSTRUCT_M_HAT(s, norms);
    M_infty = CONSTRUCT_M_INFTY(s, norms, T, scaling_non_temporal);

    %% Create frequency-dependent blocks
    N = length(w_infty);
    M_blocks = cell(1, N);
    for k = 1:N
        if w_infty(k) == 0
            M_blocks{k} = M_infty;
        else
            M_blocks{k} = M_infty / 2;
        end
    end

    %% Assemble block-diagonal matrix
    M_infty_ = blkdiag(M_hat, M_blocks{:});
end


function M_infty = CONSTRUCT_M_INFTY(s, norms, T, scaling_non_temporal)
% CONSTRUCT_M_INFTY  Build the freestream inflow energy-weight matrix.
%
%   M_infty = CONSTRUCT_M_INFTY(s, norms, T, scaling_non_temporal)
%
%   Constructs a sparse diagonal matrix that weights the freestream
%   upstream state by the Chu energy norm scaled by the shock arc length
%   and a time or reference-time factor.
%
%   Inputs:
%       s             - Solution structure
%       norms                - Logical vector [pressure, u, v, entropy, (w)]
%       T                    - Time horizon
%       scaling_non_temporal - If true, use advection time T_ref for scaling
%
%   Outputs:
%       M_infty - Sparse diagonal matrix of size (nv*Nx x nv*Nx)

    %% Extract freestream quantities
    Nx = s.mesh.Nchi;
    gamma_star = s.freestream.gamma_star * ones(Nx, 1);
    p = s.freestream.p_0 * ones(Nx, 1);
    rho = s.freestream.rho_0 * ones(Nx, 1);
    a = s.freestream.a_0 * ones(Nx, 1);

    %% Shock arc length projected normal to the freestream
    V_inf_mag = sqrt(s.freestream.rho_u_0^2 + s.freestream.rho_v_0^2);
    u_hat = s.freestream.rho_u_0 / V_inf_mag;
    v_hat = s.freestream.rho_v_0 / V_inf_mag;
    shock_faces = sub2ind(size(s.mesh.bt_area), (1:Nx)', s.shock.cell_indices(:, 1));
    proj = s.mesh.bt_x_normal(shock_faces) * u_hat + s.mesh.bt_y_normal(shock_faces) * v_hat;
    shock_arc_length = s.mesh.bt_area(shock_faces) .* abs(proj);

    %% Time scaling
    if scaling_non_temporal
        scaling = GET_T_REF(s);
    else
        scaling = T;
    end

    %% Chu weights of the upstream state
    W = CHU_NORM_WEIGHTS(rho, a, gamma_star, p, norms, numel(PERTURBATION_FIELDS(s)));
    values = W .* shock_arc_length * scaling;
    M_infty = spdiags(values(:), 0, numel(values), numel(values));
end
