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
%       norms                - 4-element logical vector [pressure, u-momentum, v-momentum, entropy]
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
%       norms                - 4-element logical vector
%       T                    - Time horizon
%       scaling_non_temporal - If true, use advection time T_ref for scaling
%
%   Outputs:
%       M_infty - Sparse diagonal matrix of size (4*Nx x 4*Nx)

    %% Extract freestream quantities
    Nx = s.mesh.Nchi;
    Ny = s.mesh.Neta;
    gamma_star = s.freestream.gamma_star * ones(Nx, 1);
    p = s.freestream.p_0 * ones(Nx, 1);
    rho = s.freestream.rho_0 * ones(Nx, 1);
    u = s.freestream.rho_u_0 / s.freestream.rho_0 * ones(Nx, 1);
    v = s.freestream.rho_v_0 / s.freestream.rho_0 * ones(Nx, 1);
    a = s.freestream.a_0 * ones(Nx, 1);
    volume = s.mesh.volume;
    flow_cells = s.shock.flow_cells;
    indices = reshape(1:Nx, Nx, 1);

    %% Compute shock arc lengths and post-shock volume
    shock_arc_length = zeros(Nx, 1);
    total_arc_length = 0;
    total_post_shock_volume = 0;
    volume_advect = zeros(Nx, 1);

    for i = 1:Nx
        shock_arc_length(i) = s.mesh.bt_area(i, s.shock.cell_indices(i, 1)) * s.mesh.bt_y_normal(i, s.shock.cell_indices(i, 1));
        total_arc_length = total_arc_length + shock_arc_length(i);
    end
    for i = 1:Nx
        for j = 1:Ny
            total_post_shock_volume = total_post_shock_volume + volume(i, j) * flow_cells(i, j);
        end
        volume_advect(i, 1) = total_post_shock_volume;
    end

    %% Determine scaling factor
    if scaling_non_temporal
        T_ref = GET_T_REF(s);
        scaling = T_ref;
    else
        scaling = T;
    end

    %% Initialize storage for sparse triplets
    M = zeros(4, Nx);
    M_m = zeros(4, Nx);
    M_n = zeros(4, Nx);

    %% Pressure norm contribution
    if norms(1) == true
        M(1, :) = rho .* a.^2 ./ (2 .* (gamma_star .* p).^2) .* shock_arc_length * scaling;
        M_n(1, :) = 0*Nx + indices;
        M_m(1, :) = 0*Nx + indices;
    end

    %% X-momentum norm contribution
    if norms(2) == true
        M(2, :, :) = rho/2 .* shock_arc_length * scaling;
        M_n(2, :, :) = 1*Nx + indices;
        M_m(2, :, :) = 1*Nx + indices;
    end

    %% Y-momentum norm contribution
    if norms(3) == true
        M(3, :, :) = rho/2 .* shock_arc_length * scaling;
        M_n(3, :, :) = 2*Nx + indices;
        M_m(3, :, :) = 2*Nx + indices;
    end

    %% Entropy norm contribution
    if norms(4) == true
        M(4, :, :) = (gamma_star - 1) .* p ./ (2 * gamma_star) .* shock_arc_length * scaling;
        M_n(4, :, :) = 3*Nx + indices;
        M_m(4, :, :) = 3*Nx + indices;
    end

    %% Assemble sparse matrix from triplets
    values = M(:);
    vector_indices_m = M_m(:);
    vector_indices_n = M_n(:);

    % Keep only non-zero index entries
    mask = (vector_indices_m ~= 0) & (vector_indices_n ~= 0);
    values = values(mask);
    vector_indices_m = vector_indices_m(mask);
    vector_indices_n = vector_indices_n(mask);

    M_infty = sparse(vector_indices_m, vector_indices_n, values, 4*Nx, 4*Nx, size(vector_indices_m, 1));
end
