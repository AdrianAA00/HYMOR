function s = SETUP_AZIMUTHAL_MODE(s, chemistry)
% SETUP_AZIMUTHAL_MODE  Prepare the linearization of an azimuthal mode.
%
%   s = SETUP_AZIMUTHAL_MODE(s, chemistry)
%
%   If s.stability_analysis.azimuthal_wavenumber = m is set (3D-axisymmetric
%   only), adds the azimuthal momentum rho_w (zero in the base flow) to the
%   state and stores in s.azimuthal what the azimuthal-mode terms need:
%       .m, .odd          - azimuthal wavenumber and whether it is odd
%       .base             - base-flow fields after the boundary update
%                           (u, v on the extended grid; p, T and the hoop
%                           stress tau_tt in the cells; the full .var struct
%                           and the shock)
%       .axis_boundaries  - boundaries lying on the symmetry axis (r = 0)
%   The linear state then reads [rho; rho_u; rho_v; rho_E; rho_w; (shock)]
%   (PERTURBATION_FIELDS). Otherwise any previous set-up is removed and the
%   regular axisymmetric (or 2D) operator is built.
%
%   Inputs:
%       s         - Solution structure (base flow) in linearization mode.
%       chemistry - Chemistry model structure.
%
%   Outputs:
%       s - Solution structure with s.azimuthal and s.var.rho_w, or without
%           them (CLEAR_AZIMUTHAL_MODE).
%
% Part of: Hypersonics Stability MATLAB Solver - Stability Analysis / Linear Stability Module

    s = CLEAR_AZIMUTHAL_MODE(s);
    if ~isfield(s.stability_analysis, 'azimuthal_wavenumber') ...
            || isempty(s.stability_analysis.azimuthal_wavenumber)
        return
    end

    m = s.stability_analysis.azimuthal_wavenumber;
    assert(s.PDE_dimension == "3D-axisymmetric", ...
        'stability_analysis.azimuthal_wavenumber requires PDE_dimension = "3D-axisymmetric"');
    assert(isscalar(m) && m >= 0 && m == round(m), ...
        'stability_analysis.azimuthal_wavenumber must be a non-negative integer');
    assert(all(s.mesh.y_corner >= -AXIS_TOLERANCE(s), 'all'), 'Azimuthal modes need the domain at y >= 0 (y = r)');

    %% Base flow after the same state update as the linearization
    s_base = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry), chemistry);
    base.var = s_base.var;
    base.u = s_base.var.rho_u ./ s_base.var.rho;
    base.v = s_base.var.rho_v ./ s_base.var.rho;
    base.p = s_base.var.p(2:end-1, 2:end-1);
    base.T = s_base.var.T(2:end-1, 2:end-1);
    [gu, gv] = FACE_GRADIENTS(s.mesh, base.u, base.v);
    base.tau_tt = HOOP_STRESS(s_base.var.mu_star(2:end-1, 2:end-1) / s.freestream.Re, gu, gv, ...
                              base.v(2:end-1, 2:end-1) ./ s.mesh.y);
    if s.shock.enabled
        base.shock_x = s_base.shock.points_x;
        base.shock_y = s_base.shock.points_y;
        base.shock_beta = s_base.shock.beta;
        base_spline = COMPUTE_BETA(s_base);   % before the symmetry override of beta(1)
        base.shock_beta_spline_1 = base_spline.shock.beta(1);
    end

    s.azimuthal.m = m;
    s.azimuthal.odd = mod(m, 2) == 1;
    s.azimuthal.base = base;
    s.azimuthal.axis_boundaries = AXIS_BOUNDARIES(s);
    s.var.rho_w = zeros(size(s.var.rho));
end

%% ========================================================================
%  Local Helpers
%  ========================================================================
function names = AXIS_BOUNDARIES(s)
% AXIS_BOUNDARIES  Symmetry boundaries whose faces all lie on r = 0.
    y = s.mesh.y_corner;
    tol = AXIS_TOLERANCE(s);
    lines = struct('chi0', y(1, :), 'chi1', y(end, :), 'eta0', y(:, 1), 'eta1', y(:, end));
    names = strings(1, 0);
    for b = ["chi0", "chi1", "eta0", "eta1"]
        if strcmp(s.boundary_conditions.("boundary_" + b).name, 'symmetry') && all(abs(lines.(b)) <= tol)
            names(end+1) = b; %#ok<AGROW>
        end
    end
end

function tol = AXIS_TOLERANCE(s)
% AXIS_TOLERANCE  Distance from r = 0 still counted as on the axis.
%   Shock-fitted meshes do not place the axis exactly at y = 0: the shock
%   spline leaves the axis end of the shock at y ~ 1e-8 R and the wall
%   corner at -1e-16. 1e-6 of the domain size is far above that and far
%   below any distance between the axis and a non-axis boundary.
    tol = 1e-6 * max(abs(s.mesh.y_corner), [], 'all');
end
