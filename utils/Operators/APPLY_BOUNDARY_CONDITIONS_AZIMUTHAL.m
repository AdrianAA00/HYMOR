function s = APPLY_BOUNDARY_CONDITIONS_AZIMUTHAL(s)
% APPLY_BOUNDARY_CONDITIONS_AZIMUTHAL  Ghost cells of an azimuthal mode.
%
%   s = APPLY_BOUNDARY_CONDITIONS_AZIMUTHAL(s)
%
%   Completes APPLY_BOUNDARY_CONDITIONS for a perturbation with azimuthal
%   wavenumber m (see FLUX_AZIMUTHAL_MODE):
%     - Sets the ghost cells of the azimuthal momentum rho_w on every
%       boundary from its condition (no-slip: w = 0, inflow: w = 0, slip or
%       symmetry plane: dw/dn = 0, outflow: extrapolation, periodic: copy).
%     - On the symmetry axis (r = 0) the ghost cell is the interior cell
%       rotated by pi, so a mode picks up (-1)^m. The ghost perturbation of
%       every variable is the m = 0 mirror image (the regular symmetry
%       condition) times (-1)^m; rho_w, odd in the m = 0 mirror, is
%       -(-1)^m times the interior value.
%   The shock boundary is set by UPDATE_SHOCK_JUMP_PROPERTIES.
%
%   Inputs:
%       s - Solution structure in linearization mode (s.azimuthal set),
%           after APPLY_BOUNDARY_CONDITIONS.
%
%   Outputs:
%       s - Solution structure with azimuthal-mode ghost cells.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    for boundary = ["chi0", "chi1", "eta0", "eta1"]
        [ghost, inner, inner2, opposite] = BOUNDARY_LINES(boundary, size(s.var.rho));
        switch s.boundary_conditions.("boundary_" + boundary).name
            case {'no_slip_adiabatic', 'no_slip_isothermal', 'lid_driven_cavity'}
                s.var.rho_w(ghost{:}) = -s.var.rho_w(inner{:});
            case {'inflow_subsonic', 'inflow_supersonic'}
                s.var.rho_w(ghost{:}) = 0;
            case {'outflow_subsonic', 'outflow_supersonic', 'outflow_NRCBC'}
                s.var.rho_w(ghost{:}) = 2 * s.var.rho_w(inner{:}) - s.var.rho_w(inner2{:});
            case 'outflow_supersonic_1st'
                s.var.rho_w(ghost{:}) = s.var.rho_w(inner{:});
            case 'periodic'
                s.var.rho_w(ghost{:}) = s.var.rho_w(opposite{:});
            case 'symmetry'
                if any(s.azimuthal.axis_boundaries == boundary)
                    s.var.rho_w(ghost{:}) = -(-1)^s.azimuthal.m * s.var.rho_w(inner{:});
                    if s.azimuthal.odd
                        s = FLIP_GHOST_PERTURBATION(s, ghost);
                    end
                else
                    s.var.rho_w(ghost{:}) = s.var.rho_w(inner{:});
                end
            case 'shock'
                % Set at the shock by UPDATE_SHOCK_JUMP_PROPERTIES
        end
    end
end

%% ========================================================================
%  Local Helpers
%  ========================================================================
function s = FLIP_GHOST_PERTURBATION(s, ghost)
% FLIP_GHOST_PERTURBATION  Reverse the sign of the ghost perturbation of the
%   fields set by the symmetry condition: q_ghost -> 2 q0_ghost - q_ghost.
    base = s.azimuthal.base.var;
    for f = ["rho", "rho_u", "rho_v", "rho_E", "p", "T", "gamma_star", "cv_star"]
        s.var.(f)(ghost{:}) = 2 * base.(f)(ghost{:}) - s.var.(f)(ghost{:});
    end
end

function [ghost, inner, inner2, opposite] = BOUNDARY_LINES(boundary, n)
% BOUNDARY_LINES  Index sets {rows, cols} of a boundary's ghost line, the two
%   adjacent interior lines and the interior line across a periodic pair.
    switch boundary
        case "chi0"
            ghost = {1, ':'};    inner = {2, ':'};    inner2 = {3, ':'};    opposite = {n(1)-1, ':'};
        case "chi1"
            ghost = {n(1), ':'}; inner = {n(1)-1, ':'}; inner2 = {n(1)-2, ':'}; opposite = {2, ':'};
        case "eta0"
            ghost = {':', 1};    inner = {':', 2};    inner2 = {':', 3};    opposite = {':', n(2)-1};
        case "eta1"
            ghost = {':', n(2)}; inner = {':', n(2)-1}; inner2 = {':', n(2)-2}; opposite = {':', 2};
    end
end
