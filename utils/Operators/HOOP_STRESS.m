function tau_tt = HOOP_STRESS(mu, gu, gv, h)
% HOOP_STRESS  Cell value of the axisymmetric hoop stress
%   tau_theta_theta = mu (2 h - 2/3 (div_xy + h)).
%
%   tau_tt = HOOP_STRESS(mu, gu, gv, h)
%
%   div_xy = du/dx + dv/dy is the mean of the four face values. For an
%   axisymmetric field h = v / r; for an azimuthal mode m (cos amplitude)
%   h = (v + m w) / r, kept as one quotient so that it stays regular on the
%   axis for odd m.
%
%   Inputs:
%       mu     - (Nchi x Neta) viscosity in the cells (non-dimensional, / Re).
%       gu, gv - Face gradients of u and v (FACE_GRADIENTS).
%       h      - (Nchi x Neta) (1/r) d_theta(w) + v / r in the cells.
%
%   Outputs:
%       tau_tt - (Nchi x Neta) hoop stress.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    div_xy = CELL_AVERAGE(gu.x_lr + gv.y_lr, gu.x_bt + gv.y_bt);
    tau_tt = mu .* (2 * h - 2/3 * (div_xy + h));
end
