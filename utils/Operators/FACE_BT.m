function f_bt = FACE_BT(f)
% FACE_BT  Mean of the two cells adjacent to each bottom-top (eta) face.
%
%   f_bt = FACE_BT(f)
%
%   Inputs:
%       f - Cell field on the extended grid (Nchi+2 x Neta+2).
%
%   Outputs:
%       f_bt - (Nchi x Neta+1) face values.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    f_bt = (f(2:end-1, 1:end-1) + f(2:end-1, 2:end)) / 2;
end
