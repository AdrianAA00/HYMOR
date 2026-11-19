function f_lr = FACE_LR(f)
% FACE_LR  Mean of the two cells adjacent to each left-right (chi) face.
%
%   f_lr = FACE_LR(f)
%
%   Inputs:
%       f - Cell field on the extended grid (Nchi+2 x Neta+2).
%
%   Outputs:
%       f_lr - (Nchi+1 x Neta) face values.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    f_lr = (f(1:end-1, 2:end-1) + f(2:end, 2:end-1)) / 2;
end
