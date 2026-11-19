function f_c = CELL_AVERAGE(f_lr, f_bt)
% CELL_AVERAGE  Cell value as the mean of its four face values.
%
%   f_c = CELL_AVERAGE(f_lr, f_bt)
%
%   Inputs:
%       f_lr - (Nchi+1 x Neta) values on the left-right faces.
%       f_bt - (Nchi x Neta+1) values on the bottom-top faces.
%
%   Outputs:
%       f_c - (Nchi x Neta) cell values.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    f_c = (f_lr(1:end-1, :) + f_lr(2:end, :) + f_bt(:, 1:end-1) + f_bt(:, 2:end)) / 4;
end
