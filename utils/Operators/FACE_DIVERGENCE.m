function div = FACE_DIVERGENCE(F_lr, F_bt, mesh)
% FACE_DIVERGENCE  Net outward face integral per unit cell volume.
%
%   div = FACE_DIVERGENCE(F_lr, F_bt, mesh)
%
%   Inputs:
%       F_lr - (Nchi+1 x Neta) flux through the left-right faces, F . n * area,
%              with the face normal pointing towards increasing chi.
%       F_bt - (Nchi x Neta+1) flux through the bottom-top faces, normal
%              towards increasing eta.
%       mesh - Mesh structure (volume).
%
%   Outputs:
%       div - (Nchi x Neta) sum of the outward fluxes divided by the volume.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    div = (F_lr(2:end, :) - F_lr(1:end-1, :) + F_bt(:, 2:end) - F_bt(:, 1:end-1)) ./ mesh.volume;
end
