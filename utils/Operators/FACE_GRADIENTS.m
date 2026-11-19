function varargout = FACE_GRADIENTS(mesh, varargin)
% FACE_GRADIENTS  Cartesian gradients of cell fields on every cell face.
%
%   [g1, g2, ...] = FACE_GRADIENTS(mesh, f1, f2, ...)
%
%   On each face, the gradient is the solution of
%       grad(f) . d = df_d    (d: centroid-to-centroid vector across the face)
%       grad(f) . t = df_t    (t: corner-to-corner vector along the face)
%   with df_d the difference of the two adjacent cell values and df_t the
%   difference of the corner values (average of the four cells around each
%   corner). The face-normal and tangential differences are consistent on
%   non-orthogonal meshes; gradients of linear fields are exact wherever
%   each corner is the mean of its four surrounding centroids (smoothly
%   varying meshes: second-order error).
%
%   Inputs:
%       mesh   - Mesh structure (x_Ext, y_Ext, x_corner, y_corner).
%       f1,... - Cell fields on the extended grid (Nchi+2 x Neta+2).
%
%   Outputs:
%       g1,... - Structs with fields
%                  .x_lr, .y_lr  (Nchi+1 x Neta)  on the left-right faces
%                  .x_bt, .y_bt  (Nchi x Neta+1)  on the bottom-top faces
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% Inverse of [dx dy; tx ty] on each face family
    lr = INVERSE_2X2( ...
        mesh.x_Ext(2:end, 2:end-1) - mesh.x_Ext(1:end-1, 2:end-1), ...
        mesh.y_Ext(2:end, 2:end-1) - mesh.y_Ext(1:end-1, 2:end-1), ...
        mesh.x_corner(:, 2:end) - mesh.x_corner(:, 1:end-1), ...
        mesh.y_corner(:, 2:end) - mesh.y_corner(:, 1:end-1));
    bt = INVERSE_2X2( ...
        mesh.x_Ext(2:end-1, 2:end) - mesh.x_Ext(2:end-1, 1:end-1), ...
        mesh.y_Ext(2:end-1, 2:end) - mesh.y_Ext(2:end-1, 1:end-1), ...
        mesh.x_corner(2:end, :) - mesh.x_corner(1:end-1, :), ...
        mesh.y_corner(2:end, :) - mesh.y_corner(1:end-1, :));

    %% Gradients of every field
    varargout = cell(1, numel(varargin));
    for k = 1:numel(varargin)
        f = varargin{k};
        f_corner = (f(1:end-1, 1:end-1) + f(2:end, 1:end-1) + f(1:end-1, 2:end) + f(2:end, 2:end)) / 4;
        df_d = f(2:end, 2:end-1) - f(1:end-1, 2:end-1);
        df_t = f_corner(:, 2:end) - f_corner(:, 1:end-1);
        g.x_lr = lr.xd .* df_d + lr.xt .* df_t;
        g.y_lr = lr.yd .* df_d + lr.yt .* df_t;
        df_d = f(2:end-1, 2:end) - f(2:end-1, 1:end-1);
        df_t = f_corner(2:end, :) - f_corner(1:end-1, :);
        g.x_bt = bt.xd .* df_d + bt.xt .* df_t;
        g.y_bt = bt.yd .* df_d + bt.yt .* df_t;
        varargout{k} = g;
    end
end

function c = INVERSE_2X2(dx, dy, tx, ty)
% INVERSE_2X2  Entries of inv([dx dy; tx ty]), element-wise.
    det = dx .* ty - dy .* tx;
    c.xd =  ty ./ det;   c.xt = -dy ./ det;
    c.yd = -tx ./ det;   c.yt =  dx ./ det;
end
