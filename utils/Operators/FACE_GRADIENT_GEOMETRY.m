function geometry = FACE_GRADIENT_GEOMETRY(mesh)
% FACE_GRADIENT_GEOMETRY  Mesh coefficients of the face-gradient reconstruction.
%
%   geometry = FACE_GRADIENT_GEOMETRY(mesh)
%
%   On each face, FACE_GRADIENTS solves [dx dy; tx ty] grad(f) = [df_d; df_t]
%   with d the centroid-to-centroid vector across the face and t the
%   corner-to-corner vector along it. This returns the entries of the
%   inverse matrices, which depend on the mesh only; GENERATE_MESH stores
%   them as mesh.face_gradient.
%
%   Inputs:
%       mesh - Mesh structure (x_Ext, y_Ext, x_corner, y_corner).
%
%   Outputs:
%       geometry - Struct with .lr and .bt, each holding the inverse entries
%                  .xd, .xt, .yd, .yt on the left-right and bottom-top faces.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    geometry.lr = INVERSE_2X2( ...
        mesh.x_Ext(2:end, 2:end-1) - mesh.x_Ext(1:end-1, 2:end-1), ...
        mesh.y_Ext(2:end, 2:end-1) - mesh.y_Ext(1:end-1, 2:end-1), ...
        mesh.x_corner(:, 2:end) - mesh.x_corner(:, 1:end-1), ...
        mesh.y_corner(:, 2:end) - mesh.y_corner(:, 1:end-1));
    geometry.bt = INVERSE_2X2( ...
        mesh.x_Ext(2:end-1, 2:end) - mesh.x_Ext(2:end-1, 1:end-1), ...
        mesh.y_Ext(2:end-1, 2:end) - mesh.y_Ext(2:end-1, 1:end-1), ...
        mesh.x_corner(2:end, :) - mesh.x_corner(1:end-1, :), ...
        mesh.y_corner(2:end, :) - mesh.y_corner(1:end-1, :));
end

function c = INVERSE_2X2(dx, dy, tx, ty)
% INVERSE_2X2  Entries of inv([dx dy; tx ty]), element-wise.
    det = dx .* ty - dy .* tx;
    c.xd =  ty ./ det;   c.xt = -dy ./ det;
    c.yd = -tx ./ det;   c.yt =  dx ./ det;
end
