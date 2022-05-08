function angles = bondAngles(BL, uv) 
% bondAngles(BL, uv)
% angles of bonds in a connected network
% 
% Parameters
% ----------
% BL : #bonds x 2 int
%   bond list, indices into uv of connected vertices
% uv : #vertices x 2 numeric
%   vertex positions in 2d
%
% Returns
% -------
% angles : #bonds x 1 float
%   bond angles in plane of vertices uv
%

% take arctangent of vertex2 - vertex1
angles = atan2(uv(BL(:, 2), 2)-uv(BL(:, 1), 2), ...
    uv(BL(:, 2), 1)-uv(BL(:, 1), 1))  ;