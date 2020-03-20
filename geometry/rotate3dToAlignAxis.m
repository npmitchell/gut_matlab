function [rot, rotx_rotz ] = rotate3dToAlignAxis(ax2alignx, varargin)
% ROTATE3DTOALIGNAXIS(ax2align, varargin)
%
% Parameters
% ----------
% ax2alignx : 1x3 float array
%   axis to rotate to the x axis
% varargin : optional, ax2alignz
%   axis to rotate to the z axis
%
% Returns 
% -------
% rot : 3x3 float array 
%   rotation matrix taking ax2alignx to the x axis, and also taking
%   ax2alignz to the z axis if supplied
% rotx_rotz : length 2 cell array of 3x3 float arrays
%   individual rotation matrices such that rot = rotz * rotx
%   ax2alignz must be supplied for this output to be returned.
% 
% NPMitchell 2020
 

% compute rotation matrix using this procedure: 
% https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
xhat = [1, 0, 0] ;
zhat = [0, 0, 1] ;
ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0] ;
RU = @(A,B) eye(3) + ssc(cross(A,B)) + ...
     ssc(cross(A,B))^2*(1-dot(A,B))/(norm(cross(A,B))^2) ;
% Normalize supplied vector
ax2alignx = ax2alignx / norm(ax2alignx) ;
% rotz aligns AP to xhat (x axis)
rotx = RU(ax2alignx, xhat) ;

if nargin > 1
    % Rotate dorsal to the z axis
    dvec = varargin{1} ;
    % Normalize the supplied vector
    dhat = dvec / norm(dvec) ;
    rotz = RU(dhat, zhat) ;
    rot = rotz * rotx  ;
    
    if nargout > 1
        rotx_rotz = {rotx, rotz} ;
    end
else
    if nargout > 1
        error('Two rotations were requested but only one rotation axis was given')
    end
    rot = rotx ;
end