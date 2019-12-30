function [totalVolume, totalArea] = meshVolumeArea(pts, tri)
% MESHVOLUMEAREA(pts, tri) compute area and volume of a triangulated mesh
%   Given a surface triangulation, compute the volume enclosed using
%   divergence theorem.
%   Assumption: Triangle nodes are ordered correctly, i.e.,computed normal is outwards
% 
% Parameters
% ----------
% pts: N x 3 float
%   3D coordinates of the mesh vertices
% tri: M x 3 integer 
%   face connectivity list indexing into pts 
% 
% Returns
% -------
% totalVolume: float
%   total volume enclosed
% totalArea : float 
%   total area of surface  
%
% Authors: K. Suresh (suresh@engr.wisc.edu) and NPMitchell 2019

p = pts' ;
t = tri' ;
% Compute the vectors d13 and d12
d13= [(p(1,t(2,:))-p(1,t(3,:))); (p(2,t(2,:))-p(2,t(3,:)));  (p(3,t(2,:))-p(3,t(3,:)))];
d12= [(p(1,t(1,:))-p(1,t(2,:))); (p(2,t(1,:))-p(2,t(2,:))); (p(3,t(1,:))-p(3,t(2,:)))];
cr = cross(d13,d12,1);  %cross-product (vectorized)
area = 0.5 * sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);  % Area of each triangle
totalArea = nansum(area);
crNorm = sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);
zMean = (p(3,t(1,:))+p(3,t(2,:))+p(3,t(3,:)))/3;
nz = -cr(3,:) ./ crNorm;  % z component of normal for each triangle
volume = area .* zMean .* nz;  % contribution of each triangle
totalVolume = nansum(volume);  %divergence theorem
