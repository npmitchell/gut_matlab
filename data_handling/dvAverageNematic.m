function [mag_ap, theta_ap] = dvAverageNematic(magnitude, theta)
%[mag_ap, theta_ap] = DVAVERAGENEMATIC(magnitude, theta)
% Given a nematic field defined on a rectilinear grid, with 
% Q(i,j) being in the ith ap position and jth dv position,
% average the nematic field over the dv positions (dimension 2)
% 
% Parameters
% ----------
% magnitude : nU x nV numeric array
%   magnitude of nematic field in 2D rectilinear grid
% theta : nU x nV float array, with values as angles in radians
%   angle (mod pi) of nematic field in 2D rectilinear grid
%
% Returns
% -------
% mag_ap : nU x 1 float array
%   average magnitude of nematic field along the ap dimension
% theta_ap : nU x 1 float array
%   average angle (mod pi) of nematic field along the ap 
%   dimension
ap_x = magnitude .* cos(2*theta) ;
ap_y = magnitude .* sin(2*theta) ;
ap_xy = [mean(ap_x, 2) , mean(ap_y, 2)] ;
mag_ap = vecnorm(ap_xy, 2, 2) ;
theta_averages = atan2(ap_xy(:, 2), ap_xy(:, 1)) ;
theta_ap = 0.5 * mod(theta_averages, 2*pi) ;