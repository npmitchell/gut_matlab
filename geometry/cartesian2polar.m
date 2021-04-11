function [rr, phi] = cartesian2polar(xy)
% Convert cartesian coordinates to polar coordinates

rr = vecnorm(xy, 2, 2) ;
phi = atan2(xy(:, 2), xy(:, 1)) ;