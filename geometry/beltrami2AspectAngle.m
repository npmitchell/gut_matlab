function [ar, theta] = beltrami2AspectAngle(mus)
%[ar, theta] = beltrami2AspectAngle(mus)
% Convert beltrami coefficient mu = a+ib into aspect ratio and elongation 
% angle 
%
ar = (1 + abs(mus)) ./ (1-abs(mus)) ;
theta = 0.5 * angle(mus) ;