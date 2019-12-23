function [ssdiff] = matchCurves(vars, curv1, curv2)
% Given two curves, interpolate the first in reference to the second and
% find the sum squared differences betweeen the two curves as a function of
% a time offset. This can be used to minimize difference between two time
% offset curves. 
%
% parameters
% ----------
% vars : a vector containing the variables to be minimized
%
% vars1 : float
%   time to shift curv1 relative to curv2's t axis to match the two curves
% vars2 : float
%   time dilation factor
% vars3: float
%   vertical offset of function
% curv1 : N x 1 float array
%   first column is timestamp, second is a value f(t)
%
% Returns 
% -------
% ssdiff : float
%   the sum of the squared differences between curv1(t2) and curv2(t2),
%   where t2 are the timestamps given in curv2
vars1 = vars(1); % time offset
vars2 = vars(2); % time stretch 
vars3 = vars(3); % vertical offset
vars4 = vars(4); % vertical stretch
% Unpack 
t1 = curv1(:, 1) ;
t2 = curv2(:, 1) ;
p1 = vars4*curv1(:, 2) + vars3 ;
p2 = curv2(:, 2) ;
toff = vars2*t1 + vars1 ;

% Interpolate curv1 at times t2
curv1intrp = interp1(toff, p1, t2, 'pchip');
ssdiff = sum((curv1intrp(t2) - p2).^2) ;

end

