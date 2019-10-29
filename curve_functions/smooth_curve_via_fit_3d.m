function [ssx, xp, yp, zp] = smooth_curve_via_fit_3d(ss, skelrs, polyorder, framelen)
%SMOOTH_CURVE_VIA_FIT Smooth a curve by fitting it to a ploynomial. 
%   Assumes reasonably uniform spacing in s for the fit (so equal weight 
%   for savgol filter).
% 
% Parameters
% ----------
% ss : the pathlength parameterization of the input curve
% skelrs : N x 3 input curve coordinates
% polyorder : int
% framelen : odd int
%
% Returns
% -------
% ssx : pathlength parameterization
% xp, yp, zp : N x 1 float arrays of the smoothed curve coords in 3d
% 
% NPMitchell 2019 

xskel = skelrs(:, 1) ;
yskel = skelrs(:, 2) ;
zskel = skelrs(:, 3) ;

% Note we ignore the variations in ds to do this fit
scx = savgol(xskel, polyorder, framelen)' ;
scy = savgol(yskel, polyorder, framelen)' ;
scz = savgol(zskel, polyorder, framelen)' ;
% sc = [scx, scy, scz] ;

% Fit the smoothed curve
xcoeffs = polyfit(ss, scx, 7) ;
ycoeffs = polyfit(ss, scy, 7) ;
zcoeffs = polyfit(ss, scz, 7) ;
ssx = linspace(min(ss), max(ss), 100) ;
xp = polyval(xcoeffs, ssx);
yp = polyval(ycoeffs, ssx);
zp = polyval(zcoeffs, ssx);

end

