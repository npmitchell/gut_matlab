function [tangent, normal, binormal] = frenetSeretFrame(ss, xp, yp, zp)
%frenetSeretFrame Compute Frenet-Seret frame for a 1d curve in 3d
%   Detailed explanation goes here

dsx = gradient(ss) ;

% First calc rate of change of curve along curve
gradc_raw = [gradient(xp'), gradient(yp'), gradient(zp')] ; 
gradc = bsxfun(@rdivide, gradc_raw, dsx(:)) ;
gradc_ds = vecnorm(gradc, 2, 2) ;
% Compute the tangent to the curve
tangent = bsxfun(@rdivide, gradc, gradc_ds(:)) ;

% Compute normal 
normal_raw = [gradient(tangent(:, 1)), ...
    gradient(tangent(:, 2)), ...
    gradient(tangent(:, 3))] ; 
normalc = bsxfun(@rdivide, normal_raw, dsx(:)) ;
normalc_ds = vecnorm(normalc, 2, 2) ;
normal = bsxfun(@rdivide, normalc, normalc_ds(:)) ;

% Compute binormal from cross product
binormalc = cross(tangent, normal) ;
binormalc_ds = vecnorm(binormalc, 2, 2) ;
binormal = bsxfun(@rdivide, binormalc, binormalc_ds(:)) ;
end

