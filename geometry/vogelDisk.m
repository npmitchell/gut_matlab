function [xy, rphi] = vogelDisk(npts)

golden_angle = pi * (3 - sqrt(5)) ;
rho = sqrt(1:npts) / sqrt(npts) ;
phi = (1:npts) * golden_angle ;

% convert to xy cartesian
xy = [rho .* cos(phi); rho .* sin(phi)]' ;

if nargout > 1
    rphi = [rho; phi]' ;
end

