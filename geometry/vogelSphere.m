function points = vogelSphere(npts)
% generate unit sphere with npts points approximately equidistant on a sphere
%

golden_angle = pi * (3 - sqrt(5)) ;
theta = golden_angle * (1:npts) ;
zcap = linspace(1 - 1.0 / npts, 1.0 / npts - 1, npts) ;
radius = sqrt(1 - zcap .* zcap) ;
points = zeros([npts, 3]) ;
points(:,1) = radius .* cos(theta) ;
points(:,2) = radius .* sin(theta) ;
points(:,3) = zcap ;