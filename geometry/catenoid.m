function [zz, rr] = catenoid(Rphi, Rz, xy)
% Parameterization of catenoid / barrel

xx = xy(:, 1) ;
yy = xy(:, 2) ;
rr = Rphi + Rz - Rz * cosh(yy/Rz) ;
zz = rr .* cos(asin( xx ./ (Rphi + Rz - Rz * cosh(yy/Rz)))) ;