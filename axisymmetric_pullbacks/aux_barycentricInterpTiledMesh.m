function [pt0, pt1, tr0, tr0_orig, tr1, tr1_orig, trisa] = aux_barycentricInterpTiledMesh(mesh0, mesh1,...
    Ysc0, xsc0, Ysc1, xsc1, ...
    use_shifty, shifty, i, x0, y0, x1, y1, x2Xpix, y2Ypix, dx2dX, dy2dY)
%AUX_BARYCENTRICINTERPTILEDMESH auxiliary function for Orbifold pipeline
% Considered abandoned at 5pm 01-02-2020
% NPMitchell 2019


% get embedded vector in R^3 for t0
% Obtain the equivalent of v2D: ie the 2D vertices: mesh0.u
mesh0x = x2Xpix(mesh0.sphi(:, 1), Ysc0, xsc0) ;
if use_shifty
    % The shift in pixels of the current frame = shifty(i)
    mesh0y = y2Ypix(mesh0.sphi(:, 2), shifty(i), Ysc0) ;
else
    mesh0y = y2Ypix(mesh0.sphi(:, 2), Ysc0) ;
end

% Create extended mesh (copied above and below)
mesh0xy = [mesh0x, mesh0y ] ;
mabove = [mesh0x, mesh0y + dy2dY(1., Ysc0)] ;
mbelow = [mesh0x, mesh0y - dy2dY(1., Ysc0)] ;
mabove2 = [mesh0x, mesh0y + dy2dY(2., Ysc0)] ;
mbelow2 = [mesh0x, mesh0y - dy2dY(2., Ysc0)] ;
m0xy = [mesh0xy; mabove; mbelow; mabove2; mbelow2] ;
% mesh faces for t0 concatenated = mf0c
mf0 = mesh0.f ;
mf0c = [mf0; mf0 + length(mesh0x); mf0 + 2 * length(mesh0x); ...
    mf0 + 3 * length(mesh0x); mf0 + 4 * length(mesh0x)] ;
tr0 = triangulation(mf0c, m0xy) ;
[t0_contain, baryc0] = pointLocation(tr0, [x0(:), Ysc0 - y0(:)]) ;    
% Interpolate the position in 3D given relative position within 2D
% triangle. Label as 'a' for t0 and 'b' for t1.
% x123(i) is the x coords of the elements of triangle t_contain(i)
vxa = mesh0.v(:, 1) ;
vya = mesh0.v(:, 2) ;
vza = mesh0.v(:, 3) ;
assert(size(vxa, 1) == size(mesh0x, 1))
% Modulo the vertex IDs: trisa are the triangle vertex IDs
tria = tr0.ConnectivityList(t0_contain, :) ;

% Make sure normals are pointing the right way
% tmp = faceNormal(tr0)
% v21 = x0(trisa(:, 2), :) - mesh0.v(trisa(:, 1), :) ;
% v31 = mesh0.v(trisa(:, 3), :) - mesh0.v(trisa(:, 1), :) ;

trisa = mod(tria, size(vxa, 1)) ;
trisa(trisa == 0) = size(vxa, 1) ;
x123a = vxa(trisa) ;
y123a = vya(trisa) ;
z123a = vza(trisa) ;
% Multiply the vertex positions by relative weights.
% Note that baryc gives weights of the three vertices of triangle
% t_contain(i) for pointLocation x0(i), y0(i)
pt0 = [sum(baryc0 .* x123a, 2), sum(baryc0 .* y123a, 2), sum(baryc0 .* z123a, 2) ] ;


% Now do next mesh
% Find xyz for matching position in t1 xy plane
% get embedded vector in R^3 for t0
% The shift in pixels of the current frame = shifty(i)
mesh1x = x2Xpix(mesh1.sphi(:, 1), Ysc1, xsc1) ;
if use_shifty
    mesh1y = y2Ypix(mesh1.sphi(:, 2), shifty(i + 1), Ysc1) ;
else
    mesh1y = y2Ypix(mesh1.sphi(:, 2), Ysc1) ;
end
% Create extended mesh (copied above and below), also shifted by shifty
mesh1xy = [mesh1x, mesh1y ] ;
mabove = [mesh1x, mesh1y + dy2dY(1., Ysc1)] ;
mbelow = [mesh1x, mesh1y - dy2dY(1., Ysc1)] ;
mabove2 = [mesh1x, mesh1y + dy2dY(2., Ysc1)] ;
mbelow2 = [mesh1x, mesh1y - dy2dY(2., Ysc1)] ;
m1xy = [mesh1xy; mabove; mbelow; mabove2; mbelow2] ;
mf1 = mesh1.f ;
mf1c = [mf1; mf1 + length(mesh1x); mf1 + 2 * length(mesh1x); ...
     mf1 + 3 * length(mesh1x); mf1 + 4 * length(mesh1x)] ;
tr1 = triangulation(mf1c, m1xy) ;
% x123(i) is the x coords of the elements of triangle t_contain(i)
[t1_contain, baryc1] = pointLocation(tr1, [x1(:), Ysc1 - y1(:)]) ;
vxb = mesh1.v(:, 1) ;
vyb = mesh1.v(:, 2) ;
vzb = mesh1.v(:, 3) ;
try
    trisb = mod(tr1.ConnectivityList(t1_contain, :), size(vxb, 1)) ;
catch
    % Plot the problem
    tmpx = x1(:);
    tmpy = y1(:) ;
    close all
    figure('visible', 'on')
    triplot(tr1)
    hold on;
    plot(x1(:), y1(:), 'g.')
    plot(tmpx(isnan(t1_contain)), tmpy(isnan(t1_contain)), 'ro')
    trisb = mod(tr1.ConnectivityList(t1_contain, :), size(vxb, 1)) ;
end
trisb(trisb == 0) = size(vxb, 1) ;
x123b = vxb(trisb) ;
y123b = vyb(trisb) ;
z123b = vzb(trisb) ;
pt1 = [sum(baryc1 .* x123b, 2), sum(baryc1 .* y123b, 2), sum(baryc1 .* z123b, 2) ] ;


tr0_orig = triangulation(mf0, mesh0xy) ;
tr1_orig = triangulation(mf1, mesh1xy) ;

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Check that pt0 is giving the same answer as in interpolation
% tileCount = [2, 2] ;
% mesh0.u = mesh0.sphi ;
% [ tm0f, tm0v2d, tm0v3d ] = tileAnnularCutMesh( mesh0, tileCount );
% tm0x = x2Xpix(tm0v2d(:, 1), Ysc0, xsc0) ;
% tm0y = y2Ypix(tm0v2d(:, 2), Ysc0) ;
% xainterp = scatteredInterpolant(tm0x, tm0y, tm0v3d(:, 1)) ;
% yainterp = scatteredInterpolant(tm0x, tm0y, tm0v3d(:, 2)) ;
% zainterp = scatteredInterpolant(tm0x, tm0y, tm0v3d(:, 3)) ;
% pt0_check = [xainterp(x0(:), y0(:)), yainterp(x0(:), y0(:)), zainterp(x0(:), y0(:))] ;
% 
% % Graphically check the quality of the interpolation
% clf; 
% scatter3(pt0(:, 1), pt0(:, 2), pt0(:, 3)); hold on;
% plot3(pt0_check(:, 1), pt0_check(:, 2), pt0_check(:, 3), 's')
% clf; 
% scatter(pt0(:, 1), pt0(:, 2)); hold on;
% plot(pt0_check(:, 1), pt0_check(:, 2), 's')
% title('explicit barycentric evaluation versus interpolation: velocity evaluation points')
% xlabel('x [pix]')
% ylabel('y [pix]')
% clf ;
% scatter3(pt0(:,1)-pt0_check(:, 1), pt0(:,2)-pt0_check(:, 2), pt0(:,3)-pt0_check(:, 3), 'markeredgealpha', 0.05)
% title('Error between barycoords and interpolation')
% xlabel('x [pix]')
% ylabel('y [pix]')
% zlabel('z [pix]')
% ylabel('y [pix]')
% pause(1)
% err = [std(pt0(:,1)-pt0_check(:, 1)), ...
%     std(pt0(:,2)-pt0_check(:, 2)), ...
%     std(pt0(:,3)-pt0_check(:, 3))] ;
% 
% err2 = [mean(abs(pt0(:,1)-pt0_check(:, 1))), ...
%     mean(abs(pt0(:,2)-pt0_check(:, 2))), ...
%     mean(abs(pt0(:,3)-pt0_check(:, 3)))] ;
% ranges = [max(pt0(:, 1)) - min(pt0(:,1)), ...
%     max(pt0(:, 2)) - min(pt0(:,2)), max(pt0(:, 3)) - min(pt0(:,3))] ;
% disp(err ./ ranges)


