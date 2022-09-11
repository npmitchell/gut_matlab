%% Compare methods of extracting centerline

RR = 0.6 ; % radius of coil
R1 = 0.5 ; % amount of radial undulation / variation along pathlength
debug = false ;

% make spiral tube
LL = 10 ;
zz = linspace(0, LL, 100)  ;
% width from each point along z to consider inside mesh
ww = 1 + R1 * cos(zz*2*pi/LL) ;
w2 = 1 * ones(size(ww)) % + R1 * cos(zz*0.5) ;
% points in 3D of the coil
rr = [ RR * cos(zz*pi/LL); RR * sin(zz*pi/LL); zz]' ;
r2 = [ RR * cos(zz(1)) - cos(zz*pi/LL); RR * sin(zz(1)) * ones(size(zz)); zz]' ;

% Build 3D volume segmentation
dx = 0.1 ;
xlin = -3:dx:3 ;
ylin = -3:dx:3 ;
zlin = -2:dx:LL+2 ;
[segx, segy, segz] = meshgrid(xlin, ylin, zlin) ;
seg = zeros(size(segx)) ;

% segment the volume
for ii = 1:length(zz)
    % include voxels near this point ii
    xpt = rr(ii, 1) ;
    ypt = rr(ii, 2) ;
    zpt = rr(ii, 3) ;
    seg = seg + ((segx-xpt).^2 + (segy-ypt).^2 + (segz-zpt).^2 < ww(ii)) ;
    xpt = r2(ii, 1) ;
    ypt = r2(ii, 2) ;
    zpt = r2(ii, 3) ;
    seg = seg + ((segx-xpt).^2 + (segy-ypt).^2 + (segz-zpt).^2 < w2(ii)) ;
end
seg = seg>0 ;

% Inspect seg
if debug
    for ii = 1:5:max(size(seg))
        if ii<= size(seg, 1) 
            subplot(1, 3, 1)
            imagesc(squeeze(seg(ii, :, :)))
            axis equal
        end
        if ii<= size(seg, 2) 
            subplot(1, 3, 2)
            imagesc(squeeze(seg( :, ii, :)))
            axis equal
        end
        if ii <= size(seg, 3)
            subplot(1, 3, 3)
            imagesc(squeeze(seg( :, :, ii)))
            axis equal
        end
        pause(1e-9)
    end
end

%% Make mesh 
mesh = isosurface(seg) ;
mesh.v = mesh.vertices ;
mesh.f = mesh.faces ;
disp('smoothing vertices')
% mesh.v = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], 0.04, 'implicit', mesh.v) ;

% Compute skeleton
[sy, sx, sz] = ind2sub(size(seg), find(bwskel(seg))) ;

%% Plot result
[colors, names] = define_colors ;
sky = colors(6, :) ;
disp('plotting resulting mesh')
clf
trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none',...
    'facecolor', sky, 'facealpha', 0.5)
axis equal
hold on;
plot3(sx, sy, sz, '.-')
lighting gouraud
camlight 