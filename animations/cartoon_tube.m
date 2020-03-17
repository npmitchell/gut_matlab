
addpath('/mnt/data/code/gut_matlab/addpath_recurse/')
addpath('/mnt/data/code/gut_matlab/plotting/')
addpath('/mnt/data/code/gut_matlab/plotting/export_fig/')
addpath('/mnt/data/code/gut_matlab/mesh_handling/')
addpath_recurse('/mnt/data/code/gptoolbox/')

x = linspace(-1, 1, 100) ;
theta = (0:0.1:(2 * pi + 0.1))' ;
xx = (ones(size(theta)) * x)' ;
% Create a dummy y array
y = linspace(0, 1, length(theta)) ;
yy = ones(size(x))' * y ;
xy = [xx(:), yy(:)] ;
faces = defineFacesRectilinearGrid(xy, length(x), length(theta)) ;
outdir = './cartoon_tube_persp/';
mkdir(outdir) ;
outdir3 = './cartoon_tube_persp_local/';
mkdir(outdir3) ;
% outdir2 = './cartoon_tube/';
% mkdir(outdir2) ;
% outdir4 = './cartoon_tube_local/';
% mkdir(outdir4) ;

ratio = 0.5 ;
dmyk = 0 ;
tmax = 0.8 ;
sigma = 0.1 ;
mu = 0. ;

% Get final colors
y = 1 - tmax * exp(-(x-mu).^2 / sigma.^2) ;
cfinal = (ones(size(theta)) * y - 1)' ;

for t = 0:0.01:tmax
    disp(['t=' num2str(t)])
    
    
    y = 1 - t * exp(-(x-mu).^2 / sigma.^2) ;
    radius = ones(size(theta)) * y ;
    
    yy = (cos(theta) * ones(size(x))) .* radius ;
    zz = (sin(theta) * ones(size(x))) .* radius ;
    cc = (radius - 1)';
    yy = yy' ;
    zz = zz' ;
    
    % Plot the tube outline
    % plot(x, y1, 'k-')
    % hold on;
    % plot(x, y2, 'k-')
    set(gcf, 'visible', 'off')
    h = trisurf(faces, xx, yy, zz, cc, 'EdgeColor', 'none') ; 
    caxis([-tmax, tmax])   
    colormap(bwr)
    axis equal
    axis off
    [AO,C,l] = apply_ambient_occlusion(h, 'Factor', 0.4, 'Samples', 1000, 'AddLights', true) ;
    camlight(90,  20) 
    
    % Now plot the flow along the tube
    % xg = imresize(xx, [20, 20]) ;
    % yg = imresize(yy, [20, 20]) ;
    % zg = imresize(zz, [20, 20]) ;
    % keep = abs(xg) > 0.2 ;
    % xg = xg(keep) ;
    % yg = yg(keep) ;
    % zg = zg(keep) ;
    % vv = -0.1 * xg ./ abs(xg) ; 
    % hold on; 
    % q  = quiver3(xg, yg, zg, vv, 0 * vv, 0 * vv, 0 ) ;
    % q.Color = [0, 162, 255];
    
    set(gca, 'Color', 'k')
    set(gcf, 'color', 'k')
    axis equal
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
    set(gca, 'ztick', [])
    view(-20, 30) ;
    export_fig([outdir sprintf('%06d', dmyk) '.png'])
    % view(2)
    % export_fig([outdir2 sprintf('%06d', dmyk) '.png'])
    clf
    
    % Same with full color
    set(gcf, 'visible', 'off')
    h = trisurf(faces, xx, yy, zz, cfinal, 'EdgeColor', 'none') ; 
    caxis([-tmax, tmax])
    
    set(gca, 'Color', 'k')
    set(gcf, 'color', 'k')
    axis equal
    axis off
    % axis off
    
    [AO,C,l] = apply_ambient_occlusion(h, 'Factor', 0.4, 'Samples', 1000, 'AddLights', true) ;
    camlight(90,  20) 
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
    set(gca, 'ztick', [])
    colormap(bwr)
    view(-20, 30) ;
    export_fig([outdir3 sprintf('%06d', dmyk) '.png'])
    % view(2)
    % export_fig([outdir4 sprintf('%06d', dmyk) '.png'])
    clf
    
    dmyk = dmyk + 1 ;
end