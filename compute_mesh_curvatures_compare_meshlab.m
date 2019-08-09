%% MEAN & GAUSSIAN CURVATURE ==============================================
% By Noah P Mitchell & Dillon Cislo 02/2019
%==========================================================================

clear; close all; clc;
addpath('./ply_codes/')
addpath('./plotting/')
addpath('./DiscreteCurvature/')
coolwarm = coolwarm();

% Prepare path for meshes ===========================================
path = '../data/48Ygal4-UAShisRFP/20170329_objFiles/';
meshes = dir([path, 'test.ply']) ;
mesh_checkgauss = dir([path, 'test_w_gausscurv.ply']) ;
mesh_checkmean = dir([path, 'test_w_meancurv.ply']) ;
outdir = fullfile(path, 'curvature/') ;
imdir = fullfile(path, 'curvature_images/') ;
clims = [-3 3] ;
check = true ;
dx = 0.2619 * 2 ;
ii = 1;
xlims = [0, 200] ;
ylims = [-50, 50] ;
zlims = [-50, 50] ;

% Create output dirs
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
if ~exist(imdir, 'dir')
    mkdir(imdir)
end

%% Compute curvature for each mesh
disp('Loading meshes')
% for ii = 1:length(meshes)
% Load the mesh
meshfn = fullfile(meshes(ii).folder, meshes(ii).name) ;
[tri, pts, vn] = ply_read(meshfn, 'tri') ;
pts = transpose(pts) * dx;
pts(:, 3) = - pts(:, 3) + mean(pts(:, 3));
tri = transpose(tri) ;

% Check Gaussian curvature
meshfn = fullfile(mesh_checkgauss(ii).folder, mesh_checkgauss(ii).name) ;
[tri_g, pts_g, vn_g, K_g, data, comments] = ply_read(meshfn, 'tri') ;
pts_g = transpose(pts_g) * dx;
pts_g(:, 3) = - pts_g(:, 3) + mean(pts_g(:, 3));
tri_g = transpose(tri_g) ;

% Check mean curvature
meshfn = fullfile(mesh_checkmean(ii).folder, mesh_checkmean(ii).name) ;
[tri_m, pts_m, vn_m, H_m, data, comments] = ply_read(meshfn, 'tri') ;
pts_m = transpose(pts_m) * dx;
pts_m(:, 3) = - pts_m(:, 3) + mean(pts_m(:, 3));
tri_m = transpose(tri_m) ;


% View Result --------------------------------------------------------
if check
    figh = figure('visible','on');
    trisurf(tri, pts(:,1), pts(:,2), pts(:,3));
    axis equal
end
%%

% Check if any faces have normals on both sides (ie are their own front
% and back)
trisort = sort(tri, 2) ;
trisort = sortrows(trisort) ;
[trisu, ia, ib] = unique(trisort, 'rows', 'stable') ;
keep_face = false(size(trisort, 1), 1);
keep_face(ia) = true;
faces_discard = find(~keep_face) ;
% Get the rows that were duplicates
offending = trisort(~keep_face, :);
tri2rm = [];
for jj=1:size(offending, 1)
    newrm = find(all(trisort == offending(jj, :), 2)) ;
    tri2rm = [tri2rm; newrm(:)] ;
end
tri2rm = unique(tri2rm) ;
trikeep = transpose(setdiff(1:size(trisort, 1), tri2rm)) ;
face = trisort(trikeep, :) ;
vertex = pts;

%% Now compute curvatures
disp('Compute the curvatures')
[ KK, Kr ] = calculate_gaussian_curvature( face, vertex );
[ Habs, Hn ] = calculate_mean_curvature( face, vertex );

% Compute sign of mean curvature using normals
% Outward pointing unit normals on the mesh are vn
% Load these normals


% The mean curvature unit normals
nn = -Hn ./ sqrt( sum( Hn.^2, 2 ) );

% The sign of the mean curvature
sgnH = sign( dot( vn, nn, 2 ) );

% The signed mean curvature
HH = sgnH .* Habs;

% % View mean curvature normals ---------------------------------------------
% trisurf(face, vertex(:,1), vertex(:,2), vertex(:,3), ...
%     'FaceColor',[0.8 0.8 1.0]);
% axis equal
% hold on
% quiver3( vertex(:,1), vertex(:,2), vertex(:,3), ...
%     nn(:,1), nn(:,2), nn(:,3), ...
%     0.5, 'Color', 'b' );


%% VIEW RESULTS =======================================================
disp('View the results')
% View the Gaussian curvature -----------------------------------------
figh = figure('visible','off');
colormap(coolwarm)
patch('Faces', face, 'Vertices', vertex, ...
    'FaceVertexCData', KK, 'FaceColor', 'flat', ...
    'edgecolor', [0.5 0.5 0.5]);

axis equal
% caxis(clims)
clim = mean(abs(KK)) + std(abs(KK)) ;
caxis([-clim, clim])
cbar = colorbar ;
cbar.Label.String = 'Gaussian curvature, {\it K}' ;
xlabel('AP position [\mum]')
zlabel('DV position [\mum]')
title('Gaussian curvature')
xlim(xlims) ;
zlim(zlims) ;   
view(0,0) 
% Save the figure
try
    namestr = split(meshes(ii).name, '_mesh') ;
    namestr = split(namestr(1), 'pointCloud_') ;
    namestr = namestr(2) ;
    name = join(['gausscurv_', namestr, '.png'], '');
    filename = fullfile(imdir, name) ;
    saveas(figh, filename{1})
    close(figh)
catch 
    namestr = namestr{1} ;
    name = join(['gausscurv_', namestr, '.png'], '');
    filename = fullfile(imdir, name) ;
    saveas(figh, filename)
    close(figh)
end

% View the MeshLab Gaussian curvature -------------------------------------
figh = figure('visible','off');
colormap(coolwarm)
patch('Faces', tri_g, 'Vertices', pts_g, ...
    'FaceVertexCData', K_g, 'FaceColor', 'flat', ...
    'edgecolor', [0.5 0.5 0.5]);

axis equal
% caxis(clims)
clim = mean(abs(K_g)) + std(abs(K_g)) ;
caxis([-clim, clim])
cbar = colorbar ;
cbar.Label.String = 'Gaussian curvature, {\it K}' ;
xlabel('AP position [\mum]')
zlabel('DV position [\mum]')
xlim(xlims) ;
zlim(zlims) ;   
title('MeshLab Gaussian curvature')
view(0,0) 
% Save the figure
try
    namestr = split(meshes(ii).name, '_mesh') ;
    namestr = split(namestr(1), 'pointCloud_') ;
    namestr = namestr(2) ;
    name = join(['gausscurv_meshlab_', namestr, '.png'], '');
    filename = fullfile(imdir, name) ;
    saveas(figh, filename{1})
    close(figh)
catch 
    namestr = namestr{1} ;
    name = join(['gausscurv_meshlab_', namestr, '.png'], '');
    filename = fullfile(imdir, name) ;
    saveas(figh, filename)
    close(figh)
end

%%
% View the Mean curvature -----------------------------------------
figh = figure('visible','off');
colormap(coolwarm)
patch('Faces', face, 'Vertices', vertex, ...
    'FaceVertexCData', HH, 'FaceColor', 'flat', ...
    'edgecolor', [0.5 0.5 0.5]);

axis equal
% caxis(clims)
clim = mean(abs(HH)) + std(abs(HH)) ;
caxis([-clim, clim])
cbar = colorbar ;
cbar.Label.String = 'Mean Curvature, {\it H}' ;
xlabel('AP position [\mum]')
zlabel('DV position [\mum]')
xlim(xlims) ;
zlim(zlims) ;   
title('Mean curvature')
view(0,0) 
% Save the figure
try
    namestr = split(meshes(ii).name, '_mesh') ;
    namestr = split(namestr(1), 'pointCloud_') ;
    namestr = namestr(2) ;
    name = join(['meancurv_', namestr, '.png'], '');
    filename = fullfile(imdir, name) ;
    saveas(figh, filename{1})
    close(figh)
catch 
    namestr = namestr{1} ;
    name = join(['meancurv_', namestr, '.png'], '');
    filename = fullfile(imdir, name) ;
    saveas(figh, filename)
    close(figh)
end

% View the MeshLab Mean curvature -------------------------------------
figh = figure('visible','off');
colormap(coolwarm)
patch('Faces', tri_m, 'Vertices', pts_m, ...
    'FaceVertexCData', H_m, 'FaceColor', 'flat', ...
    'edgecolor', [0.5 0.5 0.5]);

axis equal
% caxis(clims)
clim = mean(abs(H_m)) + std(abs(H_m)) ;
caxis([-clim, clim])
cbar = colorbar ;
cbar.Label.String = 'Mean curvature, {\it H}' ;
xlabel('AP position [\mum]')
zlabel('DV position [\mum]')
title('MeshLab mean curvature')
xlim(xlims) ;
zlim(zlims) ;   
view(0,0) 
% Save the figure
try
    namestr = split(meshes(ii).name, '_mesh') ;
    namestr = split(namestr(1), 'pointCloud_') ;
    namestr = namestr(2) ;
    name = join(['meancurv_meshlab_', namestr, '.png'], '');
    filename = fullfile(imdir, name) ;
    saveas(figh, filename{1})
    close(figh)
catch 
    namestr = namestr{1} ;
    name = join(['meancurv_meshlab_', namestr, '.png'], '');
    filename = fullfile(imdir, name) ;
    saveas(figh, filename)
    close(figh)
end

% end