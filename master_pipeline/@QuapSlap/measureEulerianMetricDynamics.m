function measureEulerianMetricDynamics(QS, options)
% measureEulerianMetricDynamics(QS, options)
%   Measure area of each circumferential slice of the QS meshes over time.
% Measures the following lengths with Eulerian (fixed) boundaries between
% lobes and folds, over time
%
%     lobe1    fold1    lobe2     fold2  lobe3   fold3  lobe4
%  <--------->|<--->|<---------->|<--->|<------>|<--->|<------>
%             ^
%             |
%          this boundary is fixed in space at a location determined by the
%          position of fold1 at its onset of foling, minus some distance
%          0.5 * options.foldW, which is taken by default to be 2% of the
%          length of the object at the onset of fold1's folding event
%
% This code can be generalized to watch for different features
% stored in QS.features instead of folds.
%
%
% NPMitchell 2020

%% Unpack QS
nU = QS.nU ;
nV = QS.nV ;
outdir = fullfile(QS.dir.mesh, 'eulerianKinematics') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
[~, ~, ~, xyzlims] = QS.getXYZLims() ;

%% Unpack options
pivMethod = 'simpleAvg' ;
writheStyle = 'Levitt' ;
overwrite = false ;
preview = false ;
foldW = 15 ;  % full width attributed to each fold, in microns

if isfield(options, 'overwrite')
    overwrite = options.overwrite  ;
end
if isfield(options, 'preview')
    preview = options.preview  ;
end
if isfield(options, 'foldW')
    foldW = options.foldW ;
end

%% Count area in each region
QS.getFeatures('folds', 'fold_onset')
folds = QS.features.folds ;
fons = QS.features.fold_onset ;

%% Create planes for separating area measurements based on folds
if size(folds, 2) == 3
    planeIds = [0*folds(fons(1), 1), ...
                folds(fons(1), 1), folds(fons(1), 1), ...
                folds(fons(2), 2), folds(fons(2), 2), ...
                folds(fons(3), 3), folds(fons(3), 3), ...
                    nU * ones(size(folds(fons(3), 1)))] ;
    labels = {'anterior lobe', 'anterior fold', ...
        'second lobe', 'middle fold', 'third lobe', ...
        'posterior fold', 'fourth lobe'} ;
else
    error('handle case of >3 or <3 folds here.')
end
planeIds = max(1, min(planeIds, nU)) ;
colors = QS.plotting.colors ;
markers = QS.plotting.markers ;

%% Convert planeIds into X positions
% Create two x=const planes for each feature (fold)
xplanes = zeros(2*length(fons), 1) ;
dmyk = 1 ;
for kk = 1:length(fons)
    fidx = fons(kk) ;
    % for this time id 
    tp = QS.xp.fileMeta.timePoints(fidx) ;
    
    % Load mesh
    fn = sprintf(QS.fullFileBase.spcutMeshSmRSC, tp) ;
    mesh = load(fn, 'spcutMeshSmRSC') ;
    mesh = mesh.spcutMeshSmRSC ;
    
    % Arrange vertices into a nU*nV*3 grid
    vtx = reshape(mesh.v, [nU, nV-1, 3]) ;
    xplanes(dmyk) = mean(vtx(planeIds(dmyk+1), :, 1)) - foldW ; 
    dmyk = dmyk + 1;
    xplanes(dmyk) = mean(vtx(planeIds(dmyk+1), :, 1)) + foldW ;
    dmyk = dmyk + 1 ;
end
disp('planeIds: ')
planeIds
disp('xplanes: ')
xplanes

%% Load face velocities from disk
QS.getVelocityAverage('vf')

%% Now find surface area of each sausage slice
fArea = zeros(length(QS.xp.fileMeta.timePoints), length(xplanes)+1) ;
for tp = QS.xp.fileMeta.timePoints
    % Obtain timepoint id for lists/cells/arrays
    tidx = QS.xp.tIdx(tp) ;
    
    disp(['t = ', num2str(tp)])
    
    % Load mesh
    fn = sprintf(QS.fullFileBase.spcutMeshSmRSC, tp) ;
    mesh = load(fn, 'spcutMeshSmRSC') ;
    mesh = mesh.spcutMeshSmRSC ;
   
    % Store vertices and faces for this timepoint
    vertices = mesh.v ;  % will append to this as more vertices added
    faces = mesh.f ;     % will append to this as more faces added
    
    % Consider
    
    % Debug test
    % mesh = read_ply_mod('/mnt/data/code/gut_matlab/mesh_handling/test_meshes/tube_simple_h1p00_R0p00_w1p00.ply') ;
    %mesh.v(:, 1) = mesh.v(:, 1) * 15  + 0.3*rand(size(mesh.v(:, 1))) ;
    %[faces, vertices] = remove_vertex_from_mesh(mesh.f, mesh.v, find(mesh.v(:, 1) > 9)) ;
    %xplanes = [5] ;
    
    % [faces, vertices] = remove_vertex_from_mesh(mesh.f, mesh.v, find(mesh.v(:, 1) > 80))
    % [faces, vertices] = remove_vertex_from_mesh(faces, vertices, find(vertices(:, 1) < 0))
    % xplanes = [30] ;
    
    tri2rm = [] ;       % will append intersected faces for removal
    addpts = [] ;       % will append new points for bisected faces
    addtri = [] ;       % will append new triangles for subsectioned faces
    flux   = zeros(length(xplanes), 1) ;
    for qq = 1:length(xplanes)
        xplane = xplanes(qq) ;
        
        % To Do: use instead:
        % [VV,FF,birth,UT,E] = (V,F,p,n,varargin)
         
        % Plane intersection with mesh
        if tidx < max(QS.xp.fileMeta.timePoints)
            [~, augv, detail] = meshPlaneIntersect(faces, vertices, xplane) ; 

            % This timepoint's velocity and relevant faces
            v0f = squeeze(QS.velocityAverage.vf(tidx, :, :)) ;
            fIntx = detail.intersected_faces ;

            % Measure the flux of tangent vector field through xplane
            [v0n, v0t] = resolveTangentNormalVector(faces(fIntx, :), mesh.v, ...
                v0f(fIntx, :)) ;
            
            % For flux, need area flux, not just velocity flux
            % Estimate perimeter of ring from radius measurement sampled
            % by weighted average of fIntx.       
            lseg = detail.line_segments ;
            dseg = vecnorm(augv(lseg(:, 1), :) ...
                           - augv(lseg(:, 2), :), 2, 2) ;
            flux(qq) = sum(v0t(:, 1) .* dseg) ;
        end
        
        % again for mesh itself
        [faces, vertices] = meshPlaneIntersect(faces, vertices, xplane) ; 
        % Add the new vertices and faces
        % addpts = cat(1, addpts, detail.added_vertices) ;
        % addtri = cat(1, addtri, detail.added_faces) ;
        % tri2rm = cat(1, tri2rm, detail.intersected_faces) ;
    end
    
    % Add the new vertices and faces
    % vertices = cat(1, vertices, addpts) ;
    % faces = cat(1, faces, addtri) ;
    % Prune triangles that were intersected. Note all vertices are kept.
    % faces(tri2rm, :) = [] ;
    
    % Check it
    if preview && mod(tp, 30) == 0
        clf
        trisurf(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), ...
             1:size(faces, 1), 'edgecolor', 'none') 
        axis equal
        view(2)
        if preview
            pause(0.0000001)
        end
        xlabel('ap position [\mum]')
        ylabel('lateral position [\mum]')
        zlabel('dv position [\mum]')
        xlim(xyzlims(1, :))
        ylim(xyzlims(2, :))
        zlim(xyzlims(3, :))
        saveas(gcf, fullfile(outdir, sprintf('cut_planes_%06d.png', tp)))
        clf
    end
    
    % Build triangulation
    tri = triangulation(faces, vertices) ;
    % cc = incenter(tri) ;
    cc = barycenter(vertices, faces) ;
    
    for qq = 1:length(xplanes)
        xplane = xplanes(qq) ;
        
        % First plane: find region behind plane
        if qq == 1
            faceInclude{qq} = find(cc(:, 1) < xplane) ;
        end
        % Find region ahead of plane, behind next plane
        if qq < length(xplanes) 
            xplaneNext = xplanes(qq + 1) ;
            faceInclude{qq + 1} = find(cc(:, 1) > xplane & ...
                                       cc(:, 1) < xplaneNext) ;            
        else        
            % Find region ahead of plane -- there is no next plane
            faceInclude{qq + 1} = find(cc(:, 1) > xplane) ;
        end        
    end
    
    % For each chunk, find the area
    for qq = 1:length(faceInclude)
        % Face areas between the planes
        fArea(tidx, qq) = sum(0.5 * doublearea(vertices, faces(faceInclude{qq}, :))) ;
    end
    
    % Store net flux into/out of each region
    for qq = 1:(length(xplanes) + 1)
        if qq == 1
            fluxes(tidx, qq) = - flux(qq) ;
        elseif qq == length(xplanes) + 1
            fluxes(tidx, qq) = flux(qq - 1) ;
        else
            % into is from prev plane, out of through next plane
            fluxes(tidx, qq) = flux(qq-1) - flux(qq) ;
        end
    end
end



%% Save it
labels = {['$x <$ ' num2str(round(xplanes(1))) '$\mu$m'], ...
    [num2str(round(xplanes(1))) '$\mu$m $<x <$ ' num2str(round(xplanes(2))) '$\mu$m'], ...
    [num2str(round(xplanes(2))) '$\mu$m $<x <$ ' num2str(round(xplanes(3))) '$\mu$m'], ...
    [num2str(round(xplanes(3))) '$\mu$m $<x <$ ' num2str(round(xplanes(4))) '$\mu$m'], ...
    [num2str(round(xplanes(4))) '$\mu$m $<x <$ ' num2str(round(xplanes(5))) '$\mu$m'], ...
    [num2str(round(xplanes(5))) '$\mu$m $<x <$ ' num2str(round(xplanes(6))) '$\mu$m'], ...
    [num2str(round(xplanes(6))) '$\mu$m $<x $']} ;

clf
t0 = QS.t0set() ;
tps = (QS.xp.fileMeta.timePoints - t0) * QS.timeInterval ;
for qq = 1:size(fArea, 2)
    if mod(qq, 2) == 0
        plot(tps, fArea(:, qq), '--')
    else
        plot(tps, fArea(:, qq), '-')
    end
    hold on;
end
legend(labels, 'location', 'eastoutside', 'Interpreter', 'Latex')
xlabel(['time [', QS.timeUnits, ']'],  'Interpreter', 'Latex')
ylabel(['area [' QS.spaceUnits '$^2$]'], 'Interpreter', 'Latex')
title('Area of mesh regions', 'Interpreter', 'Latex')
fn =  fullfile(outdir, 'areas_per_region.png') ;
disp(['saving figure: ' fn]) 
saveas(gcf, fn)
xlim([-Inf, 40])
fn =  fullfile(outdir, 'areas_per_region_zoom.png') ;
disp(['saving figure: ' fn]) 
saveas(gcf, fn) 
clf



%% Normalized plot
clf
t0Idx = QS.xp.tIdx(t0) ;
tps = (QS.xp.fileMeta.timePoints - t0) * QS.timeInterval ;
for qq = 1:size(fArea, 2)
    if mod(qq, 2) == 0
        plot(tps, fArea(:, qq) ./ fArea(t0Idx, qq), '--')
    else
        plot(tps, fArea(:, qq) ./ fArea(t0Idx, qq), '-')
    end
    hold on;
end
legend(labels, 'location', 'eastoutside', 'Interpreter', 'Latex')
xlabel(['time [', QS.timeUnits, ']'],  'Interpreter', 'Latex')
ylabel(['area $A/A_0$'], 'Interpreter', 'Latex')
title('Area of mesh regions', 'Interpreter', 'Latex')
fn = fullfile(outdir, 'areas_per_region_normalized.png') ;
disp(['saving figure: ' fn]) 
saveas(gcf, fn) 
xlim([-Inf, 40])
ylim([0.9, 1.1])
legend(labels, 'location', 'eastoutside', 'Interpreter', 'Latex')
saveas(gcf, fullfile(outdir, 'areas_per_region_normalized_zoom.png'))


%% Save fluxes as image
clf
for qq = 1:size(fArea, 2)
    if mod(qq, 2) == 0
        plot(tps, fluxes(:, qq), '--')
    else
        plot(tps, fluxes(:, qq), '-')
    end
    hold on;
end
legend(labels, 'location', 'eastOutside', 'Interpreter', 'Latex')
xlabel(['time [', QS.timeUnits, ']'],  'Interpreter', 'Latex')
ylabel(['flux [' QS.spaceUnits '/' QS.timeUnits '$^2$]'], ...
    'Interpreter', 'Latex')
title('Flux into mesh regions', 'Interpreter', 'Latex')
fn = fullfile(outdir, 'fluxes_per_region.png') ;
disp(['saving figure: ' fn]) 
saveas(gcf, fn) 
xlim([-Inf, 40])
saveas(gcf, fullfile(outdir, 'fluxes_per_region_zoom.png'))
clf


%% Save cumulative fluxes as image
clf
for qq = 1:size(fArea, 2)
    if mod(qq, 2) == 0
        plot(tps, cumsum(fluxes(:, qq)), '--')
    else
        plot(tps, cumsum(fluxes(:, qq)), '-')
    end
    hold on;
end
legend(labels, 'location', 'eastOutside', 'Interpreter', 'Latex')
xlabel(['time [', QS.timeUnits, ']'],  'Interpreter', 'Latex')
ylabel(['$\int^t \, \mathrm{d}t \, \bf{v}_{\parallel} \cdot \hat{n}$ [' QS.spaceUnits '/' QS.timeUnits '$^2$]'], ...
    'Interpreter', 'Latex')
title('Cumulative flux into mesh regions', 'Interpreter', 'Latex')
fn = fullfile(outdir, 'cumfluxes_per_region.png') ;
disp(['saving figure: ' fn]) 
saveas(gcf, fn) 
xlim([-Inf, 40])
saveas(gcf, fullfile(outdir, 'cumfluxes_per_region_zoom.png'))
clf


%% Save area and fluxes in folds only
labels = {'anterior fold area', ...
    'anterior fold flux', ...
    'middle fold area', ...
    'middle fold flux', ...
    'posterior fold area', ...
    'posterior fold flux'} ;

clf
colors = QS.plotting.colors ;
dmyk = 1 ;
for qq = 1:size(fArea, 2)
    if mod(qq, 2) == 0
        plot(tps, fArea(:, qq), '.', 'color', colors(dmyk, :))
        hold on ;
        plot(tps, cumsum(fluxes(:, qq)), '-', 'color', colors(dmyk, :))
        dmyk = dmyk + 1 ;
    end
end
legend(labels, 'location', 'eastOutside', 'Interpreter', 'Latex')
xlabel(['time [', QS.timeUnits, ']'],  'Interpreter', 'Latex')
ylabel(['$A$ [' QS.spaceUnits '$^2$], ', ...
    '$\int^t \, \mathrm{d}t \, \bf{v}_{\parallel}  \cdot \hat{n}$ [', ...
    QS.spaceUnits, '/', QS.timeUnits, '$^2$]'], ...
    'Interpreter', 'Latex')
title('areas \& fluxes: folds', 'Interpreter', 'Latex')
fn = fullfile(outdir, 'areas_fluxes_per_fold.png') ;
disp(['saving image: ' fn])
saveas(gcf, fn)
xlim([-Inf, 40])
fn = fullfile(outdir, 'areas_fluxes_per_fold_zoom.png') ;
disp(['saving image: ' fn])
saveas(gcf, fn)
clf


%% Save area and fluxes in folds only
labels = {'anterior fold area', ...
    'anterior fold flux + $A_0$', ...
    'middle fold area', ...
    'middle fold flux + $A_0$', ...
    'posterior fold area', ...
    'posterior fold flux + $A_0$'} ;

clf
colors = QS.plotting.colors ;
dmyk = 1 ;
for qq = 1:size(fArea, 2)
    if mod(qq, 2) == 0
        plot(tps, fArea(:, qq), '.', 'color', colors(dmyk, :))
        hold on;
        plot(tps, fArea(1, qq) + cumsum(fluxes(:, qq)), ...
            '-', 'color', colors(dmyk, :))
        dmyk = dmyk + 1 ;
    end
end
legend(labels, 'location', 'eastOutside', 'Interpreter', 'Latex')
xlabel(['time [', QS.timeUnits, ']'],  'Interpreter', 'Latex')
ylabel(['$A$ [' QS.spaceUnits '$^2$], ', ...
    '$\int^t \, \mathrm{d}t \, \bf{v}_{\parallel}  \cdot \hat{n}$ [', ...
    QS.spaceUnits, '/', QS.timeUnits, '$^2$]'], ...
    'Interpreter', 'Latex')
title('areas \& fluxes: folds', 'Interpreter', 'Latex')
fn = fullfile(outdir, 'areas_fluxes_per_fold_offset.png') ;
disp(['saving image: ' fn])
saveas(gcf, fn)
xlim([-Inf, 40])
fn = fullfile(outdir, 'areas_fluxes_per_fold_offset_zoom.png') ;
disp(['saving image: ' fn])
saveas(gcf, fn)
clf


%% Save area and fluxes in lobes only

labels = {'lobe 1 area', ...
    'lobe 1 flux', ...
    'lobe 2 area', ...
    'lobe 2 flux', ...
    'lobe 3 area', ...
    'lobe 3 flux', ...
    'lobe 4 area', ...
    'lobe 4 flux'} ;
clf
colors = QS.plotting.colors ;
dmyk = 1 ;
for qq = 1:size(fArea, 2)
    if mod(qq, 2) == 1
        plot(tps, fArea(:, qq), '.', 'color', colors(dmyk, :))
        hold on;
        plot(tps, cumsum(fluxes(:, qq)), '-', 'color', colors(dmyk, :))
        dmyk = dmyk + 1 ;
    end
end
legend(labels, 'location', 'eastOutside', 'Interpreter', 'Latex')
xlabel(['time [', QS.timeUnits, ']'],  'Interpreter', 'Latex')
ylabel(['$A$ [' QS.spaceUnits '$^2$], ', ...
    '$\int^t \, \mathrm{d}t \, \bf{v}_{\parallel}  \cdot \hat{n}$ [', ...
    QS.spaceUnits, '/', QS.timeUnits, '$^2$]'], ...
    'Interpreter', 'Latex')
title('areas \& fluxes: lobes', 'Interpreter', 'Latex')
fn = fullfile(outdir, 'areas_fluxes_per_lobe.png') ;
disp(['saving image: ' fn])
saveas(gcf, fn)
xlim([-Inf, 40])
fn = fullfile(outdir, 'areas_fluxes_per_lobe_zoom.png') ;
disp(['saving image: ' fn])
saveas(gcf, fn)
clf


%% Save area and fluxes + A0 in lobes only
labels = {'lobe 1 area', ...
    'lobe 1 flux + $A_0$', ...
    'lobe 2 area', ...
    'lobe 2 flux + $A_0$', ...
    'lobe 3 area', ...
    'lobe 3 flux + $A_0$', ...
    'lobe 4 area', ...
    'lobe 4 flux + $A_0$'} ;
clf
colors = QS.plotting.colors ;
dmyk = 1 ;
for qq = 1:size(fArea, 2)
    if mod(qq, 2) == 1
        plot(tps, fArea(1, qq) + cumsum(fluxes(:, qq)), ...
            '-', 'color', colors(dmyk, :))
        hold on;
        plot(tps, fArea(:, qq), '.', 'color', colors(dmyk, :))
        dmyk = dmyk + 1 ;
    end
end
legend(labels, 'location', 'eastOutside', 'Interpreter', 'Latex')
xlabel(['time [', QS.timeUnits, ']'],  'Interpreter', 'Latex')
ylabel(['$A$ [' QS.spaceUnits '$^2$], ', ...
    '$\int^t \, \mathrm{d}t \, \bf{v}_{\parallel}  \cdot \hat{n}$ [', ...
    QS.spaceUnits, '/', QS.timeUnits, '$^2$]'], ...
    'Interpreter', 'Latex')
title('areas \& fluxes: lobes', 'Interpreter', 'Latex')
fn = fullfile(outdir, 'areas_fluxes_per_lobe_offset.png') ;
disp(['saving image: ' fn])
saveas(gcf, fn)
xlim([-Inf, 40])
fn = fullfile(outdir, 'areas_fluxes_per_lobe_offset_zoom.png') ;
disp(['saving image: ' fn])
saveas(gcf, fn)
clf




disp('done with Euler Metric Dynamics')
