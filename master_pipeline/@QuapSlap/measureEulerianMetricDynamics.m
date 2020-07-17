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
% This code can be easily generalized to watch for different features
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
planeIds
xplanes

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
   
    % Consider each plane
    vertices = mesh.v ;  % will append to this as more vertices added
    faces = mesh.f ;     % will append to this as more faces added
    tri2rm = [] ;        % will append intersected faces for removal
    
    % Debug test
    % mesh = read_ply_mod('/mnt/data/code/gut_matlab/mesh_handling/test_meshes/tube_simple_h1p00_R0p00_w1p00.ply') ;
    %mesh.v(:, 1) = mesh.v(:, 1) * 15  + 0.3*rand(size(mesh.v(:, 1))) ;
    %[faces, vertices] = remove_vertex_from_mesh(mesh.f, mesh.v, find(mesh.v(:, 1) > 9)) ;
    %xplanes = [5] ;
    
    % [faces, vertices] = remove_vertex_from_mesh(mesh.f, mesh.v, find(mesh.v(:, 1) > 80))
    % [faces, vertices] = remove_vertex_from_mesh(faces, vertices, find(vertices(:, 1) < 0))
    % xplanes = [30] ;
    
    addpts = [] ;
    addtri = [] ;
    for qq = 1:length(xplanes)
        xplane = xplanes(qq) ;
        
        % To Do: use instead:
        % [VV,FF,birth,UT,E] = half_space_intersect(V,F,p,n,varargin)
         
        % Determine which faces intersect the plane
        %         |            |
        % -.---.--|---.----.---|--.-----.
        %         |            |
        %         |            |
        %        intersect1   intersect2
        %
        % xright --> the faces are right of the lower bound
        % xleft  --> the faces are left of the upper bound
        % Which side of the plane is each vertex on?
        triRight = reshape(vertices(faces, 1) > xplane, size(faces)) ;
        xright = sum(triRight, 2) ;
        % Intersecting faces are those whose xright row is not 000 or 111
        intersect1 = find(xright < 3 & xright > 0) ;
                        
        % For each face that intersects, make new triangles
        % mesh_surface.faces = mesh.f ;
        % mesh_surface.vertices = mesh.v ;
        % plane1.faces = [1 2 3; 3, 4, 1] ;
        % plane1.vertices = [xlower, -ww, -ww; ...
        %                    xlower, -ww, ww; ...
        %                    xlower, ww, ww; ...
        %                    xlower, ww, -ww] ;
        % [intMatrix, intSurface] = SurfaceIntersection(mesh_surface, plane1)
        nvtx = size(vertices, 1) ;
        nfaces = size(faces, 1) ;
        tside = [] ;
        for pp = 1:length(intersect1) 
            fidR = intersect1(pp) ;
            vIds = faces(fidR, :)' ;
            vtcs = vertices(vIds, :) ;
            vtx = vtcs(:, 1) ;
            
            % Sort so that the left point(s) is(are) first
            [vtx, sortId] = sort(vtx) ;
            vIds = vIds(sortId) ;
            vtcs = vtcs(sortId, :) ;
            ds = abs(vtx - xplane) ;
            side = (vtx - xplane > 0) ;
            
            % For each segment that traverses the plane, get new point on
            % the plane
            a1 = nvtx + size(addpts, 1) + 1 ;   % index for added point #1
            a2 = nvtx + size(addpts, 1) + 2 ;   % index for added point #2
            if all(side == [0; 0; 1])
                % seg 13
                v13 = vtcs(1, :) - vtcs(3, :) ;
                addpt1 = vtcs(3, :) + ds(3) / (ds(1) + ds(3)) * v13 ;
                % seg 23
                v23 = vtcs(2, :) - vtcs(3, :) ;
                addpt2 = vtcs(3, :) + ds(3) / (ds(2) + ds(3)) * v23 ; 
                
                % Check if addpt1 and addpt2 already are added points
                if ~isempty(addpts) && any(vecnorm(addpt1 - addpts, 2, 2) < eps)
                    a1 = nvtx + find(vecnorm(addpt1 - addpts, 2, 2) < eps) ;
                    a2 = a2 - 1;
                else
                    % This is a new point. Add it to addpts
                    addpts = cat(1, addpts, addpt1) ; 
                end
                if ~isempty(addpts) && any(vecnorm(addpt2 - addpts, 2, 2) < eps) 
                    a2 = nvtx + find(vecnorm(addpt2 - addpts, 2, 2) < eps) ;
                else
                    % This is a new point. Add it to addpts
                    addpts = cat(1, addpts, addpt2) ; 
                end
                % Add triangles
                tri1 = [vIds(3), a1, a2] ;
                tri2 = [vIds(1), a1, vIds(2)] ;
                tri3 = [vIds(2), a1, a2] ;
                % Which side are these triangles on
                % tside==1 is on right, previous mesh segment=0
                % tside = [tside 1 0 0 ]; 
                addtri = [addtri; tri1; tri2; tri3] ;
                
                % Check it
                if preview
                    vtxCheck = cat(1, vertices, addpts) ;
                    newtri = [tri1; tri2; tri3] ;
                    trisurf(newtri, vtxCheck(:, 1), vtxCheck(:, 2), ...
                        vtxCheck(:, 3)) ;
                    hold on;
                    plot3(vertices(vIds(1), 1), vertices(vIds(1), 2), ...
                        vertices(vIds(1), 3), 'ro')
                    plot3(vertices(vIds(2), 1), vertices(vIds(2), 2), ...
                        vertices(vIds(2), 3), 'go')
                    plot3(vertices(vIds(3), 1), vertices(vIds(3), 2), ...
                        vertices(vIds(3), 3), 'bo')
                    plot3(addpt1(1), addpt1(2), addpt1(3), 'ko')
                    plot3(addpt2(1), addpt2(2), addpt2(3), 'ko')
                    pause(0.01)
                end
                
            elseif all(side == [0; 1; 1])
                % seg 13
                v13 = vtcs(1, :) - vtcs(3, :) ;
                addpt1 = vtcs(3, :) + ds(3) / (ds(1) + ds(3)) * v13 ;
                % seg 23
                v21 = vtcs(2, :) - vtcs(1, :) ;
                addpt2 = vtcs(1, :) + ds(1) / (ds(2) + ds(1)) * v21 ; 
                
                % Check if addpt1 and addpt2 already are added points
                if ~isempty(addpts) && any(vecnorm(addpt1 - addpts, 2, 2) < eps)
                    a1 = nvtx + find(vecnorm(addpt1 - addpts, 2, 2) < eps) ;
                    a2 = a2 - 1;
                else
                    % This is a new point. Add it to addpts
                    addpts = cat(1, addpts, addpt1) ; 
                end
                if ~isempty(addpts) && any(vecnorm(addpt2 - addpts, 2, 2) < eps) 
                    a2 = nvtx + find(vecnorm(addpt2 - addpts, 2, 2) < eps) ;
                else
                    % This is a new point. Add it to addpts
                    addpts = cat(1, addpts, addpt2) ; 
                end
                
                % Add triangles
                tri1 = [vIds(3), a1, a2] ;
                tri2 = [vIds(2), vIds(3), a2] ;
                tri3 = [vIds(1), a1, a2] ;
                % Which side are these triangles on
                % tside==1 is right, left is tside=0
                % tside = [tside 1 1 0 ]; 
                addtri = [addtri; tri1; tri2; tri3] ;
                
                % Check it
                if preview
                    vtxCheck = cat(1, vertices, addpts) ;
                    
                    % Inspect this face alone
                    newtri = [tri1; tri2; tri3] ;
                    trisurf(newtri, vtxCheck(:, 1), vtxCheck(:, 2), ...
                        vtxCheck(:, 3)) ;
                    hold on;
                    plot3(vertices(vIds(1), 1), vertices(vIds(1), 2), ...
                        vertices(vIds(1), 3), 'ro')
                    plot3(vertices(vIds(2), 1), vertices(vIds(2), 2), ...
                        vertices(vIds(2), 3), 'go')
                    plot3(vertices(vIds(3), 1), vertices(vIds(3), 2), ...
                        vertices(vIds(3), 3), 'bo')
                    plot3(addpt1(1), addpt1(2), addpt1(3), 'ko')
                    plot3(addpt2(1), addpt2(2), addpt2(3), 'ko')
                    pause(0.01)
                end
            end
            
        end
        tri2rm = cat(1, tri2rm, intersect1) ;
    end
    vertices = cat(1, vertices, addpts) ;
    faces = cat(1, faces, addtri) ;
    
    % Prune triangles that were intersected. Note all vertices are kept.
    faces(tri2rm, :) = [] ;
    
    % Check it
    if (preview || true) && mod(tp, 10) == 0
        clf
        trisurf(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), ...
             1:size(faces, 1), 'edgecolor', 'none') 
        axis equal
        view(2)
        if preview
            pause(0.1)
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
    cc = incenter(tri) ;
    
    for qq = 1:length(xplanes)
        disp(['plane ' num2str(qq)])
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
tps = (QS.xp.fileMeta.timePoints - t0) * QS.timeinterval ;
for qq = 1:size(fArea, 2)
    if mod(qq, 2) == 0
        plot(tps, fArea(:, qq), '--')
    else
        plot(tps, fArea(:, qq), '-')
    end
    hold on;
end
legend(labels, 'location', 'best', 'Interpreter', 'Latex')
xlabel(['time [', QS.timeunits, ']'],  'Interpreter', 'Latex')
ylabel(['area [' QS.spaceunits '$^2$]'], 'Interpreter', 'Latex')
title('Area of mesh regions', 'Interpreter', 'Latex')
saveas(gcf, fullfile(outdir, 'areas_per_region.png'))
xlim([-Inf, 40])
saveas(gcf, fullfile(outdir, 'areas_per_region_zoom.png'))
clf

%% Normalized plot
clf
t0Idx = QS.xp.tIdx(t0) ;
tps = (QS.xp.fileMeta.timePoints - t0) * QS.timeinterval ;
for qq = 1:size(fArea, 2)
    if mod(qq, 2) == 0
        plot(tps, fArea(:, qq) ./ fArea(t0Idx, qq), '--')
    else
        plot(tps, fArea(:, qq) ./ fArea(t0Idx, qq), '-')
    end
    hold on;
end
legend(labels, 'location', 'best', 'Interpreter', 'Latex')
xlabel(['time [', QS.timeunits, ']'],  'Interpreter', 'Latex')
ylabel(['area $A/A_0$'], 'Interpreter', 'Latex')
title('Area of mesh regions', 'Interpreter', 'Latex')
saveas(gcf, fullfile(outdir, 'areas_per_region_normalized.png'))
xlim([-Inf, 40])
ylim([0.9, 1.1])
legend(labels, 'location', 'best', 'Interpreter', 'Latex')
saveas(gcf, fullfile(outdir, 'areas_per_region_normalized_zoom.png'))


disp('done with Euler Metric Dynamics')
