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

%% Unpack options
pivMethod = 'simpleAvg' ;
writheStyle = 'Levitt' ;
overwrite = false ;
foldW = 0.04 ;  % full width of fraction of zeta attributed to each fold

if isfield(options, 'overwrite')
    overwrite = options.overwrite  ;
end
if isfield(options, 'pivMethod')
    pivMethod = options.pivMethod ;
end
if isfield(options, 'writheStyle')
    writheStyle = options.writheStyle ;
end
if isfield(options, 'foldW')
    foldW = options.foldW ;
end

%% Count area in each region
QS.getFeatures('folds', 'fold_onset')
folds = QS.features.folds ;
fons = QS.features.fold_onset ;

%% Create planes for separating area measurements based on folds
% fold half width in integers sampling zeta axis
fhw = foldW * 0.5 * nU ;
if size(folds, 2) == 3
    planeIds = [0*folds(fons(1), 1), ...
                folds(fons(1), 1) - fhw, folds(fons(1), 1) + fhw, ...
                folds(fons(2), 2) - fhw, folds(fons(2), 2) + fhw, ...
                folds(fons(3), 3) - fhw, folds(fons(3), 3) + fhw, ...
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
    xplanes(dmyk) = mean(vtx(planeIds(dmyk+1), :, 1)) ;
    dmyk = dmyk + 1;
    xplanes(dmyk) = mean(vtx(planeIds(dmyk+1), :, 1)) ;
    dmyk = dmyk + 1 ;
    
end

%% Now find surface area of each sausage slice
for tp = QS.xp.fileMeta.timePoints
    % Load mesh
    fn = sprintf(QS.fullFileBase.spcutMeshSmRSC, tp) ;
    mesh = load(fn, 'spcutMeshSmRSC') ;
    mesh = mesh.spcutMeshSmRSC ;
    
    % Consider each plane
    for qq = 1:length(xplanes)
        xlower = xplanes(qq) ;
        xupper = xplanes(qq+1) ;
         
        % Determine which faces intersect the plane
        % xright --> the faces are right of the lower bound
        % xleft  --> the faces are left of the upper bound
        % Which side of the plane is each vertex on?
        triRight = mesh.v(mesh.f, 1) > xlower ;
        xright = sum(triRight, 2) ;
        % Intersecting faces are those whose xright row is not 000 or 111
        intersect1 = find(xright < 3 & xright > 0) ;
        
        % Same thing for upper bounding plane
        triLeft = mesh.v(mesh.f, 1) < xupper ;
        xleft = sum(triLeft, 2) ;
        % Intersecting faces are those whose xleft row is not 000 or 111
        intersect2 = find(xleft < 3 & xleft > 0) ;
        
        % Faces entirely in between have xright 111 and xleft 111
        faceWithin = find(xright == 3 & xleft == 3) ;
        
        % Face areas between the planes
        fArea = 0.5 * doublearea(mesh.f(faceWithin, :), mesh.v) ;
        
        % For each face that intersects, make new triangles
        mesh_surface.faces = mesh.f ;
        mesh_surface.vertices = mesh.v ;
        [intMatrix, intSurface] = SurfaceIntersection(mesh_surface, plane1)
        for fidR = intersect1
            % vec1 = vtx(2) - vtx(1) ;
            % vec2 = vtx(3) - vtx(1) ;
            
            fArea = 0.5 * doublearea() ;
        end
        for fidR = intersect2
            
        end
    end
end