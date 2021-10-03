function visualizeTracking3D(QS, options)
% visualizeTracking3D(QS, options)
%
%
% Parameters
% ----------
% options : struct with fields
%   timePoints : timepoints to visualize
%   method : 'segmentation' or 'nuclei'
%   subdir : 'labeled_groundTruth'
%
%

timePoints = QS.xp.fileMeta.timePoints ;
subdir = 'labeled_groundTruth';
method = 'segmentation' ;
coordSys = 'spsme' ;
if isfield(options, 'timePoints')
    timePoints = optoins.timePoints ;
end
if isfield(options, 'method')
    method = optoins.method ;
end
if isfield(options, 'subdir')
    subdir = optoins.subdir ;
end
if isfield(options, 'coordSys')
    coordSys = options.coordSys ;
end


for tidx = 1:length(timePoints)
    tp = timePoints(tidx) ;
    QS.setTime(tp) ;
    if strcmpi( method, 'segmentation')
        % load segmentation tracks
        segIm = load(fullfile(QS.dir.tracking, options.subdir, ...
                sprintf('tracks_label_%06d.mat', tp))) ;
            
        % load mesh and image
        if strcmpi(coordSys, 'spsme')
            % Get the image size and mesh
            im = imread(sprintf(QS.fullFileBase.im_sp_sme, timePoints(1))) ;
            cutMesh = QS.getCurrentSPCutMeshSm() ;
            glueMesh = QS.getCurrentSPCutMeshSmRSC() ;
            cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1)) ;
            glueMesh.u(:, 1) = glueMesh.u(:, 1) / max(glueMesh.u(:, 1)) ;
            shiftY = size(im, 1) * 0.5 ;
            doubleCovered = true ;
        else
            error('did not recognize coordSys')
        end
        
        [~, ~, ~, xyzlims] = QS.getXYZLims() ;
        nCells = max(segIm(:)) ;
        
        % polygons and centroids in 2d 
        props = regionprops(segIm, 'centroid') ;
        seg2d.cdat = struct() ;
        seg2d.cdat.centroid = zeros(max(segIm(:)), 2) ;
        for cid = 1:max(segIm(:))
            seg2d.cdat.centroid(cid, :) = props(cid).Centroid ;
        end
        polygons = polygonsFromSegmentation(segIm) ;        
        seg2d.cdat.polygons = polygons ;
        centroids = seg2d.cdat.centroid ;
        
        % Tile annular cutmesh
        [ TF, TV2D, TV3D, TVN3D ] = tileAnnularCutMesh(cutMesh, [1,1]) ;
        
        % Create 3d mapping of each polygon of the tracked segmentation
        % centroids in 3d
        cuv = QS.XY2uv( im, centroids, doubleCovered, 1, 1) ;
        [c3D, cMeshFaces] = interpolate2Dpts_3Dmesh(TF, TV2D, TV3D, cuv) ;
        for qq = 1:max(segIm)
            % polygons in 3d
            puv = QS.XY2uv( im, polygons{qq}, doubleCovered, 1, 1) ;
            [p3D, cMeshFaces] = interpolate2Dpts_3Dmesh(TF, TV2D, TV3D, puv) ;
        end
        
        allVertices
        opts = struct() ;
        opts.centroids = c3D ;
        polygonsToPatchTriangles3D(p3d, opts)
        patch('Faces',ff,'Vertices',vv,...
            'FaceVertexCData',cos(2*ang1V(:)),'FaceColor','flat', ...
            'Edgecolor', 'none');
        axis equal
        xlim(xyzlims(1, :))
        ylim(xyzlims(2, :))
        zlim(xyzlims(3, :))
        
        % Plot mesh for occlusion
        
    else
        error('handle nuclear tracking case here')
    end
    
end



