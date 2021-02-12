function generateCellSegmentation3D(QS, options)

% Unpack options
timePoints = QS.xp.fileMeta.timePoints ;

if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
end

% Setting the current timepoint clears non-timepoint segmentations
for tp = timePoints
    QS.setTime(tp)

    % Obtain the segmentation in 2d 
    seg2d = QS.getCurrentSegmentation2D(options) ;

    % obtain the current cut mesh in APDV coordinates -->
    % 2d coordinates in plane of MESH
    cutMesh = QS.getCurrentSPCutMeshSmRS() ;

    % tile annular cut mesh
    tileCount = [1,  1] ;
    [ faces, v2D, v3D ] = tileAnnularCutMesh(cutMesh, tileCount) ;
    
    % evaluation points are XY
    nVertices = length(seg2d.seg2d.Vdat) ;
    XY = zeros(nVertices, 2) ;
    for qq = 1:nVertices
        XY(qq, :) = [seg2d.seg2d.Vdat(qq).vertxcoord, seg2d.seg2d.Vdat(qq).vertycoord] ;
    end
    % Collate cell centroids
    nCells = length(seg2d.seg2d.Cdat) ;
    centroids = zeros(nCells, 2) ;
    for qq = 1:nCells
        centroids(qq, :) = seg2d.seg2d.Cdat(qq).centroid.coord ;
    end    
    
    % Check that we are interpolating over the same domain of
    % parameterization
    if strcmpi(seg2d.coordSys, 'spsme')
        Ly = size(seg2d.segIm, 1) ;
        Lx = size(seg2d.segIm, 2) ;
        doubleCovered = true ;
    else
        error(['Did not code for coordSys ' coordSys 'yet'])
    end
    uv = QS.XY2uv([Lx, Ly], XY, doubleCovered, max(cutMesh.u(:, 1)), max(cutMesh.u(:, 2))) ;
    cntrds = QS.XY2uv([Lx, Ly], centroids, doubleCovered, max(cutMesh.u(:, 1)), max(cutMesh.u(:, 2))) ;
    
    [c3d, vertexMeshFaces] = ...
        interpolate2Dpts_3Dmesh(faces, v2D, v3D, uv) ;
    
    % Get fieldfaces for centroids
    [cellCntrd, cellMeshFaces] = interpolate2Dpts_3Dmesh(faces, v2D, v3D, cntrds) ;
    fN = faceNormal(triangulation(faces, v3D)) ;
    cellNormals = fN(cellMeshFaces, :) ; 
    % Not needed: already normalized by construction
    % cellNormals = cellNormals ./ vecnorm(cellNormals, 2, 2) ;
    
    % Also get vectors point along zeta from cell centroid to lie along y
    jac = jacobian2Dto3DMesh(v2D, v3D, faces) ;
    
    %% Get aspect ratios of polygons 
    % rotate to 2d plane, with x along zeta, and y along phi
    for cid = 1:nCells
        % Obtain cell vertices in 3d
        cellVtx0 = c3d(seg2d.seg2d.Cdat(cid).nverts(:), :) ;
        cellVtx = cellVtx0 - cellCntrd(cid, :) ;
        
        zetahat = jac{cellMeshFaces(cid)} * [1, 0]' ;
        zetahat = zetahat / vecnorm(zetahat) ;
               
        rot2xy = rotate3dToAlignAxis(cellNormals(cid, :), zetahat) ;
        % Note: this one doesn't do the secondary rotation
        % rot = createRotationVector3d(cellNormals(cid, :), [0, 0, 1]) ;
        cell2d = rot2xy * cellVtx ;
        
        % Look at yz plane in which cell lives (very nearly, tangent plane)
        %   POLYGEOM( X, Y ) returns area, X centroid,
        %   Y centroid and perimeter for the planar polygon
        %   specified by vertices in vectors X and Y.
        %
        %   [ GEOM, INER, CPMO ] = POLYGEOM( X, Y ) returns
        %   area, centroid, perimeter and area moments of 
        %   inertia for the polygon.
        %   GEOM = [ area   X_cen  Y_cen  perimeter ]
        %   INER = [ Ixx    Iyy    Ixy    Iuu    Ivv    Iuv ]
        %     u,v are centroidal axes parallel to x,y axes.
        %   CPMO = [ I1     ang1   I2     ang2   J ]
        %     I1,I2 are centroidal principal moments about axes
        %         at angles ang1,ang2.
        %     ang1 and ang2 are in radians.
        %     J is centroidal polar moment.  J = I1 + I2 = Iuu + Ivv
        [ geom, iner, cpmo ] = polygeom( cell2d(:, 2), cell2d(:, 3) ) ;

        % Checking
        % assert(dot(cellNormals(cid, :), zetahat) < eps)
        subplot(2, 1, 1)
        plot3(cell2d(:, 2), cell2d(:, 3), '.-')
        subplot(2, 1, 2)
        plot3(cellVtx(:, 1), cellVtx(:, 2), cellVtx(:, 3), '.-')
        
        
    end
    
    %% Draw bonds
    Xs = zeros(4*length(c3d(:, 1)), 1) ;
    Ys = Xs ;
    Zs = Xs ;
    Us = Xs ;
    Vs = Xs ;
    Ws = Xs ;
    dmyk = 1 ;
    for qq = 1:length(seg2d.seg2d.Vdat)
        for id = seg2d.seg2d.Vdat(qq).nverts
            Xs(dmyk) = c3d(qq, 1) ;
            Ys(dmyk) = c3d(qq, 2) ; 
            Zs(dmyk) = c3d(qq, 3) ;
            Us(dmyk) = c3d(id, 1) - c3d(qq, 1) ;
            Vs(dmyk) = c3d(id, 2) - c3d(qq, 2) ; 
            Ws(dmyk) = c3d(id, 3) - c3d(qq, 3) ;
            dmyk = dmyk + 1 ;
        end
    end
    plot3(c3d(:, 1), c3d(:, 2), c3d(:, 3), '.')
    hold on;
    q = quiver3(Xs,Ys,Zs, Us, Vs, Ws, 0, 'color', [ 0.8500    0.3250    0.0980]);
    axis equal
    q.ShowArrowHead = 'off';
    [~, ~, ~, xyzlims] = QS.getXYZLims() ;
    xlim(xyzlims(1, :))
    ylim([xyzlims(2, 1), 0])
    zlim(xyzlims(3, :))
    view(0,0)
    xlabel(['ap position, [' QS.spaceUnits ']'], 'Interpreter', 'latex')
    ylabel(['lateral position, [' QS.spaceUnits ']'], 'Interpreter', 'latex')
    zlabel(['dv position, [' QS.spaceUnits ']'], 'Interpreter', 'latex')
    save(gcf, imfn) 
    
    %% Store in struct
    seg3d = QS.currentSegmentation.seg2D ;
    seg3d.Vdat.xyzrs = c3d ;
    seg3d.vertexMeshFaces = vertexMeshFaces ;
    seg3d.cellMeshFaces = cellMeshFaces ;
    seg3d.cellNormals = cellNormals ;
    
    save(outfn, 'seg3d', 'coordSys')
end


