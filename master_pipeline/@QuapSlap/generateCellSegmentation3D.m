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
    QS.getCellSegmentation2D(options)

    % obtain the current cut mesh in APDV coordinates
    QS.getCurrentSPCutMeshSmRS() ;

    % 2d coordinates in plane of MESH
    cutMesh = QS.currentMesh.spcutMeshSmRS ;

    % tile annular cut mesh
    tileCount = [1,  1] ;
    TM = tileAnnularCutMesh(cutMesh, tileCount) ;
    v2d = TM.u  ;
    v3d = TM.v  ;
    faces = TM.f ;

    % evaluation points are uv
    uv = cellsSegmentation.cells2D.Vdat.xy ;

    [c3d, fieldfaces] = ...
        interpolate2Dpts_3Dmesh(faces, v2d, v3d, uv) ;

    % Store in struct
    QS.currentSegmentation.seg3D = QS.currentSegmentation.seg2D ;
    QS.currentSegmentation.seg3D.
end