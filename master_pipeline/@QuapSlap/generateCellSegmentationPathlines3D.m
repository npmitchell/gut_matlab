function generateCellSegmentationPathlines3D(QS, options)
% generateCellSegmentation3D(QS, options)
% 
%  Load segmentation results from a timepoint t0Pathlines, advect the cell
%  polygons along pathlines from PIV, measure the segmentation properties
%  of the advected "tissue" pattern frozen in the Lagrangian frame as the
%  material deforms.
%
% NPMitchell 2021


% Unpack options
timePoints = QS.xp.fileMeta.timePoints ;
maxCellSize = Inf ;  % maximum size allowed to include a cell
overwrite = false ;
overwriteImages = true ;
useCorrected = true ;
debug = false ;
preview = false ;

[~, ~, ~, xyzlims] = QS.getXYZLims() ;
t0 = QS.t0set() ;

if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
end
if isfield(options, 'debug')
    debug = options.debug ;
end
if isfield(options, 'xyzlims')
    xyzlims = options.xyzlims ;
end
if isfield(options, 'useCorrectedSegmentation')
    useCorrected = options.useCorrectedSegmentation ;
elseif isfield(options, 'useCorrected')
    useCorrected = options.useCorrected ;
end
if isfield(options, 'segmentationPathlines')
    tmp = options.segmentationPathlines ;
    segVertexPathlines2D = tmp.segVertexPathlines2D ;
    segVertexPathlines3D = tmp.segVertexPathlines3D ;
    cellIDs = tmp.cellIDs ;
    segPathlinesPassedAsOption = true ;
else
    segPathlinesPassedAsOption = false ;
end

% Directory preparation
imdir = fullfile(QS.dir.segmentation, 'pathlines', 'images') ;
if ~exist(imdir, 'dir')
    mkdir(imdir)
end

% Prep for dicing data by lobes
features = QS.getFeatures() ;
folds = features.folds ;
nLobes = size(folds, 2) + 1 ;

%% Create synthetic pathlines for cells from t=0
% synethetic cell vertex pathlines in pullback space
if ~exist(fullfile(QS.dir.segmentation, 'pathlines'), 'dir')
    mkdir(fullfile(QS.dir.segmentation, 'pathlines'))
end
cellVertexPathlineFn = fullfile(QS.dir.segmentation, 'pathlines', ...
    sprintf('cellVertexPathlines_%06dt0.mat', t0)) ;
if ~exist(cellVertexPathlineFn, 'file') || overwrite 
    disp('Creating vertex pathlines from scratch...')
    QS.setTime(t0) ;
    if useCorrected
        seg2d = getCurrentSegmentation2DCorrected(QS) ;
        cellV0 = [] ;
        cellIDs = [] ;
        for cellID = 1:length(seg2d.seg2d.cdat.polygons)
            cVtx = seg2d.seg2d.cdat.polygons{cellID} ;
            if ~isempty(cVtx)
                cellV0 = [cellV0; [cVtx(:, 2), cVtx(:, 1)] ] ;            
                cellIDs = [cellIDs; cellID * ones(size(cVtx, 1), 1)] ;
            end
        end
    else
        seg2d = getCurrentSegmentation2D(QS) ;
        cellV0 = seg2d.seg2d.vdat.v ;
        cellIDs = 1:length(seg2d.seg2d.cdat.polygons) ;
    end
    opts = struct('preview', true) ;
    [segVertexPathlines2D, segVertexPathlines3D] = ...
        QS.samplePullbackPathlines(cellV0, opts) ;
    save(cellVertexPathlineFn, 'segVertexPathlines2D', ...
        'segVertexPathlines3D', 'cellIDs')
else
    disp('Loading vertex pathlines from segmentation advection...')
    if ~segPathlinesPassedAsOption
        load(cellVertexPathlineFn, 'segVertexPathlines2D', ...
            'segVertexPathlines3D', 'cellIDs')
    end
end

% Load cell segmentation 2d at t0
QS.setTime(t0) ;
if useCorrected
    seg02d = QS.getCurrentSegmentation2DCorrected() ;
else
    seg02d = QS.getCurrentSegmentation2D() ;
end

% Setting the current timepoint clears non-timepoint segmentations
close all
mc2t = [] ;
ms2t = [] ;
c2t_low25 = [] ;
c2t_high75 = [] ;
c2t_std = [] ;
c2t_ste = [] ;
s2t_low25 = [] ;
s2t_high75 = [] ;
s2t_std = [] ;
s2t_ste = [] ;
mean_mratio = [] ;
mean_moiratio = [] ;
median_moiratio = [] ;
mratio_low25 = [] ;
mratio_high75 = [] ;
mratio_std = [] ;
mratio_ste = [] ;
dmy = 1 ;
cos2thetaM = zeros(99, 1) ;
sin2thetaM = cos2thetaM ;
aspectM = cos2thetaM ;
nAPBins = 20 ;
mean_qc2ts = zeros(length(timePoints), nAPBins) ;
mean_qs2ts = zeros(length(timePoints), nAPBins) ;
std_qc2ts = zeros(length(timePoints), nAPBins) ;
std_qs2ts = zeros(length(timePoints), nAPBins) ;
mean_c2ts = zeros(length(timePoints), nAPBins) ;
mean_s2ts = zeros(length(timePoints), nAPBins) ;
meanQLobeAspects = zeros(nLobes, length(timePoints)) ;
meanQLobeAspectStds = zeros(nLobes, length(timePoints)) ;
meanQLobeThetas = zeros(nLobes, length(timePoints)) ;
for tp = timePoints
    sprintf(['t = ' num2str(tp)])
    QS.setTime(tp)
    tidx = QS.xp.tIdx(tp) ;
    
    outfn = fullfile(QS.dir.segmentation, 'pathlines', ...
        sprintf([QS.fileBase.segmentation3d '.mat'], QS.currentTime)) ;
    if ~exist(outfn, 'file') || overwrite 

        % Obtain the segmentation in 2d --> use reference segmentation
        seg2d = seg02d ;
        % replace vertices with advected ones from t0
        seg2d.seg2d.vdat.v = squeeze(segVertexPathlines2D(tidx, :, :)) ;
        
        % replace centroids with correct ones (ie advected in XY plane)
        if ~useCorrected
            % Polygons are indices into vertices
            nCells = length(seg2d.seg2d.Cdat) ;
            for pp = 1:nCells 
                poly = seg2d.seg2d.cdat.polygons{pp} ;
                if ~isempty(poly)
                    try
                        geom = polygeom( seg2d.seg2d.vdat.v(poly, 1), ...
                            seg2d.seg2d.vdat.v(poly, 2) ) ;
                        seg2d.seg2d.cdat.centroid(pp, :) = geom(2:3) ; 
                    catch
                        error('Could not get polygon information in 3D')
                    end
                end
            end
            
            % Flush superfluous info that was loaded from memory
            seg2d.seg2d.Cdat = rmfield(seg2d.seg2d.Cdat, 'centroid') ;
            seg2d.seg2d.Vdat = rmfield(seg2d.seg2d.Vdat, 'vertxcoord') ;
            seg2d.seg2d.Vdat = rmfield(seg2d.seg2d.Vdat, 'vertycoord') ;
        else
            % Polygons are computed from vertices indexed by cid
            nCells = length(seg2d.seg2d.cdat.polygons) ;
            
            if preview
                tmpCOMs = seg2d.seg2d.cdat.centroid ;
            end
            
            for cid = 1:nCells
                poly = seg2d.seg2d.vdat.v(cellIDs == cid, :) ;
                if ~isempty(poly)
                    geom = polygeom( poly(:, 1), poly(:, 2) ) ;
                    seg2d.seg2d.cdat.centroid(cid, :) = geom(2:3) ; 
                end
            end
            
            % THe COMs should appear suitably advected
            if preview
                plot(tmpCOMs(:, 1), tmpCOMs(:, 2), '.')
                hold on; plot(seg2d.seg2d.cdat.centroid(:, 1), ...
                    seg2d.seg2d.cdat.centroid(:, 2), 'o')
            end
        end
        
        % obtain the current cut mesh in APDV coordinates -->
        % 2d coordinates in plane of MESH
        cutMesh = QS.getCurrentSPCutMeshSmRS() ;

        % Normalize the u coordinate
        cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1)) ;

        % tile annular cut mesh
        tileCount = [1,  1] ;
        [ faces, v2D, v3D ] = tileAnnularCutMesh(cutMesh, tileCount) ;

        % Collate cell vertices as XY
        XY = seg2d.seg2d.vdat.v ;

        % Collate cell centroids
        centroids = seg2d.seg2d.cdat.centroid ;
        
        % Check that we are interpolating over the same domain of
        % parameterization
        if strcmpi(seg2d.coordSys, 'spsme')
            Ly = size(seg2d.segIm, 1) ;
            Lx = size(seg2d.segIm, 2) ;
            doubleCovered = true ;
        else
            error(['Did not code for coordSys ' coordSys 'yet'])
        end
        % uv coordinates of the cell vertices
        umax = max(cutMesh.u(:, 1)) ;
        vmax = max(cutMesh.u(:, 2)) ;
        uv = QS.XY2uv([Lx, Ly], XY, doubleCovered, umax, vmax) ;

        % Centroids in uv coordinates
        cntrds = QS.XY2uv([Lx, Ly], centroids, doubleCovered, umax, vmax) ;

        % cell vertices in 3d
        [c3d, vertexMeshFaces] = ...
            interpolate2Dpts_3Dmesh(faces, v2D, v3D, uv) ;

        % Get fieldfaces for centroids
        [cellCntrd, cellMeshFaces] = interpolate2Dpts_3Dmesh(faces, v2D, v3D, cntrds) ;
        fN = faceNormal(triangulation(faces, v3D)) ;
        cellNormals = fN(cellMeshFaces, :) ; 
        % Not needed: already normalized by construction
        % cellNormals = cellNormals ./ vecnorm(cellNormals, 2, 2) ;

        % Also get vectors point along zeta from cell centroid to lie along y
        jac2d3d = jacobian2Dto3DMesh(v2D, v3D, faces) ;
        % jac3d2d = jacobian3Dto2DMesh(v2D, v3D, faces) ;

        %% Get aspect ratios of polygons 
        % rotate to 2d plane, with x along zeta, and y along phi
        areas = nan(nCells, 1) ;
        perim = nan(nCells, 1) ;
        moment1 = nan(nCells, 1) ;
        ang1 = nan(nCells, 1) ;
        moment2 = nan(nCells, 1) ;
        ang2 = nan(nCells, 1) ;
        moinertia = nan(nCells, 2, 2) ;
        cellQ2d = {} ;

        for cid = 1:nCells
            % Obtain cell vertices in 3d
            if useCorrected
                cell2d0 = seg2d.seg2d.vdat.v(cellIDs == cid, :) ;
                cellVtx0 = c3d(cellIDs == cid, :) ;
            else
                cell2d0 = seg2d.seg2d.vdat.v(seg2d.seg2d.cdat.polygons{cid}, :) ;
                cellVtx0 = c3d(seg2d.seg2d.cdat.polygons{cid}, :) ;
            end
            cellVtx = cellVtx0 - cellCntrd(cid, :) ;

            if ~isempty(cellVtx)

                % Note: this approach is tricky since we have to map to
                % 3d and then back to 2d.
                % To figure out which direction to take to z, map vec to 3d

                % dzeta3d points towards the mapped z axis but in original 3d
                % embedding space
                dzeta3d = (jac2d3d{cellMeshFaces(cid)} * [0, 1]')' ;
                dzeta3d = dzeta3d / vecnorm(dzeta3d, 2, 2) ;

                % rotate cell to nearly yz plane
                rot2xy = rotate3dToAlignAxis(cellNormals(cid, :), dzeta3d) ;
                % Note: this one doesn't do the secondary rotation
                % rot = createRotationVector3d(cellNormals(cid, :), [0, 0, 1]) ;
                cell_quasi2d = (rot2xy * cellVtx')' ;
                
                % % Alternative method based on jacobian --> actually this 
                % % is a bit limited since it assumes small cells 
                % % (same size or smaller than the face) to get the 
                % % stretching in each dimension correct
                % jac = jac3d2d{cellMeshFaces(cid)} ;
                % dilation = jacobian2dilation(jac) ;
                % cell2d = (jac * cellVtx0')' ;
                % mapped centroid should be the same as pre-mapped centroid
                % up to rescaling in each dimension separately
                % cntrd2d = jac * cellCntrd(cid, :)' ;
                % cntrds(cid, :)

                % Check 2d cell polygon
                if debug
                    disp('debugging cell polygon in 3d and quasi2d...')
                    clf
                    % dchi3d points towards the mapped x axis in embedding space
                    dchi3d = (jac2d3d{cellMeshFaces(cid)} * [1, 0]')' ;
                    dchi3d = dchi3d / vecnorm(dchi3d, 2, 2) ;
                    
                    % Plot the cell in 3d and 2d
                    subplot(2, 2, 1)
                    plot(cell2d0(:, 1), cell2d0(:, 2), '.-');
                    axis equal
                    title('cell in 2d pullback XY', 'interpreter', 'latex')
                    subplot(2, 2, 2)

                    % Vector to transform = dzeta since this emanates from
                    % centroid just as the normal does.
                    zeta2d = (rot2xy * dzeta3d')' ;
                    try
                        assert(all(abs(zeta2d - [0, 0, 1]) < 1e-7))
                    catch
                        error('Rotation was not successful!')
                    end

                    % plot it
                    plot3(cell_quasi2d(:, 1), cell_quasi2d(:, 2), ...
                        cell_quasi2d(:, 3), '.-'); 
                    axis equal
                    title('cell projected onto tangent plane and oriented', ...
                        'interpreter', 'latex')
                    hold on; 

                    subplot(2, 2, 3)
                    plot3(cellVtx0(:, 1), cellVtx0(:, 2), cellVtx0(:, 3), ...
                        '.-')
                    hold on;
                    zplus = cellCntrd(cid, :) + dzeta3d * mean(var(cellVtx0)); 
                    xplus = cellCntrd(cid, :) + dchi3d * mean(var(cellVtx0));
                    plot3dpts([cellCntrd(cid, :); zplus], 'r-')
                    plot3dpts([cellCntrd(cid, :); xplus], 'g-')
                    title('cell in 3d with computed $(\zeta, \chi)$ coordinates', ...
                        'interpreter', 'latex')
                    axis equal
                end

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

                % Discard 3d info and compute
                if size(cellVtx0, 1) > 2
                    try
                        [ geom, iner, cpmo ] = polygeom( cell_quasi2d(:, 2), ...
                            cell_quasi2d(:, 3) ) ;
                    catch
                        error('here')
                    end
                    areas(cid) = geom(1) ;
                    perim(cid) = geom(4) ;
                    cellQ2d{cid} = cell_quasi2d ;
                    moment1(cid) = cpmo(1) ;
                    ang1(cid) = cpmo(2) ;
                    moment2(cid) = cpmo(3) ;
                    ang2(cid) = cpmo(4) ;
                    mratio(cid) = moment2(cid) / moment1(cid) ;
                    moinertia(cid, :, :) = [iner(4) -iner(6); -iner(6) iner(5)] ;
                    
                end
            else
                disp(['bad cell: ' num2str(cid)])
            end
        end
        
        %% Save results stored in struct
        seg3d = struct('vdat', struct(), 'cdat', struct(), ...
            'qualities', struct(), 'map', struct()) ;
        
        % Mesh information for mapping
        seg3d.map.f = faces ;
        seg3d.map.v2d = v2D ;
        seg3d.map.v3d = v3D ;
        
        % vertex data in 3d
        seg3d.vdat.xyzrs = c3d ;
        seg3d.vdat.uv = uv ;
        seg3d.vdat.meshFaces = vertexMeshFaces ;
        if useCorrected
            seg3d.cellIDs = cellIDs ;
        else
            seg3d.vdat.NL = seg2d.seg2d.vdat.NL ;
            seg3d.vdat.BL = seg2d.seg2d.vdat.BL ;
            seg3d.vdat.fourfold = seg2d.seg2d.vdat.fourfold ;
        end
        
        % cell data in 3d
        seg3d.cdat.centroids_uv = cntrds ;
        seg3d.cdat.centroids_3d = cellCntrd ;
        seg3d.cdat.meshFaces = cellMeshFaces ;
        seg3d.cdat.normals = cellNormals ;
        seg3d.cdat.polygons = seg2d.seg2d.cdat.polygons ;
        
        % cell qualities
        seg3d.qualities = struct() ;
        seg3d.qualities.areas = areas ;
        seg3d.qualities.perim = perim ;
        seg3d.qualities.moment1 = moment1 ;
        seg3d.qualities.moment2 = moment2 ;
        seg3d.qualities.ang1 = ang1 ;
        seg3d.qualities.ang2 = ang2 ;
        seg3d.qualities.mInertia = moinertia ;
        seg3d.qualities.cellQ2d = cellQ2d ;
        seg3d.qualities.readme = struct(...
            'areas', ['area of each cell in squared ' QS.spaceUnits], ...
            'perim', ['perimeter of each cell in ' QS.spaceUnits], ...
            'moment1', ['smaller moment of inertia (about long axis) ', ...
                'in density * squared ' QS.spaceUnits], ...
            'moment2', ['larger moment of inertia (about short axis) ', ...
                'in density * squared ' QS.spaceUnits], ...
            'mInertia', ['centroidal moment of inertia in coordSys ', ...
                'coords Iuu Iuv Ivv, with uv parallel with XY but ', ...
                'offset to centroid, in density * squared ', ...
                QS.spaceUnits], ...
            'ang1', 'angle of long axis in coordSys, in radians', ...
            'ang2', 'angle of long axis in coordSys, in radians', ...
            'cellQ2d', ['quasi-2d cell polygon from embedding space, ', ...
                    'but with cell centroid surface normal rotated ', ...
                    'to be along x axis'], ...
            'nematicTensor', ['n^T * n - 0.5 * [1, 0; 0, 1], ', ...
                    'where n is along long axis'], ...
            'nematicStrength', 'abs(sqrt(MOIEigenvalueRatio)) - 1, strength of elongation' ) ;
            
        % cell statistics 
        % find which are "good" cells to consider
        keep = find(~isnan(ang1) & (areas < maxCellSize) & ...
            moment1 > 0 & moment2 > 0) ;
        seg3d.statistics.keep = keep ;
        seg3d.statistics.maxCellSize = maxCellSize ;
        
        % which coordinate system has been used for segmentation
        coordSys = seg2d.coordSys ;
        disp(['saving segmentation in 3d to: ' outfn])
        save(outfn, 'seg3d', 'coordSys')

        %% Medians of orientation and moment ratio over TIME    
        % Compute cell statistics
        nCells = length(seg3d.qualities.areas) ;
        keep = seg3d.statistics.keep ;
        areas = seg3d.qualities.areas ;

        % Compute mean MOI
        iuu = nanmean(seg3d.qualities.mInertia(keep, 1, 1)) ;
        iuv = - nanmean(seg3d.qualities.mInertia(keep, 1, 2)) ;
        ivv = nanmean(seg3d.qualities.mInertia(keep, 2, 2)) ;
        iii = [ iuu  -iuv ;
             -iuv   ivv ];
        [ eig_vec, eig_val ] = eig(iii);
        meanMoI = iii ;
        meanMoment1 = eig_val(1,1);
        meanMoment2 = eig_val(2,2);
        meanAng1 = atan2( eig_vec(2,1), eig_vec(1,1) );
        meanAng2 = atan2( eig_vec(2,2), eig_vec(1,2) );

        % Compute mean MOI, area-weighted 
        % --> isn't MOI already area-weighted? No, it's weighted oddly. 
        weight = areas(keep) ;
        weights = weight ./ nansum(weight) ;

        % Could rescale MOIs by sqrt(determinant) so each is nematic tensor 
        %   with unit 'size'. How to do this properly? Define Q tensor for each
        %   by Q =  q (n^T n - II/2), where q = (sqrt(I_1/I_2) - 1) is the
        %   magnitude of the anisotropy

        % Compute nematic tensor for each
        mratio = seg3d.qualities.moment2 ./ seg3d.qualities.moment1 ;
        strength = zeros(nCells, 1) ;
        QQ = zeros(nCells, 2, 2) ;
        for qq = 1:nCells
            if ~isempty(intersect(keep, qq))
                tt = mod(seg3d.qualities.ang1(qq), pi) ;
                nn = [cos(tt), sin(tt)] ;
                % Create traceless symmetric matrix using unit vec
                strength(qq) = abs(sqrt(mratio(qq))) - 1 ;
                QQ(qq, :, :) = nn' * nn - 0.5 * [1, 0; 0, 1] ;
            end
        end
        % Take mean shape from nematic tensor
        meanQ = squeeze(mean(strength .* QQ, 1)) ;
        [ eig_vec, eig_val ] = eig(meanQ) ;
        meanQMoment1 = eig_val(1,1);
        meanQMoment2 = eig_val(2,2);
        meanQTheta = atan2( eig_vec(2,2), eig_vec(1,2) );

        % Weight by areas
        meanQW = squeeze(sum(weights .* strength(keep) .* QQ(keep, :, :), 1)) ;
        [ eig_vec, eig_val ] = eig(meanQW) ;
        meanQMoment1Weighted = eig_val(1,1);
        meanQMoment2Weighted = eig_val(2,2);
        meanQThetaWeighted = atan2( eig_vec(2,2), eig_vec(1,2) );

        % STORE NEMATIC INFO IN QUALITIES
        seg3d.qualities.nematicTensor = QQ ;
        seg3d.qualities.nematicStrength = strength ;

        % This is not helpful.
        % i11 = seg3d.qualities.mInertia(keep, 1, 1) ;
        % i12 = seg3d.qualities.mInertia(keep, 1, 2) ;
        % i22 = seg3d.qualities.mInertia(keep, 2, 2) ;
        % sqrtdet = zeros(length(keep), 1) ;
        % for qq = 1:length(keep)
        %     sqrtdet(qq) = sqrt(abs(det([i11(qq), -i12(qq); -i12(qq), i22(qq)]))) ;
        % end
        % iuu = nansum(weights .* i11 ./ sqrtdet) ;
        % iuv = -nansum(weights .* i12 ./ sqrtdet) ;
        % ivv = nansum(weights .* i22 ./ sqrtdet) ;
        % iiiW = [ iuu  -iuv ;
        %      -iuv   ivv ];
        % [ eig_vec, eig_val ] = eig(iiiW);
        % meanMoIWeighted = iiiW ;
        % meanMoment1Weighted = eig_val(1,1);
        % meanMoment2Weighted = eig_val(2,2);
        % meanAng1Weighted = atan2( eig_vec(2,1), eig_vec(1,1) );
        % meanAng2Weighted = atan2( eig_vec(2,2), eig_vec(1,2) );

        % Compute mean MOI, area-weighted, roll off weights for largest cells
        weight(weight > 0.5 * maxCellSize) = maxCellSize - weight(weight > 0.5 * maxCellSize) ;
        weights = weight ./ nansum(weight) ;

        meanQWB = squeeze(sum(weights .* strength(keep) .* QQ(keep, :, :), 1)) ;
        [ eig_vec, eig_val ] = eig(meanQWB) ;
        meanQMoment1WeightBounded = eig_val(1,1);
        meanQMoment2WeightBounded = eig_val(2,2);
        meanQThetaWeightBounded = atan2( eig_vec(2,2), eig_vec(1,2) );

        % iuu = nansum(weights .* i11 ./ sqrtdet) ;
        % iuv = -nansum(weights .* i12 ./ sqrtdet) ;
        % ivv = nansum(weights .* i22 ./ sqrtdet) ;
        % iii = [ iuu  -iuv ;
        %      -iuv   ivv ];
        % [ eig_vec, eig_val ] = eig(iii);
        % meanMoIWeightBounded = iii ;
        % meanMoment1WeightBounded = eig_val(1,1) ;
        % meanMoment2WeightBounded = eig_val(2,2) ;
        % meanAng1WeightBounded = atan2( eig_vec(2,1), eig_vec(1,1) );
        % meanAng2WeightBounded = atan2( eig_vec(2,2), eig_vec(1,2) );

        % Other measures
        ang1 = seg3d.qualities.ang1 ;
        mratio = seg3d.qualities.moment2 ./ seg3d.qualities.moment1 ;
        mratio_principal = mratio(keep) ;
        % c1 = sqrt(mratio_principal(:)) .* cos(2 * ang1(keep)) ;
        % s1 = sqrt(mratio_principal(:)) .* sin(2 * ang1(keep)) ;
        % ar = vecnorm([mean(c1), mean(s1)]) ;
        % theta = 0.5 * atan2(nanmean(s1), nanmean(c1)) ;
        ars = sqrt(mratio_principal(:)) ;
        cos2thetas = cos(2 * ang1(keep)) ;
        sin2thetas = sin(2 * ang1(keep)) ;

        % %% Cell statistics weighted by area
        % % Weight the mean by cell area
        % weight = areas(keep) ;
        % weights = weight ./ nansum(weight) ;
        % c2 = weights .* c1(:) ;
        % s2 = weights .* s1(:) ;
        % ar_weighted = vecnorm([sum(c2), sum(s2)]) ;
        % theta_weighted = 0.5 * atan2(sum(s2), sum(c2)) ;
        % 
        % %% Cell statistics weighted by bounded area
        % weight = areas(keep) ;
        % weight(weight > 0.5 * maxCellSize) = maxCellSize - weight(weight > 0.5 * maxCellSize) ;
        % weights = weight ./ nansum(weight) ;
        % c3 = weights .* c1(:) ;
        % s3 = weights .* s1(:) ;
        % ar_weighted_bounded = vecnorm([sum(c3), sum(s3)]) ;
        % theta_weighted_bounded = 0.5 * atan2(sum(s3), sum(c3)) ;

        % seg3d.statistics.meanAspect = ar ;
        % seg3d.statistics.meanCos2Theta =  ;
        % seg3d.statistics.meanSin2Theta = sin(2*theta) ;
        % seg3d.statistics.meanAspectWeighted = ar_weighted ;
        % seg3d.statistics.meanThetaWeighted = theta_weighted ;
        % seg3d.statistics.meanAspectBoundedWeight = ar_weighted_bounded ;
        % seg3d.statistics.meanThetaBoundedWeight = theta_weighted_bounded ;
        % 
        % seg3d.statistics.aspect25 = prctile(ars, 25.0) ;
        % seg3d.statistics.aspect75 = prctile(ars, 75.0) ;
        % seg3d.statistics.cos2theta25 = prctile(cos2thetas, 25.0) ;
        % seg3d.statistics.cos2theta75 = prctile(cos2thetas, 75.0) ;
        % seg3d.statistics.sin2theta25 = prctile(sin2thetas, 25.0) ;
        % seg3d.statistics.sin2theta75 = prctile(sin2thetas, 75.0) ;

        % The aspect ratio is related to the meanQ 
        % For perfectly aligned cells, meanQ would have a norm(Q) = 0.5, so
        % abs(norm(meanQ)) * 2 = |ar| - 1
        try
            assert(abs(sqrt(abs(4 * det(meanQ))) - 2 * norm(meanQ)) < 1e-7)
        catch
            error('here')
        end

        % Statistics by AP position
        ap_pos = seg3d.cdat.centroids_uv(keep, 1) ;
        xedges = linspace(0, 1, nAPBins + 1) ;
        [mid_ap, mean_c2t, std_c2t] = ...
            binDataMeanStdWeighted(ap_pos, mratio_principal .* cos2thetas, ...
                xedges, weights) ;
        [mid_ap, mean_s2t, std_s2t] = ...
            binDataMeanStdWeighted(ap_pos, mratio_principal .* sin2thetas, ...
                xedges, weights) ;


        %% Statistics by Lobe (between features.folds)
        nU = QS.nU ;
        fold0 = double(folds(tidx, :)) / double(nU) ;
        nLobes = length(fold0(:)) + 1 ;
        foldt = [0; fold0(:); 1] ;
        ap_pos = seg3d.cdat.centroids_uv(keep, 1) ;
        [~, lobes_Q11, lobes_std_Q11] = binDataMeanStdWeighted(ap_pos, ...
            strength(keep) .* squeeze(QQ(keep, 1, 1)), foldt, weights) ;
        [~, lobes_Q12, lobes_std_Q12] = binDataMeanStdWeighted(ap_pos, ...
            strength(keep) .* squeeze(QQ(keep, 1, 2)), foldt, weights) ;
        [~, lobes_Q21, lobes_std_Q21] = binDataMeanStdWeighted(ap_pos, ...
            strength(keep) .* squeeze(QQ(keep, 2, 1)), foldt, weights) ;
        [~, lobes_Q22, lobes_std_Q22] = binDataMeanStdWeighted(ap_pos, ...
            strength(keep) .* squeeze(QQ(keep, 2, 2)), foldt, weights) ;
        
        % standard error on the mean
        lobes_ste_Q11 = lobes_std_Q11 / sqrt(length(keep)) ;
        
        % Check that result is still traceless and symmetric
        assert(all(abs(lobes_Q11 + lobes_Q22) < 1e-7))
        assert(all(abs(lobes_Q12 - lobes_Q12) < 1e-7))

        % Collate lobe information
        meanQLobeAspect = zeros(nLobes, 1) ;
        meanQLobeTheta = zeros(nLobes, 1) ;
        meanQLobeAspectStd = zeros(nLobes, 1) ;
        meanQLobeAspectSte = zeros(nLobes, 1) ;
        meanQLobeThetaStd = zeros(nLobes, 1) ;
        meanQLobeThetaSte = zeros(nLobes, 1) ;
        for lobe = 1:nLobes
            meanQ_lobes{lobe} = [lobes_Q11(lobe), lobes_Q12(lobe); ...
                lobes_Q21(lobe), lobes_Q22(lobe)] ;
            stdQ_lobes{lobe} = [lobes_std_Q11(lobe), lobes_std_Q12(lobe); ...
                lobes_std_Q21(lobe), lobes_std_Q22(lobe)] ;
            if std_ste == 2
                steQ_lobes{lobe} = [lobes_ste_Q11(lobe), lobes_ste_Q12(lobe); ...
                    lobes_ste_Q21(lobe), lobes_ste_Q22(lobe)] ;
            end
            
            % diagonalize this lobeQ
            [ eig_vec, eig_val ] = eig(meanQ_lobes{lobe});
            try
                assert(abs(abs(eig_val(2,2)) - abs(eig_val(1,1))) < 1e-7)
            catch
                error('Something is wrong with traceless or symmetry')
            end
            meanQLobeAspect(lobe) = norm(meanQ_lobes{lobe}) * 2 + 1 ;
            meanQLobeTheta(lobe) = atan2( eig_vec(2,2), eig_vec(1,2) );

            % Uncertainty in average is given by error propagation
            % lambda = 0.5 * [trace +/- sqrt(tr^2 - 4*det)] 
            % Now, the trace is guaranteed to be zero, but not sure that means
            % unc_trace = 0. If not, then we would propagate errors to be

            unc_tr = sqrt(2 * lobes_std_Q11(lobe).^2) ;
            unc_det = sqrt(2 * (lobes_Q11(lobe) * lobes_std_Q11(lobe)).^2 ...
                + 2 * (lobes_Q12(lobe) * lobes_std_Q12(lobe)).^2) ;
            determ = abs(det(meanQ_lobes{lobe})) ;
            % We ask for the uncertainty of the positive eigenvector
            unc_lambda = 0.5 * sqrt(unc_tr.^2 + unc_det.^2 / (determ)) ; 

            % NOTE: |eigenvalue| of symm traceless matrix == norm(matrix)
            meanQLobeAspectStd(lobe) = 2 * unc_lambda ;
            meanQLobeAspectSte(lobe) = 2 * unc_lambda ;
            
            % For angle uncertainty, note that Q is traceless symmetric so
            % we can say Q = [A,B;B,-A]. 
            % Then the eigenvalues are lambda = +/- sqrt(A^2+B^2)
            % Plugging back in allows us to find the eigvects as theta(A,B)
            % so that we can get dtheta(A,B,dA,dB).
            % dtheta = Sqrt[D[theta,A]^2 dA^2 + D[theta,B]^2 dB^2]
            %        = 0.5 * sqrt[ (B^2 dA^2 + A^2 dB^2) / (A^2+B^2)^2 ]
            assert(lobes_Q11(lobe) == - lobes_Q22(lobe))
            assert(lobes_Q21(lobe) == - lobes_Q12(lobe))
            numerator1 = lobes_Q11(lobe)^2 * (lobes_std_Q12(lobe))^2 ;
            numerator2 = lobes_Q12(lobe)^2 * (lobes_std_Q11(lobe))^2 ;
            numerator = numerator1+numerator2 ;
            denominator = (lobes_Q11(lobe)^2 + lobes_Q12(lobe)^2)^2 ;
            meanQLobeThetaStd(lobe) = 0.5 * sqrt(numerator / denominator) ;

            if std_ste == 2
                error('handle here')
            end
        end

        %% Save 
        tmp = seg3d.statistics ;
        seg3d.statistics = struct() ;
        seg3d.statistics.keep = tmp.keep ;
        seg3d.statistics.maxCellSize = tmp.maxCellSize ;
        % Mean tensor stats
        seg3d.statistics.meanMoI = meanMoI ;
        seg3d.statistics.meanQ = meanQ ;
        seg3d.statistics.meanQWeighted = meanQW ;
        seg3d.statistics.meanQWeightBounded = meanQWB ;
        seg3d.statistics.meanMoment1 = meanMoment1 ;
        seg3d.statistics.meanMoment2 = meanMoment2 ;
        seg3d.statistics.meanQMoment1 = meanQMoment1 ;
        seg3d.statistics.meanQMoment2 = meanQMoment2 ;
        seg3d.statistics.meanQMoment1Weighted = meanQMoment1Weighted ;
        seg3d.statistics.meanQMoment2Weighted = meanQMoment2Weighted ;
        seg3d.statistics.meanQMoment1WeightBounded = meanQMoment1WeightBounded ;
        seg3d.statistics.meanQMoment2WeightBounded = meanQMoment2WeightBounded ;

        % mean tensor angles
        seg3d.statistics.meanAng1 = meanAng1 ;
        seg3d.statistics.meanQTheta = meanQTheta ;
        seg3d.statistics.meanQThetaWeighted = meanQThetaWeighted ;
        seg3d.statistics.meanQThetaWeightBounded = meanQThetaWeightBounded ;

        % Mean Q tensor (weightedBounded) for each lobe
        seg3d.statistics.lobes = struct() ;
        seg3d.statistics.lobes.meanQLobes = meanQ_lobes ;
        seg3d.statistics.lobes.stdQLobes = stdQ_lobes ;
        seg3d.statistics.lobes.meanQLobeAspect = meanQLobeAspect ;
        seg3d.statistics.lobes.meanQLobeAspectStd = meanQLobeAspectStd ;
        seg3d.statistics.lobes.meanQLobeTheta = meanQLobeTheta ;

        % Raw distributions
        seg3d.statistics.aspect25 = prctile(ars, 25.0) ;
        seg3d.statistics.aspect75 = prctile(ars, 75.0) ;
        seg3d.statistics.cos2theta25 = prctile(cos2thetas, 25.0) ;
        seg3d.statistics.cos2theta75 = prctile(cos2thetas, 75.0) ;
        seg3d.statistics.sin2theta25 = prctile(sin2thetas, 25.0) ;
        seg3d.statistics.sin2theta75 = prctile(sin2thetas, 75.0) ;
        seg3d.statistics.aspectStd = std(ars) ;
        seg3d.statistics.aspectSte = std(ars) / sqrt(length(keep)) ;
        seg3d.statistics.cos2thetaStd = std(cos2thetas) ;
        seg3d.statistics.cos2thetaSte = std(cos2thetas) / sqrt(length(keep)) ;
        seg3d.statistics.sin2thetaStd = std(sin2thetas) ;
        seg3d.statistics.cos2thetaSte = std(sin2thetas) / sqrt(length(keep)) ;
    

        % AP averaging
        seg3d.statistics.apBins = mid_ap ;
        seg3d.statistics.apCos2Theta = mean_c2t ;
        seg3d.statistics.apSin2Theta = mean_s2t ;
        seg3d.statistics.apCos2ThetaStd = std_c2t ;
        seg3d.statistics.apSin2ThetaStd = std_s2t ;

        disp(['Saving seg3d now with statistics to: ' outfn])
        save(outfn, 'seg3d', 'coordSys')
    else
        % seg3d = QS.loadCurrentSegmentation3D() ;
        disp(['Loading seg3d with statistics from: ' outfn])
        seg3d = load(outfn) ;
        coordSys = seg3d.coordSys ;
        seg3d = seg3d.seg3d ;
        
    end
   
    % Store for later plotting
    keep = seg3d.statistics.keep ;
    mratio = seg3d.qualities.moment2 ./ seg3d.qualities.moment1 ;
    mratio_principal = mratio(keep) ;
    ars = sqrt(mratio_principal(:)) ;
    ang1 = seg3d.qualities.ang1 ; 
    cos2thetas = cos(2 * ang1(keep)) ;
    sin2thetas = sin(2 * ang1(keep)) ;

    c2t_low25 = [c2t_low25, prctile(cos2thetas, 25.0)] ;
    c2t_high75 = [c2t_high75, prctile(cos2thetas, 75.0)] ;
    s2t_low25 = [s2t_low25, prctile(sin2thetas, 25.0)] ;
    s2t_high75 = [s2t_high75, prctile(sin2thetas, 75.0)] ;
    mc2t = [mc2t, cos(2*seg3d.statistics.meanQThetaWeightBounded)] ;
    ms2t = [ms2t, sin(2*seg3d.statistics.meanQThetaWeightBounded)] ;
    c2t_std = [c2t_std, std(cos2thetas)] ;
    s2t_std = [s2t_std, std(sin2thetas)] ;
    c2t_ste = [c2t_ste, std(cos2thetas) / sqrt(length(keep))] ;
    s2t_ste = [s2t_ste, std(sin2thetas) / sqrt(length(keep))] ;
    
    meanQLobeAspects(:, dmy) = seg3d.statistics.lobes.meanQLobeAspect ;
    meanQLobeAspectStds(:, dmy) = seg3d.statistics.lobes.meanQLobeAspectStd ;
    try
        meanQLobeAspectStes(:, dmy) = seg3d.statistics.lobes.meanQLobeAspectSte ;
    catch
        disp('no ste in stats!')
    end
    meanQLobeThetas(:, dmy) = seg3d.statistics.lobes.meanQLobeTheta ;
    mean_c2ts(dmy, :) = seg3d.statistics.apCos2Theta ;
    mean_s2ts(dmy, :) = seg3d.statistics.apSin2Theta ;
    
    mean_qc2ts(dmy, :) = seg3d.statistics.apCos2Theta ;
    mean_qs2ts(dmy, :) = seg3d.statistics.apSin2Theta ;
    std_qc2ts(dmy, :) = seg3d.statistics.apCos2ThetaStd ;
    std_qs2ts(dmy, :) = seg3d.statistics.apSin2ThetaStd ;
    try
        ste_qc2ts(dmy, :) = seg3d.statistics.apCos2ThetaSte ;
        ste_qs2ts(dmy, :) = seg3d.statistics.apSin2ThetaSte ;
    catch
        disp('no ste in stats!')
    end
        
    mean_mratio = [mean_mratio, norm(seg3d.statistics.meanQWeighted) * 2 + 1 ] ;
    median_moiratio = [median_moiratio, median(ars)] ;
    mean_moiratio = [mean_moiratio, ...
        sqrt(seg3d.statistics.meanMoment2/seg3d.statistics.meanMoment1) ] ;

    mratio_low25 = [mratio_low25, seg3d.statistics.aspect25 ] ;
    mratio_high75 = [mratio_high75, seg3d.statistics.aspect75 ] ;
    mratio_std = [mratio_std, seg3d.statistics.aspectStd] ;
    mratio_ste = [mratio_ste, seg3d.statistics.aspectSte] ;
    
    %% Plot this timepoint's segmentation in 3d
    aux_plotCellSegmentation3D(QS, tp, seg3d, imdir, overwriteImages, xyzlims, ~useCorrected)
    
    %% Plot as histogram
    edges = linspace(-1, 1, 100) ;
    edgesAR = linspace(1, 5, 100) ;
    tmp = histcounts(cos(2* seg3d.qualities.ang1), edges) ;
    cos2thetaM(:, dmy) = tmp / sum(tmp) ;
    tmp = histcounts(sin(2* seg3d.qualities.ang1), edges) ;
    sin2thetaM(:, dmy) = tmp / sum(tmp) ;
    tmp = histcounts(ars, edgesAR) ;
    aspectM(:, dmy) = tmp / sum(tmp) ;
    dmy = dmy + 1 ;
end

segSubDir = 'pathlines' ;
mid_ap = seg3d.statistics.apBins ;
aux_plotCellSegmentation3DStats




%% Compare to true segmentation  -- nematic strength and direction for each lobe 
timeList = {timePoints, timePoints(timePoints < t0 + 75 & timePoints > -30)} ;
timeStr = {'', '_tlimit'} ;
for std_ste = 1:2
    for pp = 1:2
        close all 
        fig = figure('units', 'centimeters', 'position', [0,0,figH,figH]) ;
        time2do = timeList{pp} ;
        pInds = ismember(timePoints, time2do) ;

        imfn = fullfile(QS.dir.segmentation, 'pathlines', ...
            ['cell_anisotropy_lobes_signed_COMPARE', timeStr{pp}]) ;
        % Collate results
        kk = 1;
        meanQLobeAspectsTrue = zeros(nLobes, length(time2do)) ;
        meanQLobeAspectStdsTrue = zeros(nLobes, length(time2do)) ;
        meanQLobeAspectStesTrue = zeros(nLobes, length(time2do)) ;
        meanQLobeThetasTrue = meanQLobeAspectsTrue ;
        for tp = time2do
            QS.setTime(tp) ;
            if useCorrected     
                seg3d = QS.getCurrentSegmentation3DCorrected() ; 
            else
                seg3d = QS.getCurrentSegmentation3D() ; 
            end
            meanQLobeAspectsTrue(:, kk) = seg3d.seg3d.statistics.lobes.meanQLobeAspect ;
            meanQLobeAspectStdsTrue(:, kk) = seg3d.seg3d.statistics.lobes.meanQLobeAspectStd ;
            if std_ste == 2
                meanQLobeAspectStesTrue(:, kk) = seg3d.seg3d.statistics.lobes.meanQLobeAspectSte ;
            end
            meanQLobeThetasTrue(:, kk) = seg3d.seg3d.statistics.lobes.meanQLobeTheta ;
            kk = kk + 1;
        end

        for lobe = 1:nLobes

            c2tTrue = cos(2*meanQLobeThetasTrue(lobe, :))  ;
            trueLine = 0.5*c2tTrue .* (squeeze(meanQLobeAspectsTrue(lobe, :)) - 1) ;
            trueUncs = 0.5*c2tTrue .* (squeeze(meanQLobeAspectStdsTrue(lobe, :)) - 1) ;
            

            % subplot(ceil(nLobes * 0.5), 2, lobe)
            c2t = cos(2*meanQLobeThetas(lobe, pInds))  ;
            midline = 0.5 * c2t .* (squeeze(meanQLobeAspects(lobe, pInds)) - 1) ;
            uncs = 0.5 * c2t .* (squeeze(meanQLobeAspectStds(lobe, pInds)) - 1) ;
            
            % Standard errors
            if std_ste == 2
                trueUncEs = 0.5*c2tTrue .* (squeeze(meanQLobeAspectStesTrue(lobe, :)) - 1) ;
                uncEs = 0.5 * c2t .* (squeeze(meanQLobeAspectStes(lobe, pInds)) - 1) ;
            end
            
            timestamps = timePoints - t0 ;
            if contains(QS.timeUnits, 'min')
                timestamps = timestamps / 60 ;
                timeunits = 'hr';
            else
                timeunits = QS.timeUnits ;
            end
            x2 = [timestamps, fliplr(timestamps)] ;
            % fill(x2,[midline-abs(uncs), fliplr(midline+abs(uncs))], ...
            %     colors(lobe, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
            %     'HandleVisibility', 'off');

            hs{lobe} = errorbar(midline, trueLine, trueUncs, trueUncs, uncs, uncs) ;
            hold on;
            %hs{lobe} = plot(midline, tr, '.-', 'color', colors(lobe, :)) ;

            legendentries{lobe} = ['chamber ' num2str(lobe)] ;
        end
        ylims = ylim() ;
        ylim([-max(abs(ylims)), max(abs(ylims))])

        % Mark zero line
        plot([-0.75,1], [-0.75,1], 'k--', 'HandleVisibility','off')
        % Labels
        legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        ylabel('cell anisotropy, $Q_{xx}$',   'interpreter', 'latex')
        xlabel('integrated tissue shear',   'interpreter', 'latex')
        sgtitle('endoderm orientation over time', 'interpreter', 'latex')
        axis equal
        saveas(gcf, [imfn, '.png'])
        saveas(gcf, [imfn, '.pdf'])
    end
end



%% Testing for shape characterization
% % Is the result dependent on the distribution of vertices? No! yay.
% xx = 2 * [0, 1, 1, 1, 1, 0];
% yy = [0, 0, 0.5, 2, 3, 1];
% [ geom, iner, cpmo ] = polygeom( xx, yy ) ;
% plot([xx xx(1)], [yy yy(1)], '.-')
% axis equal; hold on;
% 
% 
% xx2 = 2 * [0, 1, 1, 0];
% yy2 = [0, 0, 3, 1];
% [ geom2, iner2, cpmo2 ] = polygeom( xx2, yy2 ) ;
% plot([xx2 xx2(1)], [yy2 yy2(1)], 'o--')
% axis equal
% 
% moi = [iner(4) -iner(6); -iner(6) iner(5)]
% [tmp, tmp2] = eig(moi)
%
% areas(cid) = geom(1) ;
% perim(cid) = geom(4) ;
% pcentroid(cid) = geom(2:3) ;
% moment1(cid) = cpmo(1) ;
% ang1(cid) = cpmo(2) ;
% moment2(cid) = cpmo(3) ;
% ang2(cid) = cpmo(4) ;
% mratio(cid) = moment2(cid) / moment1(cid) ;
% moinertia(cid, :) = [iner(4) iner(6) iner(5)] ;

end


