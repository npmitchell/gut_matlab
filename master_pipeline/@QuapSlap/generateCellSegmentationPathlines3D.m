function generateCellSegmentationPathlines3D(QS, options)
% generateCellSegmentation3D(QS, options)
% 
%
% NPMitchell 2021


% Unpack options
timePoints = QS.xp.fileMeta.timePoints ;
maxCellSize = 200 ;  % maximum size allowed to include a cell
overwrite = false ;
overwriteImages = false ;
debug = false ;
[~, ~, ~, xyzlims] = QS.getXYZLims() ;
t0 = QS.t0set() ;

if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
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
    QS.setTime(t0) ;
    seg2d = getCurrentSegmentation2D(QS) ;
    cellV0 = seg2d.seg2d.vdat.v ;
    opts = struct('preview', true) ;
    [segVertexPathlines2D, segVertexPathlines3D] = ...
        QS.samplePullbackPathlines(cellV0, opts) ;
    save(cellVertexPathlineFn, 'segVertexPathlines2D', ...
        'segVertexPathlines3D')
else
    load(cellVertexPathlineFn, 'segVertexPathlines2D', ...
        'segVertexPathlines3D')
end

% Load cell segmentation 2d at t0
QS.setTime(t0) ;
seg02d = QS.getCurrentSegmentation2D() ;

% Setting the current timepoint clears non-timepoint segmentations
close all
mc2t = [] ;
ms2t = [] ;
c2t_low = [] ;
c2t_high = [] ;
s2t_low = [] ;
s2t_high = [] ;
mean_mratio = [] ;
mean_moiratio = [] ;
median_moiratio = [] ;
mratio_low = [] ;
mratio_high = [] ;
dmy = 1 ;
cos2thetaM = zeros(99, 1) ;
sin2thetaM = cos2thetaM ;
aspectM = cos2thetaM ;
nAPBins = 20 ;
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
        % replace vertices
        seg2d.seg2d.vdat.v = squeeze(segVertexPathlines2D(tidx, :, :)) ;
        % replace centroids with correct ones
        for pp = 1:length(seg2d.seg2d.cdat.polygons)
            poly = seg2d.seg2d.cdat.polygons{pp} ;
            if ~isempty(poly)
                geom = polygeom( seg2d.seg2d.vdat.v(poly, 1), ...
                    seg2d.seg2d.vdat.v(poly, 2) ) ;
                seg2d.seg2d.cdat.centroid(pp, :) = geom(2:3) ; 
            end
        end
        seg2d.seg2d.Cdat = rmfield(seg2d.seg2d.Cdat, 'centroid') ;
        seg2d.seg2d.Vdat = rmfield(seg2d.seg2d.Vdat, 'vertxcoord') ;
        seg2d.seg2d.Vdat = rmfield(seg2d.seg2d.Vdat, 'vertycoord') ;
        
        % obtain the current cut mesh in APDV coordinates -->
        % 2d coordinates in plane of MESH
        cutMesh = QS.getCurrentSPCutMeshSmRS() ;

        % Normalize the u coordinate
        cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1)) ;

        % tile annular cut mesh
        tileCount = [1,  1] ;
        [ faces, v2D, v3D ] = tileAnnularCutMesh(cutMesh, tileCount) ;

        % Collate cell vertices as XY
        nVertices = length(seg2d.seg2d.vdat) ;
        % XY = zeros(nVertices, 2) ;
        % for qq = 1:nVertices
        %     XY(qq, :) = [seg2d.seg2d.Vdat(qq).vertxcoord, seg2d.seg2d.Vdat(qq).vertycoord] ;
        % end
        XY = seg2d.seg2d.vdat.v ;

        % Collate cell centroids
        nCells = length(seg2d.seg2d.Cdat) ;
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
            cell2d0 = seg2d.seg2d.vdat.v(seg2d.seg2d.cdat.polygons{cid}, :) ;
            cellVtx0 = c3d(seg2d.seg2d.cdat.polygons{cid}, :) ;
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
                    % dchi3d points towards the mapped x axis in embedding space
                    dchi3d = (jac2d3d{cellMeshFaces(cid)} * [1, 0]')' ;
                    dchi3d = dchi3d / vecnorm(dchi3d, 2, 2) ;
                    
                    % Plot the cell in 3d and 2d
                    subplot(2, 2, 1)
                    plot(cell2d0(:, 1), cell2d0(:, 2), '.-');
                    axis equal
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
                    hold on; 

                    subplot(2, 2, 3)
                    plot3(cellVtx0(:, 1), cellVtx0(:, 2), cellVtx0(:, 3), ...
                        '.-')
                    hold on;
                    zplus = cellCntrd(cid, :) + dzeta3d * mean(var(cellVtx0)); 
                    xplus = cellCntrd(cid, :) + dchi3d * mean(var(cellVtx0));
                    plot3dpts([cellCntrd(cid, :); zplus])
                    plot3dpts([cellCntrd(cid, :); xplus])
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
                    [ geom, iner, cpmo ] = polygeom( cell_quasi2d(:, 2), ...
                        cell_quasi2d(:, 3) ) ;

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
        seg3d.vdat.NL = seg2d.seg2d.vdat.NL ;
        seg3d.vdat.BL = seg2d.seg2d.vdat.BL ;
        seg3d.vdat.fourfold = seg2d.seg2d.vdat.fourfold ;
        
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
                    'to be along x axis']) ;
            
        % cell statistics 
        % find which are "good" cells to consider
        keep = find(~isnan(ang1) & (areas < maxCellSize) & ...
            moment1 > 0 & moment2 > 0) ;
        seg3d.statistics.keep = keep ;
        seg3d.statistics.maxCellSize = maxCellSize ;
        
        % which coordinate system has been used for segmentation
        coordSys = seg2d.coordSys ;
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
        % Check that result is still traceless and symmetric
        assert(all(abs(lobes_Q11 + lobes_Q22) < 1e-7))
        assert(all(abs(lobes_Q12 - lobes_Q12) < 1e-7))

        % Collate lobe information
        meanQLobeAspect = zeros(nLobes, 1) ;
        meanQLobeTheta = zeros(nLobes, 1) ;
        meanQLobeAspectStd = zeros(nLobes, 1) ;
        for lobe = 1:nLobes
            meanQ_lobes{lobe} = [lobes_Q11(lobe), lobes_Q12(lobe); ...
                lobes_Q21(lobe), lobes_Q22(lobe)] ;
            stdQ_lobes{lobe} = [lobes_std_Q11(lobe), lobes_std_Q12(lobe); ...
                lobes_std_Q21(lobe), lobes_std_Q22(lobe)] ;

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
            unc_lambda = 0.5 * sqrt(unc_tr.^2 + unc_det.^2 / (determ)) ; 

            % NOTE: |eigenvalue| of symm traceless matrix == norm(matrix)
            meanQLobeAspectStd(lobe) = 2 * unc_lambda ;

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

        % AP averaging
        seg3d.statistics.apBins = mid_ap ;
        seg3d.statistics.apCos2Theta = mean_c2t ;
        seg3d.statistics.apSin2Theta = mean_s2t ;
        seg3d.statistics.apCos2ThetaStd = std_c2t ;
        seg3d.statistics.apSin2ThetaStd = std_s2t ;

        save(outfn, 'seg3d', 'coordSys')
    else
        % seg3d = QS.loadCurrentSegmentation3D() ;
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

    c2t_low = [c2t_low, prctile(cos2thetas, 25.0)] ;
    c2t_high = [c2t_high, prctile(cos2thetas, 75.0)] ;
    s2t_low = [s2t_low, prctile(sin2thetas, 25.0)] ;
    s2t_high = [s2t_high, prctile(sin2thetas, 75.0)] ;
    mc2t = [mc2t, cos(2*seg3d.statistics.meanQThetaWeightBounded)] ;
    ms2t = [ms2t, sin(2*seg3d.statistics.meanQThetaWeightBounded)] ;
    
    meanQLobeAspects(:, dmy) = seg3d.statistics.lobes.meanQLobeAspect ;
    meanQLobeAspectStds(:, dmy) = seg3d.statistics.lobes.meanQLobeAspectStd ;
    meanQLobeThetas(:, dmy) = seg3d.statistics.lobes.meanQLobeTheta ;
    mean_c2ts(dmy, :) = seg3d.statistics.apCos2Theta ;
    mean_s2ts(dmy, :) = seg3d.statistics.apSin2Theta ;
    mean_mratio = [mean_mratio, norm(seg3d.statistics.meanQWeighted) * 2 + 1 ] ;
    median_moiratio = [median_moiratio, median(ars)] ;
    mean_moiratio = [mean_moiratio, ...
        sqrt(seg3d.statistics.meanMoment2/seg3d.statistics.meanMoment1) ] ;

    mratio_low = [mratio_low, seg3d.statistics.aspect25 ] ;
    mratio_high = [mratio_high, seg3d.statistics.aspect75 ] ;
    
    %% Plot this timepoint's segmentation in 3d
    plotCells(QS, tp, seg3d, imdir, overwriteImages, xyzlims)
    
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

%% Define some colors
colors = define_colors() ;
bluecol = colors(1, :) ;
redcol = colors(2, :) ;
yelcol = colors(3, :) ;

%% Plot mean +/- pctile over time
imfn = fullfile(QS.dir.segmentation, 'pathlines', 'cell_anisotropy.png') ;
clf
% shade(timePoints - t0, bndlow, timePoints, bndhigh)
x2 = [timePoints - t0, fliplr(timePoints - t0)] ;
yyaxis left
fill(x2, [mratio_low, fliplr(mratio_high)], bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
hold on;
plot(timePoints - t0, mean_mratio, '.-', 'color', bluecol)
yyaxis right
fill(x2, [c2t_low, fliplr(c2t_high)], redcol, 'facealpha', 0.3, 'edgecolor', 'none');
hold on;
% shadedErrorBar(timePoints - t0, mean(y,1),std(y),'lineProps','g');
plot(timePoints - t0, mc2t, '.-', 'color', redcol)
% addaxis(timePoints - t0, ms2t, '.-', 'color', yelcol)
hold on;
fill(x2, [s2t_low, fliplr(s2t_high)], yelcol, 'facealpha', 0.3, 'edgecolor', 'none');
plot(timePoints - t0, ms2t, '.-', 'color', yelcol)

xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
yyaxis left
ylabel('aspect ratio $\sqrt{I_{1}/I_{2}}$',   'interpreter', 'latex')
yyaxis right
ylabel('nematic orientation $\cos 2\theta$, $\sin2\theta$',   'interpreter', 'latex')
title('endoderm orientation over time', 'interpreter', 'latex')
saveas(gcf, imfn)


%% Plot nematic strength and direction for each lobe 
imfn = fullfile(QS.dir.segmentation, 'pathlines', 'cell_anisotropy_lobes.png') ;
clf
for lobe = 1:nLobes
    subplot(ceil(nLobes * 0.5), 2, lobe)
    midline = squeeze(meanQLobeAspects(lobe, :)) ;
    uncs = squeeze(meanQLobeAspectStds(lobe, :)) ;
    timestamps = timePoints - t0 ;
    if contains(QS.timeUnits, 'min')
        timestamps = timestamps / 60 ;
    end
    x2 = [timestamps, fliplr(timestamps)] ;
    fill(x2,[midline-uncs, fliplr(midline+uncs)], ...
        bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
    hold on;
    plot(timestamps, midline, '.-', 'color', bluecol)
    ylim([1, Inf])
    
    yyaxis right
    % fill(x2, [c2t_low, fliplr(c2t_high)], redcol, 'facealpha', 0.3, 'edgecolor', 'none');
    hold on;
    plot(timestamps, mod(meanQLobeThetas(lobe, :), pi)/pi, '.-')
    % plot(timestamps, sin(2*meanQLobeThetas(lobe, :)), '.-')
    % 'color', redcol)
    ylim([0, 1])
    
    if mod(lobe, 2) == 1 
        yyaxis left
        ylabel('aspect ratio, $a=2||Q|| + 1$',   'interpreter', 'latex')
    else
        yyaxis right
        ylabel('nematic orientation $\theta/\pi$',   'interpreter', 'latex')
    end
    
    % Time label
    if lobe > nLobes - 2
        if contains(QS.timeUnits, 'min')
            xlabel('time [hr]', 'interpreter', 'latex')  
        else
            xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')  
        end
    end
end
sgtitle('endoderm orientation over time', 'interpreter', 'latex')
saveas(gcf, imfn)

%% Plot nematic strength and direction for each lobe 
imfn = fullfile(QS.dir.segmentation, 'pathlines', ...
    'cell_anisotropy_lobes_signed.png') ;
clf
for lobe = 1:nLobes
    % subplot(ceil(nLobes * 0.5), 2, lobe)
    c2t = cos(2*meanQLobeThetas(lobe, :))  ;
    midline = c2t .* (squeeze(meanQLobeAspects(lobe, :)) - 1) ;
    uncs = c2t .* (squeeze(meanQLobeAspectStds(lobe, :)) - 1) ;
    timestamps = timePoints - t0 ;
    if contains(QS.timeUnits, 'min')
        timestamps = timestamps / 60 ;
    end
    x2 = [timestamps, fliplr(timestamps)] ;
    fill(x2,[midline-abs(uncs), fliplr(midline+abs(uncs))], ...
        colors(lobe, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
        'HandleVisibility', 'off');
    hold on;
    hs{lobe} = plot(timestamps, midline, '.-', 'color', colors(lobe, :)) ;
    
    legendentries{lobe} = ['chamber ' num2str(lobe)] ;
end
ylims = ylim() ;
ylim([-max(abs(ylims)), max(abs(ylims))])

% Mark zero line
plot(timestamps, 0*timestamps, 'k--', 'HandleVisibility','off')
% Labels
legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
ylabel('elongation, $2||Q|| \cos 2\theta$',   'interpreter', 'latex')
if contains(QS.timeUnits, 'min')
    xlabel('time [hr]', 'interpreter', 'latex')  
else
    xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')  
end
sgtitle('endoderm orientation over time', 'interpreter', 'latex')

ylim([-2.5, 2.5])
saveas(gcf, imfn)
saveas(gcf, [imfn(1:end-3), 'pdf'])


%% Compare to true segmentation  -- nematic strength and direction for each lobe 


timeList = {timePoints, timePoints(timePoints < t0 + 75 & timePoints > -30)} ;
timeStr = {'', '_tlimit'} ;
for pp = 1:2
    time2do = timeList{pp} ;
    pInds = contains(timePoints, time2do) ;
    
    imfn = fullfile(QS.dir.segmentation, 'pathlines', ...
        ['cell_anisotropy_lobes_signed_COMPARE', timeStr{pp}, '.png']) ;
    clf
    % Collate results
    kk = 1;
    meanQLobeAspectsTrue = zeros(nLobes, length(time2do)) ;
    meanQLobeAspectStdsTrue = zeros(nLobes, length(time2do)) ;
    meanQLobeThetasTrue = meanQLobeAspectsTrue ;
    for tp = time2do
        QS.setTime(tp) ;
        seg3d = QS.getCurrentSegmentation3D() ; 
        meanQLobeAspectsTrue(:, kk) = seg3d.seg3d.statistics.lobes.meanQLobeAspect ;
        meanQLobeAspectStdsTrue(:, kk) = seg3d.seg3d.statistics.lobes.meanQLobeAspectStd ;
        meanQLobeThetasTrue(:, kk) = seg3d.seg3d.statistics.lobes.meanQLobeTheta ;
        kk = kk + 1;
    end

    for lobe = 1:nLobes

        c2tTrue = cos(2*meanQLobeThetasTrue(lobe, :))  ;
        trueLine = c2tTrue .* (squeeze(meanQLobeAspectsTrue(lobe, :)) - 1) ;
        trueUncs = c2tTrue .* (squeeze(meanQLobeAspectStdsTrue(lobe, :)) - 1) ;

        % subplot(ceil(nLobes * 0.5), 2, lobe)
        c2t = cos(2*meanQLobeThetas(lobe, pInds))  ;
        midline = c2t .* (squeeze(meanQLobeAspects(lobe, pInds)) - 1) ;
        uncs = c2t .* (squeeze(meanQLobeAspectStds(lobe, pInds)) - 1) ;
        timestamps = timePoints - t0 ;
        if contains(QS.timeUnits, 'min')
            timestamps = timestamps / 60 ;
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
    plot([-2,2], [-2,2], 'k--', 'HandleVisibility','off')
    % Labels
    legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
    ylabel('elongation, $2||Q|| \cos 2\theta$',   'interpreter', 'latex')
    if contains(QS.timeUnits, 'min')
        xlabel('time [hr]', 'interpreter', 'latex')  
    else
        xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')  
    end
    sgtitle('endoderm orientation over time', 'interpreter', 'latex')

    ylim([-2.5, 2.5])
    saveas(gcf, imfn)
    saveas(gcf, [imfn(1:end-3), 'pdf'])
end


%% Plot histograms
imfn = fullfile(QS.dir.segmentation, 'pathlines', 'cell_anisotropy_hist.png') ;
clf
colormap(cividis)
subplot(2, 2, 1)
imagesc(timePoints - t0, edges, cos2thetaM)
set(gca,'YDir','normal')
caxis([0, 0.05])
xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
ylabel('nematic orientation $\cos2\theta$',   'interpreter', 'latex')
subplot(2, 2, 2)
imagesc(timePoints - t0, edges, sin2thetaM)
set(gca,'YDir','normal')
caxis([0, 0.05])
xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
ylabel('nematic orientation $\sin2\theta$',   'interpreter', 'latex')
subplot(2, 2, 3)
imagesc(timePoints - t0, edgesAR, aspectM)
set(gca,'YDir','normal')
caxis([0, 0.05])
xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
ylabel('aspect ratio $\sqrt{I_{1}/I_{2}}$',   'interpreter', 'latex')
subplot(4, 2, 6)
cb = colorbar('location', 'south') ;
ylabel(cb, 'probability', 'interpreter', 'latex')
caxis([0, 0.05])
axis off
saveas(gcf, imfn) 

%% Aspect ratio distributions, variation of mean shape aspect
imfn = fullfile(QS.dir.segmentation, 'pathlines', 'cell_anisotropy_mratio.png') ;
clf
colors = define_colors() ;
bluecol = colors(1, :) ;
% shade(timePoints - t0, bndlow, timePoints, bndhigh)
x2 = [timePoints - t0, fliplr(timePoints - t0)] ;
fill(x2, [mratio_low, fliplr(mratio_high)], ...
    bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
hold on;
% shadedErrorBar(timePoints - t0, mean(y,1),std(y),'lineProps','g');
plot(timePoints - t0, median_moiratio, '.-')
plot(timePoints - t0, mean_moiratio, '.-')
plot(timePoints - t0, mean_mratio, '.-')
legend({'25-75\%', 'median $\sqrt{I_1/I_2}$', ...
    '$\sqrt{\lambda_1^{\langle I \rangle}/\lambda_2^{\langle I \rangle}}$', ...
    '$2||\langle Q\rangle|| + 1$'}, 'interpreter', 'latex')

xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
ylabel('aspect ratio',   'interpreter', 'latex')
title('endoderm orientation over time', 'interpreter', 'latex')
saveas(gcf, imfn)



%% Plot as a function of AP position and time (kymograph)
imfn = fullfile(QS.dir.segmentation, 'pathlines', 'ap_kymographs_c2t_s2t.png') ;
clf
subplot(1, 2, 1)
imagesc(mid_ap, timePoints-t0, medfilt2(mean_c2ts, [3, 1])) ;
caxis([-10, 10])
colormap(blueblackred)
xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
ylabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
cb = colorbar('location', 'southOutside') ;
ylabel(cb, '$\sqrt{\lambda_1/\lambda_2}\cos 2\theta$', ...
    'interpreter', 'latex')
subplot(1, 2, 2)
imagesc(mid_ap, timePoints-t0, medfilt2(mean_s2ts, [3, 1])) ;
caxis([-10, 10])
colormap(blueblackred)
xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
ylabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
cb = colorbar('location', 'southOutside') ;
ylabel(cb, '$\sqrt{\lambda_1/\lambda_2}\sin 2\theta$', ...
    'interpreter', 'latex')
sgtitle('cell anisotropy kymographs', 'interpreter', 'latex')
saveas(gcf, imfn)

%% Time derivative of filtered image AP position
imfn = fullfile(QS.dir.segmentation, 'pathlines', 'ap_kymographs_dc2t_ds2t.png') ;
clf
cfiltered = medfilt2(mean_c2ts, [3, 1]) ;
[~, dc2t] = gradient(cfiltered) ;
sfiltered = medfilt2(mean_s2ts, [3, 1]) ;
[~, ds2t] = gradient(sfiltered) ;
subplot(1, 2, 1)
imagesc(mid_ap, timePoints-t0, imgaussfilt(medfilt2(dc2t, [3, 1]), 0.5)) ;
caxis([-3, 3])
xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
ylabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
colormap(blueblackred)
cb = colorbar('location', 'southOutside') ;
ylabel(cb, '$\partial_t\left(\sqrt{\lambda_1/\lambda_2} \cos 2\theta\right)$', ...
    'interpreter', 'latex')
subplot(1, 2, 2)
imagesc(mid_ap, timePoints-t0, imgaussfilt(medfilt2(ds2t, [3, 1]), 0.5)) ;
caxis([-3, 3])
xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
ylabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
cb = colorbar('location', 'southOutside') ;
ylabel(cb, '$\partial_t\left(\sqrt{\lambda_1/\lambda_2} \sin 2\theta\right)$', ...
    'interpreter', 'latex')
sgtitle('cell anisotropy kymographs', 'interpreter', 'latex')
saveas(gcf, imfn)


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

function plotCells(QS, tp, seg3d, imdir, overwrite, xyzlims)

    %% Draw cells colored by area
    t0 = QS.t0() ;
    titlestr = ['$t=$' sprintf('%03d', tp-t0) ' ' QS.timeUnits] ;
    
    % Easiest way is to triangulate all polygons using centroids
    % This is fine if the cells are all convex
    faces = seg3d.cdat.polygons ;
    keep = seg3d.statistics.keep ;
    
    areas = seg3d.qualities.areas ; 
    ang1 = seg3d.qualities.ang1 ; 
    ang2 = seg3d.qualities.ang2 ; 
    mratio = seg3d.qualities.moment2 ./ seg3d.qualities.moment1 ;
    moinertia = seg3d.qualities.mInertia ;
    c3d = seg3d.vdat.xyzrs ;
    cellCntrd = seg3d.cdat.centroids_3d ;
    keep = seg3d.statistics.keep ;
    
    nCells = length(faces) ;
    nVertices = size(seg3d.vdat.uv, 1) ;
    dmyk = 1 ;
    ff = zeros(nCells * 7, 3) ;
    areaV = NaN * zeros(nCells * 7, 1) ;
    ang1V = areaV ;
    mratioV = areaV ;
    oparmV = areaV ;
    IxxV = areaV ;
    IxyV = areaV ;
    IyyV = areaV ;
    for cid = 1:nCells
        if ismember(cid, keep)
            face = faces{cid} ;
            if ~isempty(face)
                for vid = 1:length(face)
                    if vid < length(face)
                        addface = [face(vid), face(vid+1), nVertices + cid] ;
                    else
                        addface = [face(vid), face(1), nVertices + cid] ;
                    end
                    ff(dmyk, :) = addface ;
                    areaV(dmyk) = areas(cid) ;
                    ang1V(dmyk) = ang1(cid) ;
                    mratioV(dmyk) = mratio(cid) ;
                    IxxV(dmyk) = moinertia(cid, 1) ;
                    IxyV(dmyk) = moinertia(cid, 2) ;
                    IyyV(dmyk) = moinertia(cid, 3) ;
                    oparmV(dmyk) = (mratio(cid) - 1) * cos(2*ang1(cid)) ;
                    dmyk = dmyk + 1 ;
                end
            end
        end
    end
    ff = ff(1:dmyk-1, :) ;
    areaV = areaV(1:dmyk-1) ;
    ang1V = ang1V(1:dmyk-1) ;
    IxxV = IxxV(1:dmyk-1) ;
    IxyV = IxyV(1:dmyk-1) ;
    IyyV = IyyV(1:dmyk-1) ;
    mratioV = mratioV(1:dmyk-1) ;
    oparmV = oparmV(1:dmyk-1) ;
    
    % Extend vertices to include centroids
    vv = [c3d; cellCntrd] ;
    
    %% Color segmentation by area
    imfn = fullfile(imdir, sprintf('cellseg3d_area_%06d.png', tp)) ;
    if ~exist(imfn, 'file') || overwrite
        clf
        patch('Faces',ff,'Vertices',vv,...
            'FaceVertexCData',areaV(:),'FaceColor','flat', ...
            'Edgecolor', 'none');
        cb = colorbar ;
        ylabel(cb, ['area [' QS.spaceUnits '$^2$]'],   'interpreter', 'latex')
        caxis([0, nanmean(areas) + 3*nanstd(areas)])
        axis equal
        view(0,0)
        xlim(xyzlims(1, :))
        ylim(xyzlims(2, :))
        zlim(xyzlims(3, :))
        colormap viridis
        xlabel(['ap position [' QS.spaceUnits ']'], 'interpreter', 'latex')
        ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
        zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
        title(titlestr, 'interpreter', 'latex')
        saveas(gcf, imfn)
    end
    
    %% Color segmentation by moi ratio -- sp coord sys and principal coords
    imfns = {fullfile(imdir, sprintf('cellseg3d_mratioSP_log_%06d.png', tp)), ...
        fullfile(imdir, sprintf('cellseg3d_mratioSP_%06d.png', tp))};
    for tmp = 1:2
        if ~exist(imfns{tmp}, 'file') || overwrite 
            clf
            if tmp == 1
                patch('Faces',ff,'Vertices',vv,...
                    'FaceVertexCData',real(log10(sqrt(IyyV ./ IxxV))), ...
                    'FaceColor','flat', ...
                    'Edgecolor', 'none');
                cb = colorbar ;
                ylabel(cb, '$\log_{10} \sqrt{I_{\phi\phi}/I_{\zeta\zeta}}$',   'interpreter', 'latex')
                caxis([-1, 1])
                bbr256 = blueblackred ;
            else
                patch('Faces',ff,'Vertices',vv,...
                    'FaceVertexCData',real(sqrt(IyyV ./ IxxV)), ...
                    'FaceColor','flat', ...
                    'Edgecolor', 'none');
                cb = colorbar ;
                ylabel(cb, '$\sqrt{I_{\phi\phi}/I_{\zeta\zeta}}$',   'interpreter', 'latex')
                caxis([0, 2])
            end
            colormap(bbr256)
            axis equal
            view(0,0)
            xlim(xyzlims(1, :))
            ylim(xyzlims(2, :))
            zlim(xyzlims(3, :))
            xlabel(['ap position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            title(titlestr, 'interpreter', 'latex')
            saveas(gcf, imfns{tmp})
        end
    end
    
    %% Order parameter    
    imfn = fullfile(imdir, sprintf('cellseg3d_order_%06d.png', tp)) ;
    if ~exist(imfn, 'file') || overwrite
        patch('Faces',ff,'Vertices',vv,...
            'FaceVertexCData',cos(2*ang1V(:)),'FaceColor','flat', ...
            'Edgecolor', 'none');
        cb = colorbar ;
        ylabel(cb, '$\cos 2\theta$',   'interpreter', 'latex')
        caxis([-1, 1])
        colormap(blueblackred)
        axis equal
        view(0,0)
        xlabel(['ap position [' QS.spaceUnits ']'], 'interpreter', 'latex')
        ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
        zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
        title(titlestr, 'interpreter', 'latex')
        saveas(gcf, imfn)
    end
    
    
    %% Draw bonds
    imfn = fullfile(imdir, sprintf('cellseg3d_bonds_full_%06d.png', tp)) ;
    if ~exist(imfn, 'file') || overwrite 
        for fullID = [0, 1] 
            Xs = zeros(4*length(c3d(:, 1)), 1) ;
            Ys = Xs ;
            Zs = Xs ;
            Us = Xs ;
            Vs = Xs ;
            Ws = Xs ;
            dmyk = 1 ;
            for qq = 1:nVertices
                for id = seg3d.vdat.NL(qq, :)
                    if id > 0
                        Xs(dmyk) = c3d(qq, 1) ;
                        Ys(dmyk) = c3d(qq, 2) ; 
                        Zs(dmyk) = c3d(qq, 3) ;
                        Us(dmyk) = c3d(id, 1) - c3d(qq, 1) ;
                        Vs(dmyk) = c3d(id, 2) - c3d(qq, 2) ; 
                        Ws(dmyk) = c3d(id, 3) - c3d(qq, 3) ;
                        dmyk = dmyk + 1 ;
                    end
                end
            end
            plot3(c3d(:, 1), c3d(:, 2), c3d(:, 3), '.')
            hold on;
            q = quiver3(Xs,Ys,Zs, Us, Vs, Ws, 0, 'color', [ 0.8500    0.3250    0.0980]);
            axis equal
            q.ShowArrowHead = 'off';
            [~, ~, ~, xyzlims] = QS.getXYZLims() ;
            xlim(xyzlims(1, :))
            if fullID
                ylim(xyzlims(2, :))
                imfn = fullfile(imdir, sprintf('cellseg3d_bonds_full_%06d.png', tp)) ;
            else
                ylim([xyzlims(2, 1), 0])
                imfn = fullfile(imdir, sprintf('cellseg3d_bonds_%06d.png', tp)) ;    
            end
            zlim(xyzlims(3, :))
            view(0,0)
            xlabel(['ap position, [' QS.spaceUnits ']'], 'Interpreter', 'latex')
            ylabel(['lateral position, [' QS.spaceUnits ']'], 'Interpreter', 'latex')
            zlabel(['dv position, [' QS.spaceUnits ']'], 'Interpreter', 'latex')
            title(titlestr, 'interpreter', 'latex')
            saveas(gcf, imfn) 
            clf
        end
    end
    
    %% Statistics
    statsfn = fullfile(imdir, sprintf('stats_%06d.png', tp)) ;
    if ~exist(statsfn, 'file') || overwrite || true
        clf
        % plot(sqrt(i11), areas(keep), '.') ; hold on;
        % plot(sqrt(i22), areas(keep), '.') ; hold on;
        subplot(2, 1, 1)
        plot(areas(keep), sqrt(seg3d.qualities.moment2(keep)), '.') ; hold on;
        plot(areas(keep), sqrt(seg3d.qualities.moment1(keep)), '.') ; 
        xlabel(['area [' QS.spaceUnits '$^2$]'], 'interpreter', 'latex')
        ylabel('$\sqrt{\lambda_2}, \sqrt{\lambda_1}$', 'interpreter', 'latex')
        subplot(2, 2, 3)
        ars_tmp = sqrt(seg3d.qualities.moment2(keep)./seg3d.qualities.moment1(keep)) ;
        plot(areas(keep), ars_tmp, '.') ; 
        xlabel(['area [' QS.spaceUnits '$^2$]'], 'interpreter', 'latex')
        ylabel('$\sqrt{\lambda_2 / \lambda_1}$', 'interpreter', 'latex')
        % Fit to line to see if there is variation
        [cc, SS] = polyfit(areas(keep), ars_tmp, 1) ;
        uncs = sqrt(abs(SS.R)) / SS.df ;
        title(['$\sqrt{\lambda_2 / \lambda_1} = ($' ...
            num2str(round(cc(1), 1, 'significant')) '$\pm$' ...
            num2str(round(uncs(1,1), 1, 'significant')) '$)A + $' ...
            num2str(round(cc(2), 3, 'significant')) '$\pm$' ...
            num2str(round(uncs(2,2), 1, 'significant'))], ...
            'interpreter', 'latex')
        
        subplot(2, 2, 4)
        plot(areas(keep), ars_tmp, '.') ;
        corrs = corrcoef(areas(keep), ars_tmp) ; 
        title(['correlation = ' ...
            num2str(round(corrs(1, 2), 2, 'significant'))], ...
            'interpreter', 'latex')
        
        ylim([0.5, 5])
        xlabel(['area [' QS.spaceUnits '$^2$]'], 'interpreter', 'latex')
        ylabel('$\sqrt{\lambda_2 / \lambda_1}$', 'interpreter', 'latex')
        sgtitle(titlestr, 'interpreter', 'latex')
        saveas(gcf, statsfn)
    end
end

