function generateCellSegmentation3D(QS, options)
% generateCellSegmentation3D(QS, options)
%
% Define Q tensor for each cell by
%   Q =  q (n^T n - II/2), where q = (sqrt(I_1/I_2) - 1) is the
%   magnitude of the anisotropy
% 
%
% The Q tensor is related to the aspect ratio of cells via:
%   abs(norm(Q)) * 2 = |ar| - 1
% The factor of two arises on the LHS since the 
% norm(n'*n - Identity*0.5) = 0.5 for any unit vector n, and an isotropic
% cell will have norm(Q) = 0.  
%
% NPMitchell 2021

% Colors for cosine, sine:
tmp = define_colors ;
% CScolors = [0.90    0.55    0.55 ;
%         0.62    0.76    0.84 ] ;
CScolors = tmp([9,5], :) ;


% Unpack options
timePoints = QS.xp.fileMeta.timePoints ;
maxCellSize = Inf ;  % maximum size allowed to include a cell
overwrite = false ;
corrected = false ;
coordSys = 'spsme' ;
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
if isfield(options, 'coordSys')
    coordSys = options.coordSys ;
    coordSys = lower(erase(coordSys, '_')) ;
    % If the coordinate system is (s,phi) relaxed aspect ratio smoothed
    % extended, correct for variants of naming convention.
    if strcmpi(coordSys, 'rsme') || strcmpi(coordSys, 'rspsme')
        coordSys = 'sprsme' ;
    end
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
if isfield(options, 'correctedSegmentation')
    corrected = options.correctedSegmentation ;
end

% Directory preparation
if corrected
    segDir = fullfile(QS.dir.segmentation, 'seg3d_corrected') ;
    if ~exist(segDir, 'dir')
        mkdir(segDir)
    end
    imdir = fullfile(QS.dir.segmentation, 'seg3d_corrected', 'images') ;
else
    imdir = fullfile(QS.dir.segmentation, 'seg3d', 'images') ;
end
if ~exist(imdir, 'dir')
    mkdir(imdir)
end

% Prep for dicing data by lobes
features = QS.getFeatures() ;
folds = features.folds ;
nLobes = size(folds, 2) + 1 ;

% Setting the current timepoint clears non-timepoint segmentations
close all
mc2t = [] ;
ms2t = [] ;
c2t_low25 = [] ;
c2t_high75 = [] ;
s2t_low25 = [] ;
s2t_high75 = [] ;
c2t_std = [] ;
s2t_std = [] ;
c2t_ste = [] ;
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
meanQLobeAspects = zeros(nLobes, length(timePoints)) ;
meanQLobeAspectStds = zeros(nLobes, length(timePoints)) ;
meanQLobeThetas = zeros(nLobes, length(timePoints)) ;
for tp = timePoints
    sprintf(['t = ' num2str(tp)])
    QS.setTime(tp)
    tidx = QS.xp.tIdx(tp) ;
    
    if corrected
        outfn = sprintf(QS.fullFileBase.segmentation3dCorrected, QS.currentTime) ;
    else
        outfn = sprintf(QS.fullFileBase.segmentation3d, QS.currentTime) ;
    end
    if ~exist(outfn, 'file') || overwrite

        % Obtain the segmentation in 2d 
        if corrected
            seg2d = QS.getCurrentSegmentation2DCorrected(options) ;
        else
            seg2d = QS.getCurrentSegmentation2D(options) ;
        end
        
        % obtain the current cut mesh in APDV coordinates -->
        % 2d coordinates in plane of MESH
        if strcmp(coordSys, 'spsme') || strcmp(coordSys, 'sprsme') 
            cutMesh = QS.getCurrentSPCutMeshSmRS() ;

            % Normalize the u coordinate
            cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1)) ;
        else
            error(['Code for this coordSys: ' coordSys])
        end

        % tile annular cut mesh
        tileCount = [1,  1] ;
        [ faces, v2D, v3D ] = tileAnnularCutMesh(cutMesh, tileCount) ;

        % Collate cell vertices as XY
        % XY = zeros(nVertices, 2) ;
        % for qq = 1:nVertices
        %     XY(qq, :) = [seg2d.seg2d.Vdat(qq).vertxcoord, seg2d.seg2d.Vdat(qq).vertycoord] ;
        % end
        if corrected
            % corrected segmentation stores polygons as vertices
            polygons = seg2d.seg2d.cdat.polygons ;
            XY = [] ;
            pgonIDs = [] ;
            for pp = 1:length(polygons)
                pgon = polygons{pp} ;
                if ~isempty(pgon)
                    pgon = pgon(:, [2, 1]) ;
                    if isempty(XY)
                        XY = pgon ;
                    else
                        XY = [XY; pgon] ;
                    end
                    nn = length(pgonIDs) ;
                    pgonIDs(nn+1:nn+size(pgon, 1)) = pp ;
                end
            end
        else
            % seg2d stores a network of vertices and bonds
            XY = seg2d.seg2d.vdat.v ;
        end
        
        % Collate cell centroids
        if corrected
            centroids = seg2d.seg2d.cdat.centroid ;
            nCells = size(centroids, 1) ;
        else
            nCells = length(seg2d.seg2d.Cdat) ;
            centroids = zeros(nCells, 2) ;
            for qq = 1:nCells
                centroids(qq, :) = seg2d.seg2d.Cdat(qq).centroid.coord ;
            end    
        end
        
        % Check that we are interpolating over the same domain of
        % parameterization
        if strcmpi(seg2d.coordSys, 'spsme') || ...
            strcmpi(seg2d.coordSys, 'sprsme')
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

        if corrected
            poly3d = cell(length(polygons), 1) ;
            polygonVtxIds = cell(length(polygons), 1) ;
            cellVtx2D = cell(length(polygons), 1) ;
            for pp = 1:length(polygons)
                poly3d{pp} = c3d(pgonIDs == pp, :) ;
                polygonVtxIds{pp} = find(pgonIDs == pp) ;
                cellVtx2D{pp} = uv(polygonVtxIds{pp}, :) ;
            end

            % Obtain cell vertices in 3d
            % cell2d0 = seg2d.seg2d.cdat.polygons{cid} ;
            % cellVtx2D = uv(pgonIDs == cid, :) ;
            % cellVtx0 = c3d(pgonIDs == cid, :) ;
            [c3d, cellCntrd3d, areas, perim, moment1, ang1, ...
                moment2, ang2, moinertia, cellQ2d, ...
                cellMeshFaces, vertexMeshFaces, cellNormals ] = ...
                polygonNetwork3dMeasurements(faces, v3D, v2D, uv, polygonVtxIds, cntrds) ;
        else
            % Obtain cell vertices in 3d
            % original segmentation stores indices into vertices for
            % polygon shapes
            % cell2d0 = seg2d.seg2d.vdat.v(seg2d.seg2d.cdat.polygons{cid}, :) ;
            polygonVtxIds = seg2d.seg2d.cdat.polygons ;
            for pp = 1:length(polygonVtxIds)
                cellVtx2D{pp} = uv(polygonVtxIds{cid}, :) ;
                % cellVtx0 = c3d(seg2d.seg2d.cdat.polygons{cid}, :) ;
            end
            [c3d, cellCntrd3d, areas, perim, moment1, ang1, ...
                moment2, ang2, moinertia, cellQ2d,  cellMeshFaces, vertexMeshFaces] = ...
                polygonNetwork3dMeasurements(faces, v3D, v2D, uv, polygonVtxIds, cntrds) ;
        end
        
        %% Save results stored in struct
        if corrected
            seg2d = QS.currentSegmentation.seg2dCorrected ;
        else
            seg2d = QS.currentSegmentation.seg2d ;
        end
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
        if ~corrected
            seg3d.vdat.NL = seg2d.seg2d.vdat.NL ;
            seg3d.vdat.BL = seg2d.seg2d.vdat.BL ;
            seg3d.vdat.fourfold = seg2d.seg2d.vdat.fourfold ;
        end
        
        % cell data in 3d
        seg3d.cdat.centroids_uv = cntrds ;
        seg3d.cdat.centroids_3d = cellCntrd3d ;
        seg3d.cdat.meshFaces = cellMeshFaces ;
        seg3d.cdat.normals = cellNormals ;
        if corrected
            seg3d.cdat.polygons = poly3d ;
            seg3d.vdat.vertexPolygonMemberID = pgonIDs ;
            seg3d.cdat.polygonVertexID = cell(length(polygons), 1) ;
            dmyk = 1 ;
            for pp = 1:length(polygons)
                if length(polygons{pp}) > 1
                    seg3d.cdat.polygonVertexID{pp} = dmyk:(dmyk+length(polygons{pp})-1) ;
                    dmyk = dmyk + length(polygons{pp}) ;
                end
            end
        else
            seg3d.cdat.polygons = seg2d.seg2d.cdat.polygons ;
        end
        
        % Cell nematic tensor Q and strength of elongation qstrength
        % --> add to struct after filtering
        % mratio = seg3d.qualities.moment2 ./ seg3d.qualities.moment1 ;
        % strength = zeros(nCells, 1) ;
        % QQ = zeros(nCells, 2, 2) ;
        % for qq = 1:nCells
        %     if ~isempty(intersect(keep, qq))
        %         tt = mod(seg3d.qualities.ang1(qq), pi) ;
        %         nn = [cos(tt), sin(tt)] ;
        %         % Create traceless symmetric matrix using unit vec
        %         strength(qq) = abs(sqrt(mratio(qq))) - 1 ;
        %         QQ(qq, :, :) = nn' * nn - 0.5 * [1, 0; 0, 1] ;
        %     end
        % end
        
        % cell qualities, including Q= q*(n^T n - II/2)
        seg3d.qualities = struct() ;
        seg3d.qualities.areas = areas ;
        seg3d.qualities.perim = perim ;
        seg3d.qualities.moment1 = moment1 ;
        seg3d.qualities.moment2 = moment2 ;
        seg3d.qualities.ang1 = ang1 ;
        seg3d.qualities.ang2 = ang2 ;
        seg3d.qualities.mInertia = moinertia ;
        seg3d.qualities.cellQ2d = cellQ2d ;
        % seg3d.qualities.nematicTensor = QQ ;   % added after filtering
        % seg3d.qualities.nematicStrength = strength ; % added after filtering
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
            'ang2', 'angle of short axis in coordSys, in radians', ...
            'cellQ2d', ['quasi-2d cell polygon from embedding space, ', ...
                    'but with cell centroid surface normal rotated ', ...
                    'to be along x axis'], ...
            'nematicTensor', 'n^T * n - 0.5 * [1, 0; 0, 1], where n is along long axis', ...
            'nematicStrength', 'abs(sqrt(MOIEigenvalueRatio)) - 1, strength of elongation' ) ;
            
        % cell statistics 
        % find which are "good" cells to consider
        keep = find(~isnan(ang1) & (areas < maxCellSize) & ...
            moment1 > 0 & moment2 > 0) ;
        seg3d.statistics.keep = keep ;
        seg3d.statistics.maxCellSize = maxCellSize ;
         
        % which coordinate system has been used for segmentation
        coordSys = seg2d.coordSys ;
        save(outfn, 'seg3d', 'coordSys')
    else
        if corrected
            seg3d = QS.loadCurrentSegmentation3DCorrected() ;
        else
            seg3d = QS.loadCurrentSegmentation3D() ;
        end
        coordSys = seg3d.coordSys ;
        seg3d = seg3d.seg3d ;
        
    end
    
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
    allweights = zeros(size(areas)) ;
    allweights(keep) = weights ;
    
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
    seg3d.statistics.weights = allweights ;
    
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
    allweightsBounded = zeros(size(areas)) ;
    allweightsBounded(keep) = weights ;
    seg3d.statistics.weightsBounded = allweightsBounded ;
    
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

    c2t_low25 = [c2t_low25, prctile(cos2thetas, 25.0)] ;
    c2t_high75 = [c2t_high75, prctile(cos2thetas, 75.0)] ;
    s2t_low25 = [s2t_low25, prctile(sin2thetas, 25.0)] ;
    s2t_high75 = [s2t_high75, prctile(sin2thetas, 75.0)] ;
    mc2t = [mc2t, cos(2*meanQThetaWeightBounded)] ;
    ms2t = [ms2t, sin(2*meanQThetaWeightBounded)] ;
    c2t_std = [c2t_std, std(cos2thetas)] ;
    s2t_std = [s2t_std, std(sin2thetas)] ;
    c2t_ste = [c2t_ste, std(cos2thetas) / sqrt(length(keep))] ;
    s2t_ste = [s2t_ste, std(sin2thetas) / sqrt(length(keep))] ;
    
    % The aspect ratio is related to the meanQ 
    % For perfectly aligned cells, meanQ would have a norm(Q) = 0.5, so
    % abs(norm(meanQ)) * 2 = |ar| - 1
    try
        assert(abs(sqrt(abs(4 * det(meanQ))) - 2 * norm(meanQ)) < 1e-7)
    catch
        error('here')
    end
    % Multiplying by two here since norm gives 1/2*(strength_Q)
    mean_mratio = [mean_mratio, norm(meanQW) * 2 + 1 ] ;
    median_moiratio = [median_moiratio, median(ars)] ;
    mean_moiratio = [mean_moiratio, ...
        sqrt(meanMoment2/meanMoment1) ] ;
    mratio_low25 = [mratio_low25, prctile(ars, 25.0)] ;
    mratio_high75 = [mratio_high75, prctile(ars, 75.0) ] ;
    mratio_std = [mratio_std, std(ars)] ;
    mratio_ste = [mratio_ste, std(ars) / sqrt(length(keep))] ;
    
    % Statistics by AP position
    ap_pos = seg3d.cdat.centroids_uv(keep, 1) ;
    xedges = linspace(0, 1, nAPBins + 1) ;
    [mid_ap, mean_qc2t_ap, std_qc2t_ap] = ...
        binDataMeanStdWeighted(ap_pos, strength(keep) .* cos2thetas, ...
            xedges, weights) ;
    [mid_ap, mean_qs2t_ap, std_qs2t_ap] = ...
        binDataMeanStdWeighted(ap_pos, strength(keep) .* sin2thetas, ...
            xedges, weights) ;
        
    %% Statistics by Lobe (between features.folds)
    nU = QS.nU ;
    fold0 = double(folds(tidx, :)) / double(nU) ;
    nLobes = length(fold0(:)) + 1 ;
    foldt = [0; fold0(:); 1] ;
    ap_pos = seg3d.cdat.centroids_uv(keep, 1) ;
    
    % Statistics by Lobe -- ars.*cos(2theta), ars.*sin(2theta)
    [~, mean_qc2t_lobes, std_qc2t_lobes] = binDataMeanStdWeighted(ap_pos, ...
        strength(keep) .* cos(2 * ang1(keep)), foldt, weights) ;
    [~, mean_qs2t_lobes, std_qs2t_lobes] = binDataMeanStdWeighted(ap_pos, ...
        strength(keep) .* sin(2 * ang1(keep)), foldt, weights) ;
    
    % Other measures
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
    
    % Store for later plotting
    meanQLobeAspects(:, dmy) = meanQLobeAspect ;
    meanQLobeAspectStds(:, dmy) = meanQLobeAspectStd ;
    meanQLobeThetas(:, dmy) = meanQLobeTheta ;
    mean_qc2ts(dmy, :) = mean_qc2t_ap ;
    mean_qs2ts(dmy, :) = mean_qs2t_ap ;
    std_qc2ts(dmy, :) = std_qc2t_ap ;
    std_qs2ts(dmy, :) = std_qs2t_ap ;
    
    
    %% Save 
    tmp = seg3d.statistics ;
    seg3d.statistics = struct() ;
    seg3d.statistics.keep = tmp.keep ;
    seg3d.statistics.maxCellSize = tmp.maxCellSize ;
    seg3d.statistics.weights = tmp.weights ;
    seg3d.statistics.weights = tmp.weightsBounded ;
    
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
    seg3d.statistics.lobes.meanQCos2Theta = mean_qc2t_lobes ;
    seg3d.statistics.lobes.stdQCos2Theta = std_qc2t_lobes ;
    seg3d.statistics.lobes.meanQSin2Theta = mean_qs2t_lobes ;
    seg3d.statistics.lobes.stdQSin2Theta = std_qs2t_lobes ;
    
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
    seg3d.statistics.apCos2Theta = mean_qc2t_ap ;
    seg3d.statistics.apSin2Theta = mean_qs2t_ap ;
    seg3d.statistics.apCos2ThetaStd = std_qc2t_ap ;
    seg3d.statistics.apSin2ThetaStd = std_qs2t_ap ;
    save(outfn, 'seg3d', 'coordSys')
    
    %% Plot this timepoint's segmentation in 3d
    aux_plotCellSegmentation3D(QS, tp, seg3d, imdir, overwriteImages, xyzlims, corrected)
    
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

% FIgure size in cm
figW = 18 ;
figH = 9 ;

%% Plot mean +/- pctile over time
if corrected
    segSubDir = 'seg3d_corrected' ;
else
    segSubDir = 'seg3d' ;
end


aux_plotCellSegmentation3DStats


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

