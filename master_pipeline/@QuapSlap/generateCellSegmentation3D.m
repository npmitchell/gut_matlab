function generateCellSegmentation3D(QS, options)
% generateCellSegmentation3D(QS, options)
% 
%
% NPMitchell 2021


% Unpack options
timePoints = QS.xp.fileMeta.timePoints ;
maxCellSize = 200 ;  % maximum size allowed to include a cell
overwrite = false ;
debug = false ;
[~, ~, ~, xyzlims] = QS.getXYZLims() ;
t0 = QS.t0set() ;

if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
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
imdir = fullfile(QS.dir.segmentation, 'seg3d', 'images') ;
if ~exist(imdir, 'dir')
    mkdir(imdir)
end

% Setting the current timepoint clears non-timepoint segmentations
medop = [] ;
op_low = [] ;
op_high = [] ;
mean_mratio = [] ;
mratio_low = [] ;
mratio_high = [] ;
for tp = timePoints
    sprintf(['t = ' num2str(tp)])
    QS.setTime(tp)
    tidx = QS.xp.tIdx(tp) ;
    
    outfn = sprintf(QS.fullFileBase.segmentation3d, QS.currentTime) ;
    if ~exist(outfn, 'file') || overwrite || true

        % Obtain the segmentation in 2d 
        seg2d = QS.getCurrentSegmentation2D(options) ;

        % obtain the current cut mesh in APDV coordinates -->
        % 2d coordinates in plane of MESH
        cutMesh = QS.getCurrentSPCutMeshSmRS() ;

        % Normalize the u coordinate
        cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1)) ;

        % tile annular cut mesh
        tileCount = [1,  1] ;
        [ faces, v2D, v3D ] = tileAnnularCutMesh(cutMesh, tileCount) ;

        % Collate cell vertices as XY
        nVertices = length(seg2d.seg2d.Vdat) ;
        % XY = zeros(nVertices, 2) ;
        % for qq = 1:nVertices
        %     XY(qq, :) = [seg2d.seg2d.Vdat(qq).vertxcoord, seg2d.seg2d.Vdat(qq).vertycoord] ;
        % end
        XY = seg2d.seg2d.vdat.v ;

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
        areas = NaN * ones(nCells, 1) ;
        perim = NaN * ones(nCells, 1) ;
        moment1 = NaN * ones(nCells, 1) ;
        ang1 = NaN * ones(nCells, 1) ;
        moment2 = NaN * ones(nCells, 1) ;
        ang2 = NaN * ones(nCells, 1) ;
        moinertia = NaN * ones(nCells, 2, 2) ;
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
        
        %% Compute cell statistics
        keep = find(~isnan(ang1) & (areas < maxCellSize) & ...
            moment1 > 0 & moment2 > 0) ;
        mratio_principal = mratio(keep) ;
        c1 = sqrt(mratio_principal(:)) .* cos(2 * ang1(keep)) ;
        s1 = sqrt(mratio_principal(:)) .* sin(2 * ang1(keep)) ;
        ar = vecnorm([mean(c1), mean(s1)]) ;
        theta = 0.5 * atan2(mean(s1), mean(c1)) ;
        thetas = 0.5 * atan2(s1, c1) ;
        
        %% Cell statistics weighted by area
        % Weight the mean by cell area
        weight = areas(keep) ;
        weights = weight ./ nansum(weight) ;
        c2 = weights .* c1(:) ;
        s2 = weights .* s1(:) ;
        ar_weighted = vecnorm([sum(c2), sum(s2)]) ;
        theta_weighted = 0.5 * atan2(sum(s2), sum(c2)) ;
        
        %% Cell statistics weighted by bounded area
        weight = areas(keep) ;
        weight(weight > 0.5 * maxCellSize) = maxCellSize - weight(weight > 0.5 * maxCellSize) ;
        weights = weight ./ nansum(weight) ;
        c3 = weights .* c1(:) ;
        s3 = weights .* s1(:) ;
        ar_weighted_bounded = vecnorm([sum(c3), sum(s3)]) ;
        theta_weighted_bounded = 0.5 * atan2(sum(s3), sum(c3)) ;
        
        %% Save results stored in struct
        seg2d = QS.currentSegmentation.seg2d ;
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
        seg3d.statistics.keep = keep ;
        seg3d.statistics.maxCellSize = maxCellSize ;
        seg3d.statistics.meanAspect = ar ;
        seg3d.statistics.meanTheta = theta ;
        seg3d.statistics.meanAspectWeighted = ar_weighted ;
        seg3d.statistics.meanThetaWeighted = theta_weighted ;
        seg3d.statistics.meanAspectBoundedWeight = ar_weighted_bounded ;
        seg3d.statistics.meanThetaBoundedWeight = theta_weighted_bounded ;
        %     seg3d.statistics.aspectWeighted25 = prctile(ar_weighted, 25.0) ;
        %     seg3d.statistics.aspectWeighted75 = prctile(ar_weighted, 75.0) ;
        %     seg3d.statistics.thetaWeighted25 = prctile(theta_weighted, 25.0) ;
        %     seg3d.statistics.thetaWeighted75 = prctile(theta_weighted, 75.0) ;
        %     seg3d.statistics.aspectBoundedWeight25 = prctile(ar_weighted_bounded, 25.0) ;
        %     seg3d.statistics.aspectBoundedWeight75 = prctile(ar_weighted_bounded, 75.0) ;
        %     seg3d.statistics.thetaBoundedWeight25 = prctile(theta_weighted_bounded, 25.0) ;
        %     seg3d.statistics.thetaBoundedWeight75 = prctile(theta_weighted_bounded, 75.0) ;
    	seg3d.statistics.aspect25 = prctile(sqrt(mratio_principal(:)), 25.0) ;
    	seg3d.statistics.aspect75 = prctile(sqrt(mratio_principal(:)), 75.0) ;
    	seg3d.statistics.theta25 = prctile(thetas, 25.0) ;
    	seg3d.statistics.theta75 = prctile(thetas, 75.0) ;
        
        % which coordinate system has been used for segmentation
        coordSys = seg2d.coordSys ;
        save(outfn, 'seg3d', 'coordSys')
    else
        seg3d = QS.loadCurrentSegmentation3D() ;
        seg3d = seg3d.seg3d ;
        areas = seg3d.qualities.areas ; 
        ang1 = seg3d.qualities.ang1 ; 
        ang2 = seg3d.qualities.ang2 ; 
        mratio = seg3d.qualities.moment2 ./ seg3d.qualities.moment1 ;
        moinertia = seg3d.qualities.mInertia ;
        c3d = seg3d.vdat.xyzrs ;
        cellCntrd = seg3d.cdat.centroids_3d ;
        mAvg = seg3d.statistics.meanAspectBoundedWeight ;
        thetaAvg = seg3d.statistics.meanThetaBoundedWeight ;
        m25 = seg3d.statistics.aspectWeighted25 ;
        m75 = seg3d.statistics.aspectWeighted75 ;
        m25 = seg3d.statistics.thetaWeighted25 ;
        m75 = seg3d.statistics.thetaWeighted75 ;
        keep = seg3d.statistics.keep ;
    end
    
    %% Medians of orientation and moment ratio over TIME
    % find which are "good" cells to consider
    op_low = [op_low, seg3d.statistics.theta25] ;
    op_high = [op_high, seg3d.statistics.theta75] ;
    medop = [medop, seg3d.statistics.meanThetaBoundedWeight] ;
    
    % Moment ratio
    % Iuu = moinertia(:, 1, 1) ;
    % Ivv = moinertia(:, 2, 2) ;
    % 
    % mratio_sp = sqrt(Ivv(keep) ./ Iuu(keep)) ;
    % mratio_splow = [mratio_low, prctile(mratio_sp, 25.0)] ;
    % mratio_sphigh = [mratio_high, prctile(mratio_sp, 75.0)] ;
    
    mean_mratio = [mean_mratio, seg3d.statistics.meanAspectBoundedWeight ] ;
    mratio_low = [mratio_low, seg3d.statistics.aspect25 ] ;
    mratio_high = [mratio_high, seg3d.statistics.aspect75] ;
    
    %% Draw cells colored by area
    titlestr = ['$t=$' sprintf('%03d', tp-t0) ' ' QS.timeUnits] ;
    
    % Easiest way is to triangulate all polygons using centroids
    % This is fine if the cells are all convex
    faces = seg3d.cdat.polygons ;
    nCells = length(faces) ;
    nVertices = size(seg3d.vdat.uv, 1) ;
    dmyk = 1 ;
    areaV = NaN * zeros(nCells * 7, 1) ;
    ang1V = areaV ;
    mratioV = areaV ;
    oparmV = areaV ;
    IxxV = areaV ;
    IxyV = areaV ;
    IyyV = areaV ;
    ff = zeros(nCells * 7, 3) ;
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
end

%% Plot mean +/- pctile over time
imfn = fullfile(QS.dir.segmentation, 'seg3d', 'cell_anisotropy.png') ;
clf
colors = define_colors() ;
bluecol = colors(1, :) ;
% shade(timePoints - t0, bndlow, timePoints, bndhigh)
x2 = [timePoints - t0, fliplr(timePoints - t0)] ;
fill(x2, [mratio_low, fliplr(mratio_high)], bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
hold on;
plot(timePoints - t0, mean_mratio, '.-')
yyaxis right
fill(x2, [op_low, fliplr(op_high)], bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
hold on;
% shadedErrorBar(timePoints - t0, mean(y,1),std(y),'lineProps','g');
plot(timePoints - t0, medop, '.-')
xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
yyaxis right
ylabel('aspect ratio $\sqrt{I_{1}/I_{2}}$',   'interpreter', 'latex')
yyaxis right
ylabel('nematic orientation $\cos 2\theta$',   'interpreter', 'latex')
title('endoderm orientation over time', 'interpreter', 'latex')
saveas(gcf, imfn)



%% Plot mean +/- pctile over time
imfn = fullfile(QS.dir.segmentation, 'seg3d', 'cell_anisotropy_mratio_log.png') ;
clf
colors = define_colors() ;
bluecol = colors(1, :) ;
% shade(timePoints - t0, bndlow, timePoints, bndhigh)
x2 = [timePoints - t0, fliplr(timePoints - t0)] ;
fill(x2, [mratio_low, fliplr(mratio_high)], bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
hold on;
% shadedErrorBar(timePoints - t0, mean(y,1),std(y),'lineProps','g');
plot(timePoints - t0, mean_mratio, '.-')

xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
ylabel('$\log_{10} \sqrt{I_{\phi\phi}/I_{\zeta\zeta}}$',   'interpreter', 'latex')
title('endoderm orientation over time', 'interpreter', 'latex')
saveas(gcf, imfn)

imfn = fullfile(QS.dir.segmentation, 'seg3d', 'cell_anisotropy_mratio.png') ;
clf
colors = define_colors() ;
bluecol = colors(1, :) ;
% shade(timePoints - t0, bndlow, timePoints, bndhigh)
x2 = [timePoints - t0, fliplr(timePoints - t0)] ;
fill(x2, [10.^(mratio_low), fliplr(10.^(mratio_high))], ...
    bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
hold on;
% shadedErrorBar(timePoints - t0, mean(y,1),std(y),'lineProps','g');
plot(timePoints - t0, 10.^(mean_mratio), '.-')

xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
ylabel('$\sqrt{I_{\phi\phi}/I_{\zeta\zeta}}$',   'interpreter', 'latex')
title('endoderm orientation over time', 'interpreter', 'latex')
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
