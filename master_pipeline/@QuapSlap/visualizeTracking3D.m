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

xwidth = 16 ; % cm
ywidth = 10 ; % cm

% Glean options from texturePatch run
metafn = fullfile(QS.dir.texturePatchIm, 'overlays', 'metadat.mat') ;
load(metafn, 'metadat', 'Options')
metadat.subdir = 'overlays';

% Default options
timePoints = QS.xp.fileMeta.timePoints ;
subdir = 'labeled_groundTruth';
method = 'segmentation' ;
overwrite = false ;
preview = false ;
coordSys = 'spsme' ;
blackFigure = false ;
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'blackFigure')
    blackFigure = options.blackFigure ;
end
if isfield(options, 'method')
    method = options.method ;
end
if isfield(options, 'buffer') || isfield(options, 'bufferXYZ')
    if isfield(options, 'buffer')
        bufferXYZ = options.buffer ;
    elseif isfield(options, 'bufferXYZ')
        bufferXYZ = options.bufferXYZ ;
    end
    % interpret bufferXYZ if length-1 float or length(3) array
    if numel(bufferXYZ) == 1
        bufferXYZ = [-bufferXYZ, bufferXYZ] ;
    elseif numel(bufferXYZ) == 3
        bufferXYZ = [-bufferXYZ(1), bufferXYZ(1); ...
            -bufferXYZ(2), bufferXYZ(2); ...
            -bufferXYZ(3), bufferXYZ(3)] ;
    else
        assert(all(size(bufferXYZ) == [3, 2]))
    end
    % Get XYZlimits
    [~, ~, ~, xyzlims] = QS.getXYZLims() ;
    xyzlims = xyzlims + bufferXYZ ;
    metadat.xyzlim = xyzlims ;
else
    xyzlims = metadat.xyzlim ;
end
if isfield(options, 'subdir')
    subdir = options.subdir ;
end
if isfield(options, 'coordSys')
    coordSys = options.coordSys ;
end


            
% Plot each timepoint in 3d            
for tidx = 1:length(timePoints)
    tp = timePoints(tidx) ;
    QS.setTime(tp) ;
    if strcmpi( method, 'segmentation')
        seg2dFn = fullfile(QS.dir.tracking, subdir, sprintf('seg2d_%06d.mat', tp)) ;
        outFigFn = fullfile(QS.dir.tracking, subdir, 'embeddingFigure', sprintf('trackSeg3d_%06d.png', tp)) ;
        outImFn = fullfile(QS.dir.tracking, subdir, 'embeddingTexture', sprintf('trackSeg3d_%06d.png', tp)) ;
        if ~exist(fullfile(QS.dir.tracking, subdir, 'embeddingFigure'), 'dir')
            mkdir(fullfile(QS.dir.tracking, subdir, 'embeddingFigure'))
        end
        if ~exist(fullfile(QS.dir.tracking, subdir, 'embeddingTexture'), 'dir')
            mkdir(fullfile(QS.dir.tracking, subdir, 'embeddingTexture'))
        end
        if exist(seg2dFn, 'file') && ~overwrite
            load(seg2dFn, 'seg2d', 'segIm')
            centroids = seg2d.cdat.centroid ;
            polygons = seg2d.cdat.polygons ;
            cellVtx2D = seg2d.vdat.v ;
            cellIDs = seg2d.cdat.vertexCellIDs ;
        else
            % load segmentation tracks
            segIm = load(fullfile(QS.dir.tracking, subdir, ...
                    sprintf('tracks_label_%06d.mat', tp))) ;
            segIm = segIm.imlabel ;
            nCells = max(segIm(:)) ;

            % check it
            if preview 
                % labeledImage = bwlabel(segIm);
                tmp = label2rgb(segIm+1, 'jet', 'k', 'shuffle') ;
                imshow(tmp)
            end

            % polygons and centroids in 2d 
            props = regionprops(segIm, 'centroid') ;
            seg2d.cdat = struct() ;
            seg2d.cdat.centroid = zeros(max(segIm(:)), 2) ;
            for cid = 1:nCells
                seg2d.cdat.centroid(cid, :) = props(cid).Centroid ;
            end
            polygons = polygonsFromSegmentation(segIm) ;        
            seg2d.cdat.polygons = polygons ;
            centroids = seg2d.cdat.centroid ;

            % Unravel polygons into cellIDs and cellVtx2D
            cellVtx2D = nan(100*length(polygons), 2) ;
            cellIDs = cell(nCells, 1) ;
            dmyk = 1 ;
            for cid = 1:nCells
                poly = polygons{cid} ;
                nAdd = size(poly, 1) ;
                cellVtx2D(dmyk:dmyk+nAdd-1, :) = poly ;
                cellIDs{cid} = dmyk:(dmyk+nAdd-1) ;
                dmyk = dmyk + nAdd ;
            end
            cellVtx2D = cellVtx2D(1:dmyk-1, :) ;
            seg2d.vdat = struct() ;
            seg2d.vdat.v = cellVtx2D ;
            seg2d.cdat.vertexCellIDs = cellIDs ;

            % save processed segmentation2d
            save(seg2dFn, 'seg2d', 'segIm')
        end


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

        % Load/push/tile annular cutmesh
        disp('Loading meshes')
        cutMesh = QS.getCurrentSPCutMeshSmRS() ;
        glueMesh = QS.getCurrentSPCutMeshSmRSC() ;
        cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1) ) ;
        glueMesh.u(:, 1) = glueMesh.u(:, 1) / max(glueMesh.u(:, 1) ) ;
        
        
        % Make sure vertex normals are normalized and outward facing for
        % cutMesh, but not for glueMesh
        normal_shift = metadat.normal_shift * QS.APDV.resolution ;
        if QS.flipy
            % glueMesh vertex normals point IN, so we can shrink it a bit
            cutMesh.vn = - cutMesh.vn ;
        else
            % glueMesh vertex normals point IN, so we can shrink it a bit
            glueMesh.vn = - glueMesh.vn ;
        end
        cutMesh.vn = cutMesh.vn ./ sqrt( sum( cutMesh.vn.^2, 2 ) );
        % Normally evolve vertices
        cutMesh.v = cutMesh.v + normal_shift .* cutMesh.vn;
        glueMesh.v = glueMesh.v + normal_shift .* glueMesh.vn ;
        % Check that faces are outward facing
        trisurf(triangulation(cutMesh.f, cutMesh.v), cutMesh.vn(:, 1), ...
            'edgecolor', 'none')
        
        
        %% PUSH INTO 3D
        seg3dFn = fullfile(QS.dir.tracking, subdir, sprintf('seg3d_%06d.mat', tp)) ;
        if ~exist(seg3dFn, 'file') || overwrite 
            % tile annular cutmesh for triple covering 
            [ TF, TV2D, TV3D, TVN3D ] = tileAnnularCutMesh(cutMesh, [1,1]) ;

            % Allow tracks to disappear for this frame
            disp('Subset of the tracks that are active')
            keepBinary = ~isnan(centroids(:, 1)) ;
            keep = find(keepBinary) ;
            % check that no cells have missing x position but not missing y pos
            assert(all(keepBinary == ~isnan(centroids(:, 2)))) ;

            % Filter
            nCells = max(segIm(:)) ;
            centroids = centroids(keep, :) ;
            cellIDs = cellIDs(keep) ;

            % measurements in 3D of tracked cells
            disp('Perform measurements in 3D')
            cellVtxU = QS.XY2uv(im, cellVtx2D(:, [2, 1]), doubleCovered, 1, 1) ;
            cntrds_uv = QS.XY2uv( im, centroids, doubleCovered, 1, 1) ;

            % embed with measurements
            [c3d, cellCntrd3d, areas, perim, moment1, ang1, ...
                moment2, ang2, moinertia, cellQ2d, cellMeshFaces, vertexMeshFaces, cellNormals] = ...
                polygonNetwork3dMeasurements(TF, TV3D, TV2D, cellVtxU, cellIDs, cntrds_uv) ;

            % Uncertainty in areas is perim * resolution
            disp('Error estimation/propagation')
            TVXY = QS.uv2XY(im, TV2D, doubleCovered, 1, 1) ;
            gg = constructFundamentalForms(TF, TV3D, TVXY) ;
            cntrdSqrtDetg = zeros(length(keep), 1) ;
            perim2d = cntrdSqrtDetg ;
            for cid = 1:length(keep)
                faceID = cellMeshFaces(cid) ;
                cntrdSqrtDetg(cid) = sqrt(det(gg{faceID})) ;
                tmp = regionprops(segIm == keep(cid), 'perimeter') ;
                perim2d(cid) = tmp.Perimeter ;
            end
            unc_areas = cntrdSqrtDetg .* perim2d ;

            save(seg3dFn, 'keep', 'keepBinary', 'cellIDs', ...
                 'centroids', 'cntrds_uv', 'cellCntrd3d', ...
                'c3d', 'cellVtxU', 'areas', 'unc_areas', ...
                'moment1', 'ang1', 'moment2', 'moinertia', 'cellMeshFaces')
        else
            load(seg3dFn, 'keep', 'keepBinary', 'cellIDs', ...
                 'centroids', 'cntrds_uv', 'cellCntrd3d', ...
                'c3d', 'cellVtxU', 'areas', 'unc_areas', ...
                'moment1', 'ang1', 'moment2', 'moinertia', 'cellMeshFaces')
        end
        
        % Plot in 3d
        if ~exist(outImFn, 'file') || overwrite
            close all
            fig = figure('Visible', 'Off', 'units', 'centimeters', ...
                'position', [0,0,xwidth,ywidth]) ;

            % Draw cell faces
            opts = struct() ;
            opts.centroids = cellCntrd3d ;
            opts.vertexBasedTiling = true ;
            [ff, vv, faceMemberIDs] = ...
                polygonsToPatchTriangles3D(c3d, cellIDs, opts) ;
            patch('Faces',ff,'Vertices',vv,...
                'FaceVertexCData',areas(faceMemberIDs),'FaceColor','flat', ...
                'Edgecolor', 'none', 'linewidth', 0.01);
            hold on;

            % Draw contours of cells
            lw = 0.01 ;
            for cid = 1:length(cellIDs)
                poly = c3d(cellIDs{cid}, :) ;
                plot3(poly(:, 1), poly(:, 2), poly(:, 3), 'k-', ...
                    'linewidth', lw)
            end

            axis equal
            xlim(xyzlims(1, :))
            ylim(xyzlims(2, :))
            zlim(xyzlims(3, :))

            colormap viridis
            caxis([0, 100])
            cb = colorbar() ;
            ylabel(cb, ['area, ' QS.spaceUnits '$^{2}$]'], ...
                 'interpreter', 'latex')

            timeinterval = QS.timeInterval ;
            timeunits = QS.timeUnits ;
            t0 = QS.t0set() ;
            titlestr = ['$t = $' num2str(tp*timeinterval-t0) ' ' timeunits] ;
            if blackFigure
                title(titlestr, 'Interpreter', 'Latex', 'Color', 'white') 
            else
                title(titlestr, 'Interpreter', 'Latex', 'Color', 'k') 
            end
            xlabel('AP position [$\mu$m]', 'Interpreter', 'Latex')
            ylabel('lateral position [$\mu$m]', 'Interpreter', 'Latex')
            zlabel('DV position [$\mu$m]', 'Interpreter', 'Latex')

            set(fig, 'PaperUnits', 'centimeters');
            set(fig, 'PaperPosition', [0 0 xwidth ywidth]);
            view(-20, 20) ;

            % Plot mesh for occlusion
            trisurf(triangulation(glueMesh.f, glueMesh.v), 'faceColor', 'k') ;


            % Make background black & Make tick labels white
            if blackFigure
                set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w', 'zcol', 'w')
                set(gcf, 'InvertHardCopy', 'off');
                set(gcf, 'Color', 'k')
                set(gcf, 'color', 'k')
            else
                set(gcf, 'color', 'w')
                set(gca, 'color', 'k')
            end

            % Use export_fig instead, from plotting/export_fig/
            % saveas(fig, fullfile(figvdir, fnv))
            export_fig(outFigFn, '-nocrop', '-r200')

            pim = imread(outFigFn) ;
            size(pim)

            % Texturepatch counterpart  
            opts = metadat ;
            opts.timePoints = timePoints(tidx) ;
            opts.plot_perspective = true ;
            opts.plot_dorsal = false ;
            opts.plot_ventral = false ;
            opts.plot_left = false ;
            opts.plot_right = false ;
            opts.blackFigure = false ;
            opts.makeColorbar = true ;
            QS.plotSeriesOnSurfaceTexturePatch(opts, Options)


            figoutdir = fullfile(QS.dir.texturePatchIm , opts.subdir) ;
            figPerspDir = fullfile(figoutdir, 'perspective') ;
            timFn = fullfile(figPerspDir, sprintf('patch_persp_%06d.png', tp)) ;
            tim = imread(timFn) ;
            size(tim)

            % Add images together -- area to add is in axisROI
            % tim on outside, pim right of colorbarX
            im2 = tim ;
            axisROIFn = fullfile(QS.dir.tracking, subdir, 'embeddingROI_for_overlay.mat') ;
            if ~exist(axisROIFn, 'file')
                axisROI = roipoly(0.5*pim+0.5*tim) ;
                colorbarX = 1010 ;
                save(axisROIFn, 'axisROI', 'colorbarX')    
            else
                load(axisROIFn, 'axisROI', 'colorbarX')
            end
            % sum together inside axisROI
            sim = tim + 0.9 * pim ;
            pR = squeeze(sim(:, :, 1)) ;
            pG = squeeze(sim(:, :, 2)) ;
            pB = squeeze(sim(:, :, 3)) ;
            im2R = squeeze(im2(:,:,1)) ;
            im2G = squeeze(im2(:,:,1)) ;
            im2B = squeeze(im2(:,:,1)) ;
            im2R(axisROI) = pR(axisROI) ;
            im2G(axisROI) = pG(axisROI) ;
            im2B(axisROI) = pB(axisROI) ;
            im2 = cat(3, im2R, im2G, im2B) ;
            % pim right of colorbarX
            im2(:, colorbarX:end, :) = pim(:, colorbarX:end, :) ;

            disp(['Saving combined image to: ' outImFn])
            imwrite(im2, outImFn)

            if preview
                set(gcf, 'visible', 'on')
                imshow(im2)
                pause(3)
            end
        end
    else
        error('handle nuclear tracking case here')
    end
    
end



