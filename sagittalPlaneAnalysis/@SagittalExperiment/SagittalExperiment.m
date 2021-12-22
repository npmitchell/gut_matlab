classdef SagittalExperiment < handle
   
    properties
        currentTime = 0 
        tidx = 1 
        timePoints = []         % timestamps of the experiment (in file names)
        dt = 2 ;                % increment in time between timepoints with 
                                % indices differing by 1. For example, if
                                % timePoints are [0,1,2,4] and these are 
                                % [1,1,2] minutes apart, then timeInterval 
                                % is 1. 
        resolution              % SE.spaceUnits per pixel
        name = 'sagittal_Z%d_T%03d' ;
        imExten = '.png' ;
        labelExten = '.tiff' ;
        npts_skel = 10000 ;
        minsz_hole = 50 ;
        zplanes = []            % 
        channels = [1,2];       % Channels
        norms = [120, 250] ;    % intensity normalization for channel 1/2
        timeUnits = 'min'       % units of the timeInterval (ex 'min')
        spaceUnits = '$\mu$m'   % units of the embedding space (ex '$\mu$m')
        dir = struct()          % struct of directory paths where data lives -- to be sprintf(..., zplane)
        filename = struct('landmarks', '', ...
            'lobeFilters', '')     % struct of filename paths where data lives -- to be sprintf(..., zplane)
        zDir = struct('landmarks', '', ...
            'lobeFilters', '')         % str, directory where data lives -- evaluated at z planes
        t0                      % reference time in the experiment
        bw = cell(0)
        landmarks = cell(0) 
        folds = struct()
        pathMasks = cell(0) 
        thicknessMasks = cell(0) 
        skelIntensity = cell(0) 
        zplane
        skel = struct() ;
        thickness = cell(0) 
        plotting = struct('preview', false, ... % display intermediate results
            'colors', [])                       % color cycle for SE
        im = []                 % N-channel TIFF image
        imL = [] ;              % labeled image (RGB)
    end
    
    % Some methods are hidden from public view. These are used internally
    % to the class.
    methods (Hidden)
        function im = loadImage(~, fname, exten)
            % add default file extension if not present
            if ~contains(fname, '.')
                fname = [fname exten] ;
            end
            disp(['Loading ', fname])
            if strcmpi(fname(end-3:end), 'tif' ) || strcmpi(fname(end-4:end), 'tiff') 
                im = loadtiff(fname) ;
            else
                im = imread(fname) ;
            end
        end
    end
    
    % Public methods, accessible from outside the class and reliant on 
    % properties of the class instance
    methods
        
        function SE = SagittalExperiment(settings)
            SE.dir.data = settings.dataDir ;
            SE.channels = settings.channels ;
            SE.dt = settings.dt ;
            SE.spaceUnits = settings.spaceUnits ;
            SE.timeUnits = settings.timeUnits ;
            SE.timePoints = settings.timePoints ;
            SE.name = settings.name ;
            SE.norms = settings.norms ;
            SE.npts_skel = settings.npts_skel ;
            SE.minsz_hole = settings.minsz_hole ;
            SE.zplanes = settings.zplanes ;
            SE.resolution = settings.resolution ;
            SE.plotting.colors = define_colors ;
            if isfield(settings, 'imExten')
                SE.imExten = settings.imExten ;
            else
                SE.imExten = '.tif' ;
            end
            if isfield(settings, 'labelExten')
                SE.labelExten = settings.labelExten ; 
            else
                SE.labelExten = SE.imExten ; 
            end
            
            % Create directories
            SE.dir.zIm = fullfile(SE.dir.data, 'z%d') ; 
            zDir = SE.dir.zIm ;
            SE.dir.zLabeled = [zDir '_labeled' ] ; 
            SE.dir.zOut = [zDir '_results'] ;
            zOutDir = SE.dir.zOut ;
            SE.dir.skel = fullfile(zOutDir, 'skels') ;
            SE.dir.bw = fullfile(zOutDir, 'bw') ;
            SE.dir.landmarks = fullfile(zOutDir, 'landmarks') ;
            SE.dir.pathMasks = fullfile(zOutDir, 'pathmasks') ;
            SE.dir.thicknessMasks = fullfile(zOutDir, 'thicknessmasks') ;
            SE.dir.thickness = fullfile(zOutDir, 'thickness') ;
            SE.dir.thickness_images = fullfile(zOutDir, 'thickness_images') ;
            SE.dir.skelIntensity = fullfile(zOutDir, 'skelIntensity') ;
            SE.dir.curvature = fullfile(zOutDir, 'curvature') ;
            SE.dir.kymographs = fullfile(zOutDir, 'kymographs') ;
            SE.dir.thickness_kymographs =  ...
                fullfile(zOutDir, 'thickness_kymographs') ;
            SE.dir.lobeFilters = fullfile(zOutDir, 'lobeFilters') ;
            
            % Filenames
            SE.filename.kymo = fullfile(SE.dir.kymographs, 'kymoData.mat') ;
            SE.filename.tKymo = fullfile(SE.dir.thickness_kymographs, 'kymoThickness.mat') ;
            SE.filename.curvature = fullfile(SE.dir.curvature, 'curvature.mat') ;
            SE.filename.landmarks = fullfile(SE.dir.landmarks, ...
                    sprintf('landmarks_%04d.mat', SE.currentTime)) ;
                
            for zz = SE.zplanes
                % Make dirs
                dirs = {sprintf(SE.dir.zIm, zz), ...
                    sprintf(SE.dir.zLabeled, zz), ...
                    sprintf(SE.dir.zOut, zz), ...
                    sprintf(SE.dir.skel, zz), ...
                    sprintf(SE.dir.bw, zz), ...
                    sprintf(SE.dir.landmarks, zz), ...
                    sprintf(SE.dir.pathMasks, zz), ...
                    sprintf(SE.dir.thicknessMasks, zz), ...
                    sprintf(SE.dir.thickness, zz), ...
                    sprintf(SE.dir.thickness_images, zz), ...
                    sprintf(SE.dir.skelIntensity, zz), ...
                    sprintf(SE.dir.lobeFilters, zz)} ;
                for dd = 1:length(dirs)
                    if ~exist(dirs{dd}, 'dir')
                        mkdir(dirs{dd})
                    end
                end
            end
        end
        
        function setTime(SE, tt)
            % Set the current time of the dataset and clear current data
            % which was associated with the previously considered time
            %
            % Parameters
            % ----------
            % tt : int or float
            %   timePoint to set to be current, from available times in
            %   SE.xp.fileMeta.timePoints
            %
            if tt ~= SE.currentTime
                SE.reset() ;
            end
            SE.currentTime = tt ;
            SE.tidx = find(SE.timePoints == tt) ;
            
            % set specific filenames
            SE.filename.landmarks = fullfile(SE.zDir.landmarks, ...
                    sprintf('landmarks_%04d.mat', SE.currentTime)) ;
            SE.filename.lobeFilters = fullfile(SE.zDir.lobeFilters, ...
                    sprintf('lobeFilters_%04d.mat', SE.currentTime)) ;
        end
        
        function setZPlane(SE, zz)
            % Set the current time of the dataset and clear current data
            % which was associated with the previously considered time
            %
            % Parameters
            % ----------
            % zz : int 
            %   zplane to set to be current, from available times in
            %   SE.zplanes
            %
            if zz ~= SE.zplane
                SE.reset() ;
            end
            SE.zplane = zz ;
            SE.zDir.zIm = sprintf(SE.dir.zIm, zz) ;
            SE.zDir.zLabeled = sprintf(SE.dir.zLabeled, zz) ;
            SE.zDir.zOut = sprintf(SE.dir.zOut, zz) ;
            SE.zDir.skel = sprintf(SE.dir.skel, zz) ;
            SE.zDir.bw = sprintf(SE.dir.bw, zz) ;
            SE.zDir.landmarks = sprintf(SE.dir.landmarks, zz) ;
            SE.zDir.pathMasks = sprintf(SE.dir.pathMasks, zz) ;
            SE.zDir.thicknessMasks = sprintf(SE.dir.thicknessMasks, zz) ;
            SE.zDir.thickness = sprintf(SE.dir.thickness, zz) ;
            SE.zDir.thickness_images = sprintf(SE.dir.thickness_images, zz);
            SE.zDir.skelIntensity = sprintf(SE.dir.skelIntensity, zz) ;
            SE.zDir.curvature = sprintf(SE.dir.curvature, zz) ;
            SE.zDir.kymographs = sprintf(SE.dir.kymographs, zz) ;
            SE.zDir.thickness_kymographs = ...
                sprintf(SE.dir.thickness_kymographs, zz) ;
            SE.zDir.lobeFilters = sprintf(SE.dir.lobeFilters, zz) ;
            
            % Set specific filenames
            SE.filename.kymo = sprintf(SE.filename.kymo, zz) ;
            SE.filename.tKymo = sprintf(SE.filename.tKymo, zz) ;
            SE.filename.curvature = sprintf(SE.filename.tKymo, zz) ;                
        end
        
        function reset(SE)
            % clear current timepoint's data for SE instance
            SE.landmarks = cell(0) ;
            SE.pathMasks = cell(0) ;
            SE.skel = struct() ;
            SE.skelIntensity = struct() ;
            SE.thicknessMasks = cell(0) ;
            SE.thickness = cell(0) ;
            SE.bw = [] ;
            SE.im = [] ;
            SE.imL = [] ;
        end
        
        function t0 = t0set(SE, t0)
            % t0set(SE, t0) Set time offset to 1st fold onset 
            SE.t0 = t0 ;
        end
        
        function im = getIm(SE)
            if ~isempty(SE.im)
                im = SE.im ;
            else
                % Original image
                if sum(SE.name == '%') == 2
                    fname = fullfile(SE.zDir.zIm, sprintf(SE.name, SE.zplane, SE.currentTime)) ;
                elseif sum(SE.name == '%') == 1
                    fname = fullfile(SE.zDir.zIm, sprintf(SE.name, SE.currentTime)) ;
                else
                    error('Filename SE.name seems to have more than 2 numeric fields.')
                end
                im = SE.loadImage(fname, SE.imExten) ;
                if size(im, 3) == 2
                    imn = zeros(size(im, 1), size(im, 2), 3) ;
                    imn(:, :, 1) = double(im(:, :, 1)) / SE.norms(1) ;
                    imn(:, :, 2) = double(im(:, :, 2)) / SE.norms(2) ;
                    im = imn ;
                    clearvars im2 
                end
                SE.im = im ; 
            end
        end
        
        function imL = getImL(SE)
            if ~isempty(SE.imL)
                imL = SE.imL ;
            else
                % Labeled image
                if sum(SE.name == '%') == 2
                    fname = fullfile(SE.zDir.zLabeled, sprintf(SE.name, SE.zplane, SE.currentTime)) ;
                elseif sum(SE.name == '%') == 1
                    fname = fullfile(SE.zDir.zLabeled, sprintf(SE.name, SE.currentTime)) ;
                else
                    error('Filename SE.name seems to have more than 2 numeric fields.')
                end
                imL = SE.loadImage(fname, SE.labelExten) ;
                SE.imL = imL ;
            end
        end
        
        function bw = getBW(SE, overwrite)
            if nargin < 2
                overwrite = false ;
            end
            if ~isempty(SE.bw)
                bw = SE.bw ;
            else
                % Segment out the membrane/gut cross section
                bwFn = fullfile(SE.zDir.bw, sprintf('bw_%04d.mat', SE.currentTime)) ;
                bwFnPNG = fullfile(SE.zDir.bw, sprintf('bw_%04d.png', SE.currentTime)) ;
                if exist(bwFn, 'file') && ~overwrite
                    load(bwFn, 'bw')
                    if ~exist(bwFnPNG, 'file')                     
                        imwrite(bw, bwFnPNG)
                    end
                else
                    % Probabilities (ch1 is midgut)
                    if sum(SE.name == '%') == 1
                        pname = fullfile(SE.zDir.zIm, [sprintf(SE.name, SE.currentTime), '_Probabilities.h5']) ;    
                    elseif sum(SE.name == '%') == 2
                        pname = fullfile(SE.zDir.zIm, [sprintf(SE.name, SE.zplane, SE.currentTime), '_Probabilities.h5']) ;
                    else
                        error('SE.name seems to have more than 2 numeric fields.')
                    end
                    prob = h5read(pname, '/exported_data') ;

                    % Segment the probabilities
                    bw = squeeze(prob(1, :, :)) > 0.5 ;
                    % get largest region size and display it
                    regions = bwconncomp(bw) ;
                    regs = regionprops(regions) ;
                    % Get region sizes
                    regsizes = zeros(length(regs), 1) ;
                    for dmyq = 1:length(regs)
                        regsizes(dmyq) = regs(dmyq).Area ;
                    end
                    [maxszk, regID] = maxk(regsizes, 3) ;
                    canvas = false(size(bw)) ;
                    canvas(cell2mat(regions.PixelIdxList(regID(1)))) = true;            
                    % Check if more regions needed
                    clf
                    
                    % Invert, filter, then invert to close holes
                    canvas2 = ~canvas ;
                    canvas2 = bwareafilt(canvas2, [SE.minsz_hole, Inf]) ;
                    canvas2 = ~canvas2 ;

                    imshow(canvas2')
                    title(['segmentation, t = ' num2str(SE.currentTime)]) 
                    doneRegs = false ;
                    nextId = 2 ;
                    msg = 'More regions needed?[y/n, default=no]' ;
                    title(msg)
                    moreRegs = input(msg, 's') ;
                    while ~doneRegs
                        if ~isempty(moreRegs)
                            if strcmpi(moreRegs, 'y')
                                canvas(cell2mat(regions.PixelIdxList(regID(nextId)))) = true;
                                nextId = min(nextId + 1, length(regID)) ;
                                
                                % Invert, filter, then invert to close holes
                                canvas2 = ~canvas ;
                                canvas2 = bwareafilt(canvas2, [SE.minsz_hole, Inf]) ;
                                canvas2 = ~canvas2 ;
                                imshow(canvas2')

                                moreRegs = input('More regions needed?[y/n, default=no]', 's') ;
                            else
                                doneRegs = true ;
                            end
                        else
                            doneRegs = true ;
                        end
                    end
                    bw = bwareafilt(canvas, [maxszk(nextId-1)-1, maxszk(1)+1]) ;

                    % Invert, filter, then invert to close holes
                    wb = ~bw ;
                    wb = bwareafilt(wb, [SE.minsz_hole, Inf]) ;
                    bw = ~wb ;

                    % Save segmentation
                    save(bwFn, 'bw')
                    imwrite(bw, bwFnPNG)
                end
                
                % Attribute bw
                SE.bw = bw ;
            end
        end
        
        function landmarks = getLandmarks(SE, overwrite)
            if nargin < 2
                overwrite = false ;
            end
            if ~isempty(SE.landmarks)
                landmarks = SE.landmarks ;
            else
                % Load or compute landmarks
                lmkfn = fullfile(SE.zDir.landmarks, ...
                    sprintf('landmarks_%04d.mat', SE.currentTime)) ;
                if exist(lmkfn, 'file') && ~overwrite
                    load(lmkfn, 'landmarks')
                else
                    % Load previous landmarks for reference
                    if SE.tidx > 1
                        options = struct() ;
                        prevfn = fullfile(SE.zDir.landmarks, ...
                            sprintf('landmarks_%04d.mat', SE.timePoints(SE.tidx-1))) ;
                        tmp = load(prevfn, 'landmarks') ;
                        options.reference_landmarks = tmp.landmarks ; 
                    else
                        options = struct() ;
                    end
                    options.bgim = SE.imL ;

                    options.title = ['Identify landmarks, t=' num2str(SE.currentTime)] ;
                    landmarks = queryLandmarks(SE.bw', options) ; 
                    % Save landmarks
                    save(lmkfn, 'landmarks')
                end

                % Attribute landmarks
                SE.landmarks = landmarks ;
            end
        end
        
        function pathMasks = getPathMasks(SE, overwrite, forceSelection, prev_or_fix)
            if nargin < 2
                overwrite = false ;
            end
            if nargin < 3
                forceSelection = false ;
            end
            if nargin < 4
                prev_or_fix = 'prev' ;
            end
            if ~isempty(SE.pathMasks)
                pathMasks = SE.pathMasks ;
            else
                % Path masks for curve extraction
                pmfn = fullfile(SE.zDir.pathMasks, ...
                    sprintf('pathMasks_%04d.mat', SE.currentTime)) ;
                if exist(pmfn, 'file') && ~overwrite
                    load(pmfn, 'pathMasks')
                else
                    % Load previous pathMasks for reference
                    options = struct() ;
                    if SE.tidx > 1
                        if strcmpi(prev_or_fix, 'prev')
                            prevfn = fullfile(SE.zDir.pathMasks, ...
                                sprintf('pathMasks_%04d.mat', SE.timePoints(SE.tidx-1))) ;
                            tmp = load(prevfn, 'pathMasks') ;
                            options.reference_pathMasks = tmp.pathMasks ; 
                        elseif strcmpi(prev_or_fix, 'fix')
                            prevfn = fullfile(SE.zDir.pathMasks, ...
                                sprintf('pathMasks_%04d.mat', SE.timePoints(SE.tidx))) ;
                            tmp = load(prevfn, 'pathMasks') ;
                            options.reference_pathMasks = tmp.pathMasks ; 
                        end
                    end
                    SE.getBW() ;
                    SE.getImL() ;
                    SE.getLandmarks() ;
                    options.bgim = SE.imL ;
                    options.default_required = true ;
                    options.forceSelection = forceSelection ;

                    % Query masks for computing geodesics
                    options.title = ['Identify path mask, t=' num2str(SE.currentTime)] ;
                    pathMasks = queryPathMasks(SE.bw', SE.landmarks, options) ;
                    % Save path masks
                    save(pmfn, 'pathMasks')
                end
                
                % Attribute pathMasks
                SE.pathMasks = pathMasks ;
            end
        end
        
        function thicknessMasks = getThicknessMasks(SE, overwrite, path_or_prev)
            % Parameters
            % ----------
            % overwrite : bool (default=false)
            %   overwrite previous results
            % path_or_prev : str ('path', 'prev', or 'fix')
            %
            if nargin < 2
                overwrite = false ;
            end
            if nargin < 3
                path_or_prev = 'prev' ;
            end
            if ~isempty(SE.thicknessMasks)
                thicknessMasks = SE.thicknessMasks ;
            else
                % Thickness masks for thickness extraction via DT sampling
                tmfn = fullfile(SE.zDir.thicknessMasks, ...
                    sprintf('thicknessMasks_%04d.mat', SE.currentTime)) ;
                tmFigFn = fullfile(SE.zDir.thicknessMasks, ...
                    sprintf('thicknessMasks_%04d.png', SE.currentTime)) ;
                if exist(tmfn, 'file') && ~overwrite
                    load(tmfn, 'thicknessMasks')
                else
                    % Load previous pathMasks for reference
                    SE.getImL() ;
                    options = struct() ;
                    options.bgim = SE.imL ;
                    if SE.tidx > 1
                        if strcmpi(path_or_prev, 'prev')
                            prevfn = fullfile(SE.zDir.thicknessMasks, ...
                                sprintf('thicknessMasks_%04d.mat', SE.timePoints(SE.tidx-1))) ;
                            tmp = load(prevfn, 'thicknessMasks') ;
                            options.reference_pathMasks = tmp.thicknessMasks ; 
                        elseif strcmpi(path_or_prev, 'path')
                            prevfn = fullfile(SE.zDir.pathMasks, ...
                                sprintf('pathMasks_%04d.mat', SE.currentTime)) ;
                            tmp = load(prevfn, 'pathMasks') ;
                            options.reference_pathMasks = tmp.pathMasks ; 
                        elseif strcmpi(path_or_prev, 'fix')
                            prevfn = fullfile(SE.zDir.thicknessMasks, ...
                                sprintf('thicknessMasks_%04d.mat', SE.currentTime)) ;
                            tmp = load(prevfn, 'thicknessMasks') ;
                            options.reference_pathMasks = tmp.thicknessMasks ; 
                        end
                    end
                    options.ref_priority = true ;
                    options.bgim = SE.imL ;
                    options.default_required = true ;
                    
                    % Query masks for computing thickness
                    DTasPath = input('DT masks for thickness same as pathMasks? [Y/n]', 's') ;
                    
                    if isempty(DTasPath) || contains(lower(DTasPath), 'y')  
                        thicknessMasks = SE.getPathMasks() ;
                        % Save thickness masks
                        save(tmfn, 'thicknessMasks')
                    else
                        SE.getLandmarks() ;
                        SE.getIm() ;
                        SE.getBW() ;
                        options.title = ['Identify DT mask, t=' num2str(SE.currentTime)] ;
                        options.promptString = ['Are DT masks required? t=' num2str(SE.currentTime) ' [Y/n]'] ;
                        thicknessMasks = queryPathMasks(SE.bw', SE.landmarks, options) ;
                        % Save thickness masks
                        save(tmfn, 'thicknessMasks')
                    end

                    % Save thickness mask image
                    close all
                    figure('visible', 'off')
                    imshow(max(SE.imL(:))*0.2 + 0.8 * SE.imL); hold on;
                    for qq = 1:length(thicknessMasks)
                        tmqq = thicknessMasks{qq} ;
                        colors = define_colors(length(tmqq.mask)) ;
                        for pp = 1:length(tmqq.mask)
                            patch(tmqq.xroi{pp}, tmqq.yroi{pp}, colors(pp, :), ...
                                'FaceAlpha', 0.2, 'EdgeColor', 'none') ;
                        end
                    end
                    title(['thickness DT mask: t = ' num2str(SE.currentTime)])
                    saveas(gcf, tmFigFn)
                    close all
                end
                % Attribute thicknessMasks
                SE.thicknessMasks = thicknessMasks ;
            end
        end
        
        skel = getSkel(SE, overwrite)        
        thickness = getThickness(SE, overwrite)
        skelIntensity = getSkelIntensity(SE, overwrite)
        kymoData = getKymographs(SE, overwriteData, overwriteImages) 
        tKymoData = getThicknessKymographs(SE, overwriteData, overwriteImages) 
        curvature = getSkelCurvature(SE, overwriteData, overwriteImages)
        
        function identifyFolds(SE)
            % Identify fold locations for all timepoints
            SE.folds.v1 = zeros(length(SE.timePoints), 1) ;
            SE.folds.v2 = zeros(length(SE.timePoints), 1) ;
            SE.folds.v3 = zeros(length(SE.timePoints), 1) ;
            SE.folds.d1 = zeros(length(SE.timePoints), 1) ;
            SE.folds.d2 = zeros(length(SE.timePoints), 1) ;
            SE.folds.d3 = zeros(length(SE.timePoints), 1) ;

            for tidx = 1:length(SE.timePoints)
                % Plot the image
                imL = SE.getImL() ;
                imshow(imL); hold on;
                tidx0 = SE.tidx ;
                colors = parula(length(SE.timePoints)+1) ;
                for tt = tidx0:length(SE.timePoints)
                    SE.setTime(SE.timePoints(tt))
                    lm = SE.getLandmarks() ;
                    for qq = 1:length(lm)
                        ids = strsplit(lm{qq}.id, '/') ;
                        for pp = 1:length(ids)-1
                            if ~contains(ids{pp}, 'u')
                                scatter(lm{qq}.v(pp, 1), lm{qq}.v(pp, 2), 40, ...
                                    'markeredgecolor', colors(tt-tidx0+1, :), ...
                                    'markeredgeAlpha', 0.5) ;
                            end
                        end
                    end 
                end
                SE.setTime(SE.timePoints(tidx0))
                % Overlay skeleton
                thick = SE.getThickness() ;
                skel = SE.getSkel() ;
                skels = skel.skel_ss ;
                tum = thick.avg_thickness_um ;
                for qq = 1:length(skels)
                    scatter(skels{qq}(:, 1), skels{qq}(:, 2), 2, tum{qq})
                end
                cb = colorbar() ;
                ylabel(cb, ['thickness [' SE.spaceUnits ']'], 'interpreter', 'latex')
            
                % 
                
            end
            
            % Save results
        end
        
        function lobeFilters = getLobeFilters(SE, overwrite)
            if nargin < 2
                overwrite = false ;
            end
            
            lobeMaskFn = SE.filename.lobeFilters ;
            if exist(lobeMaskFn, 'file') && ~overwrite
                lobeFilters  = load(lobeMaskFn) ;
            else
                % Load previous masks, if not first timepoint, or load
                % existing masks on disk if first timepoint
                if SE.tidx > 1
                    prevMaskFn = fullfile(SE.zDir.lobeFilters, ...
                            sprintf('lobeFilters_%04d.mat', SE.timePoints(SE.tidx-1))) ;
                else
                    try
                        prevMaskFn = fullfile(SE.zDir.lobeFilters, ...
                            sprintf('lobeFilters_%04d.mat', SE.timePoints(SE.tidx))) ;
                    catch
                        disp('No previous lobeMask file found')
                        prevMaskFn = '' ;
                    end
                end

                % Plot the image
                imL = SE.getImL() ;
                imshow(imL); hold on;
                tidx0 = SE.tidx ;
                colors = parula(length(SE.timePoints)+1) ;
                for tt = tidx0:length(SE.timePoints)
                    SE.setTime(SE.timePoints(tt))
                    lm = SE.getLandmarks() ;
                    for qq = 1:length(lm)
                        ids = strsplit(lm{qq}.id, '/') ;
                        for pp = 1:length(ids)-1
                            if ~contains(ids{pp}, 'u')
                                scatter(lm{qq}.v(pp, 1), lm{qq}.v(pp, 2), 40, ...
                                    'markeredgecolor', colors(tt-tidx0+1, :), ...
                                    'markeredgeAlpha', 0.5) ;
                            end
                        end
                    end 
                end
                SE.setTime(SE.timePoints(tidx0))
                % Overlay skeleton
                thick = SE.getThickness() ;
                skel = SE.getSkel() ;
                skels = skel.skel_ss ;
                tum = thick.avg_thickness_um ;
                for qq = 1:length(skels)
                    scatter(skels{qq}(:, 1), skels{qq}(:, 2), 2, tum{qq})
                end
                cb = colorbar() ;
                ylabel(cb, ['thickness [' SE.spaceUnits ']'], 'interpreter', 'latex')

                % ID Lobes for this timepoint
                masks = cell(8, 1);
                polygons = cell(8, 1);
                thickness = cell(8, 1);
                thickness_avg = zeros(8, 1);
                thickness_std = zeros(8, 1);
                % try to load previous polygon
                prevMask = cell(8, 1) ;
                prevPoly = cell(8, 1) ;
                if exist(prevMaskFn, 'file')
                    prev = load(prevMaskFn) ;
                    prevMask = prev.masks ;
                    prevPoly = prev.polygons ;
                end

                for qq = 1:length(masks)
                    % show previous as overlay
                    if ~isempty(prevPoly{qq})
                        plot(prevPoly{qq}(:, 1), prevPoly{qq}(:, 2), 'w-') ;
                    end
                    msg = ['Select Lobe ' num2str(qq) ' mask'] ;
                    disp(msg)
                    title(msg)
                    [lobem, lobex, lobey] = roipoly() ;
                    lobePoly = [lobex, lobey] ;
                    if isempty(lobem)
                        % no polygon given, use prev if available
                        useprev = input('Use previous? [Y/n]', 's') ;
                        if ~contains(lower(useprev), 'n')
                            % Use previous
                            masks{qq} = prevMask{qq} ;
                            polygons{qq} = prevPoly{qq} ;                         
                        end
                    else
                        masks{qq} = lobem ;
                        polygons{qq} = lobePoly ;
                    end
                    if isempty(polygons{qq})
                        disp('no thickness data for this lobe')
                        thickness{qq} = [] ;
                        thickness_avg(qq) = NaN ;
                        thickness_std(qq) = NaN ;   
                    else
                        % Find all unique-pixel thickness measurements inside 
                        % the polygon
                        thickness{qq} = [] ;
                        for ii = 1:length(skels)
                            [~, idx] = uniqueRows(round(skels{ii})) ;
                            inside = find(inpolygon(skels{ii}(idx, 1), ...
                                skels{ii}(idx, 2), ...
                                polygons{qq}(:, 1), polygons{qq}(:, 2))) ;
                            ths = thick.avg_thickness_um{ii}(idx) ;
                            thickness{qq} = [thickness{qq} ths(inside)] ;
                        end
                        thickness_avg(qq) = mean(thickness{qq})
                        thickness_std(qq) = std(thickness{qq})
                    end
                end
                % Save lobeFilters and lobePolygons
                disp(['Saving lobeFilters to ' lobeMaskFn])
                save(lobeMaskFn, 'masks', 'polygons', ...
                    'thickness', 'thickness_avg', 'thickness_std')
                close all
                
                if nargout > 0
                    lobeFilters = struct('masks', masks, ...
                        'polygons', polygons, ...
                        'thickness', thickness, ...
                        'thickness_avg', thickness_avg, ...
                        'thickness_std', thickness_std) ;
                end
            end
        end
        
        function getLobeThicknesses(SE, overwrite)
            if nargin < 2
                overwrite = false ;
            end
            tavgs = zeros(length(SE.timePoints), 8) ;
            tstds = zeros(length(SE.timePoints), 8) ;
            for tidx = 1:length(SE.timePoints) 
                tp = SE.timePoints(tidx) ;
                SE.setTime(tp) ; 
                lobe = getLobeFilters(SE, overwrite) ;
                tavgs(tidx, :) = lobe.thickness_avg ;
                tstds(tidx, :) = lobe.thickness_std ;
            end
            
            % Plot them
            for norm = 0:1
                close all
                % First side (ventral or right)
                % colors = get(gca,'colororder') ;            
                colors = SE.plotting.colors ; 
                subplot(2, 1, 1)
                for qq = 1:4
                    lineprops = {'-', 'markerfacecolor', colors(qq, :), ...
                        'lineWidth', 2, 'color', colors(qq, :)} ;
                    if norm
                        hh = shadedErrorBar(SE.timePoints, tavgs(:, qq)/tavgs(1, qq), ...
                            tstds(:, qq)/tavgs(1, qq), 'lineprops', lineprops, ...
                            'transparent',1, 'patchSaturation', 0.1) ;
                    else
                        hh = shadedErrorBar(SE.timePoints, tavgs(:, qq), ...
                            tstds(:, qq), 'lineprops', lineprops, ...
                            'transparent',1, 'patchSaturation', 0.1) ;
                    end
                    hold on;
                    hhs(qq) = hh.patch ;
                end
                legend(hhs, ...
                    {'lobe 1', 'lobe 2', 'lobe 3', 'lobe 4'})
                if norm
                    ylabel('thickness $t/t_0$', 'interpreter', 'latex')
                else
                    ylabel(['thickness [' SE.spaceUnits ']'], 'interpreter', 'latex')
                end
                title('Thickness of lobes over time', 'interpreter', 'latex')

                
                % Other side (dorsal or left)
                subplot(2, 1, 2)
                for qq = 1:4
                    q2 = 8-qq+1 ;
                    lineprops = {'-', 'markerfacecolor', colors(qq, :), ...
                        'lineWidth', 2, 'color', colors(qq, :)} ;
                    if norm
                        shadedErrorBar(SE.timePoints, tavgs(:, q2)/tavgs(1, q2), ...
                            tstds(:, q2)/tavgs(1, q2), 'lineprops', lineprops, ...
                            'transparent',1, 'patchSaturation', 0.1)
                    else
                        shadedErrorBar(SE.timePoints, tavgs(:, q2), ...
                            tstds(:, q2), 'lineprops', lineprops, ...
                            'transparent',1, 'patchSaturation', 0.1)                        
                    end
                    hold on;
                end
                xlabel(['time [' SE.timeUnits ']'], 'interpreter', 'latex')
                if norm
                    ylabel('thickness $t/t_0$', 'interpreter', 'latex')
                    saveas(gcf, fullfile(SE.zDir.lobeFilters, 'lobeThickness_norm.png'))
                else
                    ylabel(['thickness [' SE.spaceUnits ']'], 'interpreter', 'latex')
                    saveas(gcf, fullfile(SE.zDir.lobeFilters, 'lobeThickness.png'))
                end
            end
        end
    end
    
    methods (Static)       
    end
    
end