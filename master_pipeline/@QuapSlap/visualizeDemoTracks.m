function visualizeDemoTracks(QS, Options) 
% Load demo tracks (segmented or tracked objects) and plot over pullback
% data.
%
% Parameters
% ----------
% QS : 
% Options : optional struct with fields
%
% Returns
% -------
% <none>
%
% Saves to disk
% -------------
% image sequence of RGB overlays
%
% NPMitchell 2021

coordSys = 'spsme' ;
timePoints = QS.xp.fileMeta.timePoints ;
fixedFrame = true ; % true for fixed coordinate systems
t0 = QS.t0set() ;

if isfield(Options, coordSys)
    coordSys = Options.coordSys ;
end

tracks2demo = dir(fullfile(QS.dir.tracking, 'demoTracks', 'demoTracks*.mat')) ;
for qq = 1:length(tracks2demo)
    tname = strsplit(tracks2demo(qq).name, '.mat') ;
    tname = tname{1} ;
    outDir = fullfile(QS.dir.tracking, 'demoTracks', tname) ;
    if ~exist(outDir, 'dir')
        mkdir(outDir)
    end
    tmp = load(fullfile(tracks2demo(qq).folder, tracks2demo(qq).name)) ;
    pos = tmp.positions ;
    timePointsUsed = tmp.timePointsUsed ;
    Ncells = size(pos, 2) ;    
    
    colors = distinguishable_colors(Ncells, [0,0,0;1,1,1]) ;
    
    % Get XY limits for regions
    if fixedFrame
        % Get the image size just once
        im = imread(sprintf(QS.fullFileBase.im_sp_sme, timePointsUsed(1))) ;
        bw0 = false(size(im)) ;
    end
    
    % XYlims for plots
    xmin = Inf ;
    xmax = -Inf ;
    ymin = Inf ;
    ymax = -Inf ;
    for tidx = 1:length(timePointsUsed)

        % if no fixed frame, every timepoint has different size image
        if ~fixedFrame
            im = imread(sprintf(QS.fullFileBase.im_sp_sme, timePointsUsed(1))) ;
        end

        % Get XYlim for this all cells at this timepoint
        for cellID = 1:Ncells
            cellPos = pos{tidx, cellID} ;
            [cx, cy] = ind2sub(size(im), cellPos) ;
            xmin = min(xmin, min(cx)) ;
            xmax = max(xmax, max(cx)) ;
            ymin = min(ymin, min(cy)) ;
            ymax = max(ymax, max(cy)) ;
        end
    end
    xmin = max(1, xmin - buffer) ;
    ymin = max(1, ymin - buffer) ;
    xmax = min(size(im, 1), xmax + buffer) ;
    ymax = min(size(im, 2), ymax + buffer) ;
    
    % Make video for all timepoints
    for tidx = 1:length(timePointsUsed)
        tp = timePointsUsed(tidx) ;
        % Load pullback
        if strcmpi(coordSys, 'spsme')
            im = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
        end
        
        segR = zeros(size(im)) ;
        segG = zeros(size(im)) ;
        segB = zeros(size(im)) ;
        segA = zeros(size(im)) ; % alpha channel
        for cellID = 1:Ncells
            cellPos = pos{tidx, cellID} ;
            if size(cellPos, 2) == 1
                % bw = false(size(im)) ;
                % bw(cellPos) = true ;
                segR(cellPos) = colors(cellID, 1) ;
                segG(cellPos) = colors(cellID, 2) ;
                segB(cellPos) = colors(cellID, 3) ;
                segA(cellPos) = 1 ;
            else
                error("handle XY coordinates here")
                seg(cellPos(:, 1), cellPos(:, 2), 1) = colors(cellID, 1) ;
                seg(cellPos(:, 1), cellPos(:, 2), 2) = colors(cellID, 2) ;
                seg(cellPos(:, 1), cellPos(:, 2), 3) = colors(cellID, 3) ;
            end
        end
        clf
        imRGB = cat(3, im,im,im) ;
        imRGB = imRGB(xmin:xmax, ymin:ymax, :) ;
        seg = cat(3, segR, segG, segB) ;
        seg = seg(xmin:xmax, ymin:ymax, :) ;
        
        % Plot it
        imshow(imRGB) ;
        hold on;
        sh = imshow(seg) ;
        set(sh, 'alphaData', 0.3 * segA(xmin:xmax, ymin:ymax))
        title(['$t = $' num2str(timePointsUsed(tidx) * QS.timeInterval - t0) ' ' QS.timeUnits], ...
            'interpreter', 'latex')
        set(gcf, 'Color', 'w')
        
        % Save as figure
        outfn = fullfile(outDir, sprintf([tname '_%06d.png'], tp)) ;
        export_fig(outfn, '-r150', '-nocrop')
        
        % Save as image
        FF = getframe(gca) ;
        outfn = fullfile(outDir, sprintf(['im_' tname '_%06d.png'], tp)) ;        
        imwrite(FF.cdata, outfn)
    end
end