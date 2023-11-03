% Draw curves of muscle along AP axis on pullback images
clc; close all; 
clearvars -except tubi

addpath('/mnt/data/code/gut_matlab/travelling_salesman/tsp_ga')

outdir = fullfile(tubi.dir.tracking, 'muscleStripes');
trackfn = fullfile(tubi.dir.tracking, 'muscle_tracks.mat') ;
tmp = load(trackfn) ;
tracks = tmp.tracks ;
imFileBase = tubi.fullFileBase.im_sp_sm ;

ncurves = 4 ;
timePoints = tubi.xp.fileMeta.timePoints ;

for curv = 1:ncurves 
        
    % First grab all tracks available for this curv by looking for any
    % point in a chosen polygon enclosing the curve at t0
    tidx0 = tubi.xp.tIdx(tubi.t0) ;
    % Gather all tracks at tidx0
    ids = [] ;
    for track = 1:length(tracks)
        ids = cat(1, ids, tracks{track}(tidx0, :)) ;
    end

    polyfn = fullfile(outdir, sprintf('polygon%d_t0.mat', curv)) ;
    if ~exist(polyfn, 'file')
        % Gather all tracks for this stripe
        disp('Gathering all tracks for this stripe at t=t0')
    
        % Load image
        close all
        im = imread(sprintf(imFileBase, tubi.t0)) ;
        imshow(im) ;
        hold on;
        plot(ids(:, 1), ids(:, 2), 'o')
    
        % Select tracks that are part of this curve
        h = drawpolygon;
        
        % Wait for the user to finish drawing the polygon
        wait(h);
        
        % Get the polygon's vertices
        polygonVertices = h.Position;
    
        % Save the polygon
        save(polyfn, 'polygonVertices') 
    else
        load(polyfn, 'polygonVertices')
    end
    
    % Check which points are inside the drawn polygon
    inPolygon = inpolygon(ids(:, 1), ids(:, 2), ...
        polygonVertices(:, 1), polygonVertices(:, 2));
    
    % Select the xy points that are inside the polygon
    sids = ids(inPolygon, :);
    
    % Display the selected points on the image
    clf
    im = imread(sprintf(imFileBase, tubi.t0)) ;
    imshow(im) ;
    hold on;
    plot(ids(:, 1), ids(:, 2), 'o')
    hold on;
    plot(sids(:, 1), sids(:, 2), 'ro', 'MarkerSize', 10);
    title('Tracks that are part of this stripe')
    waitfor(gcf)

    % Load the curve positions that are already defined on disk
    curvfn = fullfile(outdir, 'curvePositions.mat') ;
    if exist(curvfn, 'file')
        load(curvfn, 'curvePositions')
    else
        curvePositions = cell(length(timePoints), 1) ;
    end

    % Now add to these tracks to fill in the stripe
    allPts = {} ;
    for tidx = 1:length(timePoints)
        tp = timePoints(tidx) ;

        % Get known tracks for this timepoint
        disp(['Gathering tracks for this stripe at t=' num2str(tp)])
        ids = [] ;
        for track = 1:length(tracks)
            ids = cat(1, ids, tracks{track}(tidx, :)) ;
        end
        % Keep only selected tracks for this curve
        ids = ids(inPolygon, 1:2) ; 

        % Add more points as needed to draw curves along rostral/caudal axis
        clf
        im = imread(sprintf(imFileBase, tp)) ;
        imshow(im) ;
        hold on;
        plot(ids(:, 1), ids(:, 2), 'o')

        % Add to the current ids
        if isempty(curvePositions{tidx})
            h = drawfreehand;
            % Wait for the user to modify the curve
            wait(h);
            % Store the modified coordinates for the current image
            curvePositions{tidx} = h.Position;
            % Save updated curvePositions
            disp('saving updates...')
            save(curvfn, 'curvePositions')
        else
            % h = drawpolygon(gca, 'Position', curvePositions{tidx});
            % title(['t=' num2str(tp)])
            % pause(0.1)
        end
        


        % Adding to all points
        % ids = [ids; curvePositions{tidx}] ;
        % allPts{tidx} = curvePositions{tidx} ;

        % Connect these points in order using traveling salesman
        % userConfig = struct('xy', ids);
        % res = tsp_ga(userConfig);

        % error('here')
    end

end


