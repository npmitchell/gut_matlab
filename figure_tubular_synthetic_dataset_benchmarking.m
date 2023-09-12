% Plot benchmarking for synthetic dataset -- positions, velocities, and
% tracking via pullbacks against groundTruth

close all 
clear
clc

overwrite = false ;

%% Add paths
dataDir = '/mnt/data/tubular_test/single_coil_dynamic_fromScratch2/' ;  % This is the path where the example dataset is downloaded
origpath = '/mnt/data/code/tubular/' ; % This should be the path to tubular, like /mnt/data/code/tubular/
cd(origpath)
addpath(origpath) ;
addpath(fullfile('utility', 'addpath_recurse'))
addpath_recurse('utility')
% Here we do use gptoolbox code: add to the MATLAB path here the location
% of your gptoolbox copy:
addpath_recurse('/mnt/data/code/gptoolbox')
addpath('TexturePatch')
addpath('DECLab')
addpath('RicciFlow_MATLAB')
addpath(fullfile('utility','plotting'))
addpath(fullfile('utility','plotting'))
% go back to the data
cd(dataDir)


%%
disp('loading xp struct from disk')
load(fullfile(dataDir, 'xp.mat'), 'xp', 'opts')

%% TubULAR class instance
disp('defining TubULAR class instance (tubi= tubular instance)')
tubi = TubULAR(xp, opts) ;
disp('done defining TubULAR instance')

%% Extract 3D positions of nuclei from training from SPhi pullbacks
fn = fullfile(tubi.dir.tracking, 'tracks2d.mat') ;

if ~exist(fn, 'file') || overwrite
    fnim = tubi.fullFileBase.im_sp_sme ;
    for tp = tubi.xp.fileMeta.timePoints
        im = imread(sprintf(fnim, tp)) ;
        % extract just the green channel (nuclei)
        im = im(:, :, 2) ;

        % Segment the nuclei
        bw = imbinarize(im) ;
        s = regionprops(bw, 'centroid') ;
        centroids = cat(1,s.Centroid);
        % centroids = centroids(centroids(:, 2) < size(im, 1)*.75 ...
        %     & centroids(:, 2) > size(im, 1)*0.25, :) ;

        % Check it
        imshow(im); hold on; 
        plot(centroids(:, 1), centroids(:, 2), 'o')

        % Track the nuclei over time
        if tp > tubi.xp.fileMeta.timePoints(1)        
            % Check if we have the same number of nuclei over time
            % assert(size(centroids, 1) == nnuclei) ;

            % append the centroids to the running list
            clist = cat(1, clist, ...
                cat(2, centroids, tp*ones(size(centroids(:, 1)))));
        else
            nnuclei = size(centroids, 1) ;
            clist = cat(2, centroids, tp*ones(size(centroids(:, 1)))) ;
        end
    end

    % unscramble the tracks
    disp('unscrambling tracks using Crocker-Grier')
    tracks_CG = track(clist, 60) ;
    disp('unscrambling tracks using nearest neighbor search')
    tracks_NN = trackNearest(clist, 60) ;

    % Plot result -- CrockerGrier
    clf
    im = imread(sprintf(fnim, tubi.xp.fileMeta.timePoints(1))) ;
    im = im(:, :, 2)  ;
    im = im - min(im(:)) ;
    im = 255 - im ; 
    % im = 2 * double(im) / double(max(im(:))) * 255 ;
    % im(im > 255) = 255 ;
    % im = uint8(im) ;
    % im = im + (255 - max(im(:))) ;
    % im = uint
    imshow(im); hold on;
    tracks2d_CGCell = {} ;
    dmyk = 1 ;
    for id = 1:max(tracks_CG(:, 4))
        track_i = tracks_CG(tracks_CG(:, 4) == id, :) ;
        if track_i(1, 2) < size(im, 1)*.75 ...
             && track_i(1, 2) > size(im, 1)*0.25 
            % collate this track with others in a cell
            tracks2d_CGCell{dmyk} = track_i ;
            plot(track_i(:, 1), track_i(:, 2), '-')
            % waitforbuttonpress

            dmyk = dmyk + 1 ;
        end
    end
    set(gcf, 'color', 'w')
    mkdir(tubi.dir.tracking)
    export_fig(fullfile(tubi.dir.tracking, 'tracks_CG.png'), '-r600')
    
    % Plot result -- NearestNeighbor
    clf
    im = imread(sprintf(fnim, tubi.xp.fileMeta.timePoints(1))) ;
    im = im(:, :, 2) ;
    im = im - min(im(:)) ;
    im = 255 - im ; 
    imshow(im); hold on;
    tracks2d_NNCell = {} ;
    dmyk = 1 ;
    for id = 1:max(tracks_NN(:, 4))
        track_i = tracks_NN(tracks_NN(:, 4) == id, :) ;
        if track_i(1, 2) < size(im, 1)*.75 ...
             && track_i(1, 2) > size(im, 1)*0.25 
            % collate this track with others in a cell
            tracks2d_NNCell{dmyk} = track_i ;
            plot(track_i(:, 1), track_i(:, 2), '-')
            % waitforbuttonpress

            dmyk = dmyk + 1 ;
        end
    end
    set(gcf, 'color', 'w')
    export_fig(fullfile(tubi.dir.tracking, 'tracks_NN.png'), '-r600')

    save(fn, 'tracks2d_CGCell', 'tracks_CG', 'tracks2d_NNCell', 'tracks_NN')

else
    load(fn, 'tracks2d_CGCell', 'tracks_CG', 'tracks2d_NNCell', 'tracks_NN')
end


%% Convert to 3d tracks cell arrays
fn = fullfile(tubi.dir.tracking, 'tracks3d.mat') ;

if ~exist(fn, 'file') || overwrite 
    
    % Crocker-Grier tracking in 2D --> push to 3D
    tracks3d_CG = cell(size(tracks2d_CGCell)) ;
    for id = 1:length(tracks2d_CGCell)
        disp(['pushing forward into 3d: track ' ...
            num2str(id) '/' num2str(length(tracks2d_CGCell))])
        track_i = tracks2d_CGCell{id} ;
        first = true ;

        % extract the 3d position of the current track at each timepoint
        for tp = tubi.xp.fileMeta.timePoints
            tubi.setTime(tp);

            % get xy position
            trackpt = track_i(track_i(:, 3) == tp, 1:2) ;

            if ~isempty(trackpt)

                % convert mesh uv to pixel space
                mesh = tubi.getCurrentSPCutMeshSm ;

                umax = max(mesh.u(:, 1)) ;
                vmax = 1.0 ;
                doubleCovered = true ;
                XY = tubi.uv2XY(im, mesh.u, doubleCovered, umax, vmax) ;        
                pt3d = interpolate2Dpts_3Dmesh(mesh.f, XY, ...
                    mesh.v, trackpt) ;
                assert(all(size(pt3d) == [1 3]))

                if first
                    track3d = [pt3d, tp] ;
                    first = false ;
                else
                    track3d = cat(1, track3d, [pt3d, tp]) ;
                end
            end
        end
        tracks3d_CG{id} = track3d ;
    end
    
    % nearest neighbor tracking in 3D
    tracks3d_NNCell = cell(size(tracks2d_NNCell)) ;
    for id = 1:length(tracks2d_NNCell)
        disp(['pushing forward into 3d: track ' ...
            num2str(id) '/' num2str(length(tracks2d_NNCell))])
        track_i = tracks2d_NNCell{id} ;
        first = true ;

        % extract the 3d position of the current track at each timepoint
        for tp = tubi.xp.fileMeta.timePoints
            tubi.setTime(tp);

            % get xy position
            trackpt = track_i(track_i(:, 3) == tp, 1:2) ;

            if ~isempty(trackpt) && all(~isnan(trackpt))

                % convert mesh uv to pixel space
                mesh = tubi.getCurrentSPCutMeshSm ;

                umax = max(mesh.u(:, 1)) ;
                vmax = 1.0 ;
                doubleCovered = true ;
                XY = tubi.uv2XY(im, mesh.u, doubleCovered, umax, vmax) ;        
                pt3d = interpolate2Dpts_3Dmesh(mesh.f, XY, ...
                    mesh.v, trackpt) ;
                assert(all(size(pt3d) == [1 3]))

                if first
                    track3d = [pt3d, tp] ;
                    first = false ;
                else
                    track3d = cat(1, track3d, [pt3d, tp]) ;
                end
            end
        end
        tracks3d_NNCell{id} = track3d ;
    end
    save(fn, 'tracks3d_CGCell', 'tracks3d_NNCell')
else
    load(fn, 'tracks3d_CGCell', 'tracks3d_NNCell')
end


%% Match experiment results for Crocker-Grier (CG) to groundTruth positions
fn = fullfile(tubi.dir.tracking, 'exptTracks_groundTruth_CG.mat') ;
if ~exist(fn, 'file') || overwrite    
    % Build an array of 3d points over time from ground truth
    groundTruth = [] ;
    for tp = tubi.xp.fileMeta.timePoints
        tmp = load(fullfile(tubi.dir.data, ...
            sprintf('nuclei_positions_%06d.mat', tp))) ;
        % T x ptID x space array
        groundTruth(tp, :,:) = tmp.p3d ;
    end
    
    % point match tracks to initial point cloud
    t1 = tubi.xp.fileMeta.timePoints(1) ;
    tmp = load(fullfile(tubi.dir.data, ...
        sprintf('nuclei_positions_%06d.mat', t1))) ;
    p3d_init = tmp.p3d ;

    % experiment Tracks are T x pID x space dimensions
    nTracks = size(groundTruth, 2) ;
    nTimePoints = length(tubi.xp.fileMeta.timePoints) ;
    exptTracks_CG = zeros(nTimePoints, size(p3d_init, 1), size(p3d_init, 2)) ;

    % look for any tracks with tp=t1 and point match them to p3d_init
    clf
    % window size of expts with indices into groundTruth
    gtIDs_CG = zeros(length(tracks3d_CGCell), 1) ;   
    % window size of groundTruth with indices into Expts
    exIDs_CG = zeros(size(groundTruth, 2), 1) ;
    for id = 1:length(tracks3d_CGCell)
        track_i = tracks3d_CGCell{id} ;
        % if the track spans the whole range of timepoints, add it
        if track_i(1, 4) == t1 && numel(track_i(:, 4)) == nTimePoints
            % this is the ground truth ptID to assign tracks3d id to.
            gtIDs_CG(id) = pointMatch(track_i(1, 1:3), p3d_init) ;
            
            % index into original tracks3D cell for ith track is tID0
            exIDs_CG(gtIDs_CG(id)) = id ;
            exptTracks_CG(track_i(:, 4), gtIDs_CG(id), :) = track_i(:, 1:3) ;

            % check it
            plot3(exptTracks_CG(:, gtIDs_CG(id), 1), exptTracks_CG(:, gtIDs_CG(id), 2), ...
                exptTracks_CG(:, gtIDs_CG(id), 3), '.-');
            hold on;
            % end
        end
    end


    % check it
    for tp = tubi.xp.fileMeta.timePoints
        plot3(exptTracks_CG(tp, :, 1), exptTracks_CG(tp, :, 2), ...
            exptTracks_CG(tp, :, 3), '.');
        axis equal
        pause(0.01)
    end

    % bad tracks
    badtracks_CG = find(exptTracks_CG == 0 | isnan(exptTracks_CG) ) ;
    [r,c, l] = ind2sub(size(exptTracks_CG), badtracks_CG) ;
    badtracks_CG = unique(c) ;
    goodTracks_CG = setdiff(1:nTracks, badtracks_CG) ;

    save(fn, 'exptTracks_CG', 'groundTruth', 'gtIDs_CG',...
        'exIDs_CG', 'goodTracks_CG', 'nTimePoints', 'nTracks')
else
    load(fn, 'exptTracks_CG', 'groundTruth', 'gtIDs_CG', ...
        'exIDs_CG', 'goodTracks_CG', 'nTimePoints', 'nTracks')
end
disp('done')
    

%% Match experiment results for nearest neighbor (NN) to groundTruth positions
fn = fullfile(tubi.dir.tracking, 'exptTracks_groundTruth_NN.mat') ;
if ~exist(fn, 'file') || overwrite    
    % point match tracks to initial point cloud
    t1 = tubi.xp.fileMeta.timePoints(1) ;
    tmp = load(fullfile(tubi.dir.data, ...
        sprintf('nuclei_positions_%06d.mat', t1))) ;
    p3d_init = tmp.p3d ;

    % experiment Tracks are T x pID x space dimensions
    nTracks = size(groundTruth, 2) ;
    nTimePoints = length(tubi.xp.fileMeta.timePoints) ;
    exptTracks_NN = zeros(nTimePoints, size(p3d_init, 1), size(p3d_init, 2)) ;

    % look for any tracks with tp=t1 and point match them to p3d_init
    clf
    % window size of expts with indices into groundTruth
    gtIDs_NN = zeros(length(tracks3d_NNCell), 1) ;   
    % window size of groundTruth with indices into Expts
    exIDs_NN = zeros(size(groundTruth, 2), 1) ;
    for id = 1:length(tracks3d_NNCell)
        track_i = tracks3d_NNCell{id} ;
        % if the track spans the whole range of timepoints, add it
        if track_i(1, 4) == t1 && numel(track_i(:, 4)) == nTimePoints
            % this is the ground truth ptID to assign tracks3d id to.
            gtIDs_NN(id) = pointMatch(track_i(1, 1:3), p3d_init) ;
            
            % % check if we've already pointmatched to this one
            % if tID0(tID(id)) > 0
            %     % determine length of track
            % else
            %     add_track = true ;
            % end
            % if add_track

            % index into original tracks3D cell for ith track is tID0
            exIDs_NN(gtIDs_NN(id)) = id ;
            exptTracks_NN(track_i(:, 4), gtIDs_NN(id), :) = track_i(:, 1:3) ;

            % check it
            plot3(exptTracks_NN(:, gtIDs_NN(id), 1), exptTracks_NN(:, gtIDs_NN(id), 2), ...
                exptTracks_NN(:, gtIDs_NN(id), 3), '.-');
            hold on;
            % end
        end
    end


    % check it
    for tp = tubi.xp.fileMeta.timePoints
        plot3(exptTracks_NN(tp, :, 1), exptTracks_NN(tp, :, 2), ...
            exptTracks_NN(tp, :, 3), '.');
        axis equal
        pause(0.01)
    end

    % bad tracks
    badtracks_NN = find(exptTracks_NN == 0 | isnan(exptTracks_NN) ) ;
    [r,c, l] = ind2sub(size(exptTracks_NN), badtracks_NN) ;
    badtracks_NN = unique(c) ;
    goodTracks_NN = setdiff(1:nTracks, badtracks_NN) ;

    save(fn, 'exptTracks_NN', 'groundTruth', 'gtIDs_NN',...
        'exIDs_NN', 'goodTracks_NN')
else
    load(fn, 'exptTracks_NN', 'groundTruth', 'gtIDs_NN', ...
        'exIDs_NN', 'goodTracks_NN')
end
disp('done')

%% Extract tracks from ground Truth data in 3D
fn = fullfile(tubi.dir.tracking, 'tracks3D_groundTruth.mat') ;
if ~exist(fn, 'file') || overwrite
    first = true ;
    for tp = tubi.xp.fileMeta.timePoints
        tmp = load(fullfile(tubi.dir.data, ...
            sprintf('nuclei_positions_%06d.mat', tp))) ;
        % T x ptID x space array
        tptrack = cat(2, tmp.p3d, tp*ones(size(tmp.p3d(:, 1)))) ;
        if first 
            trueTracks = tptrack ;
            first = false; 
        else
            trueTracks = cat(1, trueTracks, tptrack) ;
        end
    end
    % Perform tracking
    tracks3d_gtCG = track(trueTracks, 6.5) ;
    tracks3d_gtNN = trackNearest(trueTracks) ;
    
    tracks3d_gtCGCell = {} ;
    for id = 1:max(tracks3d_gtCG(:, end))
        track_i = tracks3d_gtCG(tracks3d_gtCG(:, end) == id, :) ; 
        % collate this track with others in a cell
        tracks3d_gtCGCell{id} = track_i ;
    end
    
    tracks3d_gtNNCell = {} ;
    for id = 1:max(tracks3d_gtNN(:, end))
        track_i = tracks3d_gtNN(tracks3d_gtNN(:, end) == id, :) ; 
        % collate this track with others in a cell
        tracks3d_gtNNCell{id} = track_i ;
    end
    
    save(fn, 'tracks3d_gtCG', 'tracks3d_gtNN', ...
        'tracks3d_gtCGCell', 'tracks3d_gtNNCell')
else
    load(fn, 'tracks3d_gtCG', 'tracks3d_gtNN', ...
        'tracks3d_gtCGCell', 'tracks3d_gtNNCell')
end
    

% check result
% plot(tracks3D_nearest(:, 5), tracks3D_nearest(:, 4) .* ~isnan(tracks3D_nearest(:, 1)), '.-')

disp('done')

%% Match gt results for Crocker-Grier (CG) to groundTruth positions
fn = fullfile(tubi.dir.tracking, 'groundTruthTracks_groundTruth_CG.mat') ;
if ~exist(fn, 'file') || overwrite   
    
    % point match tracks to initial point cloud
    t1 = tubi.xp.fileMeta.timePoints(1) ;
    tmp = load(fullfile(tubi.dir.data, ...
        sprintf('nuclei_positions_%06d.mat', t1))) ;
    p3d_init = tmp.p3d ;

    % experiment Tracks are T x pID x space dimensions
    nTracks = size(groundTruth, 2) ;
    nTimePoints = length(tubi.xp.fileMeta.timePoints) ;
    gtTracks_CG = zeros(nTimePoints, size(p3d_init, 1), size(p3d_init, 2)) ;

    % look for any tracks with tp=t1 and point match them to p3d_init
    % window size of expts with indices into groundTruth
    gtIDs_gtCG = zeros(length(tracks3d_gtCGCell), 1) ;   
    % window size of groundTruth with indices into Expts
    exIDs_gtCG = zeros(size(groundTruth, 2), 1) ;
    for id = 1:length(tracks3d_gtCGCell)
        track_i = tracks3d_gtCGCell{id} ;
        % if the track spans the whole range of timepoints, add it
        if ~isempty(track_i)
            if track_i(1, 4) == t1 && numel(track_i(:, 4)) == nTimePoints
                % this is the ground truth ptID to assign tracks3d id to.
                gtIDs_gtCG(id) = pointMatch(track_i(1, 1:3), p3d_init) ;

                % index into original tracks3D cell for ith track is tID0
                exIDs_gtCG(gtIDs_gtCG(id)) = id ;
                gtTracks_CG(track_i(:, 4), gtIDs_gtCG(id), :) = track_i(:, 1:3) ;

                % check it
                plot3(gtTracks_CG(:, gtIDs_gtCG(id), 1), gtTracks_CG(:, gtIDs_gtCG(id), 2), ...
                    gtTracks_CG(:, gtIDs_gtCG(id), 3), '.-');
                hold on;
                % end
            end
        end
    end

    % bad tracks
    badtracks_CG = find(gtTracks_CG == 0 | isnan(gtTracks_CG) ) ;
    [r,c, l] = ind2sub(size(gtTracks_CG), badtracks_CG) ;
    badtracks_CG = unique(c) ;
    goodTracks_gtCG = setdiff(1:nTracks, badtracks_CG) ;

    save(fn, 'gtTracks_CG', 'groundTruth', 'gtIDs_gtCG',...
        'exIDs_gtCG', 'goodTracks_gtCG')
else
    load(fn, 'gtTracks_CG', 'groundTruth', 'gtIDs_gtCG', ...
        'exIDs_gtCG', 'goodTracks_gtCG')
end
disp('done')

%% Match gt results for Nearest-Neighbor (NN) to groundTruth positions
fn = fullfile(tubi.dir.tracking, 'groundTruthTracks_groundTruth_NN.mat') ;
if ~exist(fn, 'file') || overwrite   
    
    % point match tracks to initial point cloud
    t1 = tubi.xp.fileMeta.timePoints(1) ;
    tmp = load(fullfile(tubi.dir.data, ...
        sprintf('nuclei_positions_%06d.mat', t1))) ;
    p3d_init = tmp.p3d ;

    % experiment Tracks are T x pID x space dimensions
    nTracks = size(groundTruth, 2) ;
    nTimePoints = length(tubi.xp.fileMeta.timePoints) ;
    gtTracks_NN = zeros(nTimePoints, size(p3d_init, 1), size(p3d_init, 2)) ;

    % look for any tracks with tp=t1 and point match them to p3d_init
    % window size of expts with indices into groundTruth
    gtIDs_gtNN = zeros(length(tracks3d_gtNNCell), 1) ;   
    % window size of groundTruth with indices into Expts
    exIDs_NN = zeros(size(groundTruth, 2), 1) ;
    for id = 1:length(tracks3d_gtNNCell)
        track_i = tracks3d_gtNNCell{id} ;
        % if the track spans the whole range of timepoints, add it
        if track_i(1, 4) == t1 && numel(track_i(:, 4)) == nTimePoints
            % this is the ground truth ptID to assign tracks3d id to.
            gtIDs_gtNN(id) = pointMatch(track_i(1, 1:3), p3d_init) ;
            
            % index into original tracks3D cell for ith track is tID0
            exIDs_gtNN(gtIDs_gtNN(id)) = id ;
            gtTracks_NN(track_i(:, 4), gtIDs_gtNN(id), :) = track_i(:, 1:3) ;

            % check it
            plot3(gtTracks_NN(:, gtIDs_gtNN(id), 1), gtTracks_NN(:, gtIDs_gtNN(id), 2), ...
                gtTracks_NN(:, gtIDs_gtNN(id), 3), '.-');
            hold on;
            % end
        end
    end

    % bad tracks
    badtracks_NN = find(gtTracks_NN == 0 | isnan(gtTracks_NN) ) ;
    [r,c, l] = ind2sub(size(gtTracks_NN), badtracks_NN) ;
    badtracks_NN = unique(c) ;
    goodTracks_gtNN = setdiff(1:nTracks, badtracks_NN) ;

    save(fn, 'gtTracks_NN', 'groundTruth', 'gtIDs_gtNN',...
        'exIDs_gtNN', 'goodTracks_gtNN')
else
    load(fn, 'gtTracks_NN', 'groundTruth', 'gtIDs_gtNN', ...
        'exIDs_gtNN', 'goodTracks_gtNN')
end

%% Plot number of tracks
close all
figure('Position', [0, 0, 200, 300])
colors = define_colors ;

nParticles =  size(groundTruth, 2) ;

bar(1, numel(goodTracks_gtCG)/ nParticles, ...
    'facecolor', colors(1,:), 'edgecolor', 'none')
hold on;
bar(2, numel(goodTracks_gtNN)/ nParticles, ...
    'facecolor', colors(2,:), 'edgecolor', 'none')
bar(3, numel(goodTracks_CG)/ nParticles, ...
    'facecolor', colors(3,:), 'edgecolor', 'none')
bar(4, numel(goodTracks_NN)/ nParticles, ...
    'facecolor', colors(4,:), 'edgecolor', 'none')

ypos = 1.15 ;
text(1, ypos, num2str(numel(goodTracks_gtCG)), 'horizontalAlignment', 'center')
text(2, ypos, num2str(numel(goodTracks_gtNN)), 'horizontalAlignment', 'center')
text(3, ypos, num2str(numel(goodTracks_CG)), 'horizontalAlignment', 'center')
text(4, ypos, num2str(numel(goodTracks_NN)), 'horizontalAlignment', 'center')

hold on;
plot([0,5], [1,1], 'k--')
xlim([0,5])
ylim([0, 1.3])
xticks(1:4)

set(gca, 'fontsize', 8)
xticklabels({'3D tracking - CG', '3D tracking - NN', ...
    'tracking with tubular - CG', ...
    'tracking with tubular - NN'})
ylabel('fraction of cells tracked throughout entire timecourse', ...
    'fontsize', 10)
axis square
xtickangle(45)
set(gcf, 'color', 'w')
export_fig(fullfile(tubi.dir.tracking, 'number_of_tracks.pdf'), '-nocrop')

%% plot length of tracks
close all
figure('Position', [0, 0, 200, 300])
Ltracks_gtCG = zeros(length(tracks3d_gtCGCell), 1) ;
for tid = 1:max(tracks3d_gtCG(:, end))
    track_ii = tracks3d_gtCGCell{tid}(:,1) ;
    Ltracks_gtCG(tid) = numel(find(~isnan(track_ii))) ;
end
Ltracks_gtNN = zeros(length(tracks3d_gtNNCell), 1) ;
for tid = 1:length(tracks3d_gtNNCell)
    track_ii = tracks3d_gtNNCell{tid}(:,1) ;
    Ltracks_gtNN(tid) = numel(find(~isnan(track_ii))) ;
end
Ltracks_exCG = zeros(length(tracks3d_CGCell), 1) ;
for tid = 1:length(tracks3d_CGCell)
    track_ii = tracks3d_CGCell{tid}(:,1) ;
    Ltracks_exCG(tid) = numel(find(~isnan(track_ii))) ;
end
Ltracks_exNN = nan(length(tracks3d_NNCell), 1) ;
for tid = 1:length(tracks3d_NNCell)
    track_ii = tracks3d_NNCell{tid}(:,1) ;
    Ltracks_exNN(tid) = numel(find(~isnan(track_ii))) ;
end

colors = define_colors ;
% histogram for length of tracks with CG search in GroundTruth
[counts, edges] = histcounts(Ltracks_gtCG, ...
    [min(Ltracks_gtCG)-0.5:max(Ltracks_gtCG)+0.5]) ;
midy = (edges(1:end-1) + edges(2:end))/2 ;
xx = ones(numel(counts),1) ;
boxyViolinPlot(xx, midy, counts / sum(counts),...
    1, 'center', colors(1, :), 'none')
% histogram for length of tracks with nearest search in GroundTruth
[counts, edges] = histcounts(Ltracks_gtNN, ...
    (min(Ltracks_gtNN)-0.5):(max(Ltracks_gtNN)+0.5)) ;
midy = (edges(1:end-1) + edges(2:end))/2 ;
xx = ones(numel(counts),1)+1 ;
boxyViolinPlot(xx, midy, counts / sum(counts),...
    1, 'center', colors(2, :), 'none')

% histogram for length of tracks with CG search in expt
[counts, edges] = histcounts(Ltracks_exCG, ...
    (min(Ltracks_exCG)-0.5):(max(Ltracks_exCG)+0.5)) ;
midy = (edges(1:end-1) + edges(2:end))/2 ;
xx = ones(numel(counts),1)+2 ;
boxyViolinPlot(xx, midy, counts / sum(counts),...
    1, 'center', colors(3, :), 'none')

% histogram for length of tracks with nearest search in expt
[counts, edges] = histcounts(Ltracks_exNN, ...
    (min(Ltracks_exNN)-0.5):(max(Ltracks_exNN)+0.5)) ;
midy = (edges(1:end-1) + edges(2:end))/2 ;
xx = ones(numel(counts),1)+3 ;
boxyViolinPlot(xx, midy, counts / sum(counts),...
    1, 'center', colors(4, :), 'none')
hold on;
plot([0,5], nTimePoints*[1,1], 'k--')
xticks(1:4)
xlim([0,5])
set(gca, 'fontsize', 8)
xticklabels({'3D tracking - CG', '3D tracking - NN', ...
    'tracking with tubular - CG', ...
    'tracking with tubular - NN'})
ylabel('Length of tracks [timepoints]', 'fontsize', 10)
axis square 
xtickangle(45)
set(gcf, 'color', 'w')
export_fig(fullfile(tubi.dir.tracking, 'length_of_tracks.pdf'), '-nocrop')


clearvars first tptrack tmp


%% Check correlation
close all
fig = figure('position', [0,0, 400, 400]) ;
xx = exptTracks_CG(:, goodTracks_CG, 1) ;
yy = exptTracks_CG(:, goodTracks_CG, 2) ;
zz = exptTracks_CG(:, goodTracks_CG, 3) ;
XX = groundTruth(:, goodTracks_CG, 1); 
YY = groundTruth(:, goodTracks_CG, 2); 
ZZ = groundTruth(:, goodTracks_CG, 3); 
subplot(2, 2, 1)
s = scatter(XX(:), xx(:), 5, 'filled') ;
set(gca, 'fontsize', 8)
xlabel('true cell x coordinate [pixels]', 'fontsize', 10)
ylabel('measured cell x coordinate [pixels]', 'fontsize', 10)
set(s, 'markerfacealpha', 0.2)
axis square
axis equal
subplot(2, 2, 2)
s = scatter(YY(:), yy(:), 5, 'filled') ;
set(gca, 'fontsize', 8)
xlabel('true cell y coordinate [pixels]', 'fontsize', 10)
ylabel('measured cell y coordinate [pixels]', 'fontsize', 10 )
set(s, 'markerfacealpha', 0.2)
axis square
axis equal
subplot(2, 2, 3)
s = scatter(ZZ(:), zz(:), 5, 'filled') ;
set(gca, 'fontsize', 8)
xlabel('true cell z coordinate [pixels]', 'fontsize', 10)
ylabel('measured cell z coordinate [pixels]', 'fontsize', 10)
set(s, 'markerfacealpha', 0.2)
corrx = corrcoef(xx(:), XX(:)) ;
corry = corrcoef(yy(:), YY(:)) ;
corrz = corrcoef(zz(:), ZZ(:)) ;
axis square
axis equal
set(gcf, 'color', 'w')
fnout = fullfile(tubi.dir.tracking, 'position_benchmark') ;
saveas(gcf, [fnout '.pdf'])
export_fig([fnout '.png'], '-r600')

% Plot error
close all
fig = figure('position', [0,0, 200, 200]) ;

dx = xx - XX ;
dy = yy - YY ;
dz = zz - ZZ ;
colors = [     0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250] ;
plot([0,4], [0,0], 'k--', 'linewidth', 0.5)
hold on;
violin([dx(:), dy(:), dz(:)],'x', [1 2 3], ...
          'facecolor',colors,'edgecolor','none',...
          'bw',1e-3, 'mc',[],'medc', 'k')
xticks([1,2,3])
xticklabels({'x', 'y', 'z'})
xlim([0.5, 3.5])
ylim([-2,2])
yticks(-3:3)
set(gca, 'fontsize', 8)
ylabel('positional error of nuclei [pixels]', 'fontsize', 10)
xlabel('dimension', 'fontsize', 10)
set(gcf, 'color', 'w')
axis square
axis equal
fnout = fullfile(tubi.dir.tracking, 'position_error') ;
saveas(gcf, [fnout '.pdf'])
saveas(gcf, [fnout '.fig'])

     
%% Error in velocities
vx = diff(xx, 1) ;
vy = diff(yy, 1) ;
vz = diff(zz, 1) ;
VX = diff(XX, 1) ;
VY = diff(YY, 1) ;
VZ = diff(ZZ, 1) ;
dv = cat(3, vx-VX, vy-VY, vz-VZ) ;
dv_norm = sqrt((vx - VX).^2 + (vy-VY).^2 + (vz-VZ).^2) ;
Vnorm = sqrt(VX.^2 + VY.^2 + VZ.^2) ;
Vrms = sqrt(sum(Vnorm(:).^2) / numel(Vnorm)) ;
dv_rms = sqrt(sum(dv_norm.^2, 2) / size(dv_norm, 2)) ;
V_rms = sqrt(sum(Vnorm.^2, 2) / size(Vnorm, 2)) ;
readme = struct('dv', '3D measured velocity minus 3D true velocity',...
    'dv_norm', 'Norm of velocity error for each tracked nucleus over time',...
    'Vnorm', 'Norm of each true velocity vector for each tracked nucleus over time',...
    'Vrms', 'Root-mean-square true velocities',...
    'dv_rms', 'Root-mean-square velocity errors as a function of timepoint',...
    'V_rms', 'Root-mean-square true velocities as a function of timepoint') ;

save(fullfile(tubi.dir.tracking, 'velocity_error.mat'), ...
    'dv', 'dv_norm', 'Vnorm', 'Vrms', 'dv_rms', 'V_rms', ...
    'readme')

% plot it
close all
fig = figure('position', [0,0, 600, 240]) ;
subplot(1, 3, 1)
histogram(dv_norm(:) ./ Vnorm(:))
xlim([0,0.5])
set(gca, 'fontsize', 8)
axis square
xlabel('|v_{meas}-v_{true}|/|v_{true}|') %, 'interpreter', 'latex')
ylabel('counts', 'fontsize', 10)

subplot(1, 3, 2)
histogram(dv_norm(:) ./ Vrms(:))
xlim([0,0.5])
axis square
xlabel('|v_{meas}-v_{true}|/v^{RMS}_{true}') %, 'interpreter', 'latex')
ylabel('counts', 'fontsize', 10)

subplot(1, 3, 3)
set(gcf, 'color', 'w')
plot(tubi.xp.fileMeta.timePoints(1:end-1), V_rms, '.-')
hold on;
plot(tubi.xp.fileMeta.timePoints(1:end-1), dv_rms, '.-')
legend({'v^{RMS}_{true}', '|v_{meas}-v_{true}|^{RMS}'})
set(gca, 'fontsize', 8)
xlabel('timepoint', 'fontsize', 10)
ylabel('RMS velocity [pixels/timepoint]', 'fontsize', 10)
axis square
fnout = fullfile(tubi.dir.tracking, 'velocity_error') ;
saveas(gcf, [fnout '.pdf'])


%% Save example timepoint of velocities in 3d
close all
fig = figure('position', [0,0, 1000, 300]) ;
tdemo = 6 ;
subplot(1, 3, 1)
quiver3(XX(tdemo, :), YY(tdemo,:), ZZ(tdemo, :), ...
    VX(tdemo, :), VY(tdemo, :), VZ(tdemo,:), 0, 'color', colors(1, :)) 
xlabel('x [pixels]')
ylabel('y [pixels]')
zlabel('z [pixels]')
axis equal
subplot(1, 3, 2)
quiver3(xx(tdemo, :), yy(tdemo,:), zz(tdemo, :), ...
    vx(tdemo, :), vy(tdemo, :), vz(tdemo,:),...
    0, 'color', colors(2,:)) 
xlabel('x [pixels]')
ylabel('y [pixels]')
zlabel('z [pixels]')
axis equal
subplot(1, 3, 3)
%s = quiver3(XX(tdemo, :), YY(tdemo,:), ZZ(tdemo, :), ...
%    vx(tdemo, :), vy(tdemo, :), vz(tdemo,:), 'linewidth', 1.2) 

quiver3(XX(tdemo, :), YY(tdemo,:), ZZ(tdemo, :), ...
    VX(tdemo, :), VY(tdemo, :), VZ(tdemo,:), 'color', [0, 0.4470, 0.7410]) 
hold on;
quiver3(xx(tdemo, :), yy(tdemo,:), zz(tdemo, :), ...
    vx(tdemo, :), vy(tdemo, :), vz(tdemo,:), ...
    'color', [0.8500, 0.3250, 0.0980]) 
xlabel('x [pixels]')
ylabel('y [pixels]')
zlabel('z [pixels]')
axis equal
set(gcf, 'color', 'w')
fnout = fullfile(tubi.dir.tracking, 'velocity_example') ;
export_fig( [fnout '.png'], '-r600', '-nocrop')









