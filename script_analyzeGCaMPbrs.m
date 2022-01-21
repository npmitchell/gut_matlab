%% Script for analyzing high-temporal-resolution calcium pulses with brs
clear all
addpath('./addpath_recurse/')
addpath_recurse('/mnt/data/code/gut_matlab/')

rootdir = ['/mnt/data/confocal_data/gut/mef2GAL4klarUASbrsChGCaMP6s/' , ... 
    '202201061700_gcampDynamics_0p3spf_0p75pcAt50pc488_0p75pc561_getPulses/'] ;
dirs = {'1659', '1700_1sPause', ...
    '1705_noPause', '1710_noPause', ...
    '1715_about30sPause', '1720_about1sPause' } ;
addFrames = [0, 3, 0, 0, 10, 3] ;
spf = 0.3 ; % seconds per frame
pix2um = 1. / 3.0842 ; % conversion from pixels: um / pix -->  0.3258147

% Filtering options
thres = 15 ;        % minimum intensity for a peak to be included
filt = fspecial('gaussian', 7,1) ; % filter to prevent salt-n-pepper noise
edg = 3 ;  % pixels from the edge within which we discard peaks
preview = false ;
overwrite = false ; % overwrite previous results

timestamp = 0 ;
tstamp_int = 1 ;
fnList0 = {} ;
fnList1 = {} ;
for ii = 1:length(dirs)
    ddir = dirs{ii} ;
    disp(['Considering dir ' num2str(ii) ': ' dirs{ii}])
    
    % Get filenames
    addF = addFrames(ii) ;
    timestamp = timestamp + addF * spf ;
    
    trackFn = fullfile(rootdir, [dirs{ii} '_tracks.mat']) ;
    fns0 = dir(fullfile(rootdir, ddir, '*C0*.tif')) ;
    fns1 = dir(fullfile(rootdir, ddir, '*C1*.tif')) ;
    % there should be equal number of channel 0 as channel 1 images
    assert(length(fns0) == length(fns1)) 
    if ~exist(trackFn, 'file') || overwrite
        
        for jj = 1:length(fns0)
            disp(['fn jj = ' num2str(jj)])
            fn0 = fullfile(fns0(jj).folder, fns0(jj).name) ;
            im0 = imread(fn0) ;
            fn1 = fullfile(fns1(jj).folder, fns1(jj).name) ;
            im1 = imread(fn1) ;
            
            % Append to filename lists
            fnList0{tstamp_int} = fn0 ;
            fnList1{tstamp_int} = fn1 ;

            % Find peaks using Channels' sum
            % im = im0 + im1 ;
            if jj == 1 
                maxIm1 = im1 ;
                if ii == 1 
                    firstFrame = im1 ;
                end
            else
                maxIm1 = max(maxIm1, im1) ;
            end
            
            [peaks, labelM1, brightness] = FastPeakFind(im1, thres, filt , edg, 2) ;

            % get mean intensity in each peak region in channel0 (GFP)
            % se = offsetstrel('ball',4,4);
            % labelM = imdilate(labelM0, se) ;
            stats = regionprops(labelM1, im0, 'Area', 'MeanIntensity');
            calcium = [stats.Area]' .* [stats.MeanIntensity]';

            % Check the peaks
            if preview 
                imshow(im0); 
                hold on; 
                plot(peaks(:, 1), peaks(:, 2), 'ro')
                caxis([0, 100])
                figure(2)
                imshowpair(im,labelM>0,'montage')
            end

            % Collate info into matrix
            pkInfo = cat(2, peaks, ...
                tstamp_int * ones(size(peaks, 1), 1), ...
                timestamp * ones(size(peaks, 1), 1), ...
                brightness, calcium) ;

            % save peaks
            ensureDir(fullfile(fns0(jj).folder, 'peaks')) ;
            fn0 = fullfile(fns0(jj).folder, 'peaks', ...
                [fns0(jj).name(1:end-3) '.txt']) ;
            header = ['peakX, peakY, frame#, timestamp, ', ...
                '561nm emission, GCaMP6s emission' ] ;
            write_txt_with_header(fn0, pkInfo, header)

            if ii == 1 && jj == 1
                % new peak arrays for both this directory and all time
                allPeaks = pkInfo ;
                thesePeaks = pkInfo ;
            elseif jj == 1
                % new peak array for this directory
                allPeaks = [allPeaks ; pkInfo] ;
                thesePeaks = pkInfo  ;
            else
                % no new peak arrays, cat instead
                allPeaks = [allPeaks ; pkInfo] ;
                thesePeaks = [thesePeaks; pkInfo ] ;
            end

            % Update timestamp
            timestamp = timestamp + spf ;
            tstamp_int = tstamp_int + 1;
        end

        %% Mean image
        imwrite(maxIm1, fullfile(rootdir, [dirs{ii} '_maxIm1.png']))        
        
        %% Track these dots for this dir
        maxdispl = 8 ;
        param = struct() ;
        param.mem = 100 ;
        param.dim = 2 ;
        param.quiet = false ;
        param.good = 2 ;
        % tracks = track(allPeaks(:, 1:3), maxdispl, param) ;
        thesePeaks2 = thesePeaks(:, [1, 2, 4, 5, 6, 3]) ;
        tracks = track(thesePeaks2, maxdispl, param) ;

        trackCell = {} ;
        for qq = 1:max(tracks(:, end))
            trackID = find(tracks(:, end) == qq) ;
            trackq = tracks(trackID, 1:end-1) ;
            trackCell{qq} = trackq ;
        end
        peaks_sequence = thesePeaks2 ;
        readme = struct() ;
        readme.peaks_sequence = ['ordered by frame: peakX, peakY, timestamp, ', ...
                '561nm emission, GCaMP6s emission, frame#' ] ;
        readme.tracks = ['ordered by track: peakX, peakY, timestamp, ', ...
                '561nm emission, GCaMP6s emission, frame#, track#' ] ;
        readme.trackCell = ['ordered by track: for each entry, values for a single track', ...
                ' are: peakX, peakY, timestamp, ', ...
                '561nm emission, GCaMP6s emission, frame#' ] ;
        save(trackFn, ...
            'peaks_sequence', 'tracks', 'trackCell', 'readme')

        % Plot all the tracks & save figure
        clf; hold on;
        for qq = 1:length(trackCell)
            plot(trackCell{qq}(:, 1),trackCell{qq}(:, 2)) ;
        end
        axis equal
        xlim([1, size(im0, 2)])
        ylim([1, size(im0, 1)])
        saveas(gcf, fullfile(rootdir, [dirs{ii} '_tracks.png']))
    else
        load(trackFn, 'peaks_sequence') ;
        % x, y, timestamp (s), GFP, RFP, tstamp_integer
        old_tstamp_int = tstamp_int ;
        timestamp = peaks_sequence(end, 3) ;
        tstamp_int = peaks_sequence(end, end) ;
        timestamp = timestamp + spf ;
        tstamp_int = tstamp_int + 1;
        
        % Add this peak sequence to allPeaks
        thesePeaks2 = peaks_sequence(:, [1, 2, 6, 3, 4, 5]) ;
        % x, y, tstamp_integer, timestamp (s), GFP, RFP, 
        if ii == 1
            allPeaks = thesePeaks2 ;
            fns1 = dir(fullfile(rootdir, ddir, '*C1*.tif')) ;
            fn1 = fullfile(fns1(1).folder, fns1(1).name) ;
            firstFrame = imread(fn1) ;
        else
            allPeaks = [allPeaks ; thesePeaks2] ;
        end
        
        % Append to filename lists
        for jj = 1:length(fns0)
            fn0 = fullfile(fns0(jj).folder, fns0(jj).name) ;
            fn1 = fullfile(fns1(jj).folder, fns1(jj).name) ;
            fnList0{old_tstamp_int + jj - 1} = fn0 ;
            fnList1{old_tstamp_int + jj - 1} = fn1 ;
        end
    end
end

% Track dots
maxdispl = 8 ;
param = struct() ;
param.mem = 100 ;
param.dim = 2 ;
param.quiet = false ;
param.good = 40 ;
% tracks = track(allPeaks(:, 1:3), maxdispl, param) ;
allPeaks2 = allPeaks(:, [1, 2, 4, 5, 6, 3]) ;
tracks = track(allPeaks2, maxdispl, param) ;

trackCell = {} ;
for qq = 1:max(tracks(:, end))
    trackID = find(tracks(:, end) == qq) ;
    trackq = tracks(trackID, :) ;
    trackCell{qq} = trackq ;
end
peaks_sequence = allPeaks2 ;
readme = struct() ;
readme.peaks_sequence = ['ordered by frame: peakX, peakY, timestamp, ', ...
        '561nm emission, GCaMP6s emission, frame#' ] ;
readme.tracks = ['ordered by track: peakX, peakY, timestamp, ', ...
        '561nm emission, GCaMP6s emission, frame#, track#' ] ;
readme.trackCell = ['ordered by track: for each entry, values for a single track', ...
        ' are: peakX, peakY, timestamp, ', ...
        '561nm emission, GCaMP6s emission, frame#' ] ;
save(fullfile(rootdir, 'tracks.mat'), 'tracks', 'peaks_sequence', 'trackCell')

% Plot all the tracks
clf; hold on;
for qq = 1:length(trackCell)
    plot(trackCell{qq}(:, 1),trackCell{qq}(:, 2)) ;
end
axis equal
xlim([1, size(firstFrame, 2)])
ylim([1, size(firstFrame, 1)])
drawnow
set(gcf, 'color', 'k')
set(gca, 'color', 'k')
FF = getframe() ;
imout = imresize(FF.cdata, size(firstFrame)) ;
imwrite(imout, fullfile(rootdir, 'tracks.png'))



%% MASK starting frame for dots of interest
for ii = 1:length(dirs)
    ddir = dirs{ii} ;
    disp(['Considering dir ' num2str(ii) ': ' dirs{ii}])
    if ii == 1
        maxim = imread( fullfile(rootdir, [dirs{ii} '_maxIm1.png'])) ;
    else
        maxim = max(maxim, imread( fullfile(rootdir, [dirs{ii} '_maxIm1.png']))) ;
    end
end

% Get ROI to keep tracks only within region of interest
im = imread( fullfile(rootdir, 'tracks.png')) ;
im = flipud(im) ;
imMerge = cat(3, maxim + im(:,:,1), maxim + im(:,:,2), maxim + im(:,:,3)) ;


% bw = (255*3 - sum(im, 3)) > 0 ;
% RR = im(:, :, 1) ;
% GG = im(:, :, 2) ;
% BB = im(:, :, 3) ;
% RR(RR==255 & GG==255 & BB==255) = 0 ;
roiFn = fullfile(rootdir, 'roiKeep.mat') ;
if exist(roiFn, 'file')
    disp('Loading ROI from disk')
    load(roiFn, 'roi', 'roix', 'roiy', 'keep')
else
    disp('Define ROI on image as polygon')
    [roi, roix, roiy] = roipoly(imMerge) ;
    keep = [] ;
    for qq = 1:length(trackCell)
        trackq = trackCell{qq} ;
        inside = inpolygon(trackq(:, 1), trackq(:, 2), roix, roiy) ;
        if mean(inside > 0.5)
            keep = [keep , qq] ;
        end
    end
    save(roiFn, 'roi', 'roix', 'roiy', 'keep')
end

%% Assemble pruned track list as cell, storing first x positions
dmyi = 1 ;
keepTracks = [] ;  % a matrix with X Y Time GFP RFP Frame# NewTrack#
keepTrackCell = {} ; % cell array of kept tracks: X Y Time GFP RFP Frame# NewTrack#
x0 = [] ; 
times = [] ;
gcamp = [] ;
gradgcamp = [] ;
firstXY = [] ;       % for all kept peaks, mark first appearance location XY
firstFrameXY = [] ;  % for kept peaks that appear in the first frame, mark XY
firstTP_threshold = 1000 ; % seconds, cutoff for first timepoint being early enough
for qq = 1:length(trackCell)
    % only keep track if inside roi and if the first timepoint for this
    % track is early enough
    if ismember(qq, keep) 
        % find if first timepoint is before threshold
        firstTP = trackCell{qq}(1, 3) ;
        if firstTP < firstTP_threshold
            thisTrack = trackCell{qq}(:, 1:6) ;
            keepTrackCell{dmyi} = [ thisTrack, dmyi*ones(size(thisTrack, 1), 1) ] ;
            if dmyi > 1
                keepTracks = [keepTracks ; ...
                    [ thisTrack, dmyi*ones(size(thisTrack, 1), 1) ] ] ;
            else
                keepTracks = [ thisTrack, dmyi*ones(size(thisTrack, 1), 1) ] ;                
            end
            
            % get initial x positions 
            x0 = [x0; trackCell{qq}(1, 1) * ones(size(thisTrack, 1), 1)] ;
            thisTimes = thisTrack(:, 3) ;
            times = [times; thisTimes ] ;
            thisGcamp = thisTrack(:, 5) ./ thisTrack(:, 4) ;
            gcamp = [gcamp ; thisGcamp ] ;
            thisGradgcamp = gradient(thisGcamp, thisTimes) ;
            gradgcamp = [gradgcamp ; thisGradgcamp ] ;
            firstXY = [firstXY ; thisTrack(1, 1:2) ] ;
            
            % Mark the XY coords if this peak is in the first frame (t=0)
            if any(thisTrack(:, 3) == 0)
                tmp = find(thisTrack(:, 3) == 0) ;
                assert(tmp == 1)
                assert(length(tmp) == 1)
                firstFrameXY = [firstFrameXY; thisTrack(tmp, 1:2)];
            else
                firstFrameXY = [firstFrameXY; [-Inf, -Inf]];                
            end
                    
            dmyi = dmyi + 1 ;
        end
    end
end
disp('done assembling pruned list')






%% Debugging with images present -- choose a peak
close all
imshow(firstFrame)
colormap greys
% how many points to track
nn = 6 ;
% pts = ginput(3) ;
pts =  [   145.3793  129.7817 ;...
  325.3063  252.5285 ;...
   76.9915  214.0066 ] ;
pts2add = ginput(3) ;
pts = [pts ; pts2add] ;
ids = [] ;
myTracks = {} ;
for pid = 1:length(pts)
    % Could choose closest out of all peaks ever (using firstXY)
    dist2pt = vecnorm(firstXY - pts(pid, :), 2, 2) ;
    % Instead choose closest out of the first frame's peaks
    % dist2pt = vecnorm(firstFrameXY - pts(pid, :), 2, 2) ;
    
    [~, ids(pid)] = min(dist2pt) ;
    myTracks{pid} = keepTrackCell{ids(pid)} ;
end

[colors, names] = define_colors ;
colors = colors ./ vecnorm(colors, 2, 2) ;
smoothW = 25 ;

close all
clf
ax1 = subplot(2, 2, 1) ;
ax2 = subplot(2, 2, 2) ;
ax3 = subplot(2, 2, 3) ;
ax4 = subplot(2, 2, 4) ;
% initial title
figtitle = sgtitle(['$t=$'], 'interpreter', 'latex') ;
q2do = 1:100:length(fnList0) ;
q2do = [q2do, setdiff(1:50:length(fnList0), q2do)] ;
q2do = [q2do, setdiff(1:20:length(fnList0), q2do)] ;
% q2do = [q2do, setdiff(1:length(fnList0), q2do)] ;
for qq = q2do
    im0 = imread(fnList0{qq})  ;
    set(gcf, 'CurrentAxes', ax1)
    cla
    imshow(im0)
    caxis([0, 200])
    % colormap inferno
    
    % Lower image -- RFP
    im1 = imread(fnList1{qq})  ;
    set(gcf, 'CurrentAxes', ax3)
    cla
    imshow(im1)
    caxis([0, 150])
    colormap greys
    
    % clear other subplots
    set(gcf, 'CurrentAxes', ax2)
    cla
    set(gcf, 'CurrentAxes', ax4)
    cla
    
    % identify peaks being tracked
    for pid = 1:length(pts)
        % Upper image -- GFP
        set(gcf, 'CurrentAxes', ax1)
        if ismember( qq, myTracks{pid}(:, 6) )
            indx = find(myTracks{pid}(:, 6) == qq) ;
            hold on;
            plot(myTracks{pid}(indx, 1), myTracks{pid}(indx, 2), ...
                'o', 'color', colors(pid, :)) ;
        else
            disp(['pt ' num2str(ids(pid)) ' is not part of frame ' num2str(qq)])
        end
        
        % Lower image -- RFP
        set(gcf, 'CurrentAxes', ax3)
        if ismember( qq, myTracks{pid}(:, 6) )
            indx = find(myTracks{pid}(:, 6) == qq) ;
            hold on;
            plot(myTracks{pid}(indx, 1), myTracks{pid}(indx, 2), ...
                'o', 'color', colors(pid, :)) ;
        else
            disp(['pt ' num2str(ids(pid)) ' is not part of frame ' num2str(qq)])
        end
    end
    
    % curve of intensity over time
    set(gcf, 'CurrentAxes', ax2)
    for pid = 1:length(pts)
        % smooth the RFP denominator
        denom = movmean(myTracks{pid}(:, 4), smoothW)  ;
        gcampnorm = myTracks{pid}(:, 5) ./ denom ;
        plot(myTracks{pid}(:, 3), ...
            gcampnorm, ...
            'color', colors(pid, :))
        hold on;
        if ismember( qq, myTracks{pid}(:, 6) )
            indx = find(myTracks{pid}(:, 6) == qq) ;
            plot(myTracks{pid}(indx, 3), ...
                gcampnorm(indx), 'o', ...
                'color', colors(pid, :))
        end
    end
    
    
    % Lower intensity vs time
    set(gcf, 'CurrentAxes', ax4)
    for pid = 1:length(pts)
        % smooth the RFP denominator
        denom = movmean(myTracks{pid}(:, 4), smoothW)  ;
        plot(myTracks{pid}(:, 3), denom, 'color', colors(pid, :)) ;
        hold on;
        if ismember( qq, myTracks{pid}(:, 6) )
            indx = find(myTracks{pid}(:, 6) == qq) ;
            plot(myTracks{pid}(indx, 3), ...
                denom(indx), 'o', ...
                'color', colors(pid, :))
        end
    end
    
    % Update the title
    % figtitle.String = ['$t=$' num2str(myTracks{pid}(indx, 3)/60.) ' min' ] ;
    figtitle.String = ['$t=$' num2str(qq * spf/60.) ' min' ] ;
    drawnow
    
    % Save the image
    ensureDir( fullfile(rootdir, 'demo'))
    saveas(gcf, fullfile(rootdir, 'demo', sprintf('t%06d.png', qq)))
end








%% Plot intensity variations

% interpolate into heatmap
x0jitter = ( x0 + 0.01*rand(size(x0)) ) * pix2um ;
si = scatteredInterpolant(x0jitter, times, gcamp) ;

% check it 
close all
scatter(x0jitter, times / 60.0, 3, gcamp, 'filled')
set(gca, 'color', 'k')
set(gcf, 'color', 'w')
cb = colorbar ;
ylabel(cb, '$I_{\mathrm{GCaMP6s}} / I_{brs}$', 'interpreter', 'latex')
xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
ylabel('time [min]', 'interpreter', 'latex')
saveas(gcf, fullfile(rootdir, 'peaks_kymo.png')) 

% check gradient
clf
scatter(x0jitter, times / 60.0, 3, gradgcamp, 'filled')
caxis([0, max(gradgcamp(:))])
set(gca, 'color', 'k')
cb = colorbar ;
ylabel(cb, '$\partial_t \left(I_{\mathrm{GCaMP6s}} / I_{brs}\right)$', 'interpreter', 'latex')
xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
ylabel('time [min]', 'interpreter', 'latex')
caxis([0, 10])
saveas(gcf, fullfile(rootdir, 'peaks_kymo_gradient.png')) 



