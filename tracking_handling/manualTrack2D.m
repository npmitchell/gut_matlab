% function [tracks, trackGraph] = manualTrack2D(currentTracks, fileBase, timePoints, trackOutfn, nTracks)
% Manually track nTracks objects in 2D grayscale image sequence.
% May run as script.
%
% Parameters
% ----------
% currentTracks : cell array or empty
%   existing tracks to build off of
% timePoints : length #timepoints numeric array
%   timepoints in which to track
% trackOutfn : str path
%   where to periodically save the results as .mat
% nTracks : int
%   how many objects to track
% 
% Returns
% -------
% tracks : nTracks x 1 cell array of #timepoints x 2 float arrays
%   positions of each track over time. NaNs are permitted for frames with
%   no object (disappearances)
% trackGraph : digraph object
%   digraph representation of the tracking results with Nodes and Edges
%
% Todo
% ----
% incorporate graph structure for tracks
% 
% NPMitchel 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script version options
subdir = 'muscle_normalShiftn10_p05_n50_s1p00_lambda0p005_maxProj';
imDir = fullfile('./', subdir, 'muscle_imagestack_LUT') ;
trackOutfn = fullfile('./muscle_tracks.mat') ;

load('./muscle_tracks.mat', 'tracks')
timePoints = 1:60 ;
nTracks = 300 ;
fileBase = fullfile(imDir, 'Time_%06d_c1_stab_pbspsm_LUT.tif') ;


currentTracks = tracks ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Unpack inputs

%--------
% inserted here
load('./muscle_tracks.mat', 'tracks')
%--------

%%
pausetime = 0.5 ;
if isempty(tracks)
    tracks = cell(nTracks, 1) ;
end
bluecolor = [0 , 0.4470, 0.7410 ];
orange = [ 0.8500,  0.3250 , 0.0980 ];
lwidth = 3 ;
markerSize = 20 ;


% Consider each object
for ii = 1:nTracks
    
    % Build tracking cell for this object
    if length(tracks) > ii-1
        trackii = tracks{ii} ;
    else
        tracks{ii} = zeros(length(timePoints), 2) ;
        trackii = tracks{ii} ;
    end
    if isempty(trackii)
        trackii = zeros(length(timePoints), 2) ;
    end
    
    % Consider each timepoint
    tidx = 1 ;
    Xlim = [];  
    Ylim = [] ;
    recap = false ;
    while tidx < length(timePoints) + 1
        if ~any(trackii(tidx, :))
            recap = true ;
            tp = timePoints(tidx) ;
            
            % Initialize the Figure
            if ~ishandle(fig)
                fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
            end
            
            if ~exist('Xlim', 'var')
                Xlim = [];  
                Ylim = [] ;
            end
            
            [ax0, ax1, ax2, minX, maxX, minY, maxY] = plotCurrentTrack(ii, tidx, timePoints, fileBase, trackii, ...
                currentTracks, Xlim, Ylim) ;
            
            % Now acquire or get more info about nearby timepoints
            fprintf('Press any key to acquire or p to play')
            sgtitle('Press <space> to acquire, p: play, o: play from start, a: t-1, s: t+1')
            pause
            currkey=get(gcf, 'CurrentKey'); 
            
            
            while ismember(currkey, {'p', 'o', 'a', 's'}) 
                
                switch currkey
                    case {'a', 's'}
                        % Decrememt/Increment timepoint
                        if strcmpi(currkey, 'a')
                            tidx = tidx - 1;
                        elseif strcmpi(currkey, 's')
                            tidx = tidx + 1;
                        end
                        [ax0, ax1, ax2, minX, maxX, minY, maxY] = plotCurrentTrack(ii, tidx, timePoints, fileBase, trackii, ...
                        currentTracks, Xlim, Ylim) ;

                        % replay or move on to acquire
                        fprintf('Press any key to acquire or p to replay')
                        pause
                        currkey=get(gcf,'CurrentKey'); 
                    case {'p', 'o'}
                        % Play nearby timepoints to help identify cell center                
                        Xlim = get(gca, 'XLim');
                        Ylim = get(gca, 'YLim');

                        if strcmp(currkey, 'p')
                            tidx2play = max(1, tidx-4):min(tidx+2, length(timePoints)) ;
                        elseif strcmp(currkey, 'o')
                            tidx2play = 1:length(timePoints) ;
                        end
                        times2play = timePoints(tidx2play) ;
                        axes(ax2) ;
                        for kk = 1:length(times2play)
                            ttemp = times2play(kk) ;
                            im1 = imread(sprintf(fileBase, ttemp));

                            cla
                            imshow(im1) ;
                            hold on;
                            % Plot other lineages in blue
                            for trackID = 1:ii-1
                                plot(tracks{trackID}(tidx2play(kk), 1), ...
                                    tracks{trackID}(tidx2play(kk), 2), 'o', ...
                                    'markerSize', markerSize, ...
                                    'lineWidth', lwidth, 'color', bluecolor)
                            end
                            % Plot current lineage in orange
                            plot(trackii(tidx2play(kk), 1), ...
                                trackii(tidx2play(kk), 2), 'o', ...
                                'markerSize', markerSize, ...
                                'lineWidth', lwidth, 'color', orange)
                            if tidx > 1
                                set(gca, 'XLim', Xlim , 'YLim', Ylim);
                            end
                            title(['playing track ' num2str(ii) ': t=' num2str(ttemp)])

                            pause(pausetime)
                        end

                        % Re-gain original image onto axis
                        im = imread(sprintf(fileBase, tidx));
                        cla()
                        imshow(im) ;
                        title(['Track ' num2str(ii) ': t=' num2str(tp)])
                        set(gca, 'XLim', Xlim , 'YLim', Ylim);

                        % replay or move on to acquire
                        fprintf('Press any key to acquire or p to replay')
                        pause
                        currkey=get(gcf,'CurrentKey'); 
                    otherwise
                        disp('done with adjustments')
                
                end  % end of the switch
                
            end % end of the while loop
            
            % Acquire the XY coordinate for this timepoint
            msg = 'Click or <space> with crosshairs: acquire / <escape>: no detection\n' ;
            sgtitle(msg)
            disp(msg)
            [xx,yy] = ginput(1) ;
            
            
            currkey=get(gcf, 'CurrentKey'); 
            if strcmpi(currkey, 'escape') || strcmpi(currkey, 'backspace')
                trackii(tidx, :) = [NaN, NaN] ;
            else
                trackii(tidx, :) = [xx,yy] ;
            end
            
            Xlim = get(gca, 'XLim');
            Ylim = get(gca, 'YLim');
            
            % prep for next timepoint
            % Save tracks so far
            tracks{ii} = trackii ;
            disp(['saving developing tracks to : ' trackOutfn])
            save(trackOutfn, 'tracks')

            tidx = tidx + 1 ;
            
        else
            disp(['done with tidx = ' num2str(tidx)])
            tidx = tidx + 1 ;
        end
        
    end
    tracks{ii} = trackii ;
    currentTracks = tracks ;
    % Save tracks
    disp(['saving completed tracks to : ' trackOutfn])
    save(trackOutfn, 'tracks')
    
    
    %% Instant replay
    if recap
        clf
        set( fig, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        times2play = timePoints ;
        for kk = 1:length(times2play)
            ttemp = times2play(kk) ;
            im1 = imread(sprintf(fileBase, ttemp));
            imshow(im1) ;
            hold on;

            for iprev = 1:ii-1
                trackprev = currentTracks{iprev} ;
                plot(trackprev(kk, 1), trackprev(kk, 2), 'o', 'color', ...
                    bluecolor, 'markerSize', markerSize, 'lineWidth', lwidth)
            end

            tracknow = currentTracks{ii} ;
            plot(tracknow(kk, 1), tracknow(kk, 2), 'o', 'color', orange, ...
                'markerSize', markerSize, 'lineWidth', lwidth)

            title(['playing tracks 1-' num2str(ii) ': t=' num2str(ttemp)])

            pause(0.01 * pausetime)
        end
    end
    
    %% Review track
end

% Optionally convert to digraph
if nargout > 1
    trackGraph = trackingCell2Graph(tracks) ;
end

