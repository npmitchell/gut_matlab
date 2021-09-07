function [tracks, trackGraph] = manualTrack2D(currentTracks, timePoints, trackOutfn, nTracks)
% Manually track nTracks objects in 2D grayscale image sequence
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

% Unpack inputs
pausetime = 0.5 ;
tracks = currentTracks ;
if isempty(tracks)
    tracks = cell(nTracks, 1) ;
end

% Consider each object
for ii = 1:nTracks
    % Build tracking cell for this object
    trackii = currentTracks{ii} ;
    if isempty(trackii)
        trackii = zeros(length(timePoints), 2) ;
    end
    
    % Consider each timepoint
    for tidx = 1:length(timePoints)

        if ~any(trackii(tidx, :))
            tp0 = timePoints(max(1, tidx-1)) ;
            tp = timePoints(tidx) ;

            % Overlay previous and current image
            im0 = imread(sprintf(fileBase, tp0));
            im = imread(sprintf(fileBase, tp));
            rgb = cat(3, im0, im, im) ;

            ax0 = subtightplot(2, 2, 1) ;
            imshow(im0) ;
            if tidx > 1
                set(gca, 'XLim', Xlim , 'YLim', Ylim);
                hold on;
                plot(trackii(1:tidx-1, 1), trackii(1:tidx-1, 2), '-')
                plot(trackii(tidx-1, 1), trackii(tidx-1, 2), 'o')
                hold off;
            end
            title(['t=' num2str(tp0)])
            
            ax1 = subtightplot(2, 2, 2) ;
            imshow(rgb);
            if tidx > 1
                set(gca, 'XLim', Xlim , 'YLim', Ylim);
                hold on;
                plot(trackii(1:tidx-1, 1), trackii(1:tidx-1, 2), '-')
                plot(trackii(tidx-1, 1), trackii(tidx-1, 2), 'o')
                hold off;
            end
            title('merge')
            
            ax2 = subtightplot(2, 1, 2) ;
            imshow(im) ;
            if tidx > 1
                set(gca, 'XLim', Xlim , 'YLim', Ylim);
            end
            title(['Track ' num2str(ii) ': t=' num2str(tp)])
            
            % Now acquire or get more info about nearby timepoints
            fprintf('Press any key to acquire or p to play')
            pause
            currkey=get(gcf,'CurrentKey'); 
            
            
            % Play nearby timepoints to help identify cell center
            while strcmp(currkey, 'p') 
                
                Xlim = get(gca, 'XLim');
                Ylim = get(gca, 'YLim');
            
                tidx2play = max(1, tidx-4):min(tidx+2, length(timePoints)) ;
                times2play = timePoints(tidx2play) ;
                axes(ax2) ;
                for kk = 1:length(times2play)
                    ttemp = times2play(kk) ;
                    im1 = imread(sprintf(fileBase, ttemp));

                    subtightplot(2, 1, 2)
                    imshow(im1) ;
                    hold on;
                    plot(trackii(tidx2play(kk), 1), ...
                        trackii(tidx2play(kk), 2), 'o')
                    if tidx > 1
                        set(gca, 'XLim', Xlim , 'YLim', Ylim);
                    end
                    title(['playing track ' num2str(ii) ': t=' num2str(ttemp)])
                    
                    pause(pausetime)
                end
                
                % Re-gain original image onto axis
                imshow(im) ;
                title(['Track ' num2str(ii) ': t=' num2str(tp)])
                set(gca, 'XLim', Xlim , 'YLim', Ylim);
                        
                % replay or move on to acquire
                fprintf('Press any key to acquire or p to replay')
                pause
                currkey=get(gcf,'CurrentKey'); 
            end
            
            % Acquire the XY coordinate for this timepoint
            fprintf('Click on cell center or hit key for NaN \n')
            [xx,yy] = ginput(1) ;
            Xlim = get(gca, 'XLim');
            Ylim = get(gca, 'YLim');
            trackii(tidx, :) = [xx,yy] ;
            
            % prep for next timepoint
        end
    end
    tracks{ii} = trackii ;
    % Save tracks
    disp(['saving tracks to : ' trackOutfn])
    save(trackOutfn, 'tracks')
end

% Optionally convert to digraph
if nargout > 1
    trackGraph = trackingCell2Graph(tracks) ;
end

