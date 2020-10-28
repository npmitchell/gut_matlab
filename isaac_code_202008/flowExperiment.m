%% Define a class which contains methods for processing and analyzing flow data
classdef flowExperiment < handle
    % Flows are measured in microns per dt, where dt is the time
    % resolution. To compare different datasets, rescale time by dt 
    %
    properties
        
        % Properties that need to be set by the user
        dataDir
        dataSetName
        lastTimePoint
        fileSize
        trackFileName
        resolution
        dt                  % time resolution in minutes
        lobeSplit
        
        % Properties that can be deduced from set properties
        subSamplingFactor
        rawDataFileName
        probabilitiesFileName
        maskedDataFileName
        maskedMIPFileName
        ch1FileName
        ch2FileName
        ch1FlowFileName
        ch2FlowFileName
        
        % Properties that will be created through methods
        trackData
        trackVals
        fullFlowField
        binnedVelocities
        binnedVelocities2
        relativeVelocities
        relativeVelocities2
        nuclearCenters
        nuclearCenters2
        nuclearPhases
        nuclearPhases2
        normedDotProducts
        normedDotProducts2
        
    end
    
    methods
        
        function flex = flowExperiment(options)
            
            flex.dataDir = options.dataDir ;
            flex.dataSetName = options.dataSetName ;
            flex.lastTimePoint = options.lastTimePoint ;
            flex.fileSize = options.fileSize ;
            flex.trackFileName = options.trackFileName ;
            flex.resolution = options.res ;
            flex.lobeSplit = options.lobeSplit ;
            flex.dt = options.dt ;
            
            flex.subSamplingFactor = 2 ;
            flex.rawDataFileName = fullfile(flex.dataDir, [flex.dataSetName '_Ch%01d_T%02d']) ;
            flex.probabilitiesFileName = fullfile(flex.dataDir, [ flex.dataSetName '_Ch1_T%02d_Probabilities']) ;
            flex.maskedDataFileName = fullfile(flex.dataDir, [ flex.dataSetName '_Ch%01d_T%02d_masked']) ;
            flex.maskedMIPFileName = fullfile(flex.dataDir, 'Ch%01d_MIPs', [flex.dataSetName '_Ch%01d_T%02d_maskedMIP']) ;
            flex.ch1FileName = fullfile(flex.dataDir, 'Ch1_MIPs', [flex.dataSetName '_Ch1_T%02d_maskedMIP']) ;
            flex.ch2FileName = fullfile(flex.dataDir, 'Ch2_MIPs', [flex.dataSetName '_Ch2_T%02d_maskedMIP']) ;
            flex.ch1FlowFileName = fullfile(flex.dataDir, 'PIVlab_membrane.mat') ;
            flex.ch2FlowFileName = fullfile(flex.dataDir, 'PIVlab_nuclei.mat') ;
            
        end
        
        function makeMaskedMIPs(flex, zbounds)
            
            % Take raw data and probability fields from training on that
            % raw data and outputting masked MIPs
            
            cd(flex.dataDir)
            
            for cc = 1:2
                clear maskData
                surfaceArray = zeros(flex.lastTimePoint, flex.fileSize(1), flex.fileSize(2)) ;
                for tt = 1:flex.lastTimePoint % iterate over timepoints
                    
                    rawData = h5read([sprintf(flex.rawDataFileName , cc, tt), '.h5'], '/inputData') ; % importdata raw data of this particular timepoint
                    probField = h5read([sprintf(flex.probabilitiesFileName, tt), '.h5'], '/exported_data') ; % outputed prob from Ilastik corresponding to this timepoint
                    probField = squeeze(probField(1,:,:,:) ) ;
                    
                    % Mask the raw data and projects onto 2D
                    [maskedData, medianValues] = projectMaskedSurface(double(rawData), probField, 2) ;
                    
                    % maskedData is 3D TIFF * binarized/thresholded training
                    saveastiff(maskedData, [sprintf(flex.maskedDataFileName, cc, tt), '.tiff'])  ;
                    disp(['Saving ', sprintf([sprintf(flex.maskedDataFileName, cc, tt), '.tiff'])]) ;
                    
                    surfaceArray(tt, :, :) = medianValues ;
                    mask = round(squeeze(surfaceArray(tt, :, :))) ; % extract z(x,y) position from surfaceArray
                    
                    [xx,yy] = meshgrid(1:size(maskedData, 2), 1:size(maskedData, 1));
                    mask3d = false(size(maskedData)) ; % initiate a boolean 3D array that will mask the rawData
                    
                    xmax = size(maskedData, 1) ;
                    ymax = size(maskedData, 2) ;
                    zmax = size(maskedData, 3) ;
                    
                    % Set the bounds of the fucial volume in Z
                    if nargin < 2 
                        lowbnd = 0 ;
                        upbnd = zmax ;
                    else
                        lowbnd = zbounds(1) ;
                        upbnd = zbounds(2) ;
                    end
                    
                    for chadd = lowbnd:upbnd
                        % what page to set to true?
                        zind = mask + chadd ;
                        % convert to linear index
                        lin = (zind - 1) * xmax * ymax + (xx-1)*ymax + yy;
                        % clip linear indices to be inside volume
                        lin = lin(lin > 0) ;
                        lin = lin(lin < xmax * ymax * zmax) ;
                        % Set this layer to be true
                        mask3d(lin) = true ;
                    end
                    
                    maskedData= uint8(maskedData).*uint8(mask3d) ; % mask the raw data using our boolean 3d mask
                    
                    maskedData = imadjust(uint8(max(maskedData, [], 3))) ; % take the max value along the z-dim and adjust to new map
                    
                    if isfolder(fullfile(flex.dataDir, ['Ch',num2str(cc),'_MIPs']))
                        mkdir(fullfile(flex.dataDir, ['Ch',num2str(cc),'_MIPs'])) ;
                    end
                    
                    % Save this masked timepoint channel as tiff
                    saveastiff(maskedData, [sprintf(flex.maskedMIPFileName, cc, cc, tt), '.tiff'])  ;
                    close all
                end
                if cc == 1
                    save(fullfile(flex.dataDir, 'surfaceArray.m'), 'surfaceArray') ;
                end
            end
            
        end
        
        function smoothMIPs(flex)
            
            %Here we will importdata images, smooth and threshold each 
            %timepoint in projected space, then write each timepoint to
            %disk for analysis in PIVlab
            
            cd(flex.dataDir)
            
            % Channel 1 smoothing in space
            for tt = 1:flex.lastTimePoint
                if isfile([sprintf(flex.ch1FileName,tt) '_smoothed.tiff'])
                    im = loadtiff([sprintf(flex.ch1FileName,tt) '.tiff']) ;
                    ims = medfilt2(im, [2 2]) ;
                    saveastiff(ims, [sprintf(flex.ch1FileName,tt) '_smoothed.tiff']) ;
                end
            end
            
            % Channel 2 smoothing in space
            for tt = 1:flex.lastTimePoint
                if isfile([sprintf(flex.ch2FileName,tt) '_smoothed.tiff'])
                    im = loadtiff([sprintf(flex.ch2FileName,tt) '.tiff']) ;
                    ims = medfilt2(im, [2 2]) ;
                    saveastiff(ims, [sprintf(flex.ch2FileName,tt) '_smoothed.tiff']) ;
                end
            end
            
        end
        
        function flex = loadTrackData(flex)
            
            % read in the track data
            flex.trackData = uint8(squeeze(h5read( ...
                fullfile(flex.dataDir, flex.trackFileName), ...
                '/exported_data'))) ; 
            trackIdents = unique(flex.trackData) ; % all the grayscale values of tracked nuclei
            flex.trackVals = trackIdents(2:end) ;
            
        end
        
        function flex = buildFullFlow(flex)
            
            % Outputs
            % -------
            % fullflow : 4 x N x M float array
            %   fullflow(1, i, j) is the u component of the flow in channel
            %   1 at full resolution TIFF pixel (i, j) in (x,y) coords in
            %   um per dt 
            %   fullflow(2, :, :) is the v component of ch 1
            %   fullflow(3, :, :) is the u component of ch 2
            %   fullflow(4, :, :) is the v component of ch 2
            %   
            
            cd(flex.dataDir)
            
            % Preallocate the array that stores the full res flow data
            fullflow = zeros(4, flex.fileSize(1), flex.fileSize(2), flex.lastTimePoint-1) ;
            
            % importdata flow structs for the diff channels
            flowMembrane = importdata(flex.ch1FlowFileName) ;
            
            for tt = 1:(flex.lastTimePoint-1)
                xcoord1 = flowMembrane.x{tt} ;
                ycoord1 = flowMembrane.y{tt} ;
                error('todo: check what happens if you use u_filtered')
                ucomp1= flowMembrane.u_original{tt} ;
                vcomp1 = flowMembrane.v_original{tt} ;
                % perform median filtering on each slice
                ucomp1 = medfilt2(inpaintn(ucomp1)) ;
                vcomp1 = medfilt2(inpaintn(vcomp1)) ;
                for xx = 1:size(xcoord1, 2)
                    for yy = 1:size(ycoord1, 1)
                        xval = xcoord1(1, xx) ;
                        yval = ycoord1(yy, 1) ;
                        fullflow(1, yval, xval, tt) = flex.resolution*ucomp1(yy,xx) ;
                        fullflow(2, yval, xval, tt) = flex.resolution*vcomp1(yy,xx) ;
                    end
                end
            end
            
            % importdata flow struct
            flowNuclei = importdata(flex.ch2FlowFileName) ;
            
            for tt = 1:(flex.lastTimePoint-1)
                xcoord2 = flowNuclei.x{tt} ;
                ycoord2 = flowNuclei.y{tt} ;
                error('todo: check what happens if you use u_filtered')
                ucomp2 = flowNuclei.u_original{tt} ;
                vcomp2 = flowNuclei.v_original{tt} ;
                for xx = 1:size(xcoord1, 2)
                    for yy = 1:size(ycoord1, 1)
                        xval = xcoord2(1, xx) ;
                        yval = ycoord2(yy, 1) ;
                        fullflow(3, yval, xval, tt) = flex.resolution*ucomp2(yy,xx) ;
                        fullflow(4, yval, xval, tt) = flex.resolution*vcomp2(yy,xx) ;
                    end
                end
            end
            
            % filter the data in time using a normed triangle filter
            conv = [1 2 3 4 5 4 3 2 1] ; % the vector used for convolution
            hh = conv / sum(conv) ;
            hh = reshape(hh, [1,1,1,length(conv)]) ;
            fullflow = imfilter(fullflow, hh, 'replicate') ;
            
            flex.fullFlowField = fullflow ;
            
            % Save to disk
            outputFn = fullfile(flex.dataDir, 'fullflow.m') ;
            if ~isfile(outputFn)
                save(outputFn, 'fullflow')
            end
            
        end
        
        function flex = maskBinFlow(flex)
            % Mask both channels' flow fields using the probability fields 
            % for nuclear training.
            
            
            % Centroids of each nucleus, and the median velocity within
            % each nuclear region
            centers = zeros(flex.lastTimePoint-1, length(flex.trackVals), 2) ; % centers of nuclei (tt, nuclei, coord)
            binvels = zeros(4, flex.lastTimePoint-1, length(flex.trackVals)) ; % binned velocities (channel/comp, tt, nuclei)
            
            centers2 = zeros(2, flex.lastTimePoint-1, length(flex.trackVals), 2) ; % centers of nuclei (tt, nuclei, coord)
            binvels2 = zeros(2, 4, flex.lastTimePoint-1, length(flex.trackVals)) ; % binned velocities (channel/comp, tt, nuclei)
            
            normdots = zeros(flex.lastTimePoint, length(flex.trackVals)) ;
            normdots2 = zeros(2, flex.lastTimePoint, length(flex.trackVals)) ;
            
            % Load the flow field for both channels and all time
            if isempty(flex.fullFlowField)
                flex.fullFlowField = importdata([flex.dataDir 'fullflow.m']) ;
            end
            
            % Consider each timepoint (last axis)
            for tt = 1:(flex.lastTimePoint-1)
                slice = squeeze(flex.fullFlowField(:,:,:,tt)) ;
                
                % Consider each nucleus/track
                for ni = 1:length(flex.trackVals)
                    % obtain segmentation value (color) for this track
                    trackVal = flex.trackVals(ni) ;
                    
                    % Find all entries of trackData for this nucleus/track
                    % NOTE: sometimes tt was first entry, sometimes last,
                    % so check for this. Assume that time is short relative
                    % to spatial dimension
                    if size(flex.trackData, 1) < size(flex.trackData, 3)
                        mask = squeeze(flex.trackData(tt,:,:)) == trackVal ;
                    elseif size(flex.trackData, 1) > size(flex.trackData, 3)
                        mask = squeeze(flex.trackData(:,:,tt)) == trackVal ;
                    end
                    
                    if mean(mask(:)) > 0 % if the nuclei exists on this slice, find the center
                        s = regionprops(mask, 'centroid') ;
                        centers(tt, ni, :) = s.Centroid ;
                        
                        for ii = 1:4 % median the masked flow data and store it
                            slicem = squeeze(slice(ii, :, :)).*mask ;
                            slicem(slicem == 0) = NaN ;
                            if ~isnan(nanmedian(slicem(:)))
                                binvels(ii, tt, ni) = nanmedian(slicem(:)) ;
                                [~, col] = ind2sub(flex.fileSize, median(find(mask))) ;
                                
                                % Split the left and right lobes (A,P)
                                if col < flex.lobeSplit
                                    % anterior lobe
                                    centers2(1, tt, ni, :) = s.Centroid ;
                                    binvels2(1, ii, tt, ni) = nanmedian(slicem(:)) ;
                                elseif col > flex.lobeSplit
                                    % posterior lobe
                                    centers2(2, tt, ni, :) = s.Centroid ;
                                    binvels2(2, ii, tt, ni) = nanmedian(slicem(:)) ;
                                end
                            end
                        end
                        
                    end
                end
            end
            
            % Pixels which were not tracked are set to NaN
            binvels(binvels == 0) = NaN;
            binvels2(binvels2 == 0) = NaN;
            
            % Calculate relative motion by doing (u2-u1) and (v2-v1)
            % relativeVelocities2 is for left/right lobes separately
            flex.relativeVelocities = zeros(2, flex.lastTimePoint-1, length(flex.trackVals)) ;
            flex.relativeVelocities2 = zeros(2, 2, flex.lastTimePoint-1, length(flex.trackVals)) ;
            
            % Subtract the nuclei-binned channel velocities
            flex.relativeVelocities(1, :, :) = binvels(3, :, :) - binvels(1, :, :) ; % u_rel = u2-u1
            flex.relativeVelocities(2, :, :) = binvels(4, :, :) - binvels(2, :, :) ; % v_rel = v2 - v1
            flex.relativeVelocities2(:, 1, :, :) = binvels2(:, 3, :, :) - binvels2(:, 1, :, :) ; % u_rel = u2-u1
            flex.relativeVelocities2(:, 2, :, :) = binvels2(:, 4, :, :) - binvels2(:, 2, :, :) ; % v_rel = v2 - v1
            
            for ni = 1:length(flex.trackVals)
                for tt = 1:(flex.lastTimePoint-1)
                    A = [binvels(1, tt, ni) binvels(2, tt, ni)] ; % membrane vector
                    B = [binvels(3, tt, ni) binvels(4, tt, ni)] ; % nuclei vector
                    normdots(tt, ni) = dot(A, B)/norm(A)/norm(B) ;
                    for ll = 1:2
                        A = [binvels2(ll, 1, tt, ni) binvels2(ll, 2, tt, ni)] ; % membrane vector
                        B = [binvels2(ll, 3, tt, ni) binvels2(ll, 4, tt, ni)] ; % nuclei vector
                        normdots2(ll, tt, ni) = dot(A, B)/norm(A)/norm(B) ;
                    end
                end
            end
            
            % Output to disk binned velocities, relative velocities, dot
            % products, and nuclear center positions
            flex.binnedVelocities = binvels ;
            flex.binnedVelocities2 = binvels2 ;
            save(fullfile(flex.dataDir, 'binvels.m'), flex.binnedVelocities) ;
            save(fullfile(flex.dataDir, 'binvels2.m'), flex.binnedVelocities2) ;
            
            flex.relativeVelocities = relvels ;
            flex.relativeVelocities2 = relvels2 ;
            save(fullfile(flex.dataDir, 'relvels.m'), flex.relativeVelocities) ;
            save(fullfile(flex.dataDir, 'relvels2.m'), flex.relativeVelocities2) ;
            
            save(fullfile(flex.dataDir, 'normdots.m'), flex.normedDotProducts) ;
            save(fullfile(flex.dataDir, 'normdots2.m'),  flex.normedDotProducts2) ;
            flex.normedDotProducts = normdots ;
            flex.normedDotProducts2 = normdots2 ;
            
            flex.nuclearCenters = centers ;
            flex.nuclearCenters2 = centers2 ;
            save(fullfile(flex.dataDir, 'centers.m'), flex.nuclearCenters) ;
            save(fullfile(flex.dataDir, 'centers2.m'), flex.nuclearCenters2) ;
            
        end
        
        function flex = phaseMagVisualization(flex)
            %
            % Outputs
            % -------
            % phases : #timepoints x #tracks/nuclei float array
            % 	relative phases between velocities averaged over 
            %   each nucleus 
            %   >> saved to fullfile(flex.dataDir, 'phases.m')
            % phases2 : 2 x #timepoints x #tracks/nuclei float array
            %   each lobe's relative phase between velocities median 
            %   over each nucleus 
            %   >> saved to fullfile(flex.dataDir, 'phases2.m')
            
            if isempty(flex.fullFlowField)
                flex.fullFlowField = importdata(fullfile(flex.dataDir, 'fullflow.m')) ;
            end
            if isempty( flex.binnedVelocities)
                flex.binnedVelocities = importdata(fullfile(flex.dataDir, 'binvels.m')) ;
            end
            if isempty( flex.binnedVelocities2)
                flex.binnedVelocities2 = importdata(fullfile(flex.dataDir, 'binvels2.m')) ;
            end
            if isempty( flex.relativeVelocities)
                flex.relativeVelocities = importdata(fullfile(flex.dataDir, 'relvels.m')) ;
            end
            if isempty( flex.relativeVelocities2)
                flex.relativeVelocities2 = importdata(fullfile(flex.dataDir, 'relvels2.m')) ;
            end
            if isempty( flex.nuclearCenters)
                flex.nuclearCenters = importdata(fullfile(flex.dataDir, 'centers.m')) ;
            end
            if isempty( flex.nuclearCenters2)
                flex.nuclearCenters2 = importdata(fullfile(flex.dataDir, 'centers2.m')) ;
            end
            if isempty(flex.nuclearPhases)
                flex.nuclearPhases = importdata(fullfile(flex.dataDir, 'phases.m')) ;
            end
            if isempty( flex.nuclearPhases2)
                flex.nuclearPhases2 = importdata(fullfile(flex.dataDir, 'phases2.m')) ;
            end
            if isempty( flex.normedDotProducts)
                flex.normedDotProducts = importdata(fullfile(flex.dataDir, 'normdots.m')) ;
            end
            if isempty( flex.normedDotProducts2)
                flex.normedDotProducts2 = importdata(fullfile(flex.dataDir, 'normdots2.m')) ;
            end
            
            % Pre-allocate the magnitudes, phases, relative velocities
            colors = flipud(single(hsv))  ; % define color map
            phases = zeros(flex.lastTimePoint-1, length(flex.trackVals)) ;
            mags = zeros(flex.lastTimePoint-1, length(flex.trackVals)) ;
            phases2 = zeros(2, flex.lastTimePoint-1, length(flex.trackVals)) ;
            mags2 = zeros(2, flex.lastTimePoint-1, length(flex.trackVals)) ;
            
            relvelx = squeeze(flex.relativeVelocities(1,:,:)) ;
            relvely = squeeze(flex.relativeVelocities(2,:,:)) ;
            maxMag = sqrt(abs(max(relvelx(:)))^2 + abs(max(relvely(:)))^2) ; % calculate largest relative velocity magnitude
            
            for tt = 1:(flex.lastTimePoint-1)
                close all
                figure('units','normalized','outerposition',[0 0 1 1]) ;
                title(sprintf('Relative motion, t=%02d min', tt))
                ch1 = imadjust(uint8(loadtiff([sprintf(flex.ch1FileName,tt), '.tiff']))) ;  % read in membrane image for this TP
                
                % channel 1 MIP
                memRGB = cat(3, ch1, ch1, ch1) ;
                
                % add a red scale bar 5 um in length
                error('add scalebar here')
                memRGB((500:506), 2*(440:(440+round(5/flex.resolution))) , 1) = 255 ;
                memRGB((500:506), 2*(440:(440+round(5/flex.resolution))) , 2) = 0 ;
                memRGB((500:506), 2*(440:(440+round(5/flex.resolution))) , 3) = 0 ;
                
                imshow(memRGB) ;
                hold on
                xlim([0 flex.fileSize(1)]) ;
                ylim([0 flex.fileSize(2)]) ;
                cidx = squeeze(flex.nuclearCenters(tt,:,:)) ; % centers of nuclei for this TP
                % For each nucleus, find slice in this timepoint, if any
                for ni = 1:length(flex.trackVals)
                    
                    if size(flex.trackData,1) < size(flex.trackData, 3)
                        trackSlice = squeeze(flex.trackData(tt,:,:)) ; % the current tp from trackData
                    else
                        trackSlice = squeeze(flex.trackData(:,:,tt)) ; % the current tp from trackData
                    end
                    
                    % Create mask of this track in this timepoint
                    trackVal = flex.trackVals(ni) ; % the int label of the current nucleus
                    trackSlice(trackSlice ~= trackVal) = 0 ; % set all values not equal to the nuclei's label to zero
                    % transpose due to ilastik Output properties
                    trackSlice = trackSlice' ;
                    
                    % Check that this track is present in the current TP
                    if length(unique(trackSlice)) > 1
                        bound = bwboundaries(trackSlice) ; % grab the boundary of this nuclei
                        bound = bound{1} ;
                        if ~isnan(relvelx(tt, ni)) % if this relative velocity exists
                            if ~isnan(relvely(tt, ni))
                                % Magnitude of relative velocity (norm)
                                currMag = sqrt(flex.relativeVelocities(1, tt, ni)^2 + ...
                                    flex.relativeVelocities(2, tt, ni)^2) ; % calculate magnitude of current nuclei's relative velocity
                                if isnan(currMag)
                                    currMag = 0 ;
                                end
                                % divide the current relative velocity 
                                % magnitude by the maximum relative velocity magnitude
                                % (this will be a float between 0 and 1 and
                                % will later serve as the alpha for this nucleus)
                                scaledMag = currMag/maxMag ; 
                                % find the phase of this nuclei's relative velocity
                                phase = angle(flex.binnedVelocities(1, tt,ni) +...
                                    1i*flex.binnedVelocities(2, tt,ni)) - ...
                                    angle(flex.binnedVelocities(3, tt,ni) + ...
                                    1i*flex.binnedVelocities(4, tt,ni)) ; 
                                phase = mod(phase, 2*pi) ;
                                if isnan(phase)
                                    error('phase should be real number')
                                    phase = 0 ;
                                end
                                phases(tt,ni) = phase ;
                                mags(tt, ni) = scaledMag ;
                                
                                % Convert the phase to a colormap row
                                phase2row = abs(round((length(colors)/2/pi)*(phase-pi))) ; % this is the function that maps phase to a row in the colormap
                                if phase2row > 256
                                    phase2row = 256 ;
                                elseif phase2row < 1
                                    phase2row = 1 ;
                                end
                                if scaledMag > 1
                                    scaledMag = 1 ;
                                end
                                % grab the color for this nucleus using its 
                                % phase conversion
                                color = colors(phase2row, :) ; 
                                % Add the nucleus patch using boundary and
                                % color
                                patch(bound(:,1), bound(:,2), color, 'FaceAlpha', scaledMag, 'EdgeColor', 'none') ;
                            else
                                error('should not be the case: missing relvely')
                            end
                        else
                            error('should not be the case: missing relvelx')
                        end
                    end
                end
                
                % Add a legend
                hold on
                phasebar('colormap', colors, 'rad', 'location', 'se', 'size', 0.2) ;
                
                % Add velocity arrows to the axis
                mem = quiver(cidx(:,1), cidx(:,2), squeeze(flex.binnedVelocities(1, tt, :)), squeeze(flex.binnedVelocities(2, tt,:)), 1, 'LineWidth', 1, 'color', 'green') ;
                nuc = quiver(cidx(:,1), cidx(:,2), squeeze(flex.binnedVelocities(3, tt,:)), squeeze(flex.binnedVelocities(4, tt,:)), 1, 'LineWidth', 1, 'color', 'red') ;
                legend([mem, nuc] , 'Muscle', 'Endoderm', 'location', 'northeastoutside') ;
                
                % todo: scalebar
                text(500, 480 , '5 um', 'FontSize', 18, 'Color', 'black') ;
                
                if ~isfolder(fullfile(flex.dataDir, 'phasemag_visualization'))
                    mkdir(fullfile(flex.dataDir, 'phasemag_visualization')) ;
                end
                
                % Save to disk
                saveas(gcf, fullfile(flex.dataDir, ...
                    'phasemag_visualization', ...
                    [sprintf('phasemag_T%02d', tt), '.png'])) ;
            end
            close all
            
            % Same as above except now for individual lobes
            for tt = 1:(flex.lastTimePoint-1)
                maxMag2 = zeros(2) ;
                for ll = 1:2
                    relvelx = squeeze(flex.relativeVelocities2(ll, 1,:,:)) ;
                    relvely = squeeze(flex.relativeVelocities2(ll, 2,:,:)) ;
                    maxMag2(ll) = sqrt(abs(max(relvelx(:)))^2 + abs(max(relvely(:)))^2) ; % calculate largest relative velocity magnitude
                    for ni = 1:length(flex.trackVals)
                        if size(flex.trackData,1) < size(flex.trackData, 3)
                            trackSlice = squeeze(flex.trackData(tt,:,:)) ; % the current tp from trackData
                        else
                            trackSlice = squeeze(flex.trackData(:,:,tt)) ; % the current tp from trackData
                        end
                        trackVal = flex.trackVals(ni) ; % the int label of the current nucleus
                        trackSlice(trackSlice ~= trackVal) = 0 ; % set all values not equal to the nuclei's label to zero
                        trackSlice = trackSlice' ;
                        if length(unique(trackSlice)) > 1
                            if ~isnan(relvelx(tt, ni)) % if this relative velocity exists
                                if ~isnan(relvely(tt, ni))
                                    currMag2 = sqrt(flex.relativeVelocities2(ll, 1, tt, ni)^2 + flex.relativeVelocities2(ll, 2, tt, ni)^2) ; % calculate magnitude of current nuclei's relative velocity
                                    scaledMag2 = currMag2/maxMag2(ll) ; % divide the current relative velocity magnitude by the maximum relative velocity magnitude
                                    phase = angle(flex.binnedVelocities2(ll, 1, tt,ni) + 1i*flex.binnedVelocities2(ll, 2, tt,ni))-angle(flex.binnedVelocities2(ll, 3, tt,ni) + 1i*flex.binnedVelocities2(ll, 4, tt,ni)) ; % find the phase of this nuclei's relative velocity
                                    phase = mod(phase, 2*pi) ;
                                    phases2(ll, tt,ni) = phase ;
                                    mags2(ll, tt, ni) = scaledMag2 ;
                                end
                            end
                        end
                    end
                end
            end
            
            save(fullfile(flex.dataDir, 'phases.m'), 'phases') ;
            flex.nuclearPhases = phases ;
            
            save(fullfile(flex.dataDir, 'phases2.m'), 'phases2') ;
            flex.nuclearPhases = phases2 ;
            
        end
        
        function createMasterPlot(flex)
            
            % The master plot will include a 95% confidence interval
            % ellipse. Code for that was taken from Vision Dummy.
            
            if isempty(flex.fullFlowField)
                flex.fullFlowField = importdata(fullfile(flex.dataDir, 'fullflow.m')) ;
            end
            if isempty( flex.binnedVelocities)
                flex.binnedVelocities = importdata(fullfile(flex.dataDir, 'binvels.m')) ;
            end
            if isempty( flex.binnedVelocities2)
                flex.binnedVelocities2 = importdata(fullfile(flex.dataDir, 'binvels2.m')) ;
            end
            if isempty( flex.relativeVelocities)
                flex.relativeVelocities = importdata(fullfile(flex.dataDir, 'relvels.m')) ;
            end
            if isempty( flex.relativeVelocities2)
                flex.relativeVelocities2 = importdata(fullfile(flex.dataDir, 'relvels2.m')) ;
            end
            if isempty( flex.nuclearCenters)
                flex.nuclearCenters = importdata(fullfile(flex.dataDir, 'centers.m')) ;
            end
            if isempty( flex.nuclearCenters2)
                flex.nuclearCenters2 = importdata(fullfile(flex.dataDir, 'centers2.m')) ;
            end
            if isempty(flex.nuclearPhases)
                flex.nuclearPhases = importdata(fullfile(flex.dataDir, 'phases.m')) ;
            end
            if isempty( flex.nuclearPhases2)
                flex.nuclearPhases2 = importdata(fullfile(flex.dataDir, 'phases2.m')) ;
            end
            if isempty( flex.normedDotProducts)
                flex.normedDotProducts = importdata(fullfile(flex.dataDir, 'normdots.m')) ;
            end
            if isempty( flex.normedDotProducts2)
                flex.normedDotProducts2 = importdata(fullfile(flex.dataDir, 'normdots2.m')) ;
            end
            
            % Prepare figure
            close all
            figure('units','normalized','outerposition',[0 0 1 1])
            lin = linspace(-1, 1, 3) ;
            cmap = parula ;
            
            % Map phases from (0, 2pi) to (-pi, pi)
            for ii = 1:length(flex.nuclearPhases(:))
                if flex.nuclearPhases(ii) > pi
                    flex.nuclearPhases(ii) = flex.nuclearPhases(ii) - 2*pi ;
                end
            end
            
            colors = zeros(flex.lastTimePoint-1, 3) ;
            
            for tt = 1:(flex.lastTimePoint-1)
                
                color = cmap(round(length(parula)/flex.lastTimePoint*tt), :) ;
                colors(tt,:) = color ;
                
                subplot(2, 3, 1)
                title('u-comp velocities') ;
                ylabel('V_u for Membrane data [um/min]') ;
                xlabel('V_u for Nuclear data [um/min]') ;
                hold on
                uu1 = squeeze(flex.binnedVelocities(1, tt, :)) ;
                uu2 = squeeze(flex.binnedVelocities(3, tt,:)) ;
                uu1 = filloutliers(uu1, 'nearest', 'percentiles', [15 85]) ;
                uu2 = filloutliers(uu2, 'nearest', 'percentiles', [15 85]) ;
                plot(lin, lin, '--')
                scatter(uu2, uu1, 5, color, 'filled', 'markeredgecolor', 'none') ;
                xlim([-.5 .5])
                ylim([-.5 .5])
                axis equal
                
                subplot(2, 3, 2)
                title('v-comp velocities') ;
                ylabel('V_v for Membrane data [um/min]') ;
                xlabel('V_v for Nuclear data [um/min]') ;
                hold on
                vv1 = squeeze(flex.binnedVelocities(2, tt, :)) ;
                vv2 = squeeze(flex.binnedVelocities(4, tt,:)) ;
                vv1 = filloutliers(vv1, 'nearest', 'percentiles', [15 85]);
                vv2 = filloutliers(vv2, 'nearest', 'percentiles', [15 85]) ;
                plot(lin, lin, '--')
                scatter(vv2, vv1, 5, color, 'filled', 'markeredgecolor', 'none') ;
                xlim([-.5 .5])
                ylim([-.5 .5])
                axis equal
                
                subplot(2, 3, 3)
                title('Relative motion') ;
                ylabel('V_v_,_r_e_l') ;
                xlabel('V_u_,_r_e_l') ;
                hold on
                scatter(squeeze(flex.relativeVelocities(1, tt, :)), squeeze(flex.relativeVelocities(2, tt, :)), 5, color, 'filled', 'markeredgecolor', 'none') ;
                xlim([-.5 .5])
                ylim([-.5 .5])
                
            end
            
            subplot(2, 3, 4)
            N2 = zeros(flex.lastTimePoint-1, 9) ;
            for tt = 1:flex.lastTimePoint-1
                title('Relative Phase') ;
                ylabel('Count') ;
                xlabel('Phase') ;
                xlim([-pi pi]);
                hold on
                phase = flex.nuclearPhases(tt, flex.nuclearPhases(tt, :) ~= 0) ;
                [~,edges] = histcounts(phase, linspace(-pi, pi, 10));
                N2(tt, :) = histcounts(phase ,edges); % Bin using the same edges
            end
            
            ctrs = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
            b = bar(ctrs, N2, 'stacked', 'edgecolor', 'none') ;
            for k = 1:size(flex.nuclearPhases,1)
                b(k).FaceColor = colors(k,:) ;
            end
            
            subplot(2, 3, 5)
            relmags = zeros(flex.lastTimePoint-1, length(flex.trackVals)) ;
            for tt = 1:flex.lastTimePoint-1
                title('Relative Magnitude |V_m_e_m|/|V_n_u_c|') ;
                ylabel('Count') ;
                xlabel('Relative Magnitude') ;
                relmags(tt, :) = sqrt((squeeze(flex.binnedVelocities(1, tt, :))).^2 + (squeeze(flex.binnedVelocities(2, tt, :))).^2)./sqrt((squeeze(flex.binnedVelocities(3, tt, :))).^2 + (squeeze(flex.binnedVelocities(4, tt, :))).^2) ;
                hold on
                relmag = relmags(tt, relmags(tt,:) ~=0) ;
                [~,edges] = histcounts(relmag, linspace(0, 2, 10));
                N2(tt, :) = histcounts(relmag ,edges); % Bin using the same edges
            end
            ctrs = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
            b = bar(ctrs, N2, 'stacked', 'FaceColor', 'flat', 'edgecolor', 'none') ;
            for k = 1:size(flex.nuclearPhases,1)
                b(k).CData = colors(k,:) ;
            end
            
            subplot(2, 3, 6)
            for tt = 1:flex.lastTimePoint-1
                title('Normed Dot Products') ;
                ylabel('Count') ;
                xlabel('Normed Dot Product Values');
                hold on
                normdot = flex.normedDotProducts(tt, flex.normedDotProducts(tt,:) ~=0) ;
                [~,edges] = histcounts(normdot, linspace(-1, 1, 10));
                N2(tt, :) = histcounts(normdot ,edges); % Bin using the same edges
            end
            ctrs = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
            b = bar(ctrs, N2, 'stacked', 'FaceColor', 'flat', 'edgecolor', 'none') ;
            for k = 1:size(flex.nuclearPhases,1)
                b(k).CData = colors(k,:) ;
            end
            
            if ~isfolder(fullfile(flex.dataDir, 'correlation plots'))
                mkdir(fullfile(flex.dataDir, 'correlation plots')) ;
            end
            
            relvelucomp = squeeze(flex.relativeVelocities(1,:,:)) ;
            relvelvcomp = squeeze(flex.relativeVelocities(2,:,:)) ;
            data = horzcat(relvelucomp(:),relvelvcomp(:)) ;
            [r_ellipse, X0, Y0] = error_ellipse(data) ;
            subplot(2, 3, 3)
            p = plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'color', 'r','linestyle', '-') ;
            legend(p, '95% Confidence Interval') ;
            sgtitle(['Correlation Plot for ' flex.dataSetName]) ;
            saveas(gcf, fullfile(flex.dataDir, 'correlation plots', 'master_corr.png')) ;
            
            clear relmags
            
            % Two-lobe version of same plot
            for ll = 1:2
                
                close all
                figure('units','normalized','outerposition',[0 0 1 1])
                lin = linspace(-1, 1, 3) ;
                
                cmap = parula ;
                
                lobePhases = squeeze(flex.nuclearPhases2(ll, :, :)) ;
                for ii = 1:length(lobePhases(:))
                    if lobePhases(ii) > pi
                        lobePhases(ii) = lobePhases(ii) - 2*pi ;
                    end
                end
                
                for tt = 1:(flex.lastTimePoint-1)
                    
                    color = cmap(round(length(parula)/flex.lastTimePoint*tt), :) ;
                    colors(tt,:) = color ;
                    
                    subplot(2, 3, 1)
                    title(['u-comp velocities lobe ' num2str(ll)]) ;
                    ylabel('V_u for Membrane data [um/min]') ;
                    xlabel('V_u for Nuclear data [um/min]') ;
                    hold on
                    uu1 = squeeze(flex.binnedVelocities2(ll, 1, tt, :)) ;
                    uu2 = squeeze(flex.binnedVelocities2(ll, 3, tt,:)) ;
                    uu1 = filloutliers(uu1, 'nearest', 'percentiles', [15 85]) ;
                    uu2 = filloutliers(uu2, 'nearest', 'percentiles', [15 85]) ;
                    plot(lin, lin, '--')
                    scatter(uu2, uu1, 5, color, 'filled', 'markeredgecolor', 'none') ;
                    xlim([-.5 .5])
                    ylim([-.5 .5])
                    axis equal
                    
                    subplot(2, 3, 2)
                    title(['v-comp velocities lobe ' num2str(ll)]) ;
                    ylabel('V_v for Membrane data [um/min]') ;
                    xlabel('V_v for Nuclear data [um/min]') ;
                    hold on
                    vv1 = squeeze(flex.binnedVelocities2(ll, 2, tt, :)) ;
                    vv2 = squeeze(flex.binnedVelocities2(ll, 4, tt,:)) ;
                    vv1 = filloutliers(vv1, 'nearest', 'percentiles', [15 85]);
                    vv2 = filloutliers(vv2, 'nearest', 'percentiles', [15 85]) ;
                    plot(lin, lin, '--')
                    scatter(vv2, vv1, 5, color, 'filled', 'markeredgecolor', 'none') ;
                    xlim([-.5 .5])
                    ylim([-.5 .5])
                    axis equal
                    
                    subplot(2, 3, 3)
                    title(['Relative motion lobe ' num2str(ll)]) ;
                    ylabel('V_v_,_r_e_l') ;
                    xlabel('V_u_,_r_e_l') ;
                    hold on
                    scatter(squeeze(flex.relativeVelocities2(ll, 1, tt, :)), squeeze(flex.relativeVelocities2(ll, 2, tt, :)), 5, color, 'filled', 'markeredgecolor', 'none') ;
                    xlim([-.5 .5])
                    ylim([-.5 .5])
                    
                end
                
                subplot(2, 3, 4)
                for tt = 1:flex.lastTimePoint-1
                    title(['Relative Phase ' num2str(ll)]) ;
                    ylabel('Count') ;
                    xlabel('Phase') ;
                    xlim([-pi pi]);
                    hold on
                    phase = lobePhases(tt, lobePhases(tt, :) ~= 0) ;
                    [~,edges] = histcounts(phase, linspace(-pi, pi, 10));
                    N2(tt, :) = histcounts(phase ,edges); % Bin using the same edges
                end
                ctrs = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
                b = bar(ctrs, N2, 'stacked', 'edgecolor', 'none') ;
                for k = 1:size(lobePhases,1)
                    b(k).FaceColor = colors(k,:) ;
                end
                
                subplot(2, 3, 5)
                relmags2 = zeros(2, flex.lastTimePoint-1, length(flex.trackVals)) ;
                for tt = 1:flex.lastTimePoint-1
                    title(['Relative Magnitude |V_m_e_m|/|V_n_u_c| ' num2str(ll)]) ;
                    ylabel('Count') ;
                    xlabel('Relative Magnitude') ;
                    relmags2(ll, tt, :) = sqrt((squeeze(flex.binnedVelocities2(ll, 1, tt, :))).^2 + (squeeze(flex.binnedVelocities2(ll, 2, tt, :))).^2)./sqrt((squeeze(flex.binnedVelocities2(ll, 3, tt, :))).^2 + (squeeze(flex.binnedVelocities2(ll, 4, tt, :))).^2) ;
                    hold on
                    relmag = relmags2(ll, tt, relmags2(ll, tt,:) ~=0) ;
                    [~,edges] = histcounts(relmag, linspace(0, 2, 10));
                    N2(tt, :) = histcounts(relmag ,edges); % Bin using the same edges
                end
                ctrs = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
                b = bar(ctrs, N2, 'stacked', 'FaceColor', 'flat', 'edgecolor', 'none') ;
                for k = 1:size(lobePhases,1)
                    b(k).CData = colors(k,:) ;
                end
                
                subplot(2, 3, 6)
                for tt = 1:flex.lastTimePoint-1
                    title(['Normed Dot Products ' num2str(ll)]) ;
                    ylabel('Count') ;
                    xlabel('Normed Dot Product Values');
                    hold on
                    normdot = flex.normedDotProducts2(ll, tt, flex.normedDotProducts2(ll, tt,:) ~=0) ;
                    [~,edges] = histcounts(normdot, linspace(-1, 1, 10));
                    N2(tt, :) = histcounts(normdot ,edges); % Bin using the same edges
                end
                ctrs = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
                b = bar(ctrs, N2, 'stacked', 'FaceColor', 'flat', 'edgecolor', 'none') ;
                for k = 1:size(lobePhases,1)
                    b(k).CData = colors(k,:) ;
                end
                
                if isfolder(fullfile(flex.dataDir, 'correlation plots'))
                    mkdir(fullfile(flex.dataDir, 'correlation plots')) ;
                end
                
                relvelucomp = squeeze(flex.relativeVelocities2(ll, 1,:,:)) ;
                relvelvcomp = squeeze(flex.relativeVelocities2(ll, 2,:,:)) ;
                data = horzcat(relvelucomp(:),relvelvcomp(:)) ;
                [r_ellipse, X0, Y0] = error_ellipse(data) ;
                subplot(2, 3, 3)
                p = plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'color', 'r','linestyle', '-') ;
                legend(p, '95% Confidence Interval') ;
                            
                sgtitle(['Correlation Plot (lobe ' num2str(ll) ') for ' flex.dataSetName]) ;

                saveas(gcf, fullfile(flex.dataDir, 'correlation plots', ['master_corr_lobe' num2str(ll) '.png'])) ;
            end
        end
    end
    
end
