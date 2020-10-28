%% COMPARING FLOW FIELDS
% Isaac Breinyn 2020
%
% This is a script that loads PIVlab flow data (with an option to load ilastik
% training for it), and compare it to other flow data
%
% NOTE: The PIV data used in this script is formated into arrays of u and v
% component of the velocities for EACH channel. Therefore, arrays with a
% first dim of 4 will be in the order of u1, v1, u2, then v2. The letter
% corresponds to the component of the velocity vector, and the number to
% the prevelant channel.

%% Clean Matlab and Set up Global Parameters

clear; close all; clc ;

addpath_recurse('F:/Streichan/code/gut_matlab/') ;

dataDir = 'F:\Streichan\data\202003111830_mef2gal4klarUASCAAXmChHiFP_wo_63x_e4e5\2038_e5' ; % the directory where the data lives

% Assign global parameters

fn1 = '2038_e5_Ch1_T%02d_MIP' ; % file name of first flow data
fn2 = '2038_e5_Ch2_T%02d_MIP' ; % file name of second flow data

flowfn1 = 'PIVlab_membrane.mat' ;
flowfn2 = 'PIVlab_nuclei.mat' ;

fnp = [fn2, '_Probabilities.h5' ] ; % file name of probabilities that you want to load

fntrack = '/Ch2_MIPs/2038_e5_Ch2_MIPs-myData_Object Identities.h5' ;

ltp = 33 ; % "last timepoint" or the number of TPs in your dataset

res = 0.0833 ; % the data resolution in um/px

conv = [1 2 3 4 5 4 3 2 1] ; % the vector used for convolution

fsize = [512 1024] ;

cd(dataDir) % navigate to the data directory

%% Load the Tracking Data

trackData = uint8(squeeze(h5read(fullfile(dataDir, fntrack), '/exported_data'))) ; % load the track data
trackData = trackData(:,:,1:(ltp-1)) ;
trackVals = unique(trackData) ; % all the grayscale values of tracked nuclei
trackVals = trackVals(2:end) ;

%% Pre-processing for PIVlab
%Here we will load images, smooth and threshold them, then write them back
%for analysis in PIVlab

dataArray = zeros(fsize(1),fsize(2),ltp) ;

for tt = 1:ltp
    if ~exist(fullfile(dataDir,'Ch1_MIPs', [sprintf(fn1,tt) '_smoothed.tiff']))
        im = loadtiff(fullfile(dataDir,'Ch1_MIPs', [sprintf(fn1,tt) '.tiff'])) ;
        ims = medfilt2(im, [2 2]) ;
        %ims(ims<50) = 0;
        saveastiff(ims, fullfile(dataDir,'Ch1_MIPs', [sprintf(fn1,tt) '_smoothed.tiff'])) ;
    end
end

for tt = 1:ltp
    if ~exist(fullfile(dataDir,'Ch2_MIPs', [sprintf(fn2,tt) '_smoothed.tiff']))
        im = loadtiff(fullfile(dataDir,'Ch2_MIPs', [sprintf(fn2,tt) '.tiff'])) ;
        ims = medfilt2(im, [2 2]) ;
        % ims(ims<50) = 0;
        saveastiff(ims, fullfile(dataDir,'Ch2_MIPs', [sprintf(fn2,tt) '_smoothed.tiff'])) ;
    end
end

%% Build and Smooth both Channels

% load flow structs for the diff channels
flowMembrane = load(fullfile(cd, flowfn1)) ;

% Preallocate the 512x512 array that stores the full res flow data
fullflow = zeros(4, fsize(1), fsize(2), ltp-1) ;

for tt = 1:(ltp-1)
    xcoord1 = flowMembrane.x{tt} ;
    ycoord1 = flowMembrane.y{tt} ;
    ucomp1= flowMembrane.u_original{tt} ;
    vcomp1 = flowMembrane.v_original{tt} ;
    % perform median filtering on each slice
    ucomp1 = medfilt2(inpaintn(ucomp1)) ;
    vcomp1 = medfilt2(inpaintn(vcomp1)) ;
    for xx = 1:size(xcoord1, 2)
        for yy = 1:size(ycoord1, 1)
            xval = xcoord1(1, xx) ;
            yval = ycoord1(yy, 1) ;
            fullflow(1, xval, yval, tt) = res*ucomp1(yy,xx) ;
            fullflow(2, xval, yval, tt) = res*vcomp1(yy,xx) ;
        end
    end
end

% load flow struct
flowNuclei = load(fullfile(cd, flowfn2)) ;

for tt = 1:(ltp-1)
    xcoord2 = flowNuclei.x{tt} ;
    ycoord2 = flowNuclei.y{tt} ;
    ucomp2 = flowNuclei.u_original{tt} ;
    vcomp2 = flowNuclei.v_original{tt} ;
    for xx = 1:size(xcoord1, 2)
        for yy = 1:size(ycoord1, 1)
            xval = xcoord2(1, xx) ;
            yval = ycoord2(yy, 1) ;
            fullflow(3, xval, yval, tt) = res*ucomp2(yy,xx) ;
            fullflow(4, xval, yval, tt) = res*vcomp2(yy,xx) ;
        end
    end
end

% filter the data in time using a normed triangle filter
hh = conv / sum(conv) ;
hh = reshape(hh, [1,1,1,length(conv)]) ;
fullflow = imfilter(fullflow, hh, 'replicate') ;

%% Mask and Bin the Flow data for Both Channels

centers = zeros(ltp-1, length(trackVals), 2) ; % centers of nuclei (tt, nuclei, coord)
binvels = zeros(4, ltp-1, length(trackVals)) ; % binned velocities (channel/comp, tt, nuclei)

for tt = 1:(ltp-1)
    slice = squeeze(fullflow(:,:,:,tt)) ;
    for ni = 1:length(trackVals)
        trackVal = trackVals(ni) ;
        mask = squeeze(trackData(:,:,tt)) == trackVal ;
        
        if mean(mask(:)) > 0 % if the nuclei exists on this slice, find the center
            s = regionprops(mask, 'centroid') ;
            centers(tt, ni, :) = s.Centroid ;
            
            for ii = 1:4 % median the masked flow data and store it
                slicem = squeeze(slice(ii, :, :)).*mask ;
                slicem(slicem == 0) = NaN ;
                if ~(nanmedian(slicem(:)) == NaN)
                    binvels(ii, tt, ni) = nanmedian(slicem(:)) ;
                end
            end
        end
    end
end

% Calculate relative motion by doing (u2-u1) and (v2-v1)
relvels = zeros(2, ltp-1, length(trackVals)) ;

for ni = 1:length(trackVals)
    for tt = 1:(ltp-1)
        relvels(1, tt, ni) = binvels(3, tt, ni) - binvels(1, tt, ni) ; % u_rel = u2-u1
        relvels(2, tt, ni) = binvels(4, tt, ni) - binvels(2, tt, ni) ; % v_rel = v2 - v1
    end
end

%% Plot the raw Data with overlayed Binned Velocity Vectors

for tt = 20 %1:(ltp-1)
    close all
    figure('units','normalized','outerposition',[0 0 1 1])
    
    ch1 = imadjust(uint8(loadtiff(fullfile(dataDir, 'Ch1_MIPs', [sprintf(fn1,tt), '.tiff'])))) ;
    ch2 = imadjust(uint8(loadtiff(fullfile(dataDir, 'Ch2_MIPs', [sprintf(fn2,tt), '.tiff'])))) ;
    max = imadjust(uint8(ch1 + ch2)) ;
    max(max > 255) = 255 ;
    im = cat(3, ch2, ch1, max) ;
    
    imshow(im) ;
    title(['Flow Vectors for each Channel TP=' num2str(tt)])  ;
    
    cidx = squeeze(centers(tt,:,:)) ;
    
    hold on
    quiver(cidx(:,1), cidx(:,2), squeeze(binvels(1, tt, :)), squeeze(binvels(2, tt,:)), 1, 'LineWidth', 1, 'color', 'green') ;
    quiver(cidx(:,1), cidx(:,2), squeeze(binvels(3, tt,:)), squeeze(binvels(4, tt,:)), 1, 'LineWidth', 1, 'color', 'red') ;
    legend('Membrane', 'Nuclei', 'location', 'northeastoutside', 'FontSize', 5) ;
    
    saveas(gcf, fullfile(dataDir, 'data_flow_overlays', [sprintf('overlay_T%02d', tt), '.png'])) ;
end
%% Plot the average U and V components over time

correl = zeros(ltp-1, 2) ;

for tt = 1:ltp -1
    uu1 = squeeze(binvels(1, tt, :)) ;
    uu2 = squeeze(binvels(3, tt,:)) ;
    
    vv1 = squeeze(binvels(2, tt, :)) ;
    vv2 = squeeze(binvels(4, tt,:)) ;
    
    corm1 = corrcoef(uu1, uu2, 'Rows', 'pairwise') ;
    correl(tt, 1) = corm1(2) ;
    
    corm2 = corrcoef(vv1, vv2, 'Rows', 'pairwise') ;
    correl(tt, 2) = corm2(2) ;
end

lin = linspace(-1, 1, 3) ;

for tt = 1:(ltp-1) % Create a plot for each TP individually
    
    close all
    figure('units','normalized','outerposition',[0 0 1 1])
    
    subplot(2, 3, 1)
    title(['u-comp velocities TP ' num2str(tt)]) ;
    ylabel('V_u for Membrane data [um/min]') ;
    xlabel('V_u for Nuclear data [um/min]') ;
    hold on
    uu1 = squeeze(binvels(1, tt, :)) ;
    uu2 = squeeze(binvels(3, tt,:)) ;
    uu1 = filloutliers(uu1, 'nearest', 'percentiles', [15 85]) ;
    uu2 = filloutliers(uu2, 'nearest', 'percentiles', [15 85]) ;
    plot(lin, lin, '--')
    s1 = scatter(uu2, uu1, 6, 'blue', 'filled') ;
    xlim([-1 1])
    ylim([-1 1])
    
    subplot(2, 3, 2)
    title(['v-comp velocities TP ' num2str(tt)]) ;
    ylabel('V_v for Membrane data [um/min]') ;
    xlabel('V_v for Nuclear data [um/min]') ;
    hold on
    vv1 = squeeze(binvels(2, tt, :)) ;
    vv2 = squeeze(binvels(4, tt,:)) ;
    vv1 = filloutliers(vv1, 'nearest', 'percentiles', [15 85]);
    vv2 = filloutliers(vv2, 'nearest', 'percentiles', [15 85]) ;
    plot(lin, lin, '--')
    s2 = scatter(vv2, vv1, 6, 'magenta', 'filled') ;
    xlim([-1 1])
    ylim([-1 1])
    
    subplot(2, 3, 3)
    title(['Relative motion for TP ' num2str(tt)]) ;
    ylabel(['V_v,rel for TP ' num2str(tt)]) ;
    xlabel(['V_u,rel for TP ' num2str(tt)]) ;
    hold on
    s3 = scatter(squeeze(relvels(1, tt, :)), squeeze(relvels(2, tt, :)), 6, 'black', 'filled') ;
    xlim([-1 1])
    ylim([-1 1])
    
    subplot(2, 3, [4 6])
    title('Correlation')
    ylabel('Correlation') ;
    xlabel('Time [min]') ;
    hold on
    plot(correl(:,1), 'blue') ;
    plot(correl(:,2), 'magenta') ;
    plot(tt, correl(tt,1), 'bo', 'LineWidth', 2, 'MarkerSize', 16);
    plot(tt, correl(tt,2), 'ro', 'LineWidth', 2, 'MarkerSize', 16);
    legend('Channel 1', 'Channel 2') ;
    ylim([-1 1])
    
    s1.MarkerFaceAlpha = 0.9 ;
    s2.MarkerFaceAlpha = 0.9 ;
    s3.MarkerFaceAlpha = 0.9 ;
    
    saveas(gcf, fullfile(dataDir, 'correlation plots', [sprintf('corr_T%02d', tt), '.png'])) ;
    
end

close all
figure('units','normalized','outerposition',[0 0 1 1])

for tt = 1:(ltp-1)
    
    subplot(2, 3, 1)
    title('u-comp velocities') ;
    ylabel('V_u for Membrane data [um/min]') ;
    xlabel('V_u for Nuclear data [um/min]') ;
    hold on
    uu1 = squeeze(binvels(1, tt, :)) ;
    uu2 = squeeze(binvels(3, tt,:)) ;
    uu1 = filloutliers(uu1, 'nearest', 'percentiles', [15 85]) ;
    uu2 = filloutliers(uu2, 'nearest', 'percentiles', [15 85]) ;
    plot(lin, lin, '--')
    s1 = scatter(uu2, uu1, 6, 'blue', 'filled') ;
    xlim([-1 1])
    ylim([-1 1])
    
    subplot(2, 3, 2)
    title('v-comp velocities') ;
    ylabel('V_v for Membrane data [um/min]') ;
    xlabel('V_v for Nuclear data [um/min]') ;
    hold on
    vv1 = squeeze(binvels(2, tt, :)) ;
    vv2 = squeeze(binvels(4, tt,:)) ;
    vv1 = filloutliers(vv1, 'nearest', 'percentiles', [15 85]);
    vv2 = filloutliers(vv2, 'nearest', 'percentiles', [15 85]) ;
    plot(lin, lin, '--')
    s2 = scatter(vv2, vv1, 6, 'magenta', 'filled') ;
    xlim([-1 1])
    ylim([-1 1])
    
    subplot(2, 3, 3)
    title('Relative motion') ;
    ylabel('V_v,rel') ;
    xlabel('V_u,rel') ;
    hold on
    s3 = scatter(squeeze(relvels(1, tt, :)), squeeze(relvels(2, tt, :)), 6, 'black', 'filled') ;
    xlim([-1 1])
    ylim([-1 1])
    
    s1.MarkerFaceAlpha = 0.25 ;
    s2.MarkerFaceAlpha = 0.25 ;
    s3.MarkerFaceAlpha = 0.25 ;
    
end

subplot(2, 3, [4 6])
title('Correlation')
ylabel('Correlation') ;
xlabel('Time [min]') ;
hold on
plot(correl(:,1), 'blue') ;
plot(correl(:,2), 'magenta') ;
legend('Channel 1', 'Channel 2') ;
ylim([-1 1])

saveas(gcf, fullfile(dataDir, 'correlation plots', 'master_corr.png')) ;

%% Visualize phase and magnitude of velocity differences for each Nuclei

relvelx = squeeze(relvels(1,:,:)) ;
relvely = squeeze(relvels(2,:,:)) ;
maxMag = sqrt(abs(nanmax(relvelx(:)))^2 + abs(nanmax(relvely(:)))^2) ; % calculate largest relative velocity magnitude
colors = hsv ; % define color map
phases = zeros(ltp-1, length(trackVals)) ;
mags = zeros(ltp-1, length(trackVals)) ;
for tt =20 % 1:(ltp-1)
    close all
    figure('units','normalized','outerposition',[0 0 1 1]) ;
    title(sprintf('Relative motion, t=%02d min', tt))
    ch1 = imadjust(uint8(loadtiff(fullfile(dataDir, 'Ch1_MIPs', [sprintf(fn1,tt), '.tiff'])))) ; % load membrane image for this TP
    memRGB = cat(3, ch1, ch1, ch1) ;
    % add a red scale bar 5 um in length
    memRGB((500:506), 2*(440:(440+round(5/res))) , 1) = 255 ;
    memRGB((500:506), 2*(440:(440+round(5/res))) , 2) = 0 ;
    memRGB((500:506), 2*(440:(440+round(5/res))) , 3) = 0 ;
    
    hold on
    imshow(memRGB, 'InitialMagnification', 1000) ;
    xlim([0 fsize(2)])
    ylim([0 fsize(1)])
    cidx = squeeze(centers(tt,:,:)) ; % centers of nuclei for this TP
    for ni = 1:length(trackVals)
        trackSlice = trackData(:,:,tt) ; % the current tp from trackData
        trackVal = trackVals(ni) ; % the int label of the current nucleus
        trackSlice(trackSlice ~= trackVal) = 0 ; % set all values not equal to the nuclei's label to zero
        trackSlice = trackSlice' ;
        if length(unique(trackSlice)) > 1
            bound = bwboundaries(trackSlice) ; % grab the boundary of this nuclei
            bound = bound{1} ;
            if ~isnan(relvelx(tt, ni)) % if this relative velocity exists
                if ~isnan(relvely(tt, ni))
                    currMag = sqrt(relvels(1, tt, ni)^2 + relvels(2, tt, ni)^2) ; % calculate magnitude of current nuclei's relative velocity
                    scaledMag = currMag/maxMag ; % divide the current relative velocity magnitude by the maximum relative velocity magnitude
                    % (this will be a float between 0 and 1 and will later serve as the alpha for this nucleus)
                    phase = angle(relvels(1, tt,ni) + j*relvels(2, tt,ni)) ; % find the phase of this nuclei's relative velocity
                    phase = mod(phase, 2*pi) ;
                    phases(tt,ni) = phase ;
                    mags(tt, ni) = scaledMag ;
                    phase2row = abs(round((length(colors)/2/pi)*phase)) ; % this is the function that maps phase to a row in the colormap
                    color = colors(phase2row+1, :) ; % grab the color for this nucleus using its phase conversion
                    phandle = patch(bound(:,1), bound(:,2), color, 'FaceAlpha', scaledMag, 'EdgeColor', 'none') ;
                end
            end
        end
    end
    
    set(gca,'xdir','reverse')
    hold on
    phasebar(colors, 'rad', 'location', 'se', 'size', 0.2) ;
    mem = quiver(cidx(:,1), cidx(:,2), squeeze(binvels(1, tt, :)), squeeze(binvels(2, tt,:)), 1, 'LineWidth', 1, 'color', 'green') ;
    nuc = quiver(cidx(:,1), cidx(:,2), squeeze(binvels(3, tt,:)), squeeze(binvels(4, tt,:)), 1, 'LineWidth', 1, 'color', 'red') ;
    legend([mem, nuc] , 'Muscle', 'Endoderm', 'location', 'northeastoutside') ;
    text(500, 480 , '5 um', 'FontSize', 18, 'Color', 'black') ;
    if ~exist(fullfile(dataDir, 'phasemag_visualization'))
        mkdir(fullfile(dataDir, 'phasemag_visualization')) ;
    end
    saveas(gcf, fullfile(dataDir, 'phasemag_visualization', [sprintf('phasemag_%02d', tt), '.png'])) ;
end


%% Calculate Autocorrelation

lindt = round(linspace(-(ltp-1), (ltp -1)- 1, 2*(ltp-1))) ;
vdv = zeros(2, ltp-1, length(trackVals), 2*(ltp-1)) ;
Z1 = lindt' ;
Z2 = lindt' ;

for ch = 1:2
    % Assign array dims in binvels for each data channel
    if ch == 1
        ch1 = 1 ;
        ch2 = 2 ;
    elseif ch == 2
        ch1 = 3 ;
        ch2 = 4 ;
    end
    for tt = 1:(ltp-1)
        for ni = 1:length(trackVals)
            % The curent vector for this channel, tp, and nucleus
            vt(1) = binvels(ch1, tt, ni) ;
            vt(2) = binvels(ch2, tt, ni) ;
            for dt = 1:(ltp-1)
                % The vector at TP = tt + dt
                vdt(1) = binvels(ch1, dt, ni) ;
                vdt(2) = binvels(ch2, dt, ni) ;
                
                vdv(ch, tt, ni, (ltp-1) + dt - tt) = dot(vt, vdt) ; % Dot the current vector with the vector at TP = tt + dt
            end
            v = squeeze(vdv(ch, tt, ni, :)) ;
            if ch == 1
                Z1 = horzcat(Z1, v) ;
            elseif ch == 2
                Z2 = horzcat(Z2, v) ;
            end
        end
    end
end

vstats1 = binnedstats(Z1, 1) ;
vstats2 = binnedstats(Z2, 1) ;


%% Plot AutoCorrelation

close all
figure('units','normalized','outerposition',[0 0 1 1])
ylabel('V(t) $\cdot$ V(t+dt)', 'fontsize', 25, 'interpreter', 'latex') ;
xlabel('dt', 'fontsize', 25, 'interpreter', 'latex') ;
title('Autocorrelation', 'fontsize', 25, 'interpreter', 'latex') ;
hold on

%plot(lindt, vdv, '--', 'LineWidth', 1) ;