% First make the pullback image and save the scalebar info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Redo Pullbacks with time-smoothed meshes ===============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Skip if already done
disp('Create pullback using S,Phi coords with time-averaged Meshes')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tt = QS.t0set() % , QS.xp.fileMeta.timePoints]
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    tidx = QS.xp.tIdx(tt);
    
    % Load the data for the current time point ------------------------
    QS.setTime(tt) ;
    
    % OPTION 1: Keep constant luminosity throughout, modify default 
    % intensity limits.
    % if tidx == 1        
    %     % Use first timepoint's intensity limits throughout
    %     QS.setDataLimits(QS.xp.fileMeta.timePoints(1), 1.0, 99.995)
    % end
    % QS.data.adjustlow 
    % QS.data.adjusthigh
    
    % OPTION 2 : Adjust intensity to scale from timepoint to timepoint
    adjustlow = 1.00 ;         % floor for intensity adjustment
    adjusthigh = 99.9 ;        % ceil for intensity adjustment (clip)
    QS.data.adjustlow = adjustlow ;
    QS.data.adjusthigh = adjusthigh ;
    
    % Establish custom Options for MIP
    pbOptions = struct() ;
    pbOptions.overwrite = true ;
    pbOptions.numLayers = [0 0] ; % previously [7, 7] ;  % previously [5,5]
    pbOptions.layerSpacing = 0.75 ;
    pbOptions.generate_rsm = false ;
    pbOptions.generate_spsm = false ;
    pbOptions.generate_sphi = false ;
    pbOptions.generate_uvprime = false ;
    pbOptions.generate_ruvprime = false ;
    pbOptions.generate_ricci = true ;
    pbOptions.imSize = 2000 ;
    QS.data.adjustlow = adjustlow ;
    QS.data.adjusthigh = adjusthigh ;
    QS.generateCurrentPullbacks([], [], [], pbOptions) ;
end


%% Now load the image, select a region and get the scale in that region
datdir = '/mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/msls_output/gridCoords_nU0100_nV0100/PullbackImages_010step_ricci/';
imfn = fullfile(datdir, 'MAX_Time_000123_c1_stab_ricci_z11to14.tif') ;
ggfn = fullfile(datdir, 'Time_000123_c1_stab_ricci_info.mat') ;
scaleBarum= 10 ;

% get the scale at a particular point in the image
im = imread(imfn) ;
dat = load(ggfn) ;
clf
imshow(im)
XY = ginput(1) ;
idx = pointMatch(XY, dat.bcXY) ;
hold on ; 
pix2um = mean(dat.PBpix2um(idx, :))

disp(['scale is um/pix: ' num2str(pix2um)])

% Now crop image and save
xlims = min(round(xlim), size(im, 2)) ;
ylims = min(round(ylim), size(im, 1)) ;
snap = im(ylims(1):ylims(2), xlims(1):xlims(2)) ;
% imwrite(snap, fullfile(datdir, 'snap3_leftVentral_noScale.tif')) ;
imwrite(snap, fullfile(datdir, 'snap3_leftDorsal_noScale.tif')) ;
snap(15:20, 240:240+round(scaleBarum/pix2um)) = 255 ;
figure ; imshow(snap)
%imwrite(snap, fullfile(datdir, sprintf('snap3_leftVentral_wScale_%02dum.tif', scaleBarum))) ;
imwrite(snap, fullfile(datdir, sprintf('snap4_leftDorsal_wScale_%02dum.tif', scaleBarum))) ;

% Check that we are sampling in the right place on original image
figure ; 
imshow(im) ;
plot(dat.bcXY(:, 1), dat.bcXY(:, 2), '.')
plot(dat.bcXY(idx, 1), dat.bcXY(idx, 2), 'o')


