%% UVPrime coordinate manupulation and measurements

%% Create u'v' coordinate cutMeshes (conformal with minimally twisted boundaries)
options = struct() ;
options.overwrite = false ;
options.save_ims = true ;
options.preview = false ;
QS.generateUVPrimeCutMeshes(options)

%% TexturePatch image pullbacks 
for tidx = 1:length(QS.xp.fileMeta.timePoints)
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    QS.setTime(tp) ;
    pbOptions.overwrite = false ;
    pbOptions.generate_uv = false ;
    pbOptions.generate_uphi = false ;
    pbOptions.generate_relaxed = false ;
    pbOptions.generate_rsm = true ;
    pbOptions.generate_uvprime = true ;
    pbOptions.generate_ruvprime = false ;
    pbOptions.numLayers = [7, 7] ;  % previously [5,5]
    pbOptions.layerSpacing = 0.75 ;
    QS.generateCurrentPullbacks([], [], [], pbOptions) ;
end

% TILE/EXTEND (smoothed) UVPrime images in Y and resave =======================================
% Skip if already done
options = struct() ;
options.overwrite = false ;
options.coordsys = 'uvprime' ;
QS.doubleCoverPullbackImages(options)
disp('done')

%% Perform PIV on Extended uvprime coordinates for quasiconformal measurement
% Perform this in PullbackImages_XXstep_uvprime/piv/
%
% % Compute PIV in PIVLab
% % ---------------------
% % Open PIVLab
% % Select all frames in meshDir/PullbackImages_010step_sphi/smoothed_extended/
% % Select Sequencing style 1-2, 2-3, ... 
% % Image Preprocessing (used to select all, but now:)
% %  --> CHOOSE EITHER WAY: Enable CLAHE with >=20 pix
% %  --> DO NOT Enable highpass with 15 pix
% %  --> DO NOT Enable Intensity capping
% %  --> Wiener2 denoise filter with 3 pix
% %  --> DO NOT Auto constrast stretch
% % PIV settings: 
% %  --> 128 (32 step), 64 (32 step), 32 (16 step), 16 (8 step)
% %  --> Linear window deformation interpolator
% %  --> 5x repeated correlation 
% %  --> Disable auto-correlation
% % Post-processing
% %  --> Standard deviation filter: 7 stdev
% %  --> Local median filter: thres=5, eps=0.1
% %  --> Interpolate missing data
% % Export 
% %  --> File > Save > MAT file

%% Pathlines in uv' coords
options = struct() ;
options.overwrite = false ;
options.overwriteImages = false ;
QS.measureUVPrimePathlines(options)

%%
options = struct() ;
options.overwrite = false ;
options.overwriteImages = false ;
options.climit = 1 ;
options.coordSys = 'uvprime' ;
QS.measureBeltramiCoefficient(options)

%% Compare to nonlinear isothermal description
options = struct() ;
options.overwrite = true ;
options.overwriteImages = false ;
options.coordSys = 'uvprime' ;
QS.compareBeltramiToConstriction(options) ;

%% Compare to linearized description
options = struct() ;
options.overwrite = true ;
options.overwriteImages = false ;
QS.compareBeltramiToLinearizedConstriction(options) ;