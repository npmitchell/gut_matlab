%% MIDGUT MUSCLE ANALYSIS
% Isaac B. Breinyn 2020
%
% A script for analyzing the muscle layer of the Drosophila midgut after
% prepared via an ImSaNE pipeline as well as trained on using Ilastik
%
% TO DO BEFORE RUNNING: Run setup.m from ImsaneV1.2.3 and train on the
% muscle data in Ilastik. Have outputted probabilities ready.

%% Clean Matlab and set Global Parameters

% clean the terminal and all variables, close all open figures

clear; close all; clc;

% add any necessary paths

addpath_recurse('/mnt/data/code/gut_matlab/') ;

% set all Global Parameters here

dataDir = '/mnt/data/confocal_data/gut/mef2Gal4klarUASCAAXmChHiFP/202003111630_mef2gal4klarUASCAAXmChHiFP_wo_great/1644_folding_30s_2um_la8_4p5zoom_63x/' ; % working directory
fn = '1644folding_T%02d' ; % working file name
ssfactor = 2 ; % the subsampling factor of your muscle data (used later to return to origonal aspect ratio)
ltp = 33 ; % the last TP of your datasets (this will be used for iterating over your dataset)
ifg = 1   ; % the channel that corresponds to the fg of your outputted probabilities (1 or 2) (ilastik foreground)
dataChannel = 1 ; % the channel used for your raw data (1 or 2)

%% Load and Mask Data using Training

cd(dataDir)
fnp = [fn, '_Probabilities'] ; % Outputted probabilities file name
fnm = [fn, '_masked'] ; % Masked data file name
ssfArray = ones(ssfactor, ssfactor, ssfactor) ; % create cubic 3D array the size of your ssfactor in each dim to perform kronecker tensor product

for tt = 0:ltp % iterate over timepoints in the datasaet
%     rawData = loadtiff([sprintf(fn,tt), '.tif']) ; % load raw data of this particular timepoint (loadtiff courtesy of YoonOh Tak)
%     if dataChannel == 1
%         rawData = rawData(:,:,1:12) ;
%     elseif dataChannel == 2
%         rawData = rawData(:,:,13:24) ;
%     end
    
    rawData = h5read([sprintf(fn,tt), '.h5'], '/inputData') ; % load raw data of this particular timepoint 
    rawData = squeeze(rawData(dataChannel, :, :, :)) ; % assign which channel is used, and remove that dim
    rawData = permute(rawData, [3 2 1]) ; % permute the raw data array
    
    probUnitAspect = h5read([sprintf(fnp,tt), '.h5'], '/exported_data') ; % outputed prob from Ilastik corresponding to this timepoint
    probUnitAspect = squeeze(probUnitAspect(ifg, :, :, :)) ; % assign which channel is used, and remove that dim
    probUnitAspect = permute(probUnitAspect, [3 2 1]) ; % permute the prob array
    
    %probOrigAspect= superkron(probUnitAspect, ssfArray) ; % take kronecker tensor product of unit aspect ratio prob and subsampling array
    
    maskData = rawData.*probUnitAspect ; % mask the origonal data using the origonal aspect ratio probabilities from ilastik
    saveastiff(maskData, fullfile(dataDir, [sprintf(fn, tt), '_masked.tif']))  ;
    disp(['Saving ', sprintf(fnm, tt), '.tif']) ;
end