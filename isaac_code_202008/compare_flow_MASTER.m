%% COMPARING FLOW FIELDS MASTER -- single dataset
% Isaac Breinyn 2020
%
% Run setup.m before this script
%
% This is a script that loads multichannel confocal images and then (with
% user input and training along the way) can analyze the likeness of the
% flow between the two channels.
%
% Naming conventions must remain consistent across all datasets for this
% pipeline to work reliably. Access the flowExperiment.m file for
% information on how to name files.
%
% NOTE: The PIV data used in this script is formated into arrays of u and v
% component of the velocities for EACH channel. Therefore, arrays with a
% first dim of 4 will be in the order of u1, v1, u2, then v2. The letter
% corresponds to the component of the velocity vector, and the number to
% the prevelant channel.

%% Run setup.m from imsane path 
% setup

%% Clean Matlab and add Paths
clear; close all; clc ;

% add paths to where the code lives
% addpath_recurse('L:\\Streichan\\code\\gut_matlab\\') ;
% addpath_recurse('L:\\Streichan\\code\\') ;
addpath_recurse('/mnt/data/code/gut_matlab/') ;
addpath('/mnt/data/code/isaac_code/')

%% WARNING: Users should ONLY make changes to this section of the script.

%%% Assign paramaters for class flowExperiment %%%

% options.dataDir = 'L:\\Streichan\\data\\202003111630_mef2gal4klarUASCAAXmChHiFP_wo_great\\1644_folding_30s_2um_la8_4p5zoom_63x\\' ; % the directory where the data lives
options.dataDir = fullfile([filesep 'mnt'],'data','confocal_data',...
    'gut', 'relative_motion', 'mef2Gal4klarUASCAAXmChHiFP', ...
    '202003111630_mef2gal4klarUASCAAXmChHiFP_wo_great',...
    '1644_folding_30s_2um_la8_4p5zoom_63x') ; % the directory where the data lives
options.dataSetName = '1644_folding' ; % the "label" of this dataset (the prefix to the file names)
options.lastTimePoint = 32 ; % "last timepoint" or the number of TPs in your dataset
options.fileSize = [512 512] ; % the size of the file
options.trackFileName = fullfile('Ch2_MIPs','1644_folding_Ch2_T01_MIP_smoothed_stack-myData_Object Identities.h5') ; % the name of the object identities file (from ilastik)
options.res = 0.0802 ; % the data resolution in um/px
options.lobeSplit = 256 ; % define where you want to put the barrier that defines where the two lobes seperate
options.dt = 0.5 ;  % time resolution in minutes

% todo: add mirrored boolean

save(fullfile(options.dataDir, 'options.mat'), 'options')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(options.dataDir)
%% Instantiate the class flowExperiment and populate necessary properties

flex = flowExperiment(options) ;

%% The first step is to acquire data that you would like to analyze. This data must be split into two channels and seperate timepoints.
% NOTE: MAKE SURE DATA IS IN THIS ORIENTATION
%
%      D
%  A      P
%      V
%
% Open LIF in ImageJ, Flip/rotate in imageJ, save as separate tiff for each
% timepoint and channel. 
% 
% Once in the correct orientation, rescale the data to unit aspect ratio 
% and then convert into h5s for ilastik training using 
% muscle_surface_pipeline.m. Only the membrane channel (channel 1)
% should be used for ilastik training. This is the step necessary to remove
% any glial cells or non-gut structures from the FOV. 

%% ILASTIK: Train on ch1 and extract probabilities
% foreground channel (ch1 of the training) is the membrane, background (ch2
% of the training) is else. 

%% Mask data and then extract MIPs 
% Once you have the probabilities for channel 1 of your data, you can mask
% both channels and make MIPs of the data. These MIPs will then be used
% later in the pipeline.

flex.makeMaskedMIPs() ; 

% After running, you should have a collection of masked tiffs as well as
% MIPs of those masked tiffs.

%% Smooth the MIPs from last step
% This gets rid of any salt-and-pepper noise and allows for higher quality
% PIV analysis.

flex.smoothMIPs() ;

% Now that you have smoothed MIPs of both channels, there are two things
% that need to be done. 

%% Extract and load tracking data
% Load the masked MIPs of the nuclear data (channel 2) into a pixel
% classification + object identification project. Train on nuclei and
% export an object identification file. Once you have done this, load in
% the tracking data.

flex.loadTrackData() ;

%% PIVlab training
% Load smoothed MIPs of both channels into their own PIVlab sessions to
% export the flow of each channel.
% Save the PIV output as File > Export > mat File, with filenames:
% flex.ch1FlowFileName = fullfile(flex.dataDir, 'PIVlab_membrane.mat') ;
% flex.ch2FlowFileName = fullfile(flex.dataDir, 'PIVlab_nuclei.mat') ;

%% Build Full Flow Field using flow data
% Once you have exported the flow of each channel, build an array that
% contains the flow data of each channel and store it in flex

flex.buildFullFlow() ;

%% Mask and Bin the Flow Field
% Now that the flow field has been built, use the masked data to mask and
% bin the flow field into discrete nuclear segments.

flex.maskBinFlow() ;

%% Flow Phase and Magnitude Visualization
% This step is extremely important. Not only will this step
% calculate the relative phases of flow, but it will create a series of
% plots that will aid in any debugging. These plot will color the nuclei
% based on the relative phase between each channel's flow vector, as well
% as add an opacity to the nucleus that corresponds to the ratio between
% that channel's velocity vector, and the max velocity vector value of that
% dataset. This can provide information on velocity range.

flex.phaseMagVisualization() ;

%% Create Master Plots
% The final step is to organize all of the data. This step will create a
% master plot that contains a variety of information including velocity
% distributions, and interchannel similarity. This step will also create
% two similar plots, each reserved for indivudal lobes in the FOV

flex.createMasterPlot() ;




