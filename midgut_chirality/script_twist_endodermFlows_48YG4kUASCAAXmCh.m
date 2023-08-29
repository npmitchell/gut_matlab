%% EXAMPLE MASTER PIPELINE FOR TIME SERIES DATA (3D + time)
% by NPMitchell & Dillon Cislo

% This is a pipeline to analyze dynamic tube-like surfaces in 3D data.
% A tube-like surface is one that is either cylindrical in topology or
% elongated and spherical in topology. If the initial mesh is a spherical
% topology, then two 'endcaps' will be truncated in order to transform it 
% to a cylindrical topology.

%% Clear workspace ========================================================
% We start by clearing the memory and closing all figures
clear; close all; clc;
cd /mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/
zwidth =1 ;
% cd /mnt/data/48Ygal4-UAShistRFP/201904031830_great/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2/data/deconvolved_16bit/
% zwidth = 1 ;
% cd /mnt/data/handGAL4klarHandGFPhistGFP/202105072030_1mWGFP/deconvolved_16bit/
% % cd /mnt/data/mef2GAL4klarUASCAAXmChHiFP/202003151700_1p4um_0p5ms3msexp/
% zwidth = 2 ;


dataDir = cd ;

%% ADD PATHS TO THIS ENVIRONMENT ==========================================
disp('adding paths...')
origpath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(origpath))
addpath(fileparts(origpath))
addpath(genpath('../'))
addpath(genpath('../utility'))
addpath(genpath('/mnt/data/code/gptoolbox'))
addpath('../TexturePatch')
addpath('../DECLab')
addpath(fullfile('../utility','plotting'))
addpath(fullfile('../utility','plotting'))
% go back to the data
cd(dataDir)
disp('done adding paths')

if ~exist(fullfile(dataDir, 'xp.mat'), 'file')
    error('Could not load xp.mat')
else
    disp('loading xp struct from disk')
    load(fullfile(dataDir, 'xp.mat'), 'xp', 'opts')
end

%% TubULAR class instance
disp('defining TubULAR class instance (tubi= tubular instance)')
tubi = TubULAR(xp, opts) ;
disp('done defining TubULAR instance')

%% Measure the Twist in the flows 
tubi.measureTwist()

%% Some other extra functionality: Show accumulated strain as magnitude
tubi.plotPathlineStrain

%% Some other extra functionality: Visualize Lagrangian patch followed by piv Pathlines
opts = struct() ;
opts.timePoints = 96:10:206 ;
opts.demoPatchName = 'demoPatch005' ;
tubi.visualizeSegmentationPatch(opts) 

%% Some other extra functionality: Look at the rms velocity over time
[vrms, timestamps] = tubi.measureRMSvelocityOverTime ;

