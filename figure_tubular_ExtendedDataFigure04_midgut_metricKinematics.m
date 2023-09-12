%% PANELS A-C
%  See figure_tubular_Figure03_metricKinematics_midgut for panels a-c

%% PANEL D
% See
% figure_tubular_ExtendedDataFigure04_midgut_metricKinematics_correlationPlot.m
% on local machine (Mac OS in Dropbox) for panel d

%% PANEL E --> run the following

addpath('/mnt/data/code/gut_matlab/plotting/')

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

%% load xp
disp('loading xp struct from disk')
load(fullfile(dataDir, 'xp.mat'), 'xp', 'opts')

%% TubULAR class instance
disp('defining TubULAR class instance (tubi= tubular instance)')
tubi = TubULAR(xp, opts) ;
disp('done defining TubULAR instance')

%%
% Plot Pathline Kinematics
% options = struct() ;
% tubi.plotPathlineMetricKinematics(options)



mKDir = fullfile(tubi.dir.metricKinematics.root, ...
    strrep(sprintf(['' 'lambda%0.3f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d'], ...
    lambda, lambda_mesh, lambda_err, nmodes, zwidth), '.', 'p'));

%% Get fold locations for pathlines
fons = t0 - tubi.xp.fileMeta.timePoints(1) ;

%% Colormap
% bwr256 = bluewhitered(256) ;
bwr256 = colormap(brewermap(256, '*RdBu')) ;

%% Test incompressibility of the flow on the evolving surface
% We relate the normal velocities to the divergence / 2 * H.
tps = tubi.xp.fileMeta.timePoints(1:end-1) - t0;

% Output directory is inside metricKinematics dir
mKPDir = fullfile(mKDir, sprintf('pathline_%04dt0', t0Pathline)) ;
datdir = fullfile(mKPDir, 'measurements') ;

% colormap prep    
close all
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;
cmap = brewermap(512, '*RdBu') ;
close all

% Compute or load all timepoints
apKymoFn = fullfile(datdir, 'apKymographMetricKinematics.mat') ;

tmp = load(apKymoFn) ;
gdot_apM = tmp.gdot_apM ;

%% crop the matrix
edge = 5 ;
sinclude = edge:tubi.nU-edge ;
tidx0 = tubi.xp.tIdx(tubi.t0) ;
tps = (0:90)/60 ;
ss = (sinclude)/tubi.nU ;
gdot = gdot_apM(tidx0:tidx0+90, sinclude) ;

% initialize plot with kymo
close all
figure('Units', 'centimeters', 'Position', [0,0,4.4,4.4])
set(gca, 'fontsize', 7)
imagesc(ss, tps, gdot) ;
clim([-0.05, 0.05])
colormap(cmap)

% add folds
% tubi.getFeatures ;
% tubi.folds

xlim([0, 1])
ylim([0,1.5])

outputdir = fullfile(tubi.dir.data, 'TubULAR_Figures', 'ExtendedDataFigure4e') ;
if ~exist(outputdir, 'dir')
    mkdir(outputdir)
end
saveas(gcf, fullfile(outputdir, 'ExtendedDataFigure4e.pdf' ))
colorbar
saveas(gcf, fullfile(outputdir, 'ExtendedDataFigure4e_colorbar.pdf' ))

% Save data
save(fullfile(outputdir, 'ExtendedDataFigure4e.mat'), 'ss', 'tps', 'gdot')

fn = fullfile(outputdir, 'raw_data_ExtendedDataFigure4e_x.txt') ;
dat = ss ;
header = 'longitudinal position s of strain rate measurements, measured in the material frame of reference, as a fraction of total length' ;
write_txt_with_header(fn, dat, header)

fn = fullfile(outputdir, 'raw_data_ExtendedDataFigure4e_y.txt') ;
dat = tps ;
header = 'timestamps relative to constriction onset of strain rate measurements [min]' ;
write_txt_with_header(fn, dat, header)

fn = fullfile(outputdir, 'raw_data_ExtendedDataFigure4e_arealStrainRate.txt') ;
dat = gdot ;
header = 'areal strain rate 1/2 * Trace[g^{-1} dg/dt], in units of [1/min]' ;
write_txt_with_header(fn, dat, header)
