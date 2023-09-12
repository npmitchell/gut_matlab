%% Create minimal data for figure 3

addpath('/mnt/data/code/tubular/utility/plotting/')

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

%% Tangential velocity as vector phasemap

top = [0,90] ;
washout2d = 0 ; 

lambda = tubi.smoothing.lambda ;
lambda_mesh = tubi.smoothing.lambda_mesh ;
lambda_err = tubi.smoothing.lambda_err ;
nmodes = tubi.smoothing.nmodes ;
zwidth = tubi.smoothing.zwidth ;
nU = tubi.nU ;
nV = tubi.nV ;
climit = 0.2 ;
climit_veln = climit * 10 ;

outputdir = fullfile(tubi.dir.data, 'TubULAR_Figures', 'Figure03a') ;
if ~exist(outputdir, 'dir')
    mkdir(outputdir)
end

% Load piv results
disp('Obtaining raw piv to get field size')
tubi.getPIV() ;
piv = tubi.piv.raw ;
% get size of images to make
gridsz = size(tubi.piv.raw.x{1}) ;
% Define Nx1 and Mx1 float arrays for xspace and yspace
xx = tubi.piv.raw.x{tidx}(1, :) ;
yy = tubi.piv.raw.y{tidx}(:, 1) ;

% Get tangential velocities from disk
tubi.getVelocityAverage()
velstruct = tubi.velocityAverage ;
v2dsmMum = velstruct.v2dum ;
vtscale = 5 ; % um/min
qsubsample = 10 ; 
t0 = tubi.t0 ;

%% Compute or load all timepoints
clc
for tp = tubi.t0 + [30,60,90] 
    tidx = tubi.xp.tIdx(tp) ;
    close all
    disp(['t = ' num2str(tp)])
    tidx = tubi.xp.tIdx(tp) ;

    tubi.setTime(tp) ;
    v2dsmum_ii = squeeze(v2dsmMum(tidx, :, :)) ;
    vx = reshape(v2dsmum_ii(:, 1), gridsz) ;
    vy = reshape(v2dsmum_ii(:, 2), gridsz) ;

    vxb = imgaussfilt(vx, 4) ;
    vyb = imgaussfilt(vy, 4) ;

    % read pullback image
    im = imread(sprintf(tubi.fullFileBase.im_sp_sme, tp)) ;
    ylims = [0.25 * size(im, 1), 0.75 * size(im, 1)] ;
    imw = ones(size(im)) ;

    qopts = struct() ;
    qopts.qsubsample = qsubsample ;
    qopts.overlay_quiver = true ;
    qopts.qscale = 10 ;
    qopts.quiver_line_width = 1.5;
    qopts.label = '$v_t$ [$\mu$m/min]' ;
    qopts.title = ['Smoothed velocity, $v_t$: $t=$', ...
        num2str(tp - t0), tubi.timeUnits] ;
    qopts.ylim = ylims ;
    xyf = struct() ;
    xyf.x = xx ;
    xyf.y = yy ;

    % Plot the coarse-grained tangential velocity as heatmap on top of im
    [h1, h2, h3, ax, cax, pbax] =vectorFieldHeatPhaseOnImage(imw, xyf, vxb, vyb, vtscale, qopts) ;    
    outfn = fullfile(outputdir, sprintf('vtangential_%06d.fig', tp)) ;
    saveas(gcf, outfn)
    outfn = fullfile(outputdir, sprintf('vtangential_%06d_fig.png', tp)) ;
    saveas(gcf, outfn)

    set(gcf, 'CurrentAxes', ax)
    F = getframe() ;
    outfn = fullfile(outputdir, sprintf('vtangential_%06d.png', tp)) ;
    imwrite(F.cdata, outfn)

    outfn = fullfile(outputdir, sprintf('raw_data_Figure3a_tangential_velocity_snapshot%06d_vx.txt', tp)) ;
    header = 'in-plane tangential velocity field projected onto pullback coordinates, component in s direction  [um/min]. Note the velocity field double-covers the surface, such that the middle half of the field in phi is a single cover.' ;
    write_txt_with_header(outfn, vxb, header)

    outfn = fullfile(outputdir, sprintf('raw_data_Figure3a_tangential_velocity_snapshot%06d_vy.txt', tp)) ;
    header = 'in-plane tangential velocity field projected onto pullback coordinates, component in phi direction [um/min]. Note the velocity field double-covers the surface, such that the middle half of the field in phi is a single cover.' ;
    write_txt_with_header(outfn, vyb, header)

    dataPlotted = struct('imw', imw, 'xyf', xyf, 'vxb', vxb, 'vyb', vyb, 'vtscale', vtscale, 'qopts', qopts) ;
    outfn = fullfile(outputdir, sprintf('dataPlotted_%06d.mat', tp)) ;
    save(outfn, 'dataPlotted')
end

