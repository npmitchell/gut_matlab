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

outputdir = fullfile(tubi.dir.data, 'TubULAR_Figures', 'Figure1') ;
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
vnsmM = velstruct.vn ;
vtscale = 5 ; % um/min
vnscale = 2 ;
qsubsample = 10 ; 
t0 = tubi.t0 ;
alphaVal = 0.7 ; 

%% Compute or load all timepoints
invertImage = true ;
clc
for tp = tubi.t0 + [69] 
    tidx = tubi.xp.tIdx(tp) ;
    close all
    disp(['t = ' num2str(tp)])
    tidx = tubi.xp.tIdx(tp) ;

    tubi.setTime(tp) ;
    v2dsmum_ii = squeeze(v2dsmMum(tidx, :, :)) ;
    vnsm_ii = squeeze(vnsmM(tidx, :, :)) ;
    vx = reshape(v2dsmum_ii(:, 1), gridsz) ;
    vy = reshape(v2dsmum_ii(:, 2), gridsz) ;
    vn = reshape(vnsm_ii, gridsz) ;

    vxb = imgaussfilt(vx, 4) ;
    vyb = imgaussfilt(vy, 4) ;

    % read pullback image
    im = imread(sprintf(tubi.fullFileBase.im_sp_sme, tp)) ;
    ylims = [0.25 * size(im, 1), 0.75 * size(im, 1)] ;
    imw = (max(im(:))-im) * washout2d + max(im(:)) * (1-washout2d) ;


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
    close all
    fig = figure('units', 'centimeters', ...
            'outerposition', [0 0 36 18], 'visible', 'off') ;
    [h1, h2, h3, ax, cax, pbax] =vectorFieldHeatPhaseOnImage(imw, xyf, vxb, vyb, vtscale, qopts) ;    
    outfn = fullfile(outputdir, sprintf('vtangential_%06d.fig', tp)) ;
    saveas(gcf, outfn)
    outfn = fullfile(outputdir, sprintf('vtangential_%06d_fig.png', tp)) ;
    saveas(gcf, outfn)

    set(gcf, 'CurrentAxes', ax)
    F = getframe() ;
    outfn = fullfile(outputdir, sprintf('vtangential_%06d.png', tp)) ;
    imwrite(F.cdata, outfn)

    outfn = fullfile(outputdir, sprintf('raw_data_Figure1c_tangential_velocity_snapshot%06d_vx.txt', tp)) ;
    header = 'in-plane tangential velocity field projected onto pullback coordinates, component in s direction  [um/min]. Note the velocity field double-covers the surface, such that the middle half of the field in phi is a single cover.' ;
    write_txt_with_header(outfn, vxb, header)

    outfn = fullfile(outputdir, sprintf('raw_data_Figure1c_tangential_velocity_snapshot%06d_vy.txt', tp)) ;
    header = 'in-plane tangential velocity field projected onto pullback coordinates, component in phi direction [um/min]. Note the velocity field double-covers the surface, such that the middle half of the field in phi is a single cover.' ;
    write_txt_with_header(outfn, vyb, header)

    dataPlotted = struct('imw', imw, 'xyf', xyf, 'vxb', vxb, 'vyb', vyb, 'vtscale', vtscale, 'qopts', qopts) ;
    outfn = fullfile(outputdir, sprintf('raw_data_Figure1c_tangential_velocity_dataPlotted_%06d.mat', tp)) ;
    save(outfn, 'dataPlotted')


    % Normal velocity
    close all
    fig = figure('units', 'centimeters', ...
            'outerposition', [0 0 36 18], 'visible', 'off') ;
    % axes('Position', [0.1 0.05, 0.85, 0.4])
    % qopts.axPosition = [0.2 0.55, 0.85, 0.4] ;
    labelOpts = struct() ;
    % labelOpts.cbPosition = [.9 .1 .02 .3] ;
    labelOpts.label = '$v_n$ [$\mu$m/min]' ;
    labelOpts.title = ['normal tissue velocity, $v_n$'] ;
    labelOpts.cmap = bwr ;
    im = imread(sprintf(tubi.fullFileBase.im_sp_sme, tp)) ;
    if invertImage
        imw = max(im(:)) - im ;
        % imw = (max(im(:))-im) * washout2d + max(im(:)) * (1-washout2d) ;
    else
        imw = im * washout2d + max(im(:)) * (1-washout2d) ;
    end
    imw = cat(3, imw,imw,imw) ;
    [tmpx, tmpy] = meshgrid(xx, yy) ;
    tmpxy = [tmpx(:), tmpy(:)] ;
    scalarFieldOnImage(imw, tmpxy, vn, alphaVal, vnscale, ...
        labelOpts) ;
    ylim(ylims)

    % save the figs
    outfn = fullfile(outputdir, sprintf('vn_%06d.fig', tp)) ;
    saveas(gcf, outfn)
    outfn = fullfile(outputdir, sprintf('vn_%06d_fig.png', tp)) ;
    saveas(gcf, outfn)

    F = getframe() ;
    outfn = fullfile(outputdir, sprintf('vn_%06d.png', tp)) ;
    imwrite(F.cdata, outfn)

    outfn = fullfile(outputdir, sprintf('raw_data_Figure1c_normal_velocity_snapshot%06d_vx.txt', tp)) ;
    header = 'normal velocity sampled in (s,phi) pullback space  [um/min]. Note the velocity field double-covers the surface, such that the middle half of the field in phi is a single cover.' ;
    write_txt_with_header(outfn, vn, header)

    dataPlotted = struct('imw', imw, 'xy', tmpxy, 'vn', vn, 'vnscale', vnscale) ;
    outfn = fullfile(outputdir, sprintf('raw_data_Figure1c_normal_velocity_dataPlotted_%06d.mat', tp)) ;
    save(outfn, 'dataPlotted')


    % %% Other version of normal velocity
    % top = [0,90] ;
    % % cmap1 = brewermap(256, '*RdBu') ;
    % cmap = bwr ;
    % cmapname = 'bwr' ;
    % 
    % 
    % lambda = tubi.smoothing.lambda ;
    % lambda_mesh = tubi.smoothing.lambda_mesh ;
    % lambda_err = tubi.smoothing.lambda_err ;
    % nmodes = tubi.smoothing.nmodes ;
    % zwidth = tubi.smoothing.zwidth ;
    % nU = tubi.nU ;
    % nV = tubi.nV ;
    % climit = 0.2 ;
    % climit_veln = climit * 10 ;
    % 
    % mKDir = fullfile(tubi.dir.metricKinematics.root, ...
    %     strrep(sprintf(['lambda%0.3f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d'], ...
    %     lambda, lambda_mesh, lambda_err, nmodes, zwidth), '.', 'p'));
    % outdir = fullfile(mKDir, 'measurements') ;
    % 
    % close all
    % disp(['t = ' num2str(tp)])
    % tidx = tubi.xp.tIdx(tp) ;
    % 
    % tubi.setTime(tp) ;
    % mesh = tubi.getCurrentSPCutMeshSmRSC ;
    % cutMesh = tubi.getCurrentSPCutMeshSmRS ;
    % 
    % % Check for timepoint measurement on disk
    % nfn = fullfile(outdir, sprintf('veln_vertices_%06d.mat', tp)) ;
    % 
    % % Load timeseries measurements
    % load(nfn, 'veln_filt', 'veln_ap', 'veln_l', 'veln_r', 'veln_d', 'veln_v')
    % 
    % % separate 2d/3d data
    % veln2d = veln_filt ;
    % 
    % % plot veln in 2d
    % clf
    % set(gcf, 'Units', 'centimeters')
    % set(gcf, 'Position', [0,0,40,20])
    % set(gcf, 'color', 'w')
    % ss = cutMesh.u(:, 1) ;
    % ssvals = tubi.a_fixed * ss / max(ss) ;
    % trisurf(cutMesh.f, ssvals, ...
    %     cutMesh.u(:, 2), zeros(size(cutMesh.u(:,1))),...
    %     veln2d, 'edgecolor', 'none')
    % xlim([0, tubi.a_fixed])
    % ylim([0, 1]) 
    % caxis([-climit_veln, climit_veln]) 
    % axis off
    % view(top)
    % fillfig
    % F = getframe() ;
    % outfn = fullfile(outputdir, sprintf('veln2d_%06d.png', tp)) ;
    % imwrite(F.cdata, outfn)
    % outfn = fullfile(outputdir, sprintf('veln2d_%06d.txt', tp)) ;
    % header = 'normal velocity [um/min]' ;
    % write_txt_with_header(outfn, veln3d, header)

    
end


%% PANEL B
addpath('/mnt/data/code/tubular/utility/bfmatlab/')
% copy pullback images from sphi at 123 and 163, 168?
% 
% Establish texture patch options
metadat = struct() ;
metadat.reorient_faces = false ;            % set to true if some mesh normals may be inverted (requires gptoolbox if true)
metadat.normal_shift = tubi.normalShift ;   % normal push, in pixels, along normals defined in data XYZ space
metadat.texture_axis_order = [1 2 3] ;      % texture space sampling. If the surface and dataspace have axis permutation, enter that here
Options = struct() ; 
Options.PSize = 5 ;          % Psize is the linear dimension of the grid drawn on each triangular face. Set PSize > 1 for refinement of texture on each triangle of the surface triangulation. Higher numbers are slower but give more detailed images.
Options.numLayers = [0, 0];  % how many layers to MIP over/bundle into stack, as [outward, inward]
Options.layerSpacing = 2 ;   % Distance between layers over which we take MIP, in pixels, 

metadat.plot_dorsal = false ;
metadat.plot_ventral = false ;
metadat.plot_left = true ;
metadat.plot_right = false ;
metadat.plot_perspective = false ;
metadat.blackFigure = false ;
metadat.plot_time_points = [tubi.xp.tIdx(123), tubi.xp.tIdx(163), tubi.xp.tIdx(168),  tubi.xp.tIdx(173)] ;

% Plot on surface for all timepoints 
tubi.plotSeriesOnSurfaceTexturePatch(metadat, Options)


