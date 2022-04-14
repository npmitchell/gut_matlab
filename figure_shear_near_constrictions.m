%% eLife reviews for anisotropy 

datdir = '/Users/npmitchell/Dropbox/Soft_Matter/PAPER/gut_paper/figure_drafting/figure2_kinematics/seg3d_corrected' ;
fns = dir(fullfile(datdir, 'Time*.mat')) ;

% Colormap
addpath(genpath('~/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/gut_matlab/plotting/'))
colors = colormap(0.9 * cmap('I1', 'N', length(fns))) ;


%% WARMUP: UV coordinates
% Plot each timepoint
clf ;
for tidx = 1:length(fns)
    fn = fullfile(fns(tidx).folder, fns(tidx).name) ;
    seg = load(fn) ;
    seg = seg.seg3d ;
    stat = seg.statistics ;
    ap = stat.apStats.aspectWeighted.apBins ;
    c2t = stat.apStats.aspectWeighted.apCos2Theta ;
    cs = stat.apStats.aspectWeighted.apCos2ThetaStd ;
    ce = stat.apStats.aspectWeighted.apCos2ThetaSte ;
    
    % Plot the anisotropy
    lineProps = {'-','color', colors(tidx, :)} ;
    % shadedErrorBar(ap, c2t, cs, 'lineProps', lineProps)
    hold on;
    shadedErrorBar(ap, c2t, ce, 'lineProps', lineProps)
    
end

cb = colorbar ;

%% Find where cells go
% Quick instantiation of a QuapSlap object 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace ========================================================
% We start by clearing the memory and closing all figures
clear; close all; clc;
% change this path, for convenience
cd /mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data

dataDir = cd ;
meshDir = fullfile(dataDir, 'deconvolved_16bit', 'msls_output') ;

% PATHS ==================================================================
origpath = matlab.desktop.editor.getActiveFilename;
cd(fullfile(fileparts(origpath), 'master_pipeline'))
aux_paths_and_colors
cd(dataDir)

% LOAD EXISTING MASTER SETTINGS
disp('Loading masterSettings from ./masterSettings.mat')
load('./masterSettings.mat', 'masterSettings')
% Unpack existing master settings
stackResolution = masterSettings.stackResolution ;
nChannels = masterSettings.nChannels ;
channelsUsed = masterSettings.channelsUsed ;
timePoints = masterSettings.timePoints ;
ssfactor = masterSettings.ssfactor ;
% whether the data is stored inverted relative to real position
flipy = masterSettings.flipy ; 
timeInterval = masterSettings.timeInterval ;  % physical interval between timepoints
timeUnits = masterSettings.timeUnits ; % physical unit of time between timepoints
spaceUnits = masterSettings.spaceUnits ; % unit of distance of full resolution data pixels ('$\mu$m')
scale = masterSettings.scale ;      % scale for conversion to 16 bit
file32Base = masterSettings.file32Base ; 
fn = masterSettings.fn ;
fn_prestab = masterSettings.fn_prestab ;
set_preilastikaxisorder = masterSettings.set_preilastikaxisorder ;
swapZT = masterSettings.swapZT ;
t0_for_phi0 = masterSettings.t0_for_phi0 ;
nU = masterSettings.nU ;
nV = masterSettings.nV ;


% Fill in
scale = masterSettings.scale ;      % scale for conversion to 16 bit
file32Base = masterSettings.file32Base ; 
fn = masterSettings.fn ;
fn_prestab = masterSettings.fn_prestab ;
set_preilastikaxisorder = masterSettings.set_preilastikaxisorder ;

dir32bit = fullfile(dataDir, 'deconvolved_32bit') ;
dir16bit = fullfile(dataDir, 'deconvolved_16bit') ;
dir16bit_prestab = fullfile(dir16bit, 'data_pre_stabilization') ;

% I. INITIALIZE ImSAnE PROJECT ===========================================
% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.
cd(dir16bit)
dataDir = cd ;
projectDir = dataDir ;
% [ projectDir, ~, ~ ] = fileparts(matlab.desktop.editor.getActiveFilename); 
cd(projectDir);
if projectDir(end) ~= '/'
    projectDir = [projectDir '/'];
end

% Start by creating an experiment object, optionally pass on the project
% directory (otherwise it will ask), and change into the directory of the
% data.  This serves as a front-end for data loading, detection, fitting
% etc.
xp = struct() ;
% A filename base template - to be used throughout this script
fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fn, '.tif'];
fileMeta.nChannels          = nChannels;
fileMeta.timePoints         = timePoints ;
fileMeta.stackResolution    = stackResolution;
fileMeta.swapZT             = 1;

% first_tp is also required, which sets the tp to do individually.
first_tp = 1 ;
expMeta                     = struct();
expMeta.channelsUsed        = channelsUsed ;
expMeta.channelColor        = 1;
expMeta.description         = 'Drosophila gut';
expMeta.dynamicSurface      = 1;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.fitTime             = fileMeta.timePoints(first_tp);
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';

% Now set the meta data in the experiment.
xp.fileMeta = fileMeta;
xp.expMeta = expMeta;
xp.tIdx(timePoints) = 1:length(timePoints) ;

clear fileMeta expMeta
% % Set detect options ------------------------------------------------------
% xp.setDetectOptions( detectOptions );
% disp('done')
xp.detectOptions = struct() ;
xp.detectOptions.ssfactor = 4 ;
xp.detector.options.ofn_smoothply = 'mesh_%06d' ;

% QS DEFINITION
opts.meshDir = meshDir ;
opts.flipy = flipy ;
opts.timeInterval = timeInterval ;
opts.timeUnits = timeUnits ;
opts.spaceUnits = '$\mu$m' ;
opts.nV = nV ;
opts.nU = nU ;
opts.normalShift = 10 ;
opts.a_fixed = 2.0 ;
opts.adjustlow = 1.00 ;                  %  floor for intensity adjustment
opts.adjusthigh = 99.9 ;                 % ceil for intensity adjustment (clip)
% opts.adjustlow = 0 ;                  %  floor for intensity adjustment
% opts.adjusthigh = 0 ;                 % ceil for intensity adjustment (clip)
opts.phiMethod = 'curves3d' ;
opts.lambda_mesh = 0 ;
opts.lambda = 0 ;
opts.lambda_err = 0 ;
%  opts.lambda_mesh = 0.002 ;
%  opts.lambda = 0.01 ;
%  opts.lambda_err = 0.01 ;
 
disp('defining QS')

% build needed files
QS = QuapSlap(xp, opts) ;
disp('done')

%% Now use QS to find pathline xy and xyz
cellVertexPathlineFn = fullfile(QS.dir.segmentation, 'pathlines', ...
    sprintf('cellVertexPathlines_%06dt0.mat', QS.t0set())) ;
tmp = load(cellVertexPathlineFn, 'segVertexPathlines2D', ...
         'segVertexPathlines3D', 'cellIDs') ;
segP2 = segVertexPathlines2D ;
segP3 = segVertexPathlines3D ;

% find position where folds occur in pullback space
tref = 126 ;
QS.setTime(tref)
tmp =  QS.getCurrentSegmentation2DCorrected ;
xy0 = tmp.seg2d.cdat.centroid ;
[xpath, ypath] = QS.pullbackPathlines(xy0(:, 1), xy0(:,2), tref) ;

cellCntrdPathlineFn = fullfile(QS.dir.segmentation, 'pathlines', ...
    sprintf('cellCentroidPathlines_%06dt0.mat', tref)) ;
segCentroidPathlines2D = cat(3, xpath, ypath) ;
save(cellCntrdPathlineFn, 'segCentroidPathlines2D')

%% Get fold positions 
tref = 126 ;
QS.setTime(tref)
features = QS.getFeatures() ;
folds = features.folds ;
fold0 = folds(tref, :) ;
% fold0 is: 33, 57, 78




%% Aspect ratio and anisotropy along uv coordinates
% Plot each timepoint
close all
tps = [126, 166, 206] ;
% tps = 96:10:206 ;
Ntps = length(tps) ;
sz = 10 ;
alph = 0.05;
ylims1 = [-1,4] ;
ylims2 = [1,5] ;
psat = 0.3 ; 

fig1 = figure('units', 'centimeters', 'position', [0, 0, 12, 5]);
fig2 = figure('units', 'centimeters', 'position', [0, 0, 12, 5]);
colors = colormap(0.9 * cmap('I1', 'N', length(tps))) ;
for qq = 1:length(tps)
    
    % MAC:
    % fn = fullfile(fns(tidx).folder, fns(tidx).name) ;
    % seg = load(fn) ;
    
    
    tp = tps(qq) ;
    QS.setTime(tp) ;
    seg = QS.getCurrentSegmentation3DCorrected() ;
    
    % note positions of folds
    fq = folds(QS.xp.tIdx(tp), :) ;
    
    seg = seg.seg3d ;
    ar =  sqrt(seg.qualities.moment2 ./ seg.qualities.moment1) ;
    qcos2t = (1-ar) .* cos(seg.qualities.ang1*2) ;
    
    % ap binned with weights
    stat = seg.statistics ;
    weights = stat.weights ;
    qcos2t(weights==0) = 0 ;
    ap = stat.apStats.aspectWeighted.apBins ;
    c2t = -stat.apStats.aspectWeighted.apCos2Theta ;
    cs = stat.apStats.aspectWeighted.apCos2ThetaStd ;
    ce = stat.apStats.aspectWeighted.apCos2ThetaSte ;
        
    % Plot position along ap axis and oriented anitostropy
    try
        set(0, 'CurrentFigure', fig1)
    catch
        figure()
    end
    
    subplot(1,Ntps, qq)
    xx = seg.cdat.centroids_uv(:, 1)  ;
    scatter(xx, qcos2t, sz, 'filled', 'markerfacecolor', colors(qq, :), ...
        'markeredgecolor', 'none', 'markerfacealpha', alph)
    hold on;
    for ff = fq
        plot(double(ff)/QS.nU * [1,1], ylims1, 'k--' ) ;
    end
    ylim(ylims1)
    
    allData = struct() ;
    allData.scatter = struct('xx', xx, ...
        'qcos2t', qcos2t,'ar',ar, 'weights', weights, 'theta', seg.qualities.ang1) ;
    
    %% Plot the anisotropy
    lineProps = {'-','color', colors(qq, :)} ;
    % shadedErrorBar(ap, c2t, cs, 'lineProps', lineProps)
    hold on;
    shadedErrorBar(ap, c2t, cs, 'lineProps', lineProps, 'patchSaturation', psat)
    shadedErrorBar(ap, c2t, ce, 'lineProps', lineProps, 'patchSaturation', 1)
    xlabel('ap position [$s/L$]', 'interpreter', 'latex')
    ylabel('cell anisotropy $(1-a/b)\cos 2 \theta$', 'interpreter', 'latex')
    
    allData.shadedBar = struct('ap', ap, ...
        'c2t', c2t, ...
        'cs', cs, ...
        'ce', ce) ;
    
    %% Get change in anisotropy
    % ap with bins (redo for checking)
    xedges = sort([0, ap, 1]) ;
    % remove edges within fold bounds
    % e2rm = [] ;
    % for f0Id = 1:3
    %     e2rm = [e2rm, find(xedges > double(fold0(f0Id))/QS.nU-0.05+eps & ...
    %         xedges < double(fold0(f0Id))/QS.nU +0.05-eps)] ;
    %  end   
    % e2keep = setdiff(1:length(xedges), e2rm) ;
    % xedges = xedges(e2keep) ;
    % xedges = sort([xedges,...
    %     double(fold0) ./ QS.nU - 0.05, ...
    %     double(fold0) ./ QS.nU + 0.05]) ;

    yy = seg.cdat.centroids_uv(:, 2)  ;
    tmpOpts = struct('timePoints', tref:tps(qq)) ;
    disp('projecting pathlines back onto original')
    xypix = QS.uv2XY([2000,2000], [xx, yy], true, 1, 1);
    [XX,YY] = QS.pullbackPathlines(xypix(:, 1), xypix(:,2),...
        tps(qq), tmpOpts) ;
    
    % check pullback pathlines
    % for pp = 1:size(XX,1)
    %    scatter(XX(pp, :), YY(pp, :), '.')
    %    pause(0.001)
    % end
    
    x0 = XX(1, :)/2000 ;
    y0 = YY(1, :)/2000 ;
        
    [mid_ap, mean_QAc2t_ap, std_QAc2t_ap, ~, ste_QAc2t_ap] = ...
        binDataMeanStdWeighted(x0, qcos2t, ...
            xedges, ones(size(xx))) ;
    [mid_ap, mean_QAc2t_apw, std_QAc2t_apw, ~, ste_QAc2t_apw] = ...
        binDataMeanStdWeighted(x0, qcos2t, ...
            xedges, weights) ;

    % lineProps = {'-','color', colors(qq, :)} ;
    % shadedErrorBar(mid_ap, mean_QAc2t_ap, ste_Aspect_ap, 'lineProps', lineProps, ...
    %     'patchSaturation', psat)
    % shadedErrorBar(mid_ap, mean_QAc2t_ap, std_Aspect_ap, 'lineProps', lineProps, ...
    %     'patchSaturation', 1)
    % lineProps = {'-','color', 'k'} ;
    % shadedErrorBar(mid_ap, mean_QAc2t_apw, ste_Aspect_apw, 'lineProps', lineProps)
    % shadedErrorBar(mid_ap, mean_QAc2t_apw, std_Aspect_apw, 'lineProps', lineProps)
    
    allData.lagrangianParameterizationMeasures = struct('x0', x0, 'y0', y0,...
        'mid_ap', mid_ap, 'mean_QAc2t_ap', mean_QAc2t_ap, ...
        'std_QAc2t_ap', std_QAc2t_ap, ...
        'ste_QAc2t_ap', ste_QAc2t_ap, ...
        'mean_QAc2t_apw', mean_QAc2t_apw, ...
        'std_QAc2t_apw', std_QAc2t_apw, ...
        'ste_QAc2t_apw', ste_QAc2t_apw) ;
    
    
    %% Plot position along ap axis and aspect ratio (raw) anitostropy
    try
        set(0, 'CurrentFigure', fig2)
    catch
        figure()
    end
    subplot(1, Ntps, qq)
    scatter(xx, ar, sz, 'filled', 'markerfacecolor', colors(qq, :), ...
        'markeredgecolor', 'none', 'markerfacealpha', alph)
    hold on;
    arQ = sqrt(stat.apStats.aspectWeighted.apCos2Theta.^2 + stat.apStats.aspectWeighted.apCos2Theta.^2) ;
    arm = stat.apStats.mean_Aspect_ap ;
    ars = stat.apStats.std_Aspect_ap ;
    are = stat.apStats.ste_Aspect_ap ;
    lineProps = {'-','color', colors(qq, :)} ;
    shadedErrorBar(ap, arm, ars, 'lineProps', lineProps)
    shadedErrorBar(ap, arm, are, 'lineProps', lineProps)

    check_debug = false ;
    if check_debug
        % ap with bins (redo for checking)
        weights = stat.weights ;
        xedges = [0, ap, 1] ;
        ar(weights==0) = 0 ;
        [mid_ap, mean_Aspect_ap, std_Aspect_ap, ~, ste_Aspect_ap] = ...
            binDataMeanStdWeighted(xx, ar, ...
                xedges, ones(size(xx))) ;
        [mid_ap, mean_Aspect_apw, std_Aspect_apw, ~, ste_Aspect_apw] = ...
            binDataMeanStdWeighted(xx, ar, ...
                xedges, weights) ;

        lineProps = {'-','color', colors(qq, :)} ;
        shadedErrorBar(mid_ap, mean_Aspect_ap, ste_Aspect_ap, 'lineProps', lineProps)
        shadedErrorBar(mid_ap, mean_Aspect_ap, std_Aspect_ap, 'lineProps', lineProps)

        % lineProps = {'-','color', 'k'} ;
        % shadedErrorBar(mid_ap, mean_Aspect_apw, ste_Aspect_apw, 'lineProps', lineProps)
        % shadedErrorBar(mid_ap, mean_Aspect_apw, std_Aspect_apw, 'lineProps', lineProps)
    end
    
    for ff = fq
        plot(double(ff)/QS.nU * [1,1], ylims2, 'k--' ) ;
    end
    ylim(ylims2)
    xlabel('ap position [$s/L$]', 'interpreter', 'latex')
    ylabel('cell aspect ratio, $a/b$', 'interpreter', 'latex')
    
    allDatas{qq} = allData ;
end
    
set(0, 'CurrentFigure', fig1)
fn1 = fullfile(QS.dir.segmentation, ['aspect_ratio_over_time_' num2str(tps(1)) '_' ...
    num2str(tps(end)) '_Ntps' num2str(Ntps) '_eLifeResub_00.pdf']) ;
saveas(fig1, fn1)

set(0, 'CurrentFigure', fig2)
fn2 = fullfile(QS.dir.segmentation, ['aspect_ratio_over_time_' num2str(tps(1)) '_' ...
    num2str(tps(end)) '_Ntps' num2str(Ntps) '_eLifeResub_01.pdf']) ;
saveas(fig2, fn2)

fn = fullfile(QS.dir.segmentation, 'eLife_resubmission_plot.mat') ;
save(fn, 'allDatas')


%% Plot change over time
close all
mean0 = allDatas{1}.lagrangianParameterizationMeasures.mean_QAc2t_apw ;
mean1 = allDatas{end}.lagrangianParameterizationMeasures.mean_QAc2t_apw ;
std0 = allDatas{1}.lagrangianParameterizationMeasures.std_QAc2t_apw ;
std1 = allDatas{end}.lagrangianParameterizationMeasures.std_QAc2t_apw ;
ste0 = allDatas{1}.lagrangianParameterizationMeasures.ste_QAc2t_apw ;
ste1 = allDatas{end}.lagrangianParameterizationMeasures.ste_QAc2t_apw ;
stdD = sqrt(std0.^2 + std1.^2) ;
steD = sqrt(ste0.^2 + ste1.^2) ;

fig1 = figure('units', 'centimeters', 'position', [0, 0, 12, 5]);
subplot(1, Ntps, round(Ntps*0.5))
lineProps = {'-','color', 'k'} ;
hold on;
shadedErrorBar(mid_ap, mean1-mean0, stdD, 'lineProps', lineProps, 'patchSaturation', 0.1)
shadedErrorBar(mid_ap, mean1-mean0, steD, 'lineProps', lineProps, 'patchSaturation', psat)
ylims = [-3,0.5] ;
hold on;
    for ff = fold0
        plot(double(ff)/QS.nU * [1,1], ylims, 'k--' ) ;
    end
    xlabel('ap position at $t=0$, [$s/L$]', 'interpreter', 'latex')
    ylabel('change in cell anisotropy, $\Delta$', 'interpreter', 'latex')
    ylim(ylims)
    
fn2 = fullfile(QS.dir.segmentation, ['aspect_ratio_over_time_' num2str(tps(1)) '_' ...
    num2str(tps(end)) '_Ntps' num2str(Ntps) '_eLifeResub_Change.pdf']) ;
saveas(gcf, fn2)

 
ylim([-2, 0])
fn2 = fullfile(QS.dir.segmentation, ['aspect_ratio_over_time_' num2str(tps(1)) '_' ...
    num2str(tps(end)) '_Ntps' num2str(Ntps) '_eLifeResub_Change_zoom.pdf']) ;
saveas(gcf, fn2)

 