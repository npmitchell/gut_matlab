
%% Clear workspace ========================================================
% We start by clearing the memory and closing all figures
clear; close all; clc;
% change this path, for convenience
% cd /mnt/crunch/48Ygal4-UAShistRFP/201904031830_great/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2/data/
% cd /mnt/crunch/48YGal4UasLifeActRuby/201904021800_great/Time6views_60sec_1p4um_25x_1p0mW_exp0p150_3/data/
% cd /mnt/data/48YGal4UasLifeActRuby/201902201200_unusualfolds/Time6views_60sec_1p4um_25x_obis1_exp0p35_3/data/
% cd /mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data
% cd /mnt/crunch/48Ygal4UASCAAXmCherry/201903211930_great/Time6views_60sec_1p4um_25x_1p0mW_exp0p150/data/
% cd /mnt/crunch/48Ygal4UASsqhGFP/201902271940_excellent_notunpacked/Time6views_60sec_1p4um_25x_1p5mW_exp1p0_3/data/
% cd /mnt/data/mef2GAL4klarUASCAAXmChHiFP/202003151700_1p4um_0p5ms3msexp/data/
% cd /mnt/data/mef2GAL4klarUASCAAXmChHiFP/202003151700_1p4um_0p5ms3msexp/data/


% .=========.
% |  VIP10  |
% .=========.
% cd /mnt/crunch/gut/48YGal4UasLifeActRuby/201907311600_48YGal4UasLifeActRuby_60s_exp0p150_1p0mW_25x_1p4um
% cd /mnt/crunch/gut/48YGal4klarUASCAAXmChHiFP/202001221000_60sec_1p4um_25x_1mW_2mW_exp0p25_exp0p7/Time3views_1017/data/
% cd /mnt/crunch/gut/Mef2Gal4klarUASCAAXmChHiFP/202003151700_1p4um_0p5ms3msexp/Time3views_1/data/
% cd /mnt/crunch/gut/Mef2Gal4klarUASCAAXmChHiFP/202007151930_1p4um_0p5msexp/Time3views_25x_60s/data/
% cd /mnt/crunch/gut/handGAL4klarHandGFPhistGFP/202105072030_1p4um_0p2ms_1mWGFP/view4_1p4um_25x_1mW_exp0p2ms_120spf/data ;

% .============.
% |  Home Mac  |
% .============.
cd /Users/npmitchell/Desktop/HandGFP_celltracking/handGAL4klarHandGFPhistGFP/


codepath = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/' ;
ms_scriptDir = fullfile(codepath, 'morphsnakes_wrapper') ;


dataDir = cd ;

%% PATHS ==================================================================
origpath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(origpath))
aux_paths_and_colors
cd(dataDir)

%% Global options
% Decide whether to change previously stored detection Opts, if they exist
overwrite_masterSettings = false ;
overwrite_mips = false ;
overwrite_detOpts = false ;
run_full_dataset_ms = false ;
overwrite_alignAPDVOpts = false ;
overwrite_APDVCOMs = false ;
overwrite_APDVMeshAlignment = false ;
overwrite_alignedMeshIms = false ;
overwrite_centerlines = false ;
overwrite_centerlineIms = false ;
overwrite_TextureMeshOpts = false ;
overwrite_endcapOpts = false ;
overwrite_idAnomClines = false ;
overwrite_cleanCylMesh = false ;

%% DEFINE NEW MASTER SETTINGS

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
    
dir32bit = fullfile(dataDir, 'deconvolved_32bit') ;
dir16bit = fullfile(dataDir, 'deconvolved_16bit') ;
dir16bit_prestab = fullfile(dir16bit, 'data_pre_stabilization') ;

%% END OF EXPERIMENT METADATA =============================================
% =========================================================================
% =========================================================================
cd(dir16bit)

%% I. INITIALIZE ImSAnE PROJECT ===========================================
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
xp = project.Experiment(projectDir, dataDir);

% Set file and experiment meta data
%
% We assume on individual image stack for each time point, labeled by time.
%  To be able to load the stack, we need to tell the project wehre the data
%  is, what convention is assumed for the file names, available time
%  points, and the stack resolution.  Options for modules in ImSAnE are
%  organized in MATLAB structures, i.e a pair of field names and values are
%  provided for each option.
%
% The following file metadata information is required:
% * 'directory'         , the project directory (full path)
% * 'dataDir'           , the data directory (full path)
% * 'filenameFormat'    , fprintf type format spec of file name
% * 'timePoints'        , list of itmes available stored as a vector
% * 'stackResolution'   , stack resolution in microns, e.g. [0.25 0.25 1]
%
% The following file metadata information is optional:
% * 'imageSpace'        , bit depth of image, such as uint16 etc., defined
%                         in Stack class
% * 'stackSize'         , size of stack in pixels per dimension 
%                         [xSize ySize zSize]
% * 'swapZT'            , set=1 if time is 3rd dimension and z is 4th

% A filename base template - to be used throughout this script
fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fn, '.tif'];
fileMeta.nChannels          = nChannels;
fileMeta.timePoints         = timePoints ;
fileMeta.stackResolution    = stackResolution;
fileMeta.swapZT             = masterSettings.swapZT;

% Set required additional information on the experiment. A verbal data set
% description, Jitter correct by translating  the sample, which time point
% to use for fitting, etc.
%
% The following project metadata information is required:
% * 'channelsUsed'      , the channels used, e.g. [1 3] for RGB
% * 'channelColor'      , mapping from element in channels used to RGB = 123
% * 'dynamicSurface'    , Not implemented yet, future plan: boolean, false: static surface
% * 'detectorType'      , name of detector class, e.g. radielEdgeDetector
%                         ,(user threshholded), fastCylinderDetector
% * 'fitterType'        , name of fitter class
%
% The following project meta data information is optional:
% * 'description'     , string describing the data set set experiments metadata, 
%                                such as a description, and if the surface is dynamic,
%                                or requires drift correction of the sample.
% * 'jitterCorrection', Boolean, false: No fft based jitter correction 

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
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();

clear fileMeta expMeta

%% LOAD THE FIRST TIME POINT ==============================================
xp.setTime(xp.fileMeta.timePoints(1)) ;
% xp.loadTime(xp.fileMeta.timePoints(first_tp));
% xp.rescaleStackToUnitAspect();

%% SET DETECT OPTIONS =====================================================
% Must run this section for later functionality.
% Mesh extraction options
if run_full_dataset_ms
    run_full_dataset = projectDir ; 
else
    run_full_dataset = 'none' ;
end

% Load/define the surface detection parameters
msls_detOpts_fn = fullfile(projectDir, 'msls_detectOpts.mat') ;
if exist(msls_detOpts_fn, 'file') && ~overwrite_detOpts
    load(msls_detOpts_fn, 'detectOptions')
else
    channel = 1;
    foreGroundChannel = 1;
    % ssfactor = 4;
    % niter = 15 ;
    % niter0 = 15 ;
    % ofn_ply = 'mesh_' ; % mesh_apical_ms_stab_' ; 
    % ofn_ls = 'msls_' ; % mesh_apical_stab_' ;
    ofn_smoothply = 'mesh_apical_stab_' ; % mesh_apical_stab_' ;
    % lambda1 = 1 ;
    % lambda2 = 1 ;
    % exit_thres = 0.0001 ;
    % if contains(projectDir, 'LifeActRuby') || contains(projectDir, 'CAAXmCherry') 
    %     nu = 4 ;  % For volumetric
    %     smoothing = 0.20 ;
    % else
    %     nu = 0.00 ;
    %     smoothing = 0.10 ;
    % end
    % pre_nu = -4 ;
    % pre_smoothing = 0 ;
    % post_nu = 3 ;
    % post_smoothing = 1 ;
    zdim = 2 ;
    
    % for caax_great
    ssfactor=4
    niter=35
    niter0=35
    ofn_ply='mesh_ms_'
    ofn_ls='msls_'
    ofn_smoothply = 'mesh_stab_' ; % mesh_apical_stab_' ;
    ms_scriptDir='/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/'
    pre_nu=-5
    pre_smoothing=0
    lambda1=1
    lambda2=1
    exit_thres=0.000005
    smoothing=0.20
    nu=0
    post_nu=2
    post_smoothing=5
    
    init_ls_fn = 'msls_initguess.h5';
    % mlxprogram = fullfile(meshlabCodeDir, ...
    %     'laplace_refine_HCLaplace_LaplaceSPreserve_QuadEdgeCollapse60kfaces.mlx') ;
    mlxprogram = fullfile(meshlabCodeDir, ...
        'laplace_surface_rm_resample30k_reconstruct_LS3_1p2pc_ssfactor4.mlx') ;
    radius_guess = 40 ;
    center_guess = '200,75,75' ;
    dtype = 'h5' ;
    mask = 'none' ;
    prob_searchstr = '_stab_Probabilities.h5' ;
    preilastikaxisorder= set_preilastikaxisorder; ... % axis order in input to ilastik as h5s. To keep as saved coords use xyzc
    ilastikaxisorder= 'cxyz'; ... % axis order as output by ilastik probabilities h5
    imsaneaxisorder = 'xyzc'; ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc
    include_boundary_faces = true ;
    smooth_with_matlab = -1;
    
    % Name the output mesh directory ------------------------------------------
    % msls_exten = ['_prnu' strrep(strrep(num2str(pre_nu, '%d'), '.', 'p'), '-', 'n')];
    % msls_exten = [msls_exten '_prs' strrep(num2str(pre_smoothing, '%d'), '.', 'p') ];
    % msls_exten = [msls_exten '_nu' strrep(num2str(nu, '%0.2f'), '.', 'p') ];
    % msls_exten = [msls_exten '_s' strrep(num2str(smoothing, '%0.2f'), '.', 'p') ];
    % msls_exten = [msls_exten '_pn' num2str(post_nu, '%d') '_ps',...
    %     num2str(post_smoothing)];
    % msls_exten = [msls_exten '_l' num2str(lambda1) '_l' num2str(lambda2) ];
    meshDir = [projectDir 'msls_output' filesep];
    % meshDir = [meshDir msls_exten '/'] ;

    % Surface detection parameters --------------------------------------------
    detectOptions = struct( 'channel', channel, ...
        'ssfactor', ssfactor, ...
        'niter', niter,...
        'niter0', niter0, ...
        'lambda1', lambda1, ...
        'lambda2', lambda2, ...
        'nu', nu, ...
        'smoothing', smoothing, ...
        'post_nu', post_nu, ...
        'post_smoothing', post_smoothing, ...
        'exit_thres', exit_thres, ...
        'foreGroundChannel', foreGroundChannel, ...
        'fileName', sprintf( fn, xp.currentTime ), ...
        'mslsDir', meshDir, ...
        'ofn_ls', ofn_ls, ...
        'ofn_ply', ofn_ply,...
        'ms_scriptDir', ms_scriptDir, ...
        'timepoint', xp.currentTime, ...
        'zdim', zdim, ...
        'pre_nu', pre_nu, ...
        'pre_smoothing', pre_smoothing, ...
        'ofn_smoothply', ofn_smoothply, ...
        'mlxprogram', mlxprogram, ...
        'init_ls_fn', init_ls_fn, ... % set to none to load prev tp
        'run_full_dataset', run_full_dataset,... % projectDir, ... % set to 'none' for single tp
        'radius_guess', radius_guess, ...
        'dset_name', 'exported_data',...
        'save', true, ... % whether to save images of debugging output
        'center_guess', center_guess,... % xyz of the initial guess sphere ;
        'plot_mesh3d', false, ...
        'dtype', dtype,...
        'mask', mask,...
        'mesh_from_pointcloud', false, ...
        'prob_searchstr', prob_searchstr, ...
        'physicalaxisorder', imsaneaxisorder, ... 
        'preilastikaxisorder', preilastikaxisorder, ...
        'ilastikaxisorder', ilastikaxisorder, ... 
        'include_boundary_faces', include_boundary_faces, ...
        'smooth_with_matlab', smooth_with_matlab, ...
        'pythonVersion', ''); % version of python to call = '2' or '3', as string

    % save options
    if exist(msls_detOpts_fn, 'file')
        disp('Overwriting detectOptions --> renaming existing as backup')
        backupfn1 = [msls_detOpts_fn '_backup1'] ;
        if exist(backupfn1, 'file')
            backupfn2 = [msls_detOpts_fn '_backup2'] ; 
            system(['mv ' backupfn1 ' ' backupfn2])
        end
        system(['mv ' msls_detOpts_fn ' ' backupfn1])
    end
    disp('Saving detect Options to disk')
    save(msls_detOpts_fn, 'detectOptions') ;
end

% Overwrite certain parameters for script structure
detectOptions.fileName = sprintf( fn, xp.currentTime ) ;
detectOptions.run_full_dataset = run_full_dataset ;
detectOptions.ms_scriptDir = ms_scriptDir ;
meshDir = detectOptions.mslsDir ;
% These are now in QS.
% meshFileBase = [ofn_smoothply '%06d'] ;
% alignedMeshBase = [ofn_smoothply '%06d_APDV_um'] ;

% Set detect options ------------------------------------------------------
xp.setDetectOptions( detectOptions );
disp('done')

%% Define QuapSlap object
nU = masterSettings.nU ;
nV = masterSettings.nV ;
opts.meshDir = meshDir ;
opts.flipy = flipy ;
opts.timeInterval = timeInterval ;
opts.timeUnits = timeUnits ;
opts.spaceUnits = spaceUnits ;
opts.nU = nU ;
opts.nV = nV ;
opts.normalShift = 10 ;
opts.a_fixed = 2.0 ;
opts.adjustlow = 1.00 ;         % floor for intensity adjustment
opts.adjusthigh = 99.9 ;        % ceil for intensity adjustment (clip)
opts.phiMethod = 'curves3d' ;
opts.lambda_mesh = 0.002 ;
opts.lambda = 0.01 ;
opts.lambda_err = 0.01 ;
disp('defining QS')
QS = QuapSlap(xp, opts) ;
disp('done')

%% Inspect learning+tracking ground truth output

subdir = 'endoderm_normalShift05_p09_n00_s0p75_lambda0p0002'; 
outDir = fullfile(QS.dir.tracking, ...
    subdir, 'images_groundTruthTracking') ;
tracksfn = fullfile(QS.dir.im_sp_sm, subdir, 'tracks', ...
    '%05d.h5') ;
outGraphFn = fullfile(QS.dir.im_sp_sm, subdir, ...
    'tracks', 'graph0_tracking.mat') ;
rawImFileBase = fullfile(QS.dir.im_sp_sm, subdir,'endoderm_imagestack_LUT', ...
    'Time_%06d_c1_stab_pbspsm_LUT.tif') ; 

if ~exist(outGraphFn, 'file')
    Options = struct() ;
    Options.maxN_for_plot = 5000 ;
    Options.allowSplitting = false ;
    Options.allowMerging = false ;
    GG = unpackManualIlastikGroundTruthH5(tracksfn, timePoints, outDir, ...
        rawImFileBase, Options) ;
    disp(['Saving ' outGraphFn])
    save(outGraphFn, 'GG')
else
    load(outGraphFn, 'GG')
end

%% Check quality of GG
saveDir = fullfile(outDir, 'graph_result') ;
Options = struct() ;
Options.plot_generations = false ;
Options.faceAlpha = 0.6 ;
plotGraphTrackingVoronoi(GG, rawImFileBase, saveDir, Options)
saveDir = fullfile(outDir, 'graph_result_scatter') ;
plotGraphTrackingScatter(GG, rawImFileBase, saveDir, Options)

%% Manual acquisition for muscle tracks
subdir = 'muscle_normalShiftn10_p05_n50_s1p00_lambda0p005_maxProj';
imDir = fullfile(QS.dir.im_sp_sm, subdir, 'muscle_imagestack_LUT') ;
trackOutfn = fullfile(QS.dir.tracking, 'muscle', 'muscle_tracks.mat') ;

% codeDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/';
% addpath(fullfile(codeDir, 'ParhyaleCellTracker/external/subtightplot/'))
% addpath(fullfile(codeDir, 'gut_matlab/addpath_recurse/'))
% addpath_recurse('/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/gut_matlab/')
% subdir = 'muscle_normalShiftn10_p05_n50_s1p00_lambda0p005_maxProj';
% imDir = fullfile('./', subdir, 'muscle_imagestack_LUT') ;
% trackOutfn = fullfile('./muscle_tracks.mat') ;

load(fullfile(QS.dir.tracking, 'muscle_tracks.mat'), 'tracks')
timePoints = 1:60 ;
nTracks = 300 ;
fileBase = fullfile(imDir, 'Time_%06d_c1_stab_pbspsm_LUT.tif') ;
manualTrack2D(tracks, fileBase, timePoints, trackOutfn, nTracks) ;


%% Muscle add tracks
subdir = 'muscle';
imDir = fullfile(QS.dir.im_sp_sm, subdir) ;
trackOutfn = fullfile(QS.dir.tracking, 'muscle_tracks.mat') ;
load(fullfile(QS.dir.tracking, 'muscle_tracks.mat'), 'tracks')
timePoints = 1:60 ;
nTracks = 300 ;
fileBase = fullfile(imDir, 'Time_%06d_c1_stab_pbspsm.tif') ;
manualTrack2D(tracks, fileBase, timePoints, trackOutfn, nTracks) ;

%% Endoderm add tracks

trackOutfn = fullfile(QS.dir.tracking, 'endoderm_tracks_correction.mat') ;
load(trackOutfn, 'tracks')
timePoints = 1:60 ;
imDir = fullfile(QS.dir.im_sp_sm, 'endoderm') ;
fileBase = fullfile(imDir, 'Time_%06d_c1_stab_pbspsm.tif') ;

tracks2Add = length(tracks)+1:length(tracks) + 10 ;
tracks = manualTrack2D(tracks, fileBase, timePoints, trackOutfn, tracks2Add) ;

%% Open manual tracking gui
addpath_recurse('/mnt/data/code/ParhyaleCellTracker/')
close all
Options = struct() ;
Options.allowOverlaps = true  ;
Options.drawGenerations = false ;
Options.colorByLineage = true ;
Options.infoIndex = 1 ;
[Gout, divStruct] = parhyale_master_gui(GG, rawImFileBase, Options) ;

%% Convert digraph to cell arrays of tracks
outTrackFn = fullfile(QS.dir.tracking, 'endoderm_tracks.mat') ;
tracks = trackingGraph2Cell(GG, timePoints) ;
save(outGraphFn, 'tracks', 'GG')

%% Manual correction of endoderm tracks
subdir = 'endoderm' ;
imDir = fullfile(QS.dir.im_sp_sm, subdir) ;
timePoints = 1:60 ;
fileBase = fullfile(imDir, QS.fileBase.im_sp_sm) ;
trackOutfn = fullfile(QS.dir.tracking, 'endoderm_tracks_correction.mat') ;
load(trackOutfn, 'tracks') ;
tracks2Add = 590 ;
[newTracks, newG] = manualCorrectTracks2D(tracks, fileBase, timePoints, trackOutfn, tracks2Add) ;


%% Instant replay of tracks
if recap
    clf
    set( fig, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    times2play = timePoints ;
    for kk = 1:length(times2play)
        ttemp = times2play(kk) ;
        im1 = imread(sprintf(fileBase, ttemp));
        imshow(im1) ;
        hold on;

        for iprev = 1:ii-1
            trackprev = currentTracks{iprev} ;
            plot(trackprev(kk, 1), trackprev(kk, 2), 'o', 'color', ...
                bluecolor, 'markerSize', markerSize, 'lineWidth', lwidth)
        end

        tracknow = currentTracks{ii} ;
        plot(tracknow(kk, 1), tracknow(kk, 2), 'o', 'color', orange, ...
            'markerSize', markerSize, 'lineWidth', lwidth)

        title(['playing tracks 1-' num2str(ii) ': t=' num2str(ttemp)])

        pause(0.01 * pausetime)
    end
end

%% Relative motion
Options = struct() ;
Options.max_dv = Inf ;
Options.max_du = 50 ;
subdir1 = 'endoderm' ;
subdir2 = 'muscle' ;
track1fn = fullfile(QS.dir.tracking, 'endoderm_tracks_correction.mat') ;
track2fn = fullfile(QS.dir.tracking, 'muscle_tracks.mat') ;
QS.measureRelativeMotionTracks(track1fn, track2fn)


function measureRelativeMotionTracks(QS, track1fn, track2fn)
%
% Options : struct with fields (optional)
%   timePoints : nTimePoints x 1 numeric array
%       timepoints of the tracks, so that track{ii}(tidx, :) is the
%       XY position of track ii at time timePoints(tidx). Note this cannot 
%       be some arbitrary subset of the experiment timepoints
%       QS.xp.fileMeta.timePoints but must instead match the true track 
%       info!
%   max_du : float in units of QS.spaceUnits (optional, default=Inf)
%       maximum possible physical displacement between nuclei, for 
%       filtering out unreasonable results
%
    relMotionFn = fullfile(QS.dir.tracking, 'relative_motion_tracks.mat') ;
    relMotionStartFn = fullfile(QS.dir.tracking, 'relative_motion_trackStart.png') ;

    tmp1 = load(track1fn) ;
    tmp2 = load(track2fn) ;

    if ~isfield(tmp1, 'tracks') && isfield(tmp1, 'GG')
        tracks1 = trackingGraph2Cell(tmp1.GG) ;
    else
        tracks1 = tmp1.tracks ;
    end
    if ~isfield(tmp2, 'tracks') && isfield(tmp2, 'GG')
        tracks2 = trackingGraph2Cell(tmp2.GG) ;
    else
        tracks2 = tmp2.tracks ;
    end

    n1 = length(tracks1) ;
    n2 = length(tracks2) ;

    % mode for distance
    distMode = 'geodesic' ; %'geodesics' 'fastEuclidean' ;
    preview = false ;

    % Options handling
    doubleCovered = false ;
    umax = 1.0 ;
    vmax = 1.0 ;
    coordSys = 'spsmrs' ;
    max_du = Inf ;
    timePoints = QS.xp.fileMeta.timePoints ;
    if isfield(Options, 'doubleCovered')
        doubleCovered = Options.doubleCovered ;
    end
    if isfield(Options, 'umax')
        umax = Options.umax ;
    end
    if isfield(Options, 'vmax')
        vmax = Options.vmax ;
    end
    if isfield(Options, 'max_du')
        max_du = Options.max_du ;
    end
    if isfield(Options, 'coordSys')
        coordSys = Options.coordSys ;
    end
    if isfield(Options, 'timePoints')
        timePoints = Options.timePoints ;
    end

    % Colormap and shorthand
    nTracks = length(tracks2) ;
    colors = jetshuffle(nTracks) ;
    nTimePoints = length(timePoints) ;

    % Pre-allocate positions of starting points in tracks of tracks1 (endoderm)
    nearby = zeros(n1, 2) ;
    for trackID = 1:n1
        nearby(trackID, :) = tracks1{trackID}(1, 1:2) ;
    end

    overwrite = false ;
    if exist(relMotionFn, 'file') && ~overwrite
        load(relMotionFn, 'dusEuclidean', 'dusGeodesic', ...
            'tracks1', 'tracks2', 'pairIDs', 'U0s', 'V0s', ...
            'geodesicPaths', 'ptBarycenters', 'ptFaceLocations', ...
            'euclideanDistanceTraveledU', 'euclideanDistanceTraveledV', ...
            'v3d_u', 'v3d_v', 'nSaved') 
        % track1 3 & track2 191: t = 17 geodesic distance
        if exist(relMotionFn, 'file')
            movefile(relMotionFn, [relMotionFn '_backup'])
        end

    else
        % For each track in tracks2, find initially nearest track in tracks1
        dusGeodesic = zeros(nTracks, nTimePoints) ;
        pairIDs = zeros(nTracks, 1) ;
        U0s = zeros(nTracks, 2) ;
        V0s = zeros(nTracks, 2) ;
        geodesicPaths = cell(nTracks, 1) ;
        ptBarycenters = cell(nTracks, 1) ;
        ptFaceLocations = cell(nTracks, 1) ;
        euclideanDistanceTraveledU = dusGeodesic ;
        euclideanDistanceTraveledV = dusGeodesic ;
        v3d_u = zeros(nTracks, nTimePoints, 3) ;
        v3d_v = zeros(nTracks, nTimePoints, 3) ;
        nSaved = 1 ;
    end

    for ii = nSaved:nTracks
        disp(['matching track ' num2str(ii)])

        % Get starting positions u0 and v0 for layer1 and 2
        V0 = tracks2{ii}(1, 1:2) ;
        V0s(ii, :) = V0 ;
        VV = tracks2{ii}(:, 1:2) ;

        % Look for initially nearby nuclei in pullback space
        dists = vecnorm(nearby - V0, 2, 2) ;
        [~, minID] = min(dists) ;
        UU = tracks1{minID}(:, 1:2) ;
        U0s(ii, :) = UU(1, 1:2) ;

        % Project into 3D
        if strcmpi(coordSys, 'spsm') || strcmpi(coordSys, 'spsmrs')
            im = imread(fullfile(QS.dir.im_sp_sm, subdir1, ...
                sprintf(QS.fileBase.im_sp_sm, timePoints(1)))) ;
            im2 = imread(fullfile(QS.dir.im_sp_sm, subdir2, ...
                sprintf(QS.fileBase.im_sp_sm, timePoints(1)))) ;
            assert(all(size(im) == size(im2)))
        end
        uu = QS.XY2uv(im, UU, doubleCovered, umax, vmax) ;
        vv = QS.XY2uv(im, VV, doubleCovered, umax, vmax) ;

        % For each timepoint, project into 3d for that mesh
        geodesics = cell(nTimePoints, 1) ;
        ptBarys = cell(nTimePoints, 1) ;
        ptFaces = cell(nTimePoints, 1) ;
        duEuclidean = zeros(nTimePoints, 1) ;
        duGeodesic = zeros(nTimePoints, 1) ;
        Deltau = duGeodesic ;
        Deltav = duGeodesic ;
        u3ds = zeros(nTimePoints, 3) ;
        v3ds = zeros(nTimePoints, 3) ;
        for tidx = 1:nTimePoints
            tp = timePoints(tidx) ;
            QS.setTime(tp) 

            if ~any(isnan(uu(tidx, :))) && ~any(isnan(vv(tidx, :))) && ...
                all(uu(tidx, :) > 0) && all(vv(tidx, :) > 0)
                [u3d, fieldfacesU, ~, barycU] = QS.uv2APDV(uu(tidx, :), coordSys) ;
                [v3d, fieldfacesV, ~, barycV] = QS.uv2APDV(vv(tidx, :), coordSys) ;
                
                if tidx == 1 
                    u3d0 = u3d ;
                    v3d0 = v3d ;
                    fieldfacesU0 = fieldfacesU ;
                    fieldfacesV0 = fieldfacesV ;
                end

                progressString = ['trackA ' num2str(ii) ' & trackB ' num2str(minID)...
                        ': t = ' num2str(timePoints(tidx)) ...
                        ' | U=[' sprintf('%0.0f,%0.0f', UU(tidx, 1), UU(tidx, 2)) ...
                        '], V=[' sprintf('%0.0f,%0.0f', VV(tidx, 1), VV(tidx, 2)) ']'] ;

                % Euclidean distance in 3d     
                duEuclidean(tidx) = vecnorm(u3d-v3d, 2, 2) ;
                
                if fieldfacesU == fieldfacesV 
                    % Points are very close --> on the same face: 
                    % Euclidean distance in 3D is equal to the geodesic 
                    % distance      
                    disp(progressString)
                    duGeodesic(tidx) = duEuclidean(tidx) ;
                    ptFaces(tidx, :) = [fieldfacesU, fieldfacesV] ;
                else
                    % Measure distance of these tracks over surface as geodesic
                    disp([progressString ' geodesic distance'])
                    mesh = QS.getCurrentSPCutMeshSmRSC() ;
                    nvtx = size(mesh.v, 1) ;

                    if tidx == 17 && ii == 3 && minID == 191
                        pause
                    end

                    [geodesicP, pointLocations] = surfaceGeodesicPairs( mesh.f, ...
                            mesh.v, [nvtx+1, nvtx+2], [mesh.v; u3d; v3d] ) ;
                    ptFace = [ pointLocations(nvtx+1).face, ...
                        pointLocations(nvtx+2).face ] ;
                    ptBary = [ pointLocations(nvtx+1).barycentricCoordinates; ...
                        pointLocations(nvtx+2).barycentricCoordinates ] ;
                    pp = geodesicP{1} ;
                    duGeodesic(tidx) = sum(vecnorm(diff(pp), 2, 2)) ;

                    geodesics{tidx} = pp ;
                    ptBarys(tidx, :, :) = ptBary ;
                    ptFaces(tidx, :) = ptFace ;
                    
                    assert(all(barycU == ptBary(1, :)))
                    assert(all(barycV == ptBary(2, :)))
                    assert(all(ptFace == [fieldfacesU, fieldfacesV] ))
                end

                % Compare to distance traveled
                % if fieldfacesU == fieldfacesU0 || strcmpi(distMode, 'fastEuclidean')
                Deltau(tidx) = vecnorm(u3d - u3d0, 2, 2) ;
                % if fieldfacesV == fieldfacesV0 || strcmpi(distMode, 'fastEuclidean')
                Deltav(tidx) = vecnorm(v3d - v3d0, 2, 2) ;

                % Keep 3d positions
                u3ds(tidx, :) = u3d ;
                v3ds(tidx, :) = v3d ;
            else
                disp('NaN for one track!')
                duGeodesic(tidx) = NaN ;
                Deltau(tidx) = NaN ;
                Deltav(tidx) = NaN ;
                u3ds(tidx, :) = NaN ;
                v3ds(tidx, :) = NaN ;
            end
        end
        
        % Look for initially nearby nuclei in pullback space
        dists = vecnorm(nearby - V0, 2, 2) ;
        [~, minID] = min(dists) ;
        UU = tracks1{minID}(:, 1:2) ;
        U0s(ii, :) = UU(1, 1:2) ;


        % Store all these vectors
        geodesicPaths{ii} = geodesics ;
        ptBarycenters{ii} = ptBarys ;
        ptFaceLocations{ii} = ptFaces ;
        euclideanDistanceTraveledU(ii, :) = Deltau ;
        euclideanDistanceTraveledV(ii, :) = Deltav ;
        dusEuclidean(ii, :) = duEuclidean ;
        dusGeodesic(ii, :) = duGeodesic ;
        pairIDs(ii) = minID ;
        v3d_u(ii, :, :) = u3ds ;
        v3d_v(ii, :, :) = v3ds ;
        nSaved = ii ;
        disp(['Saving track pairs [' num2str(ii) '/' num2str(nTracks) '] to ' relMotionFn])
        save(relMotionFn, 'dusEuclidean', 'dusGeodesic', ...
            'tracks1', 'tracks2', 'pairIDs', 'U0s', 'V0s', ...
            'geodesicPaths', 'ptBarycenters', 'ptFaceLocations', ...
            'euclideanDistanceTraveledU', 'euclideanDistanceTraveledV', ...
            'v3d_u', 'v3d_v', 'nSaved') 
        
        % xlim, ylim
        Xlim = U0s(ii, 1) + [-200, 200] ;
        Ylim = U0s(ii, 2) + [-200, 200] ;

        % Plot starting correspondences
        if preview 
            clf
            sz = 200 ;
            % Axis 1
            ax0 = subtightplot(2, 2, 1) ;
            imshow(0.25*im);
            hold on;
            scatter(nearby(:, 1), nearby(:, 2), sz, [0, 0.447, 0.741], 'filled')
            scatter(U0s(:, 1), U0s(:, 2), sz, colors, 'filled', 'markeredgecolor', 'k')
            title([subdir1 ' positions'])
            xlim(Xlim)
            ylim(Ylim)
            % Axis 2
            ax1 = subtightplot(2, 2, 2) ;
            imshow(im2);
            hold on;
            scatter(V0s(:, 1), V0s(:, 2), sz*4, colors, 's', 'linewidth', 5)
            title([subdir2 ' positions'])
            xlim(Xlim)
            ylim(Ylim)
            % Axis 3
            ax2 = subtightplot(2, 1, 2) ;
            imshow(min(0.25*im+im2, 255));
            hold on;
            scatter(nearby(:, 1), nearby(:, 2), sz, [0, 0.447, 0.741], 'filled')
            scatter(V0s(:, 1), V0s(:, 2), sz*4, colors, 's', 'linewidth', 5)
            scatter(U0s(:, 1), U0s(:, 2), sz, colors, 'filled', 'markeredgecolor', 'k')
            for jj = 1:nTracks
                plot([U0s(jj, 1), V0s(jj, 1)], [U0s(jj, 2), V0s(jj, 2)], '-', ...
                    'color', colors(jj, :))
            end
            hold off;
            xlim(Xlim)
            ylim(Ylim)
            pause(1)
        end
    end
    
    % Store in-plane positions too, as array
    Uall = nan(nTracks, nTimePoints, 2) ;
    Vall = nan(nTracks, nTimePoints, 2) ;
    phaseRaw = nan(nTracks, nTimePoints) ;
    magXYRaw = nan(nTracks, nTimePoints) ;
    phaseRelative = nan(nTracks, nTimePoints) ;
    magXYRelative = nan(nTracks, nTimePoints) ;
    for ii = 1:nTracks
        % Get starting positions u0 and v0 for layer1 and 2
        Uall(ii, :, :) = tracks1{pairIDs(ii)}(:, 1:2) ;
        Vall(ii, :, :) = tracks2{ii}(:, 1:2) ;
        
        dY = Vall(ii, :, 2) - Uall(ii, :, 2) ;
        dX = Vall(ii, :, 1) - Uall(ii, :, 1) ;
        
        phaseRaw(ii, :) = atan2(dY, dX) ;
        magXYRaw(ii, :) = vecnorm(Vall(ii, :, :) - Uall(ii, :, :), 2, 3) ;
        
        dY0 = Vall(ii, 1, 2) - Uall(ii, 1, 2) ;
        dX0 = Vall(ii, 1, 1) - Uall(ii, 1, 1) ;
        
        phaseRelative(ii, :) = atan2(dY - dY0, dX - dX0) ;
        magXYRelative(ii, :) = sqrt((dY - dY0).^2 + (dX - dX0).^2) ;
        
    end
    
    % Convert to approximate inplane distance
    % Convert to approximate distance using metric
    bcmetric
    for tidx = 1:nTimePoints
        QS.setTime(timePoints(tidx)) ;
        mesh = QS.getCurrentSPCutMeshSmRSC() ;
        % [gcell, bcell] = constructFundamentalForms(mesh.f, mesh.v, mesh.u) ;
        [~, v0t, v0t2d] = ...
            resolveTangentNormalVelocities(faces, vertices, v0s, fieldfaces, ...
                mesh.u) ;

        magUVRelative = metric
    end

    
    % First consider raw relative distances, then displacements from
    % starting positions
    figDir = QS.dir.tracking; 
    for analysisMode = 1:2
        % Raw or relative
        if analysisMode == 1
            % raw displacement vectors
            duG = dusGeodesic ;
            duE = dusEuclidean ;
            phaseV = phaseRaw ;
            exten = '_raw' ;
        else
            % relative to starting points
            duG = dusGeodesic - dusGeodesic(:, 1) ;
            duE = dusEuclidean - dusEuclidean(:, 1) ;
            phaseV = phaseRelative ;
            exten = '' ;
        end
        
        
        %% 1D plot: curves for distance over time
        close all
        hf = figure('Position', [100 100 800 340], 'units', 'centimeters');
        h1 = subplot(1, 3, 1) ;
        plot(timePoints, duG)
        ylabel(['relative motion of nuclei between layers [' QS.spaceUnits ']'], ...
            'interpreter', 'latex')
        xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
        if analysisMode == 1
            title('Geodesic separation')
        else
            title('Change in geodesic separation')
        end
        ylimsRel = ylim ;

        h1 = subplot(1, 3, 2) ;
        plot(timePoints, duE)
        ylabel(['relative motion of nuclei between layers [' QS.spaceUnits ']'], ...
            'interpreter', 'latex')
        xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
        if analysisMode == 1
            title('Euclidean separation in 3D')
        else
            title('Change in Euclidean separation in 3D')
        end
        ylimsRel = ylim ;
        
        h2 = subplot(1, 3, 3) ;
        plot(timePoints, euclideanDistanceTraveledU')
        ylabel(['distance traveled by endoderm nuclei [' QS.spaceUnits ']'], ...
            'interpreter', 'latex')
        xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
        title('Euclidean distance traveled by endoderm in 3D')
        ylimsInp = ylim ;

        % set axis limits to match
        ylims = [min(ylimsRel(1), ylimsInp(1)), max(ylimsRel(2), ylimsInp(2))] ;
        ylim(ylims)
        axes(h1)
        ylim(ylims)
        saveas(gcf, fullfile(figDir, ['relative_motion_' distMode '_individual_curves' exten '.pdf']))
        
        
        %% Polar/2d scatter plot over time
        figure ;
        if analysisMode == 1
            duV = duG ;
        else
            duV = duE ;
        end
        
        alphaVal = 0.1 ;
        for ii = 1:nTracks
            lastIdx = find(isnan(duG(ii, :))) ;
            if isempty(lastIdx)
                xx = abs(duV(ii, :)) .* cos(phaseV(ii, :)) ;
                yy = abs(duV(ii, :)) .* sin(phaseV(ii, :)) ;
                aColor = 1:nTimePoints ;
            else
                xx = abs(duV(ii, 1:lastIdx-1)) .* cos(phaseV(ii, 1:lastIdx-1)) ;
                yy = abs(duV(ii, 1:lastIdx-1)) .* sin(phaseV(ii, 1:lastIdx-1)) ;
                aColor = 1:lastIdx-1 ;
            end
            mfig = patch([xx NaN], [yy NaN], [aColor nTimePoints]);
            % alphamap = 0.2 * ones(length(xx), 1) ;
            % alphamap(end) = 0 ;
            set(mfig,'FaceColor','none','EdgeColor','flat',...
               'LineWidth',2, ...
               'FaceVertexAlphaData',alphaVal,...
               'EdgeAlpha',alphaVal)
           pause
        end
        % Take means for each timepoint
        meanX = mean(du
        meanY = 
        xlabel('relative ap displacement')
        ylabel('relative dv displacement')
        saveas(gcf, fullfile(figDir, ['relative_motion_' distMode '_inPlane' exten '.pdf']))

        %% Histogram version
        fewTidx = 1:15:nTimePoints ;
        tColors = viridis(length(fewTidx)) ;
        qq = 1 ;
        close all
        for tidx = fewTidx
            polarhistogram(phaseV(:, tidx), 15, ...
                'EdgeColor', 'none', 'FaceColor', tColors(qq, :), ...
                'FaceAlpha',.2)
            hold on;
            pause(0.2)
            qq = qq + 1 ;
        end    
    end
        
    %% Plot all starting correspondences
    clf
    sz = 20 ;
    % Axis 1
    ax0 = subtightplot(2, 2, 1) ;
    imshow(0.25*im);
    hold on;
    scatter(nearby(:, 1), nearby(:, 2), sz*0.2, [0, 0.447, 0.741])
    scatter(U0s(:, 1), U0s(:, 2), sz, colors, 'filled')
    title([subdir1 ' positions'])
    % Axis 2
    ax1 = subtightplot(2, 2, 2) ;
    imshow(im2);
    hold on;
    scatter(V0s(:, 1), V0s(:, 2), sz*2, colors, 's')
    title([subdir2 ' positions'])
    % Axis 3
    ax2 = subtightplot(2, 1, 2) ;
    imshow(min(0.25*im+im2, 255));
    hold on;
    scatter(U0s(:, 1), U0s(:, 2), sz, colors, 'filled')
    scatter(V0s(:, 1), V0s(:, 2), sz*2, colors, 's')
    for jj = 1:nTracks
        plot([U0s(jj, 1), V0s(jj, 1)], [U0s(jj, 2), V0s(jj, 2)], '-', ...
            'color', colors(jj, :))
    end
    hold off;
    saveas(gcf, relMotionStartFn)
    
    
    %% Filter bad points
    if ~isnan(max_du) && isfinite(max_du)
        dusGeodesic(dusGeodesic(:) > max_du) = NaN ;
    end
    %% Parameterize by time -- mean+/- std
    figDir = QS.dir.tracking ;
    du_tp = dusGeodesic - dusGeodesic(:, 1) .* ones(size(dusGeodesic)) ;
    means_du = mean(du_tp, 1, 'omitnan') ;
    stds_du = std(du_tp, 1, 'omitnan') ;

    means_dd = mean(euclideanDistanceTraveledU, 1, 'omitnan') ;
    stds_dd = std(euclideanDistanceTraveledU, 1, 'omitnan') ;

    % Handle minutes to hours conversion
    if contains(lower(QS.timeUnits), 'min')
        timeConversion = 1/60 ;
        xlabelString = 'time [hr]' ;
    else
        xlabelString = ['time [' QS.timeUnits ']'] ;
    end

    % Plot it
    close all
    bigRed = [ 0.62    0.76    0.84 ];
    bigBlue = [ 0.90    0.55    0.55] ;
    big3 = [ 0.89    0.10    0.11] ;
    big4 = [ 0.12    0.47    0.70] ;
    lineProps1 = {'-','color', bigRed} ;
    lineProps2 = {'-','color', bigBlue} ;

    hf = figure('Position', [100 100 300 200], 'units', 'centimeters');
    h1 = shadedErrorBar(timePoints * QS.timeInterval * timeConversion, ...
        means_du, stds_du, 'lineProps', lineProps1) ;
    hold on;
    h2 = shadedErrorBar(timePoints * QS.timeInterval * timeConversion, ...
        means_dd, stds_dd, 'lineProps', lineProps2) ;
    legend('motion', 'relative motion')
    ylabelString = ['displacement [' QS.spaceUnits ']'] ;
    xlabel(xlabelString, 'interpreter', 'latex')
    ylabel(ylabelString, 'interpreter', 'latex')
    ylim([0, max(means_dd + stds_dd)])
    saveas(gcf, fullfile(figDir, ['relative_motion_' distMode '.pdf']))


    %% Plot in 3d using initial du
    du0 = dusGeodesic(1, :) ;
    [du0, indx] = sort(du0) ;
    assert(length(du0) == nTimePoints)
    for qq = 1:nTimePoints
        curv = dusGeodesic(indx(qq), :) ;
        deltau = euclideanDistanceTraveledU(indx(qq), :) ;
        plot3(du0(qq)* ones(size(curv)), timePoints, curv, '.-', 'color', bluecolor)
        plot3(du0(qq)* ones(size(curv)), timePoints, deltau, '.-', 'color', orange)
        hold on;
    end
    pause
    
    %% Plot phase of displacement at some timepoints
    
    
end

%% Plot the track pairs in 3d colored by pairID
% for hand dataset, load metadat options for texturepatch:
% 
Options = struct() ; 
Options.layerLabel = 'muscle' ;
QS.plotRelativeMotionTracks(Options)


%% Ensure all pairIDs are good tracks in endoderm layer

subdir = 'endoderm' ;
imDir = fullfile(QS.dir.im_sp_sm, subdir) ;
timePoints = 1:60 ;
fileBase = fullfile(imDir, QS.fileBase.im_sp_sm) ;
trackOutfn = fullfile(QS.dir.tracking, 'endoderm_tracks_correction.mat') ;
load(trackOutfn, 'tracks') ;
tracks2Correct = pairIDs(pairIDs > 100 & pairIDs < 541) ;
tracks2Correct =  tracks2Correct(18:end)
tracks2Correct = unique(tracks2Correct) ;
[newTracks, newG] = manualCorrectTracks2D(tracks, fileBase, timePoints, trackOutfn, tracks2Correct) ;

% 104: tp 47-end:NaN
% finished 1-12 of tracks2Correct











%% Assign Basic Lineage IDs ===============================================

%--------------------------------------------------------------------------
% Determine the 'progenitor' cells which have no parents in the tracking
% graph structure
%--------------------------------------------------------------------------

% Find the in-degree of each node
inDeg = indegree(G);

assert(isequal(unique(inDeg), [0; 1]), 'Invalid tracking structure');

progCells = find(inDeg == 0);

clear inDeg

%--------------------------------------------------------------------------
% Add a lineage ID field to each node of the tracking graph
%--------------------------------------------------------------------------

lineageID = -ones(size(G.Nodes,1),1);
lineageID(progCells) = 1:numel(progCells);

G.Nodes.LineageID = lineageID;

clear lineageID

%--------------------------------------------------------------------------
% Update the lineage ID field of each node using a depth-first search
%--------------------------------------------------------------------------

for i = 1:numel(progCells)
    
    % Get the IDs of all nodes descended from the current progenitor cell
    descNodes = dfsearch(G, progCells(i), ...
        { 'discovernode', 'edgetonew' } );
    descNodes = descNodes.Node;
    descNodes = descNodes( ~isnan(descNodes) );
    
    descNodes( ismember(descNodes, progCells) ) = [];
    
    % Udpate the lineage ID
    G.Nodes(descNodes,:).LineageID = repmat(i, numel(descNodes), 1);
    
end

clear descNodes

assert( ~any(G.Nodes.LineageID < 0), 'Lineage improperly assigned!' );




