
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
cd /mnt/crunch/gut/handGAL4klarHandGFPhistGFP/202105072030_1p4um_0p2ms_1mWGFP/view4_1p4um_25x_1mW_exp0p2ms_120spf/data ;

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
    % ofn_smoothply = 'mesh_apical_stab_' ; % mesh_apical_stab_' ;
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
        'ofn_smoothply', ofn_smoothply, ...
        'pre_nu', pre_nu, ...
        'pre_smoothing', pre_smoothing, ...
        'mlxprogram', mlxprogram, ...
        'init_ls_fn', init_ls_fn, ... % set to none to load prev tp
        'run_full_dataset', run_full_dataset,... % projectDir, ... % set to 'none' for single tp
        'radius_guess', radius_guess, ...
        'dset_name', 'exported_data',...
        'center_guess', center_guess,... % xyz of the initial guess sphere ;
        'save', true, ... % whether to save images of debugging output
        'plot_mesh3d', false, ...
        'dtype', dtype,...
        'mask', mask,...
        'mesh_from_pointcloud', false, ...
        'prob_searchstr', prob_searchstr, ...
        'preilastikaxisorder', preilastikaxisorder, ... 
        'ilastikaxisorder', ilastikaxisorder, ... 
        'physicalaxisorder', imsaneaxisorder, ... 
        'include_boundary_faces', include_boundary_faces, ...
        'smooth_with_matlab', smooth_with_matlab) ;

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

    
    
%% Open manual tracking gui
addpath_recurse('/mnt/data/code/ParhyaleCellTracker/')
close all
Options = struct() ;
Options.allowOverlaps = true  ;
Options.drawGenerations = false ;
Options.colorByLineage = true ;
Options.infoIndex = 1 ;
[Gout, divStruct] = parhyale_master_gui(GG, rawImFileBase, Options) ;



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




