%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quick instantiation of a QuapSlap object 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace ========================================================
% We start by clearing the memory and closing all figures
clear; close all; clc;
% change this path, for convenience
% cd /mnt/crunch/48Ygal4-UAShistRFP/201904031830_great/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2/data/
% cd /mnt/crunch/48YGal4UasLifeActRuby/201904021800_great/Time6views_60sec_1p4um_25x_1p0mW_exp0p150_3/data/
% cd /mnt/data/48YGal4UasLifeActRuby/201902201200_unusualfolds/Time6views_60sec_1p4um_25x_obis1_exp0p35_3/data/
cd /mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data
% .=========.
% |  VIP10  |
% .=========.
% cd /mnt/crunch/gut/48YGal4UasLifeActRuby/201907311600_48YGal4UasLifeActRuby_60s_exp0p150_1p0mW_25x_1p4um
% cd /mnt/crunch/gut/48YGal4klarUASCAAXmChHiFP/202001221000_60sec_1p4um_25x_1mW_2mW_exp0p25_exp0p7/Time3views_1017/data/
% cd /mnt/crunch/gut/Mef2Gal4klarUASCAAXmChHiFP/202003151700_1p4um_0p5ms3msexp/Time3views_1/data/
dataDir = cd ;

% PATHS ==================================================================
origpath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(origpath))
aux_paths_and_colors
cd(dataDir)

% Global options
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
overwrite_cutMesh = false ;
overwrite_spcutMesh = false ;
overwrite_pullbacks = false ;

% DEFINE NEW MASTER SETTINGS
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

% Try loading
timeInterval = masterSettings.timeInterval ;  % physical interval between timepoints
timeUnits = masterSettings.timeUnits ; % physical unit of time between timepoints        
spaceUnits = masterSettings.spaceUnits ;  % microns as $\mu$m


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
fileMeta.swapZT             = 1;

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

% SET DETECT OPTIONS =====================================================
% Must run this section for later functionality.
% Load/define the surface detection parameters
msls_detOpts_fn = fullfile(projectDir, 'msls_detectOpts.mat') ;
load(msls_detOpts_fn, 'detectOptions')

% Overwrite certain parameters for script structure
detectOptions.fileName = sprintf( fn, xp.currentTime ) ;
detectOptions.run_full_dataset = false ;
detectOptions.ms_scriptDir = ms_scriptDir ;
meshDir = detectOptions.mslsDir ;
% These are now in QS.
% meshFileBase = [ofn_smoothply '%06d'] ;
% alignedMeshBase = [ofn_smoothply '%06d_APDV_um'] ;

% Set detect options ------------------------------------------------------
xp.setDetectOptions( detectOptions );
disp('done')

%% QS DEFINITION
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
QS = QuapSlap(xp, opts) ;
disp('done')
