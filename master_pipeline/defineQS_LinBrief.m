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
%cd /mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data
 cd /home/yuzhenglin/membraneproject/caaxdataset
% .=========.
% |  VIP10  |
% .=========.
% cd /mnt/crunch/gut/48YGal4UasLifeActRuby/201907311600_48YGal4UasLifeActRuby_60s_exp0p150_1p0mW_25x_1p4um
% cd /mnt/crunch/gut/48YGal4klarUASCAAXmChHiFP/202001221000_60sec_1p4um_25x_1mW_2mW_exp0p25_exp0p7/Time3views_1017/data/
% cd /mnt/crunch/gut/Mef2Gal4klarUASCAAXmChHiFP/202003151700_1p4um_0p5ms3msexp/Time3views_1/data/
dataDir = cd ;
meshDir = fullfile(dataDir, 'msls_output')

% PATHS ==================================================================
origpath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(origpath))
aux_paths_and_colors
cd(dataDir)

% DEFINE NEW MASTER SETTINGS
overwrite_masterSettings = false ;
if overwrite_masterSettings || ~exist('./masterSettings.mat', 'file')
    % Metadata about the experiment
    stackResolution = [.2619 .2619 .2619] ;
    nChannels = 1 ;
    channelsUsed = 1 ;
    timePoints = 10:263 ;
    ssfactor = 4 ;
    % whether the data is stored inverted relative to real position
    flipy = true ; 
    timeInterval = 1 ;  % physical interval between timepoints
    timeUnits = 'min' ; % physical unit of time between timepoints
    spaceUnits = '$\mu$m' ; % physical unit of time between timepoints
    scale = 0.04 ;      % scale for conversion to 16 bit
    % file32Base = 'TP%d_Ch0_Ill0_Ang0,45,90,135,180,225,270,315.tif'; 
    file32Base = 'TP%d_Ch0_Ill0_Ang0,60,120,180,240,300.tif'; 
    % file32Base = 'TP%d_Ch0_Ill0_Ang0,60,120,180,240,300.tif'; 
    fn = 'Time_%06d_c1_stab';
    fn_prestab = 'Time_%06d_c1.tif';
    set_preilastikaxisorder = 'xyzc' ;
    swapZT = 1 ;
    masterSettings = struct('stackResolution', stackResolution, ...
        'nChannels', nChannels, ...
        'channelsUsed', channelsUsed, ...
        'timePoints', timePoints, ...
        'ssfactor', ssfactor, ...
        'flipy', flipy, ...
        'timeInterval', timeInterval, ...
        'timeUnits', timeUnits, ...
        'spaceUnits', timeUnits, ...
        'scale', scale, ...
        'file32Base', file32Base, ...
        'fn', fn,...
        'fn_prestab', fn_prestab, ...
        'swapZT', swapZT, ...
        'set_preilastikaxisorder', set_preilastikaxisorder, ...
        't0_for_phi0', 110, ... % 40 for mef2 single channel, 110 for CAAX excellent
        'nU', 100, ...  % 150 for mef2 data with posterior midgut loop
        'nV', 100); 
    disp('Saving masterSettings to ./masterSettings.mat')
    if exist('./masterSettings.mat', 'file')
        ui = input('This will overwrite the masterSettings. Proceed (Y/n)?', 's') ;
        if ~isempty(ui) && (strcmp(ui(1), 'Y') || strcmp(ui(1), 'y'))
            save('./masterSettings.mat', 'masterSettings')
            loadMaster = false ;
        else
            disp('Loading masterSettings from disk instead of overwriting')
            loadMaster = true ;
        end
    else
        save('./masterSettings.mat', 'masterSettings')
        loadMaster = false ;
    end
else
    loadMaster = true ;
end


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

clear fileMeta expMeta
 
% % SET DETECT OPTIONS =====================================================
% % Must run this section for later functionality.
% % Load/define the surface detection parameters
% msls_detOpts_fn = fullfile(projectDir, 'msls_detectOpts.mat') ;
% load(msls_detOpts_fn, 'detectOptions')
% 
% % Overwrite certain parameters for script structure
% detectOptions.fileName = sprintf( fn, xp.currentTime ) ;
% detectOptions.run_full_dataset = false ;
% detectOptions.ms_scriptDir = ms_scriptDir ;
% meshDir = detectOptions.mslsDir ;
% % These are now in QS.
% % meshFileBase = [ofn_smoothply '%06d'] ;
% % alignedMeshBase = [ofn_smoothply '%06d_APDV_um'] ;
% 
% % Set detect options ------------------------------------------------------
% xp.setDetectOptions( detectOptions );
% disp('done')
xp.detectOptions = struct() ;
xp.detectOptions.ssfactor = 4 ;
xp.detector.options.ofn_smoothply = 'mesh_%06d' ;

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

%% Cell Segmentation
% options = struct() ;
% options.overwrite = true;
% options.cellSize = 50 ;
% options.strelRadius = 0 ;
% options.gaussKernel = 1 ;
% options.heightMinimum = 2.2 ;
% options.coordSys= 'sprsme';
% options.timePoints = 213:15:258;% 93:15:200 ;
% QS.generateCellSegmentation2D(options) 

% Now manual corrections in GIMP
options = struct() ;
options.overwrite = true ;
options.overwriteImages = true;
options.timePoints = 96:1:205;
QS.processCorrectedCellSegmentation2D(options) 

% Generate
options.overwrite = true ;
options.correctedSegmentation = true;
QS.generateCellSegmentation3D(options) 

