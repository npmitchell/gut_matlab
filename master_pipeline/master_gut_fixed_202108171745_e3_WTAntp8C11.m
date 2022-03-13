%% MASTER GUT PULLBACK -- FIXED SAMPLE
% NPMitchell
%
% This is a pipeline to take the surface of the fixed Drosophila gut and
% map to the plane

% We start by clearing the memory and closing all figures
clear; close all; clc;

% temporary path def
cd /mnt/data/antibodies_lightsheet/48YGAL4CAAXmCh_antp8C11_1t50/202108171745_e3_16a_1p2um_0p3ms0p5ms_3mWGFP_20mWRFP/
dataDir = pwd;

%% IMSANE SETUP FOR DETECTOR
% 
% cd('/mnt/data/code/imsane_for_git/imsane/')
% Run setup
% cd(projectDir)

%% PATHS ==================================================================
origpath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(origpath))
aux_paths_and_colors
cd(dataDir)

%% DEFINE NEW MASTER SETTINGS
overwrite_masterSettings = true ;
if overwrite_masterSettings || ~exist('./masterSettings.mat', 'file')
    % Metadata about the experiment
    stackResolution = [.2619 .2619 .2619] ;
    nChannels = 2 ;
    channelsUsed = [1,2] ;
    timePoints = 1 ;
    ssfactor = 4 ;
    % whether the data is stored inverted relative to real position
    flipy = true ; 
    timeInterval = 1 ;  % physical interval between timepoints
    timeUnits = 'min' ; % physical unit of time between timepoints
    scale = [1.0 0.25] ;      % scale for conversion to 16 bit
    file32Base = 'TP%d_Ch%d_Ill0_Ang0,45,90,135,180,225,270,315.tif'; 
    % file32Base = 'TP%d_Ch0_Ill0_Ang0,60,120,180,240,300.tif'; 
    fn = 'fixed_sample_c1';
    spaceUnits = '$\mu$m';  % microns as $\mu$m
    swapZT = 1 ;
    nU = 100 ;
    nV = 100 ;
    
    set_preilastikaxisorder = 'xyzc' ;
    masterSettings = struct('stackResolution', stackResolution, ...
        'nChannels', nChannels, ...
        'channelsUsed', channelsUsed, ...
        'timePoints', timePoints, ...
        'ssfactor', ssfactor, ...
        'flipy', flipy, ...
        'timeInterval', timeInterval, ...
        'timeUnits', timeUnits, ...
        'spaceUnits', spaceUnits, ...
        'scale', scale, ...
        'file32Base', file32Base, ...
        'fn', fn,...
        'swapZT', swapZT, ...
        'nU', nU, ...
        'nV', nV, ...
        'set_preilastikaxisorder', set_preilastikaxisorder); 
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

if loadMaster
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
    
    % Try loading
    timeInterval = masterSettings.timeInterval ;  % physical interval between timepoints
    timeUnits = masterSettings.timeUnits ; % physical unit of time between timepoints        
    spaceUnits = masterSettings.spaceUnits ;  % microns as $\mu$m
    
    
    % Fill in
    swapZT = masterSettings.swapZT ;
    nU = masterSettings.nU ;
    nV = masterSettings.nV ;


    scale = masterSettings.scale ;      % scale for conversion to 16 bit
    file32Base = masterSettings.file32Base ; 
    fn = masterSettings.fn ;
    fn_prestab = masterSettings.fn_prestab ;
    set_preilastikaxisorder = masterSettings.set_preilastikaxisorder ;
end
dir16bit = dataDir ; % fullfile(dataDir, 'deconvolved_16bit') ;

%% INITIALIZE ImSAnE PROJECT ==============================================
% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.

dataDir    =  cd; 
projectDir = dataDir ;
% [ projectDir, ~, ~ ] = fileparts(matlab.desktop.editor.getActiveFilename); 

%% CREATE EXPERIMENT
% Start by creating an experiment object, optionally pass on the project
% directory (otherwise it will ask), and change into the directory of the
% data.  This serves as a front-end for data loading, detection, fitting
% etc.
xp = project.Experiment(projectDir, dataDir);

% Set file and experiment meta data
% Set required additional information on the files.
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
%
% * 'imageSpace'        , bit depth of image, such as uint16 etc., defined
%                         in Stack class
% * 'stackSize'         , size of stack in pixels per dimension 
%                         [xSize ySize zSize]
% * 'swapZT'            , set=1 if time is 3rd dimension and z is 4th

% A filename base template - to be used throughout this script
% the 32 bit fn
fn = 'TP0_Ch%d_Ill0_Ang0,45,90,135,180,225,270,315' ;
% the 16 bit fn
file16name = 'fixed_sample_c%d' ;                   

fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fn, '.tif'];
fileMeta.nChannels          = 2;
fileMeta.timePoints         = timePoints ;
fileMeta.stackResolution    = [.2619 .2619 .2619];
fileMeta.swapZT             = swapZT ;

% Set required additional information on the experiment. A verbal data set
% description, Jitter correct by translating  the sample, which time point
% to use for fitting, etc.
%
% The following project metadata information is required:
%
% * 'channelsUsed'      , the channels used, e.g. [1 3] for RGB
% * 'channelColor'      , mapping from element in channels used to RGB = 123
% * 'dynamicSurface'    , Not implemented yet, future plan: boolean, false: static surface
% * 'detectorType'      , name of detector class, e.g. radielEdgeDetector
%                         ,(user threshholded), fastCylinderDetector
% * 'fitterType'        , name of fitter class
%
% The following project meta data information is optional:
%
% * 'description'     , string describing the data set set experiments metadata, 
%                                such as a description, and if the surface is dynamic,
%                                or requires drift correction of the sample.
% * 'jitterCorrection', Boolean, false: No fft based jitter correction 

% first_tp is also required, which sets the tp to do individually.
first_tp = 1 ;
expMeta                     = struct();
expMeta.channelsUsed        = channelsUsed;
expMeta.channelColor        = [1, 2];
expMeta.description         = 'Drosophila gut';
expMeta.dynamicSurface      = 0;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.fitTime             = fileMeta.timePoints(first_tp);
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';

%% 32 to 16 BIT CONVERSION
% Check that data is 16 bit. If not, convert to 16bit
pc2use = 99.99;
scalemethod = 'user-defined' ; % use either 'user-defined' or 'prctile'

assert(length(scale) == length(expMeta.channelsUsed)) 

for channel_check = expMeta.channelsUsed
    scalemax = scale(channel_check) ;
    fullFileName = [sprintf(fn, channel_check) '.tif'] ;
    info = imfinfo(fullFileName) ;
    full16fn = [sprintf(file16name, channel_check) '.tif'] ;
    bitDepth = info.BitDepth ;

    if (bitDepth == 32) && ~isfile(full16fn)

        disp([fullFileName ' is not 16bit, converting to ' full16fn])

        % Note that imread only loads a single frame
        % A = imread(fullFileName) ;
        % scalemin = double(min(A(:))) ;
        % scalemax = double(max(A(:))) ;
        disp('')
        disp('Reading 32 bit file to convert...')
        A = readSingleTiff(fullFileName) ;
        tmpA = A(:) ;
        disp('')
        disp('Computing scalemin, scalemax')

        % Optional step here to figure out what the cutoff
        % intensity should be
        % tmpA_no_ouliers = tmpA(tmpA < pcntile(tmpA, 99)) ;
        % thisstd = std(tmpA_no_ouliers) ;
        % check it using histogram(tmpA)
        thismedian = median(tmpA) ;

        %goodmedian = 2559.00;
        %worstmedian = 420.00;
        %range2correct = goodmedian - worstmedian ;
        %normal_pc2use = 99.9999 ;
        %worstcase_pc2use = 99.99 ;
        %diffpc = normal_pc2use - worstcase_pc2use ;
        %pc2use = normal_pc2use + diffpc * (thismedian - goodmedian) / range2correct ;
        %pc2use = max(worstcase_pc2use, pc2use) ;
        %pc2use = min(normal_pc2use, pc2use) ;
        chanpositionstart = strfind(fullFileName,'Ch');
        chanposition = fullFileName(chanpositionstart+2);
        chanposition = str2num(chanposition);
        disp('determining prctile')
        if strcmp('scalemethod', 'prctile')
            scalemax = double(prctile(tmpA, pc2use)) ;
        else
            disp('using user supplied scalemax')
        end
            
        scalemin = double(min(tmpA(tmpA > 0))) ;
        disp(['scalemax = ', num2str(scalemax)])
        disp(['scalemin = ', num2str(scalemin)])

        disp('Showing slice of scaled image')
        % Note to self to check scale:
        close all
        imagesc(squeeze(A(:, 300, :)))
        title('Checking scale of image')
        colorbar
        waitfor(gcf)
        histogram(A(:))
        ylim([0,10])
        title('Histogram of values')
        waitfor(gcf)

        % data = readSingleTiff(fullFileName);
        im2 = mat2gray(A,[scalemin scalemax]);
        im2 = uint16(2^16*im2);
        imSize = size(im2);
        
        % Check scale:
        imagesc(squeeze(im2(:, 300, :)))
        title('Checking scale of image. Close image to continue')
        axis equal
        colorbar
        waitfor(gcf)

        % Save the 16 bit image
        disp(['Saving 16bit volume to ' full16fn])
        imwrite(im2(:,:,1),full16fn,'tiff','Compression','none');
        for z = 2 : imSize(3)
            imwrite(im2(:,:,z),full16fn,'tiff','Compression','none','WriteMode','append');
        end
        disp('done saving 16bit volume')


    elseif isfile(full16fn)
        % the data is already 16 bit, so we're good
        fullFileName = [sprintf(fn, channel_check) '.tif'] ;
        disp([fullFileName ' has been converted to 16bit: ' full16fn ])

    else
        disp('File is 16bit.')
    end
end
fn = file16name ;

%% Collate multiple colors 
overwrite = false ;
fnCombined = 'fixed_sample_combined' ;
fileNameIn = fullfile(dir16bit, [fn '.tif']) ;
fileNameOut = fullfile(dir16bit, [fnCombined '.tif']) ;
if ~exist(fileNameOut, 'file') || overwrite
    collateColors(fileNameIn, fileNameOut, [], channelsUsed) ; 
end
file16name = fnCombined ;

%% INSTANTIATE EXPERIMENT CLASS
% Now set the meta data in the experiment.
fileMeta.filenameFormat = [ file16name '.tif' ] ;
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();

%% LOAD THE FIRST TIME POINT ==============================================
xp.loadTime(xp.fileMeta.timePoints(first_tp));
xp.rescaleStackToUnitAspect();

%% DETECT THE SURFACE =====================================================
% Surface detection parameters --------------------------------------------
% Must run this section for later functionality.
% Mesh extraction options
run_full_dataset = 'none' ;
overwrite_detOpts = true ;

mlxprogram = 'surface_rm_resample20k_reconstruct_LS3_1p2pc_ssfactor4.mlx';
% Mesh marching options
normal_step = -10;

% Load/define the surface detection parameters
msls_detOpts_fn = fullfile(projectDir, 'msls_detectOpts.mat') ;
if exist(msls_detOpts_fn, 'file') && ~overwrite_detOpts
    load(msls_detOpts_fn, 'detectOptions')
else
    % Define the surface detection parameters
    channel = 2;
    foreGroundChannel = 2;
    niter = 50 ;
    niter0 = 50 ;
    ofn_smoothply = 'mesh_' ;
    ofn_ply = 'mesh_ms_' ; 
    ofn_ls = 'msls_' ;
    ms_scriptDir = '/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper' ;
    lambda1 = 1 ;
    lambda2 = 1 ;
    exit_thres = 0.00001 ;
    smoothing = 4 ;
    nu = 0.5 ;
    pre_nu = 2 ;
    pre_smoothing = 1 ;
    post_nu = 2;
    post_smoothing = 8 ;
    radius_guess = 10 ;
    center_guess = 'empty_string' ;
    zdim = 3 ;
    init_ls_fn =  'mesh_initguess' ;
    dtype = 'h5';
    mask = 'none' ;
    prob_searchstr = '_Probabilities.h5';
    preilastikaxisorder = 'xyzc';
    ilastikaxisorder = 'xyzc';
    imsaneaxisorder = 'xyzc'; 
    include_boundary_faces = true ;
    smooth_with_matlab = -1;

    % Name the output mesh directory ------------------------------------------
    if projectDir(end) ~= filesep
        projectDir = [projectDir filesep];
    end
    meshDir = fullfile(projectDir, 'msls_output');

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
        'mlxprogram', fullfile('/mnt/data/code/meshlab_codes', mlxprogram), ...
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
        'pythonVersion', '2'); % version of python to call = '2' or '3', as string
    
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

% Set detect options ------------------------------------------------------
xp.setDetectOptions( detectOptions );

% clear msls_exten imwriteOptions saveDir
% clear channel foreGroundChannel
% clear niter niter0 lambda1 lambda2
% clear exit_thres smoothing nu
% clear post_nu post_smoothing

%% CREATE THE SUBSAMPLED H5 FILE FOR INPUT TO ILASTIK =====================
% skip if already done
for t = xp.fileMeta.timePoints
    h5fn = fullfile(projectDir, [sprintf(sprintf(fn, t)) '.h5']) ;
    if ~exist(h5fn, 'file')
        disp(['Did not find file: ', fullfile(projectDir, [sprintf(fn, t) '.h5'])])

        % Only load and rescale if multiple timepoints/channels
        % if length(xp.expMeta.channelsUsed) > 1 || length(xp.fileMeta.timePoints) > 1
        %     xp.loadTime(t);
        %     xp.rescaleStackToUnitAspect();
        % end
        % make a copy of the detectOptions and change the fileName
        detectOpts2 = detectOptions ;
        detectOpts2.fileName = sprintf( fn, xp.currentTime ) ;
        xp.setDetectOptions( detectOpts2 );
        xp.detector.prepareIlastik(xp.stack);
        disp(['done outputting downsampled data h5: tp=' num2str(t) ' for surface detection'])
    else
        disp(['h5 ' num2str(t) ' was already output, skipping...'])
    end
end    
disp('Open with ilastik if not already done')



%% TRAIN DATA IN ILASTIK TO IDENTIFY APICAL/YOLK ==========================
% open ilastik, train until probabilities and uncertainty are satisfactory

%% Create MorphoSnakesLevelSet from the Probabilities from ilastik ========
xp.detectSurface();
fileMeta = xp.fileMeta ;


%%% Example python run

% python /mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/run_morphsnakes.py \
%  -i fixed_sample_c1_Probabilities.h5 -o \
%  /mnt/data/antibodies_lightsheet/202002201800_w48YGal4klar_DAPI_abdA100_exp2_2mW_25x_1p4um_ms568/Time4views_60sec_1p4um_25x_2mW_exp2/data/msls_output \
%  -prenu -8 -presmooth 1 -ofn_ply mesh_ms_000000.ply -ofn_ls msls_000000.h5 -l1 1 -l2 1 -nu 0.5 -postnu 1 -channel 1 -smooth 4 \
%  -postsmooth 4 -exit 0.000001000 -channel 1 -dtype h5 -permute xyzc -ss 4 -include_boundary_faces -rad0 10 \
%  -init_ls /mnt/data/antibodies_lightsheet/202002201800_w48YGal4klar_DAPI_abdA100_exp2_2mW_25x_1p4um_ms568/Time4views_60sec_1p4um_25x_2mW_exp2/data/msls_output/mesh_initguess.h5 \
%  -n 20 -save
 
command = ['python /mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/run_morphsnakes.py  ' ...
     '-i fixed_sample_c1_Probabilities.h5 -o  ' ...
     dataDir '   ' ...
     '-prenu -6 -presmooth 1 -ofn_ply msls_ms_000000.ply -ofn_ls msls_000000.h5   ' ...
     '-l1 1 -l2 1 -nu 0.5 -postnu 2 -channel 1 -smooth 4  -postsmooth 20  ' ...
     '-exit 0.000001000 -channel 1 -dtype h5 -permute xyzc  ' ...
     '-ss 4 -include_boundary_faces -rad0 10   ' ...
     '-init_ls ' dataDir '/msls_output/msls_initguess.h5  '...
     '-n 15 -save'] 

     %  python /mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/run_morphsnakes.py \
%     -i fixed_sample_c1_Probabilities.h5 -o  \
%     /mnt/data/antibodies_lightsheet/202002201800_w48YGal4klar_DAPI_abdA100_exp2_2mW_25x_1p4um_ms568/Time4views_60sec_1p4um_25x_2mW_exp2/data/msls_output \
%     -prenu -6 -presmooth 1 -ofn_ply msls_ms_000000.ply -ofn_ls msls_000000.h5  \
%     -l1 1 -l2 1 -nu 0.5 -postnu 2 -channel 1 -smooth 4  -postsmooth 20 \
%     -exit 0.000001000 -channel 1 -dtype h5 -permute xyzc \
%     -ss 4 -include_boundary_faces -rad0 10  \
%     -init_ls /mnt/data/antibodies_lightsheet/202002201800_w48YGal4klar_DAPI_abdA100_exp2_2mW_25x_1p4um_ms568/Time4views_60sec_1p4um_25x_2mW_exp2/data/msls_output/mesh_initguess.h5 \
%     -n 15


%% Adjust the mesh to push outward a bit
% ONLY DO THIS ONCE
preSmoothMeshfn = [QS.fullFileBase.mesh(1:end-4) '_presmooth.ply'] ;
if ~exist(preSmoothMeshfn, 'file')
    mesh = read_ply_mod(sprintf(QS.fullFileBase.mesh, timePoints(1))) ;
    amesh = read_ply_mod(sprintf(QS.fullFileBase.alignedMesh, timePoints(1))) ;
    tmp = mesh ;
    tmp.v = laplacian_smooth(tmp.v ,tmp.f, "cotan", [], 0.002) ;
    tmp.vn = per_vertex_normals(tmp.v, tmp.f) ;
    % note: adding extra thickness near the posterior
    tmp.v = tmp.v + 10 * tmp.vn + (amesh.v(:, 1)/ 25) .* tmp.vn;
    trisurf(triangulation(tmp.f, tmp.v), 'edgecolor', 'none')
    axis equal
    % resave as smoothed mesh
    plywrite_with_normals(sprintf(preSmoothMeshfn, timePoints(1)), mesh.f, mesh.v, mesh.vn) ;
    plywrite_with_normals(sprintf(QS.fullFileBase.mesh, timePoints(1)), tmp.f, tmp.v, tmp.vn) ;
end

%% Define QuapSlap object
opts.meshDir = meshDir ;
opts.ilastikOutputAxisOrder = 'xyzc';  % used in APDV alignment
opts.flipy = flipy ;
opts.timeInterval = timeInterval ;
opts.timeUnits = timeUnits ;
opts.spaceUnits = spaceUnits ;
opts.nV = nV ;
opts.nU = nU ;
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


%% APDV ilastik training
% Train on anterior (A), posterior (P), background (B), and 
% dorsal anterior (D) location in different iLastik channels.
% Default order is 1-A, 2-P, 3-bg, 4-D. 
% Perform on pre-stabilized H5s, NOT on stabilized H5s. 
% anteriorChannel, posteriorChannel, and dorsalChannel specify the iLastik
% training channel that is used for each specification.
% Name the h5 file output from iLastik as ..._Probabilities_APDVcoords.h5
% with dataset name /exported_data.
% Train for anterior dorsal (D) only at the first time point, because
% that's the only one that's used. 
% Dorsal anterior for the gut is at the fused site where additional 
% 48YGAL4-expressing muscle-like cells form a seam, against heart tube.
% Posterior is at the rear of the yolk, where the endoderm closes, for 
% apical surface training. 
% Anterior is at the junction of the midgut with the foregut.

%% 3. align_meshes_APDV based on first timepoint APD training
% NOTE: later timepoints are irrelevant here, since we will use anterior
% and posterior training from earlier to extract centerlines.
% Skip if already done.

overwrite_APDVCOMs = false;
overwrite_APDVMeshAlignment = false ;
overwrite_alignAPDVOpts = false;
overwrite_alignedMeshIms = false ;

% Try to load the results to test if we must do calculation 
try
    [rot, trans] = QS.getRotTrans() ;
    [xyzlim_raw, xyzlim, xyzlim_um, xyzlim_buff] = getXYZLims(QS) ;
    assert(exist(QS.fileName.startendPt, 'file') ~= 0)
    redo_alignmesh = false ;
    
    % check all timepoints
    tidx = 1; 
    while ~redo_alignmesh && tidx <= length(xp.fileMeta.timePoints)
        tt = xp.fileMeta.timePoints(tidx) ;
        searchfn = sprintf(QS.fullFileBase.alignedMesh, tt);
        if ~strcmp(searchfn(end-3:end), '.ply')
            searchfn = [searchfn '.ply'] ;
        end
        if ~exist(searchfn, 'file')
            redo_alignmesh = true ;
            disp(['did not find aligned mesh for tp = ' num2str(tt)])
        end
        tidx = tidx + 1;
    end
    
catch
    redo_alignmesh = true ;
end

if redo_alignmesh || overwrite_APDVMeshAlignment || overwrite_APDVCOMs
    disp('Calculating/Overwriting APDVMeshAlignment')
    % OUTPUTS OF THIS SECTION
    % -------
    % xyzlim.txt 
    %   xyzlimits of raw meshes in units of full resolution pixels (ie not
    %   downsampled)
    % xyzlim_APDV.txt 
    %   xyzlimits of rotated and translated meshes in units of full resolution 
    %   pixels (ie not downsampled)
    % xyzlim_APDV_um.txt 
    %   xyz limits of rotated and translated meshes (APDV coord sys) in microns
    % rotation_APDV.txt
    %   rotation matrix to align mesh to APDV frame
    % translation_APDV.txt
    %   translation vector to align mesh to APDV frame
    % xyzlim.txt 
    %   raw bounding box in original frame (not rotated), in full res pixels
    % xyzlim_APDV.txt
    %   bounding box in rotated frame, in full resolution pixels
    % xyzlim_APDV_um.txt
    %   bounding box in rotated frame, in microns
    % apdv_coms_rs.h5
    %   Centers of mass for A,aux_paths_and_colors P, and D in microns in rotated APDV coord system
    % 
    % Notes
    % -----
    % vertices are in units of pixels (at full resolution)
    % To take mesh to rotated + translated mesh in physical units, apply:
    %         xs = mesh.vertex.z ;
    %         ys = mesh.vertex.y ;
    %         zs = mesh.vertex.x ;
    %         vtx_rs = (rot * vtx' + trans)' * resolution
    %         
    save_figs = true ;  % save images of cntrline, etc, along the way
    preview = false ;  % display intermediate results, for debugging
    dorsal_thres = 0.5 ;  % threshold for extracting Dorsal probability cloud 
    buffer = 5 ;  % extra space in meshgrid of centerline extraction, to ensure mesh contained in volume
    plot_buffer = 40;  % extra space in plots, in um
    weight = 0.1;  % for speedup of centerline extraction. Larger is less precise
    normal_step = 0.5 ;  % how far to move normally from ptmatched vtx if a/pcom is not inside mesh
    eps = 0.01 ;  % value for DT outside of mesh in centerline extraction
    meshorder = 'xyz' ;  % ordering of axes in loaded mesh wrt iLastik output
    anteriorChannel = 1;  % which channel of APD training is anterior
    posteriorChannel = 2;  % which channel of APD training is posterior 
    dorsalChannel = 4 ;  % which channel of APD training is dorsal

    clearvars opts
    optsfn = QS.fileName.apdv_options ;
    if exist(optsfn, 'file') && ~overwrite_alignAPDVOpts
        disp('Loading options from disk')
        load(optsfn, 'alignAPDVOpts', 'apdvOpts')
    else
        disp('No alignAPDV_Opts on disk or overwriting, defining')
        apdvOpts.smwindow = 0 ;
        apdvOpts.buffer = buffer ;  
        apdvOpts.plot_buffer = plot_buffer ;
        apdvOpts.anteriorChannel = anteriorChannel ;
        apdvOpts.posteriorChannel = posteriorChannel ;
        apdvOpts.tpref = timePoints(1) ;
        % Can specify separate APDV coordinate system from centerline
        apdvOpts.aProbFileName = QS.fullFileBase.apCenterlineProb ;
        apdvOpts.pProbFileName = QS.fullFileBase.apCenterlineProb ;
        
        
        alignAPDVOpts.weight = weight ;
        alignAPDVOpts.normal_step = normal_step ;
        alignAPDVOpts.eps = eps ;
        alignAPDVOpts.meshorder = meshorder ;
        alignAPDVOpts.anteriorChannel = anteriorChannel ;
        alignAPDVOpts.posteriorChannel = posteriorChannel ;
        alignAPDVOpts.dorsalChannel = dorsalChannel ;
        alignAPDVOpts.dorsal_thres = dorsal_thres ;
        alignAPDVOpts.tref = timePoints(1) ;       
        alignAPDVOpts.fn = fn ;  % filename base
        alignAPDVOpts.aProbFileName = QS.fullFileBase.apdProb ;
        alignAPDVOpts.pProbFileName = QS.fullFileBase.apdProb ;

        % save the options
        save(optsfn, 'alignAPDVOpts', 'apdvOpts')
    end

    % Add script-instance-specific options
    apdvOpts.overwrite = overwrite_APDVCOMs ;  % recompute APDV coms from training
    apdvOpts.check_slices = false ;
    apdvOpts.preview = preview ;
    apdvOpts.preview_com = false ;
        
    %% Compute APDV coordinate system-- this is where rot is defined
    QS.computeAPDVCoords(alignAPDVOpts) ;
    
    % Compute the APD COMs 
    % apdvOpts.aProbFileName = QS.fullFileBase.apdProb ; % filename pattern for the apdv training probabilities
    % apdvOpts.pProbFileName = QS.fullFileBase.apdProb ;
    apdvOpts.posteriorChannel = 2 ;
    [acom_sm, pcom_sm] = QS.computeAPDCOMs(apdvOpts) ;
    
    % Align the meshes APDV & plot them
    alignAPDVOpts.overwrite_ims = overwrite_alignedMeshIms ;  % overwrite images even if centerlines are not overwritten
    alignAPDVOpts.overwrite = overwrite_APDVCOMs || overwrite_APDVMeshAlignment ; % recompute APDV rotation, translation
    QS.alignMeshesAPDV(alignAPDVOpts) ;
    
    
    % NOTE: normals should point INWARD
    tmp = read_ply_mod(sprintf(QS.fullFileBase.alignedMesh, 1)) ;
    trisurf(triangulation(tmp.f, tmp.v), tmp.vn(:, 1), 'edgecolor', 'none')
    axis equal
    title('Do normals point inward?')
    
else
    % Display APDV COMS over time
    acom_sm = h5read(QS.fileName.apdv, '/acom_sm') ;
    pcom_sm = h5read(QS.fileName.apdv, '/pcom_sm') ;
    acoms = h5read(QS.fileName.apdv, '/acom') ;
    pcoms = h5read(QS.fileName.apdv, '/pcom') ;
    dcom = dlmread(QS.fileName.dcom) ;
    for tidx = 1:length(timePoints)
        tp = timePoints(tidx) ;
        % Plot the APDV points
        clf
        plot3(acom_sm(tidx, 1), acom_sm(tidx, 2), acom_sm(tidx, 3), 'ro')
        hold on;
        plot3(acoms(tidx, 1), acoms(tidx, 2), acoms(tidx, 3), 'r.')
        plot3(pcom_sm(tidx, 1), pcom_sm(tidx, 2), pcom_sm(tidx, 3), 'b^')
        plot3(pcoms(tidx, 1), pcoms(tidx, 2), pcoms(tidx, 3), 'b.')
        plot3(dcom(1, 1), dcom(1, 2), dcom(1, 3), 'cs')
        axis equal
        title(['t = ', num2str(tp)]) 
        pause(0.01)
        if tp > 135
            pause(1)
        end
    end
    disp('Already done')
end
disp('done')
clearvars normal_step 


%% Fix vertex normals in alignedMeshes (hack)
% for tt = QS.xp.fileMeta.timePoints ;
%     alignedmeshfn = fullfile(QS.dir.alignedMesh, ...
%         sprintf([QS.fileBase.alignedMesh '.ply'], tt)) ;
%     mm = read_ply_mod(alignedmeshfn) ;
%     xyzrs = mm.v ;
%     vn_rs = mm.vn ;
%     vn_rs(:, 1) = -vn_rs(:, 1) ;
%     vn_rs(:, 3) = -vn_rs(:, 3) ;
%     
%     % Save it
%     if overwrite || ~exist(alignedmeshfn, 'file')
%         disp('Saving the aligned mesh...')
%         disp([' --> ' alignedmeshfn])
%         plywrite_with_normals(alignedmeshfn, mm.f, xyzrs, vn_rs)
%     end
% end

%% MAKE MASKED DATA FOR PRETTY VIDEO ======================================
% Skip if already done
% Generate masks for isolating intensity data of the gut alone, for making
% pretty videos without fiducial markers and other spurious fluorescent
% bits. 
% Save output h5s trained on stabilized h5s from iLastik as 
%   -->   <QS.fileBase.name>_Probabilities_mask3d.h5
% and 
%   -->   <QS.fileBase.name>_Probabilities_maskDorsal.h5
%QS.generateMaskedData()

%% MAKE ORIENTED MASKED DATA FOR PRETTY VIDEO =============================
% Skip if already done
%QS.alignMaskedDataAPDV()

%% PLOT ALL TEXTURED MESHES IN 3D =========================================
% Skip if already done
overwrite = false ;
overwrite_TextureMeshOpts= true ;

% Get limits and create output dir
% Establish texture patch options
metafn = fullfile(QS.dir.texturePatchIm, 'metadat.mat') ;
if ~exist(metafn, 'file') || overwrite_TextureMeshOpts
    [~,~,~,xyzbuff] = QS.getXYZLims() ;
    xyzbuff(:, 1) = xyzbuff(:, 1) - 20 ; 
    xyzbuff(:, 2) = xyzbuff(:, 2) + 20 ; 
    % Define & Save metadata
    metadat.xyzlim = xyzbuff ;                  % xyzlimits
    metadat.reorient_faces = false ;            % if some normals are inverted
    metadat.normal_shift = 15; %QS.normalShift ;             % normal push, in pixels, along normals defined in data XYZ space
    metadat.texture_axis_order = QS.data.axisOrder ;    % texture space sampling
    metadat.smoothing_lambda = 0.002 ;
    
    % Psize is the linear dimension of the grid drawn on each triangular face
    Options = struct() ;
    Options.PSize = 5 ;
    Options.EdgeColor = 'none';
    QS.getRotTrans() ;
    Options.Rotation = QS.APDV.rot ;
    Options.Translation = QS.APDV.trans ;
    Options.Dilation = QS.APDV.resolution ;
    Options.numLayers = [20, 5];  % at layerSpacing=2, numLayers=2 marches ~0.5 um 
    Options.layerSpacing = 0.75 ;
    Options.smoothIter = 1000 ;
    Options.smoothLambda = 0.004 ;
    
    % Save it
    disp('Saving metadat')
    save(metafn, 'metadat', 'Options')
else
    load(metafn, 'metadat', 'Options')
end

% Use first timepoint's intensity limits throughout
% QS.setDataLimits(QS.xp.fileMeta.timePoints(1), 1.0, 99.95)
% For mef2Gal4CAAX
% QS.data.adjustlow = 1000 ;
% QS.data.adjusthigh = 65535 ;

% For UbxMutantAntp8C11 e2 202108161910_e2
QS.data.adjustlow = 0 ;
QS.data.adjusthigh = 99.99 ; 

%% Plot on surface for all TP 
options = metadat ;
options.overwrite = false ;
options.plot_dorsal = false ;
options.plot_ventral = false ;
options.plot_right = false ;
options.plot_left = false ;
options.plot_perspective = true ;
options.channel = [1] ; % if empty, plot all channels

QS.plotSeriesOnSurfaceTexturePatch(options, Options)

%% Now plot meshes for overlay
amprettyOpts = options ;
options.plot_dorsal = true ;
options.plot_ventral = true ;
options.plot_right = false ;
options.plot_left = true ;
options.plot_perspective = true ;
amprettyOpts.normal_shift = -options.normal_shift * QS.APDV.resolution ;
amprettyOpts.smoothing_lambda = 2e-3 ;
amprettyOpts.smoothIter =1000 ;
QS.plotAlignedMeshesPretty(amprettyOpts)

%% EXTRACT CENTERLINES
% Skip if already done
% Note: these just need to be 'reasonable' centerlines for topological
% checks on the orbifold cuts.
exponent = 1.0 ;
res = 4.0 ; 
cntrlineOpts.overwrite = overwrite_centerlines ;     % overwrite previous results
cntrlineOpts.overwrite_ims = overwrite_centerlineIms ;     % overwrite previous results
cntrlineOpts.weight = 0.6 ;              % for speedup of centerline extraction. Larger is less precise
cntrlineOpts.exponent = exponent ;       % how heavily to scale distance transform for speed through voxel
cntrlineOpts.res = res ;                 % resolution of distance tranform grid in which to compute centerlines
cntrlineOpts.preview = false ;           % preview intermediate results
cntrlineOpts.reorient_faces = false ;    % not needed for our well-constructed meshes
cntrlineOpts.dilation = 0 ;              % how many voxels to dilate the segmentation inside/outside before path computation
cntrlineOpts.skipErrors = true ;         % if no path found, skip timepoint
cntrlineOpts.epsilon = 5e-7 ;            % small value for weight in outside region
% Note: this takes about 400s per timepoint for res=2.0
%
QS.extractCenterlineSeries(cntrlineOpts)
disp('done with centerlines')

%% Fix flip in Y for centerlines
% aux_fix_flip

%% Identify anomalies in centerline data
% Skip if already done
idOptions.ssr_thres = 15 ;  % distance of sum squared residuals in um as threshold
idOptions.overwrite = overwrite_idAnomClines ;
QS.generateCleanCntrlines(idOptions) ;


%% Cylinder cut mesh
% Skip if already done
if overwrite_endcapOpts || ~exist(QS.fileName.endcapOptions, 'file')
    
    % THIS WORKS
    % endcapOpts = struct( 'adist_thres', 35, ...  % 20, distance threshold for cutting off anterior in pix
    %             'pdist_thres', 18, ...  % 15-20, distance threshold for cutting off posterior in pix
    %             'aOffset', [-10, -1, -1], ...
    %             'aOffsetRate', [-0.12, -0.08, 0], ...
    %             'aDistRate', [5/120, 120], ... %; -10/60, 60], ...
    %             'aCapMethod', 'ball', ...
    %             'pCapMethod', 'ball') ;   
    
    % TRY OUT CONICAL THRESHOLD
    % 55 will work if 52 is too small.
    endcapOpts = struct( 'adist_thres', 20, ...  % 20, distance threshold for cutting off anterior in pix
        'adist_thres2', 22, ...
        'pdist_thres', 16, ...  % 15-20, distance threshold for cutting off posterior in pix
        'aOffset', [16, -1, -1], ...
        'aOffset2', [-5, -1, -1], ...
        'aOffsetRate', [0.02, -0.09, 0.00], ...
        'aOffsetRate2', [0.0, 0.0, 0.00], ...
        'aDistRate', [0.6/140, 140; -0.2/20, 20; -0.3/8, 8], ...
        'aCapMethod', 'ballCone', ...
        'pCapMethod', 'ball') ;
    QS.setEndcapOptions(endcapOpts) ;
    % Save the options to disk
    QS.saveEndcapOptions() ;
else
    % load endcapOpts
    QS.loadEndcapOptions() ;
    endcapOpts = QS.endcapOptions ;
end

clearvars methodOpts
methodOpts.overwrite = overwrite_endcapOpts ;  % recompute sliced endcaps
methodOpts.save_figs = true ;   % save images of cntrline, etc, along the way
methodOpts.preview = false  ;     % display intermediate results
QS.sliceMeshEndcaps(endcapOpts, methodOpts) ;

%% Clean Cylinder Meshes
% May skip if already done
cleanCylOptions.overwrite = overwrite_cleanCylMesh ;
cleanCylOptions.save_ims = true ;
QS.cleanCylMeshes(cleanCylOptions)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORBIFOLD -> begin populating Qs.dir.mesh/gridCoords_nUXXXX_nVXXXX/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
% Overwriting options
overwrite_pullbacks = false ;
overwrite_cutMesh = false ;
overwrite_spcutMesh = false ;
% Plotting params
washout2d = 0.5 ;

%% Iterate Through Time Points to Create CutMeshes ========================
% Skip if already done
% outcutfn = fullfile(cutFolder, 'cutPaths_%06d.txt') ;

% First make cutMeshes
for tt = xp.fileMeta.timePoints(1:end)
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    tidx = xp.tIdx(tt);
    
    % Load the data for the current time point ------------------------
    QS.setTime(tt) ;
    
    %----------------------------------------------------------------------
    % Create the Cut Mesh
    %----------------------------------------------------------------------
    cutMeshfn = sprintf(QS.fullFileBase.cutMesh, tt) ;
    cutPathfn = sprintf(QS.fullFileBase.cutPath, tt) ;
    if ~exist(cutMeshfn, 'file') || ~exist(cutPathfn, 'file') || overwrite_cutMesh
        if exist(cutMeshfn, 'file')
            disp('Overwriting cutMesh...') ;
        else
            disp('cutMesh not found on disk. Generating cutMesh... ');
        end
        options = struct() ;
        options.preview = false ;
        QS.generateCurrentCutMesh(options)
        disp('Saving cutP image')
        % Plot the cutPath (cutP) in 3D
        QS.plotCutPath(QS.currentMesh.cutMesh, QS.currentMesh.cutPath)
        % compute_pullback = true ;
    else
        fprintf('Loading Cut Mesh from disk... ');
        QS.loadCurrentCutMesh()
        % compute_pullback = ~isempty(QS.currentMesh.cutPath) ;
        
        cutfn = sprintf( fullfile(fullfile(QS.dir.cutMesh, 'images'), ...
            [QS.fileBase.name, '_cut.png']), tt ) ;
        if ~exist(cutfn, 'file')
            QS.plotCutPath(QS.currentMesh.cutMesh, QS.currentMesh.cutPath) ;
        end
    end
end

%% Iterate Through Time Points to Create Pullbacks ========================
% Now make Pullbacks
for tt = xp.fileMeta.timePoints(1:end)
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    tidx = xp.tIdx(tt);
    
    % Load the data for the current time point ------------------------
    % HACK: adjustment for antpGal4 CAAX
    QS.data.adjustlow = [80, 215] ;
    QS.data.adjusthigh = [22000, 10000] ;
    QS.setTime(tt) ;
    
    %----------------------------------------------------------------------
    % Create the Cut Mesh
    %----------------------------------------------------------------------
    cutMeshfn = sprintf(QS.fullFileBase.cutMesh, tt) ;
    cutPathfn = sprintf(QS.fullFileBase.cutPath, tt) ;
    if ~exist(cutMeshfn, 'file') || ~exist(cutPathfn, 'file') || overwrite_cutMesh
        if exist(cutMeshfn, 'file')
            disp('Overwriting cutMesh...') ;
        else
            disp('cutMesh not found on disk. Generating cutMesh... ');
        end
        options = struct() ;
        options.preview = false ;
        QS.generateCurrentCutMesh(options)
        disp('Saving cutP image')
        % Plot the cutPath (cutP) in 3D
        QS.plotCutPath(QS.currentMesh.cutMesh, QS.currentMesh.cutPath)
        compute_pullback = true ;
    else
        fprintf('Loading Cut Mesh from disk... ');
        QS.loadCurrentCutMesh()
        compute_pullback = ~isempty(QS.currentMesh.cutPath) ;
        
        cutfn = sprintf( fullfile(fullfile(QS.dir.cutMesh, 'images'), ...
            [QS.fileBase.name, '_cut.png']), tt ) ;
        if ~exist(cutfn, 'file')
            QS.plotCutPath(QS.currentMesh.cutMesh, QS.currentMesh.cutPath) ;
        end
    end
    
    spcutMeshOptions.overwrite = overwrite_spcutMesh ;
    spcutMeshOptions.save_phi0patch = false ;
    spcutMeshOptions.iterative_phi0 = true ;
    spcutMeshOptions.smoothingMethod = 'none' ;
    
    QS.plotting.preview = false ;
    QS.generateCurrentSPCutMesh([], spcutMeshOptions) ;
    
    % Compute the pullback if the cutMesh is ok
    if compute_pullback 
        pbOptions.overwrite = true ;
        pbOptions.generate_uv = true ;
        pbOptions.generate_uphi = false ;
        pbOptions.generate_relaxed = true ;
        pbOptions.numLayers = {[1,1], [0,45]} ;
        pbOptions.layerSpacing = 0.8 ;
        pbOptions.falseColors = {[1,0,0], [0,1,1]} ;
        
        % For AntpGAL4 CAAX
        pbOptions.smoothIter = 1000 ;
        pbOptions.smoothLambda = 0.004 ;
        QS.generateCurrentPullbacks([], [], [], pbOptions) ;
    else
        disp('Skipping computation of pullback')
    end
    clear Options IV
        
    % Save SMArr2D (vertex positions in the 2D pullback) -----------------
    % disp(['Saving meshStack to disk: ' mstckfn])
    % save(mstckfn, 'meshStack') ;
    % 
    % %% Save SMArr2D (vertex positions in the 2D pullback) -----------------
    % disp(['Saving spmeshStack to disk: ' spmstckfn])
    % if generate_sphi_coord
    %     save(spmstckfn, 'spmeshStack') ;
    % end
end
clearvars t cutP spcutMesh spvn3d ss pp uv tileCount slin plin plotfn
clearvars IVloaded IV uphi sphi
disp('Done with generating spcutMeshes and cutMeshes')
