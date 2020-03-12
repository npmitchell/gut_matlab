%% MASTER GUT PULLBACK -- FIXED SAMPLE
% NPMitchell
%
% This is a pipeline to take the surface of the fixed Drosophila gut and
% map to the plane

% temporary path def
cd /mnt/data/antibodies_lightsheet/202002201800_w48YGal4klar_DAPI_abdA100_exp2_2mW_25x_1p4um_ms568/Time4views_60sec_1p4um_25x_2mW_exp2/data

%% IMSANE SETUP FOR DETECTOR
% 
% cd('/mnt/data/code/imsane_for_git/imsane/')
% Run setup
% cd(projectDir)

%% INITIALIZE ImSAnE PROJECT ==============================================
%
% We start by clearing the memory and closing all figures
clear; close all; clc;

% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.

dataDir    =  cd; 
projectDir = dataDir ;
% [ projectDir, ~, ~ ] = fileparts(matlab.desktop.editor.getActiveFilename); 

%% ADD PATHS
addpath('/mnt/data/code/gut_matlab/addpath_recurse/')
addpath_recurse('/mnt/crunch/djcislo/MATLAB/CGAL_Code/')
addpath_recurse('/mnt/data/code/gptoolbox/')
addpath_recurse('/mnt/data/code/imsane_for_git/imsane/')

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
fileMeta.nChannels          = 1;
fileMeta.timePoints         = 0 ;
fileMeta.stackResolution    = [.2619 .2619 .2619];
fileMeta.swapZT             = 1;

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
expMeta.channelsUsed        = 0;
expMeta.channelColor        = 0;
expMeta.description         = 'Apical membrane in Drosophila gut';
expMeta.dynamicSurface      = 0;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.fitTime             = fileMeta.timePoints(first_tp);
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';

%% 32 to 16 BIT CONVERSION
% Check that data is 16 bit. If not, convert to 16bit
pc2use = 99.999;
for channel_check = expMeta.channelsUsed
    fullFileName = [sprintf(fn, channel_check) '.tif'] ;
    info = imfinfo(fullFileName) ;
    full16fn = [sprintf(file16name, channel_check) '.tif'] ;
    bitDepth = info.BitDepth ;

    if (bitDepth == 32) && ~isfile(full16fn)

        disp([fullFileName ' is not 16bit, converting...'])

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
        scalemax = double(prctile(tmpA, pc2use)) ;
        scalemin = double(min(tmpA(tmpA > 0))) ;

        % Note to self to check scale:
        % imagesc(squeeze(A(:, 300, :)) ; 
        % histogram(A(:))

        data = readSingleTiff(fullFileName);
        im2 = mat2gray(data,[scalemin scalemax]);
        im2 = uint16(2^16*im2);
        imSize = size(im2);

        disp(['Saving 16bit volume to ' full16fn])
        for z = 1 : imSize(3)
            imwrite(im2(:,:,z),full16fn,'tiff','Compression','none','WriteMode','append');
        end
        disp('done saving 16bit volume')

        fn = file16name ;
        swapZT = 1 ;

    elseif isfile(full16fn)
        % the data is already 16 bit, so we're good
        fullFileName = [sprintf(fn, channel_check) '.tif'] ;
        disp([fullFileName ' has been converted.'])

        fn = file16name ;
        swapZT = 1 ;
    else
        disp('File is 16bit.')
    end
end

%% INSTANTIATE EXPERIMENT CLASS
% Now set the meta data in the experiment.
fileMeta.filenameFormat = [ fn '.tif' ] ;
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();

%%%%%%%%%%%%%%%%%%%%%%%%%
% Other options
%%%%%%%%%%%%%%%%%%%%%%%%%
mlxprogram = 'surface_rm_resample20k_reconstruct_LS3_1p2pc_ssfactor4.mlx';
msls_axis_order = 'yxzc';
% Mesh marching options
normal_step = 10;

% Define the surface detection parameters
channel = 2;
foreGroundChannel = 2;
ssfactor = 4;
niter = 25 ;
niter0 = 115 ;
ofn_smoothply = 'mesh_' ;
ofn_ply = 'mesh_ms_' ; 
ofn_ls = 'msls_' ;
ms_scriptDir = '/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper' ;
lambda1 = 1 ;
lambda2 = 1 ;
exit_thres = 0.00001 ;
smoothing = 0.1 ;
nu = 0.0 ;
pre_nu = -5 ;
pre_smoothing = 1 ;
post_nu = 2;
post_smoothing = 4 ;
radius_guess = 10 ;
center_guess = 'empty_string' ;

% Name the output mesh directory ------------------------------------------
if projectDir(end) ~= filesep
    projectDir = [projectDir filesep];
end
mslsDir = fullfile(projectDir, 'msls_output');

%% LOAD THE FIRST TIME POINT ==============================================
xp.loadTime(xp.fileMeta.timePoints(first_tp));
xp.rescaleStackToUnitAspect();

%% DETECT THE SURFACE =====================================================
% Surface detection parameters --------------------------------------------
detectOptions = struct('channel', 1, ...
            'ssfactor', ssfactor,... % subsampling factor: downsampling of raw data
            'niter', 100, ... % how many iterations before exit if no convergence
            'niter0', 100, ... % how many iterations before exit if no convergence for first timepoint
            'lambda1', 1, ...  % lambda1/lambda2 decides weight of inclusion/exclusion of interior/exterior
            'lambda2', 1, ...  % lambda1/lambda2 decides weight of inclusion/exclusion of interior/exterior
            'nu', nu, ... % float: how many pressure (dilation/erosion) steps per iteration
            'smoothing', smoothing,... % float: how many smoothing steps per iteration (can be <1)
            'post_nu', post_nu, ... % how many iterations to dilate (if positive) or erode (if negative) after convergence
            'post_smoothing', post_smoothing,... % how many iterations of smoothing after convergence
            'exit_thres', 1e-6, ... % convergence threshold: maximum difference between subsequent level sets upon which to exit algorithm ('close enough')
            'foreGroundChannel',foreGroundChannel, ... % the index of the first dimension of the 4d input data (if 4d)
            'fileName', sprintf( fn, xp.currentTime ), ... % the filename of h5 to train on
            'mslsDir', mslsDir, ...  % the directory for all output data/images
            'ofn_ls', ofn_ls, ...  % the output filename for level sets
            'ofn_ply', ofn_ply, ... % the output filename for PLY files
            'ms_scriptDir', ms_scriptDir, ... % the directory containing run_morphsnakes.py
            'timepoint', 0, ... % which timepoint in the data to consider
            'zdim', 2, ... % Which dimension is the z dimension
            'pre_nu', pre_nu, ... % number of dilation/erosion passes for positive/negative values
            'pre_smoothing', pre_smoothing, ... % number of smoothing passes before running MS
            'ofn_smoothply', ofn_smoothply,... % the output file name (not including path directory)
            'mlxprogram', mlxprogram, ... % the name of the mlx program to use to smooth the results. Note that if mesh_from_pointcloud==true, should take obj as input and mesh as output.
            'init_ls_fn', 'mesh_initguess', ... % the name of the initial level set to load, if any
            'run_full_dataset', false, ... % run MS on a time series, not just one file
            'radius_guess', radius_guess, ... % radius of the initial guess sphere
            'dset_name', 'exported_data', ... % the name of the dataset to load from h5
            'save', true, ... % whether to save intermediate results
            'center_guess', 'empty_string', ... % xyz of the initial guess sphere ;
            'plot_mesh3d', false, ...  % if save is true, plot intermediate results in 3d 
            'dtype', 'h5', ... % h5 or npy: use hdf5 or numpy file format for input and output ls
            'mask', 'none', ... % filename for mask to apply before running MS
            'mesh_from_pointcloud', false, ... % use a pointcloud from the marching cubes algorithm rather than a mesh to create smoothed mesh
            'prob_searchstr', '_Probabilities.h5', ... % if dataset mode, what string to seek for loading all probabilities in data directory (glob datadir/*searchstr)
            'imsaneaxisorder', 'xyzc', ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc
            'preilastikaxisorder', 'xyzc', ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
            'ilastikaxisorder', 'xyzc', ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
            'include_boundary_faces', true) ; % keep faces along the boundaries of the data volume if true

% Set detect options ------------------------------------------------------
xp.setDetectOptions( detectOptions );

% clear msls_exten imwriteOptions saveDir
% clear channel foreGroundChannel
% clear niter niter0 lambda1 lambda2
% clear exit_thres smoothing nu
% clear post_nu post_smoothing

%% CREATE THE SUBSAMPLED H5 FILE FOR INPUT TO ILASTIK =====================
% skip if already done
for c = xp.expMeta.channelsUsed
    for t = xp.fileMeta.timePoints
        if ~exist(fullfile(projectDir, [sprintf(sprintf(fn, t), c) '.h5']), 'file')
            disp(['Did not find file: ', fullfile(projectDir, [sprintf(fn, t) '.h5'])])
            xp.loadTime(t);
            xp.rescaleStackToUnitAspect();
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
end    
disp('Open with ilastik if not already done')

%% TRAIN DATA IN ILASTIK TO IDENTIFY APICAL/YOLK ==========================
% open ilastik, train until probabilities and uncertainty are satisfactory

%% Create MorphoSnakesLevelSet from the Probabilities from ilastik ========
xp.detectSurface();
fileMeta = xp.fileMeta ;
