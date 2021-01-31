%% MASTER GUT PULLBACK -- FIXED SAMPLE
% NPMitchell
%
% This is a pipeline to take the surface of the fixed Drosophila gut and
% map to the plane

% We start by clearing the memory and closing all figures
clear; close all; clc;

% temporary path def
cd /mnt/data/antibodies_lightsheet/202002201800_w48YGal4klar_DAPI_abdA100_exp2_2mW_25x_1p4um_ms568/Time4views_60sec_1p4um_25x_2mW_exp2/data

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
    nChannels = 1 ;
    channelsUsed = 1 ;
    timePoints = 0 ;
    ssfactor = 4 ;
    % whether the data is stored inverted relative to real position
    flipy = true ; 
    timeInterval = 1 ;  % physical interval between timepoints
    timeUnits = 'min' ; % physical unit of time between timepoints
    scale = 0.02 ;      % scale for conversion to 16 bit
    file32Base = 'TP%d_Ch1_Ill0_Ang0,45,90,135,180,225,270,315.tif'; 
    % file32Base = 'TP%d_Ch0_Ill0_Ang0,60,120,180,240,300.tif'; 
    fn = 'fixed_sample_c1';
    spaceUnits = '$\mu$m';  % microns as $\mu$m
    
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
    t0_for_phi0 = masterSettings.t0_for_phi0 ;
    nU = masterSettings.nU ;
    nV = masterSettings.nV ;


    scale = masterSettings.scale ;      % scale for conversion to 16 bit
    file32Base = masterSettings.file32Base ; 
    fn = masterSettings.fn ;
    fn_prestab = masterSettings.fn_prestab ;
    set_preilastikaxisorder = masterSettings.set_preilastikaxisorder ;
end
dir16bit = fullfile(dataDir, 'deconvolved_16bit') ;

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
fileMeta.nChannels          = 1;
fileMeta.timePoints         = 1 ;
fileMeta.stackResolution    = [.2619 .2619 .2619];
fileMeta.swapZT             = 0;

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
expMeta.channelsUsed        = 1;
expMeta.channelColor        = 1;
expMeta.description         = 'Drosophila gut';
expMeta.dynamicSurface      = 0;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.fitTime             = fileMeta.timePoints(first_tp);
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';

%% 32 to 16 BIT CONVERSION
% Check that data is 16 bit. If not, convert to 16bit
pc2use = 99.9;
scalemethod = 'user-defined' ; % use either 'user-defined' or 'prctile'
scalemax = 0.05 ;

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
        waitfor(gcf)
        histogram(A(:))
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

        fn = file16name ;

    elseif isfile(full16fn)
        % the data is already 16 bit, so we're good
        fullFileName = [sprintf(fn, channel_check) '.tif'] ;
        disp([fullFileName ' has been converted to 16bit: ' full16fn ])

        fn = file16name ;
    else
        disp('File is 16bit.')
    end
end

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
msls_axis_order = 'yxzc';
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

    % Name the output mesh directory ------------------------------------------
    if projectDir(end) ~= filesep
        projectDir = [projectDir filesep];
    end
    mslsDir = fullfile(projectDir, 'msls_output');

    detectOptions = struct('channel', 1, ...
            'ssfactor', ssfactor,... % subsampling factor: downsampling of raw data
            'niter', niter, ... % how many iterations before exit if no convergence
            'niter0', niter0, ... % how many iterations before exit if no convergence for first timepoint
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
            'physicalaxisorder', 'xyzc', ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc
            'preilastikaxisorder', 'xyzc', ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
            'ilastikaxisorder', 'xyzc', ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
            'include_boundary_faces', true,... % keep faces along the boundaries of the data volume if true
            'smooth_with_matlab', -1) ;
    
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
h5fn = fullfile(projectDir, [sprintf(sprintf(fn, t), c) '.h5']) ;
for c = xp.expMeta.channelsUsed
    for t = xp.fileMeta.timePoints
        if ~exist(h5fn, 'file')
            disp(['Did not find file: ', fullfile(projectDir, [sprintf(fn, t) '.h5'])])
            
            % Only load and rescale if multiple timepoints/channels
            if length(xp.expMeta.channelsUsed) > 1 || length(xp.fileMeta.timePoints) > 1
                xp.loadTime(t);
                xp.rescaleStackToUnitAspect();
            end
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


%%% Example python run
% python /mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/run_morphsnakes.py \
%  -i fixed_sample_c1_Probabilities.h5 -o \
%  /mnt/data/antibodies_lightsheet/202002201800_w48YGal4klar_DAPI_abdA100_exp2_2mW_25x_1p4um_ms568/Time4views_60sec_1p4um_25x_2mW_exp2/data/msls_output \
%  -prenu -8 -presmooth 1 -ofn_ply mesh_ms_000000.ply -ofn_ls msls_000000.h5 -l1 1 -l2 1 -nu 0.5 -postnu 1 -channel 1 -smooth 4 \
%  -postsmooth 4 -exit 0.000001000 -channel 1 -dtype h5 -permute xyzc -ss 4 -include_boundary_faces -rad0 10 \
%  -init_ls /mnt/data/antibodies_lightsheet/202002201800_w48YGal4klar_DAPI_abdA100_exp2_2mW_25x_1p4um_ms568/Time4views_60sec_1p4um_25x_2mW_exp2/data/msls_output/mesh_initguess.h5 \
%  -n 20 -save
 
 python /mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/run_morphsnakes.py \
    -i fixed_sample_c1_Probabilities.h5 -o  \
    /mnt/data/antibodies_lightsheet/202002201800_w48YGal4klar_DAPI_abdA100_exp2_2mW_25x_1p4um_ms568/Time4views_60sec_1p4um_25x_2mW_exp2/data/msls_output \
    -prenu -6 -presmooth 1 -ofn_ply msls_ms_000000.ply -ofn_ls msls_000000.h5  \
    -l1 1 -l2 1 -nu 0.5 -postnu 2 -channel 1 -smooth 4  -postsmooth 20 \
    -exit 0.000001000 -channel 1 -dtype h5 -permute xyzc \
    -ss 4 -include_boundary_faces -rad0 10  \
    -init_ls /mnt/data/antibodies_lightsheet/202002201800_w48YGal4klar_DAPI_abdA100_exp2_2mW_25x_1p4um_ms568/Time4views_60sec_1p4um_25x_2mW_exp2/data/msls_output/mesh_initguess.h5 \
    -n 15

%% Overwriting/Computing clean cylinderMesh
% Load the cylinder mesh
cylmeshfn = sprintf( cylinderMeshBase, t ) ;
mesh = read_ply_mod( cylmeshfn );
% Clean the mesh, including averaging normals onto vertices using
% angle-weighted scheme
mesh = cleanCylMesh(mesh) ;   
[adIDx, pdIDx] = aux_adjust_dIDx(mesh, t, dpFile, ADBase, PDBase, cylinderMeshCleanBase, outadIDxfn, outpdIDxfn, xp) ;

%% Save the 3d cut mesh with new indices
% This is saving the cylinder meshes with no ears. Also adIDx
% is saved in h5file.
plywrite_with_normals(mesh3dfn, mesh.f, mesh.v, mesh.vn)
% Save adIDx with new indices
save_to_h5(outadIDxfn, ['/' sprintf('%06d', t) ], adIDx, ['adIDx for t=' num2str(t) ' already exists'])
% Save pdIDx with new indices
save_to_h5(outpdIDxfn, ['/' sprintf('%06d', t) ], pdIDx, ['pdIDx for t=' num2str(t) ' already exists'])
disp('done with cylindermesh cleaning')

 % View results --------------------------------------------------------
mesh3dfigfn = sprintf( cylinderMeshCleanFigBase, t ) ;
if save_ims
    aux_plot_cleanCylMesh
end

%% Cut the mesh
cutOptions.method = 'fastest' ;
disp(['Cutting mesh using method ' cutOptions.method])
cutMesh = cylinderCutMesh( mesh.f, mesh.v, mesh.vn, adIDx, pdIDx, cutOptions );
cutP = cutMesh.pathPairs(:, 1) ;
adIDx = cutP(1) ;
pdIDx = cutP(end) ;
prevTw = twist(mesh.v(cutP, :), cntrlines{t}) ;
compute_pullback = true ;


disp(['Evolving mesh along normal shift for pullback images: shift=' num2str(normal_shift)])
% Displace normally ---------------------------------------------------
cutMesh.v = cutMesh.v + cutMesh.vn * normal_shift ;

% todo !!!!
% if 
% else
%     % Load the cutMesh
%     load(sprintf(cutMeshBase, t), 'cutMesh')
% end

% View results --------------------------------------------------------
% 
% patch( 'Faces', cutMesh.f, 'Vertices', cutMesh.u, ...
%     'FaceVertexCData', cutMesh.v(:,3), 'FaceColor', 'interp', ...
%     'EdgeColor', 'k' );
% 
% hold on
% 
% cornerColors = [ 1 0 0; 1 0 1; 0 1 1; 0 1 0 ];
% corners = [ adIDx cutMesh.pathPairs(1,2), ...
%     cutMesh.pathPairs(end,1) pdIDx ];
% 
% scatter( cutMesh.u( corners, 1 ), cutMesh.u( corners, 2 ), [], ...
%     cornerColors, 'filled' );
% 
% hold off
% 
% axis equal
% 
% clear cornerColors corners

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate s,phi coord system for rotated,scaled mesh (rs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Establishing s,phi coord system\n');

if ~exist(sprintf(spcutMeshBase, t), 'file') || overwrite_spcutMesh
    if overwrite_spcutMesh
        disp('Overwriting spcutMesh...')
    else
        disp('spcutMesh not on disk. Generating ...')
    end

    % Transform from u,v coordinates to s, phi coordinates
    % [scoords, phicoords] = generateSPhiFromUV();

    %----------------------------------------------------------------------
    % Generate tiled orbifold triangulation
    %----------------------------------------------------------------------
    tileCount = [1 1];  % how many above, how many below
    cutMeshrs = cutMesh;
    % Rotate and translate TV3D
    cutMeshrs.v = ((rot * cutMesh.v')' + trans) * resolution ;
    cutMeshrs.vn = (rot * cutMesh.vn')' ;
    [ ~, ~, TV3D, TVN3D ] = tileAnnularCutMesh( cutMesh, tileCount );
    [ TF, TV2D, TV3Drs ] = tileAnnularCutMesh( cutMeshrs, tileCount );

    %----------------------------------------------------------------------
    % Calculate abbreviated centerline from cutMesh boundaries
    %----------------------------------------------------------------------
    % Load centerline in raw units
    % cntrfn = sprintf(cntrsFileName, t) ;
    % cline = dlmread(cntrfn, ',') ;
    cline = cntrlines_rs{t} ; 
    ss = cline(:, 1) ;
    cline = cline(:, 2:end) ;
    disp('Finding relevant segment of centerline')
    [cseg, acID, pcID, bdLeft, bdRight] = centerlineSegmentFromCutMesh(cline, TF, TV2D, TV3Drs) ;

    %----------------------------------------------------------------------
    % Generate surface curves of constant s
    %----------------------------------------------------------------------
    % For lines of constant phi
    disp('Creating crude uv curves with du=const to define uspace by ds(u)')
    % Make grid
    eps = 1e-14 ;
    uspace0 = linspace( eps, cutMesh.umax - eps, nU )' ;
    vspace = linspace( eps, 1-eps, nV )' ;

    disp('Casting crude (equal dU) points into 3D...')
    crude_ringpath_ss = ringpathsGridSampling(uspace0, vspace, TF, TV2D, TV3Drs) ;

    % Resample crude_ringpath_ds made from uspace0 (equal du, not equal ds_3D in u direction)
    [uspace, eq_ringpath_ss] = equidistantSampling1D(linspace(0, 1, nU)', crude_ringpath_ss, nU, 'linear') ;
    % ensure that uspace is nU x 1, not 1 x nU
    uspace = reshape(uspace, [nU, 1]) ; 
    % hedge the first and last point to avoid NaNs
    eps = 1e-13 ;
    uspace(1) = uspace(1) + eps ;
    uspace(end) = uspace(end) - eps ;
    clearvars dsuphi curves3d uspace0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Casting resampled points into 3D (approx equal ds_3D in u dir, but variable ds_3D in v dir)...')
    % NOTE: first dimension indexes u, second indexes v
    curves3d = zeros(nU, nV, 3) ;  % in units of um
    for kk = 1:nU
        if mod(kk, 50) == 0
            disp(['u = ' num2str(kk / nU)])
        end
        uv = [cutMesh.umax * uspace(kk) * ones(size(vspace)), vspace] ;
        curves3d(kk, :, :) = interpolate2Dpts_3Dmesh(TF, TV2D, TV3Drs, uv) ;
    end 

    % Check the 3d curves 
    if preview
        figure ; hold on;
        for kk = 1:nU
            plot3(curves3d(kk, :, 1), curves3d(kk, :, 2), curves3d(kk, :, 3), '.') 
        end
        axis equal
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Compute s(u) and radius(u) for "uniform"--> evenly sample each DV hoop (0,1) so ds_3D=const \n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Resample at evenly spaced dphi in embedding space (rs, in um)
    fprintf('Resampling curves...\n')
    c3d_dsv = zeros(size(curves3d)) ;  % in units of um
    for i=1:nU
        % Note: no need to add the first point to the curve
        % since the endpoints already match exactly in 3d and
        % curvspace gives a curve with points on either
        % endpoint (corresponding to the same 3d location).
        c3d_dsv(i, :, :) = resampleCurvReplaceNaNs(squeeze(curves3d(i, :, :)), nV, true) ;
        if vecnorm(squeeze(c3d_dsv(i, 1, :)) - squeeze(c3d_dsv(i, end, :))) > 1e-7
            error('endpoints do not join! Exiting')
        end

        % Visualization for Troubleshooting:
        % triplot(TF, TV2D(:, 1), TV2D(:, 2))
        % hold on;
        % plot(uv(:, 1), uv(:, 2), '.')
    end

    % Check the 3d curves 
    if preview
        figure ; hold on;
        for kk = 1:nU
            plot3(c3d_dsv(kk, :, 1), c3d_dsv(kk, :, 2), c3d_dsv(kk, :, 3), '.') 
        end
        axis equal
    end

    fprintf('Finding s(u) and r(u) of resampled "uniform" c3ds [uniform ds in V dir]...\n')
    % mcline is the resampled centerline, with mss
    % avgpts is the raw Nx3 averaged hoops, with avgpts_ss
    [mss, mcline, radii_from_mean_uniform_rs, avgpts_ss, avgpts] = srFromDVCurves(c3d_dsv) ;

    % Used to find radius using original centerline
    % [ssv, radii, avgpts, cids] = srFromDVCurvesGivenCenterline(ss, cline, c3ds) ;
    % Could operate just on the centerline segment
    cseg_ss = ss(acID:pcID) ;
    % [ssv, radii, avgpts, cids] = srFromDVCurves(cseg_ss, cseg, c3ds) ;
    % 
    % Adjust the centerline indices to index into the full
    % centerline. Note that cseg_ss already does this for ss.
    % cids = cids + acID ;

    % Plot new centerline
    aux_plot_clineDVhoop(avgpts, avgpts_ss, cseg, cline, cseg_ss, curves3d, xyzlim, clineDVhoopFigBase, t)

    % Optional: clean curve with polynomial and point match
    % avgpts onto cleaned curve. Skipping for later.

    % Compute ringpath_ss, the mean distance traveled from one
    % line of constant u to the next
    disp('Computing ringpath_ss in "uniform" resampling (equal ds along DV)...')
    % The distance from one hoop to another is the
    % difference in position from (u_i, v_i) to (u_{i+1}, v_i).
    dsuphi = reshape(vecnorm(diff(c3d_dsv), 2, 3), [nU-1, nV]) ;
    ringpath_ds = nanmean(dsuphi, 2) ;
    ringpath_ss = cumsum([0; ringpath_ds]) ;
    clearvars dsuphi ringpath_ds

    % Save new centerline in rotated translated units
    fn = sprintf(clineDVhoopBase, t) ;
    disp(['Saving new centerline to ' fn])
    save(fn, 'mss', 'mcline', 'avgpts', 'avgpts_ss')

    % Note: radii_from_mean_uniform_rs is the radius of 
    % interpolated hoops, not the actual points

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Done making new centerline using uniformly sampled hoops\n') ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Create new3d, the regridded pts at UV, moved to sphi')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    onesUV = ones(nU, nV) ;
    uu = uspace * cutMesh.umax .* onesUV ;
    vv = (vspace .* onesUV')' ;
    uv = [uu(:), vv(:)] ;
    % Note: here interpolate uv in the TV2D coord system, then
    % use uphi as the actual 2D coordinates for these vertices
    % NOTE: unlike curves3d, new3d is NOT rotated/translated/scaled
    new3d = interpolate2Dpts_3Dmesh(TF, TV2D, TV3D, uv) ;

    IVloaded = false ;
    if t == xp.fileMeta.timePoints(1)
        % Store for next timepoint
        phiv = (vspace .* ones(nU, nV))' ;
        phi0s = zeros(size(uspace)) ;
        phi0_fit = phi0s ;
    else
        % Load previous sphi vertices in 3d 
        plotfn = sprintf(phi0fitBase, t, 0);
        if strcmp(phi_method, '3dcurves')
            % Load the previous spcutMesh and call it prev3d_sphi
            % Also note the previous spcutMesh pullback image's fn
            tmp = load(sprintf(spcutMeshBase, ...
                xp.fileMeta.timePoints(tidx-1)), 'spcutMesh') ;
            prevf = tmp.spcutMesh.f ;
            prev3d_sphi = reshape(tmp.spcutMesh.v, [nU, nV, 3]) ; 
            imfn_sp_prev = sprintf( ...
                fullfile([imFolder_sp, '/', fileNameBase, '.tif']), ...
                xp.fileMeta.timePoints(tidx-1) ) ;

            % fit the shifts in the y direction
            dmyk = 0 ;
            phi0_fit = zeros(size(uspace)) ;
            phi0s = zeros(size(uspace)) ;
            phi0_fit_kk = 1 ; % for first pass                
            phiv_kk = (vspace .* ones(nU, nV))' ;
            ensureDir([sphiDir, '/phi0_correction/'])
            while any(phi0_fit_kk > 0.002) && dmyk < 6
                disp(['Iteration ' num2str(dmyk)])
                plotfn = sprintf(phi0fitBase, t, dmyk);

                % Will we save check pullbacks to preview the algo?
                if save_phi0patch
                    patchImFn = sprintf( ...
                        fullfile(sphiDir, 'phi0_correction', [fileNameBase, '_prephi0_' num2str(dmyk) '.tif']), ...
                        xp.fileMeta.timePoints(tidx-1) )  ;
                    geomImFn = sprintf( ...
                        fullfile(sphiDir, 'phi0_correction', ['3d' fileNameBase '_prephi0_' num2str(dmyk) '.tif']), ...
                        xp.fileMeta.timePoints(tidx-1) )  ;

                    % Load the intensity data for this timepoint
                    if ~IVloaded
                        % (3D data for coloring mesh pullback)
                        xp.loadTime(t);
                        xp.rescaleStackToUnitAspect();

                        % Raw stack data
                        IV = xp.stack.image.apply();
                        IV = imadjustn(IV{1});         
                        IVloaded = true ;
                    end

                    % Texture patch options
                    Options.PSize = 5;
                    Options.EdgeColor = 'none';
                    % Texture image options
                    Options.imSize = ceil( 1000 .* [ 1 a_fixed ] );
                    Options.yLim = [0 1];

                    % Roll options into a struct
                    patchOpts.patchImFn = patchImFn ;
                    patchOpts.imfn_sp_prev = imfn_sp_prev ;
                    patchOpts.IV = IV ;
                    patchOpts.ringpath_ss = ringpath_ss ;
                    patchOpts.Options = Options ;
                    patchOpts.v3d = new3d ;
                else
                    patchOpts = [] ;
                end

                % Minimize difference in DV hoop positions wrt
                % previous pullback mesh                        
                [phi0_fit_kk, phi0s_kk] = fitPhiOffsetsFromPrevMesh(TF, TV2D, TV3D, ...
                    uspace * cutMesh.umax, phiv_kk, prev3d_sphi, -0.45, 0.45, ...
                    save_ims, plotfn, save_phi0patch, preview, patchOpts) ;

                % Update the result
                dmyk = dmyk + 1;
                phi0_fit = phi0_fit + phi0_fit_kk ;
                phi0s = phi0s + phi0s_kk ;
                phiv_kk = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;


                % plot mesh colored by the phase phi 
                % previous timepoint
                xtmp = prev3d_sphi(:, :, 1) ;
                ytmp = prev3d_sphi(:, :, 2) ;
                ztmp = prev3d_sphi(:, :, 3) ;
                phitmp = (vspace .* ones(nU, nV))' ;
                colormap parula ;
                % cmap = parula ;
                % colors = cmap(max(1, uint8(colortmp(:) * length(parula))), :) ;
                trisurf(prevf, xtmp(:), ytmp(:), ztmp(:), phitmp(:), ...
                    'FaceColor', 'interp',...
                    'EdgeColor', 'none', 'FaceAlpha', 0.25)
                axis equal
                % freezeColors

                % before fitting
                hold on;
                pe0 = find(phitmp(:) < 1e-4 | phitmp(:) > 0.99) ;
                plot3(new3d(pe0, 1), new3d(pe0, 2), ...
                    new3d(pe0, 3), '.')
                % colormap copper
                % trimesh(prevf, new3d(inds, 1), new3d(:, 2), new3d(:, 3),...
                %     phitmp(:), 'FaceColor', 'interp', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
                % freezeColors

                % after fitting
                hold on;                   
                pekk = find(mod(phiv_kk(:), 1) < 1e-4 | mod(phiv_kk(:), 1) > 0.99) ;
                plot3(new3d(pekk, 1), new3d(pekk, 2), ...
                    new3d(pekk, 3), '^')                        
                % colormap summer
                % trimesh(prevf, new3d(:, 1), new3d(:, 2), new3d(:, 3),...
                %     mod(phiv_kk(:), 1), 'FaceColor', 'interp', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
                % freezeColors
                xlabel('x [\mum]')
                ylabel('y [\mum]')
                zlabel('z [\mum]')
                view(2)
                saveas(gcf, geomImFn)
            end

        elseif strcmp(phi_method, 'texture')
            imfn_sp_prev = sprintf( ...
                fullfile([imFolder_sp, '/', fileNameBase, '.tif']), ...
                xp.fileMeta.timePoints(tidx-1) ) ;

            % Load the intensity data            
            % Load 3D data for coloring mesh pullback
            xp.loadTime(t);
            xp.rescaleStackToUnitAspect();

            % Raw stack data
            IV = xp.stack.image.apply();
            IV = imadjustn(IV{1});         
            IVloaded = true ;

            % Texture patch options
            Options.PSize = 5;
            Options.EdgeColor = 'none';
            % Texture image options
            Options.imSize = ceil( 1000 .* [ 1 a_fixed ] );
            Options.yLim = [0 1];

            % fit the shifts in the y direction
            % todo: save uncorrected patchIms,
            % could try tiling twice...
            dmyk = 0 ;
            phi0_fit = zeros(size(uspace)) ;
            phi0s = zeros(size(uspace)) ;
            phi0_fit_kk = 1 ; % for first pass                
            phiv_kk = (vspace .* ones(nU, nV))' ;
            ensureDir([sphiDir, '/phi0_correction/'])
            while any(phi0_fit_kk > 0.002) && dmyk < 6
                disp(['Iteration ' num2str(dmyk)])
                plotfn = sprintf(phi0fitBase, t, dmyk);
                patchImFn = sprintf( ...
                    fullfile([sphiDir, '/phi0_correction/', fileNameBase, '_prephi0_' num2str(dmyk) '.tif']), ...
                    xp.fileMeta.timePoints(tidx-1) )  ;
                [phi0_fit_kk, phi0s_kk] = fitPhiOffsetsFromPrevPullback(IV, ...
                    new3d, cutMesh.umax, uspace * cutMesh.umax, phiv_kk, ...
                    ringpath_ss, imfn_sp_prev, lowerboundy, upperboundy, ...
                    save_ims, plotfn, Options, ...
                    step_phi0tile, width_phi0tile, potential_sigmay, 'integer', ...
                    patchImFn) ;

                % Update the result
                dmyk = dmyk + 1;
                phi0_fit = phi0_fit + phi0_fit_kk ;
                phi0s = phi0s + phi0s_kk ;
                phiv_kk = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;
            end
        else
            error("Could not recognize phi_method: must be 'texture' or '3dcurves'")
        end
        close all

        % Store to save at this timepoint
        phiv = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;
    end

    % NOTE: We have coordinates u,phiv that we associate with
    % the 3d coordinates already mapped to uv
    uphi = [uu(:), phiv(:)] ;

    % plot(uphi(:, 1), uphi(:, 2), '.')
    % xlabel('u')
    % ylabel('\phi')
    % waitfor(gcf)

    % Recompute radii_from_mean_uniform_rs as radii_from_avgpts 
    % NOTE: all radius calculations done in microns, not pixels
    sphi3d_rs = ((rot * new3d')' + trans) * resolution ;
    radii_from_avgpts = zeros(size(sphi3d_rs, 1), size(sphi3d_rs, 2)) ;
    for jj = 1:nU
        % Consider this hoop
        hoop = squeeze(sphi3d_rs(jj, :, :)) ;
        radii_from_avgpts(jj, :) = vecnorm(hoop - avgpts(jj, :), 2, 2) ;
    end

    % Triangulate the sphigrid and store as its own cutMesh
    % sphiv = zeros(nU, nV, 2) ;
    % sphiv(:, :, 1) = sv ;
    % sphiv(:, :, 2) = phiv ;
    sv = ringpath_ss .* onesUV ;
    % % Triangulate the mesh (topology is already known):
    % tmptri = defineFacesRectilinearGrid(sp, nU, nV) ;
    % % Old version did not assume topology as given:
    % tmptri = delaunay(sv(:), phiv(:)) ;
    % disp('orienting faces of delaunay triangulation (s,phi)')
    % tmptri = bfs_orient( tmptri );

    % Define path pairs for tiling the (s,phi) cut mesh
    spcutP1 = 1:nU;
    spcutP2 = nU*nV - fliplr(0:(nU-1)) ;
    spcutMesh.pathPairs = [ spcutP1', spcutP2' ];

    % Check to see if any members of pathPairs connect to
    % non-Nearest Neighbors. Not necessary now that we assume
    % known gridded mesh topology
    % cleantri = cleanBoundaryPath2D(tmptri, [sv(:), phiv(:)], spcutMesh.pathPairs(:), true) ;

    spcutMesh.f = defineFacesRectilinearGrid(uv, nU, nV) ;
    spcutMesh.nU = nU ;
    spcutMesh.nV = nV ;
    % First resampling
    spcutMesh.v0 = new3d ;
    % spcutMesh.vrs0 = ((rot * new3d')' + trans) * resolution ;
    % Define normals based on the original mesh normals
    spvn03d = interpolate2Dpts_3Dmesh(TF, TV2D, TVN3D, uphi) ;
    spvn03d = spvn03d ./ vecnorm(spvn03d, 2, 2) ;
    spcutMesh.vn0 = spvn03d ;
    spcutMesh.sphi0 = [sv(:), phiv(:)] ;
    spcutMesh.uphi0 = uphi ;
    % Note: uv has no direct relation with cutMesh, just a grid
    % for utility and reference, but it does have unequal 
    % spacing in u in anticipation of building sphi0 as a 
    % near perfect grid.
    spcutMesh.uv = uv ;  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SECOND RESAMPLING
    % Make a new grid
    slin = linspace(0, max(spcutMesh.sphi0(:, 1)), nU) ;
    plin = linspace(0, 1, nV) ;
    [ss, pp] = meshgrid(slin, plin) ;
    % Push the endpoints on each boundary in by epsilon to
    % avoid NaNs
    eps = 1e-14 ;
    ss(:, 1) = eps ;
    ss(:, end) = ss(:, end) - eps ;
    % Transpose so that x increases with increasing index first
    ss = ss' ;
    pp = pp' ;
    sp = [ss(:), pp(:)] ;

    % Tile the spcutMesh
    tileCount = [2, 2] ;
    spcutMesh.u = spcutMesh.sphi0 ;
    spcutMesh.v = spcutMesh.v0 ;
    spcutMesh.vn = spcutMesh.vn0 ;
    [ faces, v2d, v3d, vn3d ] = tileAnnularCutMesh( spcutMesh, tileCount );
    spcutMesh = rmfield(spcutMesh, 'u') ;
    spcutMesh = rmfield(spcutMesh, 'v') ;
    spcutMesh = rmfield(spcutMesh, 'vn') ;
    spv3d = interpolate2Dpts_3Dmesh(faces, v2d, v3d, sp) ;
    % check the pts
    % plot3(spv3d(:, 1), spv3d(:, 2), spv3d(:, 3))  

    % also interpolate the normals
    spvn3d = interpolate2Dpts_3Dmesh(faces, v2d, vn3d, sp) ;
    spvn3d = spvn3d ./ vecnorm(spvn3d, 2, 2) ;

    % Define new faces for second rectilinear resampling
    % NOTE: not necessary since we already defined the topology
    % from the guess [sv(:), phiv(:)] stored as spcutMesh.sphi0
    % spcutMesh.f = defineFacesRectilinearGrid(sp, nU, nV) ;
    spcutMesh.sphi = sp ;
    spcutMesh.v = spv3d ;
    spcutMesh.vn = spvn3d ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    spcutMesh.ringpath_ss = ringpath_ss ;
    spcutMesh.radii_from_mean_uniform_rs = radii_from_mean_uniform_rs ;  % from uniform DV sampling
    spcutMesh.radii_from_avgpts = radii_from_avgpts ;
    spcutMesh.mss = mss ;       % from uniform DV sampling, also stored in centerline
    spcutMesh.mcline = mcline ; % from uniform DV sampling, also stored in centerline
    spcutMesh.avgpts = avgpts ; % from uniform DV sampling, also stored in centerline
    spcutMesh.avgpts_ss = avgpts_ss ; % from uniform sampling, also stored in centerline

    % Define optimal isoareal Affine dilation factor in s
    % tmp = spcutMesh.sphi ;
    % tmp(:, 1) = tmp(:, 1) / max(tmp(:, 1)) ;
    % arsp = minimizeIsoarealAffineEnergy( spcutMesh.f, spcutMesh.v, tmp );
    % clearvars tmp
    spcutMesh.ar = cutMesh.ar ;

    % todo: check that u coords have not shifted upon
    % redefinition of sphi0 -> sphi

    % Save s,phi and their 3D embedding
    spcutMesh.phi0s = phi0s ;
    spcutMesh.phi0_fit = phi0_fit ;
    save(sprintf(spcutMeshBase, t), 'spcutMesh') ;
else
    disp('Loading spcutMesh from disk...')
    load(sprintf(spcutMeshBase, t), 'spcutMesh') ;
    IVloaded = false ;

    % Load new centerline
    fn = sprintf(clineDVhoopBase, t) ;
    disp(['Loading new centerline from ' fn])
    load(fn, 'mss', 'mcline', 'avgpts', 'avgpts_ss')
end
fprintf('Done with generating S,Phi coords \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Create pullback using S,Phi coords \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------
% Generate Output Image Files
%--------------------------------------------------------------
imfn = sprintf( fullfile([imFolder, '/', fileNameBase, '.tif']), t ); 
imfn_r = sprintf( fullfile([imFolder_r, '/', fileNameBase, '.tif']), t ) ;
imfn_sp = sprintf( fullfile([imFolder_sp, '/', fileNameBase, '.tif']), t ) ;
imfn_up = sprintf( fullfile([imFolder_up, '/', fileNameBase, '.tif']), t ) ;
pullbacks_exist1 = exist(imfn, 'file') && exist(imfn_r, 'file') ;
pullbacks_exist2 = exist(imfn_sp, 'file') && (exist(imfn_up, 'file') || ~generate_uphi_coord) ;
if (~pullbacks_exist1 || ~pullbacks_exist2 || overwrite_pullbacks) && ~IVloaded
    % Load 3D data for coloring mesh pullback
    xp.loadTime(t);
    xp.rescaleStackToUnitAspect();

    % Raw stack data
    IV = xp.stack.image.apply();
    IV = imadjustn(IV{1});
end

if ~exist(imfn_sp, 'file') || overwrite_pullbacks
    fprintf(['Generating SP output image: ' imfn_sp]);
    % Assigning field spcutMesh.u to be [s, phi] (ringpath
    % and azimuthal angle)
    spcutMesh.u = spcutMesh.sphi ;
    aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_sp)
    spcutMesh = rmfield(spcutMesh, 'u') ;
end

if (~exist(imfn_up, 'file') || overwrite_pullbacks) && generate_uphi_coord
    fprintf(['Generating uphi output image: ' imfn_up]);
    % Assigning field spcutMesh.u to be [s, phi] (ringpath
    % and azimuthal angle)
    spcutMesh.u = spcutMesh.uphi ;
    aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_up)
    spcutMesh = rmfield(spcutMesh, 'u') ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Output Image File -- regular UV coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(imfn, 'file') || overwrite_pullbacks
    % Generate output image in uv
    fprintf(['Generating output image: ' imfn]);
    aux_generate_orbifold(cutMesh, a_fixed, IV, imfn)
else
    disp('Skipping pullback image generation since exists')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save relaxed image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if ~exist(imfn_r, 'file')
    disp('Generating relaxed image for sphi coords...')
    spcutMesh.u = spcutMesh.sphi ;
    aux_generate_orbifold(spcutMesh, spcutMesh.ar, IV, imfn_r)
    spcutMesh = rmfield(spcutMesh, 'u') ;
end
