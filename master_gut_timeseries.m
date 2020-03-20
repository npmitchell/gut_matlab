%%
% initialize
%%
%% GUT_PIPELINE_20190415: GUT CONFORMAL MAPPING PIPELINE
% by Dillon Cislo & NPMitchell
%
% This is a pipeline to take the surface of the growing Drosophila gut and
% conformally map patches of it to the unit disk

%% INITIALIZE ImSAnE PROJECT ==============================================
%
% We start by clearing the memory and closing all figures
clear; close all; clc;
addpath_recurse('/mnt/crunch/djcislo/MATLAB/CGAL_Code/')

% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.
dataDir    =  cd; 
[ projectDir, ~, ~ ] = fileparts(matlab.desktop.editor.getActiveFilename); 
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
% Set required additional information on the files.
%
% We assume on individual image stack for each time point, labeled by time.
%  To be able to load the stack, we need to tell the project wehre the data
%  is, what convention is assumed for the file names, available time
%  points, and the stack resolution.  Options for modules in ImSAnE are
%  organized in MATLAB structures, i.e a pair of field names and values are
%  provided for each option.
%
% The following file metadata information is required:
%
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
fn = 'Time_%06d_c1_stab';

fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fn, '.tif'];
fileMeta.nChannels          = 1;
fileMeta.timePoints         = 11:169;
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
expMeta.channelsUsed        = 1;
expMeta.channelColor        = 1;
expMeta.description         = 'Apical membrane in Drosophila gut';
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
xp.loadTime(xp.fileMeta.timePoints(first_tp));
xp.rescaleStackToUnitAspect();

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Other options
%%%%%%%%%%%%%%%%%%%%%%%%%
% Warning: difference from CAAX
mlxprogram = 'surface_rm_resample20k_reconstruct_LS3_ssfactor4.mlx';
msls_axis_order = 'yxzc';
% Mesh marching options
normal_step = 10;
run_full_dataset = projectDir ; % 'none' ;

% Load/define the surface detection parameters
msls_detOpts_fn = fullfile(projectDir, 'msls_detectOpts.mat') ;
if exist(msls_detOpts_fn, 'file')
    load(detectOpts)
else
    channel = 1;
    foreGroundChannel = 1;
    ssfactor = 4;
    niter = 25 ;
    niter0 = 25 ;
    ofn_ply = 'mesh_apical_ms_stab_' ; 
    ofn_ls = 'msls_apical_stab_' ;
    ofn_smoothply = 'mesh_apical_stab_' ;
    ms_scriptDir = '/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/' ;
    lambda1 = 1 ;
    lambda2 = 1 ;
    exit_thres = 0.000001 ;
    smoothing = 0.10 ;
    nu = 0.00 ;
    pre_nu = -5 ;
    pre_smoothing = 1 ;
    post_nu = 2;
    post_smoothing = 4 ;
    zdim = 2 ;
    init_ls_fn = 'msls_apical_stab_init.h5';
    mlxprogram = './surface_rm_resample20k_reconstruct_LS3_1p2pc_ssfactor4.mlx' ;
    radius_guess = 40 ;
    center_guess = '100,100,100' ;
    dtype = 'h5' ;
    mask = 'none' ;
    prob_searchstr = '_stab_Probabilities.h5' ;
    imsaneaxisorder = 'xyzc'; ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc
    preilastikaxisorder= 'xyzc'; ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
    ilastikaxisorder= 'xyzc'; ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
    include_boundary_faces = true ;
    
    % save is
    save(msls_detOpts_fn, 'channel', ...
    'foreGroundChannel', ...
    'ssfactor', ...
    'niter', ...
    'niter0', ...
    'ofn_ply', ...
    'ofn_ls', ...
    'ofn_smoothply', ...
    'ms_scriptDir', ...
    'lambda1', ...
    'lambda2', ...
    'exit_thres', ...
    'smoothing', ...
    'nu', ...
    'pre_nu', ...
    'pre_smoothing', ...
    'post_nu', ...
    'post_smoothing', ...
    'zdim', ...
    'init_ls_fn', ...
    'mlxprogram', ...
    'radius_guess', ...
    'center_guess', ...
    'dtype', ...
    'mask', ...
    'prob_searchstr', ...
    'imsaneaxisorder', ...
    'preilastikaxisorder', ...
    'ilastikaxisorder', ...
    'include_boundary_faces') ;
end

% Name the output mesh directory ------------------------------------------
msls_exten = ['_prnu' strrep(strrep(num2str(pre_nu, '%d'), '.', 'p'), '-', 'n')];
msls_exten = [msls_exten '_prs' strrep(num2str(pre_smoothing, '%d'), '.', 'p') ];
msls_exten = [msls_exten '_nu' strrep(num2str(nu, '%0.2f'), '.', 'p') ];
msls_exten = [msls_exten '_s' strrep(num2str(smoothing, '%0.2f'), '.', 'p') ];
msls_exten = [msls_exten '_pn' num2str(post_nu, '%d') '_ps',...
    num2str(post_smoothing)];
msls_exten = [msls_exten '_l' num2str(lambda1) '_l' num2str(lambda2) ];
mslsDir = [projectDir 'msls_output'];
mslsDir = [mslsDir msls_exten '/'] ;

% The dimension to use to grab extremal seeds
seeddim = 3;

% Onion Options
nLayers = 3 ;  % nLayers must be an odd int
layerDistance = 5 ;  % layerDistance is in pix
sigma = 10 ;  % Sigma smooths
makeIP = 'MIP' ;  % SIP, MIP are options for makeIP
IPonly = false ;
onionOpts = struct('nLayers', nLayers, 'layerDistance', layerDistance,...
                   'sigma', sigma, 'makeIP', makeIP, 'IPonly', IPonly);
% SOI saving options
imwriteOptions = {'tif'};
soiDir = fullfile(projectDir, ['gut_apical_conformal_msls' msls_exten]);
soiDir = [soiDir '_' num2str(nLayers) 'layer/']; 
soi_save_options = struct('dir',soiDir,'imwriteOptions',{imwriteOptions},...
                    'make8bit',false);
                
%% DETECT THE SURFACE =====================================================
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
    'mslsDir', mslsDir, ...
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
    'center_guess', center_guess,... % xyz of the initial guess sphere ;
    'save', false, ... % whether to save images of debugging output
    'plot_mesh3d', false, ...
    'dtype', dtype,...
    'mask', mask,...
    'mesh_from_pointcloud', false, ...
    'prob_searchstr', prob_searchstr, ...
    'imsaneaxisorder', imsaneaxisorder, ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc
    'preilastikaxisorder', preilastikaxisorder, ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
    'ilastikaxisorder', ilastikaxisorder, ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
    'include_boundary_faces', include_boundary_faces) ;


% Set detect options ------------------------------------------------------
xp.setDetectOptions( detectOptions );

%% CREATE THE SUBSAMPLED H5 FILE FOR INPUT TO ILASTIK =====================
% skip if already done

for t = xp.fileMeta.timePoints
    
    if ~exist(fullfile(projectDir, [sprintf(fn, t) '.h5']), 'file')
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
disp('Open with ilastik if not already done')


%if ~exist([projectDir sprintf(fn, xp.currentTime) '.h5'], 'file')
 %   xp.detector.prepareIlastik(xp.stack);
  %  disp('done outputting downsampled data h5 for surface detection')
%else
 %   disp('h5 was already output, skipping...')
%end
%disp('Open with ilastik if not already done')

%% TRAIN DATA IN ILASTIK TO IDENTIFY APICAL/YOLK ==========================
% open ilastik, train until probabilities and uncertainty are satisfactory

%% Create MorphoSnakesLevelSet from the Probabilities from ilastik ========
fileMeta = xp.fileMeta ;
if strcmp(detectOptions.run_full_dataset, projectDir)
    disp(['Running dataset mode'])
    xp.detectSurface();
else
    % Morphosnakes for all remaining timepoints INDIVIDUALLY ==============
    for tp = fileMeta.timePoints(100:end)
        xp.setTime(tp);
        % xp.rescaleStackToUnitAspect();
        detectOptions.timepoint = xp.currentTime ;
        detectOptions.fileName = sprintf( fn, xp.currentTime );
        detectOptions.init_ls_fn = 'none' ;
        xp.setDetectOptions( detectOptions );
        xp.detectSurface();
    end
end

%% APDV ilastik training
% Train on anterior (A), posterior (P), and 
% dorsal anterior (D) points in different iLastik channels.
% anteriorChannel, posteriorChannel, and dorsalChannel specify the iLastik
% training channel that is used for each specification.
% Name the h5 file output from iLastik as ..._Probabilities_apcenterline.h5
% Train for anterior dorsal (D) only at the first time point, because
% that's the only one that's used.

%% 
% 3. align_meshes_APDV
%%
cd(mslsDir)
%
% OUTPUTS OF THIS SECTION
% -------
% xyzlim.txt 
%   xyzlimits of raw meshes in units of full resolution pixels (ie not
%   downsampled)
% xyzlim_APDV.txt 
%   xyzlimits of rotated and translated meshes in units of full resolution 
%   pixels (ie not downsampled)
% xyzlim_APDV_um.txt 
%   xyz limits of rotated and translated meshes in microns
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
%   Centers of mass for A, P, and D in microns in rotated APDV coord system
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
overwrite = true ;  % recompute APDV rotation, translation
overwrite_apdvcoms = false ;  % recompute APDV coms from training
save_figs = true ;  % save images of cntrline, etc, along the way
overwrite_ims = true ;  % overwrite images even if centerlines are not overwritten
preview = false ;  % display intermediate results, for debugging
resolution = xp.fileMeta.stackResolution(1) ;  % um per pixel for full resolution (not subsampled)
dorsal_thres = 0.5 ;  % threshold for extracting Dorsal probability cloud 
buffer = 5 ;  % extra space in meshgrid of centerline extraction, to ensure mesh contained in volume
plot_buffer = 40;  % extra space in plots, in um
% ssfator == subsampling factor for the h5s used to train for mesh/acom/pcom/dcom
weight = 0.1;  % for speedup of centerline extraction. Larger is less precise
normal_step = 0.5 ;  % how far to move normally from ptmatched vtx if a/pcom is not inside mesh
eps = 0.01 ;  % value for DT outside of mesh in centerline extraction
meshorder = 'zyx' ;  % ordering of axes in loaded mesh wrt iLastik output
anteriorChannel = 1;  % which channel of APD training is anterior
posteriorChannel = 2;  % which channel of APD training is posterior 
dorsalChannel = 4 ;  % which channel of APD training is dorsal
axorder = [2, 1, 3] ;  % axis order for APD training output

opts.overwrite = overwrite ;
opts.overwrite_ims = overwrite_ims ;
opts.resolution = resolution ;
opts.dorsal_thres = dorsal_thres ;
opts.buffer = buffer ;
opts.plot_buffer = plot_buffer ;
opts.ssfactor = ssfactor ;
opts.weight = weight ;
opts.normal_step = normal_step ;
opts.eps = eps ;
opts.meshorder = meshorder ;
opts.anteriorChannel = anteriorChannel ;
opts.posteriorChannel = posteriorChannel ;
opts.dorsalChannel = dorsalChannel ;
opts.axorder = axorder ;
opts.preview = preview ;
opts.meshdir = mslsDir ;

alignMeshesAPDV(opts)