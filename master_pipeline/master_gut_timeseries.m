%% GUT MASTER PIPELINE FOR TIME SERIES DATA
% by NPMitchell & Dillon Cislo
%
% This is a pipeline to take the surface of the growing Drosophila gut and
% conformally map patches of it to the unit disk

% CONVENTIONS
% -----------
% xp.stack.image.apply() gives imageLength, imageWidth, imageDepth as XYZ.
% This forms our definition for XYZ. However, this is mirrored wrt the
% physical coordinates in the lab frame. 
% Note that this is different from bfGetReader(fn), which gives yxz -- ie
%   imageWidth, imageLength, imageDepth. In Experiment.m, y and x are 
% swapped to account for this discrepancy after using bfGetReader to build
% the empty array.
%
% Python will read h5s with dimensions backward wrt MATLAB, so if h5s are 
% saved as cxyz in MATLAB, morphsnakes/python will read as zyxc. 
% preilastikaxisorder converts from xyzc to what is given for h5. 
% Typically, ilasik puts the channel at the start in MATLAB (end in python)
% Therefore, set axisorder = 'cxyz' since MATLAB probability h5s are saved 
% as cxyz.
%
% Meshes from integralDetector are stored with normals 'outward'/basal.
% It is important to preserve the triangle orientation for texturepatch,
% but the mesh XY coordinates are flipped with respect to the physical lab
% coordinates. Therefore, we set flipy = true, and the aligned meshes have 
% apical normals. The APDV meshes are black in meshlab when viewed from the
% outside.  flipy indicates that the APDV coordinate system is mirrored 
% across XZ wrt data coordinate system XYZ. 

%% Global options
% We start by clearing the memory and closing all figures
clear; close all; clc;
% change this path, for convenience
% cd /mnt/crunch/48Ygal4-UAShistRFP/201904031830_great/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2/data/deconvolved_16bit
dataDir = cd ;

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
overwrite_cutMesh = true ;

%% DEFINE NEW MASTER SETTINGS
if overwrite_masterSettings || ~exist('./masterSettings.mat', 'file')
    % Metadata about the experiment
    stackResolution = [.2619 .2619 .2619] ;
    nChannels = 1 ;
    timePoints = 0:169 ;
    ssfactor = 4 ;
    % whether the data is stored inverted relative to real position
    flipy = true ; 
    timeinterval = 1 ;  % physical interval between timepoints
    timeunits = 'min' ; % physical unit of time between timepoints
    scale = 0.02 ;      % scale for conversion to 16 bit
    file32Base = 'TP%d_Ch0_Ill0_Ang0,45,90,135,180,225,270,315.tif'; 
    % file32Base = 'TP%d_Ch0_Ill0_Ang0,60,120,180,240,300.tif'; 
    fn = 'Time_%06d_c1_stab';
    fn_prestab = 'Time_%06d_c1.tif';
    set_preilastikaxisorder = 'xyzc' ;
    masterSettings = struct('stackResolution', stackResolution, ...
        'nChannels', nChannels, ...
        'timePoints', timePoints, ...
        'ssfactor', ssfactor, ...
        'flipy', flipy, ...
        'timeinterval', timeinterval, ...
        'timeunits', timeunits, ...
        'scale', scale, ...
        'file32Base', file32Base, ...
        'fn', fn,...
        'fn_prestab', fn_prestab, ...
        'set_preilastikaxisorder', set_preilastikaxisorder); 
    disp('Saving masterSettings to ./masterSettings.mat')
    if exist('./masterSettings.mat', 'file')
        ui = userinput('This will overwrite the masterSettings. Proceed (Y/n)?');
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
    timePoints = masterSettings.timePoints ;
    ssfactor = masterSettings.ssfactor ;
    % whether the data is stored inverted relative to real position
    flipy = masterSettings.flipy ; 
    timeinterval = masterSettings.timeinterval ;  % physical interval between timepoints
    timeunits = masterSettings.timeunits ; % physical unit of time between timepoints
    scale = masterSettings.scale ;      % scale for conversion to 16 bit
    file32Base = masterSettings.file32Base ; 
    fn = masterSettings.fn ;
    fn_prestab = masterSettings.fn_prestab ;
    set_preilastikaxisorder = masterSettings.set_preilastikaxisorder ;
end
dir32bit = fullfile(dataDir, 'deconvolved_32bit') ;
dir16bit = fullfile(dataDir, 'deconvolved_16bit') ;
dir16bit_prestab = fullfile(dir16bit, 'data_pre_stabilization') ;

%% PATHS ==================================================================
disp('Adding paths')
ms_scriptDir = '/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/' ;
codepath = '/mnt/data/code/' ;
gutpath = fullfile(codepath, 'gut_matlab') ;
meshlabCodeDir = fullfile(codepath, 'meshlab_codes') ;
addpath_recurse('/mnt/crunch/djcislo/MATLAB/CGAL_Code/')
addpath_recurse(fullfile(codepath, 'gptoolbox'))
addpath(fullfile(gutpath, 'master_pipeline'))
addpath(fullfile(gutpath, 'data_handling'))
addpath_recurse(fullfile(gutpath, 'mesh_handling'))
addpath(fullfile(gutpath, 'basics'))
addpath(fullfile(gutpath, 'distanceTransform'))
addpath_recurse(fullfile(gutpath, 'plotting'))
addpath(fullfile(gutpath, 'savgol'))
addpath_recurse(fullfile(gutpath, 'toolbox_fast_marching/'));
addpath(genpath(fullfile(gutpath, 'euclidean_orbifolds')));
addpath(genpath(fullfile(gutpath, 'TexturePatch')));
addpath_recurse(fullfile(gutpath, ['axisymmetric_pullbacks' filesep])) ;
addpath_recurse(fullfile(gutpath, ['plotting' filesep])) ;
addpath_recurse(fullfile(gutpath, ['h5_handling' filesep])) ;
addpath_recurse(fullfile(gutpath, ['curve_functions' filesep])) ;
addpath(fullfile(gutpath, ['ExtPhaseCorrelation' filesep])) ;
addpath(fullfile(gutpath, 'savgol')) ;
addpath(fullfile(codepath, 'DEC')) ;
addpath(fullfile(codepath, 'TexturePatch_for_git', 'TexturePatch')) ;
% addpath(genpath('/mnt/crunch/djcislo/MATLAB/TexturePatch'));
disp('done')

%% Define some colors
[colors, color_names] = define_colors() ;
blue = [0 0.4470 0.7410] ;
orange = [0.8500 0.3250 0.0980] ;
yellow = [0.9290, 0.6940, 0.1250] ;
purple = [0.4940, 0.1840, 0.5560] ;
green = [0.4660, 0.6740, 0.1880] ;
sky = [0.3010, 0.7450, 0.9330] ;
red = [0.6350, 0.0780, 0.1840] ;
brick = [0.800000 0.250000 0.330000] ;
light_green =[0.560000 0.930000 0.560000] ;
light_gray = [0.830000 0.830000 0.830000] ;
bwr = diverging_cmap([0:0.01:1], 1, 2) ;
disp('done adding colors')

%% END OF EXPERIMENT METADATA =============================================
% =========================================================================
% =========================================================================
%% -IIV. make MIPs for 32bit images
mipDir = fullfile(dir32bit, 'mips32bit') ;
Options.overwrite_mips = overwrite_mips ;
Options.scale = scale ;
makeMips(timePoints, dir32bit, file32Base, mipDir, Options)

%%  -IV. convert 32 to 16bit images
convert32to16bit(timePoints, scale, dir32bit, dir16bit_prestab,...
    file32Base, fn_prestab)

%% Rename stab to prestab
% fns = fullfile('./deconvolved_16bit/Time*stab')
% for qq=1:length(fns)
%     command = ['mv ' fullfile(fns.folder, fns.name) fullfile(dir16bit, fns.name
% end

%% -III. make MIPs for 16bit images
mipDir = fullfile(dir16bit_prestab, 'mips') ;
Options.overwrite_mips = overwrite_mips ;
makeMips(timePoints, dir16bit_prestab, fn_prestab, mipDir, Options)

%%  -II. stabilizeImagesCorrect.m
% name of directory to check the stabilization of mips
mips_stab_check = fullfile(mipDir, 'stab_check') ;
mipoutdir = fullfile(mipDir, 'mips_stab') ;
im_intensity = 1 ; % 0.01 ;
imref_intensity = 1 ; % 0.005 ;
% Choose bit depth as typename
typename = 'uint16' ;
% Give file names for I/O
fileNameIn = fullfile(dir16bit_prestab, fn_prestab) ;
fileNameOut = fullfile(dir16bit, [fn '.tif']) ;
rgbName = [fn '.png'] ;
typename = 'uint16' ; 
overwrite_mips = false ;
overwrite_tiffs = false ;
stabilizeImages(fileNameIn, fileNameOut, rgbName, typename, ...
    timePoints, timePoints, timePoints(1), ...
    im_intensity, imref_intensity, ...
    mipDir, mipoutdir, mips_stab_check, overwrite_mips, overwrite_tiffs)

%%   -I. master_gut_timeseries_prestab_for_training.m
cd(dir16bit)
dataDir = cd ;
makeH5SeriesPrestabForTraining(fn_prestab(1:end-4), nChannels, timePoints, ...
    ssfactor, stackResolution, dir16bit_prestab, set_preilastikaxisorder)
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
xp.setTime(xp.fileMeta.timePoints(1)) ;
% xp.loadTime(xp.fileMeta.timePoints(first_tp));
% xp.rescaleStackToUnitAspect();

%% SET DETECT OPTIONS =====================================================
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
    ssfactor = 4;
    niter = 15 ;
    niter0 = 15 ;
    ofn_ply = 'mesh_apical_ms_stab_' ; 
    ofn_ls = 'msls_apical_stab_' ;
    ofn_smoothply = 'mesh_apical_stab_' ;
    lambda1 = 1 ;
    lambda2 = 1 ;
    exit_thres = 0.0001 ;
    smoothing = 0.10 ;
    nu = 0.00 ;
    pre_nu = -4 ;
    pre_smoothing = 1 ;
    post_nu = 2;
    post_smoothing = 4 ;
    zdim = 2 ;
    init_ls_fn = 'msls_initguess.h5';
    % mlxprogram = fullfile(meshlabCodeDir, ...
    %     'laplace_refine_HCLaplace_LaplaceSPreserve_QuadEdgeCollapse60kfaces.mlx') ;
    mlxprogram = fullfile(meshlabCodeDir, ...
        'surface_rm_resample20k_reconstruct_LS3_1p2pc_ssfactor4.mlx') ;
    radius_guess = 40 ;
    center_guess = '100,100,100' ;
    dtype = 'h5' ;
    mask = 'none' ;
    prob_searchstr = '_stab_Probabilities.h5' ;
    preilastikaxisorder= set_preilastikaxisorder; ... % axis order in input to ilastik as h5s. To keep as saved coords use xyzc
    ilastikaxisorder= 'cxyz'; ... % axis order as output by ilastik probabilities h5
    imsaneaxisorder = 'xyzc'; ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc
    include_boundary_faces = true ;
    
    % Name the output mesh directory ------------------------------------------
    % msls_exten = ['_prnu' strrep(strrep(num2str(pre_nu, '%d'), '.', 'p'), '-', 'n')];
    % msls_exten = [msls_exten '_prs' strrep(num2str(pre_smoothing, '%d'), '.', 'p') ];
    % msls_exten = [msls_exten '_nu' strrep(num2str(nu, '%0.2f'), '.', 'p') ];
    % msls_exten = [msls_exten '_s' strrep(num2str(smoothing, '%0.2f'), '.', 'p') ];
    % msls_exten = [msls_exten '_pn' num2str(post_nu, '%d') '_ps',...
    %     num2str(post_smoothing)];
    % msls_exten = [msls_exten '_l' num2str(lambda1) '_l' num2str(lambda2) ];
    meshDir = [projectDir 'msls_output'];
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
        'ofn_smoothply', ofn_smoothply, ...
        'ms_scriptDir', ms_scriptDir, ...
        'timepoint', xp.currentTime, ...
        'zdim', zdim, ...
        'pre_nu', pre_nu, ...
        'pre_smoothing', pre_smoothing, ...
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
        'preilastikaxisorder', preilastikaxisorder, ... 
        'ilastikaxisorder', ilastikaxisorder, ... 
        'physicalaxisorder', imsaneaxisorder, ... 
        'include_boundary_faces', include_boundary_faces, ...
        'smooth_with_matlab', true) ;

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

%% CREATE THE SUBSAMPLED H5 FILE FOR INPUT TO ILASTIK =====================
% skip if already done

for tt = xp.fileMeta.timePoints
    if ~exist(fullfile(projectDir, 'stabilized_h5s', [sprintf(fn, tt) '.h5']), 'file')
        disp(['Did not find file: ', fullfile(projectDir, 'stabilized_h5s', [sprintf(fn, tt) '.h5'])])
        xp.loadTime(tt);
        xp.rescaleStackToUnitAspect();
        % make a copy of the detectOptions and change the fileName
        detectOpts2 = detectOptions ;
        detectOpts2.fileName = sprintf( fn, xp.currentTime ) ;
        xp.setDetectOptions( detectOpts2 );
        xp.detector.prepareIlastik(xp.stack);
        disp(['done outputting downsampled data h5: tp=' num2str(tt) ' for surface detection'])
    else
        disp(['h5 ' num2str(tt) ' was already output, skipping...'])
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

%% TRAIN NON-STABILIZED DATA IN ILASTIK TO IDENTIFY APICAL/YOLK ===========
% open ilastik, train pre-stab h5s until probabilities and uncertainty are 
% satisfactory, then run on stab images.

%% Note that old h5's used to be different order. To convert, do
% if false
%     tmp = h5read(fullfile(meshDir, init_ls_fn), '/implicit_levelset');
%     % OLD ORDER: yxz for implicit levelset
%     tmp2 = permute(tmp, [3,2,1]);
%     h5create(fullfile(meshDir, init_ls_fn), '/implicit_levelset', size(tmp2), 'Datatype', 'int8')
%     h5write(fullfile(meshDir, init_ls_fn), '/implicit_levelset', tmp2)
% end

%% Create MorphoSnakesLevelSet from the Probabilities from ilastik ========
% Now detect all surfaces
if strcmp(detectOptions.run_full_dataset, projectDir)
    assert(run_full_dataset_ms)
    disp('Running dataset mode')
    xp.detectSurface();
else
    assert(~run_full_dataset_ms)
    assert(strcmp(detectOptions.run_full_dataset, 'none'))
    % Morphosnakes for all remaining timepoints INDIVIDUALLY ==============
    for tp = xp.fileMeta.timePoints
        try
            xp.setTime(tp);
            % xp.loadTime(tp) ;
            % xp.rescaleStackToUnitAspect();

            % make a copy of the detectOptions and change the fileName
            detectOpts2 = detectOptions ;
            detectOpts2.timepoint = xp.currentTime ;
            detectOpts2.fileName = sprintf( fn, xp.currentTime );
            xp.setDetectOptions( detectOpts2 );
            xp.detectSurface();
            % For next time, use the output mesh as an initial mesh
            detectOpts2.init_ls_fn = 'none' ;
        catch
            disp('Could not create mesh -- skipping for now')
            % On next timepoint, use the tp previous to current time
            detectOptions.init_ls_fn = [detectOptions.ofn_ls, ...
                    num2str(tp - 1, '%06d' ) '.' detectOptions.dtype]
        end
    end
end

%% Inspect a single mesh
tp = 50 ;
meshfn = fullfile(meshDir, sprintf([ meshFileBase '.ply'], tp)) ;    
xp.loadTime(tp)
xp.rescaleStackToUnitAspect() ;
IV = xp.stack.image.apply() ;
IV = IV{1} ;
mesh = read_ply_mod(meshfn) ;
for leaf=1:40:size(IV, 1)
    inds = find(abs(mesh.v(:, 1) - leaf) < 5) ;
    imshow(imadjust(squeeze(IV(leaf, :, :))'))
    if any(inds)
        hold on;
        plot(mesh.v(inds, 2), mesh.v(inds, 3), 'co')
    end
    pause(1)
    clf
end
close all

%% Inspect all meshes in 3D

% Make an output directory for the quick-and-dirty inspection
outputdir = fullfile(meshDir, 'quick_mesh_inspect') ;
for tp = xp.fileMeta.timePoints
    % Load the mesh
    meshfn = fullfile(meshDir, sprintf([ meshFileBase '.ply'], tp)) ;    
    mesh = read_ply_mod(meshfn) ;
    % Plot the mesh in 3d. Color here by Y coordinate
    trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
        mesh.v(:, 2), 'edgecolor', 'none', 'Facealpha', 0.5)
    saveas(gcf, fullfile(outputdir, sprintf('inspect_%04d.png', tp)))
    close all
end

%% adjust all meshes by 0.5 --- now this is done in integralDetector
% for tp = fileMeta.timePoints
%     disp(num2str(tp))
%     meshfn = fullfile(meshDir, sprintf([ofn_smoothply '%06d.ply'], tp)) ;
%     outfn = fullfile(projectDir, sprintf([ofn_smoothply '%06d.ply'], tp)) ;
%     if ~exist(outfn, 'file')
%         mesh = read_ply_mod(meshfn) ;
%         mesh.v = mesh.v + 0.5 * [1 1 1] ;
%         plywrite_with_normals(outfn, mesh.f, mesh.v, mesh.vn)
%     else
%         disp('skip')
%     end
%     
%     meshfn = fullfile(meshDir, sprintf([ofn_ply '%06d.ply'], tp)) ;
%     outfn = fullfile(projectDir, sprintf([ofn_ply '%06d.ply'], tp)) ;
%     if ~exist(outfn, 'file')
%         mesh = read_ply_mod(meshfn) ;
%         mesh.v = mesh.v + 0.5 * [1 1 1] ;
%         plywrite(outfn, mesh.f, mesh.v)
%     else
%         disp('skip')
%     end
% end
% disp('done')


%% Check that all have been created
for tp = xp.fileMeta.timePoints
    meshfn = fullfile(meshDir, sprintf([ofn_smoothply '%06d.ply'], tp)) ;
       
    % Check that smoothply exists
    if ~exist(meshfn, 'file') 
        xp.setTime(tp);
        % xp.loadTime(tp) ;
        % xp.rescaleStackToUnitAspect();

        % make a copy of the detectOptions and change the fileName
        detectOpts2 = detectOptions ;
        detectOpts2.timepoint = xp.currentTime ;
        detectOpts2.fileName = sprintf( fn, xp.currentTime );
        if tp > xp.fileMeta.timePoints(1)
            detectOpts2.init_ls_fn = 'none' ;
        end
        xp.setDetectOptions( detectOpts2 );
        xp.detectSurface();
    else
        disp('exists')
    end
end
disp('done')

clearvars lambda1 lambda2 ilastikaxisorder mlxprogram ms_scriptDir
clearvars msls_exten msls_detOpts_fn niter niter0 nu 

%% Define QuapSlap object
opts.meshDir = meshDir ;
opts.flipy = flipy ;
opts.timeinterval = timeinterval ;
opts.timeunits = timeunits ;
opts.nV = 100 ;
opts.nU = 100 ;
opts.normalShift = 10 ;
opts.a_fixed = 2.0 ;
opts.adjustlow = 1.00 ;                  %  floor for intensity adjustment
opts.adjusthigh = 99.9 ;                 % ceil for intensity adjustment (clip)
opts.phiMethod = 'curves3d' ;
disp('defining QS')
QS = QuapSlap(xp, opts) ;
disp('done')

%% APDV ilastik training
% Train on anterior (A), posterior (P), and 
% dorsal anterior (D) points in different iLastik channels.
% Perform on pre-stabilized H5s, NOT on stabilized H5s. 
% anteriorChannel, posteriorChannel, and dorsalChannel specify the iLastik
% training channel that is used for each specification.
% Name the h5 file output from iLastik as ..._Probabilities_apcenterline.h5
% Train for anterior dorsal (D) only at the first time point, because
% that's the only one that's used.

%% 3. align_meshes_APDV or load transformations if already done
% startendptH5FileName = fullfile(centerlineDir, 'startendpt.h5') ;

% Try to load the results to test if we must do calculation 
try
    [rot, trans] = QS.getRotTrans() ;
    [xyzlim_raw, xyzlim, xyzlim_um] = getXYZLims(QS) ;
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
    axorder = [1, 2, 3] ;  % axis order for APD training output

    clearvars opts
    optsfn = fullfile(projectDir, 'alignAPDV_Opts.mat') ;
    if exist(optsfn, 'file') && ~overwrite_alignAPDVOpts
        disp('Loading options from disk')
        QS.load(optsfn, 'alignAPDVOpts')
    else
        disp('No alignAPDV_Opts on disk or overwriting, defining')
        apdvOpts.smwindow = 30 ;
        apdvOpts.dorsal_thres = dorsal_thres ;
        apdvOpts.buffer = buffer ;    error('here')
        apdvOpts.plot_buffer = plot_buffer ;
        alignAPDVOpts.weight = weight ;
        alignAPDVOpts.normal_step = normal_step ;
        alignAPDVOpts.eps = eps ;
        alignAPDVOpts.meshorder = meshorder ;
        alignAPDVOpts.anteriorChannel = anteriorChannel ;
        alignAPDVOpts.posteriorChannel = posteriorChannel ;
        alignAPDVOpts.dorsalChannel = dorsalChannel ;
        alignAPDVOpts.axorder = axorder ;
        % filename pattern for the apdv training probabilities
        alignAPDVOpts.apdProbFileName = fullfile(projectDir, 'stabilized_h5s', ...
            [fn '_Probabilities_apcenterline.h5']) ;
        % alignAPDVOpts.fn = fn ;  % filename base

        % save the optionsc
        save(optsfn, 'alignAPDVOpts')
    end

    % Add script-instance-specific options
    alignAPDVOpts.overwrite = overwrite_APDVCOMs ;  % recompute APDV coms from training
    alignAPDVOpts.check_slices = false ;
    alignAPDVOpts.preview = preview ;
    alignAPDVOpts.preview_com = false ;
    alignAPDVOpts.apdvoutdir = fullfile(meshDir, 'centerline') ;
    
    % Compute the APD COMs
    [acom_sm, pcom_sm] = computeAPDCOMs(QS, alignAPDVOpts) ;
    
    % Align the meshes APDV & plot them
    opts.overwrite_ims = overwrite_alignedMeshIms ;  % overwrite images even if centerlines are not overwritten
    opts.overwrite = overwrite_APDVCOMs || overwrite_APDVMeshAlignment ; % recompute APDV rotation, translation
    [rot, trans, xyzlim_raw, xyzlim, xyzlim_um] = QS.alignMeshesAPDV(acom_sm, pcom_sm, opts) ;
else
    disp('Already done')
end
disp('done')
clearvars normal_step 

%% MAKE MASKED DATA FOR PRETTY VIDEO ==================================
QS.generateMaskedData()

%% MAKE ORIENTED MASKED DATA FOR PRETTY VIDEO ==================================
QS.alignMaskedDataAPDV()

%% PLOT ALL TEXTURED MESHES IN 3D =========================================
% Get limits and create output dir
% Name output directory
figoutdir = fullfile(meshDir, 'images_texturepatch') ;
figddir = fullfile(figoutdir, 'dorsal') ;
figvdir = fullfile(figoutdir, 'ventral') ;
figlat1dir = fullfile(figoutdir, 'lateral1') ;
figlat2dir = fullfile(figoutdir, 'lateral2') ;
dirs = {figoutdir, figddir, figvdir, figlat1dir, figlat2dir} ;
for i = 1:length(dirs)
    if ~exist(dirs{i}, 'dir')
        mkdir(dirs{i}) ;
    end
end
% Establish texture patch options
plot_buffer = 20 ;
xyzbuff = xyzlim_um ;
xyzbuff(1, 1) = xyzbuff(1, 1) - plot_buffer ;
xyzbuff(1, 2) = xyzbuff(1, 2) + plot_buffer ;
% Make other axis limits same range as x range
aratio = 0.61803398875 ;  % 1 / golden ratio
yzrange = aratio * (xyzbuff(1, 2) - xyzbuff(1, 1)) ;
ymid = 0.0 ;
zmid = 0.5 * (xyzbuff(3, 1) + xyzbuff(3, 2)) ;
xyzbuff(2, :) = ymid + [-0.5 * yzrange, 0.5 * yzrange] ;
xyzbuff(3, :) = zmid + [-0.5 * yzrange, 0.5 * yzrange] ;

meshFileName = fullfile(meshDir, [ meshFileBase '.ply']) ;
figdirs = {figddir, figvdir, figlat1dir, figlat2dir} ;
metafn = fullfile(figoutdir, 'metadat.mat') ;

if ~exist(metafn, 'file') || overwrite_TextureMeshOpts
    % Define & Save metadata
    metadat.normal_shift = 10 ;                 % pixels
    metadat.xyzlim = xyzbuff ;                  % xyzlimits
    metadat.flipy = flipy ;                     % invert lateral dimension
    metadat.texture_axis_order = [1, 2, 3] ;    % texture space sampling
    metadat.reorient_faces = false ;            % if some normals are inverted
    metadat.timeinterval = timeinterval ;       % physical increment btween timepts
    metadat.timeunits = timeunits ;             % physical unit of time increment
    metadat.t0 = xp.fileMeta.timePoints(1) ;    % time offset (first timepoint --> t=0)
    metadat.adjustlow = 1.00 ;                  %  floor for intensity adjustment
    metadat.adjusthigh = 99.5 ;                 % ceil for intensity adjustment (clip)
    % Psize is the linear dimension of the grid drawn on each triangular face
    Options.PSize = 5;
    Options.EdgeColor = 'none';
    Options.Rotation = rot ;
    Options.Translation = trans ;
    Options.Dilation = resolution ;
    Options.numLayers = [1, -1];  % at layerSpacing 2, 2 marches ~0.5 um 
    Options.layerSpacing = 2 ;
    % Save it
    save(metafn, 'metadat', 'Options')
else
    load(metafn, 'metadat', 'Options')
end

% Plot on surface for all TP
plotSeriesOnSurfaceTexturePatch(xp, meshFileName, figdirs, overwrite, ...
    metadat, Options)
clearvars Options xyzbuff 

%% EXTRACT CENTERLINES
% Note: these just need to be 'reasonable' centerlines for topological
% checks on the orbifold cuts.
exponent = 1.0 ;
res = 3.0 ; 
cntrlineOpts.overwrite = overwrite_centerlines ;     % overwrite previous results
cntrlineOpts.overwrite_ims = overwrite_centerlineIms ;     % overwrite previous results
cntrlineOpts.weight = 0.1;            % for speedup of centerline extraction. Larger is less precise
cntrlineOpts.exponent = exponent ;            % how heavily to scale distance transform for speed through voxel
cntrlineOpts.res = res ;               % resolution of distance tranform grid in which to compute centerlines
cntrlineOpts.preview = false ;           % preview intermediate results
cntrlineOpts.reorient_faces = false ;    % not needed for our well-constructed meshes
cntrlineOpts.dilation = 1 ;              % how many voxels to dilate the segmentation inside/outside before path computation

% Note: this takes about 400s per timepoint for res=2.0
%
QS.extractCenterlineSeries(cntrlineOpts)
disp('done with centerlines')

%% Fix flip in Y for centerlines
% aux_fix_flip

%% Surface area and volume calc

%% Cylinder cut mesh
if overwrite_endcapOpts
    endcapOpts = ...
        struct( 'adist_thres', 20, ...  % distance threshold for cutting off anterior in pix
                'pdist_thres', 25);     % distance threshold for cutting off posterior in pix
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
methodOpts.preview = false ;     % display intermediate results
sliceMeshEndcaps(QS, endcapOpts, methodOpts) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORBIFOLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
% Overwriting options
overwrite_pullbacks = true ;
overwrite_cutMesh = false ;
overwrite_spcutMesh = true ;
overwrite_writhe = false ;
overwrite_SmRSIms = false ;
overwrite_spcutMeshSm = false ;
overwrite_folds = false ;
overwrite_lobedynamics = false ;
overwrite_foldims = false ;
overwrite_lobeims = false ;
overwrite_spcutMesh_smoothradii = false ;
overwrite_piv = false ;
% Other options for what to do
save_ims = true ;
nCurves_yjitter = 100 ;
nCurves_sphicoord = 1000 ;
% Plotting params
preview = false ;
washout2d = 0.5 ;
washout3d = 0.5 ;
% Parameters for cutMesh creation
maxJitter = 100 ;
maxTwChange = 0.15 ;
nsegs4path = 5 ;
a_fixed = 2 ;
% Parameters for spcutMesh creation
normal_shift = 10 ;
maxJitter = 100 ;
maxTwChange = 0.15 ;
% for phi0 calculation via texture matching
lowerboundy = -350 ;
upperboundy = 350 ;
step_phi0tile = 25 ;
width_phi0tile = 150 ;
potential_sigmay = 350 ;


%% Identify anomalies in centerline data
idOptions.ssr_thres = 15 ;  % distance of sum squared residuals in um as threshold
idOptions.overwrite = overwrite_idAnomClines ;
QS.generateCleanCntrlines(idOptions) ;

%% Clean Cylinder Meshes
cleanCylOptions.overwrite = overwrite_cleanCylMesh ;
cleanCylOptions.save_ims = true ;
QS.cleanCylMeshes(cleanCylOptions)
    
%% Iterate Through Time Points to Create Pullbacks ========================
% outcutfn = fullfile(cutFolder, 'cutPaths_%06d.txt') ;
for tt = xp.fileMeta.timePoints(142:end)
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
        error('here')
        if exist(cutMeshfn, 'file')
            disp('Overwriting cutMesh...') ;
        else
            disp('cutMesh not found on disk. Generating cutMesh... ');
        end
        QS.generateCurrentCutMesh()
        disp('Saving cutP image')
        % Plot the cutPath (cutP) in 3D
        QS.plotCutPath(QS.currentMesh.cutMesh, QS.currentMesh.cutPath)
        compute_pullback = true ;
    else
        fprintf('Loading Cut Mesh from disk... ');
        QS.loadCurrentCutMesh()
        compute_pullback = ~isempty(QS.currentMesh.cutPath) ;
    end

    spcutMeshOptions.overwrite = overwrite_spcutMesh ;
    spcutMeshOptions.save_phi0patch = true ;
    spcutMeshOptions.iterative_phi0 = true ;
    spcutMeshOptions.smoothingMethod = 'none' ;
    QS.plotting.preview = true ;
    QS.generateCurrentSPCutMesh([], spcutMeshOptions) ;
    
    % Compute the pullback if the cutMesh is ok
    if compute_pullback
        pbOptions.overwrite = overwrite_pullbacks ;
        pbOptions.generate_uv = false ;
        pbOptions.generate_uphi = false ;
        pbOptions.generate_relaxed = false ;
        QS.generateCurrentPullbacks([], [], pbOptions) ;
    else
        disp('Skipping computation of pullback')
    end
    clear Options IV
    error('here')

    % %% Save SMArr2D (vertex positions in the 2D pullback) -----------------
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

%% Preview results ========================================================
check = false ;
if check
    aux_preview_results
end

%% Phase correlation to get x shift


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TILE IMAGES IN Y AND RESAVE ============================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imDirs = {imFolder, imFolder_sp, imFolder_up} ;
imDirs_e = {imFolder_e, imFolder_sp_e, imFolder_up_e} ;
ntiles = 50 ;   %   The number of bins in each dimension for histogram equilization for a
                %   square original image. That is, the extended image will have (a_fixed *
                %   ntiles, 2 * ntiles) bins in (x,y).
for qq = 1:2
    options.histeq = true ;
    options.a_fixed = a_fixed ;
    options.ntiles = ntiles ;
    options.overwrite = overwrite_pullbacks;
    extendImages(imDirs{qq}, imDirs_e{qq}, fileNameBase, options)
    disp(['done ensuring extended tiffs for ' imDirs{qq} ' in ' imDirs_e{qq}])
end
disp('done')

%% Estimate cell density

%% FIND THE FOLDS SEPARATING COMPARTMENTS =================================
% First compute using the avgpts (DVhoop means)
guess123 = [0.2, 0.5, 0.8] ;
max_wander = 20 ; % max amount that DV hoop id can wander per timepoint
disp('Identifying lobes...')
foldfn = fullfile(lobeDir, ['fold_locations_sphi' dvexten '_avgpts.mat']) ;
if exist(foldfn, 'file') && ~overwrite_folds
    disp('Loading lobes')
    % Save the fold locations as a mat file
    load(foldfn, 'ssfold', 'folds', 'ssfold_frac', 'ssmax', 'fold_onset', ...
        'rssfold', 'rssfold_frac', 'rssmax', 'rmax')
    
else
    [folds, ssfold, ssfold_frac, ssmax, rmax, fold_onset] = identifyLobes(xp.fileMeta.timePoints,...
            spcutMeshBase, guess123, max_wander, preview, 'avgpts') ;
    
    % Compute ringpath pathlength for results found using centerline
    disp('Converting folds to ringpath_ss locations...')
    [rssfold, rssfold_frac, rssmax] = rssFromFoldID(folds, xp.fileMeta.timePoints, spcutMeshBase) ;

    % Save the fold locations as a mat file
    save(foldfn, 'rssfold', 'rssfold_frac', 'rssmax', 'rmax', ...
        'ssfold', 'folds', 'ssfold_frac', 'ssmax', 'fold_onset')
end
clearvars guess123 maxwander

% Plot results as both avgpts and ringpath distances
fold_ofn = dir(fullfile(lobeDir, ['radii_folds' dvexten '_avgpts*.png'])) ;
if (length(fold_ofn) == length(timePoints)) || overwrite_foldims
    if save_ims
        disp('Plotting ss folds...')
        aux_plot_folds(folds, ssfold, ssfold_frac, ssmax, rmax, nU, ...
            xp.fileMeta.timePoints, lobeDir, dvexten, spcutMeshBase, ...
            'avgpts', overwrite_foldims)
    end
end
fold_ofn = dir(fullfile(lobeDir, ['radii_folds' dvexten '_avgpts*.png'])) ;
if (length(fold_ofn) == length(timePoints)) || overwrite_foldims
    if save_ims
        disp('Plotting rss folds...')
        aux_plot_folds(folds, rssfold, rssfold_frac, rssmax, rmax, nU, ...
            xp.fileMeta.timePoints, lobeDir, dvexten, spcutMeshBase, ...
            'ringpath', overwrite_foldims)
    end
end
disp('done')

%% RECOMPUTE WRITHE OF MEANCURVE CENTERLINES ==============================
% First compute using the avgpts (DVhoop means)
disp('Computing/Loading writhe...')
omit_endpts = 4 ;
wrfn = fullfile(lobeDir, ['writhe_sphi' dvexten '_avgpts.mat']) ;
if ~exist(wrfn, 'file') || overwrite_writhe
    [Wr, Wr_density, dWr, Length_t, clines_resampled] = ...
        aux_compute_writhe(clineDVhoopBase, xp.fileMeta.timePoints, true, omit_endpts, preview) ;
    
    % Save the fold locations as a mat file
    save(wrfn, 'Wr', 'Wr_density', 'dWr', 'Length_t', 'clines_resampled')
else
    load(wrfn, 'Wr', 'Wr_density', 'dWr', 'Length_t', 'clines_resampled')
end
Wr_style = 'Levitt' ;
tmpfn = fullfile(writheDir, ['writhe_' Wr_style '_vs_time_comparison_DVhoop.png']) ;
if ~exist(tmpfn, 'file') || overwrite_writhe
    % Compute ringpath pathlength for results found using centerline
    area_volume_fn = fullfile(meshDir, 'surfacearea_volume_stab.mat') ;
    aux_plot_writhe(xp.fileMeta.timePoints, clines_resampled, ...
        Wr, Wr_density, dWr, Length_t, writheDir, area_volume_fn, ...
        fold_onset, Wr_style, xyzlim, clineDVhoopBase, ...
        cylinderMeshCleanBase, rot, trans, resolution, omit_endpts, false)
    
end
disp('done')

%% Compute surface area and volume for each compartment ===================
lobe_dynamics_fn = fullfile(lobeDir, ['lobe_dynamics' dvexten '.mat']) ;
if exist(lobe_dynamics_fn, 'file') && ~overwrite_lobedynamics
    % Load length, surface area, and volume dynamics for each lobe
    disp('Loading lobe length, area, and volume...')
    load(lobe_dynamics_fn, 'length_lobes', 'area_lobes', 'volume_lobes')
    tp = xp.fileMeta.timePoints - min(fold_onset) ;
else
    [length_lobes, area_lobes, volume_lobes] = aux_compute_lobe_dynamics(folds, ssfold, ssmax, lobeDir, xp.fileMeta.timePoints, ...
        spcutMeshBase, nV, nU, rot, trans, resolution, xyzlim, colors, save_ims, overwrite_lobeims) ;
    % Save surface area and volume dynamics for each lobe
    save(fullfile(lobeDir, ['lobe_dynamics' dvexten '.mat']), ...
        'length_lobes', 'area_lobes', 'volume_lobes')
end
clearvars fold_ofn

%% plot length, area, and volume for each lobe ============================
lobe_dynamics_figfn = fullfile(lobeDir, ['lobe_dynamics' dvexten '.png']) ;
fig1exist = exist(lobe_dynamics_figfn, 'file') ;
% scaled version of same plot
lobe_dynamics_figfn = fullfile(lobeDir, ['lobe_dynamics' dvexten '_scaled.png']) ;
fig2exist = exist(lobe_dynamics_figfn, 'file') ;
if save_ims && (~fig2exist || ~fig2exist || overwrite_lobeims)
    disp('Plotting lobe dynamics...')
    aux_plot_lobe_dynamics(length_lobes, area_lobes, volume_lobes, ...
            timePoints, fold_onset, colors, lobe_dynamics_figfn)
else
    disp('Skipping lobe dynamics plot since it exists...')
end

%% Plot motion of avgpts at folds in yz plane over time ===================
aux_plot_avgptcline_lobes(folds, fold_onset, lobeDir, dvexten, save_ims, ...
    overwrite_lobeims, tp, timePoints, spcutMeshBase, clineDVhoopBase)
disp('done')

%% Plot motion of DVhoop at folds in yz plane over time ===================
aux_plot_constriction_DVhoops(folds, fold_onset, foldHoopImDir, dvexten, save_ims, ...
    overwrite_lobeims, tp, timePoints, spcutMeshBase, alignedMeshFullBase, ...
    normal_shift, rot, trans, resolution, colors, xyzlim)
disp('done')

%% SMOOTH MEAN CENTERLINE RADIUS ==========================================
% todo: rework this section so that it takes place after smoothing meshes
aux_smooth_avgptcline_radius_before_mesh_smoothing(overwrite_spcutMesh_smoothradii, ...
    timePoints, spcutMeshBase, clineDVhoopBase, radiusImDir, ...
    rot, trans, resolution, xyzlim, nU, nV)

%% Smooth the sphi grid meshes in time ====================================
% Check if all already exist
redo_meshsmooth = overwrite_spcutMeshSm ;
qq = 1 ;
disp('Checking if smoothed meshes already exist on file')
while ~redo_meshsmooth && qq < length(xp.fileMeta.timePoints)
    tt = xp.fileMeta.timePoints(qq) ;
    % Check that all meshes are saved--if not, declare that we must compute
    smfn = sprintf(spcutMeshSmBase, tt) ;
    smrsfn = sprintf(spcutMeshSmRSBase, tt) ;
    smrscfn = sprintf(spcutMeshSmRSCBase, tt) ;
    if ~exist(smfn, 'file') || ~exist(smrsfn, 'file') || ~exist(smrscfn, 'file')
        redo_meshsmooth = true ;
    end
    qq = qq + 1 ;
end

if redo_meshsmooth
    disp('Mesh smoothing does not exist on file, computing')
    % Load all spcutMesh objects 
    timePoints = xp.fileMeta.timePoints ;
    vM = zeros(length(timePoints), nU*nV, 3);
    nsmM = zeros(length(timePoints), nU*nV, 3) ;
    for i = 1:length(timePoints)
        tt = timePoints(i) ;
        % Load the spcutMesh for this timepoint
        disp(['Loading spcutMesh from disk... [t = ' num2str(tt) ']'])
        load(sprintf(spcutMeshBase, tt), 'spcutMesh') ;
        vM(i, :, :) = spcutMesh.v ;
    end
    disp('built v3d matrix')
    % Filter in time axis
    disp('Building tripulse filter equivalent to tripuls(-0.5:0.1:0.5)')
    tripulse = 0:0.2:1 ;
    tripulse = [tripulse, fliplr(tripulse(1:end-1))] ;
    tripulse = tripulse ./ sum(tripulse(:)) ;
    tripulse = reshape(tripulse, [length(tripulse), 1]) ;
    % linfilt = 0.1 * ones(10, 1, 1) ;
    % ellipsoid = fspecial3('ellipsoid', [5, 1, 1]) ;
    v3dsmM = imfilter(vM, tripulse, 'replicate') ;
    % vsmM = permute(vsmM, [2,1,3]) ;
    % nsmM = permute(nsmM, [2,1,3]) ;

    % Alternative is to use Gaussian filter but seems no padding option here:
    % smoothdata(vM, 1, 'gaussian', 10)
    close all

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save the smoothed meshes, then smoothed/rotated/scaled meshes 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for qq = 1:length(xp.fileMeta.timePoints)
        tt = xp.fileMeta.timePoints(qq) ;

        % First build normals from smoothed vertices as facenormals
        vqq = squeeze(v3dsmM(qq, :, :)) ;
        nqq = per_vertex_normals(vqq, spcutMesh.f, 'Weighting', 'angle') ;
        % pack it into a struct mesh
        mesh.v = vqq ;
        mesh.vn = nqq ;
        mesh.f = spcutMesh.f ;

        % Average normals with neighboring normals
        disp('Averaging normals with neighboring normals')
        nqq = average_normals_with_neighboring_vertices(mesh, 0.5) ;
        nsmM(qq, :, :) = nqq ;

        % Next save the mesh
        smfn = sprintf(spcutMeshSmBase, tt) ;
        smrsfn = sprintf(spcutMeshSmRSBase, tt) ;
        smrscfn = sprintf(spcutMeshSmRSCBase, tt) ;
        if ~exist(smfn, 'file') || ~exist(smrsfn, 'file') || ~exist(smrscfn, 'file') || overwrite_spcutMeshSm
            vqq = squeeze(v3dsmM(qq, :, :)) ;
            nqq = squeeze(nsmM(qq, :, :)) ;
            nsmM(qq, :, :) = nqq ./ vecnorm(nqq, 2, 2) ;

            % rotate and scale
            nqqrs = (rot * nqq')' ;
            vqqrs = ((rot * vqq')' + trans) * resolution;

            spcutMeshSm.f = spcutMesh.f ;
            spcutMeshSm.v = vqq ;
            spcutMeshSm.vn = nqq ;
            spcutMeshSm.u = spcutMesh.sphi ;
            spcutMeshSm.nU = spcutMesh.nU ;
            spcutMeshSm.nV = spcutMesh.nV ;
            spcutMeshSm.pathPairs = spcutMesh.pathPairs ;

            % Resave s,phi and their 3D embedding
            disp(['Saving ' sprintf(spcutMeshSmBase, tt)])
            save(sprintf(spcutMeshSmBase, tt), 'spcutMeshSm') ;

            % Also save rotated and scaled (RS) copy of the time-smoothed mesh
            spcutMeshSmRS = spcutMeshSm ;
            spcutMeshSmRS.v = vqqrs ;
            spcutMeshSmRS.vn = nqqrs ;
            % Resave s,phi and their 3D embedding
            save(sprintf(spcutMeshSmRSBase, tt), 'spcutMeshSmRS') ;
            clearvars vqq vqqrs

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CLOSED & GLUED SMOOTHED MESHES
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % To close the mesh, do the following:
            tmp = spcutMeshSmRS ;
            tmp.u = spcutMesh.sphi ;
            spcutMeshSmRSC = glueCylinderCutMeshSeam(tmp) ;
            save(sprintf(spcutMeshSmRSCBase, tt), 'spcutMeshSmRSC') ;

            % check it
            % fig = figure ;
            % % triplot(spcutMeshSmRSC.f, spcutMeshSmRSC.u(:, 1), spcutMeshSmRSC.u(:, 2))
            % trisurf(spcutMeshSmRSC.f, spcutMeshSmRSC.v(:, 1),...
            %     spcutMeshSmRSC.v(:, 2), spcutMeshSmRSC.v(:, 3))
            % To fully close the mesh use:
            % anewpt = mean(spcutMeshSmRS.v(1:nV:end, :), 1)
            % pnewpt = mean(spcutMeshSmRS.v(nU:nV:end, :), 1)
            % waitfor(fig)
        end
    end
else
    disp('Mesh smoothing already exists on file, loading...')
    timePoints = xp.fileMeta.timePoints ;
    v3dsmM = zeros(length(timePoints), nU*nV, 3);
    nsmM = zeros(length(timePoints), nU*nV, 3) ;
    
    % Load each mesh into v3dsmM and nsmM    
    for qq = 1:length(xp.fileMeta.timePoints)
        tt = xp.fileMeta.timePoints(qq) ;
        load(sprintf(spcutMeshSmBase, tt), 'spcutMeshSm') ;
        v3dsmM(qq, :, :) = spcutMeshSm.v ;
        nsmM(qq, :, :) = spcutMeshSm.vn ;
    end
end
disp('done smoothing meshes in time')

%% Plot the time-smoothed meshes
pdir = ensureDir(fullfile(sphiSmRSPhiImDir, 'perspective')) ;
ddir = ensureDir(fullfile(sphiSmRSPhiImDir, 'dorsal')) ;
vdir = ensureDir(fullfile(sphiSmRSPhiImDir, 'ventral')) ;
ldir = ensureDir(fullfile(sphiSmRSPhiImDir, 'latL')) ;
rdir = ensureDir(fullfile(sphiSmRSPhiImDir, 'latR')) ;
for qq = 1:length(timePoints)
    tt = timePoints(qq) ;
    disp(['checking figures for time-smoothed meshes: t = ' num2str(tt)])
    % prep directories
    fig1fn = fullfile(sphiSmRSImDir, [sprintf('%04d', tt ) '.png']) ;
    fp0fn = fullfile(pdir, [sprintf('%04d', tt ) '.png']) ;
    fp1fn = fullfile(ddir, ['dorsal_' sprintf('%04d', tt ) '.png']) ;
    fp2fn = fullfile(vdir, ['ventral_' sprintf('%04d', tt ) '.png']) ;
    fp3fn = fullfile(ldir, ['latL_' sprintf('%04d', tt ) '.png']) ;
    fp4fn = fullfile(rdir, ['latR_' sprintf('%04d', tt ) '.png']) ;
    e0 = ~exist(fig1fn, 'file') ;
    e1 = ~exist(fp0fn, 'file') ;
    e2 = ~exist(fp1fn, 'file') || ~exist(fp2fn, 'file') ;
    e3 = ~exist(fp3fn, 'file') || ~exist(fp4fn, 'file') ;
    if e0 || e1 || e2 || e3 || overwrite_SmRSIms
        disp(['Consider smoothed sphi gridmesh for time ' num2str(tt)])
        % Load the spcutMesh for this timepoint
        vqq = squeeze(v3dsmM(qq, :, :)) ;
        nqq = squeeze(nsmM(qq, :, :)) ;

        % rotate and scale
        nqqrs = (rot * nqq')' ;
        vqq = ((rot * vqq')' + trans) * resolution;
    end
    
    % Plot embedding RS colored by normal vector in y direction
    if e0 || overwrite_SmRSIms
        fig = figure('visible', 'off') ;
        % Color figure with normals taken via faceNormals  
        trisurf(spcutMeshSm.f, vqq(:, 1), vqq(:, 2), vqq(:, 3), ...
            nqqrs(:, 2), 'EdgeColor', 'none')
        axis equal
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        xlabel('x [\mum]')
        ylabel('y [\mum]')
        zlabel('z [\mum]')
        title(['Smoothed mesh embedding, t=' num2str(tt)])
        disp(['Saving figure: ' fig1fn])
        saveas(fig, fig1fn)
        close all
    end
    
    % Plot embedding colored by phi
    if e1 || e2 || e3 || overwrite_SmRSIms 
        fig = figure('visible', 'off') ;
        % Plot embedding in 3d color coding (color=phi)
        pc = (1:nV) .* ones(nU, nV) / nV ;
        trisurf(spcutMeshSm.f, vqq(:, 1), vqq(:, 2), vqq(:, 3), ...
            pc(:), 'EdgeColor', 'none')
        axis equal
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        xlabel('x [\mum]')
        ylabel('y [\mum]')
        zlabel('z [\mum]')
        c = colorbar() ;
        c.Label.String = '\phi / 2\pi' ;
        title(['Smoothed mesh embedding, t=' num2str(tt)])
        disp(['Saving figure: ' fp0fn])
        saveas(fig, fp0fn)
        view(0, 90)
        saveas(fig, fp1fn)
        view(0, 270)
        saveas(fig, fp2fn)
        view(0, 0)
        saveas(fig, fp3fn)
        view(0, 180)
        saveas(fig, fp4fn)
        close all
    end
end
disp('done')
clearvars fig vqq nqqrs e0 e1 e2 e3 pdir ddir vdir rdir ldir

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE MEAN AND GAUSSIAN CURVATURES OF SMOOTHED MESHES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
getderivatives=0;
disp('Computing/Loading Curvatures...')
for tt = xp.fileMeta.timePoints
    Kfns = {fullfile(KSmDir, 'latleft', sprintf('gausscurv_latleft_%06d.png', tt)), ...
        fullfile(KSmDir, 'dorsal', sprintf('gausscurv_dorsal_%06d.png', tt)), ...
        fullfile(KSmDir, 'latright', sprintf('gausscurv_latright_%06d.png', tt)), ...
        fullfile(KSmDir, 'ventral', sprintf('gausscurv_ventral_%06d.png', tt)) };
    Hfns = {fullfile(HSmDir, 'latleft', sprintf('meancurv_latleft_%06d.png', tt)), ...
        fullfile(HSmDir, 'dorsal', sprintf('meancurv_dorsal_%06d.png', tt)), ...
        fullfile(HSmDir, 'latright', sprintf('meancurv_latright_%06d.png', tt)),...
        fullfile(HSmDir, 'ventral', sprintf('meancurv_ventral_%06d.png', tt))} ;
    KHfn = fullfile(KHSmDir, sprintf('gauss_mean_curvature_%06d.mat', tt)) ;
    if ~exist(KHfn, 'file') || overwrite_curvatures
        % load glued / closed mesh for curvature computation
        load(sprintf(spcutMeshSmRSCBase, tt), 'spcutMeshSmRSC') ;
        fv.faces = spcutMeshSmRSC.f ;
        fv.vertices = spcutMeshSmRSC.v ;
        PrincipalCurvatures = GetCurvatures( fv, getderivatives);
        gaussCurv = (PrincipalCurvatures(1,:).*PrincipalCurvatures(2,:))';
        meanCurv = -0.5 * (PrincipalCurvatures(1,:) + PrincipalCurvatures(2,:))';

        opts.outfn = Kfns ;
        opts.label = 'Gaussian curvature, $K$' ;
        opts.sscale = 0.01 ;
        opts.cbarPosition = [.85 .333 .02 .333] ;
        opts.xlim = xyzlim(1, :) ;
        opts.ylim = xyzlim(2, :) ;
        opts.zlim = xyzlim(3, :) ;
        scalarFieldOnSurface(fv.faces, fv.vertices, gaussCurv, opts)
        opts.outfn = Hfns ;
        opts.sscale = 0.3 ;
        opts.label = 'mean curvature, $H$' ;
        scalarFieldOnSurface(fv.faces, fv.vertices, meanCurv, opts)
        % Save the curvatures as a mat file
        save(KHfn, 'meanCurv', 'gaussCurv')
    else
        disp('Curvatures already computed & saved on disk')
    end
end
error('here after curvature calc')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Redo Pullbacks with time-smoothed meshes ===============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Create pullback using S,Phi coords with time-averaged Meshes')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for qq = 1:length(xp.fileMeta.timePoints)
    tt = xp.fileMeta.timePoints(qq) ;
    disp(['checking pullbacks for time-smoothed meshes: t = ' num2str(tt)])
    
    % Load time-smoothed mesh
    load(sprintf(spcutMeshSmBase, tt), 'spcutMeshSm') ;
    
    %--------------------------------------------------------------
    % Generate Output Image File
    %--------------------------------------------------------------
    imfn_spsm = sprintf( fullfile( imFolder_spsm, [fileNameBase, '.tif']), tt ) ;
    imfn_rsm  = sprintf( fullfile( imFolder_rsm, [fileNameBase, '.tif']), tt ) ;
    pullbacks_exist1 = exist(imfn_spsm, 'file') ;
    pullbacks_exist2 = exist(imfn_rsm, 'file') ;
    if ~pullbacks_exist1 || ~pullbacks_exist2 || overwrite_pullbacks
        QS.getCurrentData()
        IV = QS.currentData.IV ;
    end

    if ~exist(imfn_spsm, 'file') || overwrite_pullbacks
        fprintf(['Generating SP output image for sm mesh: ' imfn_spsm]);
        % Assigning field spcutMesh.u to be [s, phi] (ringpath
        % and azimuthal angle)
        aux_generate_orbifold( spcutMeshSm, a_fixed, IV, imfn_spsm)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save relaxed image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if ~exist(imfn_rsm, 'file') || overwrite_pullbacks
        disp('Generating relaxed image for sm mesh...')
        tmp = spcutMesh.sphi ;
        tmp(:, 1) = tmp(:, 1) / max(tmp(:, 1)) ;
        arspsm = minimizeIsoarealAffineEnergy( spcutMeshSm.f, spcutMeshSm.v, tmp );
        aux_generate_orbifold(spcutMeshSm, arspsm, IV, imfn_rsm)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TILE SMOOTHED IMAGES IN Y AND RESAVE ===================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imDirs = {imFolder_spsm, imFolder_rsm, } ;
imDirs_e = {imFolder_spsm_e, imFolder_rsm_e, } ;
ntiles = 50 ;   %   The number of bins in each dimension for histogram equilization for a
                %   square original image. That is, the extended image will have (a_fixed *
                %   ntiles, 2 * ntiles) bins in (x,y).
for qq = 1:length(imDirs)
    options.histeq = true ;
    options.a_fixed = a_fixed ;
    options.ntiles = ntiles ;
    options.overwrite = overwrite_pullbacks;
    extendImages(imDirs{qq}, imDirs_e{qq}, fileNameBase, options)
    disp(['done ensuring extended tiffs for ' imDirs{qq} ' in ' imDirs_e{qq}])
end
clearvars qq ntiles imDirs_e imDirs
disp('done')

%% CREATE PULLBACK STACKS =================================================
fprintf('Create pullback nstack using S,Phi coords with time-averaged Meshes \n');
% Load options
if ~exist(optionfn, 'file') || overwriteSmSPCutMeshStackOptions
    spcutMeshSmStackOptions.layer_spacing = 1 / resolution ; % pixel resolution matches xy
    spcutMeshSmStackOptions.n_outward = 10 ;
    spcutMeshSmStackOptions.n_inward = 8 ;
else
    load(optionfn, 'smSPCutMeshStackOptions')
end
spcutMeshSmStackOptions.overwrite = overwrite_spcutMeshSmStackOptions ;
QS.generateSPCutMeshSmStack(spcutMeshSmStackOptions)

%% TRAIN IN ILASTIK ON MIDGUT TISSUE TO GET THICKNESS
% Read thickness training output

%% MEASURE THICKNESS ======================================================
% Measure thickness from images of imFolder_spsm_e2 (extendedLUT_smoothed)
QS.measureThickness(thicknessOptions)

%% DUMP OR LOAD HERE [break]
clearvars fold_ofn fig1exist fig2exist id idx mcline nsmM prevcline tmp
clearvars fig1exist fig1fn fig2exist fp1fn fp2fn fp3fn fp4fn IV 
clearvars n_inward n_outward TV2D TV3D 
clearvars layer_spacing 
dumpfn = fullfile(meshDir, 'orbifold_dump_before_piv.mat') ;
save(dumpfn)
load(dumpfn)
clearvars dumpfn

