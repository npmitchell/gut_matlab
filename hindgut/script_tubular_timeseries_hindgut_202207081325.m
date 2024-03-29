%% EXAMPLE MASTER PIPELINE FOR TIME SERIES DATA (3D + time)
% by NPMitchell & Dillon Cislo

% This is a pipeline to analyze dynamic tube-like surfaces in 3D data.
% A tube-like surface is one that is either cylindrical in topology or
% elongated and spherical in topology. If the initial mesh is a spherical
% topology, then two 'endcaps' will be truncated in order to transform it 
% to a cylindrical topology.

% Here we analyze the hindgut dataset bynGAL4klarCAAXmChHGFP
% 202207081325_1p4um_1mW1mW_4views_0p25ms_120s

%% Clear workspace ========================================================
% We start by clearing the memory and closing all figures
clear; close all; clc;
cd /mnt/crunch/bynGAL4klarCAAXmChHGFP/202207081325_1p4um_1mW1mW_4views_0p25ms_120s/
% cd /mnt/crunch/bynGAL4klarCAAXmChHGFP/202207131248_1p4um_2p5mW1mW_4views_0p25ms_60s/deconvolved_16bit
zwidth =1 ;
% cd /mnt/data/48Ygal4-UAShistRFP/201904031830_great/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2/data/deconvolved_16bit/
% zwidth = 1 ;
% cd /mnt/data/handGAL4klarHandGFPhistGFP/202105072030_1mWGFP/deconvolved_16bit/
% % cd /mnt/data/mef2GAL4klarUASCAAXmChHiFP/202003151700_1p4um_0p5ms3msexp/
% zwidth = 2 ;


dataDir = cd ;

%% ADD PATHS TO THIS ENVIRONMENT ==========================================
origpath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(origpath))
% addpath(fileparts(origpath))
addpath(genpath('/mnt/data/code/tubular/'))
addpath(genpath('/mnt/data/code/tubular/utility'))
addpath(genpath('/mnt/data/code/gptoolbox'))
addpath('/mnt/data/code/tubular/TexturePatch')
addpath('/mnt/data/code/tubular/DECLab')
addpath(fullfile('/mnt/data/code/tubular/utility','plotting'))
addpath(fullfile('/mnt/data/code/tubular/utility','plotting'))
% go back to the data
cd(dataDir)




%% DEFINE NEW MASTER SETTINGS
if ~exist('./masterSettings.mat', 'file') 
    % Metadata about the experiment
    stackResolution = [.2619 .2619 .2619] ;  % resolution in spaceUnits per pixel
    nChannels = 2 ;             % how many channels is the data (ex 2 for GFP + RFP)
    channelsUsed = [1 2] ;          % which channels are used for analysis
    timePoints = 0:92;       % timepoints to include in the analysis
    ssfactor = 4 ;              % subsampling factor
    flipy = false ;             % whether the data is stored inverted relative to real position in lab frame
    timeInterval = 2 ;          % physical interval between timepoints
    timeUnits = 'min' ;         % physical unit of time between timepoints
    spaceUnits = [char(956) 'm'] ;     % physical unit of time between timepoints
    scale = [0.04 0.2] ;      % scale for conversion to 16 bit, one value for each channel
    file32Base = 'TP%d_Ch%d_Ill0_Ang0,45,90,135,180,225,270,315.tif'; 
    % file32Base = 'TP%d_Ch%d_Ill0_Ang0,60,120,180,240,300.tif'; 
    fn = 'Time_%06d_c%d_stab';
    fn_prestab = 'Time_%06d_c%d.tif';
    set_preilastikaxisorder = 'xyzc' ; % data axis order for subsampled h5 data (ilastik input)
    swapZT = 0 ;                % whether to swap the z and t dimensions
    masterSettings = struct('stackResolution', stackResolution, ...
        'scale', scale, ...
        'file32Base', 'TP%d_Ch%d_Ill0_Ang0,45,90,135,180,225,270,315.tif', ...
        'fn_prestab', 'Time_%06d_c%d.tif', ...
        'nChannels', nChannels, ...
        'channelsUsed', channelsUsed, ...
        'timePoints', timePoints, ...
        'ssfactor', ssfactor, ...
        'flipy', flipy, ...
        'timeInterval', timeInterval, ...
        'timeUnits', timeUnits, ...
        'spaceUnits', spaceUnits, ...
        'fn', fn,...
        'swapZT', swapZT, ...
        'set_preilastikaxisorder', set_preilastikaxisorder, ...
        'nU', 100, ...  
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
    timeInterval = masterSettings.timeInterval ;  % physical interval between timepoints
    timeUnits = masterSettings.timeUnits ; % physical unit of time between timepoints
    spaceUnits = masterSettings.spaceUnits ; % unit of distance of full resolution data pixels ('$\mu$m')
    fn = masterSettings.fn ;
    fn_prestab = masterSettings.fn_prestab ;
    set_preilastikaxisorder = masterSettings.set_preilastikaxisorder ;
    swapZT = masterSettings.swapZT ;
    nU = masterSettings.nU ;
    nV = masterSettings.nV ;
    scale = masterSettings.scale ;
end



%dir32bit = fullfile(dataDir, 'deconvolved_32bit') ;
dir32bit = fullfile(dataDir, 'deconvolved12iter') ;
dir16bit = fullfile(dataDir, 'deconvolved_16bit') ;
dir16bit_prestab = fullfile(dir16bit, 'data_pre_stabilization') ;

%% END OF EXPERIMENT METADATA =============================================
% =========================================================================
% =========================================================================
addpath(genpath('/mnt/data/code/gut_matlab/plotting/'))
addpath(genpath('/mnt/data/code/gut_matlab/tiff_handling/'))
% 
% % Skip if already done
% convert32to16bit(timePoints, scale, dir32bit, dir16bit_prestab,...
%     file32Base, fn_prestab, channelsUsed)
% 
% %% -III. make MIPs for 16bit images
% % Skip if already done
% mipDir = fullfile(dir16bit_prestab, 'mips') ;
% Options.overwrite_mips = false ;
% Options.scale = -1 ; % do NOT rescale intensities during intensity projection
% Options.channels = [1 2] ;
% makeMips(timePoints, dir16bit_prestab, fn_prestab, mipDir, Options)
% 
% % -III. make saggital MIPs for 16bit images
% % Skip if already done
% subMipDir = fullfile(dir16bit_prestab, 'substack') ;
% Options.overwrite_mips = false ;
% Options.scale = -1 ; % do NOT rescale intensities during intensity projection
% Options.overlayColors = true ;
% Options.width = 5 ;
% Options.channels = [1 2] ;
% Options.overwrite_overlays = false ;
% makeSubStackMips(timePoints, dir16bit_prestab, fn_prestab, subMipDir, Options)
% 
% %% Collate multiple colors into one TIFF pre-stabilization
% fileNameIn = fullfile(dir16bit_prestab, fn_prestab) ;
% fileNameOut = fullfile(dir16bit_prestab, 'Time_%06d.tif') ;
% collateColors(fileNameIn, fileNameOut, timePoints, channelsUsed) ; 
% 
% % % Rename stab to prestab -- hack
% % % fns = fullfile('./deconvolved_16bit/Time*stab')
% % % for qq=1:length(fns)
% % %     command = ['mv ' fullfile(fns.folder, fns.name) fullfile(dir16bit, fns.name
% % % end
% 
% % -III. make MIPs for 16bit images
% % Skip if already done
% mipDir = fullfile(dir16bit_prestab, 'mips') ;
% Options.overwrite_mips = false ;
% Options.scale = -1 ; % do NOT rescale intensities during intensity projection
% makeMips(timePoints, dir16bit_prestab, fn_prestab, mipDir, Options)
% 
% %%  -II. stabilize images, based on script stabilizeImagesCorrect.m
% % Skip if already done
% % name of directory to check the stabilization of mips
% mips_stab_check = fullfile(mipDir, 'stab_check') ;
% mipoutdir = fullfile(mipDir, 'mips_stab') ;
% % Choose bit depth as typename
% typename = 'uint16' ;
% % Give file names for I/O
% fileNameIn = fullfile(dir16bit_prestab, fn_prestab) ;
% fileNameOut = fullfile(dir16bit, [fn '.tif']) ;
% rgbName = [fn '.png'] ;
% typename = 'uint16' ; 
% Options = struct(); 
% Options.overwrite_mips = false ;
% Options.overwrite_tiffs = false ;
% Options.im_intensity = 1 ; % 0.01 ;
% Options.imref_intensity = 1 ; % 0.005 ;
% Options.stabChannel = 1 ; 
% stabilizeImages(fileNameIn, fileNameOut, rgbName, typename, ...
%     timePoints, timePoints, timePoints(1), ...
%     mipDir, mipoutdir, mips_stab_check, Options)
% Options.stabChannel = 2 ; 
% stabilizeImages(fileNameIn, fileNameOut, rgbName, typename, ...
%     timePoints, timePoints, timePoints(1), ...
%     mipDir, mipoutdir, mips_stab_check, Options)
% 
% 
% %% Collate multiple colors into one TIFF post-stabilization
% fileNameIn = fullfile(dir16bit, 'Time_%06d_c%d_stab.tif') ;
% fileNameOut = fullfile(dir16bit, 'Time_%06d.tif') ;
% collateColors(fileNameIn, fileNameOut, timePoints, channelsUsed) ; 
% 
% %%   -I. master_gut_timeseries_prestab_for_training.m
% % Skip if already done
% % cd(dir16bit)
% % dataDir = cd ;
% % masterSettings.dir16bit_prestab = dir16bit_prestab ;
% % addpath('/mnt/data/code/gut_matlab/master_pipeline/')
% % makeH5SeriesForTraining(masterSettings)

%% 
cd(dir16bit)
dataDir = dir16bit ;
fn = 'Time_%06d' ;


%% Define settings for experiment
overwrite = false ;
if ~exist(fullfile(dataDir, 'xp.mat'), 'file') || overwrite
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % PART 1: Define the metadata for the project
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd(dir16bit)
    projectDir = dataDir ;

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
    fileMeta.dataDir            = dir16bit;
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


    % % SET DETECTION OPTIONS ==================================================
    % Load/define the surface detection parameters
    msls_detOpts_fn = fullfile(projectDir, 'msls_detectOpts.mat') ;
    if exist(msls_detOpts_fn, 'file') && false
        disp('loading detectOptions')
        load(msls_detOpts_fn, 'detectOptions')
    else
        outputfilename_ply='mesh_ms_' ;
        outputfilename_ls='msls_' ;
        outputfilename_smoothply = 'mesh_' ;
        ms_scriptDir='/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/' ;   
        init_ls_fn = 'msls_initguess' ;
        meshlabCodeDir = '/mnt/data/code/meshlab_codes/';
        mlxprogram = fullfile(meshlabCodeDir, ...
            'laplace_surface_rm_resample30k_reconstruct_LS3_1p2pc_ssfactor4.mlx') ;
        prob_searchstr = '_stab_Probabilities.h5' ;
        preilastikaxisorder = set_preilastikaxisorder; ... % axis order in input to ilastik as h5s. To keep as saved coords use xyzc
        ilastikaxisorder= 'yxzc'; % default='cxyz'; ... % axis order as output by ilastik probabilities h5
        imsaneaxisorder = 'xyzc'; ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc

        % Name the output mesh directory --------------------------------------
        mslsDir = [fullfile(projectDir, 'tubular_output') filesep];

        % Surface detection parameters ----------------------------------------
        detectOptions = struct( 'channel', 1, ...
            'ssfactor', 4, ...
            'niter', 35,...
            'niter0', 160, ...
            'pre_pressure', -5, ...
            'pre_tension', 0, ...
            'pressure', 0, ...
            'tension', 0.5, ...
            'post_pressure', 2, ...
            'post_tension', 3, ...
            'exit_thres', 1e-7, ...
            'foreGroundChannel', 2, ...
            'fileName', fn, ...
            'mslsDir', mslsDir, ...
            'ofn_ls', outputfilename_ls, ...
            'ofn_ply', outputfilename_ply,...
            'ms_scriptDir', ms_scriptDir, ...
            'timepoint', timePoints(1), ...
            'zdim', 2, ...
            'ofn_smoothply', outputfilename_smoothply, ...
            'mlxprogram', mlxprogram, ...
            'init_ls_fn', init_ls_fn, ... % set to none to load prev tp
            'run_full_dataset', projectDir,... % projectDir, ... % set to 'none' for single tp
            'radius_guess', 40, ...
            'dset_name', 'exported_data',...
            'center_guess', '200,75,75',... % xyz of the initial guess sphere ;
            'save', true, ... % whether to save images of debugging output
            'plot_mesh3d', false, ...
            'dtype', 'h5',...
            'mask', 'none',...
            'mesh_from_pointcloud', false, ...
            'prob_searchstr', prob_searchstr, ...
            'preilastikaxisorder', preilastikaxisorder, ... 
            'ilastikaxisorder', ilastikaxisorder, ... 
            'physicalaxisorder', imsaneaxisorder, ... 
            'include_boundary_faces', true, ...
            'smooth_with_matlab', -1) ; % set this to >0 to use matlab laplacian filter instead of meshlab

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
    mslsDir = detectOptions.mslsDir ;

    % % Define Experiment as struct
    xp = struct('fileMeta', fileMeta, ...
        'expMeta', expMeta, 'detectOptions', detectOptions) ;
    disp('done')
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % PART 2: TubULAR -- surface parameterization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Now we have 3d data volumes and surfaces. Define a TubULAR object. 
    % To visualize data on these surfaces and compute how these surfaces deform
    % we now define TubULAR object.
    nU = masterSettings.nU ;
    nV = masterSettings.nV ;
    opts = struct() ;
    opts.meshDir = mslsDir ;        % Directory where meshes reside
    opts.flipy = flipy ;            % Set to true if data volume axes are inverted in chirality wrt physical lab coordinates
    opts.timeInterval = timeInterval ; % Spacing between adjacent timepoints in units of timeUnits 
    opts.timeUnits = timeUnits ;    % units of time, so that adjacent timepoints are timeUnits * timeInterval apart
    opts.spaceUnits = spaceUnits ;  % Units of space in LaTeX, for ex '$mu$m' for micron
    opts.nU = nU ;                  % How many points along the longitudinal axis to sample surface
    opts.nV = nV ;                  % How many points along the circumferential axis to sample surface
    opts.normalShift = 10 ;         % Additional dilation acting on surface for texture mapping
    opts.a_fixed = 2.0 ;            % Fixed aspect ratio of pullback images. Setting to 1.0 is most conformal mapping option.
    opts.adjustlow = 1.00 ;         % floor for intensity adjustment
    opts.adjusthigh = 99.9 ;        % ceil for intensity adjustment (clip)
    opts.phiMethod = 'curves3d' ;   % Method for following surface in surface-Lagrangian mapping [(s,phi) coordinates]
    opts.lambda_mesh = 0.00 ;       % Smoothing applied to the mesh before DEC measurements
    opts.lambda = 0.0 ;             % Smoothing applied to computed values on the surface
    opts.lambda_err = 0.0 ;         % Additional smoothing parameter, optional
    opts.zwidth = zwidth ;
    opts.nmodes = 7 ;
    % opts.t0 = xp.fileMeta.timePoints(1) ;   % reference timepoint used to define surface-Lagrangian and Lagrangian measurements
    % opts.t0 = 123 ;
    % opts.t0 = 37 ;
    % opts.t0 = 1 ;

    disp('saving xp struct and opts to disk')
    save(fullfile(dataDir, 'xp.mat'), 'xp', 'opts')
else
    disp('loading xp struct from disk')
    load(fullfile(dataDir, 'xp.mat'), 'xp', 'opts')
end

%% TubULAR class instance
disp('defining TubULAR class instance (tubi= tubular instance)')
tubi = TubULAR(xp, opts) ;
disp('done defining TubULAR instance')

%% 
tubi.prepareIlastik()


%% 
tubi.xp.detectOptions.niter0 = 50 ;
tubi.xp.detectOptions.post_pressure = 2 ;
tubi.xp.detectOptions.pressure = 0.2 ;
tubi.xp.detectOptions.tension = 0.1 ;
tubi.xp.detectOptions.radius_guess = 15 ;
tubi.xp.detectOptions.center_guess = 'click' ;
tubi.xp.detectOptions.previewIsosurface = false ;
tubi.getMeshes()

%% Obtain APDV coordinates of the surface. 
% There are two options for obtaining these coordinates. 
%   1. Automatically determine A and P by the extremal points of the
%   surface mesh along the elongated axis of the mesh, and define DV as
%   pointing perpendicular to this.
%   2. Train in iLastik for an anterior spot in 3d (A), a posterior spot in
%   3d (P), and a spot which is dorsal to the line connecting A and P. Any
%   dorsal point is fine, as long as it points dorsal to the AP axis
%   defined by A and P. See picture below.
%
%     example:
%                 D
%                 |
%                 |
%       A -----------------P
%
% Here we use option 2. We must prepare APDV ilastik training first outside
% MATLAB.
% Train on anterior (A), posterior (P), background (B), and 
% dorsal anterior (D) location in different iLastik channels by having a
% blob centered on a point that you wish to identify as each label (in 3D).
% anteriorChannel, posteriorChannel, and dorsalChannel specify the iLastik
% training channel that is used for each specification.
% Name the h5 file output from iLastik as ..._Probabilities_apcenterline.h5
% Training for dorsal (D) is only needed at the reference time point, t0,
% because that's the only one that's used. 
%
% A dorsal blob for the gut is marked at the site where the gut closes,
% with 48YGAL4-expressing cells form a seam.
% Posterior is at the rear of the yolk, where the endoderm closes, for 
% apical surface training.
% Anterior is at the junction of the midgut with the foregut.

%% Define global orientation frame (for viewing in canonical frame)
% Compute APDV coordinate system
alignAPDVOpts = struct() ;
alignAPDVOpts.overwrite = true ;
% alignAPDVOpts.use_APDVAxes = [2 3] ;
% alignAPDVOpts.flipAP = true ;
% alignAPDVOpts.flipDV = true ;
alignAPDVOpts.use_iLastik = true ;
alignAPDVOpts.ilastikOutputAxisOrder = 'yxzc' ;
tubi.computeAPDVCoords(alignAPDVOpts) ;
disp('done')

%% Select the endcaps for the centerline computation (A and P) and a point
% along which we will form a branch cut for mapping to the plane (D).
apdvOpts = struct() ;
apdvOpts.overwrite = true ;
apdvOpts.use_iLastik=true ;
apdvOpts.ilastikOutputAxisOrder = 'yxzc' ;
[apts_sm, ppts_sm] = tubi.computeAPDpoints(apdvOpts) ;

%% Align the meshes in the APDV global frame & plot them
alignAPDVOpts.overwrite = true ;
tubi.alignMeshesAPDV(alignAPDVOpts) ;

disp('done')

%% PLOT ALL TEXTURED MESHES IN 3D (OPTIONAL: this is SLOW) ================
% Establish texture patch options
metadat = struct() ;
metadat.reorient_faces = false ;            % set to true if some mesh normals may be inverted (requires gptoolbox if true)
metadat.normal_shift = tubi.normalShift ;   % normal push, in pixels, along normals defined in data XYZ space
metadat.texture_axis_order = [2 1 3] ;      % texture space sampling. If the surface and dataspace have axis permutation, enter that here
Options.PSize = 5 ;          % Psize is the linear dimension of the grid drawn on each triangular face. Set PSize > 1 for refinement of texture on each triangle of the surface triangulation. Higher numbers are slower but give more detailed images.
Options.numLayers = [0, 5];  % how many layers to MIP over/bundle into stack, as [outward, inward]
Options.layerSpacing = 1.2 ;   % Distance between layers over which we take MIP, in pixels, 

Options.smoothIter = 3 ;
metadat.smoothing_lambda = 0.1 ;

metadat.perspective_angle = [-25,40] ;
metadat.overwrite = true ;
Options.falseColors = [ 0 1 1; 1 0 0] ;
metadat.plot_ventral = false ;
metadat.plot_dorsal = false ;
metadat.plot_lateral1 = false ;
metadat.plot_lateral2 = false ;
% Plot on surface for all timepoints 
tubi.plotSeriesOnSurfaceTexturePatch(metadat, Options)

%% EXTRACT CENTERLINES
% Note: these just need to be 'reasonable' centerlines for topological
% checks on the orbifold cuts. Therefore, use as large a resolution ('res')
% as possible that still forms a centerline passing through the mesh
% surface, since the centerline computed here is just for constraining the 
% mapping to the plane.
cntrlineOpts.overwrite = false ;         % overwrite previous results
cntrlineOpts.overwrite_ims = false ;     % overwrite previous results
cntrlineOpts.weight = 0.1;               % for speedup of centerline extraction. Larger is less precise
cntrlineOpts.exponent = 1.0 ;            % how heavily to scale distance transform for speed through voxel
cntrlineOpts.res = 4.0 ;                 % resolution of distance tranform grid in which to compute centerlines
cntrlineOpts.preview = false ;           % preview intermediate results
cntrlineOpts.reorient_faces = false ;    % not needed for our well-constructed meshes
cntrlineOpts.dilation = 0 ;              % how many voxels to dilate the segmentation inside/outside before path computation
% Note: this can take about 400s per timepoint for res=2.0, so use as big a 
%   res value as possible.
%
tubi.generateFastMarchingCenterlines(cntrlineOpts)
disp('done with centerlines')

%% Identify anomalies in centerline data
idOptions.ssr_thres = 15 ;  % distance of sum squared residuals in um as threshold for removing spurious centerlines
tubi.cleanFastMarchingCenterlines(idOptions) ;
disp('done with cleaning up centerlines')

%% Cylinder cut mesh --> transforms a topological sphere into a topological cylinder
% Look for options on disk. If not saved, define options.
if ~exist(tubi.fileName.endcapOptions, 'file')
    endcapOpts = struct( 'adist_thres', 20, ...  % 20, distance threshold for cutting off anterior in pix
                'pdist_thres', 14, ...  % 15-20, distance threshold for cutting off posterior in pix
                'tref', tubi.t0) ;  % reference timepoint at which time dorsal-most endcap vertices are defined
    tubi.setEndcapOptions(endcapOpts) ;
    % Save the options to disk
    tubi.saveEndcapOptions() ;
else
    % load endcapOpts
    tubi.loadEndcapOptions() ;
    endcapOpts = tubi.endcapOptions ;
end

methodOpts.overwrite = true ;
methodOpts.save_figs = true ;   % save images of cutMeshes along the way
methodOpts.preview = false  ;     % display intermediate results
tubi.sliceMeshEndcaps(endcapOpts, methodOpts) ;

%% Clean Cylinder Meshes
% This removes "ears" from the endcaps of the tubular meshes (cylindrical
% meshes)
cleanCylOptions = struct() ;
cleanCylOptions.overwrite = false ;
tubi.cleanCylMeshes(cleanCylOptions)
disp('done cleaning cylinder meshes')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORBIFOLD -> begin populating tubi.dir.mesh/gridCoords_nUXXXX_nVXXXX/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overwrite = false ;
% Iterate Through Time Points to Create Pullbacks ========================
for tt = tubi.xp.fileMeta.timePoints
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    tidx = tubi.xp.tIdx(tt);
    
    % Load the data for the current time point ------------------------
    tubi.setTime(tt) ;
    
    %----------------------------------------------------------------------
    % Create the Cut Mesh
    %----------------------------------------------------------------------
    cutMeshfn = sprintf(tubi.fullFileBase.cutMesh, tt) ;
    cutPathfn = sprintf(tubi.fullFileBase.cutPath, tt) ;
    if ~exist(cutMeshfn, 'file') || ~exist(cutPathfn, 'file') || overwrite
        if exist(cutMeshfn, 'file')
            disp('Overwriting cutMesh...') ;
        else
            disp('cutMesh not found on disk. Generating cutMesh... ');
        end
        options = struct() ;
        tubi.generateCurrentCutMesh(options)
        disp('Saving cutP image')
        % Plot the cutPath (cutP) in 3D
        tubi.plotCutPath(tubi.currentMesh.cutMesh, tubi.currentMesh.cutPath)
        compute_pullback = true ;
    else
        compute_pullback = false ;
    end
    
    tubi.getCurrentUVCutMesh() ;
    
    spcutMeshOptions = struct() ;
    spcutMeshOptions.t0_for_phi0 = tubi.t0set() ;  % which timepoint do we define corners of pullback map
    spcutMeshOptions.save_phi0patch = false ;
    spcutMeshOptions.iterative_phi0 = false ;
    spcutMeshOptions.smoothingMethod = 'none' ;
    tubi.plotting.preview = false ;
    tubi.generateCurrentSPCutMesh([], spcutMeshOptions) ;
    
    % Compute the pullback if the cutMesh is ok
    if compute_pullback || ~exist(sprintf(tubi.fullFileBase.im_sp, tt), 'file')
        pbOptions = struct() ;
        tubi.generateCurrentPullbacks([], [], [], pbOptions) ;
    else
        disp('Skipping computation of pullback')
    end
        
end
disp('Done with generating spcutMeshes and cutMeshes')

%% Inspect coordinate system charts using (s,phi) coordinate system ('sp')
options = struct() ;
options.coordSys = 'uv' ;
tubi.coordinateSystemDemo(options)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 3: Further refinement of dynamic meshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Smooth the sphi grid meshes in time ====================================
options = struct() ;
options.overwrite = overwrite ;
options.width = 4 ;  % width of kernel, in #timepoints, to use in smoothing meshes
tubi.smoothDynamicSPhiMeshes(options) ;

%% Plot the time-smoothed meshes
tubi.plotSPCutMeshSmRS(options) ;

% Inspect coordinate system charts using smoothed meshes
options = struct() ;
options.coordSys = 'spsm' ;
tubi.coordinateSystemDemo(options)

%% Redo Pullbacks with time-smoothed meshes ===============================
disp('Create pullback using S,Phi coords with time-averaged Meshes')
for tt = tubi.xp.fileMeta.timePoints
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    tidx = tubi.xp.tIdx(tt);
    
    % Load the data for the current time point ------------------------
    tubi.setTime(tt) ;
    
    % Establish custom Options for MIP --> choose which pullbacks to use
    pbOptions = struct() ;
    pbOptions.numLayers = [0 0] ; % how many onion layers over which to take MIP
    pbOptions.generate_spsm = true ;
    pbOptions.generate_sp = false ;
    pbOptions.overwrite = false ;
    tubi.generateCurrentPullbacks([], [], [], pbOptions) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 4: Computation of tissue deformation, with in-plane and out-of-plane flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TILE/EXTEND SMOOTHED IMAGES IN Y AND RESAVE ============================
% Skip if already done
options = struct() ;
options.coordsys = 'spsm' ;
tubi.doubleCoverPullbackImages(options)
disp('done')

%% PERFORM PIV ON PULLBACK MIPS ===========================================
% % Compute PIV either with built-in phase correlation or in PIVLab
options = struct() ;
options.overwrite = true ;
tubi.measurePIV2d(options) ;

%% Measure velocities =====================================================
disp('Making map from pixel to xyz to compute velocities in 3d for smoothed meshes...')
options = struct() ;
options.show_v3d_on_data = false ;
tubi.measurePIV3d(options) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lagrangian dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pullback pathline time averaging of velocities
options = struct() ;
tubi.timeAverageVelocities(options)
% Velocity plots for pathline time averaging 
options.plot_vxyz = false ;
options.invertImage = true ;
options.averagingStyle = 'Lagrangian'; 
tubi.plotTimeAvgVelocities(options)
% Divergence and Curl (Helmholtz-Hodge) for Lagrangian
options = struct() ;
options.averagingStyle = 'Lagrangian' ;
options.lambda = 0 ;
options.lambda_mesh = 0 ; 
tubi.helmholtzHodge(options) ;

% Compressibility & kinematics for Lagrangian
options = struct() ;
tubi.measureMetricKinematics(options)

%% Metric Kinematics Kymographs & Correlations -- Bandwidth Filtered
options = struct() ;
tubi.plotMetricKinematics(options)

%% Pullback pathlines connecting Lagrangian grids
options = struct() ;
tubi.measurePullbackPathlines(options)

%% Pullback pathline texturepatching (PIV pathline) ======================
disp('Create pullback using pullback pathline coords')
tidx0 = tubi.xp.tIdx(tubi.t0) ;
tp2do = tubi.xp.fileMeta.timePoints(tidx0:10:end) ;
tp2do = [tp2do, setdiff(tubi.xp.fileMeta.timePoints, tp2do)] ;
for tt = tp2do
    disp(['PB Pathline texturepatch: NOW PROCESSING TIME POINT ', num2str(tt)]);
    tidx = tubi.xp.tIdx(tt);
    
    % Load the data for the current time point ------------------------
    tubi.setTime(tt) ;
    
    % Establish custom Options for MIP --> choose which pullbacks to use
    pbOptions = struct() ;
    pbOptions.numLayers = [0 0] ; % how many onion layers over which to take MIP
    pbOptions.generate_spsm = false ;
    pbOptions.generate_sp = false ;
    pbOptions.overwrite = false ;
    pbOptions.generate_pivPathline = true ;
    tubi.generateCurrentPullbacks([], [], [], pbOptions) ;
end

%% Query velocities along pathlines
options = struct() ;
tubi.measurePathlineVelocities(options)
% plot the pathline velocities 
options = struct() ;
options.gridTopology = 'triangulated' ;
tubi.plotPathlineVelocities(options)

% Measure Pathline Kinematics
options = struct() ;
tubi.measurePathlineMetricKinematics(options)

% Plot Pathline Kinematics
options = struct() ;
tubi.plotPathlineMetricKinematics(options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create ricci mesh at t0 to measure Beltrami coefficient in pathlines
options = struct() ;
options.climit = 1 ;
options.coordSys = 'ricci' ;
tubi.measureBeltramiCoefficient(options) ;

%% Strain rate (epsilon = 1/2 (djvi+divj) -vn bij)
options = struct() ;
tubi.measureStrainRate(options) 

%% Plot time-averaged strain rates in 3d on mesh
options = struct() ;
tubi.plotStrainRate3DFiltered(options) 

%% Kymograph strain rates
options = struct() ;
options.clim_trace = 0.05 ;
options.clim_deviatoric = 0.05 ;
tubi.plotStrainRate(options)

% Measure strain rate along pathlines
options = struct() ;
options.overwriteImages = false ;
options.plot_dzdp = false ;
tubi.measurePathlineStrainRate(options)

%% Measure divergence and out-of-plane deformation along pathlines
tubi.measurePathlineMetricKinematics()

% Pathline strain rate plots
options = struct() ;
options.climit = 0.05 ;
options.climitWide = 1.0 ;
tubi.plotPathlineStrainRate(options)

%% Measure strain along pathlines -- note this is from pathlines, not integrating rates
options = struct() ;
options.plot_dzdp = false ;
options.climitInitial = 0.05 ;
options.climitRamp = 0.01 ;
options.climitRatio = 1 ;
tubi.measurePathlineStrain(options)
tubi.plotPathlineStrain(options)



%% PCA decomposition
pcaTypes = {'vnVector', 'v3d', 'vt', 'H2vn', 'vnScalar', 'divv', 'gdot'} ;
% pcaTypes = {'H2vn', 'vnScalar', 'divv', 'gdot'} ;
options = struct('overwrite', true, ...
    'overwriteImages', true) ;
options.pcaTypes = pcaTypes ;
% options.meshStyles = 'sphi' ;
tubi.spaceUnits = [char(181) 'm'] ;
tubi.getPCAoverTime(options)

%% Laplace-Beltrami Spectral (LBS) decomposition
close all; clc;

% lbsTypes = {'vnVector', 'v3d', 'vt', 'H2vn', 'vnScalar', 'divv', 'gdot'} ;
lbsTypes = {'H2vn', 'vnScalar', 'divv', 'gdot'} ;
options = struct('overwrite', true, ...
    'overwriteImages', true) ;
options.lbsTypes = lbsTypes ;
% options.meshStyles = 'sphi' ;
tubi.spaceUnits = [char(181) 'm'] ;
tubi.getLBSoverTime(options)

%% Some other extra functionality: Show accumulated strain as magnitude
tubi.plotPathlineStrain

%% Some other extra functionality: Visualize Lagrangian patch followed by piv Pathlines
opts = struct() ;
opts.timePoints = 96:10:206 ;
opts.demoPatchName = 'demoPatch005' ;
tubi.visualizeSegmentationPatch(opts) 

%% Some other extra functionality: Look at the rms velocity over time
[vrms, timestamps] = tubi.measureRMSvelocityOverTime ;

