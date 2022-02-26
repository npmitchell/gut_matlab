%% Drosophila Pipeline 009
%
% Added display of celllayer in phase finder.
% Added phase clicker option to select ventral furrow on pullbacks.
%
% Requirements
% TRI_NL_BL_handling from github: npmitchell/gut_matlab
% flipThroughStackFindLayer.m
%
%% Setting up ImSAnE
%
% You may need to go into the imsane directory and run setup.m
%

%% Initialize ImSAnE project
%
% We start by clearing the memory and closing all figures.
%

% clear all ;
close all ;
clc ;

%% Add paths
addpath('/tigress/mfl2/Matlab Scripts for Tigress/bfmatlab')
addpath_recurse('/tigress/mfl2/Matlab Scripts for Tigress/gptoolbox')
addpath('/tigress/mfl2/Matlab Scripts for Tigress/TRI_NL_BL_handling')
activefn = matlab.desktop.editor.getActiveFilename;
thispath = fileparts(activefn) ;
addpath(thispath) ;

% activefn = mfilename('fullpath') ;
% thisfn = mfilename() ;
% afn_split = strsplit(activefn, thisfn) ;
% addpath(afn_split{1}) ;

%% Experiment Options
%
% Set options for the experiment. The options are listed in sections based
% on what each one is used for.
%

%
% General options for the experiment.
%
% * 'overwrite'      , the boolean whether to overwrite previous pullbacks
% * 'visualize'      , the boolean whether to show figures
% * 'fn'             , the file name convention for fused data
% * 'zDim'           , the long axis of the embryo
% * 'channels'       , the array containing what channels there are
% * 'primaryChannel' , the primary channels for surface generation
%

overwrite      = false ;  % WARNING! Use with caution
if overwrite
    disp('WARNING: Overwriting data is set to true!')
end
visualize      = false ;
fn             = 'TP0_Ch%d_Ill0_Ang0,45,90,135,180,225,270,315' ;
fn0            = fn ;
zDim           = 1 ;
channels       = 1:3 ;
primaryChannel = 3 ;

% Save original definition of channels (channelstocheck will be removed)
channels_orig = channels ;


%
% Options for image conversion.
%
% * 'file16name' , the name of the converted 16bit files
%

file16name = 'Time_000000_ch%d' ;

%
% Options for the file metadata.
%
% * 'timePointsava'   , the time points available
% * 'stackResolution' , the stack resolution in microns
% * 'swapZT'          , the option to swap the third and fourth dimension
%

timePointsava   = channels ;
stackResolution = [.25 .25 .25] ;
swapZT          = 0 ;

%
% Options for experiment metadata.
%
% * 'channelsUsed' , the channels used
% * 'channelColor' , the mapping from element in channels used to RGB = 123
% * 'description'  , the description of the data
%

channelsUsed = 1 ;
channelColor = 1 ;
description  = 'Stained embryo' ;

%
% Options for detecting the initial point cloud.
%
% * 'guessChannel'             , the integer specifying which channel(s) to
%                                use for detection
% * 'guessSigma'               , the standard deviation of the Gaussian
%                                filter in pixels
% * 'guessssfactor'            , the integer specifying the degree of data
%                                subsampling
% * 'guessnBins'               , number of radial bins to determine the
%                                point cloud in
% * 'guessrmRadialOutliers'    , the threshold to remove radial outliers
% * 'guessrmIntensityOutliers' , the threshold to remove intensity outliers
%

guessChannel             = 1 ;
guessSigma               = 2 ;
guessssfactor            = 8 ;
guessnBins               = 120 ;
guessrmRadialOutliers    = 1.2 ;
guessrmIntensityOutliers = 1.2 ;

%
% Options for visuAutomatically go through all subdirectories to find data to detect and
% generate SOI. (Hierarchy is taken to be two directories deep.)
% labelDirNames{i} '_' dataDirNames{j} '_' 

%
% * 'plot2Ddimval'   , the dimension along which to take cross section
% * 'plot2Dval'      , the slice to use for plotting
% * 'plot3Dssfactor' , the factor to reduce the number of points shown
%

plot2Ddimval   = 'z' ;
plot2Dval      = 350 ;
plot3Dssfactor = 4 ;

%
% Options for coarse fitting the surface.
%
% * 'Rfit'     , the degree of polynomial fit for radius
% * 'Xfit'     , the degree of polynomial fit for x
% * 'Yfit'     , the degree of polynomial fit for y
% * 'efit'     , the degree of polynomial fit for ellipse eccentricity
% * 'phasefit' , the degree of polynom
        %ial fit for ellipse orientation
%

Rfit     = 6 ;
Xfit     = 4 ;
Yfit     = 4 ;
efit     = 3 ;
phasefit = 0 ;

%
% Options for morphsnakes surface detection.
%
% * 'ms_scriptDir'      , the location of the morphsnakes repo
% * 'channel'           , 
% * 'foreGroundChannel' , 
% * 'niter'             , 
% * 'niter0'            , 
% * 'lambda1'           , 
% * 'lambda2'           , 
% * 'nu'                , 
% * 'smoothing'         , 
% * 'post_nu'           , 
% * 'post_smoothing'    , 
% * 'exit_thres'        , 
% * 'pre_nu'            ,
% * 'pre_smoothing'     , 
% * 'radius_guess'      , 
% * 'clip'              , the boolean wheter to clip the intensity of the
%                         raw data
% * 'save_intermediate' , the boolean whether to save intermediate steps
% * 'mlxprogram'        , the location of the mesh script to use
%

ms_scriptDir      = '/tigress/mfl2/code/morphsnakes_wrapper/' ;
channel           = -1 ;
foreGroundChannel = 1 ;
niter             = 40 ;
niter0            = 40 ;
lambda1           = 1 ; %LAMDA 1 AND 2 ARE RESPECTIVELY THE PENALTIES FOR 
                        %ACCEPTING INCORRECT PIXELS AND MISSING CORRECT PIXELS
                        %PIXELS. IF THEY ARE EQUAL, YOU ARE GIVING A SIMILAR
                        %PENALTY TO BOTH.  CHANGE THIS RATIO AND YOU CHANGE 
                        %WHAT THE PROGRAM IS WILLING TO ACCEPT. I.E
                        %INCREASING LAMDA 1 MAKES THE PROGRAM MORE
                        %CONSERVATIVE
lambda2           = 2 ;
nu                = 0.1 ;
smoothing         = 5 ;  %MAKE THIS HIGHER IF THE PLY IS ESCAPING AND FINDING THE EDGE 
post_nu           = -5 ;
post_smoothing    = 1 ;
exit_thres        = 0.0001 ;
pre_nu            = 0 ;
pre_smoothing     = 0 ;
radius_guess      = 30 ;
clip              = 70000.0 ;
save_intermediate = true ;  % true --> saves output to show step-by-step of morphsnakes
mlxprogram        = '/tigress/mfl2/code/meshlab_codes/surface_rm_resample2k_reconstruct_LS3_ssfactor4_tiger.mlx' ;

%
% Options for loading the mesh from morphsnakes.
%
% * 'normal_step_for_mesh' , the number of pixels to nomally evolve
%

normal_step_for_mesh = -27 ;

%
% Options for setting the orientation of the system.
%
% * 'ventralpresent' , the boolean whether the ventral furrow is present
% * 'celllayer'      , the underestimated guess of the thickness of the
%                      embryo cell layer in pixels
% * 'celllayerstep'  , the number of pixels to increase cell layer guess by
%                      each iteration
% * 'min_thres_wo'   , the minimum threshold for x and y eigenvalues to use
%                      if the ventral furrow is not present
%

ventralpresent = false ;
celllayer      = 10 ;
celllayerstep  = 3 ;
min_thres_wo   = 0.01 ;

%
% Options for inspecting morphsnakes fit.
%
% * 'mesh2Ddim'     , the dimension along which to visualize
% * 'mesh2Dzval'    , the value of the slice to visualize
% * 'mesh2Dnoalign' , 
%

mesh2Ddim     = 'z' ;
mesh2Dzval    = 320 ;
mesh2Dnoalign = true ;

%
% Options for the pullback styles.
%
% * 'cylinderone'     , the boolean whether to compute cylinder one
% * 'cylindertwo'     , the boolean whether to compute cylinder two
% * 'cylinderoneprop' , the boolean whether to compute cylinder one proper
% * 'cylindertwoprop' , the boolean whether to compute cylinder two proper

cylinderone     = true ;
cylindertwo     = true ;
cylinderoneprop = false ;
cylindertwoprop = false ;

%
% Options for the onion layers.
%
% * 'onionnLayers'       , the number of onion layers is an odd integer
% * 'onionlayerDistance' , the onion layer distance in pixels
% * 'onionsigma'         , the sigma used for smoothing onion layers
% * 'onionmakeIP'        , the options are SIP, MIP, or both
% * 'onionIPonly'        , 
%

onionnLayers       = 49 ;
onionlayerDistance = 3 ;
onionsigma         = 0 ;
onionmakeIP        = 'both' ;
onionIPonly        = false ;

%
% Options for saving to disc.
%
% * 'soiDir'         , the directory to save the pullbacks to
% * 'imwriteOptions' , the image writing options
% * 'make8bit'       , the boolean whether to save files as 8bit
%

soiDir         = 'cylinder_soi_c%d' ;
imwriteOptions = {'tif'} ;
make8bit       = false ;

%
% Options for clicking the phase on the pullbacks
%
% * 'phaseclicker' , the boolean whether to check pullbacks and click on
%                    phase of ventral furrow
% * 'channelclick' , the channel to use for clicking on ventral furrow
% * 'numcylinders' , an array of which cylinders were generated
%

phaseclicker = true ;
channelclick = 2 ;
numcylinders = [1 2] ;
% Change layerspec_channel is the channel used for layer specification to 
% create MIPs
layerspec_channel = 1 ;

%% Automatic Pullback Generation
%  DEFINE GLBAL PARAMS
% Automatically go through all subdirectories to find data to detect and
% generate SOI. (Hierarchy is taken to be two directories deep.)
%labelDirNames{i} '_' dataDirNames{j} '_' 

% Remove the primary channel from the list of channels
channelstocheck = channels ;
channels = channels(channels ~= primaryChannel) ;

% Automatically detect all the labeltype folders (for ex, 'Eve_Runt_Neuro')
labelfolders = dir ;
labelDirNames = {labelfolders([labelfolders.isdir]).name} ;
labelDirNames = labelDirNames(~ismember(labelDirNames,{'.','..'})) ;

% Save a file name dummy
filedummy = fn ;
%% BIG LOOP START ALL DIF STAINS
for i = 1 : length(labelDirNames)
    
    % Automatically detect all the data folders
    cd(labelDirNames{i})
    datafolders = dir ;
    dataDirNames = {datafolders([datafolders.isdir]).name} ;
    dataDirNames = dataDirNames(~ismember(dataDirNames,{'.','..'})) ;
    
    for j = 1 : length(dataDirNames) % LOOP THROUGH IN. EXPS.
        
        cd(dataDirNames{j})
        dataDir = cd ;
        projectDir = cd ;
        projectIDstr = split(projectDir, '/') ;
        projectIDstr = projectIDstr{end} ;
        
        % Define a file to look for as indicator that this is done
        parametersfile = 'pipeparameters.txt' ;
        
        % disp([labelDirNames{i} '_' dataDirNames{j} '_' ])
        % error('exiting here')
        
        if ~isfile(parametersfile) || overwrite % IF PARAMS FILE IS MISSING THEN LETS DO THAT EXP
            disp(['Making pullbacks for ' projectIDstr])
            for channel_check = channelstocheck
                %% Check that data is 16 bit. If not, convert to 16bit
                fullFileName = [sprintf(filedummy, channel_check) '.tif'] ;
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
                    if chanposition == 1
                        pc2use = 99.9999;
                    elseif chanposition == 2
                        pc2use = 99.99;
                    elseif chanposition == 3
                        pc2use = 99.9999;
                    end
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
                    fullFileName = [sprintf(filedummy, channel_check) '.tif'] ;
                    disp([fullFileName ' has been converted.'])
                    
                    fn = file16name ;
                    swapZT = 1 ;
                else
                    disp('File is 16bit.')
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Create an Experiment class instance called xp
            % Start by creating an experiment object, optionally pass on
            % the project directory (otherwise it will ask). This serves as
            % a frontend for data loading, detection, fitting, etc.
            %
            fn = file16name;
            xp = project.Experiment(projectDir, dataDir) ;
            
            %
            % Set required additional information on the files.
            %
            % We assume one individual image stack for each time point,
            % labeled by time. To be able to load the stack, we need to
            % tell the project where the data is, what convention is
            % assumed for the file names, available time points and the
            % stack resolution. Options to modules in ImSAnE are organised
            % in matlab structures, that is a pair of field name and value
            % are provided for each option.
            %
            
            fileMeta                 = struct() ;
            fileMeta.dataDir         = dataDir ;
            fileMeta.filenameFormat  = [fn '.tif'] ;
            fileMeta.timePoints      = timePointsava ;
            fileMeta.stackResolution = stackResolution ;
            fileMeta.swapZT          = swapZT ;
            
            %
            % Set required additional information on the experiment. A verbal data set
            % description, Jitter correct by translating the sample, which time point to 
            % use for fitting, etc.
            %
            % The following project metadata information is required 
            %
            % * 'channelsUsed'   , the channels used, e.g. [1 3] for RGB
            % * 'channelColor'   , mapping from element in channels used to RGB = 123
            %

            expMeta                  = struct() ;
            expMeta.channelsUsed     = channelsUsed ;
            expMeta.channelColor     = channelColor ;
            expMeta.description      = description ;
            expMeta.dynamicSurface   = 0 ;
            expMeta.jitterCorrection = 0 ;
            expMeta.fitTime          = fileMeta.timePoints(timePointsava == primaryChannel) ; 
            expMeta.detectorType     = 'surfaceDetection.fastCylinderDetector' ;
            expMeta.fitterType       = 'surfaceFitting.spherelikeFitter' ;
            
            %
            % Now set the meta data in the experiment.
            %
            
            xp.setFileMeta(fileMeta) ;
            xp.setExpMeta(expMeta) ;
            
            %
            % initNew() reads the stack size from the first available time
            % point, then initializes the fitter and deterctor and creates
            % fitOptions and detectOptions based on their defaults.
            %
            
            xp.initNew() ;
            
            %% Load the 16bit data for this experiment
            % loadTime sets currentTime, loads the stack, and resets the
            % detector and fitter with the appropriate options for that
            % time. Optionally rescale the stack to unit aspect ratio.
            %

            xp.loadTime(primaryChannel) ;
            xp.rescaleStackToUnitAspect() ;
            
            %% Prepare for ROUGH POLYNOMIAL FIT surface detection
            % Create the struct of detect options.
            %
            
            myDetectOpts = struct('channel', guessChannel, ...
                'sigma', guessSigma, ...
                'ssfactor', guessssfactor, ...
                'nBins', guessnBins, ...
                'rmRadialOutliers', guessrmRadialOutliers, ...
                'rmIntensityOutliers', guessrmIntensityOutliers, ...
                'zDim', zDim) ;
            
            %
            % Set the detect options in the project
            %
            
            xp.setDetectOptions(myDetectOpts) ;
            
            %
            % Calling detectSurface runs the surface detector and creates
            % the detector.pointCloud object.
            %
            xp.detectSurface() ;
            
            %% Visualize the first POINTCLOUD that was detected via INTENSITY THRESHOLDING
            % Inspect the point cloud over a cross section in the data.
            %
            if visualize
                % visualize method 1
                inspectOptions = struct('dimension', plot2Ddimval, ...
                    'value', plot2Dval, ...
                    'pointCloud', 'r', 'thickness', 3) ;
                figure ;
                xp.detector.inspectQuality(inspectOptions, xp.stack) ;
            
                % visualize method 2
                % Plot the points on the surface as a 3D point cloud.
                %
                xp.detector.pointCloud.inspect(plot3Dssfactor) ;
            end
            
            %% FIT POINTCLOUD DETECTED ABOVE TO A POLYNOMIAL
            % IE. THIS IS THE ROUGH POLYINOMIAL PART
            % Fit the surface coarsly to prepare levelset.
            %
            % We fit the surface in a cylindrical basis, with radius,
            % eccentricity, centre of mass, and ellipse orientation as
            % slowly varying polynomials of the z axis.
            %
            % The fitter behaves similar to the detector, but the fitting
            % of the surface is done by calling the function fitSurface()
            % in the experiment class. This function calls
            % fitter.fitSurface(detector.pointCloud), but also stores the
            % resulting fittedParam in the experiment class.
            %
            % Initially we fit the surface with a coarse set of values, to
            % obtain the initial levelset for morphsnakes.
            %

            fitOptions = struct('R', Rfit, ...
                'X', Xfit, ...
                'Y', Yfit, ...
                'e', efit, ...
                'phase', phasefit, ...
                'path', fullfile(projectDir, 'debugOutput')) ;
            xp.setFitOptions(fitOptions) ;
            try
                xp.fitSurface() ;
                surface_fitting_success = true;
            catch 
                surface_fitting_success = false ;
            end
            
            % Only continue if there was a surface that was properly fit
            if surface_fitting_success
                %
                % Now use the polynomial fit of z of ellipses to get inside/outside

                %% SAVE THE POLYNOMIAL FIT AS A LEVELSET GUESS H5 FOR MSLS (FILLED IN LEGO) 
                % If initial guess does not exist, create it.
                %
                initls_h5fn = sprintf('Time_%06d_c3_levelset.h5', xp.currentTime) ;
                if ~exist(initls_h5fn, 'file')
                    % Create initial level set from fit
                    X0 = @(Z) polyval(xp.fittedParam.pX, double(Z), xp.fittedParam.SX, xp.fittedParam.muX) ;
                    Y0 = @(Z) polyval(xp.fittedParam.pY, double(Z), xp.fittedParam.SY, xp.fittedParam.muY) ;
                    rsq = @(Z) polyval(xp.fittedParam.pR, double(Z), xp.fittedParam.SR, xp.fittedParam.muR) ;
                    ee = @(Z) polyval(xp.fittedParam.pe, double(Z), xp.fittedParam.Se, xp.fittedParam.mue) ;
                    % phi = @(Z) polyval(xp.fittedParam.pphase, double(Z), xp.fittedParam.Sphase, xp.fittedParam.muphase) ;
                    % Build the level set from the subsampled data.
                    opts = xp.detector.options ; 
                    data = xp.stack.image.apply{opts.channel} ;
                    data = data(1:opts.ssfactor:end, 1:opts.ssfactor:end, 1:opts.ssfactor:end) ;
                    ls = 0 * data ;
                    zvals = (1:size(data, zDim)) ;
                    if zDim == 1
                        xyinds = [3, 2] ;
                    elseif zDim == 2
                        error('Decide the order of the axes here.')
                    end
                    xarr = (1:size(data, xyinds(1))) * opts.ssfactor ;
                    yarr = (1:size(data, xyinds(2))) * opts.ssfactor ;
                    [xgrid, ygrid] = meshgrid(xarr, yarr) ;
                    for zind = zvals
                        zval = zind * opts.ssfactor ;
                        if mod(zval, 100) == 0
                            disp(['zval = ', num2str(zval)])
                        end
                        cosp = cos(0) ;
                        sinp = sin(0) ;
                        % If phase of the ellipse is necessary uncomment below.
                        %
                        % cosp = cos(phi(zval)) ;
                        % sinp = sin(phi(zval)) ;
                        %
                        % Get the xy values inside the polynomial.
                        rsqz = rsq(zval) ;
                        rbigsqz = rsqz / (1 - ee(zval)^2 )^2 ;
                        X0z = X0(zval) ;
                        Y0z = Y0(zval) ; 
                        dd = (((xgrid - X0z) * cosp + (ygrid - Y0z) * sinp).^2 / rsqz + ((xgrid - X0z) * sinp - (ygrid - Y0z) * cosp).^2 / rbigsqz) ;
                        ls(zind,:,:) = (dd > 0) & (dd < 1) ;
                    end
                    % Write the levelset to disk.
                    h5create(initls_h5fn,'/implicit_levelset', size(ls))
                    h5write(initls_h5fn, '/implicit_levelset', ls)
                else
                    disp('initial implicit_levelset found on disk')
                end

                %% PREPARE FOR MORPHING THE GUESS INTO LEGO AND SMOOTHED SURFACES
                % i.E prepare for filled in lego to go to lego surface &
                % smoothed surface
                % Initialize morphsnakes detectors.
                %

                expMeta.detectorType = 'surfaceDetection.integralDetector_rawdata' ;
                expMeta.fitterType = 'surfaceFitting.cylinderMeshWrapper' ;
                xp.setExpMeta(expMeta) ;
                xp.initNew()

                %
                % If you initialize a new project with xp.initNew(),
                % xp.currentTime returns 1 rather than the time point that we
                % are working with, which is why we ned to use xp.setTime().
                %

                xp.setTime(primaryChannel) ;

                %
                % File naming for morphsnakes.
                %

                ssfactor = xp.detector.options.ssfactor ;
                ofn_ply = 'mesh_apical_ms_' ; % this is the name of lego surface 
                ofn_ls = 'msls_apical_ms_' ; %this is the filled in lego h5
                ofn_smoothply = 'mesh_apical_ms_rep_' ; % this is the smoothed surfae after MeshLab

                %
                % Name the output mesh directory
                %

                msls_exten = ['_nu' strrep(num2str(nu, '%0.2f'), '.', 'p') ] ;
                msls_exten = [msls_exten '_s' num2str(smoothing, '%d') ] ;
                msls_exten = [msls_exten '_pn' num2str(post_nu, '%d') '_ps', num2str(post_smoothing)] ;
                msls_exten = [msls_exten '_l' num2str(lambda1) '_l' num2str(lambda2) ] ;
                msls_exten = strrep(msls_exten, '-', 'm') ;
                if projectDir(end) ~= '/'
                    projectDir = [projectDir '/'] ;
                end
                mslsDir = [projectDir 'msls_output'] ;
                mslsDir = [mslsDir msls_exten '/'] ;
                if ~exist(mslsDir,'dir')
                    mkdir(mslsDir)
                end

                %
                % Set the detection options for morphsnakes.
                %

                msDetectOpts = struct('channel', channel, ...
                    'ssfactor', ssfactor, ...
                    'foreGroundChannel', foreGroundChannel, ...
                    'niter', niter, ...
                    'niter0', niter0, ...
                    'lambda1', lambda1, ...
                    'lambda2', lambda2, ...
                    'nu', nu, ...
                    'smoothing', smoothing, ...
                    'post_nu', post_nu, ...
                    'post_smoothing', post_smoothing, ...
                    'exit_thres', exit_thres, ...
                    'fileName', sprintf(fn, xp.currentTime), ...
                    'mslsDir', mslsDir, ...
                    'ofn_ls', ofn_ls, ...
                    'ofn_ply', ofn_ply,...
                    'ms_scriptDir', ms_scriptDir, ...
                    'timepoint', xp.currentTime, ...
                    'zdim', zDim, ...
                    'pre_nu', pre_nu, ...
                    'pre_smoothing', pre_smoothing, ...
                    'ofn_smoothply', ofn_smoothply, ...
                    'mlxprogram', mlxprogram, ...
                    'init_ls_fn', initls_h5fn, ...
                    'run_full_dataset', 'none', ...
                    'radius_guess', radius_guess, ...
                    'dset_name', 'inputData', ...
                    'clip', clip, ...
                    'save', save_intermediate, ...
                    'center_guess', 'empty_string', ...
                    'plot_mesh3d', false, ...
                    'mask', 'none') ;

                %
                % Set the morphsnakes detect options.
                %

                xp.setDetectOptions(msDetectOpts) ;

                %
                % Move the initial level set to the msls dir.
                %

                if exist(initls_h5fn, 'file')
                    system(['mv ' initls_h5fn ' ' fullfile(mslsDir, initls_h5fn)]) ;
                end

                %
                % Create downsampled dataset as hdf5
                %

                if ~exist([projectDir sprintf(fn, xp.currentTime) '.h5'], 'file')
                    xp.detector.prepareDownsampled(xp.stack) ;
                    disp('Done outputting downsampled data h5 for surface detection.')
                else
                    disp('The h5 file was already output.')
                end

                disp('Downsampling to h5 is done.')

                %% NOW FIND THE LEGO AND SMOOTHED SURFACES FROM THE POLYNOMIAL GUESS H5
                % Detect the surface of the embryo.
                %
                xp.detectSurface()

                %% LOAD THE OUTPUT SMOOTHED SURFACE
                % Read in the mesh file.
                %channelstocheck(1)

                mesh_outfn = [ofn_smoothply, num2str(fileMeta.timePoints(fileMeta.timePoints == xp.currentTime), '%06d'), '.ply'] ;
                outputMesh = fullfile(mslsDir, mesh_outfn) ;
                mesh = read_ply_mod(outputMesh) ;
                msls_axis_order = 'zxyc' ;

                %
                % There is a rotation implicit here by permuting axes
                %

                if strcmp(msls_axis_order, 'zxyc')
                    mesh.v  = mesh.v(:, [2,3,1]) * ssfactor ;
                    % mesh.vn = mesh.vn(:, [2,3,1]) ;
                end

                %% SMOOTH THE NORMALS OF THE LOADED SMOOTHED SURFACE
                % THIS IS IMPORTANT FOR DOING SHIFTS OR ONION LAYERS BECAUSE IT
                % KEEPS THE SURFACE MOVING TOGETHER
                % Compute smoothed normals
                % 
                mesh.f = bfs_orient( mesh.f );
                % Compute vertex normals by weighting the incident triangles by
                % their incident angles to the vertex
                % This uses gptoolbox.
                mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                % Average normals with neighboring normals
                disp('Averaging normals with neighboring normals')
                mesh.vn = average_normals_with_neighboring_vertices(mesh, 0.5) ;

                %
                % Make sure vertex normals are normalized.
                %
                mesh.vn = mesh.vn ./ sqrt( sum( mesh.vn.^2, 2 ) ) ;

                %% MAKE ONE GLOBAL SHIFT ALONG THE NORMALS--> datalayer0
                % Normally evolve vertices.
                % I DEFINED NORMAL STEP FOR MESH WAY UP AT THE TOP BUT NOW IM
                % ACTUALLY USING IT HERE
                mesh.v = mesh.v + normal_step_for_mesh .* mesh.vn ;

                %% Inspect the surface
                % inspectOptions = struct('dimension', 'z','value', 250, 'pointCloud', 'b') ;
                % xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);
                im = xp.stack.image.apply() ;
                im = im{1} ;
                yval = 350 ;
                buff = 5 ;
                fig = figure('visible', 'on') ;
                imshow(squeeze(im(:, yval, :)))
                hold on
                xyslice = mesh.v(mesh.v(:, zDim) < (yval + buff) & mesh.v(:, zDim) > (yval - buff), :) ;
                plot(xyslice(:, 3), xyslice(:, 2), '.')
                %pause
                fn = fullfile(mslsDir, 'meshpreview.png') ;
                disp(['Saving a preview of the mesh fit to ' fn])
                saveas(fig, fn)

                %% Ventral furrow orientation --> this has been removed
                % If ventral furrow not present, uses the eigenvectors
                phase = 0 ;
                min_thres = min_thres_wo ;

                % Determine the orientation of the region of interest.
                %
                % Note that the first argument of pc.determineROI is the
                % margin around the point cloud.
                %
                points    = mesh.v(sum(isnan(mesh.v),2)==0,:);
                pc        = surfaceDetection.PointCloud(points);
                pc.determineROI(1, min_thres) ;
                translation = pc.ROI.translation ;
                rotation    = pc.ROI.rotation ;
                % phase       = phase - acos(rotation(3,1)) ;

                %% PREPARE FOR CYLINDER MESH WRAPPER
                % This is actually an unwrapper cartographically unwrapping
                % the surface.
                % Set the fit options for he surfacefitter
                % (cylinderMeshWrapper)
                %
                fitOptions = struct('chartSeeds', [], ...
                    'diskSeeds', [], ...
                    'phase', phase, ...
                    'transitionWidth', 0, ...
                    'fixAxis', 1, ...
                    'rotation', rotation, ...
                    'translation', translation, ...
                    'fixResolution', 1, ...
                    'resolution', []) ;
                xp.setFitOptions(fitOptions) ;
                xp.fitSurface(mesh) ;


                %% VISUALIZE THE POINTCLOUD OF THE SURFACE AGAIN (NOW IT BELONGS TO THE FITTER,ie cylinderMeshWrapper)
                % Inspect fit and point cloud over a cross section in the data.
                %
                if visualize
                    inspectOptions = struct('dimension', mesh2Ddim, ...
                        'value', mesh2Dzval, ...
                        'pointCloud', 'b', ...
                        'noalign', mesh2Dnoalign) ;
                    figure ;
                    xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack) ;
                end

                %
                % Inspect the whole mesh in three-dimensions. Note that the y
                % axis should be the AP axis of the embryo.
                %
                if visualize
                    figure ;
                    xp.fitter.inspectMesh() ;
                    xlabel('x') ;
                    ylabel('y') ;
                    zlabel('z') ;
                end

                %% PREPARE FOR PULLBACKS
                % Define the pullback styles.
                %
                xp.fitter.setDesiredChart('cylinder1', cylinderone) ;
                xp.fitter.setDesiredChart('cylinder2', cylindertwo) ;
                xp.fitter.setDesiredChart('cylinder1_proper', cylinderoneprop) ;
                xp.fitter.setDesiredChart('cylinder2_proper', cylindertwoprop) ;
                xp.generateSOI() ;
                disp('Done generating SOI for channel.')

                %
                % Set the options for the onion projections.
                %
                onionOpts = struct('nLayers', onionnLayers, ...
                    'layerDistance', onionlayerDistance, ...
                    'sigma', onionsigma, ...
                    'makeIP', onionmakeIP, ...
                    'IPonly', onionIPonly) ;

                %% GENERATE PULLBACKS IN RAM
                % Pass the region of interest and the current time to pull back
                % the stack in the desired charts. This generates the data
                % fields containing the pullback.
                %
                xp.SOI.pullbackStack(xp.stack, [], xp.currentTime, onionOpts) ;
                disp('Done pulling back onions.')

                %
                % Now we extract the data field from the surface of interest at
                % the current time, which is the time of the fit.
                %
                fitOptions    = xp.fitter.fitOptions ;
                data          = xp.SOI.getField('data_MIP') ;
                data          = data(xp.tIdx(xp.currentTime)) ;
                type          = 'cylinder' ;
                patchName     = 'cylinder2_index' ;
                transformName = 'cylinder2' ;
                pb = data.getPatch(patchName).getTransform(transformName).apply{1} ;
                if visualize
                    figure, imshow(pb',[],'InitialMagnification',66) ;
                    figure ; imagesc(pb') ; colorbar ;
                end

                %% SAVE PULLBACKS TO DISK
                % Save the pullbacks to disc.
                %
                saveDir = fullfile(projectDir, sprintf(soiDir, xp.currentTime)) ;
                options = struct('dir', saveDir, ...
                    'imwriteOptions', {imwriteOptions}, ...
                    'make8bit', make8bit) ;
                xp.SOI.save(options)

                %% SAVE A TXT FILE SAYING WHAT WE'VE ACCOMPLISHED
                % Save a plain text file containing the options used in the
                % experiment.
                %

                optionstring = ['guessChannel = ' num2str(guessChannel), ...
                    '\nguessSigma = ' num2str(guessSigma), ...
                    '\nguessssfactor = ' num2str(guessssfactor), ...
                    '\nguessnBins = ' num2str(guessnBins), ...
                    '\nguessrmRadialOutliers = ' num2str(guessrmRadialOutliers), ...
                    '\nguessrmIntensityOutliers = ' num2str(guessrmIntensityOutliers), ...
                    '\nzDim = ' num2str(zDim), ...
                    '\nRfit = ' num2str(Rfit), ...
                    '\nXfit = ' num2str(Xfit), ...
                    '\nYfit = ' num2str(Yfit), ...
                    '\nefit = ' num2str(efit), ...
                    '\nphasefit = ' num2str(phasefit), ...
                    '\nchannel = ' num2str(channel), ...
                    '\nforeGroundChannel = ' num2str(foreGroundChannel), ...
                    '\nniter = ' num2str(niter), ...
                    '\nniter0 = ' num2str(niter0), ...
                    '\nmlx_program = ' mlxprogram, ...
                    '\nlambda1 = ' num2str(lambda1), ...
                    '\nlambda2 = ' num2str(lambda2), ...
                    '\nexit_thres = ' num2str(exit_thres), ...
                    '\nsmoothing = ' num2str(smoothing), ...
                    '\nnu = ' num2str(nu), ...
                    '\npre_nu = ' num2str(pre_nu), ...
                    '\npre_smoothing = ' num2str(pre_smoothing), ...
                    '\npost_nu = ' num2str(post_nu), ...
                    '\npost_smoothing = ' num2str(post_smoothing), ...
                    '\nradius_guess = ' num2str(radius_guess), ...
                    '\nclip = ' num2str(clip), ...
                    '\ncelllayer = ' num2str(celllayer), ...
                    '\ncelllayerstep = ' num2str(celllayerstep), ...
                    '\nonionnLayers = ' num2str(onionnLayers), ...
                    '\nonionlayerDistance = ' num2str(onionlayerDistance), ...
                    '\nonionsigma = ' num2str(onionsigma), ...
                    '\nnormal_step_for_mesh = ', num2str(normal_step_for_mesh)] ;

                fparam = fopen(parametersfile, 'wt') ;
                fprintf(fparam, optionstring) ;
                fclose(fparam) ;

                %% NOW DO IT AGAIN FOR OTHER CHANNELS USING THE SURFACE FOUND ALREADY
                % 
                for b = channels

                    disp(['Processing channel ' num2str(b)]) ;

                    %
                    % Load in channel b and rescale to unit aspect ratio.
                    %

                    xp.loadTime(b) ;
                    xp.rescaleStackToUnitAspect() ;

                    %
                    % Use the same mesh as detected before to generate
                    % pullbacks for the other channels.
                    %

                    points = mesh.v ;
                    pc = surfaceDetection.PointCloud(points) ;
                    pc.determineROI(1, min_thres) ;
                    rotation = pc.ROI.rotation ;
                    translation = pc.ROI.translation ;

                    %
                    % Set the fit options for surfacefitter
                    %

                    fitOptions = struct('chartSeeds', [], ...
                        'diskSeeds', [], ...
                        'phase', phase, ...
                        'transitionWidth', 0, ...
                        'fixAxis', 1, ...
                        'rotation', rotation, ...
                        'translation', translation, ...
                        'fixResolution', 1, ...
                        'resolution', []) ;
                    xp.setFitOptions(fitOptions) ;
                    xp.fitSurface(mesh) ;

                    xp.fitter.setDesiredChart('cylinder1', cylinderone) ;
                    xp.fitter.setDesiredChart('cylinder2', cylindertwo) ;
                    xp.fitter.setDesiredChart('cylinder1_proper', cylinderoneprop) ;
                    xp.fitter.setDesiredChart('cylinder2_proper', cylindertwoprop) ;
                    xp.generateSOI() ;
                    disp('Done generating SOI.')

                    %
                    % Set the options for the onion projections.
                    %

                    onionOpts = struct('nLayers', onionnLayers, ...
                        'layerDistance', onionlayerDistance, ...
                        'sigma', onionsigma, ...
                        'makeIP', onionmakeIP, ...
                        'IPonly', onionIPonly) ;

                    %
                    % Pass the region of interest and the current time to pull back
                    % the stack in the desired charts. This generates the data
                    % fields containing the pullback.
                    %

                    xp.SOI.pullbackStack(xp.stack, [], xp.currentTime, onionOpts) ;
                    disp('Done pulling back onions')

                    %
                    % Save to disc
                    %

                    saveDir = fullfile(projectDir, sprintf(soiDir, xp.currentTime)) ;
                    options = struct('dir', saveDir, ...
                        'imwriteOptions', {imwriteOptions}, ...
                        'make8bit', make8bit) ;
                    xp.SOI.save(options)
                end
            else
            end
        else
            disp(['Data in ' projectDir ' has already been processed.'])
        end
        fn = filedummy ;
        cd ..
    end
    cd ..
end

