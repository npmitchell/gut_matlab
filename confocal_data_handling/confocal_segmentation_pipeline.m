%% Segment out a Monge-form "slab" [z0(x,y), z1(x,y)] from confocal data 
% NPMitchell 2020
%
% This is a pipeline to segment confocal data & take MIPs

% temporary path def
% cd /mnt/data/confocal_data/gut/2020/mef2GAL4klarUASsqhGFPUASCAAXmCh/202010191649_mef2GAL4klarUASsqhGFPUASCAAXmCh_40x1p6x_240spf_0p4um_5to10pc3p5to7pc_oilImmersol/

dataDir = '/mnt/data/optogenetics_confocal/' ;
dataDir = [dataDir 'antpGAL4/huygens_deconvolution_20210325/'] ;
dataDir = [dataDir '202103181730_antpG4OCRLGap43mCh_40x1p6x_5mpf_4pc3pc_to_12pc9pc_600ns_lav3_DC/'] ;
cd(dataDir)
gutDir = '/mnt/data/code/gut_matlab/' ;
addpath(fullfile(gutDir, 'addpath_recurse'))
addpath_recurse(gutDir)

%% INITIALIZE ImSAnE PROJECT ==============================================
%
% We start by clearing the memory and closing all figures
clear; close all; clc;

% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.

dataDir    =  cd; 
projectDir = dataDir ;

% A filename base template - to be used throughout this script
% the 32 bit fn
fn = '' ;
% the 16 bit fn
file16name = 'antpOCRLgap43_T%03d' ;     
um2pix = 5.6284 ;
resolution = [1/um2pix, 1/um2pix, 1.4] ;
timepoints = 1:38 ;

%% Join data into stacks
fns0 = dir('./splitChannels/*ch00.tif') ;
fns1 = dir('./splitChannels/*ch01.tif') ;
for tidx = 1:length(fns0)
    fn0 = fullfile(fns0(tidx).folder, fns0(tidx).name) ;
    fn1 = fullfile(fns1(tidx).folder, fns1(tidx).name) ;
    fn2 = sprintf([file16name '.tif'], tidx) ;
    
    if ~exist(fn2, 'file')

        im0 = readTiff4D(fn0, 1) ;
        im1 = readTiff4D(fn1, 1) ;
        im = cat(4, im0, im1) ;

        writeTiff5D(permute(im, [1, 2, 4, 3]), fn2)

        % % check it
        % for qq = 1:size(im0, 3)
        %   imagesc(squeeze(im0(:, :, qq))) ;
        %   pause(0.01)
        % end
    end
end

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
% The following file metadata information is requireL:
% * 'directory'         , the project directory (full path)
% * 'dataDir'           , the data directory (full path)
% * 'filenameFormat'    , fprintf type format spec of file name
% * 'timePoints'        , list of times available stored as a vector
% * 'stackResolution'   , stack resolution in microns, e.g. [0.25 0.25 1]
%
% The following file metadata information is optional:
%
% * 'imageSpace'        , bit depth of image, such as uint16 etc., defined
%                         in Stack class
% * 'stackSize'         , size of stack in pixels per dimension 
%                         [xSize ySize zSize]
% * 'swapZT'            , set=1 if time is 3rd dimension and z is 4th
              

fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fn, '.tif'];
fileMeta.nChannels          = 2;
fileMeta.timePoints         = timepoints;
fileMeta.stackResolution    = resolution ; % the px resolution (found in the .lif; 4 dec places)
fileMeta.swapZT             = 0;

% Set required additional information on the experiment. A verbal data set
% description, Jitter correct by translating  the sample, which time point
% to use for fitting, etc.
%
% The following project metadata information is requireL:
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
expMeta.channelsUsed        = [1,2];
expMeta.channelColor        = 1;
expMeta.description         = 'Drosophila gut';
expMeta.dynamicSurface      = 0;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.fitTime             = fileMeta.timePoints(first_tp);
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';

% %% 32 to 16 BIT CONVERSION
% % Check that data is 16 bit. If not, convert to 16bit
% pc2use = 99.9;
% scalemethod = 'user-defined' ; % use either 'user-defined' or 'prctile'
% scalemax = 0.05 ;
% 
% for channel_check = 1: fileMeta.timePoints(end) % last tp to convert to 16 bit
%     fullFileName = [sprintf(fn, channel_check) '.tiff'] ;
%     info = imfinfo(fullFileName) ;
%     full16fn = [sprintf(file16name, channel_check) '.tiff'] ;
%     bitDepth = info.BitDepth ;
% 
%     if (bitDepth == 32) && ~isfile(full16fn)
% 
%         disp([fullFileName ' is not 16bit, converting...'])
% 
%         % Note that imread only loads a single frame
%         % A = imread(fullFileName) ;
%         % scalemin = double(min(A(:))) ;
%         % scalemax = double(max(A(:))) ;
%         disp('')
%         disp('Reading 32 bit file to convert...')
%         A = readSingleTiff(fullFileName) ;
%         tmpA = A(:) ;
%         disp('')
%         disp('Computing scalemin, scalemax')
% 
%         % Optional step here to figure out what the cutoff
%         % intensity should be
%         % tmpA_no_ouliers = tmpA(tmpA < pcntile(tmpA, 99)) ;
%         % thisstd = std(tmpA_no_ouliers) ;
%         % check it using histogram(tmpA)
%         thismedian = median(tmpA) ;
% 
%         %goodmedian = 2559.00;
%         %worstmedian = 420.00;
%         %range2correct = goodmedian - worstmedian ;
%         %normal_pc2use = 99.9999 ;
%         %worstcase_pc2use = 99.99 ;
%         %diffpc = normal_pc2use - worstcase_pc2use ;
%         %pc2use = normal_pc2use + diffpc * (thismedian - goodmedian) / range2correct ;
%         %pc2use = max(worstcase_pc2use, pc2use) ;
%         %pc2use = min(normal_pc2use, pc2use) ;
%         chanpositionstart = strfind(fullFileName,'Ch');
%         chanposition = fullFileName(chanpositionstart+2);
%         chanposition = str2num(chanposition);
%         disp('determining prctile')
%         if strcmp('scalemethod', 'prctile')
%             scalemax = double(prctile(tmpA, pc2use)) ;
%         else
%             disp('using user supplied scalemax')
%         end
%             
%         scalemin = double(min(tmpA(tmpA > 0))) ;
%         disp(['scalemax = ', num2str(scalemax)])
%         disp(['scalemin = ', num2str(scalemin)])
% 
%         disp('Showing slice of scaled image')
%         % Note to self to check scale:
%         close all
%         imagesc(squeeze(A(:, 300, :)))
%         title('Checking scale of image')
%         waitfor(gcf)
%         histogram(A(:))
%         title('Histogram of values')
%         waitfor(gcf)
% 
%         % data = readSingleTiff(fullFileName);
%         im2 = mat2gray(A,[scalemin scalemax]);
%         im2 = uint16(2^16*im2);
%         imSize = size(im2);
%         
%         % Check scale:
%         imagesc(squeeze(im2(:, 300, :)))
%         title('Checking scale of image. Close image to continue')
%         axis equal
%         colorbar
%         waitfor(gcf)
% 
%         % Save the 16 bit image
%         disp(['Saving 16bit volume to ' full16fn])
%         imwrite(im2(:,:,1),full16fn,'tiff','Compression','none');
%         for z = 2 : imSize(3)
%             imwrite(im2(:,:,z),full16fn,'tiff','Compression','none','WriteMode','append');
%         end
%         disp('done saving 16bit volume')
% 
%         fn = file16name ;
% 
%     elseif isfile(full16fn)
%         % the data is already 16 bit, so we're good
%         fullFileName = [sprintf(fn, channel_check) '.tiff'] ;
%         disp([fullFileName ' has been converted.'])
% 
%         fn = file16name ;
%     else
%         disp('File is 16bit.')
%     end
% end

%% INSTANTIATE EXPERIMENT CLASS
% Now set the meta data in the experiment.
fileMeta.filenameFormat = [ file16name '.tif' ] ;
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();

%%%%%%%%%%%%%%%%%%%%%%%%%
% Other options
%%%%%%%%%%%%%%%%%%%%%%%%%
mlxprogram = 'surface_rm_resample20k_reconstruct_LS3_1p2pc_ssfactor4.mlx';
msls_axis_order = 'xyzc';
% Mesh marching options
normal_step = 10;

% Define the surface detection parameters
channel = 2;
foreGroundChannel = 2;
ssfactor = 2;
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
            'exit_thres', 1e-6, ... % convergence thresholL: maximum difference between subsequent level sets upon which to exit algorithm ('close enough')
            'foreGroundChannel',foreGroundChannel, ... % the index of the first dimension of the 4d input data (if 4d)
            'fileName', sprintf( fn, xp.currentTime ), ... % the filename of h5 to train on
            'mslsDir', mslsDir, ...  % the directory for all output data/images
            'ofn_ls', ofn_ls, ...  % the output filename for level sets
            'ofn_ply', ofn_ply, ... % the output filename for PLY files
            'ms_scriptDir', ms_scriptDir, ... % the directory containing run_morphsnakes.py
            'timepoint', 0, ... % which timepoint in the data to consider
            'zdim', 3, ... % Which dimension is the z dimension
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
            'physicalaxisorder', 'yxzc', ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc
            'preilastikaxisorder', 'xyzc', ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
            'ilastikaxisorder', 'xyzc', ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
            'include_boundary_faces', true, ... % keep faces along the boundaries of the data volume if true
            'smooth_with_matlab', -1); % if <0, use meshlab. If >0, smooth the mesh after marching cubes mesh creation using matlab instead of mlxprogram, with diffusion parameter lambda = this value. If =0, no smoothing.
            
        
% Set detect options ------------------------------------------------------
xp.setDetectOptions( detectOptions );

% clear msls_exten imwriteOptions saveDir
% clear channel foreGroundChannel
% clear niter niter0 lambda1 lambda2
% clear exit_thres smoothing nu
% clear post_nu post_smoothing



%% Adjust LUT to be approx constant over time & Isotropic resolution
% Use Fiji BioFormats Exporter to export each channel and each timepoint to
% a separate TIFF
adjustlow_pctile = 0;
adjusthigh_pctile = 99.5 ;
for tt = xp.fileMeta.timePoints
    name_out = sprintf('./combined_T%02d.tif', tt) ;
    if ~exist(name_out, 'file')
        xp.loadTime(tt)
        xp.rescaleStackToUnitAspect();
        IV = xp.stack.image.apply() ;
        % Save a TIFF of the 2color data
        for ii = 1:length(IV)
            IVii = IV{ii} ;
            vlo = double(prctile( IVii(:) , adjustlow_pctile )) / double(max(IVii(:))) ;
            vhi = double(prctile( IVii(:) , adjusthigh_pctile )) / double(max(IVii(:))) ;
            disp(['--> ', num2str(vlo), ', ', num2str(vhi), ...
                ' for ', num2str(adjustlow), '/', num2str(adjusthigh)])
            IV{ii} = imadjustn(IVii, [double(vlo); double(vhi)]) ;
        end
        IV1 = IV{1} ;
        IV2 = IV{2} ;
        im = uint16(zeros(size(IV1, 1), size(IV1,2), 2, size(IV1,3), 1)) ;
        im(:, :, 1, :, 1) = reshape(uint16(IV1), [size(IV1, 1), size(IV1, 2), 1, size(IV1, 3), 1]) ;
        im(:, :, 2, :, 1) = reshape(uint16(IV2), [size(IV1, 1), size(IV1, 2), 1, size(IV1, 3), 1]) ;
        dimensionOrder = 'XYZCT' ;   % must be in {'XYZCT', 'XYZTC', 'XYCTZ', 'XYCZT', 'XYTCZ', 'XYTZC'}           
        % bfsave(im, name_out, 'dimensionOrder', dimensionOrder, 'BigTiff', false)

        writeTiff5D(im, name_out)
    end
end

%% INSTANTIATE EXPERIMENT CLASS AGAIN WITH COMBINED STACKS
% Now set the meta data in the experiment.
fileMeta.filenameFormat = 'combined_T%02d.tif' ;
fileMeta.stackResolution = min(fileMeta.stackResolution) * [1,1,1] ;
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();


%% CREATE THE SUBSAMPLED H5 FILE FOR INPUT TO ILASTIK =====================
% skip if already done
for t = xp.fileMeta.timePoints(1:end)
    if ~exist(fullfile(projectDir, [(sprintf(file16name, t)) '.h5']), 'file')
        disp(['Did not find file: ', fullfile(projectDir, [sprintf(file16name, t) '.h5'])])
        
        % Only load and rescale if multiple timepoints/channels
        % if length(xp.fileMeta.timePoints) > 1 && t > xp.fileMeta.timePoints(1)
        xp.loadTime(t);
        xp.rescaleStackToUnitAspect();
        % end
        
        % make a copy of the detectOptions and change the fileName
        detectOpts2 = detectOptions ;
        detectOpts2.fileName = sprintf( file16name, xp.currentTime ) ;
        xp.setDetectOptions( detectOpts2 );
        xp.detector.prepareIlastik(xp.stack);
        disp(['done outputting downsampled data h5: tp=' num2str(t) ' for surface detection'])
    else
        disp(['h5 ' num2str(t) ' was already output, skipping...'])
    end
end
disp('Open with ilastik if not already done')

%% TRAIN DATA IN ILASTIK TO IDENTIFY SURFACE ==========================
% % open ilastik, train until probabilities and uncertainty are satisfactory
% 

%% MASK TIFFs with subsampled Probabilities
probChannel = 1 ;  % which channel to use as mask
clip = 0.0 ;
outdirs = {'./mips/foreground', './mips/background', './mips/all_channels'} ;
for qq = 1:length(outdirs)
    if ~exist(outdirs{qq}, 'dir')
        mkdir(outdirs{qq})
    end
end
for t = xp.fileMeta.timePoints
    % TIFF filename output
    outfnF = sprintf([ file16name '_foreground.tif'], t) ;
    outfnB = sprintf([ file16name '_background.tif'], t) ;
    tif_exist = exist(outfnF, 'file') && exist(outfnB, 'file') ;
     
    if ~tif_exist
        % Load the data first to set time
        xp.loadTime(t);
        xp.rescaleStackToUnitAspect();
        IV = xp.stack.image.apply() ;

        % Probabilities filename input
        probfn = sprintf(file16name, xp.currentTime) ;

        probField = h5read([sprintf(probfn, tt), '_Probabilities.h5'], '/exported_data') ; % outputed prob from Ilastik corresponding to this timepoint
        probField = squeeze(probField(:,:,:,probChannel) ) ;
        % Full scale resolution now
        probFull = superkron(probField, ones(ssfactor, ssfactor, ssfactor)) ;

        % Preallocate the masked data arrays
        datmF = uint16(zeros(size(IV{1},1), size(IV{1}, 2), length(IV), size(IV{1}, 3), 1)) ;
        datmB = datmF ;

        for ch = 1:length(IV)
            dat = IV{ch} ;
            if any(size(dat) ~= size(probFull))
                % resize probFull a bit -- missing row or column
                smallDim = find(size(probFull) < size(dat)) ;
                % Prob array is too small in some dimensions
                if ~isempty(smallDim)
                    for dim = smallDim
                        firstID = size(probFull, dim) ;
                        lastID = size(dat, dim) ;
                        % Build larger array with more rows/cols/slices
                        for ind = firstID:lastID
                            if dim == 1
                                probFull(ind, :, :) = probFull(lastID, :, :) ;
                            elseif dim == 2
                                probFull(:, ind, :) = probFull(:, lastID, :) ;
                            elseif dim == 3
                                probFull(:, :, ind) = probFull(:, :, lastID) ;
                            end
                        end
                    end
                end

                % Trim probabilities if too big
                bigDim = find(size(probFull) > size(dat)) ;
                % Prob array is too big in some dimensions
                if ~isempty(bigDim)
                    for dim = bigDim
                        if dim == 1
                            probFull = probFull(1:size(dat, dim), :, :) ;
                        elseif dim == 2
                            probFull = probFull(:, 1:size(dat, dim), :) ;
                        elseif dim == 3
                            probFull = probFull(:, :, 1:size(dat, dim)) ;
                        else
                            error('What dimension is this?')
                        end
                    end
                end
            end

            % Now, clip probabilities field (do this after resizing so that
            % can resize just once instead of resizing both F & B)
            probF = probFull ;
            probF(probF < clip) = 0 ;
            probF = single(probF - clip) / single(1 - clip) ;

            probB = 1 - probFull ;
            probB(probB < clip) = 0 ;
            probB = single(probB - clip) / single(1 - clip) ;

            datmF(:, :, ch, :, 1) = uint16(probF .* single(dat)) ;
            datmB(:, :, ch, :, 1) = uint16(probB .* single(dat)) ;

        end

        % Save the masked data as TIFF
        writeTiff5D(uint16(datmF), outfnF)
        writeTiff5D(uint16(datmB), outfnB)
    end    
end

%% Get MIPS normalization curves
for tidx = 1:length(xp.fileMeta.timePoints)
    t = xp.fileMeta.timePoints(tidx) ;
    disp(['t = ' num2str(t)])
    
    % mip output fn
    outMIPF = fullfile('./mips', 'foreground', sprintf(['mip_' file16name '_foreground.mat'], t)) ;
    outMIPB = fullfile('./mips', 'background', sprintf(['mip_' file16name '_background.mat'], t)) ;
    
    if ~exist(outMIPF, 'file') || ~exist(outMIPB, 'file') 
        % TIFF filename output
        outfnF = sprintf([ file16name '_foreground.tif'], t) ;
        outfnB = sprintf([ file16name '_background.tif'], t) ;
        % Load the tiffs
        datmF = readTiff4D(outfnF, 2) ;
        datmB = readTiff4D(outfnB, 2) ;

        % Make mip of this channel with designated colors
        for ch = 1:length(datmF)
            imF{ch} = squeeze(max(squeeze(datmF{ch}), [], 3)) ;
            imB{ch} = squeeze(max(squeeze(datmB{ch}), [], 3)) ;
        end
        
        % Save the mip to mat file
        save(outMIPF, 'imF')
        save(outMIPB, 'imB')
    end
end
    
%% Load/compute statistics of each mip
for tidx = 1:length(xp.fileMeta.timePoints)
    t = xp.fileMeta.timePoints(tidx) ;
    
    % mip output fn
    mipfnF = fullfile('./mips', 'foreground_mat', sprintf(['mip_' file16name '_foreground.mat'], t)) ;
    mipfnB = fullfile('./mips', 'background_mat', sprintf(['mip_' file16name '_background.mat'], t)) ;
    
    % Load the tiff
    load(mipfnF, 'imF')
    load(mipfnB, 'imB')
    if tidx == 1
        mF = zeros(length(xp.fileMeta.timePoints), 2) ;
        mB = zeros(length(xp.fileMeta.timePoints), 2) ;
        sF = zeros(length(xp.fileMeta.timePoints), 2) ;
        sB = zeros(length(xp.fileMeta.timePoints), 2) ;
    end
    
    % Take stats
    for ch = 1:length(imF)
        mF(tidx, ch) = mean(double(imF{ch}(:))) ;
        sF(tidx, ch) = std(double(imF{ch}(:)))  ;
        mB(tidx, ch) = mean(double(imB{ch}(:))) ;
        sB(tidx, ch) = std(double(imB{ch}(:)))  ;
    end
end
% Save as plot -- FOREGROUND
clf
colors = [0,1,1; ...  % cyan
          1,0,0];     % red
for ch = 1:size(mF, 2)
    lineProps = {'-','color', colors(ch, :)} ;
    hs{ch} = shadedErrorBar(xp.fileMeta.timePoints, mF(:,ch), sF(:,ch), 'lineProps', lineProps) ;
end
xlabel(['time ' timeUnits])
ylabel('Mean intensity')
title('Midgut signal')
legend({'sqh', 'membrane'})
saveas(gcf, './mips/means_stdevs_foreground.png')

% BACKGROUND
clf
for ch = 1:size(mB, 2)
    lineProps = {'-','color', colors(ch, :)} ;
    hs{ch} = shadedErrorBar(xp.fileMeta.timePoints, mB(:,ch), sB(:,ch), 'lineProps', lineProps) ;
end
xlabel(['time ' timeUnits])
ylabel('Mean intensity')
title('Background signal')
legend({'sqh', 'membrane'})
saveas(gcf, './mips/means_stdevs_background.png')

% Save means and stdevs
save('./mips/means_stdevs.mat', 'mF', 'sF', 'mB', 'sB')

%% MAKE MIPS
invert = false
for tidx = 1:length(xp.fileMeta.timePoints)
    t = xp.fileMeta.timePoints(tidx) ;
    disp(['t = ', num2str(t)])
    
    % TIFF filename output
    outfnF = sprintf([ file16name '_foreground.tif'], t) ;
    outfnB = sprintf([ file16name '_background.tif'], t) ;
    mipMatF = fullfile('./mips', 'foreground_mat', sprintf(['mip_' file16name '_foreground.mat'], t)) ;
    mipMatB = fullfile('./mips', 'background_mat', sprintf(['mip_' file16name '_background.mat'], t)) ;
    outMIPF = fullfile('./mips', 'foreground', sprintf(['mip_' file16name '_foreground.tif'], t)) ;
    outMIPB = fullfile('./mips', 'background', sprintf(['mip_' file16name '_background.tif'], t)) ;
    totalfn = fullfile('./mips', 'all_channels', sprintf(['mip_all_' file16name '.tif'], t)) ;
    
    mip_exist = exist(outMIPF, 'file') && exist(outMIPB, 'file') ;
    if ~mip_exist || ~exist(totalfn, 'file') || true
               
        % Load the tiffs
        load(mipMatF, 'imF') ;
        load(mipMatB, 'imB') ;
        
        % Make mip of this channel with designated colors
        % imF{ch} = squeeze(max(squeeze(datmF(:, :, ch, :, 1)), [], 3)) ;
        % imB{ch} = squeeze(max(squeeze(datmB(:, :, ch, :, 1)), [], 3)) ;

        for ch = 1:length(imF)
            imF{ch} = double(imF{ch}) ./ double(mF(tidx,ch) + 2 * sF(tidx,ch)) ;
            imB{ch} = double(imB{ch}) ./ double(mF(tidx,ch) + 2 * sF(tidx,ch)) ;
%             imF{ch} = double(imF{ch}) ./ double(max(imF{ch}(:))) ;
%             imB{ch} = double(imB{ch}) ./ double(max(imB{ch}(:))) ;
        end
        
        % Write MIP as RGB with falsecolors
        imoutF = uint8(zeros(size(imF{1}, 1), size(imF{1}, 2), 3)) ;
        imoutB = uint8(zeros(size(imF{1}, 1), size(imF{1}, 2), 3)) ;
        maxI = 255 ;
        for ch = 1:length(IV)
            if invert
                disp('inverting')
                colors = [1, 0, 0; ...  % cyan
                          0, 1, 1];     % red
                imoutF(:, :, 1) = imoutF(:, :, 1) + uint8(maxI * (1-imF{ch}) * colors(ch, 1)) ;
                imoutF(:, :, 2) = imoutF(:, :, 2) + uint8(maxI * (1-imF{ch}) * colors(ch, 2)) ;
                imoutF(:, :, 3) = imoutF(:, :, 3) + uint8(maxI * (1-imF{ch}) * colors(ch, 3)) ;
                imoutB(:, :, 1) = imoutB(:, :, 1) + uint8(maxI * (1-imB{ch}) * colors(ch, 1)) ;
                imoutB(:, :, 2) = imoutB(:, :, 2) + uint8(maxI * (1-imB{ch}) * colors(ch, 2)) ;
                imoutB(:, :, 3) = imoutB(:, :, 3) + uint8(maxI * (1-imB{ch}) * colors(ch, 3)) ;
            else
                colors = [0, 1, 1; ...  % cyan
                          1, 0, 0];     % red
                imoutF(:, :, 1) = imoutF(:, :, 1) + uint8(maxI * imF{ch} * colors(ch, 1)) ;
                imoutF(:, :, 2) = imoutF(:, :, 2) + uint8(maxI * imF{ch} * colors(ch, 2)) ;
                imoutF(:, :, 3) = imoutF(:, :, 3) + uint8(maxI * imF{ch} * colors(ch, 3)) ;
                imoutB(:, :, 1) = imoutB(:, :, 1) + uint8(maxI * imB{ch} * colors(ch, 1)) ;
                imoutB(:, :, 2) = imoutB(:, :, 2) + uint8(maxI * imB{ch} * colors(ch, 2)) ;
                imoutB(:, :, 3) = imoutB(:, :, 3) + uint8(maxI * imB{ch} * colors(ch, 3)) ;
            end
        end
        imwrite(imoutF, outMIPF)
        imwrite(imoutB, outMIPB)
        
        % Write total
        colors = [0.8500, 0.3250, 0.0980 ; % red
                 0.9290, 0.6940, 0.1250 ; % yellow
                 0.0000, 0.4470, 0.7410 ; % blue
                 0.4940, 0.1840, 0.5560 ] ; % purple
        imout = uint8(zeros(size(imF{1}, 1), size(imF{1}, 2), 3)) ;
        maxI = 255 ;
        for ch = 1:length(IV)
            imout(:, :, 1) = imout(:, :, 1) + uint8(maxI * imF{ch} * colors(ch, 1)) ;
            imout(:, :, 2) = imout(:, :, 2) + uint8(maxI * imF{ch} * colors(ch, 2)) ;
            imout(:, :, 3) = imout(:, :, 3) + uint8(maxI * imF{ch} * colors(ch, 3)) ;
            imout(:, :, 1) = imout(:, :, 1) + uint8(maxI * imB{ch} * colors(ch+2, 1)) ;
            imout(:, :, 2) = imout(:, :, 2) + uint8(maxI * imB{ch} * colors(ch+2, 2)) ;
            imout(:, :, 3) = imout(:, :, 3) + uint8(maxI * imB{ch} * colors(ch+2, 3)) ;
        end
        imwrite(imout, totalfn)
    end
end

%% INSTANTIATE EXPERIMENT CLASS AGAIN WITH MASKED STACKS
% Now set the meta data in the experiment.
fileMeta.filenameFormat = [file16name '_foreground.tif'] ;
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();

%% Take ratio of two channels
timeInterval = 4 ;
timeUnits = 'min' ;
chRGB = [3, 1] ;
clip = 2 ;
clipDir = fullfile('./diffs', ['clip' num2str(clip)]) ;
if ~exist(clipDir, 'dir')
    mkdir(clipDir)
end
for t = xp.fileMeta.timePoints
    
    % TIFF filename output
    dFfn0 = fullfile(clipDir, sprintf(['mip_diff_' file16name '.tif'], t)) ;
    if ~exist(dFfn0, 'file')
        % Load the data for this timepoint
        % xp.loadTime(t);
        % xp.rescaleStackToUnitAspect();
        % IV = xp.stack.image.apply() ;
        % 
        % mean1 = double(mean(IV{1}(IV{1} > 0))) ;
        % mean2 = double(mean(IV{2}(IV{2} > 0))) ;
        % dF = double(IV{1})/mean1 - double(IV{2})/mean2 ;
        % dFmip = mean(dF, 3) ;
        
        % OPTION 2: Load the mips and subtract
        mipfn = fullfile('./mips', sprintf(['mip_' file16name '_foreground.tif'], t)) ;
        im = imread(mipfn) ;
        im1 = squeeze(im(:, :, chRGB(1))) ;
        im2 = squeeze(im(:, :, chRGB(2))) ;
        mean1 = double(mean(im1(im1 > 0))) ;
        mean2 = double(mean(im2(im2 > 0))) ;
        im1 = double(im1) / mean1 ;
        im2 = double(im2) / mean2 ;
        dFmip = im1 - im2 ;
        
        % Optional: flip the image
        dFmip = fliplr(dFmip) ;        
        imagesc(dFmip)

        % Formatting and output
        caxis([-1, 1])
        colormap blueblackred
        axis off
        axis equal
        title(['t = ' num2str(t * timeInterval) ' ' timeUnits])
        cbar = colorbar() ;
        ylabel(cbar, '$I_{sqh}- I_{mem}$', 'Interpreter', 'latex')
        caxis([-clip, clip])
        % caxis(max(abs(dFmip(:))) * [-1,1])
        saveas(gcf, dFfn0)
    end
end


%% Take ratio of two channels, averaged over DV
chRGB = [3, 1] ;
clip = 2 ;
modN = 5 ;
colors = parula(ceil(length(xp.fileMeta.timePoints) / modN)) ;
% TIFF filename output
dFfn0 = fullfile('./mips', ['diffDV_compare.png']) ;
dFfn1 = fullfile('./mips', ['diffDV_substract.png']) ;
kk = 1 ;
if ~exist(dFfn0, 'file') || true
    for tidx = 1:40
        t = xp.fileMeta.timePoints(tidx) ;
        disp(['t = ', num2str(t)])
        % OPTION 2: Load the mips and subtract
        mipMatF = fullfile('./mips', 'foreground_mat', sprintf(['mip_' file16name '_foreground.mat'], t)) ;
        load(mipMatF, 'imF') ;
        im1 = squeeze(imF{1}) ;
        im2 = squeeze(imF{2}) ;
        mean1 = double(mean(im1(im1 > 0))) ;
        mean2 = double(mean(im2(im2 > 0))) ;
        im1 = double(im1) / mean1 ;
        im2 = double(im2) / mean2 ;
        midband = round(0.3*size(im1, 1)):round(0.66*size(im1, 1));
        I1 = sum(im1(midband, :), 1) ;
        I2 = sum(im2(midband, :), 1) ;
        dFmip = I1 - I2 ;
        
        % Optional: flip the image
        I1 = fliplr(I1) ;
        I2 = fliplr(I2) ;
        dFnet = fliplr(dFnet) ;
        
        if t == xp.fileMeta.timePoints(1)
            I1net = I1 ;
            I2net = I2 ;
            dFnet = dFmip ;
        end
        I1net = I1net + I1 ;        
        I2net = I2net + I2 ;   
        dFnet = dFnet + dFmip ;
        
        if mod(tidx-1, modN) == 0
            if kk == 1
                Io1 = double(mean(I1net)) ;
                Io2 = double(mean(I2net)) ;
            end
            % Formatting and output
            xL = linspace(0, 1, length(I1)) ;
            
            figure(1)
            subplot(2,1,1)
            plot(xL, double(I1net) / Io1, 'color', colors(kk, :))
            hold on;
            subplot(2,1,2)
            plot(xL, double(I2net) / Io2, 'color', colors(kk, :))    
            hold on;
            
            figure(2)
            plot(xL, double(dFnet) / Io1, 'color', colors(kk, :))
            hold on
            I1net = 0*I1net ;
            I2net = 0*I2net ;
            kk = kk + 1;
        end
    end
    
    figure(1)
    subplot(2, 1, 1)
    ylabel('$I_{sqh}$ [au]', 'Interpreter', 'latex')
    subplot(2, 1, 2)
    ylabel('$I_{mem}$ [au]', 'Interpreter', 'latex')
    xlabel('AP position, $x/L$', 'interpreter', 'latex')
    saveas(gcf, dFfn0)
    
    figure(2)
    % legend({'sqh', 'membrane'})
    ylabel('$I_{sqh} - I_{mem}$ [au]', 'Interpreter', 'latex')
    xlabel('AP position, $x/L$', 'interpreter', 'latex')
    saveas(gcf, dFfn1)
end

%% Show cross section
% % Now set the meta data in the experiment.
fileMeta.filenameFormat = '202010191649_T%02d_foreground.tif' ;
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();

slices = 220:230 ;
clearvars im
for t = xp.fileMeta.timePoints
    % TIFF filename output
    outdir = ['./mips_sliceY' num2str(slices(1)) '_' num2str(slices(end))] ;
    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end
    outmip = fullfile(outdir, sprintf(['mip_' file16name '_sliceY.tif'], t)) ;
    
    mip_exist = exist(outmip, 'file') ;
    if ~mip_exist 
        % Load the data first to set time
        xp.loadTime(t);
        xp.rescaleStackToUnitAspect();
        IV = xp.stack.image.apply() ;
        
        for ch = 1:length(IV)
            % Make mip of this channel with designated colors
            im{ch} = squeeze(max(squeeze(IV{ch}(slices, :, :)), [], 1)) ;
            im{ch} = double(im{ch}) ./ double(max(im{ch}(:))) ;
        end
        
        % Write MIP as RGB with falsecolors
        colors = [0, 1, 1; ...  % cyan
                  1, 0, 0];     % red
        imout = uint8(zeros(size(im{1}, 1), size(im{1}, 2), 3)) ;
        maxI = 255 ;
        for ch = 1:length(IV)
            imout(:, :, 1) = imout(:, :, 1) + uint8(maxI * im{ch} * colors(ch, 1)) ;
            imout(:, :, 2) = imout(:, :, 2) + uint8(maxI * im{ch} * colors(ch, 2)) ;
            imout(:, :, 3) = imout(:, :, 3) + uint8(maxI * im{ch} * colors(ch, 3)) ;
        end
        imwrite(permute(imout, [2, 1, 3]), outmip)
    end
end

%% COMPARING FLOW FIELDS MASTER -- single dataset
% Isaac Breinyn 2020
%
% Run setup.m before this script
%
% This is a script that loads multichannel confocal images and then (with
% user input and training along the way) can analyze the likeness of the
% flow between the two channels.
%
% Naming conventions must remain consistent across all datasets for this
% pipeline to work reliably. Access the flowExperiment.m file for
% information on how to name files.
%
% NOTE: The PIV data used in this script is formated into arrays of u and v
% component of the velocities for EACH channel. Therefore, arrays with a
% first dim of 4 will be in the order of u1, v1, u2, then v2. The letter
% corresponds to the component of the velocity vector, and the number to
% the prevelant channel.

%% Run setup.m from imsane path 
% setup

%% Clean Matlab and add Paths
clear; close all; clc ;

% add paths to where the code lives
% addpath_recurse('L:\\Streichan\\code\\gut_matlab\\') ;
% addpath_recurse('L:\\Streichan\\code\\') ;
addpath_recurse('/mnt/data/code/gut_matlab/') ;
addpath('/mnt/data/code/isaac_code/')

%% WARNING: Users should ONLY make changes to this section of the script.

%%% Assign paramaters for class flowExperiment %%%

% options.dataDir = 'L:\\Streichan\\data\\202003111630_mef2gal4klarUASCAAXmChHiFP_wo_great\\1644_folding_30s_2um_la8_4p5zoom_63x\\' ; % the directory where the data lives
options.dataDir = fullfile([filesep 'mnt'],'data','confocal_data',...
    'gut', 'relative_motion', 'mef2Gal4klarUASCAAXmChHiFP', ...
    '202003111630_mef2gal4klarUASCAAXmChHiFP_wo_great',...
    '1644_folding_30s_2um_la8_4p5zoom_63x') ; % the directory where the data lives
options.dataSetName = '1644_folding' ; % the "label" of this dataset (the prefix to the file names)
options.lastTimePoint = 32 ; % "last timepoint" or the number of TPs in your dataset
options.fileSize = [512 512] ; % the size of the file
options.trackFileName = fullfile('Ch2_MIPs','1644_folding_Ch2_T01_MIP_smoothed_stack-myData_Object Identities.h5') ; % the name of the object identities file (from ilastik)
options.res = 0.0802 ; % the data resolution in um/px
options.lobeSplit = 256 ; % define where you want to put the barrier that defines where the two lobes seperate
options.dt = 0.5 ;  % time resolution in minutes

% todo: add mirrored boolean

save(fullfile(options.dataDir, 'options.mat'), 'options')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(options.dataDir)
%% Instantiate the class flowExperiment and populate necessary properties

flex = flowExperiment(options) ;

%% The first step is to acquire data that you would like to analyze. This data must be split into two channels and seperate timepoints.
% NOTE: MAKE SURE DATA IS IN THIS ORIENTATION
%
%      D
%  A      P
%      V
%
% Open LIF in ImageJ, Flip/rotate in imageJ, save as separate tiff for each
% timepoint and channel. 
% 
% Once in the correct orientation, rescale the data to unit aspect ratio 
% and then convert into h5s for ilastik training using 
% muscle_surface_pipeline.m. Only the membrane channel (channel 1)
% should be used for ilastik training. This is the step necessary to remove
% any glial cells or non-gut structures from the FOV. 

%% ILASTIK: Train on ch1 and extract probabilities
% foreground channel (ch1 of the training) is the membrane, background (ch2
% of the training) is else. 

%% Mask data and then extract MIPs 
% Once you have the probabilities for channel 1 of your data, you can mask
% both channels and make MIPs of the data. These MIPs will then be used
% later in the pipeline.

flex.makeMaskedMIPs() ; 

% After running, you should have a collection of masked tiffs as well as
% MIPs of those masked tiffs.

%% Smooth the MIPs from last step
% This gets rid of any salt-and-pepper noise and allows for higher quality
% PIV analysis.

flex.smoothMIPs() ;
