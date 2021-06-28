%% Segment out a Monge-form "slab" [z0(x,y), z1(x,y)] from confocal data 
% NPMitchell 2020
%
% This is a pipeline to segment confocal data & take MIPs

% temporary path def
% cd /mnt/data/confocal_data/gut/2020/mef2GAL4klarUASsqhGFPUASCAAXmCh/202010191649_mef2GAL4klarUASsqhGFPUASCAAXmCh_40x1p6x_240spf_0p4um_5to10pc3p5to7pc_oilImmersol/

if ismac
    gutdir = '~/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/' ;
    textureDir = '~/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/TexturePatch/' ;
else
    gutdir = '/mnt/data/code/gut_matlab/' ;
    textureDir = '/mnt/data/code/TexturePatch_for_git/TexturePatch/' ;
end
addpath(fullfile(gutdir, 'addpath_recurse/'))
addpath_recurse(gutdir)
addpath(textureDir)

% dataDir = '/mnt/data/optogenetics_confocal/' ;
% dataDir = [dataDir 'antpGAL4/huygens_deconvolution_withKlar/'] ;

% WT
if ismac
    dataDir = '~/Desktop/gut/antibodies/antp8C11_1to50_ms568_1to500_PFAM_water/tif_stacks/' ;
else
    dataDir = '/mnt/data/antibodies/antp/antp8C11_1to50_ms568_1to500_PFAM_water/tif_stacks/';
end

cd(dataDir)

% We start by clearing the memory and closing all figures
clear; close all; clc;

%% INITIALIZE ImSAnE PROJECT ==============================================
%
% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.

dataDir    =  cd; 
projectDir = dataDir ;

% A filename base template - to be used throughout this script
% the 32 bit fn
fn = '' ;
% the 16 bit fn
% file16name = 'antpOCRLgap43_T%03d' ;     
file16name = '20200225_antp8C11_1to50_ms568_1to500_PFAM_2um_e*.tif' ;
fns = dir(fullfile(dataDir, file16name)) ;

% 20200225 resolution
pix2um = 0.455 ;
um2pix = 1 / pix2um ;
resolution = [pix2um, pix2um, 0.75] ;
timepoints = 1 ;

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
fileMeta.timePoints         = 1 ;
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
expMeta.description         = 'antp stains';
expMeta.dynamicSurface      = 0;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.fitTime             = fileMeta.timePoints(first_tp);
% expMeta.detectorType = 'surfaceDetection.planarEdgeDetector';
% expMeta.fitterType = 'surfaceFitting.tpsFitter';
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';
expMeta.fitterType          = 'surfaceFitting.tpsFitter';

%% INSTANTIATE EXPERIMENT CLASS
for ii = 1:length(fns)
    
    preview = true ;
    overwrite = false ;

    % Output files for this h5 file / embryo 
    outpngFn = ['./texturePatches/layer0_' fns(ii).name(1:end-4) '_c%01d.png'] ;
    outPfn1 = sprintf(outpngFn, 1) ;
    outPfn2 = sprintf(outpngFn, 2) ;
    outtifFn = ['./texturePatches/layer0_' fns(ii).name(1:end-4) '_c%01d.tif'] ;
    outTfn1 = sprintf(outtifFn, 1) ;
    outTfn2 = sprintf(outtifFn, 2) ;
    if ~exist(outPfn1, 'file') || ~exist(outPfn2, 'file') || ...
            ~exist(outTfn2, 'file') || ~exist(outTfn2, 'file') || overwrite

        % Now set the meta data in the experiment.
        fileMeta.filenameFormat = fns(ii).name ;
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
        ssfactor = 4;
        niter = 25 ;
        niter0 = 115 ;
        ofn_smoothply = 'mesh_' ;
        ofn_ply = 'mesh_ms_' ; 
        ofn_ls = 'msls_' ;
        meshlabCodeDir = '/mnt/data/code/meshlab_codes/' ;
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
        %xp.loadTime(xp.fileMeta.timePoints(first_tp));
        %xp.rescaleStackToUnitAspect();

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
                    'fileName', fns(ii).name(1:end-4), ... % the filename of h5 to train on
                    'mslsDir', mslsDir, ...  % the directory for all output data/images
                    'ofn_ls', fns(ii).name(1:end-4), ...  % the output filename for level sets
                    'ofn_ply', fns(ii).name(1:end-4), ... % the output filename for PLY files
                    'ms_scriptDir', ms_scriptDir, ... % the directory containing run_morphsnakes.py
                    'timepoint', 0, ... % which timepoint in the data to consider
                    'zdim', 3, ... % Which dimension is the z dimension
                    'pre_nu', pre_nu, ... % number of dilation/erosion passes for positive/negative values
                    'pre_smoothing', pre_smoothing, ... % number of smoothing passes before running MS
                    'ofn_smoothply', [fns(ii).name(1:end-4) '_smooth'],... % the output file name (not including path directory)
                    'mlxprogram', fullfile(meshlabCodeDir, mlxprogram), ... % the name of the mlx program to use to smooth the results. Note that if mesh_from_pointcloud==true, should take obj as input and mesh as output.
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
                    'smooth_with_matlab', -1, ... % if <0, use meshlab. If >0, smooth the mesh after marching cubes mesh creation using matlab instead of mlxprogram, with diffusion parameter lambda = this value. If =0, no smoothing.
                    'pythonVersion', '2') ; 

        % Set detect options ------------------------------------------------------
        xp.setDetectOptions( detectOptions );

        % clear msls_exten imwriteOptions saveDir
        % clear channel foreGroundChannel
        % clear niter niter0 lambda1 lambda2
        % clear exit_thres smoothing nu
        % clear post_nu post_smoothing

        %% Adjust LUT to be approx constant over time & Isotropic resolution
        convert2IsotropicCombined = false ;
        if convert2IsotropicCombined
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
                    for dmyk = 1:length(IV)
                        IVii = IV{dmyk} ;
                        vlo = double(prctile( IVii(:) , adjustlow_pctile )) / double(max(IVii(:))) ;
                        vhi = double(prctile( IVii(:) , adjusthigh_pctile )) / double(max(IVii(:))) ;
                        disp(['--> ', num2str(vlo), ', ', num2str(vhi), ...
                            ' for ', num2str(adjustlow), '/', num2str(adjusthigh)])
                        IV{dmyk} = imadjustn(IVii, [double(vlo); double(vhi)]) ;
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
            fileMeta.filenameFormat = 'wt_T%03d.tif' ;
            fileMeta.stackResolution = min(fileMeta.stackResolution) * [1,1,1] ;
            xp.setFileMeta(fileMeta);
            xp.setExpMeta(expMeta);
            xp.initNew();
        end


        %% CREATE THE SUBSAMPLED H5 FILE FOR INPUT TO ILASTIK =====================
        % skip if already done
        xp.setFileMeta(fileMeta);
        xp.setExpMeta(expMeta);
        xp.initNew();

        h5fn = fullfile(projectDir, [fns(ii).name(1:end-4) '.h5']) ;
        if ~exist(h5fn, 'file')
            disp(['Did not find file: ', h5fn])

            % Only load and rescale if multiple timepoints/channels
            % if length(xp.fileMeta.timePoints) > 1 && t > xp.fileMeta.timePoints(1)
            xp.loadTime(1);
            xp.rescaleStackToUnitAspect();
            % end

            % make a copy of the detectOptions and change the fileName
            detectOpts2 = detectOptions ;
            detectOpts2.fileName = fns(ii).name(1:end-4) ;
            xp.setDetectOptions( detectOpts2 );
            xp.detector.prepareIlastik(xp.stack);
            disp(['done outputting downsampled data h5: im#' num2str(ii) ' for surface detection'])
        else
            disp(['h5 im#' num2str(ii) ' was already output, skipping...'])
        end
        disp('Open with ilastik if not already done')


        %% TRAIN DATA IN ILASTIK TO IDENTIFY SURFACE ==========================
        % % open ilastik, train until probabilities and uncertainty are satisfactory
        % 

        %% Create MorphSnakesLevelSet from the Probabilities from ilastik ========
        % Skip if already done
        % Now detect all surfaces
        run_full_dataset_ms = false ;
        detectOptions.run_full_dataset = 'none' ;  % projectDir ; % 'none' ;  % override here

        % Morphosnakes for all remaining timepoints INDIVIDUALLY ==============
        xp.loadTime(1) ;
        xp.rescaleStackToUnitAspect();
        xp.detectSurface();
        error('here')

        %% Continue to mesh smoothing
        xp.detectSurface();

        %% Load mesh, subsample, triangulate, Texturepatch
        IV = xp.stack.image.apply() ;
        sz1 = size(IV{1}, 1) ;
        sz2 = size(IV{1}, 2) ;
        sz3 = size(IV{1}, 3) ;
        nU = round(150 * sz1 / sz2) ;
        nV = 150 ;
        uminmax = [0, sz1] ;
        vminmax = [0, sz2] ;
        %%
        layerWidth = 5 ; % 1 for RFP
        subsample = 10 ;
        zOffset = 5 ; % 20 ;
        axisOrder = 'yxz' ;

        axisOrder = erase(axisOrder, 'c') ;
        zdim = find(axisOrder== 'z') ;

        %% Onion stacks
        nPos = 15 ;
        nNeg = 10 ;

        layerWidth = 5 ;
        midLayerOffset = -4 ;
        lam = 0.01 ;  % 1 ; 

        outtifFn = ['./texturePatches/layer0_' fns(ii).name(1:end-4) '_c%01d.tif'] ;
        if ~exist(sprintf(outtifFn, 2), 'file') || overwrite
            try
                IV = xp.stack.image.apply() ;
            catch
                xp.loadTime(1);
                xp.rescaleStackToUnitAspect();
                IV = xp.stack.image.apply() ;
            end
            sz1 = size(IV{1}, 1) ;
            sz2 = size(IV{1}, 2) ;
            sz3 = size(IV{1}, 3) ;
            if contains(lower(axisOrder), 'xyz')
                szX = sz1 ; szY = sz2 ; szZ = sz3 ;
            elseif contains(lower(axisOrder), 'yxz')
                szX = sz2 ; szY = sz1 ; szZ = sz3 ;
            else
                error('handle here')
            end

            % Could load original mesh
            meshfn = fullfile(mslsDir, [fns(ii).name(1:end-4) '_smooth000000.ply']) ;
            mesh0 = read_ply_mod(meshfn) ;
            
            
            % Filter a bit
            
            % for e1:
            % zOffset = 5 ;
            % gridOptions = struct('uminmax', uminmax, 'vminmax', vminmax, ...
            %     'nU', nU, 'nV', nV, 'takeMedianMaxMin', 'interpolate', ...
            %     'szX', szX, 'szY', szY, 'zOffset', zOffset) ;
            % planeOptions = struct('maxIsPlane', false, 'thresPlane', 0.1, 'z0', 0) ;
            
            % For e2
            zOffset = 0 ;
            gridOptions = struct('uminmax', uminmax, 'vminmax', vminmax, ...
                'nU', nU, 'nV', nV, 'takeMedianMaxMin', 'interpolate', ...
                'szX', szX, 'szY', szY, 'zOffset', zOffset, ...
                'erosionRadius', 2, 'dilationRadius', 0) ;
            planeOptions = struct('maxIsPlane', true, 'thresPlane', 1, 'z0', sz3) ;
            [mesh, out_of_frame] = medianFilterMongeMesh(mesh0, gridOptions, planeOptions, preview) ;

            % Fix vertices that are off mesh or near the mesh boundary
            disp('smoothing mesh...')
            tri = triangulation(mesh.f, mesh.v) ;
            bndId = neighbors(tri, out_of_frame); % find(zr > max(medz(:))-5)) ;
            bndId = unique(bndId(~isnan(bndId))) ;
            bndId = tri.ConnectivityList(bndId, :) ;
            bndId = unique(bndId(:)) ;
            vsm0 = mesh.v ;
            mesh.v = laplacian_smooth(mesh.v, mesh.f, 'cotan', bndId, lam, ...
                'implicit', mesh.v, 1000);

            if preview
                disp('previewing smoothed mesh')
                clf
                tri = triangulation(mesh.f, mesh.v) ;
                trisurf(tri, 'edgecolor', 'none')
                view(2); title('smoothed triangulation of rescaled resampling')
                axis equal
            end

            m2d = mesh ;
            m2d.v = mesh.v(:, [1,2]) ;
            mesh.vn = zeros(size(mesh.v)) ;
            mesh.vn(:, zdim) = 1 ;
            % mesh.vn(out_of_frame, 3) = 0 ;
            mesh.v = mesh.v ;
            if preview
                trisurf(triangulation(mesh.f, mesh.v), mesh.vn(:, 3), 'edgecolor', 'none')
                title('normal vectors')
            end

            % TexturePatch
            Opts = struct() ;
            if strcmpi(axisOrder, 'xyz')
                Opts.imSize = [sz2, sz1] ;
            elseif strcmpi(axisOrder, 'yxz')
                Opts.imSize = [sz1, sz2] ;
            else
                error('handle here')
            end
            Opts.numLayers = [nPos, nNeg] ;
            Opts.layerSpacing = 1 ;
            Opts.vertexNormal = mesh.vn ;
            Opts.extrapolationMethod = 'nearest' ; % 'nearest' ;'none'; 
            if contains(lower(axisOrder), 'zyx')
                [ patchIm, imref, zeroID, MIP, SIP ] = ...
                    texture_patch_to_image(mesh.f, m2d.v, mesh.f, ...
                    mesh.v(:, [2,1,3]), IV, Opts) ;
            elseif contains(lower(axisOrder), 'yxz')
                IVtmp = IV ;
                IVtmp{1} = permute(IVtmp{1}, [2,1,3]) ;
                IVtmp{2} = permute(IVtmp{2}, [2,1,3]) ;
                [ patchIm, ~, zeroID, MIP, SIP ] = ...
                    texture_patch_to_image(mesh.f, m2d.v(:, [2,1]), mesh.f, ...
                    mesh.v(:, [2,1,3]), IVtmp, Opts) ;   
                patchIm = permute(patchIm, [2,1,3,4]) ;
                MIP = permute(MIP, [2,1,3]);
                SIP = permute(SIP, [2,1,3]) ;
                imagesc(permute(MIP, [2,1,3])); pause(1)
            elseif contains(lower(axisOrder), 'xyz')
                IVtmp = IV ;
                IVtmp{1} = permute(IVtmp{1}, [2,1,3]) ;
                IVtmp{2} = permute(IVtmp{2}, [2,1,3]) ;
                [ patchIm, imref, zeroID, MIP, SIP ] = ...
                    texture_patch_to_image(mesh.f, m2d.v(:, [2,1]), mesh.f, ...
                    mesh.v(:, [1,2,3]), IVtmp, Opts) ;            
            else    
                error('handle here')
            end
            if preview
                imshow(permute(squeeze(max(patchIm, [], 4)), [2, 1, 3]))
            end

            if preview
                figure(1)
                for page = 1:size(patchIm, 4)
                    imshow(mat2gray(squeeze(patchIm(:, :, 2, page)), [0, 1])')
                    title(['page #' num2str(page)])
                    pause(0.05)
                end
                % Show reference page
                imshow(mat2gray(squeeze(patchIm(:, :, 2, zeroID+midLayerOffset)), [0, 1])')
            end

            % Look at figure
            if preview
                figure(2)
                trisurf(triangulation(mesh.f, mesh.v), mesh.v(:, 3), 'edgecolor', 'none')
                view(2); axis equal
                pause(0.01)
            end

            % imshow(patchIm(:,:,2, zeroID))
            % caxis([0, 0.2])

            if ~exist('./texturePatches/', 'dir')
                mkdir('./texturePatches/')
            end
            ch1 = squeeze(patchIm(:, :, 1, :)) ;
            ch2 = squeeze(patchIm(:, :, 2, :)) ;

            ch1 = uint8(mat2gray(ch1, [min(ch1(:)), max(ch1(:))]) * 2^8) ;
            ch2 = uint8(mat2gray(ch2, [min(ch2(:)), max(ch2(:))]) * 2^8) ;

            ch1 = permute(ch1, [2, 1, 3]) ;
            ch2 = permute(ch2, [2, 1, 3]) ;

            chs = {reshape(ch1, [size(ch1, 1), size(ch1, 2), 1, size(ch1, 3)]), ...
                reshape(ch2, [size(ch2, 1), size(ch2, 2), 1, size(ch2, 3)])} ;

            % Save image stack
            for ch = 1:2
                outfn = sprintf(outtifFn, ch) ;
                writeTiff5D(chs{ch}, outfn, 8) ;
            end

            % Save mesh
            save(['./texturePatches/mesh_' fns(ii).name(1:end-4) '.mat'], ...
                'mesh', 'm2d', 'Opts', 'axisOrder') ;
        else
            disp(['already on disk: image#' num2str(ii)])
        end

        outpngFn = ['./texturePatches/layer0_' fns(ii).name(1:end-4) '_c%01d.png'] ;
        outfn1 = sprintf(outpngFn, 1) ;
        outfn2 = sprintf(outpngFn, 2) ;
        if ~exist(outfn1, 'file') || ~exist(outfn2, 'file')
            imfn1 = sprintf(outtifFn, 1) ;
            imfn2 = sprintf(outtifFn, 2) ;
            ims1 = loadtiff(imfn1) ;
            minLayer = max(nNeg-layerWidth+midLayerOffset, 1) ;
            maxLayer = min(nNeg+layerWidth+midLayerOffset, size(ims1, 3)) ;
            im01 = squeeze(max(ims1(:, :, minLayer:maxLayer), [], 3)) ;
            im01 = uint8(255 * mat2gray(im01, double([min(im01(:)), max(im01(:))]))) ;
            imwrite(im01, outfn1)
            ims2 = loadtiff(imfn2) ;
            im02 = squeeze(max(ims2(:, :, minLayer:maxLayer), [], 3)) ;
            im02 = uint8(255 * mat2gray(im02, double([min(im02(:)), max(im02(:))]))) ;
            imwrite(im02, outfn2)
        end
    end
    
end

%% AP variation in expression -- quantification
% To Quantify ap variation of antp expression, rotate/flip image so that 
% anterior is left, dorsal is up and save as
% ./texturePatches_aligned/layer0...png
afns = dir('./texturePatches_aligned/layer0*c2.png') ;
for qq = 1:length(afns)
    im = imread(fullfile(afns(qq).folder, afns(qq).name)) ;
    for pp = 1:2
        imshow(im)
        roi = drawrectangle ;
        roi = roi.Position ;
        x0 = roi(1) ;
        x1 = roi(1) + roi(3) ;
        y0 = roi(2) ;
        y1 = roi(2) + roi(4) ;
        
        save(fullfile(afns(qq).folder, 'rois', ...
            sprintf([afns(qq).name(1:end-4) 'rois_%d.mat'], pp)), 'roi') ;
        
        cutim = im(y0:y1, x0:x1) ;
        imshow(cutim)
        cutim(cutim < 10) = NaN ;
        plot(nanmedian(cutim, 2))
        sum(cutim, 3)
    end
end
