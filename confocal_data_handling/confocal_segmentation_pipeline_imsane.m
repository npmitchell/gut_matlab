%% Segment out a surface from confocal data using planar detector
% NPMitchell 2020
%
% This is a pipeline to segment confocal data & take MIPs

doc surfaceDetection.planarDetector

%% Initialize the project
%
% We start by creating an experiment object, which holds this metadata and 
% provides a frontend for a number of tasks such as loading the raw data.
% To construct the experiment object we need to pass dataDir, the directory 
% containing the raw data and projectDir, the directory where results of 
% the script will be saved.

clear all; close all;

% [scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
% dataDir = fullfile(scriptPath, 'rawData');
% projectDir = fullfile(scriptPath, 'projectFiles');
dataDir = cd ;
projectDir = cd ;

%% For ease of use, convert TIFFs to full resolution first
%
xp = project.Experiment(projectDir, dataDir);

% A filename base template - to be used throughout this script
fn = '202010191649_T%02d' ;                   

fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fn, '.tif'];
fileMeta.nChannels          = 2;
fileMeta.timePoints         = 0:66;
fileMeta.stackResolution    = [0.22674292866082602, 0.22674292866082602, 0.4] ; % the px resolution (found in the .lif; 4 dec places)
fileMeta.swapZT             = 0;

expMeta = struct();
expMeta.description = 'Fixed wing 2hAPF [Ecad, actin, widefield]';
expMeta.channelsUsed = [1 2];
expMeta.channelColor = [1 2];
expMeta.dynamicSurface = false;
expMeta.jitterCorrection = false;
expMeta.detectorType = 'surfaceDetection.planarEdgeDetector';
expMeta.fitterType = 'surfaceFitting.tpsFitter';

xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();

% Make some filler options for now
detectOptions = xp.detector.defaultOptions;
% Calling detectSurface runs the surface detector and creates the point
% cloud in detector.pointCloud.
xp.setDetectOptions(detectOptions);

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

%% Load a time point from the data
%
% Now that we have set up the project, we can load a time point.
% loadTime sets xp.currentTime, loads the stack into xp.stack 
% and resets detector and fitter with the default options for that time.
%
% We rescale to unit aspect ratio to make detection work better later and
% visualize in the right proportions.

xp.loadTime(0);
xp.rescaleStackToUnitAspect();

%% Now use planar detector

% planarDetector.detectSurface detects the surface as the position of the 
% maximal Gaussian z-derivative in some direction, i.e. the position of the
% largest intensity jump along some direction and smoothened over some
% scale.
%
% A number of detection options directly affect detection:
%
% * sigma :     Width of the Gaussian z-derivative.
% * channels :  Channels (summed) to use for detection.
% * zdir :      Dimension corresponding to z, minus flips direction.
% Flipping the direction can sometimes improve detection.
%
% Then there are options which filter the result and can be modified
% without redetecting:
%
% * maxIthresh:     Throw out points with MIP dimmer than this.
% * summedIthresh:  Throw out points with SIP dimmer than this.
% * sigZoutliers:   Remove height outliers after all other masks.
% * scaleZoutliers: Spatial scale of outlier removal.
%
% scaleZoutliers is the linear size of a region over which the
% distribution of height is computed, sigZoutliers is then a cutoff in
% units of standard deviation of this distribution to remove misdetected
% points far above or below the other points in the region.

detectOptions = xp.detector.defaultOptions;
detectOptions.sigma = 1;
detectOptions.zdir = 3;
detectOptions.maxIthresh = 0.1; 
detectOptions.sigZoutliers = 2; 
detectOptions.scaleZoutliers = 30; 

% Calling detectSurface runs the surface detector and creates the point
% cloud in detector.pointCloud.

xp.setDetectOptions(detectOptions);
xp.detectSurface();

% Different from the other detectors, the detected surface is
% represented not only by a PointCloud object but also by an image
% surfaceMatrix, containing z values for each xy.
% Looking at this height map masked by the filters specified in
% detectOptions one can judge how well the surface was detected.

imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [],...
                                            'InitialMagnification', 40);

%% 
% One can then find better filter parameters without redetecting the
% surface by changing the second block of options in detectOptions and 
% calling resetMask and applyMasks. 
% 
% xp.detector.resetMask();
% 
% detectOptions.maxIthresh = 0.1; 
% detectOptions.sigZoutliers = 2; 
% detectOptions.scaleZoutliers = 30; 
% 
% xp.detector.setOptions(detectOptions);    
% xp.detector.applyMasks();
% 
% imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [],...
%                                             'InitialMagnification', 40);

%%
% We can also inspect a point cloud cross section over the data with
% detector.inspectQuality. In the pointCloud option, 'c' specifies the 
% color cyan.

inspectOptions= struct('dimension', 'x', 'value', 200, 'pointCloud', 'c');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%% 
% Or we can look at the point cloud in 3d, with some subsampling factor.
ssfactor = 50;
xp.detector.pointCloud.inspect(ssfactor);


%% Fit the surface for the disc proper cells
%
% By detecting the largest intensity jump along z for each x,y in the
% E-cad channel and filtering out local outliers we have found the apical
% surface of the disc proper cells. We can now fit a smooth surface
% representation to that.
%
% tpsFitter fits the pointcloud using a thin plate spline fit. It has the
% following options:
%
% * gridSize:     Size of grid on which to generate fitted surface
%               default [50 50], full size takes long.
% * smoothing:    TPS smoothing parameter (default 1000).

fitOptions = struct('smoothing', 500, 'gridSize', [100 100]);
xp.setFitOptions(fitOptions);
xp.fitSurface();

%%
% We can visualize the result on a cross section with
% fitter.inspectQuality.

xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%%
% The detector picks up the edge of the E-cad signal but the best read out
% goes solidly through it so we want to move the surface down a little. For
% this we use zEvolve, with a shift specified in pixels.

shift = 12;
xp.zEvolve(shift);

xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%%
% We now generate the Surface Of Interest. The charts to be generated are 
% specified in xp.fitter.charts. In this case there is only one, called
% 'xy'. 

xp.generateSOI();

%% Pull back the data to the surface
% 
% We pull back the data to the SOI using pullbackStack.

xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime);

%%
% To look at the pullback, we call the data field of the SOI at the right
% time, and get a particular patch from that with getPatch. A patch is a 
% part of a surface. In this case, there is only one called xy_index.
% Then we get the data in some patch in a particular coordinate system with
% getTransform. In this case there is only one coordinate system: xy.
% What we get is an object not only holding the image data but also
% metadata and methods to manipulate it. The actual data is obtained by
% calling the method apply. This returns a cell array with entries for each
% channel.

% xp.tIdx converts the time into an index in a list of time points
tidx = xp.tIdx(0);

% the first channel is Ecad
channel = 1;

discProperPatch = xp.SOI.data(tidx).getPatch('xy_index');
discProperImage = discProperPatch.getTransform('xy').apply{channel};
figure, imshow(discProperImage, [], 'InitialMagnification', 50);

%% Save the result
%
% Finally we save the SOI using SOI.save. We set the following options:
%
% * dir:            The directory to save the SOI to.
% * imwriteOptions: Pullbacks are saved to image files using imwrite, we
% can pass options to change file format, compression etc. For example we
% could change this option to
% imwriteOptions = {'jp2', 'Mode', 'lossless'}; 
% * make8bit:       Often absolute intensities don't matter and 8 bit offers
% a large enough dynamic range. This options rescales the lookup table and
% converts to 8 bit before saving.

imwriteOptions = {'tif', 'Compression', 'deflate'};
savedir = fullfile(scriptPath, 'discProperApicalSOI');

options = struct(   'dir',              savedir,...
                    'imwriteOptions',   {imwriteOptions},...
                    'make8bit',         true);
xp.SOI.save(options)



