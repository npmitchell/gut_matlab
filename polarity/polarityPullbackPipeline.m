%% Polarity Pullback Pipeline =============================================
% Measure anisotropy via Radon transform on pullback stack. Here, use
% 'axisymmetric' pullbacks of the growing Drosophila midgut
%
% Execute from the projectDir, where the data is. 
% By NPMitchell 2019
%==========================================================================

clear; close all; clc;

%% Options
% seriestype defines what pullback we are considering for the Radon transf.
% e: extended, es: extended_shifted, rs: relaxed_shifted, 
seriestype = 'e'; 

%% Parameters
overwrite = true ;
save_ims = true ;
normal_shift = 10 ;
a_fixed = 2 ;
patch_width = 30 ;
preview = false ;
washout2d = 0.5 ;
washout3d = 0.5 ;
colorwheel_position = [.8 .01 .15 .15] ;

%% Add paths
% Add some necessary code to the path (ImSAnE should also be setup!) ------
addpath(genpath('/mnt/crunch/djcislo/MATLAB/euclidean_orbifolds'));
addpath(genpath('/mnt/data/code/gptoolbox'));
addpath(genpath('/mnt/data/code/gut_matlab/TexturePatch'));
addpath(genpath('/mnt/data/code/gut_matlab/polarity'));
addpath(genpath('/mnt/data/code/gut_matlab/PeakFinding'));
addpath_recurse('/mnt/data/code/imsaneV1.2.3/external/') ;
addpath_recurse('/mnt/data/code/gut_matlab/plotting/') ;
% addpath(genpath('/mnt/crunch/djcislo/MATLAB/TexturePatch'));

%% Define some colors
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

%% Initialize ImSAnE Project ==============================================

% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.
dataDir = pwd ;

projectDir = dataDir ;
% [ projectDir, ~, ~ ] = fileparts(matlab.desktop.editor.getActiveFilename); 
cd( projectDir );

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
fileNameBase = 'Time_%06d_c1_stab';

fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fileNameBase, '.tif'];
fileMeta.nChannels          = 1;
fileMeta.timePoints         = 110:263;
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

expMeta                     = struct();
expMeta.channelsUsed        = 1;
expMeta.channelColor        = 1;
expMeta.description         = 'Apical membrane in Drosophila gut';
expMeta.dynamicSurface      = 1;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';

% Now set the meta data in the experiment.
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();

clear fileMeta expMeta


%% Initialize Some Directory Definitions ==================================

% The top level data directory
meshDir = fullfile(dataDir, 'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1') ;

% The file name base for the full meshes
fullMeshBase = fullfile( meshDir, 'mesh_apical_stab_%06d.ply' );

% The file name base for the cylinder meshes
cylinderMeshBase = fullfile( meshDir, ...
    'cylindercut/mesh_apical_stab_%06d_cylindercut.ply' );

% The file constaing the AD/PD points
dpFile = fullfile( meshDir, ...
    'cylindercut/ap_boundary_dorsalpts.h5' );

% The dataset name base for the AD points
ADBase = '/mesh_apical_stab_%06d/adorsal';

% The dataset name base for the PD points
PDBase = '/mesh_apical_stab_%06d/pdorsal';

% The folder where the pullback images will be saved
nshift = strrep(sprintf('%03d', normal_shift), '-', 'n') ;
imFolder = fullfile(meshDir, ['PullbackImages_' nshift 'step'] ) ;
imFolder_e = [imFolder '_extended' filesep] ; % debug
imFolder_r = [imFolder '_relaxed' filesep] ; % debug
imFolder_re = [imFolder '_relaxed_extended' filesep] ;
pivDir = fullfile(meshDir, 'piv') ;
polDir = fullfile(meshDir, ['polarity' filesep 'radon' ]) ;
imFolder_es = [imFolder '_extended_shifted' filesep] ;
% The extensile scale factor in x for relaxing the mesh
arfn = fullfile(imFolder_r, 'ar_scalefactors.h5') ;

tomake = {polDir} ;
for ii = 1:length(tomake)
    dir2make = tomake{ii} ;
    if ~exist( dir2make, 'dir' )
        mkdir(dir2make);
    end
end

%% Load Pullback Mesh Stack ===============================================
% Check if cutmeshes already saved
mstckfn = fullfile(meshDir, 'meshStack_orbifold.mat') ; % debug
if exist(mstckfn, 'file') 
    load(mstckfn)
else
    msg = ['Did not find ', mstckfn] ;
    msg = [msg '--> Run Generate_Axisymmetric_Pullbacks_Orbifold.m first'];
    error(msg)
end
disp('done loading meshStack')

%% Get timestamps for the images with pullbacks (from extended imFolder)
fns = dir(fullfile(imFolder_e, '*.tif')) ;
npiv = length(fns) ;
time = zeros(npiv, 1);
meshidx = zeros(npiv, 1) ;
for ii=1:npiv
    tmp = split(fns(ii).name, '.tif') ;
    tmp = split(tmp{1}, '_') ;
    time(ii) = str2double(tmp{2}) ;
    meshidx(ii) = time(ii) - time(1) + 1 ;
end
dt = diff(time) ;
if length(time) < 1
    error('Could not load timepoint IDs')
end
disp('done building dt') 

%% PREPARE FOR RADON ======================================================
if contains(seriestype, 's')
    shiftyfn = fullfile(pivDir, 'shifty.mat') ;
    if exist(shiftyfn, 'file')
        load(shiftyfn)
    else
        disp('Loading pivresults_orbifold_pass0.mat ...')
        load(fullfile(pivDir, 'pivresults_orbifold_pass0.mat'))

        %% Subtract off the mean flow in y for each frame ========================= 
        meanv = zeros(length(v_filtered), 1) ;
        for ii=1:length(v_filtered)
            tmp = v_filtered{ii} ;
            meanv(ii) = mean(tmp(:)) ;
        end
        dy = round(meanv) ;
        shifty = [0; -cumsum(dy)] ;
        disp('computed shifty')
        save(shiftyfn, 'shifty')
    end
    disp('done loading/computing shifty')
else
    shifty = zeros(length(fns), 1);
    disp(['no shifty since seriestype is ' seriestype])
end

%% Load stretch values (width of mesh.urelax / width of mesh.u)
% Try to load ar from h5
ar = zeros(length(meshStack), 1) ;
for kk = 1:length(time)
    t = time(kk) ;
    try
        ar(t - time(1)) = h5read(arfn, ['/' sprintf('%06d', t)]) ;
    catch
        disp(['Could not find time ' num2str(t) ' in arfn'])
    end
end

if ~any(ar)
    disp('Could not load from h5! Loading from meshStack')
    ar = zeros(length(meshStack), 1) ;
    for ii = 1:length(meshStack)
        if isfield(meshStack{ii}, 'urelax')
            ar(ii) = max(meshStack{ii}.urelax(:, 1)) ;
        end
    end
    if all(ar == 0)
        error('Could not load ar (relaxed scalefactors)')
    end
end
disp('done loading ar (relaxed scalefactors)')


%% Get image sizes of images for gridding and mapping =====================
% Obtain the size in x and y of images to consider
if strcmp(seriestype, 'es')
    fns = dir(strrep(fullfile([imFolder_es, '/', fileNameBase, '.tif']), '%06d', '*')) ;
elseif strcmp(seriestype, 're')
    fns = dir(strrep(fullfile([imFolder_re, '/', fileNameBase, '.tif']), '%06d', '*')) ;
elseif strcmp(seriestype, 'e')
    fns = dir(strrep(fullfile([imFolder_e, '/', fileNameBase, '.tif']), '%06d', '*')) ;
end
% Write the images to disk and get their sizes while being written 
imsizes = zeros(length(fns), 2) ;
for ii=1:length(fns)
    outfn = fullfile(fns(ii).folder, fns(ii).name) ;
    % Note that shifting was done as follows:
    % disp(['Reading ' fns(i).name])
    % fileName = split(fns(i).name, '.tif') ;
    % fileName = fileName{1} ;
    % im = imread(fullfile(fns(i).folder, fns(i).name)) ;
    % im = circshift(im, shifty(i), 1) ;
    % imsizes(i, :) = size(im) ;
    % imwrite( im, outfn, 'TIFF' );
    im = imread(outfn) ;
    imsizes(ii, :) = size(im) ;
end
disp('done with reading image sizes')

%% RADON TRANSFORM ========================================================
if strcmp(seriestype, 'es')
    fns = dir(strrep(fullfile([imFolder_es, '/', fileNameBase, '.tif']), '%06d', '*')) ;
    stretch = true ;
elseif strcmp(seriestype, 're')
    fns = dir(strrep(fullfile([imFolder_re, '/', fileNameBase, '.tif']), '%06d', '*')) ;
    stretch = false ;
elseif strcmp(seriestype, 'e')
    fns = dir(strrep(fullfile([imFolder_e, '/', fileNameBase, '.tif']), '%06d', '*')) ;
    stretch = true ;
end

% Iterate over all images
nemsz = 15 ;  % size of nematic bar to draw in images
w = 80 ; % patch_width ;
step = 30 ;
for ii=1:length(fns)
    tidx = time(ii) ; % timestamp as integer
    timestr = sprintf('%03d', tidx) ;
    disp(['t = ' timestr])
    % Grab filename
    fileName = split(fns(ii).name, '.tif') ;
    fileName = fileName{1} ;
    im = imread(fullfile(fns(ii).folder, fns(ii).name)) ;
    image_max = max(im(:)) ;
    
    % disp('TRANSPOSING IMAGE for debug')
    % im = im' ;
    % disp('Rotatin gimage for debug')
    % im = imrotate(im, 90) ;
    
    % Get scale of image
    xsc0 = imsizes(ii, 1) ;
    ysc0 = imsizes(ii, 2) ;
    
    % Chop up the image into little chunks
    xx = w:step:(xsc0 - w) ;
    yy = w:step:(ysc0 - w) ;
    
    % Load radon for this timept or compute it
    % Use two different methods for comparison
    for res = 2
        disp(['Computing/loading polarity using algorithm ' sprintf('%01d', res)])
        options.res = res ;
        wstepstr = ['_w' sprintf('%04d', w) '_step' sprintf('%04d', step)] ;
        fn = ['polarity_' seriestype wstepstr '_res' sprintf('%01d', res)] ;
        radonfn = fullfile(polDir, fn) ;
        
        % Check for the results of this method (res = 1,2)
        poldir2d_res = fullfile(polDir, ['polarity2d_' seriestype '_res' sprintf('%01d', res) wstepstr ]) ;
        if ~exist(poldir2d_res, 'dir')
            mkdir(poldir2d_res)
        end
        if stretch        
            poldir2d_res_stretch = fullfile(polDir, ['polarity2d_' seriestype '_res' sprintf('%01d', res) '_stretchcorr' wstepstr]) ;
            if ~exist(poldir2d_res_stretch, 'dir')
                mkdir(poldir2d_res_stretch)
            end
        end
        poldir2d_res_color = [poldir2d_res '_color'] ;
        if ~exist(poldir2d_res_color, 'dir')
            mkdir(poldir2d_res_color)
        end
        
        % Figure out if we load the data from disk or compute it
        load_from_disk = false ;
        if exist(radonfn, 'file')
            try 
                h5read(radonfn, ['/' fileName '/angles_smoothcorr']) ;
                load_from_disk = true && ~overwrite;
            catch
            end
        end
        
        % Load the results if they exist, otherwise compute them
        if load_from_disk 
            disp('Loading results from h5...')
            angles = h5read(radonfn, ['/' fileName '/angles']) ;
            magnitudes = h5read(radonfn, ['/' fileName '/magnitudes']) ;
            angles_smooth = h5read(radonfn, ['/' fileName '/angles_smooth']) ;
            mag_smooth = h5read(radonfn, ['/' fileName '/mag_smooth']) ;
            xx = h5read(radonfn, ['/' fileName '/xx']) ;
            yy = h5read(radonfn, ['/' fileName '/yy']) ;
        else
            disp('Computing radon peaks...')
            % Preallocate
            angles = zeros(length(xx), length(yy)) ;
            magnitudes = zeros(length(xx), length(yy)) ;

            % Compute radon transform for each little chunk
            for j = 1:length(xx)
                for k = 1:length(yy)
                    xi = xx(j);
                    yi = yy(k) ;
                    xmin = max(1, xi - w) ;
                    xmax = min(xsc0, xi + w) ;
                    ymin = max(1, yi - w) ;
                    ymax = min(ysc0, yi + w) ;
                    chunk = im(xmin:xmax, ymin:ymax) ;
                    [xsz, ysz] = size(chunk) ;
                    xcenter = xsz * 0.5 ;
                    ycenter = ysz * 0.5 ;

                    % Mask out a circle from the patch
                    % create a xygrid
                    [xc, yc] = meshgrid(1:xsz, 1:ysz) ;
                    dist = (xc - xcenter) .^2 + (yc - ycenter) .^2 ;
                    mask = dist' < w^2 ;
                    chunk = chunk .* uint8(mask) ;
                    % imagesc(chunk)

                    % Compute radon as a function of angle
                    [angle, magnitude, results] = extractRadonNematic(chunk, options) ;
                    
                    % Store angle and magnitude of this patch in array
                    angles(j, k) = angle ;
                    magnitudes(j, k) = magnitude ;
                end        
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% Save data to h5 file
            try
                h5create(radonfn, ['/' fileName '/angles'], size(angles)) ;
            catch
                disp('angles already exists')
            end
            try
                h5create(radonfn, ['/' fileName '/magnitudes'], size(angles)) ;
            catch
                disp('magnitudes already exists')
            end
            try
                h5create(radonfn, ['/' fileName '/xx'], size(xx)) ;
            catch
                disp('xx already exists')
            end
            try
                h5create(radonfn, ['/' fileName '/yy'], size(yy)) ;
            catch
                disp('yy already exists')
            end
            
            h5write(radonfn, ['/' fileName '/angles'], angles) ;
            h5write(radonfn, ['/' fileName '/magnitudes'], magnitudes) ;
            h5write(radonfn, ['/' fileName '/xx'], xx) ;
            h5write(radonfn, ['/' fileName '/yy'], yy) ;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Stretch the vectors if the image is not relaxed
            if stretch
                scale_factor = ar(tidx - time(1) + 1) ;
                angles_corrected = atan2(sin(angles), scale_factor * cos(angles)) ;
                % save it
                try
                    h5create(radonfn, ['/' fileName '/angles_corrected'], size(angles_corrected)) ;
                catch
                    disp('angles_corrected already exists')
                end
                h5write(radonfn, ['/' fileName '/angles_corrected'], angles_corrected) ;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            %% Smooth the nematic 
            filt = (fspecial('gaussian', 5, 1)); %if needed modify the filter according to the expected peaks sizes
            % First smooth the image
            % Note: multiply by 2 to give continuity around 2pi
            nxs = conv2(cos(2 * angles), filt, 'same') ;
            nys = conv2(sin(2 * angles), filt, 'same') ;
            mag_smooth = conv2(magnitudes, filt, 'same') ;
            
            % Convert back to nematic order from polar
            angles_smooth = atan2(nys, nxs) * 0.5 ;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Save the smoothed fields
            try
                h5create(radonfn, ['/' fileName '/angles_smooth'], size(angles_smooth)) ;
            catch
                disp('angles_smooth already exists')
            end
            
            try
                h5create(radonfn, ['/' fileName '/mag_smooth'], size(mag_smooth)) ;
            catch
                disp('mag_smooth already exists')
            end            
            h5write(radonfn, ['/' fileName '/angles_smooth'], angles_smooth) ;
            h5write(radonfn, ['/' fileName '/mag_smooth'], mag_smooth) ;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Stretch the vectors if the image is not relaxed
            if stretch
                % First smooth the image
                % Note: multiply by 2 to give continuity around 2pi
                nxs = conv2(cos(2 * angles_corrected), filt, 'same') ;
                nys = conv2(sin(2 * angles_corrected), filt, 'same') ;
                
                % Convert back to nematic order from polar
                angles_smoothcorr = atan2(nys, nxs) * 0.5 ;

                % save it
                try
                    h5create(radonfn, ['/' fileName '/angles_smoothcorr'], size(angles_smoothcorr)) ;
                catch
                    disp('angles_smoothcorr already exists')
                end
                h5write(radonfn, ['/' fileName '/angles_smoothcorr'], angles_smoothcorr) ;
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Normalize magnitudes  
        magnitudes = magnitudes / median(magnitudes(:)) ;
        mag_smooth = mag_smooth / median(mag_smooth(:)) ;
        
        % Load the image filename for saving after putting nematic on top
        fileName = split(fns(ii).name, '.tif') ;
        fileName = fileName{1} ;
        
        % Save an image overlay
        outimfn = fullfile(poldir2d_res, [fileName '.png']) ;
        if ~exist(outimfn, 'file') || overwrite
            figure('units', 'normalized', ...
                'outerposition', [0 0 1 1], 'visible', 'off')
            
            imshow(im * washout2d + image_max * (1-washout2d))
            xlims = xlim ;
            ylims = ylim ;
            hold on
            % Transpose everything
            xv = nemsz * magnitudes .* cos(angles) ;
            yv = nemsz * magnitudes .* sin(angles) ;
            xvt = xv';
            yvt = yv';
            [x0, y0] = meshgrid(xx, yy) ;
            x0q = x0 - 0.5 * xvt ;
            y0q = y0 - 0.5 * yvt ;
            % quiver(x0q(:), y0q(:), xv(:), yv(:), 0, 'ShowArrowHead', 'off') ;
            % scatter(y0(:), x0(:), 'r.')
            quiver(y0q(:), x0q(:), yvt(:), xvt(:), 0, 'ShowArrowHead', 'off') ;
            axis equal
            xlim(xlims) ;
            ylim(ylims) ;
            
            % Extract image from figure axes
            patchIm = getframe(gca);
            % print('-dpng','-r300', outimfn)
            disp(['Saving image to ' outimfn]) 
            imwrite( patchIm.cdata, outimfn );
            close all
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Save image with corrected angles if stretched
        if stretch
            outimfn = fullfile(poldir2d_res_stretch, [fileName '_corr.png']) ;
            if ~exist(outimfn) || overwrite
                disp('Saving stretch corrected image')
                figure('units', 'normalized', ...
                    'outerposition', [0 0 1 1], 'visible', 'off')
                imshow(im * washout2d + image_max * (1-washout2d))
                xlims = xlim ;
                ylims = ylim ;
                hold on
                [x0, y0] = meshgrid(xx, yy) ;
                xv = nemsz * mag_smooth .* cos(angles_smoothcorr) ;
                yv = nemsz * mag_smooth .* sin(angles_smoothcorr) ;
                xvt = xv';
                yvt = yv';
                x0q = x0 - 0.5 * xvt ;
                y0q = y0 - 0.5 * yvt ;
                % quiver(x0q(:), y0q(:), xv(:), yv(:), 0, 'ShowArrowHead', 'off') ;
                % scatter(y0(:), x0(:), 'r.')
                quiver(y0q(:), x0q(:), yvt(:), xvt(:), 0, 'ShowArrowHead', 'off') ;
                axis equal
                xlim(xlims) ;
                ylim(ylims) ;
                % Extract image from figure axes
                patchIm = getframe(gca);
                % print('-dpng','-r300', outimfn)
                disp(['Saving image to ' outimfn]) 
                imwrite( patchIm.cdata, outimfn );
                close all
            end
        end
        
        %% Save colorplots (heatmap)    
        close all
        outimfn = fullfile(poldir2d_res_color, [fileName '.png']) ;
        % colors for colorplots
        [colors, ~] = define_colors(3) ;
        yellow = colors(3, :) ;
        cmap = interpolate_3color_cmap(0:0.01:1, 5, yellow, 4, false) ;
        if ~exist(outimfn, 'file') || overwrite
            disp('Saving colored image')
            
            figure('units', 'normalized', ...
                'outerposition', [0 0 1 1], 'visible', 'off')
            if stretch
                cos2t = cos(2 * angles_smoothcorr) ;
            else
                cos2t = cos(2 * angles_smooth) ;
            end
            options.alpha = magnitudes / max(magnitudes(:)) ;
            % options.alpha = 1;
            hh = heatmap_on_alphaimage(im, xx, yy, cos2t, options) ;
            
            axis equal
            axis off 
            colormap(gca, cmap)
            
            % Colorwheel as grid
            ax3 = axes('Position', colorwheel_position) ;
            [qq, pp] = meshgrid(-1:0.01:1, -1:0.01:1) ;
            % Note purposeful transposition here!
            color = cos(2 * atan2(qq, pp)) * 0.5  + 0.5 ;
            radius = vecnorm([qq(:), pp(:)]') ;
            alpha = reshape(radius, size(qq)) ;
            alpha(radius > 1) = 0 ;
            h = imagesc(qq(:), pp(:), color) ;
            set(h, 'AlphaData', alpha) ;
            colormap(gca, cmap)
            axis equal
            axis off
            title({'Radon-based', 'anisotropy'}, 'Fontweight', 'normal')

            % Save the image
            disp(['Saving figure ' outimfn])
            saveas(gcf, outimfn)
        end
        
    end
end
disp('done with radon peaks')

%% Plot average over time of cos(2theta) and sin(2theta)

% Use two different methods for comparison
% todo: also try https://www.mathworks.com/help/images/ref/imhmax.html
for res = 1:2
    resstr = ['_res' sprintf('%01d', res)] ;
    cos2t_medians = zeros(length(ar), 1) ;
    sin2t_medians = zeros(length(ar), 1) ;
    cos2t_unc = zeros(length(ar), 1) ;
    sin2t_unc = zeros(length(ar), 1) ;
    % Creates edges for histogram plot
    bin_step = 0.1 ;
    edges = -1:bin_step:1 ;
    midpts = edges(1:end-1) + bin_step * 0.5 ;
    % Pre-allocate for binning
    ccount = zeros(length(ar), length(edges) - 1) ;
    scount = zeros(length(ar), length(edges) - 1) ;
    for ii=1:length(fns)
        tidx = time(ii) ; % timestamp as integer
        indx = tidx - time(1) + 1 ; % index of this timepoint in full_time
        timestr = sprintf('%03d', tidx) ;
        disp(['t = ' timestr])
        % Grab filename
        fileName = split(fns(ii).name, '.tif') ;
        fileName = fileName{1} ;

        % Load radon for this timept or compute it
        disp(['Computing/loading polarity using algorithm ' sprintf('%01d', res)])
        options.res = res ;
        fn = ['polarity_' seriestype '_w' sprintf('%04d', w) ...
            '_step' sprintf('%04d', step) resstr '.h5' ] ;
        radonfn = fullfile(polDir, fn) ;

        % Load results for this fn
        disp('Loading results from h5...')
        angles_smooth = h5read(radonfn, ['/' fileName '/angles_smooth']) ;
        if stretch
            angles_smoothcorr = h5read(radonfn, ['/' fileName '/angles_smoothcorr']) ;
        end
        mag_smooth = h5read(radonfn, ['/' fileName '/mag_smooth']) ;
        xx = h5read(radonfn, ['/' fileName '/xx']) ;
        yy = h5read(radonfn, ['/' fileName '/yy']) ;

        % convert to cos2t, sin2t
        cos2t = cos(2 * angles_smooth) ;
        sin2t = sin(2 * angles_smooth) ;
        
        % Here just take simple medians 
        cos2t_medians(indx) = nanmedian(cos2t(cos2t ~= 0)) ;
        tmp = nanstd(cos2t(cos2t ~= 0)) ;
        cos2t_unc(indx) = tmp ; % / sqrt(length(cos2t(cos2t ~= 0))) ;

        sin2t_medians(indx) = nanmedian(sin2t(sin2t ~= 0)) ;
        tmp = nanstd(sin2t(sin2t ~= 0)) ;
        sin2t_unc(indx) = tmp ; % / sqrt(length(sin2t(sin2t ~= 0))) ;
        
        % Also take histogram
        ccount(indx, :) =  histcounts(cos2t(:), edges) ;
        scount(indx, :) =  histcounts(sin2t(:), edges) ;
    end

    cos2t_medians(cos2t_medians == 0 ) = NaN ;
    sin2t_medians(sin2t_medians == 0 ) = NaN ;

    % errorbar(aniso_medians, aniso_std)
    close all
    tx = 1:length(ar) ;
    clower = cos2t_medians - cos2t_unc ;
    cupper = cos2t_medians + cos2t_unc ;
    % sin(2t) * l
    slower = sin2t_medians - sin2t_unc ;
    supper = sin2t_medians + sin2t_unc ;

    colors = define_colors();
    blue = colors(1, :)  ;
    red = uint8(colors(2, :) * 255) ;
    % lightblue = [149 / 255, 208 / 255, 252 / 255] ;
    goodidx = find(~isnan(clower)) ;
    txg = tx(goodidx) ;
    clower = clower(goodidx) ;
    cupper = cupper(goodidx) ;
    slower = slower(goodidx) ;
    supper = supper(goodidx) ;
    cf = fill([txg, fliplr(txg)], [clower', fliplr(cupper')], blue, 'LineStyle', 'none') ;
    set(cf, 'facealpha', .1)
    hold on;
    sf = fill([txg, fliplr(txg)], [slower', fliplr(supper')], red, 'LineStyle', 'none') ;
    set(sf, 'facealpha', .1)
    ch = plot(tx, cos2t_medians, '.', 'Color', blue) ;
    sh = plot(tx, sin2t_medians, '.', 'Color', red) ;
    xlabel('time [min]')
    ylabel('$\langle \cos(2\theta) \rangle$, $\langle \sin(2\theta) \rangle$', 'interpreter', 'latex')
    title('Tissue anisotropy via Radon transform')
    % xlim([0, max(time)]) ;
    legend([ch, sh], {'$\langle \cos(2\theta) \rangle$', ...
        '$\langle \sin(2\theta) \rangle$'}, 'location', 'best', 'Interpreter', 'Latex')
    % ylim([min(aniso_medians - aniso_unc), max(aniso_medians + aniso_unc)])
    zeroh = plot([0, length(ar)], [0, 0], 'k--') ;
    set(get(get(zeroh,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off'); % Exclude line from legend

    % save it
    outfn_time = fullfile(polDir, ['radon_anisotropy' resstr '.png']) ;
    disp(['Saving figure to ' outfn_time])
    saveas(gcf, outfn_time)
    close all

    % Save data
    datfn = fullfile(polDir, ['radon_anisotropy_timeseq' resstr '.mat']) ;
    disp(['Saving time sequence data to ' datfn])
    save(datfn, 'cos2t_medians', 'cos2t_unc', 'sin2t_medians', 'sin2t_unc')
    
    % Save as histogram over time 
    % Histogram for Cos(2t)
    ccount(ccount == 0) = NaN ;
    ccount = fillmissing(ccount, 'nearest');
    imagesc(time - time(1), midpts, ccount')
    xlabel('time [min]')
    ylabel('$\langle \cos(2\theta) \rangle$', 'interpreter', 'latex')
    title('Tissue anisotropy via Radon transform')
    ax = gca;
    ax.YDir = 'normal' ;
    xlim([-Inf, max(time) - time(1)])
    % save it
    outfn_time = fullfile(polDir, ['radon_anisotropy' resstr '_histcos.png']) ;
    disp(['Saving figure to ' outfn_time])
    saveas(gcf, outfn_time)
    
    % Histogram for Sin(2t)
    scount(scount == 0) = NaN ;
    scount = fillmissing(scount, 'nearest');
    imagesc(time - time(1), midpts, scount')
    xlabel('time [min]')
    ylabel('$\langle \sin(2\theta) \rangle$', 'interpreter', 'latex')
    title('Tissue anisotropy via Radon transform')
    ax = gca;
    ax.YDir = 'normal' ;
    xlim([-Inf, max(time) - time(1)])
    % save it
    outfn_time = fullfile(polDir, ['radon_anisotropy' resstr '_histsin.png']) ;
    disp(['Saving figure to ' outfn_time])
    saveas(gcf, outfn_time)
    
    error('break')
end
error('break')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAKE MAP FROM PIXEL TO XYZ =============================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Making map from pixel to xyz to compute orientation in 3d...')
% Get position of embedding points associated with velocities 
% Compute size of domain (0,1), (0, 1)
fns = dir(strrep(fullfile([imFolder_es, '/', fileNameBase, '.tif']), '%06d', '*')) ;
im = imread(fullfile(fns(ii).folder, fns(ii).name)) ;
% size of extended image
esize = size(im) ;
% map extended image size to (0, 1), (-0.5, 1.5)
xesz = esize(1) ;
yesz = esize(2) ;
% map from pixel y to network y
fy = @(y) 2.0 * y / yesz - 0.5 ;
pix2u = @(x) 2.0 * x / yesz ;
% map from network y to pixel y
% Note that for some reason we need to flip the image
u2pix = @(x, ysc) ysc / 2.0 * x ;
y2pix = @(y, h, ysc) ysc - u2pix(y + 0.5, ysc) + h ;

pol3dfn = fullfile(polDir, 'polarity3D.mat') ;

if exist(pol3dfn, 'file')
    load(pol3dfn)
else
    pol3d = cell(length(fns), 1) ;
    
    %% Iterate over all images with flow fields
    for ii=1:length(fns) - 1
        timestr = sprintf('%03d', time(ii)) ;
        disp(['t = ' timestr])

        % Get scale of image
        ysc0 = imsizes(ii, 2) ;
        ysc1 = imsizes(ii+1, 2) ;

        mesh0 = meshStack{meshidx(ii)} ;    
        mesh1 = meshStack{meshidx(ii + 1)} ;
        assert(time(ii) + dt(ii) == time(ii+ 1))

        % Load the positions of the velocity vectors in pixels
        x0 = x{ii} ;
        y0 = y{ii} ;
        uu = u_filtered{ii} ;
        vv = v_filtered{ii} ; 
        
        % get embedded vector in R^3 for t0
        % Obtain the equivalent of v2D: ie the 2D vertices: mesh0.u
        mesh0x = x2pix(mesh0.u(:, 1), ysc0) ;
        mesh0y = y2pix(mesh0.u(:, 2), ysc0) ;

        % Create extended mesh (copied above and below)
        meshxy = [mesh0x, mesh0y ] ;
        mabove = [mesh0x, mesh0y + x2pix(1., ysc0)] ;
        mbelow = [mesh0x, mesh0y - x2pix(1., ysc0)] ;
        mabove2 = [mesh0x, mesh0y + 2 * x2pix(1., ysc0)] ;
        mbelow2 = [mesh0x, mesh0y - 2 * x2pix(1., ysc0)] ;
        m0xy = [meshxy; mabove; mbelow; mabove2; mbelow2] ;
        % mesh faces for t0 concatenated = mf0c
        mf0 = mesh0.f ;
        mf0c = [mf0; mf0 + length(mesh0x); mf0 + 2 * length(mesh0x); ...
            mf0 + 3 * length(mesh0x); mf0 + 4 * length(mesh0x)] ;
        tr0 = triangulation(mf0c, m0xy) ;
        [t0_contain, baryc0] = pointLocation(tr0, [x0(:), y0(:)]) ;    
        % Interpolate the position in 3D given relative position within 2D
        % triangle.
        % x123(i) is the x coords of the elements of triangle t_contain(i)
        vxa = mesh0.v(:, 1) ;
        vya = mesh0.v(:, 2) ;
        vza = mesh0.v(:, 3) ;
        assert(size(vxa, 1) == size(mesh0x, 1))
        % Modulo the vertex IDs: trisa are the triangle vertex IDs
        tria = tr0.ConnectivityList(t0_contain, :) ;
        
        % Make sure normals are pointing the right way
        % tmp = faceNormal(tr0)
        % v21 = x0(trisa(:, 2), :) - mesh0.v(trisa(:, 1), :) ;
        % v31 = mesh0.v(trisa(:, 3), :) - mesh0.v(trisa(:, 1), :) ;
        
        trisa = mod(tria, size(vxa, 1)) ;
        trisa(trisa == 0) = size(vxa, 1) ;
        x123a = vxa(trisa) ;
        y123a = vya(trisa) ;
        z123a = vza(trisa) ;
        % Multiply the vertex positions by relative weights.
        % Note that baryc gives weights of the three vertices of triangle
        % t_contain(i) for pointLocation x0(i), y0(i)
        pt0 = [sum(baryc0 .* x123a, 2), sum(baryc0 .* y123a, 2), sum(baryc0 .* z123a, 2) ] ;

        %% Visualize the polarity in 3d ---------------------------------------
        if save_ims
            close all
            fig = figure('Visible', 'Off') ;
            quiver3(pt0(:, 1), pt0(:, 2), pt0(:, 3), ...
                v0(:, 1), v0(:, 2), v0(:, 3), 0)
            axis equal
            title(['t=' timestr]) 
            saveas(gcf, fullfile(polarOutDir, ['piv3d_' timestr '.png']))
            close all
        end
        
        %% Obtain normal & tangential components of velocity --------------
        % Use normal vectors defined on every face by averaging from
        % vertices
        normals1 = mesh0.vn(trisa(:, 1), :) ;
        normals2 = mesh0.vn(trisa(:, 2), :) ;
        normals3 = mesh0.vn(trisa(:, 3), :) ;
        normals = normals1 + normals2 + normals3 ;
        normals = normals ./ sqrt(sum(normals.^2, 2)) ;
        
        % Project the local triangle into the tangent plane as defined by
        % the smoothed normals
        % project each triangle containing a velocity evaluation pt onto 
        % the local smoothed tangent plane
        % triangle vector 1/2/3 projected:
        tv1 = mesh0.v(trisa(:, 1), :) ;
        tv2 = mesh0.v(trisa(:, 2), :) ;
        tv3 = mesh0.v(trisa(:, 3), :) ;
        % triangle vertex positions (vertices 1,2,3 of each face)
        tv1proj = tv1 - dot(tv1 - pt0, normals, 2) .* normals;
        tv2proj = tv2 - dot(tv2 - pt0, normals, 2) .* normals ;
        tv3proj = tv3 - dot(tv3 - pt0, normals, 2) .* normals ;
        % x positions of projected triangles
        x123ap = [tv1proj(:, 1), tv2proj(:, 1), tv3proj(:, 1)] ;
        y123ap = [tv1proj(:, 2), tv2proj(:, 2), tv3proj(:, 2)] ;
        z123ap = [tv1proj(:, 3), tv2proj(:, 3), tv3proj(:, 3)] ;
        
        %% Compute the tangential polarity nematic in plane
        u21 = tv2proj - tv1proj ;
        u31 = tv3proj - tv1proj ;
        w21 = m0xy(tria(:, 2), :) - m0xy(tria(:, 1), :) ;
        w31 = m0xy(tria(:, 3), :) - m0xy(tria(:, 1), :) ;
        
        % Build jacobian for each triangle jjac(triangle index, :, :)
        jac = zeros(size(tv2proj, 1), 2, 3) ;
        jac(:, 1, 1) = w21(:, 1) ./ u21(:, 1) + w31(:, 1) ./ u31(:, 1) ;
        jac(:, 1, 2) = w21(:, 1) ./ u21(:, 2) + w31(:, 1) ./ u31(:, 2) ;
        jac(:, 1, 3) = w21(:, 1) ./ u21(:, 3) + w31(:, 1) ./ u31(:, 3) ;
        jac(:, 2, 1) = w21(:, 2) ./ u21(:, 1) + w31(:, 2) ./ u31(:, 1) ;
        jac(:, 2, 2) = w21(:, 2) ./ u21(:, 2) + w31(:, 2) ./ u31(:, 2) ;
        jac(:, 2, 3) = w21(:, 2) ./ u21(:, 3) + w31(:, 2) ./ u31(:, 3) ;
        
        % Build jacobian for each triangle jac(triangle index, :, :)
        jjac = zeros(size(tv2proj, 1), 3, 2) ;
        jjac(:, 1, 1) = u21(:, 1) ./ w21(:, 1) + u31(:, 1) ./ w31(:, 1) ;
        jjac(:, 1, 2) = u21(:, 1) ./ w21(:, 2) + u31(:, 1) ./ w31(:, 2) ;
        jjac(:, 2, 1) = u21(:, 2) ./ w21(:, 1) + u31(:, 2) ./ w31(:, 1) ;
        jjac(:, 2, 2) = u21(:, 2) ./ w21(:, 2) + u31(:, 2) ./ w31(:, 2) ;
        jjac(:, 3, 1) = u21(:, 3) ./ w21(:, 1) + u31(:, 3) ./ w31(:, 1) ;
        jjac(:, 3, 2) = u21(:, 3) ./ w21(:, 2) + u31(:, 3) ./ w31(:, 2) ;
        
        % I have checked that det(jac * jac') = det(jjac * jjac') for a
        % triangle. 
        
        % Compute 2D veclocities and metric tensor, also dilation
        v0t2d = zeros(size(v0, 1), 2) ;
        g_ab = zeros(size(v0, 1), 2, 2) ;
        dilation = zeros(size(v0, 1), 1) ;
        for qq = 1:size(v0, 1)
            qjac = squeeze(jac(qq, :, :)) ; 
            v0t2d(qq, :) = qjac * v0(qq, :)' ;
            g_ab(qq, :, :) = qjac * qjac' ;
            dilation(qq) = sqrt(det(qjac * qjac')) ;
        end
        
        %% Save the results in datstruct ----------------------------------
        datstruct.pt0 = pt0 ;
        datstruct.pt1 = pt1 ;
        datstruct.v0 = v0 ;
        datstruct.v0n = v0n ;
        datstruct.v0t = v0t ;
        datstruct.v0t2d = v0t2d ;
        pol3d{ii} = datstruct ;

        %% Draw 2D polarity
        if save_ims        
            pdir2d = fullfile(polDir, 'vt2d') ;
            if ~exist(pdir2d, 'dir')
                mkdir(pdir2d)
            end
            polcolordir2d = fullfile(polDir, 'polarity_heatmap') ;
            if ~exist(pheatdir2d, 'dir')
                mkdir(pheatdir2d)
            end
            % Load the image to put flow on top
            fileName = split(fns(ii).name, '.tif') ;
            fileName = fileName{1} ;
            im = imread(fullfile(fns(ii).folder, fns(ii).name)) ;
            figure('units', 'normalized', ...
                'outerposition', [0 0 1 1], 'visible', 'off')
            imshow(im * washout2d + image_max * (1-washout2d))
            xlims = xlim ;
            ylims = ylim ;
            hold on
            quiver(x0(:), y0(:), polv(:, 1), polv(:, 2), 0) ;
            axis equal
            xlim(xlims) ;
            ylim(ylims) ;
            % Extract image from figure axes
            patchIm = getframe(gca);
            outimfn = fullfile(pdir2d, [fileName '.png']) ;
            % print('-dpng','-r300', outimfn)
            imwrite( patchIm.cdata, outimfn );
            
            
            %% Now draw orientation as heatmap
            close all; clear alpha 
            image = mat2gray(im, [0, 256]) ;
            options.caxis = [-10 10] ;
            options.cmap = bwr ;
            v0grid = reshape(v0n, size(x0)) ;
            xfield = x0(1, :)' ;
            yfield = y0(:, 1) ;
            
            % Create the figure
            figure('units', 'normalized', ...
                'outerposition', [0 0 1 1], 'visible', 'off')
            c_handle = heatmap_on_alphaimage(image, xfield, yfield, v0grid, options) ;
            axis equal
            % Extract image from figure axes
            patchIm = getframe(gca);
            % Write figure to file
            outimfn = fullfile(pheatdir2d, [fileName '.png']) ;
            % print('-dpng','-r300', outimfn)
            imwrite( patchIm.cdata, outimfn );

        end
        
        % %% Draw 3D flows
        % if save_ims 
        % flowdir3d = fullfile(pivDir, 'vt3d') ;
        % if ~exist(flowdir3d, 'dir')
        %     mkdir(flowdir3d)
        % end
        % 
        % % Load the image to put flow on top
        % fileName = split(fns(i).name, '.tif') ;
        % fileName = fileName{1} ;
        % im = imread(fullfile(fns(i).folder, fns(i).name)) ;
        % imshow(im * washout2d + max(im) * (1-washout2d))
        % xlims = xlim ;
        % ylims = ylim ;
        % hold on
        % quiver(x0(:), y0(:), v0t2d(:, 1), v0t2d(:, 2), 0) ;
        % axis equal
        % xlim(xlims) ;
        % ylim(ylims) ;
        % outimfn = fullfile(flowdir3d, [fileName '.png']) ;
        % print('-dpng','-r300', outimfn)
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Check that baryc gives weights of the three vertices of triangle
        % t_contain(i) for pointLocation x0(i), y0(i)
        if preview
            close all
            fig = figure('Visible', 'On')
            % triplot(tr0, 'Color', blue, 'LineWidth', 0.0001) 
            hold on
            for j = 1:3
                tritest = tr0.ConnectivityList(t0_contain(j), :); 
                btest = baryc0(j, :);
                triangle = [tritest, tritest(1) ] ;
                plot(m0xy(triangle, 1), m0xy(triangle, 2), 'g.-')
                plot(x0(j), y0(j), 'o')
                xfind = sum(btest' .* m0xy(tritest, 1)) ;
                yfind = sum(btest' .* m0xy(tritest, 2)) ;
                plot(xfind, yfind, '^', 'MarkerSize', 10)        
                axis equal
            end

            % Draw other connections
            quiver(x0, y0, uu, vv, 0)
            plot(x1, y1, 'ko')
            for j = 1:3
                tritest = tr1.ConnectivityList(t1_contain(j), :); 
                btest = baryc1(j, :);
                triangle = [tritest, tritest(1) ] ;
                plot(m1xy(triangle, 1), m1xy(triangle, 2), 'r.-')
                plot(x1(j), y1(j), 'kx')
                xfind = sum(btest' .* m1xy(tritest, 1)) ;
                yfind = sum(btest' .* m1xy(tritest, 2)) ;
                plot(xfind, yfind, '^', 'MarkerSize', 10)        
                axis equal
            end

            waitfor(fig)
            close all        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        clear x0 y0 uu vv x1 y1
        clear pt0 v0 mesh0x mesh0y 
        clear pt1 v1 mesh1x mesh1y
        clear meshxy meshabove meshbelow 
        clear m0xy mf0 mf0c tr0 t0_contain baryc0
        clear m1xy mf1 mf1c tr1 t1_contain baryc1
        clear vxa vya vza x123a y123a z123a
        clear vxb vyb vzb x123b y123b z123b
        clear datstruct
    end 
    save(pol3dfn, 'piv3d') ;
end

