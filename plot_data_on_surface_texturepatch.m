%% Plot surfaces with intensity data
% Plot data on a surface using texturepatch for a dataset.
%
% By Dillon Cislo & NPMitchell
%==========================================================================


%% INITIALIZE ImSAnE PROJECT ==============================================
%
% We start by clearing the memory and closing all figures
clear; close all; clc;

% We then add some necessary code to the path
addpath(genpath('/mnt/crunch/djcislo/MATLAB/euclidean_orbifolds'));
addpath(genpath('/mnt/data/code/gptoolbox'));
addpath(genpath('/mnt/crunch/djcislo/MATLAB/TexturePatch'));
addpath_recurse('/mnt/data/code/gut_matlab/geometry/')

% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.
dataDir = [ pwd filesep ] ;
% dataDir = [ '/mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/', ...
%     'Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/' ];
meshDir = fullfile(dataDir, ['msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1' filesep]);
projectDir = [ pwd filesep ] ;
% [ projectDir, ~, ~ ] = fileparts(matlab.desktop.editor.getActiveFilename); 
% cd( projectDir );

%% DEFINE MESH STACK ======================================================
meshFileNameBase = fullfile( meshDir, 'mesh_apical_stab_%06d.ply' );
msls_axis_order = 'yxzc';
normal_step = 10;

%% CREATE EXPERIMENT ======================================================
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

%% LOAD MESH STACK ========================================================

% % Cell arrays to hold the mesh struct attributes
% v = cell(length(xp.fileMeta.timePoints),1);
% f = cell(length(xp.fileMeta.timePoints),1);
% vn = cell(length(xp.fileMeta.timePoints),1);
% 
% for t = xp.fileMeta.timePoints
%     disp(['Loading mesh ' num2str(t) ] )
%     
%     % Convert into timepoint ID
%     tidx = xp.tIdx(t);
%  
%     % Specfiy the mesh file to load
%     meshFile = sprintf( meshFileNameBase, t );
%       
%     % Read in the mesh file -----------------------------------------------
%     mesh = read_ply_mod( meshFile );
%     
%     if strcmp( msls_axis_order, 'yxzc' )
%         mesh.v = mesh.v( :, [3,2,1] ) ;
%         mesh.vn = mesh.vn( :, [3,2,1] ) ;
%     end
%     
%     % Make sure vertex normals are normalized
%     mesh.vn = mesh.vn ./ sqrt( sum( mesh.vn.^2, 2 ) );
%     
%     % Normally evolve vertices
%     mesh.v = mesh.v + normal_step .* mesh.vn;
%     
%     % Re-orient faces
%     mesh.f = reorient_facets( mesh.v, mesh.f );
%      
%     v{tidx} = mesh.v;
%     f{tidx} = mesh.f;
%     vn{tidx} = mesh.vn;
%     
% end
% disp('done loading meshes')
% 
% % Create a stack of all meshes
% % meshStack = struct( 'v', v, 'f', f, 'vn', vn );
% 
% clear meshFileNameBase meshFile
% clear msls_axis_order 


%% Determine camera angles from rotation matrix
rotname = fullfile(meshDir, 'rotation_APDV') ;
transname = fullfile(meshDir, 'translation_APDV') ;
xyzlimname = fullfile(meshDir, 'xyzlim_APDV_um') ;
% Load transformations
% Load the rotation matrix
rot = importdata([rotname '.txt']) ;
% Load the translation to put anterior to origin
trans = importdata([transname '.txt']) ;
% Load plotting limits
xyzlim = dlmread([xyzlimname '.txt'], ',', 1, 0) ;
% camrot = inv(rot) ;
% [yaw, pitch, roll] = yawpitchroll_from_rotationmatrix(camrot) ;

%% PLOT THE TEXTURED MESH IN 3D FOR FIRST TIME POINT ======================
do_singletp = false ;
if do_singletp
    % Load the first timepoint
    xp.loadTime(xp.fileMeta.timePoints(1));
    xp.rescaleStackToUnitAspect();
    
    % Specfiy the mesh file to load
    meshFile = sprintf( meshFileNameBase, t );
      
    % Read in the mesh file -----------------------------------------------
    mesh = read_ply_mod( meshFile );
    
    if strcmp( msls_axis_order, 'yxzc' )
        mesh.v = mesh.v( :, [3,2,1] ) ;
        mesh.vn = mesh.vn( :, [3,2,1] ) ;
    end
    
    % Make sure vertex normals are normalized
    mesh.vn = mesh.vn ./ sqrt( sum( mesh.vn.^2, 2 ) );
    % Normally evolve vertices
    mesh.v = mesh.v + normal_step .* mesh.vn;
    % Re-orient faces
    mesh.f = reorient_facets( mesh.v, mesh.f );
    
    % Psize is the linear dimension of the grid drawn on each triangular face
    Options.PSize = 5;
    Options.EdgeColor = 'none';
    Options.Rotation = rot ;

    IV = xp.stack.image.apply();
    IV = imadjustn(IV{1});

    % First args are physical vertices, then texture faces (same number as 
    % physical faces, but connectivity can be different), then texture
    % vertices, which can also be different. The vertices are in row, column, 
    % page format, so x and y get flipped. IV is the texture volume.
    % Options.PSize 
    texture_patch_3d( mesh.f, mesh.v, ...
        mesh.f, mesh.v(:, [2 1 3]), IV, Options );

    axis equal

    colormap bone

    clear Options IV
end

%% PLOT ALL TEXTURED MESHES IN 3D =========================================
% Get limits and create output dir
% Name output directory
nstepexten = strrep([sprintf('%04d', normal_step) 'step'], '-', 'n') ;
figoutdir = [fullfile(meshDir, ['images_patchorbifolds' nstepexten]) filesep];
if ~exist(figoutdir, 'dir')
    mkdir(figoutdir) ;
end
figddir = [fullfile(figoutdir, 'dorsal') filesep];
figvdir = [fullfile(figoutdir, 'ventral') filesep];
figlat1dir = [fullfile(figoutdir, 'lateral1') filesep];
figlat2dir = [fullfile(figoutdir, 'lateral2') filesep];
dirs = {figddir, figvdir, figlat1dir, figlat2dir} ;
for i = 1:length(dirs)
    if ~exist(dirs{i}, 'dir')
        mkdir(dirs{i}) ;
    end
end

buffer = 20 ;
resolution = min(xp.fileMeta.stackResolution) ;  % in um
xmin = xyzlim(1, 1); xmax = xyzlim(1, 2) ;
ymin = xyzlim(2, 1); ymax = xyzlim(2, 2) ;
zmin = xyzlim(3, 1); zmax = xyzlim(3, 2) ;
xmin = xmin - buffer ;
ymin = ymin - buffer ;
zmin = zmin - buffer ;
xmax = xmax + buffer ;
ymax = ymax + buffer ;
zmax = zmax + buffer ;


%% Now draw for all TPs
% Psize is the linear dimension of the grid drawn on each triangular face
Options.PSize = 5;
Options.EdgeColor = 'none';
Options.Rotation = rot ;
Options.Translation = trans ;
Options.Dilation = resolution ;
% figure parameters
xwidth = 16 ; % cm
ywidth = 10 ; % cm

for t = xp.fileMeta.timePoints
    tic 
    % Get the data in 3d for this timepoint
    % if xp.currentTime ~= t
    tidx = xp.tIdx(t) ;
    xp.loadTime(t) ;
    % end
    disp(['Applying image...'])
    IV = xp.stack.image.apply();
    IV = imadjustn(IV{1});

    % Specfiy the mesh file to load
    meshFile = sprintf( meshFileNameBase, t );
      
    % Read in the mesh file -----------------------------------------------
    mesh = read_ply_mod( meshFile );
    
    if strcmp( msls_axis_order, 'yxzc' )
        mesh.v = mesh.v( :, [3,2,1] ) ;
        mesh.vn = mesh.vn( :, [3,2,1] ) ;
    end
    
    % Make sure vertex normals are normalized
    mesh.vn = mesh.vn ./ sqrt( sum( mesh.vn.^2, 2 ) );
    % Normally evolve vertices
    mesh.v = mesh.v + normal_step .* mesh.vn;
    % Re-orient faces
    mesh.f = reorient_facets( mesh.v, mesh.f );
    
    % First args are physical vertices, then texture faces (same number as 
    % physical faces, but connectivity can be different), then texture
    % vertices, which can also be different. The vertices are in row, column, 
    % page format, so x and y get flipped. IV is the texture volume.
    % Options.PSize 
    % mesh.f = f{tidx} ;
    % mesh.v = v{tidx} ;
    % mesh.vn = vn{tidx} ;
    
    fig = figure('Visible', 'Off') ;
    % set(gcf, 'Visible', 'Off') 
    disp(['creating texture patch ' num2str(t, '%06d')])
    texture_patch_3d( mesh.f, mesh.v, ...
        mesh.f, mesh.v(:, [2 1 3]), IV, Options );
    
    % format the figure
    disp('formatting figure...')
    axis equal
    colormap bone
    xlim([xmin, xmax])
    ylim([ymin, ymax])
    zlim([zmin, zmax])
    title(['$t = $' num2str(t) ' min'], 'Interpreter', 'Latex', 'Color', 'white') 
    xlabel('AP position [$\mu$m]', 'Interpreter', 'Latex', 'Color', 'white')
    ylabel('DV position [$\mu$m]', 'Interpreter', 'Latex', 'Color', 'white')
    zlabel('lateral position [$\mu$m]', 'Interpreter', 'Latex', 'Color', 'white')
    
    % Rotate the camera angle using rotation and translation 
    % camorbit(theta, phi)
    % camup() --> ax.CameraUpVector = [sin(45) cos(45) 1]
    % camroll(dtheta)
    % This is an active rotation
    % rotate(hSurface,direction,25)
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperPosition', [0 0 xwidth ywidth]);
    
    % Make background black
    set(gca,'Color','k')
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'Color', 'k')
        
    % Make tick labels white
    labeltype = {'XTickLabel', 'YTickLabel', 'ZTickLabel'} ;
    for q = 1:2
        % get the current tick labeks
        ticklabels = get(gca, labeltype{q});
        % prepend a color for each tick label
        ticklabels_new = cell(size(ticklabels));
        for i = 1:length(ticklabels)
            ticklabels_new{i} = ['\color{white} ' ticklabels{i}];
        end
        % set the tick labels
        set(gca, labeltype{q}, ticklabels_new);
    end
    
    % Capture all four views
    disp(['saving figure...' num2str(t, '%06d')])
    view(0, 90)
    fnd = ['patchorbifold_dorsal_' num2str(t, '%06d') '.png'] ;
    saveas(fig, fullfile(figddir, fnd))
    view(0, 270)
    fnv = ['patchorbifold_ventral_' num2str(t, '%06d') '.png'] ;
    saveas(fig, fullfile(figvdir, fnv))
    
    % Lateral views
    view(0, 0)
    % make white z tick labels
    q = 3;
    % get the current tick labeks
    ticklabels = get(gca, labeltype{q});
    % prepend a color for each tick label
    ticklabels_new = cell(size(ticklabels));
    for i = 1:length(ticklabels)
        ticklabels_new{i} = ['\color{white} ' ticklabels{i}];
    end
    % set the tick labels
    set(gca, labeltype{q}, ticklabels_new);
    
    fnv = ['patchorbifold_lateral1_' num2str(t, '%06d') '.png'] ;
    saveas(fig, fullfile(figlat1dir, fnv))
    view(0, 180)
    fnv = ['patchorbifold_lateral2_' num2str(t, '%06d') '.png'] ;
    saveas(fig, fullfile(figlat2dir, fnv))
    close all
    toc
end

clear Options IV