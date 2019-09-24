%% AXISYMMETRIC_PULLBACK ==================================================
% This is a test of the efficacy of a potiential method for creating
% 'axisymmetric' pullbacks of the growing Drosophila midgut
%
% By Dillon Cislo and Noah Mitchell
%==========================================================================

clear; close all; clc;

% Add some necessary code to the path (ImSAnE should also be setup!) ------
addpath(genpath('/mnt/crunch/djcislo/MATLAB/euclidean_orbifolds'));
addpath(genpath('/mnt/data/code/gptoolbox'));
addpath(genpath('/mnt/crunch/djcislo/MATLAB/TexturePatch'));

%% Initialize ImSAnE Project ==============================================

% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.
dataDir = [ '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/', ...
    'Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/' ];

[ projectDir, ~, ~ ] = fileparts(matlab.desktop.editor.getActiveFilename); 
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
meshDir = [ dataDir ...
    'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1_20190908' ];

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
imFolder = fullfile(projectDir, 'PullbackImages' );
if ~exist( imFolder, 'dir' )
    mkdir(imFolder);
end

%% Iterate Through Time Points to Create Pullbacks ========================

for t = xp.fileMeta.timePoints(2)%xp.fileMeta.timePoints(1:50)
    
    disp(['NOW PROCESSING TIME POINT ', num2str(t)]);
    
    tidx = xp.tIdx(t);
    
    % Load the data for the current time point ----------------------------
    xp.loadTime(t);
    xp.rescaleStackToUnitAspect();

    % Load the cylinder mesh
    mesh = read_ply_mod( sprintf( cylinderMeshBase, t ) );
    
    % Consistently orient mesh faces
    mesh.f = bfs_orient( mesh.f );
    
    % Load he AD/PD vertex IDs
    adIDx = h5read( dpFile, sprintf( ADBase, t ) );
    pdIDx = h5read( dpFile, sprintf( PDBase, t ) );

    % View results --------------------------------------------------------
    % trisurf( triangulation( mesh.f, mesh.v ) );
    % hold on
    % scatter3( mesh.v(adIDx,1), mesh.v(adIDx,2), mesh.v(adIDx,3), ...
    %     'filled', 'r' );
    % scatter3( mesh.v(pdIDx,1), mesh.v(pdIDx,2), mesh.v(pdIDx,3), ...
    %     'filled', 'c' );
    % hold off
    % axis equal

    % Create the cut mesh -------------------------------------------------
    fprintf('Generating Cut Mesh... ');
    [ cutMesh, adIDx, pdIDx ] = ...
        cylinderCutMesh( mesh.f, mesh.v, mesh.vn, adIDx, pdIDx );
    fprintf('Done\n');

    % View results --------------------------------------------------------
    % P = cutMesh.pathPairs(:,1);
    % 
    % trisurf( triangulation( mesh.f, mesh.v ) );
    %
    % hold on
    %
    % line( mesh.v(P,1), mesh.v(P,2), mesh.v(P,3), ...
    %     'Color', 'c', 'LineWidth',2);
    % 
    % scatter3( mesh.v(adIDx,1), mesh.v(adIDx,2), mesh.v(adIDx,3), ...
    %    'filled', 'r' );
    % scatter3( mesh.v(pdIDx,1), mesh.v(pdIDx,2), mesh.v(pdIDx,3), ...
    %     'filled', 'm' );
    % 
    % hold off
    % 
    % axis equal
    % 
    % clear P

    % Generate pullback to rectangular domain ----------------------------------
    % The surface parameterization algorithm (optionally) takes four vertex IDs
    % as input to specify the corners of the square parameterization domain.
    % Maddeningly, the order in which these points are specified to not seem to
    % effect the output. For consistency, we perform a post-hoc correction so
    % that the final output has the following geometric ordering
    %
    %   (AD1)-------(PD1)
    %     |           |
    %     |           |
    %     |           |
    %     |           |
    %   (AD2)-------(PD2)
    %
    % Note that the pathPairs variable has the following columns:
    %   [ ( AD1 -> PD1 ), ( AD2 -> PD2 ) ]
    %--------------------------------------------------------------------------
    
    fprintf('Generating Pullback... ');
    
    % Specify surface parameterization settings ---------------------------
    param = struct();
    param.fixedShape = 2; % Square shaped boundary
    
    ad2IDx = cutMesh.pathPairs(1,2); % The duplicate AD point ID
    pd2IDx = cutMesh.pathPairs(end,2); % The duplicate PD point ID
    
    % The square domain corners
    param.corners = [ adIDx ad2IDx pd2IDx pdIDx ];
    
    % Calculate the surface pullback --------------------------------------
    [ v2D, ~ ] = surface_parameterization( cutMesh.f, cutMesh.v, param );
    
    % Shuffle the boundary to the proper orientation ----------------------
    
    % Find the points associated with the corners
    cornerIDx = pointMatch( [ 0 1; 0 0; 1 0; 1 1 ], v2D );
    
    % Create a copy of the points on the range [-1 1]x[-1 1]
    v2D0 = 2 .* v2D - 1;
    
    % Move the orignal AD point to the (0,1)-position
    v2D0 = complex(v2D0(:,1), v2D0(:,2));
    v2D0 = exp( -1i * pi * (find(cornerIDx == adIDx)-1) / 2 ) .* v2D0;
    v2D0 = [ real(v2D0), imag(v2D0) ];
    
    % Check for proper vertex ordering
    cornerIDx = circshift( cornerIDx, -(find(cornerIDx == adIDx)-1) )';
    
    if cornerIDx(3) ~= pd2IDx
        warning( 'Invalid cut mesh boundary ordering!' );
    end
    
    if cornerIDx(2) ~= ad2IDx
        v2D0 = -[ v2D0(:,2), v2D0(:,1) ];
    end
    
    % Map the updated vertices to the domain [0 1]x[0 1]
    v2D = 0.5 .* (v2D0 + 1);
    
    clear cornerIDx v2D0 ad2IDx pd2IDx
    
    % Fix the aspect ratio of the domain of parameterization to make the
    % mapping as 'isometric' as possible ----------------------------------
    
    if tidx == 1
        a = minimizeIsoarealAffineEnergy( cutMesh.f, cutMesh.v, v2D );
    end
    
    v2D = [ a .* v2D(:,1), v2D(:,2) ];
    
    fprintf('Done\n');
    
    % View results --------------------------------------------------------
    % triplot( triangulation( cutMesh.f, v2D ), 'k' );
    % 
    % hold on
    % 
    % cornerColors = [ 1 0 0; 1 0 1; 0 1 1; 0 1 0 ];
    % scatter( v2D( param.corners, 1 ), v2D( param.corners, 2 ), [], ...
    %     cornerColors, 'filled' );
    %
    % hold off
    %
    % axis equal
    
    % Generate output image file ------------------------------------------
%     fprintf('Generating output image... ');
%     
%     % Texture patch options
%     Options.PSize = 5;
%     Options.EdgeColor = 'none';
%     
%     % Raw stack data
%     IV = xp.stack.image.apply();
%     IV = imadjustn(IV{1});
%     
%     % Make a full screen image
%     figure('units', 'normalized', ...
%         'outerposition', [0 0 1 1], 'visible', 'off')
%     
%     % Plot texture patch cut mesh in 2D
%     texture_patch_3d( cutMesh.f, v2D, ...
%         cutMesh.f, cutMesh.v(:, [2 1 3]), IV, Options );
%     
%     axis equal
%     
%     % Format axes
%     xlim([0 a]); ylim([0 1]);
%     set(gca, 'xtick', []);
%     set(gca, 'ytick', []);
%     
%     colormap bone
%     
%     % Extract image from figure axes
%     patchIm = getframe(gca);
%     patchIm = rgb2gray(patchIm.cdata);
%     
%     % Close open figures
%     close all
%     
%     % Write figure to file
%     imwrite( patchIm, ...
%         sprintf( fullfile([imFolder, '/', fileNameBase, '.tif']), t ), ...
%         'TIFF' );
%     
%     clear Options IV
%     
%     fprintf('Done\n');

end

%% Show 3D Texture Patch --------------------------------------------------

% Texture patch options
Options.PSize = 5;
Options.EdgeColor = 'none';

% Raw stack data
IV = xp.stack.image.apply();
IV = imadjustn(IV{1});

% Make a full screen image
figure('units', 'normalized', 'outerposition', [0 0 1 1] )

% Plot texture patch cut mesh in 2D
texture_patch_3d( mesh.f, mesh.v, ...
    mesh.f, mesh.v(:, [2 1 3]), IV, Options );

axis equal

colormap bone

clear Options IV