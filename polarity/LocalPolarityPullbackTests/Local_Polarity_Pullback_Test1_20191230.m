%% Local_Polarity_Pullback_Test1_12302019 =================================
% A script to prototype a pipeline that measures the tissue orientation
% polarity of the growing Drosopohila midgut by acting a Radon transform on
% pullbacks of dynamically sized patches of the 3D gut configuration.
% Orientation is then collated via parallel transport into a single frame.
%
% by Dillon Cislo and Noah Mitchell 12/30/2019
%==========================================================================

clear; close all; clc;

% Add necessary directories to path ---------------------------------------
addpath(genpath('/mnt/crunch/djcislo/MATLAB/euclidean_orbifolds'));
addpath(genpath('/mnt/data/code/gptoolbox'));
addpath(genpath('/mnt/data/code/gut_matlab/'));

%% Initialize ImSAnE Project ==============================================

% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.
dataDir = ['/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/' ...
    'Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/'];

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
    'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1_20200129' ];

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


%% ========================================================================
% LOAD MESH AND CREATE TILED PULLBACK IMAGES
%
% This is just regurgitated code from the axisymmetric pullback pipelines
% so I don't have to navigate the directories of the various saved
% pullbacks used in the active polarity pullback pipeline
%==========================================================================
%==========================================================================
%==========================================================================

t = xp.fileMeta.timePoints(1);

disp(['NOW PROCESSING TIME POINT ', num2str(t)]);

tidx = xp.tIdx(t);

%--------------------------------------------------------------------------
% Load Data for Current Time Point
%--------------------------------------------------------------------------

xp.loadTime(t);
xp.rescaleStackToUnitAspect();

%% Load the cylinder mesh
t = xp.fileMeta.timePoints(1);
mesh = read_ply_mod( sprintf( cylinderMeshBase, t ) );

% Consistently orient mesh faces
mesh.f = bfs_orient( mesh.f );

% Load the AD/PD vertex IDs
adIDx = h5read( dpFile, sprintf( ADBase, t ) );
pdIDx = h5read( dpFile, sprintf( PDBase, t ) );

ad3D = mesh.v( adIDx, : );
pd3D = mesh.v( pdIDx, : );

% Perform isotropic remeshing of the surfce if desired --------------------
[ ff, vv, ~, vvn ] = isotropic_remeshing( mesh.f, mesh.v, 15, 20, 0 );
mesh.f = ff; mesh.v = vv; mesh.vn = vvn;

% Clip the ears of the triangulation and update the AD/PD points if
% necessary ---------------------------------------------------------------
[ ff, ~, vv, newVIDx ] = clipMesh( mesh.f, mesh.v );
mesh.f = ff; mesh.v = vv; mesh.vn = mesh.vn(newVIDx, :);

adIDx = pointMatch( ad3D, mesh.v );
pdIDx = pointMatch( pd3D, mesh.v );

clear ff vv ad3D pd3D vvn newVIDx

% View results ------------------------------------------------------------
trisurf( triangulation( mesh.f, mesh.v ) );
hold on
scatter3( mesh.v(adIDx,1), mesh.v(adIDx,2), mesh.v(adIDx,3), ...
    'filled', 'r' );
scatter3( mesh.v(pdIDx,1), mesh.v(pdIDx,2), mesh.v(pdIDx,3), ...
    'filled', 'c' );
hold off
axis equal

%% ------------------------------------------------------------------------
% Create the Cut Mesh
%--------------------------------------------------------------------------

fprintf('Generating Cut Mesh... ');

cutOptions = struct();
cutMesh = ...
    cylinderCutMesh( mesh.f, mesh.v, mesh.vn, adIDx, pdIDx, cutOptions );

clear  cutOptions

fprintf('Done\n');

% View results ------------------------------------------------------------
P = cutMesh.pathPairs(:,1);

trisurf( triangulation( mesh.f, mesh.v ) );

hold on

line( mesh.v(P,1), mesh.v(P,2), mesh.v(P,3), ...
    'Color', 'c', 'LineWidth',2);

scatter3( mesh.v(adIDx,1), mesh.v(adIDx,2), mesh.v(adIDx,3), ...
    'filled', 'r' );
scatter3( mesh.v(pdIDx,1), mesh.v(pdIDx,2), mesh.v(pdIDx,3), ...
    'filled', 'm' );

hold off

axis equal

clear P

%% ------------------------------------------------------------------------
% Generate Pullback to Annular Orbifold Domain
%--------------------------------------------------------------------------

fprintf('Generating Pullback... ');

cutMesh = flattenAnnulus( cutMesh );

if tidx == 1
    a = minimizeIsoarealAffineEnergy( cutMesh.f, cutMesh.v, cutMesh.u );
end

cutMesh.u = [ a .* cutMesh.u(:,1), cutMesh.u(:,2) ];

fprintf('Done\n');

% View results ------------------------------------------------------------

patch( 'Faces', cutMesh.f, 'Vertices', cutMesh.u, ...
    'FaceVertexCData', cutMesh.v(:,3), 'FaceColor', 'interp', ...
    'EdgeColor', 'k' );

hold on

cornerColors = [ 1 0 0; 1 0 1; 0 1 1; 0 1 0 ];
corners = [ adIDx cutMesh.pathPairs(1,2), ...
    cutMesh.pathPairs(end,2) pdIDx ];

scatter( cutMesh.u( corners, 1 ), cutMesh.u( corners, 2 ), [], ...
    cornerColors, 'filled' );

hold off

axis equal

clear cornerColors corners

%% ========================================================================
% BEGIN RADON TRANSFORM PIPELINE
%
% This is where the actual new code begins
% =========================================================================
% =========================================================================
% =========================================================================

tesMesh = generateFPSSTessellation( mesh.f, mesh.v, ...
    'NumPoints', 1000, 'Overlap', 0.1 );

%% View Results -----------------------------------------------------------

% A color for each Voronoi submesh
cellColors = distinguishable_colors(length(tesMesh), [1 0 0]);

inCell = cell( size(mesh.f, 1), 1 );
faceColors = zeros( size(mesh.f, 1), 3 );
for fID = 1:size(mesh.f, 1)
    
    progressbar(fID, size(mesh.f,1));
    
    curInCell = [];
    for j = 1:length(tesMesh)
        
        if ismember( fID, tesMesh(j).FinFullMesh )
            curInCell = [ curInCell j ];
        end
        
    end
    
    inCell{fID} = curInCell;
    faceColors(fID, :) = mean(cellColors(curInCell, :), 1);
    
end

clear curInCell
%% 

cIDx = reshape( [ tesMesh.Center ], 2, length(tesMesh) ).';
cIDx = cIDx(:,1);

cXYZ = barycentricToCartesian( triangulation( mesh.f, mesh.v ), ...
    cIDx, ones( numel(cIDx), 3 ) ./ 3 );

patch('Faces', mesh.f, 'Vertices', mesh.v, ...
    'FaceVertexCData', faceColors, ...
    'FaceColor', 'flat', 'EdgeColor', 'none', ...
    'SpecularStrength', 0.1, 'DiffuseStrength', 0.1, ...
    'AmbientStrength', 0.8 );
hold on
scatter3( cXYZ(:,1), cXYZ(:,2), cXYZ(:,3), ...
    'filled', 'r' )
hold off
axis equal tight
camlight

%clear cIDx cXYZ
    

