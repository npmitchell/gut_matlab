% MASK3D crops out a thin layer around a given mesh, +/- some specified
% thickness, from a 3d volume and plots it. 
% Note that an alternative approach is to use the level sets, dilate and
% erode them, then multiply that by the data. 
%
% Makes use of the SOImask3D.m function of imsane.

%% INITIALIZE ImSAnE PROJECT ==============================================
%
% We start by clearing the memory and closing all figures
clear; close all; clc;
addpath_recurse('/mnt/crunch/djcislo/MATLAB/CGAL_Code/')
addpath_recurse('/mnt/data/code/gptoolbox/')

% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.
dataDir    =  cd; 
[ projectDir, ~, ~ ] = fileparts(matlab.desktop.editor.getActiveFilename); 
cd(projectDir);

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
fn = 'Time_%06d_c1_stab';

fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fn, '.tif'];
fileMeta.nChannels          = 1;
fileMeta.timePoints         = 110:190;
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

%%%%%%%%%%%%%%%%%%%%%%%%%
% Other options
%%%%%%%%%%%%%%%%%%%%%%%%%
mlxprogram = 'surface_rm_resample20k_reconstruct_LS3_1p2pc_ssfactor4.mlx';
msls_axis_order = 'yxzc';
% Mesh marching options
normal_step = 12;

% Define the surface detection parameters
channel = 1;
foreGroundChannel = 1;
ssfactor = 4;
niter = 25 ;
niter0 = 25 ;
ofn_ply = 'mesh_apical_ms_stab_' ; 
ofn_ls = 'msls_apical_stab_' ;
ms_scriptDir = '/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/' ;
lambda1 = 1 ;
lambda2 = 1 ;
exit_thres = 0.000001 ;
smoothing = 0.1 ;
nu = 0.0 ;
pre_nu = -5 ;
pre_smoothing = 1 ;
post_nu = 2;
post_smoothing = 4 ;

% Name the output mesh directory ------------------------------------------
msls_exten = ['_prnu' strrep(strrep(num2str(pre_nu, '%d'), '.', 'p'), '-', 'n')];
msls_exten = [msls_exten '_prs' strrep(num2str(pre_smoothing, '%d'), '.', 'p') ];
msls_exten = [msls_exten '_nu' strrep(num2str(nu, '%0.2f'), '.', 'p') ];
msls_exten = [msls_exten '_s' strrep(num2str(smoothing, '%0.2f'), '.', 'p') ];
msls_exten = [msls_exten '_pn' num2str(post_nu, '%d') '_ps',...
    num2str(post_smoothing)];
msls_exten = [msls_exten '_l' num2str(lambda1) '_l' num2str(lambda2) ];
if projectDir(end) ~= '/'
    projectDir = [projectDir '/'];
end
mslsDir = [projectDir 'msls_output'];
mslsDir = [mslsDir msls_exten '/'] ;

% The dimension to use to grab extremal seeds
seeddim = 3;

% Onion Options
nLayers = 75 ;  % nLayers must be an odd int
layerDistance = 20 ;  % layerDistance is in pix
sigma = 10 ;  % Sigma smooths
makeIP = 'MIP' ;  % SIP, MIP are options for makeIP
IPonly = false ;
onionOpts = struct('nLayers', nLayers, 'layerDistance', layerDistance,...
                   'sigma', sigma, 'makeIP', makeIP, 'IPonly', IPonly);

%% LOAD THE FIRST TIME POINT ==============================================
xp.loadTime(xp.fileMeta.timePoints(first_tp));
xp.rescaleStackToUnitAspect();

%% LOAD MESH FROM FILES ===================================================

% DEFINE MESHLAB SCRIPT FOR CLEANING UP THE PLY ==========================
% add to path to meshlabserver if needed by making a symbolic link
% sudo ln -s /Applications/meshlab.app/Contents/MacOS/meshlabserver /usr/bin/meshlabserver
msls_mesh_outfn = [ofn_ply, ...
    num2str(fileMeta.timePoints(first_tp), '%06d') '.ply'];
PCfile = fullfile(mslsDir, msls_mesh_outfn);
mesh_outfn = ['mesh_apical_', num2str(fileMeta.timePoints(first_tp), '%06d'), '.ply'];
outputMesh = fullfile(mslsDir, mesh_outfn);
meshlabScript = fullfile(projectDir, mlxprogram);

% LOAD MESH FROM FILES ===================================================
% Cell arrays to hold the mesh struct attributes
v = cell(length(xp.fileMeta.timePoints),1);
f = cell(length(xp.fileMeta.timePoints),1);
vn = cell(length(xp.fileMeta.timePoints),1);

for t = xp.fileMeta.timePoints
    % Convert into timepoint ID
    tidx = xp.tIdx(t);
    
    % Clean up mesh file for this timepoint using MeshLab -----------------
    msls_mesh_outfn = [ ofn_ply, ...
        num2str( xp.fileMeta.timePoints(tidx), '%06d' ), '.ply' ];
    PCfile = fullfile( mslsDir, msls_mesh_outfn );
    mesh_outfn = ['mesh_apical_', ...
        num2str(xp.fileMeta.timePoints(tidx), '%06d'), '.ply'];
    outputMesh = fullfile(mslsDir, mesh_outfn);
    
    if ~exist( outputMesh, 'file')
        meshlabScript = fullfile(projectDir, mlxprogram);
        command = ['meshlabserver -i ' PCfile ' -o ' outputMesh, ...
                   ' -s ' meshlabScript ' -om vn'];
        % Either copy the command to the clipboard
        clipboard('copy', command);
        % or else run it on the system
        disp(['running ' command])
        system(command) 
    else
        disp(['t=', num2str(tidx) ': smoothed mesh file found, loading...'])
    end
      
    % Read in the mesh file -----------------------------------------------
    mesh = read_ply_mod( outputMesh );
    
    if strcmp( msls_axis_order, 'yxzc' )
        % Multiply by ssfactor to return to original scale
        mesh.v = mesh.v( :, [3,2,1] ) ;
        mesh.vn = mesh.vn( :, [3,2,1] ) ;
    end
    
    % Make sure vertex normals are normalized
    mesh.vn = mesh.vn ./ sqrt( sum( mesh.vn.^2, 2 ) );
    
    % Normally evolve vertices
    mesh.v = mesh.v + normal_step .* mesh.vn;
    
    % Re-orient faces
    mesh.f = reorient_facets( mesh.v, mesh.f );
    
    v{tidx} = mesh.v;
    f{tidx} = mesh.f;
    vn{tidx} = mesh.vn;
    
end

% Create a stack of all meshes
meshStack = struct( 'v', v, 'f', f, 'vn', vn );

clear mesh_outfn outputMesh mesh
clear v f vn ssfactor
disp('done loading meshStack')

%% Now make each mask =====================================================
for ii = xp.fileMeta.timePoints
    mesh = {'v', meshStack{ii}} ;
    % SOImask = SOImask3D(mesh, onionOpts, maskSize)
    % SOImask:  binary mask
    % mesh:         structure with field f, v, vn for faces, vertices,
    %               vertex normals
    % onionOpts:    onion options as provided to
    %               SurfaceOfInterest.pullbackStack
    % maskSize:     required dimensions of mask
    % requires VOXELIZE from matlab file exchange
    SOImask = SOImask3D(mesh, onionOpts, maskSize) ;
    size(SOImask)
    
    % Apply to the data
    xp.loadTime(xp.fileMeta.timePoints(first_tp));
    xp.rescaleStackToUnitAspect();
    tp_masked = SOImask .* xp.stack.apply() ;
end


