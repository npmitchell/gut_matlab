%% Plot surfaces with intensity data
% This version allows a mesh refinement but doesn't use texture_patch 
% USE THE TEXTUREPATCH VERSION INSTEAD!
% Parameters
movies_no_frames = false


%% INITIALIZE ImSAnE PROJECT ==============================================
%
cd('/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/')
cd('Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/')

% We start by clearing the memory and closing all figures
clear; close all; clc;
addpath_recurse('/mnt/crunch/djcislo/MATLAB/CGAL_Code/')
addpath_recurse('/mnt/data/code/gptoolbox/')
addpath_recurse('/mnt/data/code/gut_matlab/')

% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.
dataDir    =  cd; 
projectDir = pwd ; 
% [ projectDir, ~, ~ ] = fileparts(matlab.desktop.editor.getActiveFilename); 
cd(projectDir);
xp = project.Experiment(projectDir, dataDir);
% A filename base template - to be used throughout this script
fn = 'Time_%06d_c1_stab';

fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fn, '.tif'];
fileMeta.nChannels          = 1;
fileMeta.timePoints         = 110:245;
fileMeta.stackResolution    = [.2619 .2619 .2619];
fileMeta.swapZT             = 1;

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
normal_step = 10;

% Define the surface detection parameters
channel = 1;
foreGroundChannel = 1;
ssfactor = 4;
niter = 25 ;
niter0 = 115 ;
ofn_ply = 'mesh_apical_ms_stab_' ; 
ofn_ls = 'msls_apical_stab_' ;
ofn_smoothply = 'mesh_apical_stab_' ;
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
layerDistance = 5 ;  % layerDistance is in pix
sigma = 10 ;  % Sigma smooths
makeIP = 'MIP' ;  % SIP, MIP are options for makeIP
IPonly = false ;
onionOpts = struct('nLayers', nLayers, 'layerDistance', layerDistance,...
                   'sigma', sigma, 'makeIP', makeIP, 'IPonly', IPonly);

%% Go through all smoothed meshes to get XYZLimits
% LOAD MESH FROM FILES ===================================================
% Cell arrays to hold the mesh struct attributes
v = cell(length(xp.fileMeta.timePoints),1);
f = cell(length(xp.fileMeta.timePoints),1);
vn = cell(length(xp.fileMeta.timePoints),1);
msls_axis_order = 'xyzc' ;
normal_step = 2 ;

for t = xp.fileMeta.timePoints
    disp(['loading mesh for t= ', num2str(t)])
    % Convert into timepoint ID
    tidx = xp.tIdx(t);
    mesh_outfn = [ofn_smoothply, ...
        num2str(xp.fileMeta.timePoints(tidx), '%06d'), '.ply'];
    outputMesh = fullfile(mslsDir, mesh_outfn);
    
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
    if t == xp.fileMeta.timePoints(1)
        xmin = min(mesh.v(:,1)) ;
        ymin = min(mesh.v(:,2)) ;
        zmin = min(mesh.v(:,3)) ;
        xmax = max(mesh.v(:,1)) ;
        ymax = max(mesh.v(:,2)) ;
        zmax = max(mesh.v(:,3)) ;
    else
        xmin = min(min(mesh.v(:,1)), xmin) ;
        ymin = min(min(mesh.v(:,2)), ymin) ;
        zmin = min(min(mesh.v(:,3)), zmin) ;
        xmax = max(max(mesh.v(:,1)), xmax) ;
        ymax = max(max(mesh.v(:,2)), ymax) ;
        zmax = max(max(mesh.v(:,3)), zmax) ;
    end
end
% Set the limits
xLimits = [xmin xmax] ;
yLimits = [ymin ymax] ;
zLimits = [zmin zmax] ;

%% Formatting for video
fig = figure('units','normalized','outerposition',[0 0 1 1]);
ax = gca;
ax.NextPlot = 'replaceChildren';
if movies_no_frames
    % create videos
    F1(length(xp.fileMeta.timePoints)) = struct('cdata', [], 'colormap', [] );
    F2(length(xp.fileMeta.timePoints)) = struct('cdata', [], 'colormap', [] );
    F3(length(xp.fileMeta.timePoints)) = struct('cdata', [], 'colormap', [] );
    F4(length(xp.fileMeta.timePoints)) = struct('cdata', [], 'colormap', [] );
else
    % create frames
    outfigdir = fullfile(mslsDir, ') ;
end

% Refine the mesh
meshfn_master = 'mesh_apical_%06d.ply' ;
prec = '1p00' ;
mlxprogram = ['/mnt/data/code/meshlab_codes/refineLS3Loop_' prec 'world.mlx'] ;


for T = 1:length(xp.fileMeta.timePoints)
    disp(['capturing tp = ' num2str(T)])
    
    % create patch object
    meshfn = sprintf(meshfn_master, xp.fileMeta.timePoints(T)) ;
    meshfn = fullfile(mslsDir, meshfn) ;
    rmesh = read_ply_mod(meshfn) ;
    clear mm
    mm.faces = rmesh.f ;
    mm.vertices = rmesh.v ;

    % Get the PLY
    rmeshfn = [meshfn(1:end-4) '_refine' prec 'wu.ply'] ;
    if ~exist(rmeshfn, 'file')
        % using Meshlab script instead of refinepatch()
        command = ['meshlabserver -i ' meshfn ' -o ' rmeshfn, ...
            ' -s ' mlxprogram ' -om vn'];
        % Either copy the command to the clipboard
        clipboard('copy', command);
        % or else run it on the system
        disp(['running ' command])
        system(command)
    end
    rmesh = read_ply_mod(rmeshfn) ;

    % LOAD THE FIRST TIME POINT ===========================================
    % obtain the image from xp
    xp.loadTime(xp.fileMeta.timePoints(T));
    xp.rescaleStackToUnitAspect();

    tmp = xp.stack.image.apply();
    im = tmp{1} ;
    clip = 65536 ;
    im = im * 3 ;
    im(im > clip) = clip ;
    % im = im / clip * 65536 ;

    % Obtain the mesh indices for pt correspondence with the image data
    rv = round(rmesh.v) ;
    rvInd = sub2ind( size(im), rv(:,2), rv(:,3), rv(:,1) );
    colors = im(rvInd) ;

    % colors is Nx3 RGB or Nx1 intensity array
    patch( 'Faces', rmesh.f, 'Vertices', rmesh.v, ...
           'FaceVertexCData', colors, ...
           'FaceColor', 'interp', 'EdgeColor', 'interp' );
    colormap(gray);
    axis equal
    box off
    set(gca,'Color','k')
    set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
    set(gca,'XtickLabel',[], 'YtickLabel', [], 'ZtickLabel', []);
    xlim(xLimits); ylim(yLimits); zlim(zLimits);
    
    % Now set the camera angle
    % load translation + rotation
    % set camera angle using view()
    view()
    
    if movie_no_frames
        % capture the entire figure as a frame
        drawnow
        F1(T) = getframe(fig);
        cla
    else
        saveas(fig, 
        cla
        
    end
end

close(fig);

%% Write to Video File ====================================================
if movie_no_frames 
    v = VideoWriter(fullfile(mslsDir, 'view01.avi'));
    v.Quality = 100;
    v.FrameRate = 10;
    open(v);

    for T = 1:length(xp.fileMeta.timePoints)

        writeVideo(v, F1(T));
        disp(['wrote frame ' num2str(T) ])

    end

    close(v);
end


% patchT from 2d image to a patch object
