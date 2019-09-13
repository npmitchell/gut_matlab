%% GUT_PIPELINE_20190415: GUT CONFORMAL MAPPING PIPELINE
% by Dillon Cislo & NPMitchell
%
% This is a pipeline to take the surface of the growing Drosophila gut and
% conformally map patches of it to the unit disk

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
normal_step = 10;

% Define the surface detection parameters
channel = 1;
foreGroundChannel = 1;
ssfactor = 4;
niter = 25 ;
niter0 = 115 ;
ofn_ply = 'mesh_apical_ms_' ; 
ofn_ls = 'msls_apical_' ;
ms_scriptDir = '/mnt/data/code/gut_python/' ;
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
% SOI saving options
imwriteOptions = {'tif'};
soiDir = fullfile(projectDir, ['gut_apical_conformal_msls' msls_exten]);
soiDir = [soiDir '_' num2str(nLayers) 'layer_voronoipatches/']; 
soi_save_options = struct('dir',soiDir,'imwriteOptions',{imwriteOptions},...
                    'make8bit',false);

%% LOAD THE FIRST TIME POINT ==============================================
xp.loadTime(xp.fileMeta.timePoints(first_tp));
xp.rescaleStackToUnitAspect();

%% DETECT THE SURFACE =====================================================
% Surface detection parameters --------------------------------------------
detectOptions = struct( 'channel', channel, 'ssfactor', ssfactor, ...
    'niter', niter, ...
    'niter0', niter0, ...
    'lambda1', lambda1, ...
    'lambda2', lambda2, ...
    'nu', nu, ...
    'smoothing', smoothing, ...
    'post_nu', post_nu, ...
    'post_smoothing', post_smoothing, ...
    'exit_thres', exit_thres, ...
    'foreGroundChannel', foreGroundChannel, ...
    'fileName', sprintf( fn, xp.currentTime ), ...
    'mslsDir', mslsDir, 'ofn_ls', ofn_ls, 'ofn_ply', ofn_ply,...
    'ms_scriptDir', ms_scriptDir, ...
    'timepoint', xp.currentTime, ...
    'zdim', 2, ...
    'pre_nu', pre_nu,...
    'pre_smoothing', pre_smoothing,...
    'ofn_smoothply', 'mesh_apical_',...
    'mlxprogram', './surface_rm_resample20k_reconstruct_LS3_1p2pc_ssfactor4.mlx',...
    'init_ls_fn', 'none', ...
    'run_full_dataset', false, ...
    'radius_guess', -1, ...
    'dset_name', 'exported_data');

% Set detect options ------------------------------------------------------
xp.setDetectOptions( detectOptions );

% clear msls_exten imwriteOptions saveDir
% clear channel foreGroundChannel
% clear niter niter0 lambda1 lambda2
% clear exit_thres smoothing nu
% clear post_nu post_smoothing

%% CREATE THE SUBSAMPLED H5 FILE FOR INPUT TO ILASTIK =====================
% skip if already done
if ~exist([projectDir sprintf(fn, xp.currentTime) '.h5'], 'file')
    xp.detector.prepareIlastik(xp.stack);
    disp('done outputting downsampled data h5 for surface detection')
else
    disp('h5 was already output, skipping...')
end
disp('Open with ilastik if not already done')

%% pre processing: Make subsampled h5 files for ilastik if not already done
open_fbar = true ;
if xp.expMeta.dynamicSurface == 1
    disp('preprocessing for ilastik...')
    for t = 1:length(xp.fileMeta.timePoints)  
        time = xp.fileMeta.timePoints(t) ;
        % load the data
        if ~exist(sprintf([fn '.h5'], time), 'file') 
            disp(['preprocessing for ilastik: t=', num2str(time)])
            xp.loadTime(xp.fileMeta.timePoints(t));
            xp.rescaleStackToUnitAspect();
            if xp.expMeta.dynamicSurface 
                detectOptions.fileName = sprintf(fn, xp.currentTime) ;
                xp.setDetectOptions(detectOptions);  
                xp.detector.prepareIlastik(xp.stack);
            end
        else
            if open_fbar 
                tmp = sprintf([fn '.h5'], time) ;
                msg = strrep(['Already prepared ' tmp], '_', '\_') ;
                fbar = waitbar(t/length(xp.fileMeta.timePoints), msg) ;
                open_fbar = false ;
            else
                tmp = sprintf([fn '.h5'], time) ;
                msg = strrep(['Already prepared ' tmp], '_', '\_') ;
                waitbar(t / length(xp.fileMeta.timePoints), fbar, msg) 
            end
        end
    end
    close(fbar)
    disp('done')
end

%% TRAIN DATA IN ILASTIK TO IDENTIFY APICAL/YOLK ==========================
% open ilastik, train until probabilities and uncertainty are satisfactory

%% Create MorphoSnakesLevelSet from the Probabilities from ilastik ========
xp.detectSurface();
fileMeta = xp.fileMeta ;

%% Morphosnakes for all remaining timepoints ==============================
% % So do below: The above should do the following:
% for tp = fileMeta.timePoints(2:end)
%     xp.currentTime = tp ;
%     detectOptions.timepoint = xp.currentTime ;
%     detectOptions.fileName = sprintf( fn, xp.currentTime );
% end

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

%% SET ATLAS PARAMETERS FOR THE FIRST TIME POINT ==========================

% Load data for the first time point --------------------------------------
t = xp.fileMeta.timePoints(1); tidx = xp.tIdx(t);
xp.loadTime(t); xp.rescaleStackToUnitAspect();

% mesh object holds the whole mesh (not cut up into patches)
mesh = meshStack(tidx);

%% Plot to Define the seed points for the fitter --------------------------
ptsz = 10 ;
close all
s = scatter(mesh.v(:, 1), mesh.v(:, 2), ptsz, mesh.v(:, 2), 'filled');
alpha(s, 0.1)
xlabel('x'); ylabel('y')
axis('equal')
figure()
s = scatter(mesh.v(:, 1), mesh.v(:, 3), ptsz, mesh.v(:, 2), 'filled') ;
alpha(s, 0.1)
xlabel('x'); ylabel('z')

%% now define seeds
close all
% VorSeedsXinit = [[400, 700, 250]; [250, 700, 530]] ;
% num_vorpatches = size(VorSeedsXinit, 1);

% Alternately, grab seeds to define voronoi regions
% vorseeds = [0, 0] ;
% vorseeds(1) = find( mesh.v(:, seeddim) == min(mesh.v(:, seeddim)) );
% vorseeds(2) = find( mesh.v(:, seeddim) == max(mesh.v(:, seeddim)) );
% VorSeedsXinit = mesh.v( vorseeds, : );
VorSeedsXinit = [];
num_vorpatches = size(VorSeedsXinit, 1) ;

% Seeds for patch regions
% % First fold DV
seeds = zeros(2,3);
seeds(1, :) = [100, 1000, 300] ; 
seeds(2, :) = [550, 1000, 300] ; 
% seeds(1, :) = [370, 700, 100] ;
% seeds(2, :) = [250, 700, 450] ;
% % % Front fold DV
% seeds(3, :) = [400, 900, 50] ;
% seeds(4, :) = [250, 900, 470] ;
% % % Back fold DV
% seeds(5, :) = [270, 500, 200] ;
% seeds(6, :) = [380, 500, 480] ;
% % First fold Lat
% seeds(7, :) = [100, 550, 450];
% seeds(8, :) = [520, 550, 400];
% % Second fold Lat
% seeds(9, :) = [50, 350, 400];
% seeds(10, :) = [550, 350, 300];
% % Third fold Lat
% seeds(11, :) = [100, 830, 450];
% seeds(12, :) = [550, 830, 450];
diskSeedsXinit = seeds;
num_diskpatches = length(seeds) ;

diskSeeds = pointMatch(diskSeedsXinit, mesh.v);
VorSeeds = pointMatch(VorSeedsXinit, mesh.v);

fitOptions = struct('VorSeeds', VorSeeds, 'transitionWidth', 20,...
                    'diskSeeds', diskSeeds, 'diskRadius', 70, ...
                    'makeTMaps', false);
                
xp.setFitOptions(fitOptions);
disp('done defining seeds')

%% Create atlas ------------------------------------------------------------
% Hold off from generating SOI here.
% Creates the submeshes, indexes vertices of the larger mesh
xp.fitSurface(mesh);

% Keep track of seed propagation for later time points --------------------
VorSeedsX = cell( length(xp.fileMeta.timePoints), 1 );
diskSeedsX = cell( length(xp.fileMeta.timePoints), 1 );

VorSeedsX{1} = mesh.v(VorSeeds, :);
diskSeedsX{1} = mesh.v(diskSeeds,: );

% clear seeds diskSeeds VorSeeds diskSeedsXinit VorSeedsXinit

% The total number of charts per time point
nCharts = numel(xp.fitter.fittedParam.submeshes);

% A cell array to hold the 3D submeshes created during the fitting process
% Indices are timepoint, submesh: ie allMeshes(1,3) is timepoint 1, 
% submesh 3. Returns a struct: struct( 'v', v, 'f', f, 'vn', vn, 'u', empty_for_now? );
allMeshes = cell( length(xp.fileMeta.timePoints), nCharts );
disp('done creating atlas')

%% Check it
meshcheckDir = [soiDir 'submeshes/'] ;
if ~exist(meshcheckDir, 'dir')
    mkdir(meshcheckDir)
end
outfn = fullfile(meshcheckDir, ['submeshes' num2str(t, '%04d') '.png']) ;
% Inspect the whole mesh first
% xp.fitter.inspectMesh ;
% trisurf(triangulation(meshStack(tidx).f, meshStack(tidx).v), 'FaceColor', 'k', 'EdgeColor', 'k')  ;
trisurf(triangulation(meshStack(tidx).f, meshStack(tidx).v), ...
    'FaceVertexCData', zeros(size(meshStack(tidx).v,1),3), ...
    'FaceColor', 'flat', 'EdgeColor', 'k')  ;
% hold on
% trisurf(triangulation(allMeshes{tidx}.f, allMeshes{tidx}.v), ...
%     'FaceVertexCData', ones(size(allMeshes{tidx}.v,1),3) .* [1, 0, 0], ...
%     'FaceColor', 'flat', 'EdgeColor', 'r')  ;

xp.fitter.inspectMesh(1:nCharts, true, true) ;
title(['Submeshes for t=' num2str(t)])
view([180, 0])
saveas(gcf, outfn) ;
% xp.fitter.inspectMesh(1:6);
% figure()
% xp.fitter.inspectMesh(7:12);

%% Populate allMeshes =====================================================
% Point in domain of parametrization to hold fixed next round
fixedPtU = cell( nCharts, 1 );  
% ---> would be used in spectral conformal projection algorithm
fixedPtX = cell( length(xp.fileMeta.timePoints), 1 );
tmpX = cell( nCharts, 1 );

for k = 1:nCharts
    
    subm = xp.fitter.fittedParam.submeshes{k};
    
    % fixedPtUtmp = sqrt(mean(subm.u{1}.^2)/2);
    % fixedPtIdx = pointMatch(fixedPtUtmp, subm.u{1});
    % tmpX{k} = subm.v(fixedPtIdx,:);
    % fixedPtU{k} = subm.u{1}(fixedPtIdx,:);
    
    allMeshes{ 1, k } = subm;
    
end

% fixedPtX{tidx} = tmpX;
clear fixedPtUtmp fixedPtIdx subm tmpX
disp('done with generating allMeshes for first timepoint')

%% LOOP FOR DYNAMIC ATLAS CREATION ========================================
% We build an atlas by dividing the mesh into overlapping submeshes and
% mapping those submeshes into the plane.

for t = xp.fileMeta.timePoints(2:end)
    
    disp(['Now processing time point ', num2str(t)]);
    
    tidx = xp.tIdx(t);
    
    % Load data -----------------------------------------------------------
    xp.loadTime(t);
    xp.rescaleStackToUnitAspect();
    
    mesh = meshStack(tidx);
    
    % Initialize fitter ---------------------------------------------------
    fitOptions.VorSeeds = pointMatch(VorSeedsX{tidx-1}, mesh.v);
    fitOptions.diskSeeds = pointMatch(diskSeedsX{tidx-1}, mesh.v);
    % fitOptions.fixedPtX = fixedPtX{tidx-1};
    % fitOptions.fixedPtU = fixedPtU;
    
    xp.setFitOptions(fitOptions);
    xp.fitSurface(mesh); 
    
    % Find and propagate seeds points for the next time point -------------
    VorSeedsX{tidx} = mesh.v( fitOptions.VorSeeds, : );
    diskSeedsX{tidx} = mesh.v( fitOptions.diskSeeds, : );
    
    % tmpX = cell( nCharts, 1 );
    for k = 1:nCharts
        
        subm = xp.fitter.fittedParam.submeshes{k};
        % V = subm.v;
        % fixedPtIdx = pointMatch(fixedPtX{k}, V);         
        % tmpX{k} = V(fixedPtIdx,:);
        
        allMeshes{tidx, k} = subm;
        
    end
    
    % Save image of the submeshes
    close all
    outfn = fullfile(meshcheckDir, ['submeshes' num2str(t, '%04d') '.png']) ;
    xp.fitter.inspectMesh(1:nCharts) ;
    title(['Submeshes for t=' num2str(t)])
    view([180, 0])
    saveas(gcf, outfn) ;

    % If there were a fixed point inside disk, define here, but no need if 
    % boundary if fixed in Dirichlet energy minimization
    % fixedPtX{tidx} = tmpX;
    
end
clear mesh tmpX V fixedPtIdx subm
disp('done with dynamic atlas creation')

%% SURFACE PARAMETERIZATION (WARNING: SLOW) ===============================

SMArr2D = cell( length(xp.fileMeta.timePoints), nCharts );

% For conformal mapping to the disk (default)
param = struct();

for t = xp.fileMeta.timePoints
    
    tidx = xp.tIdx(t);
    disp(['Now processing time point ', num2str(t)]);
    
    for i = 1:nCharts
        
        disp(['Now processing submesh ', num2str(i)]);
        
        V3D = allMeshes{ tidx, i }.v;
        F = allMeshes{ tidx, i }.f;
        
        % If param is empty struct, use conformal mapping to disk using
        % Dirichlet energy minimization. The image has diameter 1 at (0.5, 0.5)
        [ v2D, ~ ] = surface_parameterization( F, V3D, param );
        
        % Translate/scale output to the unit disk (diameter 2 at (0,0))
        SMArr2D{ tidx, i } = 2.*(v2D - 0.5);
        
    end
    
end

clear v2D F V3D param;
disp('done with surface parameterization')

%% CHECK FOR DUPLICATE POINTS =============================================

for t = xp.fileMeta.timePoints
    tidx = xp.tIdx(t);
    
    for i = 1:nCharts
        % If two vertices mapped to same (u,v) coords
        [dupl, ~] = find_duplicate_rows(SMArr2D{tidx,i});
        if ~isempty(dupl)
            disp(['Submesh ', num2str(i), ...
                ' of timepoint ', num2str(t), ...
                ' contains duplicate points']);
        end
        
    end
    
end

clear dupl
disp('done checking for duplicates')

%% CHECK FOR MESH SELF-INTERSECTIONS ======================================

%**************************************************************************
% *************************************************************************
% SKIP THIS FOR NOW - GETTING SEG-FAULTS AND HAVE TO CHECK WHY
% *************************************************************************
%**************************************************************************
% 
% for t = xp.fileMeta.timePoints
%     tidx = xp.tIdx(t);
%     
%     for i = 1:nCharts
%         
%         F = allMeshes{ tidx, i };
%         VV = SMArr2D{tidx, i};
%         
%         [badMesh, ~] = mesh_self_intersection_3d( F, ...
%             [ VV, zeros(size(VV,1),1) ] );
%         if badMesh
%             disp(['Submesh ', num2str(i), ...
%                 ' of timepoint ', num2str(t), ...
%                 ' contains self-intersections!']);
%         end
%         
%     end
%     
% end
% 
% clear F VV badMesh

%% CLIP BOUNDARY POINTS TO THE UNIT DISK ==================================

for t = xp.fileMeta.timePoints
    
    tidx = xp.tIdx(t);
    
    for i = 1:nCharts
        
        v2D = SMArr2D{tidx,i};
        v2DC = complex( v2D(:,1), v2D(:,2) );
        
        outIDx = ( abs( v2DC ) > 1 );
        
        v2DC( outIDx ) = exp( 1i .* angle( v2DC( outIDx ) ) );
        
        SMArr2D{ tidx, i } = [ real(v2DC), imag(v2DC) ];
        
    end
    
end

clear v2D v2DC outIDx
disp('done clipping')

%% CORRECT FOR POSSIBLE REFLECTIONS =======================================
% This fixes possible reflections
for t = xp.fileMeta.timePoints
    tidx = xp.tIdx(t);
    
    for i = 1:nCharts

        % The face unit normals of the planar triangulation
        fN = faceNormal( triangulation( allMeshes{tidx,i}, ...
            SMArr2D{tidx,i} ) );
        
        % The sign of the z-coordinate of the unit normals
        pm = sign( fN(:,3) );
        
        if ( numel( unique( pm ) ) > 1 )
            
            % This will only be true if the faces are not consistently
            % ordered which shouldn't be possible
            warning( ['Invalid face ordering in submesh ', num2str(i), ...
                ' of timepoint ', num2str(t) ] );
            
        elseif ( unique( pm ) < 0 )
            
            % Correct the reflection if the face unit normal is pointing in
            % the negative z-direction
            disp( ['Correcting reflection in submesh ', num2str(i), ...
                ' of timepoint ', num2str(t) ]) ;
            
            uv = SMArr2D{tidx, i};
            SMArr2D{tidx, i} = [ uv(:,1), -uv(:,2) ];
           
        end
        
    end
    
end

clear pm uv fN
disp('done correcting for reflections')

%% FIND THE 'ISOAREAL' MOBIUS MAPPING OF THE FIRST TIME POINT =============

for i = 1:nCharts
    
    disp(['Now processing chart number ', num2str(i), ...
        ' of the first time point']);
    
    w = SMArr2D{1,i};
    
    a = minimizeIsoarealMobiusEnergy( allMeshes{1,i}.f, ...
        allMeshes{1,i}.v, w );
    a = complex( a(1), a(2) );
    
    w = complex( w(:,1), w(:,2) );
    w = ( w - a ) ./ ( 1 - conj(a) .* w );
    
    SMArr2D{1,i} = [ real(w), imag(w) ];
    
end

clear w a
disp('done with isoareal mobius mapping')
% todo: populateSOI here based on first TP, repopulate in the loop below.

%% Visualize the results of Mapping ---------------------------------------

% newFaces = cell(length(xp.fileMeta.timePoints),1);
% i = 1;
% for t = xp.fileMeta.timePoints
%     tidx = xp.tIdx(t);
%     newFaces{tidx} = allMeshes{tidx,i}.f;
% end
% GraphPlanarMesh4D( SMArr2D, newFaces, [], xp.fileMeta );
% clear newFaces;

%% CONFORMALLY ALIGN SURFACES (WARNING: SLOW) =============================
% Choose subset of vertices at time t and point match them in 3D to another 
% subset at time t+1. Find mobius transf of disk that most closely aligns 
% the (u, v) coords of the paired points.
 
for t = xp.fileMeta.timePoints(2:end)
    
    disp(['Now aligning time point ', num2str(t)]);
    
    tidx = xp.tIdx(t);
    
    for i = 1:nCharts
        
        disp(['Now aligning chart ', num2str(i)]);
        
        disp('Generating point sets');
        [M, B, eps] = generate_fpss( ...
            allMeshes{tidx,i}.f, allMeshes{tidx,i}.v, SMArr2D{tidx,i}, ...
            allMeshes{tidx-1,i}.f, allMeshes{tidx-1,i}.v, SMArr2D{tidx-1,i}, ...
            'NumPnts', 100, 'Method', 'Random');
        
        disp('Finding optimal Mobius transformation');
        [z, theta] = optimal_mobius_search( M, B, ...
            'Radii', eps, 'Display', false);
        
        disp('Applying Mobius transformation')
        w = SMArr2D{tidx, i}; w = complex(w(:,1), w(:,2));
        
        w = exp(1i * theta) .* ( (w - z) ./ (1 - conj(z) .* w ) );
        
        SMArr2D{tidx, i} = [real(w), imag(w)];
        
    end
    
    
end

clear M B eps z theta w;
disp('done conformally aligning surfaces')

%% SET ATLAS PARAMETERS FOR THE FIRST TIME POINT ==========================

% Load data for the first time point --------------------------------------
t = xp.fileMeta.timePoints(1); tidx = xp.tIdx(t);
xp.loadTime(t); xp.rescaleStackToUnitAspect();

mesh = meshStack(tidx);

% Initialize fitter -------------------------------------------------------
fitOptions.VorSeeds = pointMatch(VorSeedsX{tidx}, mesh.v);
fitOptions.diskSeeds = pointMatch(diskSeedsX{tidx}, mesh.v);
fitOptions.fixedPtX = fixedPtX{tidx};
fitOptions.fixedPtU = fixedPtU;
                
xp.setFitOptions(fitOptions);
xp.fitSurface(mesh);

% Set the pullbacks of the submeshes to the plane
for i = 1:nCharts
    
    % Check that the generated submeshes match previous results
    % (Unnecessary step - figure out a better way to do this)
    assert( isequal( xp.fitter.fittedParam.submeshes{i}.v, ...
        allMeshes{tidx,i}.v ) && ...
        isequal( xp.fitter.fittedParam.submeshes{i}.f, ...
        allMeshes{tidx,i}.f ), ...
        ['Abnormal result for submesh ', num2str(i), ...
        ' of time point ', num2str(t)] );
    
    % This line sets the uv coordinates for the mesh
    xp.fitter.setPullBack(SMArr2D{tidx, i}, i);
    
end
disp('done setting atlas params for first tp')

%% GENERATE THE ATLAS =====================================================
% New behavior for ImSAnE: supply PopulateSOI. This line generates an empty
% SOI with all attributes other than stuff done by populateSOI <- check what? 
xp.generateSOI('PopulateSOI', false);

%% POPULATE THE SOI =======================================================
% New behavior for ImSAnE: repopulateSOI uses supplied mesh 
% in xp.fitter.fittedParam to generate pullback mapping at current time.
xp.fitter.repopulateSOI(xp.SOI, xp.currentTime);

% Pullback the stack to the desired charts                
xp.SOI.pullbackStack(xp.stack, [], xp.currentTime, onionOpts);

% Set the allFit object to store the information about the current fitter
% Now, xp.allFit is a struct attribute containing the fitter for each time 
% point, and we set it for the first time point now.
xp.setAllFit();
disp('done populating the SOI for the first timepoint')

%%  LOOP FOR DYNAMIC ATLAS GENERATION (WARNING: SLOW) =====================
% NOTE: come here only after generating the SOI for the first time point!!!

for t = xp.fileMeta.timePoints(2:end)
    
    disp(['Now processing time point ', num2str(t)]);
    
    tidx = xp.tIdx(t);
    
    % Load data -----------------------------------------------------------
    xp.loadTime(t);
    xp.rescaleStackToUnitAspect();
    
    mesh = meshStack(tidx);
    
    % Initialize fitter ---------------------------------------------------
    fitOptions.VorSeeds = pointMatch(VorSeedsX{tidx-1}, mesh.v);
    fitOptions.diskSeeds = pointMatch(diskSeedsX{tidx-1}, mesh.v);
    fitOptions.fixedPtX = fixedPtX{tidx-1};
    fitOptions.fixedPtU = fixedPtU;
    
    xp.setFitOptions(fitOptions);
    xp.fitSurface(mesh); 
    
    % Set the pullbacks of the submeshes to the plane
    for i = 1:nCharts
        
        % Check that the generated submeshes match previous results
        % (Unnecessary step - figure out a better way to do this)
        assert( isequal( xp.fitter.fittedParam.submeshes{i}.v, ...
            allMeshes{tidx,i}.v ) && ...
            isequal( xp.fitter.fittedParam.submeshes{i}.f, ...
            allMeshes{tidx,i}.f ), ...
            ['Abnormal result for submesh ', num2str(i), ...
            ' of time point ', num2str(t)] );
        
        xp.fitter.setPullBack(SMArr2D{tidx, i}, i);
        
    end
    
    % Populate the current SOI
    xp.fitter.repopulateSOI(xp.SOI, xp.currentTime);
    
    % Pullback the stack to the desired charts
    xp.SOI.pullbackStack(xp.stack, [], xp.currentTime, onionOpts);
    
    % Sets the fitter xp.allFit(timepoint) gives that timepoint's fitter
    % xp.allFit(timepoint).fittedParam.submeshes.u{1} will be uv coords
    % xp.allFit(timepoint).fittedParam.submeshes.v will be xyz coords
    xp.setAllFit();
    
end
disp('done with dynamic atlas generation loop')

%% Visualize the Atlas in 2D ==============================================

% The maximum intensity projection of each region at the current time
dataField = xp.SOI.getField('data');
tidx = xp.tIdx(xp.currentTime);
data = dataField(1);

% NOTE: for more information on the organization of the data structures,
% see the supplemental information of the ImSAnE manuscript.
%
% The color version of the pullback, stored to be used as texture for the
% 3D rendering in the next code block

color = cell( numel(data.patches), 1 );

figure,
for i = 1:numel(data.patches)
    
    % the two channels
    R = mat2gray(data.patches{i}.apply{1});
	%G = mat2gray(data.patches{i}.apply{2});
    
    % a little bit of manual adjustment of the lookup table
    R = imadjust(R, [0 0.5]);
    %G = imadjust(G, [0 0.8]);
    
    % concatenate to make color image
    color{i} = cat(3,R,R,R);
    
    % make the background white 
    color{i}(color{i}==0) = 1;
    
    % show the map
    subplot(ceil(numel(data.patches)/2), 2, i)
    imshow(permute(color{i}, [2 1 3]),[],'InitialMagnification',66)
    %imshow(permute(color{i}, [2 1 3]),[])
end

clear R
    
%% Visualize the Atlas in 3D
%
% We take the submeshes and displace them along the average normal by a
% distance set by separation

separation = [140, 0, 600];

figure
hold on;

% The pieces

for i = 1:numel(data.patches)
    
    % the 3D coordinates of a region
    X = xp.SOI.embedding(tidx).patches{i}.apply;
    
    % the mean normal direction
    d = mean(xp.fitter.fittedParam.submeshes{i}.vn);
    d = d./norm(d);
    d = separation(i)*d;
    
    % show each region in 3D displaced along its normal by separation
    surf(X{1} + d(1),X{2} + d(2),X{3} + d(3), 'FaceColor','texturemap'); % color{i},
    %pause
end

axis equal;
axis off;
shading flat
set(gcf,'color','w');
view(60, -60); % set viewing angle

hold off;

%% Save the surface of interest to disc
%
% Here we save the SOI using SOI.save. We set the following options:
%
% * dir:            The directory to save the SOI to.
% * imwriteOptions: Pullbacks are saved to image files using imwrite, we
% can pass options to change file format, compression etc. For example we
% could change this option to
% imwriteOptions = {'jp2', 'Mode', 'lossless'}; 
% * make8bit:       Often absolute intensities don't matter and 8 bit offers
% a large enough dynamic range. This options rescales the lookup table and
% converts to 8 bit before saving.

imwriteOptions = {'tif'};
% soiDir = fullfile(projectDir, 'bad_test');
options = struct('dir', soiDir, 'imwriteOptions', {imwriteOptions},...
                    'make8bit', false);
xp.SOI.save(options)
SOI = xp.SOI ;
disp('done saving SOI')


%% INSPECT SUBMESHES IF NOT ALREADY DONE ==================================
% Plot the collection of submeshes on the surface at each tp.
% 
% for t = xp.fileMeta.timePoints(1:end)
%     disp(['Now reviewing time point ', num2str(t)]);
%     tidx = xp.tIdx(t);
%     
%     % Load data -----------------------------------------------------------
%     xp.loadTime(t);
%     xp.rescaleStackToUnitAspect();
%     mesh = meshStack(tidx);
%     
%     % Save image of the submeshes
%     close all
%     outfn = fullfile(meshcheckDir, ['submeshes' num2str(t, '%04d') '.png']) ;
%     xp.fitter.inspectMesh(1:length(seeds)) ;
%     title(['Submeshes for t=' num2str(t)])
%     saveas(gcf, outfn) ;    
% end

%% Create projection images (vorpatches) ==================================
view_ad = [-100 25] ;  % azimuth and declination of view
em = xp.SOI.embedding ;
soiDir_3dsurf = [soiDir 'apical_surface_projection/'] ;
if ~exist(soiDir_3dsurf, 'dir')
    mkdir(soiDir_3dsurf)
end
scrsz = get(0,'ScreenSize');
for tp = 1:length(fileMeta.timePoints)
    disp(['surface projection: t=' num2str(tp)])
    time = fileMeta.timePoints(tp) ;
    close all;
    figure('Position', [1 scrsz(4)/1 scrsz(3)/1 scrsz(4)/1],...
        'Color',[0 0 0], 'visible', 'off')
    hold off
    for patch = 1:num_vorpatches 
        egrids = em(tp).patches{patch}.apply() ;
        image = xp.SOI.data(tp).patches{patch}.apply() ;
        % cmp = [217,217,215]/255;
        % h = trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),...
        %     'EdgeColor','none','FaceColor',cmp, 'FaceLighting','phong');
        % im = double(xp.SOI.data(time).getPatch(patchName).apply{1});
        % egrids = xp.SOI.embedding(time).patches{1}.apply() ;
        im = double(image{1}) / 65536. ;
        h = surf(-egrids{1}, egrids{2}, -egrids{3}, im);
        hold on
        set(h, 'AmbientStrength', .8)
    end
    axis equal
    shading interp 
    grid off
    colormap gray
    % set(gca,'CLim',[0. 0.5]);
    view(view_ad)
    
    axis off
    % axis([-1000 1000 -400 400 -400 400])
    % g = light('position',[0  1 0]);
    % set(g,'Color',[1 1 1]/2);
    set(gcf,'Color',[0 0 0]);
    
    saveas(gcf, [soiDir_3dsurf sprintf('Time_%06d.png', time)])
end
close all;

%% Create projection images (other patches) ===============================
view_ad = [-100 25] ;  % azimuth and declination of view
em = xp.SOI.embedding ;
soiDir_3dsurf = [soiDir 'apical_surface_projection_disk/'] ;
if ~exist(soiDir_3dsurf, 'dir')
    mkdir(soiDir_3dsurf)
end
scrsz = get(0,'ScreenSize');
for tp = 1:length(fileMeta.timePoints)
    disp(['surface projection: t=' num2str(tp)])
    time = fileMeta.timePoints(tp) ;
    close all;
    figure('Position', [1 scrsz(4)/1 scrsz(3)/1 scrsz(4)/1],...
        'Color',[0 0 0], 'visible', 'off')
    hold off
    for patch = (num_vorpatches + 1):(num_vorpatches + num_diskpatches)
        egrids = em(tp).patches{patch}.apply() ;
        image = xp.SOI.data(tp).patches{patch}.apply() ;
        % cmp = [217,217,215]/255;
        % h = trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),...
        %     'EdgeColor','none','FaceColor',cmp, 'FaceLighting','phong');
        % im = double(xp.SOI.data(time).getPatch(patchName).apply{1});
        % egrids = xp.SOI.embedding(time).patches{1}.apply() ;
        im = double(image{1}) / 65536. ;
        h = surf(-egrids{1}, egrids{2}, -egrids{3}, im);
        hold on
        set(h, 'AmbientStrength', .8)
    end
    axis equal
    shading interp 
    grid off
    colormap gray
    % set(gca,'CLim',[0. 0.5]);
    view(view_ad)
    
    axis off
    % axis([-1000 1000 -400 400 -400 400])
    % g = light('position',[0  1 0]);
    % set(g,'Color',[1 1 1]/2);
    set(gcf,'Color',[0 0 0]);
    
    saveas(gcf, [soiDir_3dsurf sprintf('Time_%06d.png', time)])
end
close all;
disp('done with projection images')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here, we add quasiconformal pipeline to minimize differences.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write IO to put stacks in stacked tifs. Put all stacked tifs in a 
% subdir of projectDir
% Train in ilastik to find the surface to be analyzed.
% Multiply the Probability arrays by the image arrays, sum over z 
% to get new pullbacks at nice surface.
% Binarize: try simple thresholding with Gaussian blur. 
% If not, use Ilastik to segment cells. 
% Plug in the quasiconformal mapping code > in MATLAB. 

%% Load the SOI if already saved here =====================================
SOI = surfaceAnalysis.SurfaceOfInterest(soiDir) ;

%% Collate tifs from each projected layer into a 3D stack =================
addpath('/mnt/data/code/gut_matlab/saveastiff/')
fieldsDir = [soiDir '/fields/'] ;
idx = -floor(nLayers / 2):floor(nLayers / 2) ;
layer_strs = strrep(strrep(compose('data_layer_p%d', idx), ...
    'p-', 'm'), '_layer_p0', '') ;

% Iterate over conformal indices
cmpname = 'cmp_1_1_T'; 
stackname = 'stack_1_1_T';
% Make onion stack for each timepoint
for tt=1:length(xp.fileMeta.timePoints)
    tp = xp.fileMeta.timePoints(tt);
    disp(['considering tp = ' num2str(tp)])
    % Make onion stack for each patch
    for ii = 1:length(SOI.data(tt).patches)
        disp(['path ' num2str(ii)])
        % Build a cell array of image paths
        tiffn = [sprintf('/conformal_%d_index/conformal_%d/', ii, ii), ...
                 cmpname, sprintf('%04d.tif', tp)] ;
        impaths = repmat({''}, length(layer_strs), 1) ;
        for jj = 1:length(layer_strs)
            impaths{jj} = [fieldsDir layer_strs{jj} tiffn] ;
            % impaths{end + 1} = [fieldsDir layer_strs{jj} tiffn] ;
        end
        % read images into cell array
        Im = cellfun(@imread, impaths,'uni',false);
        % Concatenate images into single 3D matrix along third dimension
        myImage = cat(3, Im{:});
        % Create output directory
        odir = sprintf('%s/stacks_conformal_%d_index_%d/', ...
            fieldsDir, ii, ii);
        if ~exist(odir, 'dir')
            mkdir(odir)
        end
        outfn = [odir sprintf('%s%04d.tif',stackname,tp)];
        options.overwrite = true;
        saveastiff(myImage, outfn, options) ;
    end
end

%% Train onion stacks in ilastik to find the surface to be analyzed =======

%% Multiply the output of ilastik (Probability arrays for onions) by the 
% image arrays, sum over z to get new pullbacks at nice surface.

%% Prepare for quasiconformal timepoint matching ==========================
% First run ilastik training on the pullbacks using three colors: 
% membrane, bg, and fuzzy data
addpath_recurse('/mnt/data/code/gut_matlab/tissueAnalysisSuite/')
ilpmaster = [projectDir 'gut_apical_cylinder_msls_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1_3layer/' ] ;
ilastiksegFolder1 = [soiDir 'fields/data_ilastikseg/conformal1/'] ;
ilastiksegFolder2 = [soiDir 'fields/data_ilastikseg/conformal2/'] ;

%% First compute order parameter cos(2*theta)
% add paths for segmentation
addpath('/mnt/data/code/gut_matlab/tissueAnalysisSuite/')

% Create a label matrix
very_far = 150 ;  % distance between points that can't be a cell, in pixels
mode = 0;  % Toggle for ilastik version control
disp('loading h5 data in ilastik segmentation folder...')
mem1 = load.ilastikh5( ilastiksegFolder1, mode ) ; 
mem2 = load.ilastikh5( ilastiksegFolder2, mode ) ; 

% Could get bonds via the following
% disp('segmenting the data...')
% L = seg.memWS(mem, 50, 0, 1, 3.5) ;
% % Set bond=0 and clear_border = 1
% [L, Struct] = seg.generate_structs(L, 0, 1, 0, very_far);
% disp('done with initial segmentation')
% vdat = Struct.Vdat ;
% % Collate the vertices into array
% xv = zeros(length(vdat), 1) ;
% yv = zeros(length(vdat), 1) ;
% for i=1:length(vdat)
%     xv(i) = vdat(i).vertxcoord ;
%     yv(i) = vdat(i).vertycoord ;
% end
% xy = [xv, yv] ;
% % Collate the lines into vectors
% xv = zeros(length(vdat), 1) ;
% yv = zeros(length(vdat), 1) ;
% % Use TRIsm2NL() from npm

