%% GENERATE_AXISYMMETRIC_PULLBACK ==========================================
% Pipeline for creating 'axisymmetric' pullbacks of the growing Drosophila 
% midgut about a centerline (previously computed). Uses orbifold method to
% guarantee conformality except at x=0 and x=L
% Run from anywhere, just specify dataDir below.
%
% Inputs
% ------
% meshes as topological cylinders, in meshDir/cylindercut/
% meshDir/translation_APDV.txt
% meshDir/rotation_APDV.txt
% meshDir/xyzlim_APDV_um.txt
% 
% Returns
% -------
% meshDir/meshStack.mat
% PullbackImages_###step/ (images)
% 
% By Dillon Cislo and Noah Mitchell
%==========================================================================

clear; close all; clc;

%% Parameters
overwrite = false ;
resave_ims = false ;
save_ims = true ;
normal_shift = 10 ;
a_fixed = 2 ;
preview = false ;
washout2d = 0.5 ;
washout3d = 0.5 ;

%% Add paths
% Add some necessary code to the path (ImSAnE should also be setup!) ------
addpath(genpath('/mnt/crunch/djcislo/MATLAB/euclidean_orbifolds'));
addpath(genpath('/mnt/data/code/gptoolbox'));
addpath(genpath('/mnt/data/code/gut_matlab/TexturePatch'));
addpath('/mnt/data/code/gut_matlab/') ;
addpath_recurse('/mnt/data/code/gut_matlab/axisymmetric_pullbacks/') ;
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
dataDir = [ '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/', ...
    'Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/' ];

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
meshDir = [ dataDir ...
    'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1' ];

% The dataset name base for the AD points
ADBase = '/mesh_apical_stab_%06d/adorsal';

% The dataset name base for the PD points
PDBase = '/mesh_apical_stab_%06d/pdorsal';

% The folder where the pullback images will be saved
nshift = strrep(sprintf('%03d', normal_shift), '-', 'n') ;
imFolder = fullfile(meshDir, ['PullbackImages_' nshift 'step'] ) ;
imFolder_e = [imFolder '_extended'] ;
imFolder_r = [imFolder '_relaxed'] ;
imFolder_re = [imFolder '_relaxed_extended'] ;
pivDir = fullfile(meshDir, 'piv') ;
cutFolder = fullfile(meshDir, 'cutMesh') ;
cutMeshImagesDir = fullfile(cutFolder, 'images') ;
cylCutDir = fullfile(meshDir, 'cylindercut') ;
cylCutMeshOutDir = fullfile(cylCutDir, 'cleaned') ;

% The file name base for the cylinder meshes
cylinderMeshBase = fullfile( cylCutDir, ...
    'mesh_apical_stab_%06d_cylindercut.ply' );
cylinderMeshCleanBase = fullfile( cylCutMeshOutDir, ...
    'mesh_apical_stab_%06d_cylindercut_clean.ply' );

% The file constaing the AD/PD points
dpFile = fullfile( cylCutDir, 'ap_boundary_dorsalpts.h5' );

tomake = {imFolder, imFolder_e, imFolder_r, imFolder_re,...
    pivDir, cutFolder, cutMeshImagesDir, cylCutMeshOutDir} ;
for i = 1:length(tomake)
    dir2make = tomake{i} ;
    if ~exist( dir2make, 'dir' )
        mkdir(dir2make);
    end
end

%% Load rotation, translation, resolution
rot = dlmread(fullfile(meshDir, 'rotation_APDV.txt')) ;
xyzlim = dlmread(fullfile(meshDir, 'xyzlim_APDV_um.txt'), ',', 1, 0) ;
trans = dlmread(fullfile(meshDir, 'translation_APDV.txt'));
resolution = dlmread(fullfile(meshDir, 'resolution.txt'), ',', 1, 0) ;
if resolution(:) == resolution(1)
    resolution = resolution(1) ;
else
    error('Have not handled case for anisotropic resolution')
end

%% Iterate Through Time Points to Create Pullbacks ========================
% Check if cutmeshes already saved
mstckfn = fullfile(meshDir, 'meshStack_orbifold.mat') ;
outcutfn = fullfile(cutFolder, 'cutPaths.h5') ;
outadIDxfn = fullfile(cylCutMeshOutDir, 'adIDx.h5') ;
outpdIDxfn = fullfile(cylCutMeshOutDir, 'pdIDx.h5') ;
if exist(mstckfn, 'file') && ~overwrite
    % The cutPaths.h5 is loading here
    disp(['Loading meshStack: ' mstckfn])
    load(mstckfn)
    
    if resave_ims
        for t = xp.fileMeta.timePoints %xp.fileMeta.timePoints(1:50)
            % Plot the cutMesh3D
            % Load the mesh
            disp(['NOW PROCESSING TIME POINT ', num2str(t)]);
            tidx = xp.tIdx(t);

            % Load the cylinder mesh
            mesh = read_ply_mod( sprintf( cylinderMeshCleanBase, t ) );
            % Consistently orient mesh faces
            mesh.f = bfs_orient( mesh.f );

            % Load the cut from h5 file
            cutP = h5read(outcutfn, ['/' sprintf('%06d', t)]) ;

            %% Plot the cutPath (cutP) in 3D
            disp('Plotting cut...')
            xyzrs = ((rot * mesh.v')' + trans) * resolution ;
            fig = figure('Visible', 'Off')  ;
            fig.PaperUnits = 'centimeters';
            fig.PaperPosition = [0 0 12 12];
            sh = trimesh(mesh.f, ...
                xyzrs(:, 1), xyzrs(:,2), xyzrs(:, 3), xyzrs(:, 2), ...
                'FaceColor', 'interp', 'edgecolor', 'none', 'FaceAlpha', 0.3) ;
            hold on;
            ph = plot3(xyzrs(cutP, 1), xyzrs(cutP, 2), xyzrs(cutP, 3), 'k-', 'LineWidth', 3) ;
            axis equal
            xlim(xyzlim(1, :)); 
            ylim(xyzlim(2, :)); 
            zlim(xyzlim(3, :));
            title(['t=' sprintf('%04d', t)]) ;
            xlabel('x [$\mu$m]', 'Interpreter', 'Latex') ;
            ylabel('y [$\mu$m]', 'Interpreter', 'Latex') ;
            zlabel('z [$\mu$m]', 'Interpreter', 'Latex') ;
            cutfn = sprintf( fullfile(cutMeshImagesDir, [fileNameBase, '_cut.png']), t ) ;
            saveas(fig, cutfn)
            close all 
        end
    end
else
    meshStack = cell( length(xp.fileMeta.timePoints), 1 );

    for t = xp.fileMeta.timePoints %xp.fileMeta.timePoints(1:50)

        disp(['NOW PROCESSING TIME POINT ', num2str(t)]);

        tidx = xp.tIdx(t);

        % Load the data for the current time point ------------------------
        xp.setTime(t) ;
        
        % Load or compute clean cylindrical mesh
        mesh3dfn =  sprintf( cylinderMeshCleanBase, t ) ;
        outadIDxfn = fullfile(cylCutMeshOutDir, 'apIDx.h5') ;
        outpdIDxfn = fullfile(cylCutMeshOutDir, 'pdIDx.h5') ;
        if ~exist(mesh3dfn, 'file') || overwrite
            % Load the cylinder mesh
            mesh = read_ply_mod( sprintf( cylinderMeshBase, t ) );

            % Consistently orient mesh faces
            disp('orienting faces')
            mesh.f = bfs_orient( mesh.f );

            % Load the AD/PD vertex IDs
            disp('Loading ADPD vertex IDs...')
            adIDx = h5read( dpFile, sprintf( ADBase, t ) );
            pdIDx = h5read( dpFile, sprintf( PDBase, t ) );

            % Clip the ears of the triangulation and update the AD/PD points if
            % necessary -----------------------------------------------------------
            ad3D = mesh.v( adIDx, : );
            pd3D = mesh.v( pdIDx, : );
            disp('Clipping ears...') ;
            [ ff, ~, vv, newind] = clip_mesh_mod( mesh.f, mesh.v );

            mesh.f = ff; mesh.v = vv;
            % Remove zeros from newind in order to assign the new vtx normals
            % 
            newind = newind(newind > 0) ;
            mesh.vn = mesh.vn(newind, :) ;

            adIDx = pointMatch( ad3D, mesh.v );
            pdIDx = pointMatch( pd3D, mesh.v );

            clear ff vv ad3D pd3D
            
            %% Save the 3d cut mesh with new indices
            plywrite_with_normals(mesh3dfn, mesh.f, mesh.v, mesh.vn)
            % adIDx
            try
                h5create(outadIDxfn, ['/' sprintf('%06d', t) ], size(adIDx)) ;
            catch
                disp(['adIDx for t=' num2str(t) ' already exists'])
            end
            h5write(outadIDxfn, ['/' sprintf('%06d', t)], adIDx) ;
            % pdIDx
            try
                h5create(outpdIDxfn, ['/' sprintf('%06d', t) ], size(pdIDx)) ;
            catch
                disp(['adIDx for t=' num2str(t) ' already exists'])
            end
            h5write(outpdIDxfn, ['/' sprintf('%06d', t)], pdIDx) ;
            
            disp('done with cylindermesh cleaning')
        else
            mesh = read_ply_mod(mesh3dfn) ;
            adIDx = h5read(outadIDxfn, ['/' sprintf('%06d', t)]) ;
            pdIDx = h5read(outpdIDxfn, ['/' sprintf('%06d', t)]) ;
        end
        
        % View results --------------------------------------------------------
        % trisurf( triangulation( mesh.f, mesh.v ) );
        % hold on
        % scatter3( mesh.v(adIDx,1), mesh.v(adIDx,2), mesh.v(adIDx,3), ...
        %     'filled', 'r' );
        % scatter3( mesh.v(pdIDx,1), mesh.v(pdIDx,2), mesh.v(pdIDx,3), ...
        %     'filled', 'c' );
        % hold off
        % axis equal

        %----------------------------------------------------------------------
        % Create the Cut Mesh
        %----------------------------------------------------------------------
        fprintf('Generating Cut Mesh... ');
        cutMeshfn = fullfile(cutFolder, [fileNameBase, '_cutMesh.mat']) ;
        cutMeshfn = sprintf(cutMeshfn, t) ;
        if ~exist(cutMeshfn, 'file') || overwrite
            try
                [ cutMesh, adIDx, pdIDx, cutP ] = ...
                    cylinderCutMesh( mesh.f, mesh.v, mesh.vn, adIDx, pdIDx );
                compute_pullback = true ;
            catch
                disp('Could not cut this timepoint: Input mesh probably NOT a cylinder')
                cutP = [] ;
                compute_pullback = false ;
            end
            fprintf('Done with generating CutMesh\n');

            % Save the cut in h5 file
            try
                h5create(outcutfn, ['/' sprintf('%06d', t) ], size(cutP)) ;
            catch
                disp(['cut for t=' num2str(t) ' already exists'])
            end
            h5write(outcutfn, ['/' sprintf('%06d', t)], cutP) ;
            
            % Save cutMesh
            save(cutMeshfn, 'cutMesh', 'adIDx', 'pdIDx')

            disp('done with plotting & saving cut')
        else
            load(cutMeshfn) 
            cutP = h5read(outcutfn, ['/' sprintf('%06d', t)]) ;
            compute_pullback = ~isempty(cutP) ;
        end
        
        %% Plot the cutPath (cutP) in 3D
        if compute_pullback && save_ims
            disp('Plotting cut...')
            xyzrs = ((rot * mesh.v')' + trans) * resolution ;
            fig = figure('Visible', 'Off')  ;
            fig.PaperUnits = 'centimeters';
            fig.PaperPosition = [0 0 12 12];
            sh = trimesh(mesh.f, ...
                xyzrs(:, 1), xyzrs(:,2), xyzrs(:, 3), xyzrs(:, 2), ...
                'FaceColor', 'interp', 'edgecolor', 'none', 'FaceAlpha', 0.3) ;
            hold on;
            ph = plot3(xyzrs(cutP, 1), xyzrs(cutP, 2), xyzrs(cutP, 3), 'k-', 'LineWidth', 3) ;
            axis equal
            xlim(xyzlim(1, :)); 
            ylim(xyzlim(2, :)); 
            zlim(xyzlim(3, :));
            title(['t=' sprintf('%04d', t)]) ;
            xlabel('x [$\mu$m]', 'Interpreter', 'Latex') ;
            ylabel('y [$\mu$m]', 'Interpreter', 'Latex') ;
            zlabel('z [$\mu$m]', 'Interpreter', 'Latex') ;
            cutfn = sprintf( fullfile(cutMeshImagesDir, [fileNameBase, '_cut.png']), t ) ;
            saveas(fig, cutfn)
            close all 
        end
        
        %% Compute the pullback if the cutMesh is ok
        if compute_pullback
            % Displace normally ---------------------------------------------------
            cutMesh.v = cutMesh.v + cutMesh.vn * normal_shift ;

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
            %     'filled', 'r' );
            % scatter3( mesh.v(pdIDx,1), mesh.v(pdIDx,2), mesh.v(pdIDx,3), ...
            %     'filled', 'm' );
            % 
            % hold off
            % 
            % axis equal
            % 
            % clear P

            %----------------------------------------------------------------------
            % Generate Pullback to Annular Orbifold Domain
            %----------------------------------------------------------------------
            fprintf('Generating Pullback... ');
            cutMesh = flattenAnnulus( cutMesh );

            % Find lateral scaling that minimizes spring network energy
            ar = minimizeIsoarealAffineEnergy( cutMesh.f, cutMesh.v, cutMesh.u );
            % Assign scaling based on options: either a0 or a_fixed
            if tidx == 1 && ~a_fixed
                a_fixed = ar ;
            end      
            a = a_fixed ;
      
            % Scale the x axis by a or ar
            uvtx = cutMesh.u ;
            cutMesh.u = [ a .* uvtx(:,1), uvtx(:,2) ];
            cutMesh.urelax = [ ar .* uvtx(:,1), uvtx(:,2) ];
            fprintf('Done\n');

            % View results --------------------------------------------------------
            % 
            % patch( 'Faces', cutMesh.f, 'Vertices', cutMesh.u, ...
            %     'FaceVertexCData', cutMesh.v(:,3), 'FaceColor', 'interp', ...
            %     'EdgeColor', 'k' );
            % 
            % hold on
            % 
            % cornerColors = [ 1 0 0; 1 0 1; 0 1 1; 0 1 0 ];
            % corners = [ adIDx cutMesh.pathPairs(1,2), ...
            %     cutMesh.pathPairs(end,1) pdIDx ];
            % 
            % scatter( cutMesh.u( corners, 1 ), cutMesh.u( corners, 2 ), [], ...
            %     cornerColors, 'filled' );
            % 
            % hold off
            % 
            % axis equal
            % 
            % clear cornerColors corners

            %----------------------------------------------------------------------
            % Generate Output Image File
            %----------------------------------------------------------------------
            imfn = sprintf( fullfile([imFolder, '/', fileNameBase, '.tif']), t ); 
            if ~exist(imfn, 'file')
                fprintf(['Generating output image: ' imfn]);

                % Texture patch options
                Options.PSize = 5;
                Options.EdgeColor = 'none';

                % Load 3D data for coloring mesh pullback
                xp.loadTime(t);
                xp.rescaleStackToUnitAspect();

                % Raw stack data
                IV = xp.stack.image.apply();
                IV = imadjustn(IV{1});

                
                %% OLD VERSION
                % % Make a full screen image
                % figure('units', 'normalized', ...
                %     'outerposition', [0 0 1 1], 'visible', 'off')
                % 
                % % Plot texture patch cut mesh in 2D
                % texture_patch_3d( cutMesh.f, cutMesh.u, ...
                %     cutMesh.f, cutMesh.v(:, [2 1 3]), IV, Options );
                % 
                % hold on
                % 
                % % Plot the upper tile
                % texture_patch_3d( cutMesh.f, [ cutMesh.u(:,1), cutMesh.u(:,2)+1 ], ...
                %     cutMesh.f, cutMesh.v(:, [2 1 3]), IV, Options );
                % 
                % % Plot the lower tile
                % texture_patch_3d( cutMesh.f, [ cutMesh.u(:,1), cutMesh.u(:,2)-1 ], ...
                %     cutMesh.f, cutMesh.v(:, [2 1 3]), IV, Options );
                % 
                % hold off
                % 
                % axis equal
                % 
                % % Format axes
                % xlim([0 a]); ylim([0 1]);
                % set(gca, 'xtick', []);
                % set(gca, 'ytick', []);
                % 
                % colormap gray
                % 
                % % Extract image from figure axes
                % patchIm = getframe(gca);
                % patchIm = rgb2gray(patchIm.cdata);

                
                %% NEW Version

                fprintf('Generating output image... ');

                % Generate Tiled Orbifold Triangulation ------------------------------
                tileCount = [1 1];  % how many above, how many below
                [ TF, TV2D, TV3D ] = tileAnnularCutMesh( cutMesh, tileCount );

                % View Results -------------------------------------------------------
                % patch( 'Faces', TF, 'Vertices', TV2D, 'FaceVertexCData', ...
                %     TV3D(:,3), 'FaceColor', 'interp', 'EdgeColor', 'k' );
                % axis equal

                % Texture image options
                Options.imSize = ceil( 1000 .* [ 1 1 ] );
                Options.yLim = [0 1];

                % Raw stack data
                IV = xp.stack.image.apply();
                IV = imadjustn(IV{1});

                % profile on
                % Create texture image
                patchIm = texture_patch_to_image( TF, TV2D, TF, TV3D(:, [2 1 3]), ...
                    IV, Options );
                % profile viewer

                fprintf('Done\n');

                % View results ------------------------------------------------------------
                % imshow( patchIm );
                % set( gca, 'YDir', 'Normal' );
                
                % Write figure to file
                imwrite( patchIm, imfn, 'TIFF' );

                % Close open figures
                close all
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Save relaxed image
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            imfn_r = sprintf( fullfile([imFolder_r, '/', fileNameBase, '.tif']), t ) ;
            if ~exist(imfn_r, 'file')
                disp('Generating relaxed image...')

                % Generate Tiled Orbifold Triangulation ------------------------------
                tileCount = [2 2];  % how many above, how many below
                [ TF, TV2D, TV3D ] = tileAnnularCutMesh( cutMesh, tileCount );

                % View Results -------------------------------------------------------
                % patch( 'Faces', TF, 'Vertices', TV2D, 'FaceVertexCData', ...
                %     TV3D(:,3), 'FaceColor', 'interp', 'EdgeColor', 'k' );
                % axis equal

                % Texture image options
                Options.imSize = ceil( 1000 .* [ 1 ar ] );
                Options.yLim = [0 1];

                % Raw stack data
                IV = xp.stack.image.apply();
                IV = imadjustn(IV{1});

                % profile on
                % Create texture image
                patchIm = texture_patch_to_image( TF, TV2D, TF, TV3D(:, [2 1 3]), ...
                    IV, Options );
                % profile viewer

                fprintf('Done\n');
                
                
                % Write figure to file
                imwrite( patchIm, ...
                    sprintf( fullfile([imFolder_r, '/', fileNameBase, '.tif']), t ), ...
                    'TIFF' );            

                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % % Save extended relaxed image
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                % disp('Generating relaxed, extended image...')          
                % % Format axes
                % xlim([0 ar]); ylim([-0.5 1.5]);
                % 
                % % Extract image from figure axes
                % patchIm_e = getframe(gca);
                % patchIm_e = rgb2gray(patchIm_e.cdata);
                % 
                % % Write figure to file
                % imwrite( patchIm_e, ...
                %     sprintf( fullfile([imFolder_re, '/', fileNameBase, '.tif']), t ), ...
                %     'TIFF' );
            
                % Close open figures
                close all
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Save submesh array. Each cell element contains all the 
            % submeshes for that TP, which in this case is just one.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            meshStack{tidx} = cutMesh ;
                        
            clear Options IV

            fprintf('Done\n');
        end
    end

    %% Save SMArr2D (vertex positions in the 2D pullback) -----------------
    disp(['Saving meshStack to disk: ' mstckfn])
    save(mstckfn, 'meshStack') ;
end


%% Preview results
check = false ;
if check
    %% Show 3D Texture Patch (Optional) ---------------------------------------

    % Texture patch options
    Options.PSize = 5 ;
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

    %% Show 2D Texture Patch --------------------------------------------------
    % Texture patch options
    Options.PSize = 5;
    Options.EdgeColor = 'none';

    % Raw stack data
    IV = xp.stack.image.apply();
    IV = imadjustn(IV{1});

    % Make a full screen image
    figure('units', 'normalized', 'outerposition', [0 0 1 1])

    % Plot texture patch cut mesh in 2D
    texture_patch_3d( cutMesh.f, cutMesh.u, ...
        cutMesh.f, cutMesh.v(:, [2 1 3]), IV, Options );

    axis equal

    % Format axes
    % xlim([0 1]); ylim([0 1]);
    % set(gca, 'xtick', []);
    % set(gca, 'ytick', []);

    colormap bone

    clear Options IV

    %% Show Tiled 2D Texture Patch --------------------------------------------

    % Texture patch options
    Options.PSize = 5;
    Options.EdgeColor = 'none';

    % Raw stack data
    IV = xp.stack.image.apply();
    IV = imadjustn(IV{1});

    % Make a full screen image
    figure('units', 'normalized', 'outerposition', [0 0 1 1])

    % Plot texture patch cut mesh in 2D
    texture_patch_3d( cutMesh.f, cutMesh.u, ...
        cutMesh.f, cutMesh.v(:, [2 1 3]), IV, Options );

    hold on

    texture_patch_3d( cutMesh.f, [ cutMesh.u(:,1), cutMesh.u(:,2)+1 ], ...
        cutMesh.f, cutMesh.v(:, [2 1 3]), IV, Options );

    texture_patch_3d( cutMesh.f, [ cutMesh.u(:,1), cutMesh.u(:,2)-1 ], ...
        cutMesh.f, cutMesh.v(:, [2 1 3]), IV, Options );

    axis equal

    % Format axes
    xlim([0 a]); ylim([0 1]);
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);

    colormap bone
    clear Options IV
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TILE IMAGES IN Y AND RESAVE ============================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fns = dir(strrep(fullfile([imFolder, '/', fileNameBase, '.tif']), '%06d', '*')) ;

% Get original image size
im = imread(fullfile(fns(1).folder, fns(1).name)) ;
halfsize = round(0.5 * size(im, 1)) ;
osize = size(im) ;

for i=1:length(fns)
    if ~exist(fullfile(imFolder_e, fns(i).name), 'file')
        disp(['Reading ' fns(i).name])
        fileName = split(fns(i).name, '.tif') ;
        fileName = fileName{1} ;
        im = imread(fullfile(fns(i).folder, fns(i).name)) ;

        % im2 is as follows:
        % [ im(end-halfsize) ]
        % [     ...          ]
        % [    im(end)       ]
        % [     im(1)        ]
        % [     ...          ]
        % [    im(end)       ]
        % [     im(1)        ]
        % [     ...          ]
        % [  im(halfsize)    ]
        im2 = uint8(zeros(size(im, 1) + 2 * halfsize, size(im, 2))) ;
        im2(1:halfsize, :) = im(end-halfsize + 1:end, :);
        im2(halfsize + 1:halfsize + size(im, 1), :) = im ;
        im2(halfsize + size(im, 1) + 1:end, :) = im(1:halfsize, :);
        imwrite( im2, fullfile(imFolder_e, fns(i).name), 'TIFF' );
    else
        disp('already exists')
    end
    
end
disp('done writing extended tiffs')

%% PLOT RADIUS ============================================================


error('breaking here')
%% PERFORM PIV ============================================================
% Select all frames in PullbackImages_extended_shifted/ 
% Select Sequencing style 1-2, 2-3, ... 
% Image Preprocessing > Select All
% PIV settings: 128 (32 step), 64 (16 step) for two passes
% Post-processing: vector validation: Select velocity limits, 5 stdev
% filter, local median filter: thres=2, eps=0.1
% with default values.
% File > Save > Export all frames as pivresults.mat
disp('Loading pivresults_orbifold_pass0.mat ...')
load(fullfile(pivDir, 'pivresults_orbifold_pass0.mat'))

%% Get timestamps for the images with pullbacks
fns = dir(fullfile(imFolder_e, '*.tif')) ;
npiv = length(fns) ;
time = zeros(npiv, 1);
meshidx = zeros(npiv, 1) ;
for i=1:npiv
    tmp = split(fns(i).name, '.tif') ;
    tmp = split(tmp{1}, '_') ;
    time(i) = str2double(tmp{2}) ;
    meshidx(i) = time(i) - time(1) + 1 ;
end
dt = diff(time) ;
disp('done building dt') 

%% Subtract off the mean flow in y for each frame ========================= 
meanv = zeros(length(v_filtered), 1) ;
for i=1:length(v_filtered)
    tmp = v_filtered{i} ;
    meanv(i) = mean(tmp(:)) ;
end
dy = round(meanv) ;
shifty = [0; -cumsum(dy)] ;
disp('computed shifty')
save(fullfile(pivDir, 'shifty.mat'), 'shifty')

%% Shift each frame by shifty and resave ==================================
fns = dir(strrep(fullfile([imFolder_e, '/', fileNameBase, '.tif']), '%06d', '*')) ;
imFolder_es = [imFolder '_extended_shifted' filesep] ;
if ~exist(imFolder_es, 'dir')
    mkdir(imFolder_es) ;
end

% Write the images to disk and get their sizes while being written 
disp('Resaving shifted images if necessary...')
imsizes = zeros(length(fns), 2) ;
for i=1:length(fns)
    outfn = fullfile(imFolder_es, fns(i).name) ;
    if ~exist(outfn, 'file')
        disp(['Reading ' fns(i).name])
        fileName = split(fns(i).name, '.tif') ;
        fileName = fileName{1} ;
        im = imread(fullfile(fns(i).folder, fns(i).name)) ;
        im = circshift(im, shifty(i), 1) ;
        imsizes(i, :) = size(im) ;
        imwrite( im, outfn, 'TIFF' );
    else
        im = imread(outfn) ;
        imsizes(i, :) = size(im) ;
    end
end
disp('done with shifted images')

%% PERFORM PIV ON SHIFTED FRAMES ==========================================
% Select all frames in PullbackImages_extended_shifted/ 
% Select Sequencing style 1-2, 2-3, ... 
% Load settings: piv_set_pass1.mat
% Image Preprocessing > Select All
% PIV settings: 128 (32 step), 64 (16 step), 32 (16 step) for three passes
disp('Loading pivresults_orbifold_pass1.mat ...')
load(fullfile(pivDir, 'pivresults_orbifold_pass1.mat'))


%% MAKE MAP FROM PIXEL TO XYZ =============================================
disp('Making map from pixel to xyz to compute velocities in 3d...')
% Get position of embedding points associated with velocities 
pivOutDir = fullfile(pivDir, 'images') ;
if ~exist(pivOutDir, 'dir')
    mkdir(pivOutDir)
end

% Compute size of domain (0,1), (0, 1)
fns = dir(strrep(fullfile([imFolder_es, '/', fileNameBase, '.tif']), '%06d', '*')) ;
im = imread(fullfile(fns(i).folder, fns(i).name)) ;
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

piv3dfn = fullfile(pivDir, 'piv3d.mat') ;

if exist(piv3dfn, 'file')
    load(piv3dfn)
else
    piv3d = cell(length(fns), 1) ;
    
    %% Iterate over all images with flow fields
    for i=1:length(fns) - 1
        timestr = sprintf('%03d', time(i)) ;
        disp(['t = ' timestr])

        % Get scale of image
        ysc0 = imsizes(i, 2) ;
        ysc1 = imsizes(i+1, 2) ;

        mesh0 = meshStack{meshidx(i)} ;    
        mesh1 = meshStack{meshidx(i + 1)} ;
        assert(time(i) + dt(i) == time(i+ 1))

        % Load the positions of the velocity vectors in pixels
        x0 = x{i} ;
        y0 = y{i} ;
        uu = u_filtered{i} ;
        vv = v_filtered{i} ; 
        % Get position in next timepoint in pixels (in xy plane)
        % Clip the x position to the size of the image, and wrap the y position
        x1 = x0 + uu ;
        x1 = max(x1, 0) ;
        x1 = min(x1, imsizes(i+1, 1)) ;
        y1 = mod(y0 + vv, ysc1) ;

        % get embedded vector in R^3 for t0
        % Obtain the equivalent of v2D: ie the 2D vertices: mesh0.u
        % The shift in pixels of the current frame = shifty(i)
        mesh0x = u2pix(mesh0.u(:, 1), ysc0) ;
        mesh0y = y2pix(mesh0.u(:, 2), shifty(i), ysc0) ;

        % Create extended mesh (copied above and below), also shifted by shifty
        meshxy = [mesh0x, mesh0y ] ;
        mabove = [mesh0x, mesh0y + u2pix(1., ysc0)] ;
        mbelow = [mesh0x, mesh0y - u2pix(1., ysc0)] ;
        mabove2 = [mesh0x, mesh0y + 2 * u2pix(1., ysc0)] ;
        mbelow2 = [mesh0x, mesh0y - 2 * u2pix(1., ysc0)] ;
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

        if save_ims
            % Check out the triangulation    
            close all
            fig = figure('Visible', 'Off');
            tr0_orig = triangulation(mf0, meshxy) ;
            triplot(tr0, 'Color', blue)
            hold on
            triplot(tr0_orig, 'Color', orange)
            axis equal
            title('Triangulation in pullback space')
            ylim([0 yesz])
            xlabel('x [pix]')
            ylabel('y [pix]')
            saveas(fig, fullfile(pivOutDir, ['tri_' timestr '.png']))
            if preview
                set(fig, 'Visible', 'On')
                waitfor(fig)
            end
            close all

            triplot(tr0_orig, 'Color', orange)
        end

        % Find xyz for matching position in t1 xy plane
        % get embedded vector in R^3 for t0
        % The shift in pixels of the current frame = shifty(i)
        mesh1x = u2pix(mesh1.u(:, 1), ysc1) ;
        mesh1y = y2pix(mesh1.u(:, 2), shifty(i + 1), ysc1) ;

        % Create extended mesh (copied above and below), also shifted by shifty
        meshxy = [mesh1x, mesh1y ] ;
        mabove = [mesh1x, mesh1y + u2pix(1., ysc1)] ;
        mbelow = [mesh1x, mesh1y - u2pix(1., ysc1)] ;
        mabove2 = [mesh1x, mesh1y + u2pix(2., ysc1)] ;
        mbelow2 = [mesh1x, mesh1y - u2pix(2., ysc1)] ;
        m1xy = [meshxy; mabove; mbelow; mabove2; mbelow2] ;
        mf1 = mesh1.f ;
        mf1c = [mf1; mf1 + length(mesh1x); mf1 + 2 * length(mesh1x); ...
             mf1 + 3 * length(mesh1x); mf1 + 4 * length(mesh1x)] ;
        tr1 = triangulation(mf1c, m1xy) ;
        % x123(i) is the x coords of the elements of triangle t_contain(i)
        [t1_contain, baryc1] = pointLocation(tr1, [x1(:), y1(:)]) ;
        vxb = mesh1.v(:, 1) ;
        vyb = mesh1.v(:, 2) ;
        vzb = mesh1.v(:, 3) ;
        trisb = mod(tr1.ConnectivityList(t1_contain, :), size(vxb, 1)) ;
        trisb(trisb == 0) = size(vxb, 1) ;
        x123b = vxb(trisb) ;
        y123b = vyb(trisb) ;
        z123b = vzb(trisb) ;
        pt1 = [sum(baryc1 .* x123b, 2), sum(baryc1 .* y123b, 2), sum(baryc1 .* z123b, 2) ] ;
        
        %% Take the difference to get the velocity field ------------------
        v0 = (pt1 - pt0) / dt(i) ;
        
        %% Visualize the flow in 3d ---------------------------------------
        if save_ims
            close all
            fig = figure('Visible', 'Off') ;
            quiver3(pt0(:, 1), pt0(:, 2), pt0(:, 3), ...
                v0(:, 1), v0(:, 2), v0(:, 3), 0)
            axis equal
            title(['t=' timestr]) 
            saveas(gcf, fullfile(pivOutDir, ['piv3d_' timestr '.png']))
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
        
        % Compare with finding normals from vertex positions cross product
        v21 = mesh0.v(trisa(:, 2), :) - mesh0.v(trisa(:, 1), :) ;
        v31 = mesh0.v(trisa(:, 3), :) - mesh0.v(trisa(:, 1), :) ;
        nx = v21(:, 2) .* v31(:, 3) - v21(:, 3) .* v31(:, 2) ;
        ny = v21(:, 3) .* v31(:, 1) - v21(:, 1) .* v31(:, 3) ;
        nz = v21(:, 1) .* v31(:, 2) - v21(:, 2) .* v31(:, 1) ;
        norms = [nx, ny, nz] ;
        norms = norms ./ sqrt(sum(norms.^2, 2)) ;
        
        % check normals
        if preview
            fig = figure ;
            h1 = plot(norms(:, 1), normals(:, 1), '.') ;
            hold on
            h2 = plot(norms(:, 2), normals(:, 2), '.') ;
            h3 = plot(norms(:, 3), normals(:, 3), '.') ;
            legend('nx', 'ny', 'nz')
            xlabel('triangle normals')
            ylabel('averaged vertex normals')
            waitfor(fig)
        end
        
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
        
        % Check normals
        if preview
            % Look at fractional change
            fig = figure ;
            plot(sqrt(sum((tv3proj - tv1).^2, 2)) ./ sqrt(sum(tv1.^2, 2)))
            xlabel('vertex index')
            ylabel('$|v_p - v| / |v|$', 'Interpreter', 'Latex')
            title('Change when projecting triangles to smoothed tangent planes')
            waitfor(fig)
            
            fig = figure('Visible', 'On')
            quiver3(pt0(:, 1), pt0(:, 2), pt0(:, 3), ...
                normals(:, 1), normals(:, 2), normals(:, 3), 'Color', blue)
            hold on
            quiver3(pt0(:, 1), pt0(:, 2), pt0(:, 3), ...
                norms(:, 1), norms(:, 2), norms(:, 3), 'Color', orange)
            waitfor(fig)
            fig = figure;
            hist(sum(normals.^2, 2))
            waitfor(fig)
        end
        
        %% Take dot product of flow fields with normals
        v0n = dot(normals, v0, 2) ;
        % Subtract the normal velocity to obtain tangential velocity
        v0t = v0 - v0n .* normals ;
        
        
        %% Compute the tangential velocities in plane
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
        piv3d{i} = datstruct ;

        %% Draw 2D flows
        if save_ims        
            vtdir2d = fullfile(pivDir, 'vt2d') ;
            vndir2d = fullfile(pivDir, 'vn2d') ;
            if ~exist(vtdir2d, 'dir')
                mkdir(vtdir2d)
            end
            if ~exist(vndir2d, 'dir')
                mkdir(vndir2d)
            end
            
            % Load the image to put flow on top
            fileName = split(fns(i).name, '.tif') ;
            fileName = fileName{1} ;
            im = imread(fullfile(fns(i).folder, fns(i).name)) ;
            figure('units', 'normalized', ...
                'outerposition', [0 0 1 1], 'visible', 'off')
            imshow(im * washout2d + max(im) * (1-washout2d))
            xlims = xlim ;
            ylims = ylim ;
            hold on
            v0t2dsc = v0t2d ./ dilation ;
            quiver(x0(:), y0(:), v0t2dsc(:, 1), v0t2dsc(:, 2), 0) ;
            axis equal
            xlim(xlims) ;
            ylim(ylims) ;
            % Extract image from figure axes
            patchIm = getframe(gca);
            outimfn = fullfile(vtdir2d, [fileName '.png']) ;
            % print('-dpng','-r300', outimfn)
            imwrite( patchIm.cdata, outimfn );
            
            
            %% Now draw normal flow as heatmap
            close all; clear alpha 
            % image = im * washout2d + max(im) * (1-washout2d) ;
            image = mat2gray(im, [0, 256]) ;
            % xlims = xlim ;
            % ylims = ylim ;
            % hold on
            % % normalflow = pcolor(x0, y0, reshape(v0n, size(x0)));
            % normalflow = scatter(x0(:), y0(:), 10, v0n, 'filled') ;
            % axis equal
            % colormap( normalflow, bwr);
            % colorbar( normalflow);
            % set( normalflow, 'AlphaData', alpha );   
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
            outimfn = fullfile(vndir2d, [fileName '.png']) ;
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
    save(piv3dfn, 'piv3d') ;
end


%% Rotate and translate to aligned meshes =================================

%% Averaging in time ======================================================
% Make correspondences between faces in t_i, t_{i-1}, and t_{i+1} 

