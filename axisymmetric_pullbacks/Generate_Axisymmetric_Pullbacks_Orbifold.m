%% GENERATE_AXISYMMETRIC_PULLBACK ==========================================
% Pipeline for creating 'axisymmetric' pullbacks of the growing Drosophila 
% midgut about a centerline (previously computed). Uses orbifold method to
% guarantee conformality except at x=0 and x=L
% Run from anywhere, just specify dataDir below.
%
% Uses dorsal points for constructing curve loaded from: 
%    dpFile = fullfile( cylCutDir, 'ap_boundary_dorsalpts.h5' )
% which are saved in slice_mesh_endcaps.m
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

%% ------------------------------------------------------------------------
% Run setup.m from imsane
% -------------------------------------------------------------------------

%% Parameters
overwrite_meshStack = false ;
overwrite_pullbacks = false ;
overwrite_cleanCylMesh = false ;
overwrite_cutMesh = false ;
overwrite_spcutMesh = false ;
generate_sphi_coord = true ;
overwrite_foldims = false ;
overwrite_lobeims = false ;
overwrite_spcutMesh_smoothradii = false ;
resave_ims = false ;
save_ims = true ;
nsegs4path = 5 ;
nV = 100 ;
nU = 100 ;
dvexten = sprintf('_nU%04d_nV%04d', nU, nV) ;
nCurves_yjitter = 100 ;
nCurves_sphicoord = 1000 ;
normal_shift = 10 ;
a_fixed = 2 ;
preview = false ;
washout2d = 0.5 ;
washout3d = 0.5 ;

%% Add paths
% Add some necessary code to the path (ImSAnE should also be setup!) ------
if strcmp(computer, 'MACI64')
    rootdir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/'; 
    gutpath = fullfile(rootdir, 'gut_matlab/') ;
    addpath(genpath(fullfile(gutpath, 'euclidean_orbifolds')));
    addpath(genpath(fullfile(gutpath, 'gptoolbox')));
    addpath_recurse(fullfile(rootdir, 'imsane/external/')) ;
    
    % Setup a working directory for the project, where extracted surfaces,
    % metadata and debugging output will be stored.  Also specifiy the
    % directory containing the data.
    dataDir = ['/Volumes/Pal/48YGal4UASCAAXmCh/201902072000_excellent/', ...
        'Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/'] ;
    cd(dataDir)
else
    gutpath = '/mnt/data/code/gut_matlab/' ;
    addpath(genpath('/mnt/crunch/djcislo/MATLAB/euclidean_orbifolds'));
    addpath(genpath('/mnt/data/code/gptoolbox'));
    addpath_recurse('/mnt/data/code/imsaneV1.2.3/external/') ;
    
    % Setup a working directory for the project, where extracted surfaces,
    % metadata and debugging output will be stored.  Also specifiy the
    % directory containing the data.
    dataDir = [ '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/', ...
        'Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/' ];
    cd(dataDir)
end

addpath(genpath(fullfile(gutpath, 'TexturePatch')));
addpath(gutpath) ;
addpath_recurse(fullfile(gutpath, 'axisymmetric_pullbacks/')) ;
addpath_recurse(fullfile(gutpath, 'plotting/')) ;
addpath_recurse(fullfile(gutpath, 'mesh_handling/')) ;
addpath_recurse(fullfile(gutpath, 'h5_handling/')) ;
addpath_recurse(fullfile(gutpath, 'curve_functions/')) ;
addpath(fullfile(gutpath, 'savgol')) ;
% addpath(genpath('/mnt/crunch/djcislo/MATLAB/TexturePatch'));

%% Define some colors
[colors, color_names] = define_colors() ;
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
cylCutMeshOutImDir = fullfile(cylCutMeshOutDir, 'images') ;
centerlineDir = fullfile(meshDir, 'centerline') ;
centerlineBase = fullfile(centerlineDir, 'mesh_apical_stab_%06d_centerline_exp1p0_res1p0.txt') ;
cntrsFileName = fullfile(centerlineDir, 'mesh_apical_stab_%06d_centerline_scaled_exp1p0_res1p0.txt') ;

% centerline from DV hoops
clineDVhoopDir = fullfile(centerlineDir, ['centerline_from_DVhoops' dvexten]) ;
clineDVhoopBase = fullfile(clineDVhoopDir, 'centerline_from_DVhoops_%06d.mat');
clineDVhoopImDir = fullfile(clineDVhoopDir, 'images') ;
clineDVhoopFigBase = fullfile(clineDVhoopImDir, 'clineDVhoop_%06d.png') ;
radiusDir = fullfile(meshDir, ['radiusDVhoops' dvexten]) ;
radiusImDir = fullfile(radiusDir, 'images') ;

% The file name base for the cylinder meshes
cylinderMeshBase = fullfile( cylCutDir, ...
    'mesh_apical_stab_%06d_cylindercut.ply' );
cylinderMeshCleanBase = fullfile( cylCutMeshOutDir, ...
    'mesh_apical_stab_%06d_cylindercut_clean.ply' );
cylinderMeshCleanFigBase = fullfile( cylCutMeshOutImDir, ...
    'mesh_apical_stab_%06d_cylindercut_clean.png' );

% Lobe identification paths
lobeDir = fullfile(meshDir, 'lobes') ;

% Define cutMesh directories
sphiDir = fullfile(meshDir, ['sphi_cutMesh' dvexten]) ;
sphiBase = fullfile(sphiDir, 'mesh_apical_stab_%06d_spcutMesh.mat') ;
imFolder_sp = [imFolder '_sphi' dvexten] ;
imFolder_sp_e = [imFolder '_sphi' dvexten '_extended'] ;
imFolder_up = [imFolder '_uphi' dvexten] ;
imFolder_up_e = [imFolder '_uphi' dvexten '_extended'] ;
phi0fitBase = fullfile(sphiDir, 'phi0s_%06d.png') ; 

% The file containg the AD/PD points
dpFile = fullfile( cylCutDir, 'ap_boundary_dorsalpts.h5' );

tomake = {imFolder, imFolder_e, imFolder_r, imFolder_re,...
    pivDir, cutFolder, cutMeshImagesDir, cylCutMeshOutDir,...
    cylCutMeshOutImDir, clineDVhoopDir, clineDVhoopImDir, ...
    sphiDir, imFolder_sp, imFolder_sp_e, imFolder_sp, imFolder_sp_e,  ...
    lobeDir, radiusDir, radiusImDir} ;
for i = 1:length(tomake)
    dir2make = tomake{i} ;
    if ~exist( dir2make, 'dir' )
        mkdir(dir2make);
    end
end

%% Load rotation, translation, resolution
rot = dlmread(fullfile(meshDir, 'rotation_APDV.txt')) ;
buff = 20 ;
xyzlim = dlmread(fullfile(meshDir, 'xyzlim_APDV_um.txt'), ',', 1, 0) ;
xyzlim(:, 1) = xyzlim(:, 1) - buff ;
xyzlim(:, 2) = xyzlim(:, 2) + buff ;
xyzlim_APDV = xyzlim ;
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
spmstckfn = fullfile(meshDir, 'spmeshStack_orbifold.mat') ;
outcutfn = fullfile(cutFolder, 'cutPaths_%06d.txt') ;
outadIDxfn = fullfile(cylCutMeshOutDir, 'adIDx.h5') ;
outpdIDxfn = fullfile(cylCutMeshOutDir, 'pdIDx.h5') ;
if exist(mstckfn, 'file') && ~overwrite_meshStack
    % The cutPaths.h5 is loading here
    disp(['Loading meshStack: ' mstckfn])
    load(mstckfn)
    
    if resave_ims
        for t = xp.fileMeta.timePoints
            aux_resave_cutpath_figs
        end
    end
else
    disp('meshStack is not on disk or is to be overwritten, compute...')
    meshStack = cell( length(xp.fileMeta.timePoints), 1 );
    spmeshStack = cell( length(xp.fileMeta.timePoints), 1 );
    for t = xp.fileMeta.timePoints
        disp(['NOW PROCESSING TIME POINT ', num2str(t)]);
        tidx = xp.tIdx(t);

        % Load the data for the current time point ------------------------
        xp.setTime(t) ;
        % Load or compute clean cylindrical mesh
        mesh3dfn =  sprintf( cylinderMeshCleanBase, t ) ;
        if ~exist(mesh3dfn, 'file') || overwrite_cleanCylMesh
            disp('Overwriting/Computing clean cylinderMesh')
            % Load the cylinder mesh
            cylmeshfn = sprintf( cylinderMeshBase, t ) ;
            mesh = read_ply_mod( cylmeshfn );
            mesh = cleanCylMesh(mesh) ;   
            [adIDx, pdIDx] = aux_adjust_dIDx(mesh, t, dpFile, ADBase, PDBase, cylinderMeshCleanBase, outadIDxfn, outpdIDxfn, xp) ;
                        
            %% Save the 3d cut mesh with new indices
            % This is saving the cylinder meshes with no ears. Also adIDx
            % is saved in h5file.
            plywrite_with_normals(mesh3dfn, mesh.f, mesh.v, mesh.vn)
            % Save adIDx with new indices
            save_to_h5(outadIDxfn, ['/' sprintf('%06d', t) ], adIDx, ['adIDx for t=' num2str(t) ' already exists'])
            % Save pdIDx with new indices
            save_to_h5(outpdIDxfn, ['/' sprintf('%06d', t) ], pdIDx, ['pdIDx for t=' num2str(t) ' already exists'])
            disp('done with cylindermesh cleaning')
            
             % View results --------------------------------------------------------
            mesh3dfigfn = sprintf( cylinderMeshCleanFigBase, t ) ;
            if save_ims
                aux_plot_cleanCylMesh
            end
        else
            mesh = read_ply_mod(mesh3dfn) ;
            adIDx = h5read(outadIDxfn, ['/' sprintf('%06d', t)]) ;
            pdIDx = h5read(outpdIDxfn, ['/' sprintf('%06d', t)]) ;

            % View results --------------------------------------------------------
            mesh3dfigfn = sprintf( cylinderMeshCleanFigBase, t ) ;
            if (~exist(mesh3dfigfn, 'file') || overwrite_cleanCylMesh) && save_ims
                aux_plot_cleanCylMesh
            end
        end

        %----------------------------------------------------------------------
        % Create the Cut Mesh
        %----------------------------------------------------------------------
        cutMeshfn = fullfile(cutFolder, [fileNameBase, '_cutMesh.mat']) ;
        cutMeshfn = sprintf(cutMeshfn, t) ;
        if ~exist(cutMeshfn, 'file') || overwrite_cutMesh
            if overwrite_cutMesh
                fprintf('Generating to overwrite cutMesh...')
            else
                fprintf('cutMesh not saved. Generating cutMesh... ');
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            centerline = dlmread(sprintf(centerlineBase, t)) ;
            % try geodesic if first timepoint
            if t == xp.fileMeta.timePoints(1)
                cutOptions.method = 'fastest' ;
                disp(['Cutting mesh using method ' cutOptions.method])
                cutMesh = cylinderCutMesh( mesh.f, mesh.v, mesh.vn, adIDx, pdIDx, cutOptions );
                cutP = cutMesh.pathPairs(:, 1) ;
                adIDx = cutP(1) ;
                pdIDx = cutP(end) ;
                prevTw = twist(mesh.v(cutP, :), centerline) ;
                compute_pullback = true ;
            else 
                % If a previous Twist is not held in RAM, compute it
                % if ~exist('prevTw', 'var')
                % Load previous mesh and previous cutP
                prevcylmeshfn = sprintf( cylinderMeshCleanBase, t-1) ;
                prevmesh = read_ply_mod( prevcylmeshfn ); 
                prevcutP = dlmread(sprintf(outcutfn, t-1), ',', 1, 0) ;
                previousP = prevmesh.v(prevcutP, :) ;
                % Load previous centerline in raw units
                prevcntrfn = sprintf(centerlineBase, t-1) ;
                prevcline = dlmread(prevcntrfn, ',') ;
                % Compute Twist for this previous timepoint
                prevTw = twist(previousP, prevcline) ;
                % end
                [cutMesh, adIDx, pdIDx, cutP, prevTw] = generateCutMeshFixedTwist(mesh, adIDx, pdIDx, centerline, nsegs4path, prevTw, outcutfn, cylinderMeshCleanBase, t) ;
                compute_pullback = true ;                
            end
            
            % Store this path for the next one to be nearby
            % Limit the number of segments to nsegs4path
            % previousP = cutMesh.v(cutP, :) ;
            % pstep = round(length(cutP) / nsegs4path ) ;
            % previousP = previousP(1:pstep:end, :) ;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('Done with generating initial 3D CutMesh with cutPath\n');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Save the cutPath to txt file
            header = 'cutP (path of cut), indexing into vertices' ;
            write_txt_with_header(sprintf(outcutfn, t), cutP, header)  
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Generate pullback to rectangular domain ---------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            fprintf('Relaxing network via Affine transformation... ');
            cutMesh = flattenAnnulus( cutMesh );

            % Find lateral scaling that minimizes spring network energy
            ar = minimizeIsoarealAffineEnergy( cutMesh.f, cutMesh.v, cutMesh.u );
            % Assign scaling based on options: either a0 or a_fixed
            % if tidx == 1 && ~a_fixed
            %     a_fixed = ar ;
            % end      
            % a = a_fixed ;
      
            % Scale the x axis by a or ar
            uvtx = cutMesh.u ;
            cutMesh.u = [ a .* uvtx(:,1), uvtx(:,2) ];
            cutMesh.urelax = [ ar .* uvtx(:,1), uvtx(:,2) ];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('Done flattening cutMesh. Now saving.\n');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Save cutMesh
            save(cutMeshfn, 'cutMesh', 'adIDx', 'pdIDx', 'cutP')
            disp('done with plotting & saving cut')
        else
            fprintf('Loading Cut Mesh from disk... ');
            load(cutMeshfn) 
            cutP = dlmread(sprintf(outcutfn, t), ',', 1, 0) ;
            compute_pullback = ~isempty(cutP) ;
        end
        
        %% Plot the cutPath (cutP) in 3D
        if save_ims && overwrite_cutMesh 
            disp('Saving cutP image')
            aux_plot_cutP
        end
        
        %% Compute the pullback if the cutMesh is ok
        if compute_pullback
            disp(['Evolving mesh along normal shift for pullback images: shift=' num2str(normal_shift)])
            % Displace normally ---------------------------------------------------
            cutMesh.v = cutMesh.v + cutMesh.vn * normal_shift ;
            
            % todo !!!!
            % if 
            % else
            %     % Load the cutMesh
            %     load(sprintf(cutMeshBase, t), 'cutMesh')
            % end
            
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Generate s,phi coord system for rotated,scaled mesh (rs)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('Establishing s,phi coord system\n');

            if ~exist(sprintf(sphiBase, t), 'file') || overwrite_spcutMesh
                if overwrite_spcutMesh
                    disp('Overwriting spcutMesh...')
                else
                    disp('spcutMesh not on disk. Generating ...')
                end

                % Transform from u,v coordinates to s, phi coordinates
                % [scoords, phicoords] = generateSPhiFromUV();

                %----------------------------------------------------------------------
                % Generate tiled orbifold triangulation
                %----------------------------------------------------------------------
                tileCount = [1 1];  % how many above, how many below
                cutMeshrs = cutMesh;
                % Rotate and translate TV3D
                cutMeshrs.v = ((rot * cutMesh.v')' + trans) * resolution ;
                cutMeshrs.vn = (rot * cutMesh.vn')' ;
                [ ~, ~, TV3D ] = tileAnnularCutMesh( cutMesh, tileCount );
                [ TF, TV2D, TV3Drs ] = tileAnnularCutMesh( cutMeshrs, tileCount );

                %----------------------------------------------------------------------
                % Calculate abbreviated centerline from cutMesh boundaries
                %----------------------------------------------------------------------
                % Load centerline in raw units
                cntrfn = sprintf(cntrsFileName, t) ;
                cline = dlmread(cntrfn, ',') ;
                ss = cline(:, 1) ;
                cline = cline(:, 2:end) ;

                % Check it
                % trisurf(triangulation(TF, TV3D), 'EdgeColor', 'none', 'FaceAlpha', 0.3)
                % plot3(cline(:, 1), cline(:, 3), cline(:, 2), 'k-')
                % set(gcf, 'visible', 'on')

                disp('Finding relevant segment of centerline')
                [cseg, acID, pcID, bdLeft, bdRight] = centerlineSegmentFromCutMesh(cline, TF, TV2D, TV3Drs) ;

                %----------------------------------------------------------------------
                % Generate surface curves of constant s
                %----------------------------------------------------------------------
                % For lines of constant phi
                disp('Creating constant y curves')
                % Make grid
                eps = 1e-14 ;
                uspace = linspace( eps, max(TV2D(:, 1)) - eps, nU )' ;
                vspace = linspace( eps, 1-eps, nV )' ;

                disp('Casting points into 3D...')
                % NOTE: first dimension indexes u, second indexes v
                curves3d = zeros(nU, nV, 3) ;
                for kk = 1:nU
                    if mod(kk, 50) == 0
                        disp(['u = ' num2str(kk / nU)])
                    end
                    uv = [uspace(kk) * ones(size(vspace)), vspace] ;
                    curves3d(kk, :, :) = interpolate2Dpts_3Dmesh(TF, TV2D, TV3Drs, uv) ;
                end 

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf('Compute s(u) and radius(u), each (0,1) \n');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % todo: sample at evenly spaced dphi in embedding space
                fprintf('Resampling curves...\n')
                c3ds = zeros(size(curves3d)) ;
                for i=1:nU
                    % Note: no need to add the first point to the curve
                    % since the endpoints already match exactly in 3d and
                    % curvspace gives a curve with points on either
                    % endpoint (corresponding to the same 3d location).
                    c3ds(i, :, :) = resampleCurvReplaceNaNs(squeeze(curves3d(i, :, :)), nV, true) ;
                    if vecnorm(squeeze(c3ds(i, 1, :)) - squeeze(c3ds(i, end, :))) > 1e-7
                        error('endpoints do not join! Exiting')
                    end

                    % Visualization for Troubleshooting:
                    % triplot(TF, TV2D(:, 1), TV2D(:, 2))
                    % hold on;
                    % plot(uv(:, 1), uv(:, 2), '.')
                end

                fprintf('Finding s(u) and r(u) of resampled c3ds...\n')
                % mcline is the resampled centerline, with mss
                % avgpts is the raw Nx3 averaged hoops, with avgpts_ss
                [mss, mcline, radii_from_mean_uniform_rs, avgpts_ss, avgpts] = srFromDVCurves(c3ds) ;
                
                % Used to find radius using original centerline
                % [ssv, radii, avgpts, cids] = srFromDVCurvesGivenCenterline(ss, cline, c3ds) ;
                % Could operate just on the centerline segment
                cseg_ss = ss(acID:pcID) ;
                % [ssv, radii, avgpts, cids] = srFromDVCurves(cseg_ss, cseg, c3ds) ;
                % 
                % Adjust the centerline indices to index into the full
                % centerline. Note that cseg_ss already does this for ss.
                % cids = cids + acID ;

                % Plot new centerline
                aux_plot_clineDVhoop(avgpts, avgpts_ss, cseg, cline, cseg_ss, curves3d, xyzlim, clineDVhoopFigBase, t)

                % Optional: clean curve with polynomial and point match
                % avgpts onto cleaned curve. Skipping for later.

                % Compute ringpath_ss, the mean distance traveled from one
                % line of constant u to the next
                disp('Computing ringpath_ss...')
                % The distance from one hoop to another is the
                % difference in position from (u_i, v_i) to (u_{i+1}, v_i).
                dsuphi = reshape([vecnorm(diff(c3ds), 2, 2); 0], [nU, nV]) ;
                ringpath_ds = nanmean(dsuphi(1:(end-1), :), 2) * resolution ;
                ringpath_ss = cumsum([0; ringpath_ds]) ;
                clearvars dsuphi ringpath_ds
                
                % Save new centerline in rotated translated units
                fn = sprintf(clineDVhoopBase, t) ;
                disp(['Saving new centerline to ' fn])
                save(fn, 'mss', 'mcline', 'avgpts', 'avgpts_ss')
                
                % Note: radii_from_mean_uniform_rs is the radius of 
                % interpolated hoops, not the actual points
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf('Done with constructing new centerline\n') ;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if t == xp.fileMeta.timePoints(1)
                    % Store for next timepoint
                    phiv = (vspace .* ones(nU, nV))' ;
                    phi0s = zeros(size(uspace)) ;
                    phi0_fit = phi0s ;
                else
                    % Load previous sphi vertices in 3d if not in RAM
                    % if ~exist('prev3d_sphi', 'var')
                    tmp = load(sprintf(sphiBase, t-1), 'spcutMesh') ;
                    prev3d_sphi = reshape(tmp.spcutMesh.v, [nU, nV, 3]) ; 
                    
                    plotfn = sprintf(phi0fitBase, t);
                    [phi0_fit, phi0s] = fitPhiOffsetsFromPrevMesh(TF, TV2D, TV3D, ...
                        uspace, vspace, prev3d_sphi, -0.25, 0.25, save_ims, plotfn) ;
                    % Store for next timepoint
                    phiv = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;
                end
                disp('done computing phi0s')

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                onesUV = ones(nU, nV) ;
                uu = uspace .* onesUV ;
                vv = (vspace .* onesUV')' ;
                uphi = [uu(:), phiv(:)] ;
                uv = [uu(:), vv(:)] ;
                new3d = interpolate2Dpts_3Dmesh(TF, TV2D, TV3D, uphi) ;
                % Express the coordinates as a grid
                % new3drs = interpolate2Dpts_3Dmesh(TF, TV2D, TV3Drs, uphi) ;
                prev3d_sphi = reshape(new3d, [nU, nV, 3]) ; 
                
                % Recompute radii_from_mean_uniform_rs as radii_from_avgpts 
                % NOTE: all radius calculations done in microns, not pixels
                sphi3d_rs = ((rot * prev3d_sphi')' + trans) * resolution ;
                radii_from_avgpts = zeros(size(sphi3d_rs, 1), size(sphi3d_rs, 2)) ;
                for jj = 1:nU
                    % Consider this hoop
                    hoop = squeeze(sphi3d_rs(jj, :, :)) ;
                    radii_from_avgpts(jj, :) = vecnorm(hoop - avgpts(jj, :), 2, 2) ;
                end
                
                % Triangulate the sphigrid and store as its own cutMesh
                % sphiv = zeros(nU, nV, 2) ;
                % sphiv(:, :, 1) = sv ;
                % sphiv(:, :, 2) = phiv ;
                sv = ringpath_ss .* onesUV ;
                % Triangulate the mesh
                tmptri = delaunay(sv(:), phiv(:)) ;
                disp('orienting faces of delaunay triangulation (s,phi)')
                tmptri = bfs_orient( tmptri );

                % Define path pairs for tiling the (s,phi) cut mesh
                spcutP1 = 1:nU;
                spcutP2 = nU*nV - fliplr(0:(nU-1)) ;
                spcutMesh.pathPairs = [ spcutP1', spcutP2' ];

                % Check to see if any members of pathPairs connect to
                % non-Nearest Neighbors.
                cleantri = cleanBoundaryPath2D(tmptri, [sv(:), phiv(:)], spcutMesh.pathPairs(:), true) ;

                spcutMesh.f = cleantri ;
                spcutMesh.v = new3d ;
                % spcutMesh.vrs = ((rot * new3d')' + trans) * resolution ;
                spcutMesh.nU = nU ;
                spcutMesh.nV = nV ;
                % Define normals based on the original mesh normals
                % todo: spcutMesh.vn =  ;
                Fnx = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), cutMesh.vn(:, 1)) ;
                Fny = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), cutMesh.vn(:, 2)) ;
                Fnz = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), cutMesh.vn(:, 3)) ;
                spcutMesh.vn = zeros(size(spcutMesh.v)) ;
                spcutMesh.vn(:, 1) = Fnx(uphi(:, 1), uphi(:, 2)) ;
                spcutMesh.vn(:, 2) = Fny(uphi(:, 1), uphi(:, 2)) ;
                spcutMesh.vn(:, 3) = Fnz(uphi(:, 1), uphi(:, 2)) ;
                spcutMesh.sphi = [sv(:), phiv(:)] ;
                spcutMesh.uphi = uphi ;
                % Note: uv has no direct relation with cutMesh, just a grid for utility
                spcutMesh.uv = uv ;  
                spcutMesh.ringpath_ss = ringpath_ss ;
                spcutMesh.radii_from_mean_uniform_rs = radii_from_mean_uniform_rs ;
                spcutMesh.radii_from_avgpts = radii_from_avgpts ;
                spcutMesh.mss = mss ;       % from uniform sampling, also stored in centerline
                spcutMesh.mcline = mcline ; % from uniform sampling, also stored in centerline
                spcutMesh.avgpts = avgpts ; % from uniform sampling, also stored in centerline
                spcutMesh.avgpts_ss = avgpts_ss ; % from uniform sampling, also stored in centerline

                % Save s,phi and their 3D embedding
                spcutMesh.phi0s = phi0s ;
                spcutMesh.phi0_fit = phi0_fit ;
                save(sprintf(sphiBase, t), 'spcutMesh') ;
            else
                disp('Loading spcutMesh from disk...')
                load(sprintf(sphiBase, t), 'spcutMesh') ;

                % Load new centerline
                fn = sprintf(clineDVhoopBase, t) ;
                disp(['Loading new centerline from ' fn])
                load(fn, 'mss', 'mcline', 'avgpts', 'avgpts_ss')
            end
            fprintf('Done with generating S,Phi coords \n');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('Create pullback using S,Phi coords \n');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %--------------------------------------------------------------
            % Generate Output Image File
            %--------------------------------------------------------------
            imfn = sprintf( fullfile([imFolder, '/', fileNameBase, '.tif']), t ); 
            imfn_r = sprintf( fullfile([imFolder_r, '/', fileNameBase, '.tif']), t ) ;
            imfn_sp = sprintf( fullfile([imFolder_sp, '/', fileNameBase, '.tif']), t ) ;
            imfn_up = sprintf( fullfile([imFolder_up, '/', fileNameBase, '.tif']), t ) ;
            pullbacks_exist1 = exist(imfn, 'file') && exist(imfn_r, 'file') ;
            pullbacks_exist2 = exist(imfn_sp, 'file') && exist(imfn_up, 'file') ;
            if ~pullbacks_exist1 || ~pullbacks_exist2 || overwrite_pullbacks
                % Load 3D data for coloring mesh pullback
                xp.loadTime(t);
                xp.rescaleStackToUnitAspect();

                % Raw stack data
                IV = xp.stack.image.apply();
                IV = imadjustn(IV{1});
            end
            
            if ~exist(imfn_sp, 'file') || overwrite_pullbacks
                fprintf(['Generating SP output image: ' imfn_sp]);
                % Assigning field spcutMesh.u to be [s, phi] (ringpath
                % and azimuthal angle)
                spcutMesh.u = spcutMesh.sphi ;
                aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_sp)
                spcutMesh = rmfield(spcutMesh, 'u') ;
            end
            
            if ~exist(imfn_up, 'file') || overwrite_pullbacks
                fprintf(['Generating uphi output image: ' imfn_up]);
                % Assigning field spcutMesh.u to be [s, phi] (ringpath
                % and azimuthal angle)
                spcutMesh.u = spcutMesh.uphi ;
                aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_up)
                spcutMesh = rmfield(spcutMesh, 'u') ;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Generate Output Image File -- regular UV coordinates
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~exist(imfn, 'file') || overwrite_pullbacks
                % Generate output image in uv
                fprintf(['Generating output image: ' imfn]);
                aux_generate_orbifold(cutMesh, a_fixed, IV, imfn)
            else
                disp('Skipping pullback image generation since exists')
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Save relaxed image
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            if ~exist(imfn_r, 'file')
                disp('Generating relaxed image...')
                aux_generate_orbifold(cutMesh, ar, IV, imfn_r)
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Save submesh array. Each cell element contains all the 
            % submeshes for that TP, which in this case is just one.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            meshStack{tidx} = cutMesh ;
            if generate_sphi_coord
                spmeshStack{tidx} = spcutMesh ;
            end
            
            clear Options IV

            fprintf('Done\n');
        else
            disp('Skipping computation of pullback')
        end
    end

    %% Save SMArr2D (vertex positions in the 2D pullback) -----------------
    disp(['Saving meshStack to disk: ' mstckfn])
    save(mstckfn, 'meshStack') ;
    
    %% Save SMArr2D (vertex positions in the 2D pullback) -----------------
    disp(['Saving spmeshStack to disk: ' spmstckfn])
    if generate_sphi_coord
        save(spmstckfn, 'spmeshStack') ;
    end
end

%% Preview results
check = false ;
if check
    aux_preview_results
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TILE IMAGES IN Y AND RESAVE ============================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imDirs = {imFolder, imFolder_sp, imFolder_up} ;
imDirs_e = {imFolder_e, imFolder_sp_e, imFolder_up_e} ;
for qq = 1:2
    extendImages(imDirs{qq}, imDirs_e{qq}, fileNameBase)
    disp(['done ensuring extended tiffs for ' directory ' in ' direc_e])
end
disp('done')

%% FIND THE FOLDS SEPARATING COMPARTMENTS =================================
% First compute using the avgpts (DVhoop means)
guess123 = [0.2, 0.5, 0.8] ;
max_wander = 20 ;
disp('Identifying lobes...')
foldfn = fullfile(lobeDir, ['fold_locations_sphi' dvexten '_avgpts.mat']) ;
if exist(foldfn, 'file')
    % Save the fold locations as a mat file
    load(foldfn, 'ssfold', 'folds', 'ssfold_frac', 'ssmax', 'fold_onset')
else
    [folds, ssfold, ssfold_frac, ssmax, rmax, fold_onset] = identifyLobes(xp.fileMeta.timePoints,...
            sphiBase, guess123, max_wander, preview, 'avgpts') ;
    
    % Compute ringpath pathlength for results found using centerline
    disp('Converting folds to ringpath_ss locations...')
    [rssfold, rssfold_frac, rssmax] = rssFromFoldID(folds, xp.fileMeta.timePoints, sphiBase) ;

    % Save the fold locations as a mat file
    save(foldfn, 'rssfold', 'rssfold_frac', 'rssmax', ...
        'ssfold', 'folds', 'ssfold_frac', 'ssmax', 'fold_onset')

    % Plot results as both avgpts and ringpath distances
    if save_ims
        disp('Plotting ss folds...')
        aux_plot_folds(folds, ssfold, ssfold_frac, ssmax, rmax, nU, ...
            xp.fileMeta.timePoints, lobeDir, dvexten, sphiBase, ...
            'avgpts', overwrite_foldims)
        disp('Plotting rss folds...')
        aux_plot_folds(folds, rssfold, rssfold_frac, rssmax, rmax, nU, ...
            xp.fileMeta.timePoints, lobeDir, dvexten, sphiBase, ...
            'ringpath', overwrite_foldims)
    end
end
disp('done')

%% Compute surface area and volume for each compartment
lobe_dynamics_fn = fullfile(lobeDir, ['lobe_dynamics' dvexten '.mat']) ;
if exist(lobe_dynamics_fn, 'file')
    % Load length, surface area, and volume dynamics for each lobe
    disp('Loading lobe length, area, and volume...')
    load(lobe_dynamics_fn, 'length_lobes', 'area_lobes', 'volume_lobes')
    tp = xp.fileMeta.timePoints - min(fold_onset) ;
else
    aux_compute_lobe_dynamics(ssfold, ssmax, lobeDir, timePoints, ...
        sphiBase, rot, trans, xyzlim, colors) 
    % Save surface area and volume dynamics for each lobe
    save(fullfile(lobeDir, ['lobe_dynamics' dvexten '.mat']), ...
        'length_lobes', 'area_lobes', 'volume_lobes')
end


%% plot length, area, and volume for each lobe
lobe_dynamics_figfn = fullfile(lobeDir, ['lobe_dynamics' dvexten '.png']) ;
if save_ims && (~exist(lobe_dynamics_figfn, 'file') || overwrite_lobefigs)
    close all;
    fig = figure('visible', 'off');
    fh = cell(1, 4) ;
    tp = timePoints - min(fold_onset) ;
    for lobe = 1:4
        % Length
        subplot(3,2,1)
        plot(tp, length_lobes(:, lobe), '.', 'Color', colors(lobe, :)) ;
        hold on

        % Area
        subplot(3,2,3)
        fh{lobe} = plot(tp, area_lobes(:, lobe), '.', 'Color', colors(lobe, :)) ;
        hold on

        % Volume
        subplot(3,2,5)
        plot(tp, volume_lobes(:, lobe), '.', 'Color', colors(lobe, :)) ;
        hold on
    end

    subplot(3, 2, 1)
    xlim([min(tp), max(tp)])
    ylabel('Length [\mum]')

    subplot(3, 2, 3)
    xlim([min(tp), max(tp)])
    ylabel('Area [\mum^2]')
    legend({'\color[rgb]{ 0,0.4470,0.7410} lobe 1', ...
        '\color[rgb]{0.8500,0.3250,0.0980} lobe 2', ...
        '\color[rgb]{0.9290,0.6940,0.1250} lobe 3', ...
        '\color[rgb]{0.4940,0.1840,0.5560} lobe 4'}, 'Position', [0.55 0.4 0.1 0.2])

    subplot(3, 2, 5)
    xlim([min(tp), max(tp)])
    ylabel('Volume [\mum^3]')

    xlabel('time [min]')
    disp(['Saving summary to ' lobe_dynamics_figfn])
    saveas(fig, lobe_dynamics_figfn)
end

%% Plot motion of avgpts at lobes in yz plane over time
fold_dynamics_figfn = fullfile(lobeDir, ['constriction_dynamics' dvexten '.png']) ;
if save_ims && (~exist(fold_dynamics_figfn, 'file') || overwrite_lobefigs)
    f1pts = zeros(length(tp), 3) ;
    f2pts = zeros(length(tp), 3) ;
    f3pts = zeros(length(tp), 3) ;

    disp('Loading centerline position of each fold...')
    for kk = 1:length(xp.fileMeta.timePoints)
        % Translate to which timestamp
        t = xp.fileMeta.timePoints(kk) ;
        load(sprintf(sphiBase, t), 'spcutMesh') ;

        % Load the centerline too
        avgpts = spcutMesh.avgpts ;

        % rename the fold indices (in U)
        f1 = folds(kk, 1) ;
        f2 = folds(kk, 2) ;
        f3 = folds(kk, 3) ;

        % store distance from x axis of folds
        f1pts(kk, :) = avgpts(f1, :) ;
        f2pts(kk, :) = avgpts(f2, :) ;
        f3pts(kk, :) = avgpts(f3, :) ;
    end
    close all
    alph = 0.2 ;
    cmap = colormap ;
    fig = figure('visible', 'off'); 
    hold on;
    for qq = 1:length(tp)
        t = xp.fileMeta.timePoints(qq) ;
        % Load the centerline
        fn = sprintf(clineDVhoopBase, t) ;
        load(fn, 'avgpts')

        color = cat(2, cmap(uint8(max(1, qq * length(cmap)/length(tp))), :), alph);
        % plot([f1pts(qq, 2), f2pts(qq, 2), f3pts(qq, 2)], ...
        %     [f1pts(qq, 3), f2pts(qq, 3), f3pts(qq, 3)], '-', 'color', color)
        plot(avgpts(:, 2), avgpts(:, 3), '-', 'color', color)
    end
    sz = 40 ;
    msz = 20 ;
    scatter(f1pts(:, 2), f1pts(:, 3), sz, tp, 'o'); hold on;
    scatter(f2pts(:, 2), f2pts(:, 3), sz, tp, 's');
    scatter(f3pts(:, 2), f3pts(:, 3), sz, tp, '^');
    idx = zeros(max(xp.fileMeta.timePoints), 1) ;
    idx(xp.fileMeta.timePoints) = 1:length(tp) ;
    plot(f1pts(idx(fold_onset(1)), 2), f1pts(idx(fold_onset(1)), 3), 'ko', 'markersize', msz); 
    plot(f2pts(idx(fold_onset(2)), 2), f2pts(idx(fold_onset(2)), 3), 'ks', 'markersize', msz);
    plot(f3pts(idx(fold_onset(3)), 2), f3pts(idx(fold_onset(3)), 3), 'k^', 'markersize', msz);
    axis equal
    xlabel('y [\mum]')
    ylabel('z [\mum]')
    title('Constriction dynamics')
    cb = colorbar() ;
    cb.Label.String = 'time [min]' ;
    saveas(fig, fold_dynamics_figfn) ;
    close all
end

%% SMOOTH MEAN CENTERLINE RADIUS ==========================================
mclineM = zeros(length(xp.fileMeta.timePoints), 1000, 3) ;
for kk=1:length(xp.fileMeta.timePoints)
    t = xp.fileMeta.timePoints(kk) ;
    % Load mcline
    load(sprintf(clineDVhoopBase, t), 'mcline') ;
    mclineM(kk, :, :) = mcline ;
end

% Time-average the centerline
mcline_sm = movmean(mclineM,5,1) ;

% Optional: polynomial fit to smooth

% Re-compute radii from smoothed centerlines by pointmatching avgpts onto
% smoothed centerline and plot them
radImDir0 = fullfile(radiusImDir, 'view0') ;
radImDirD = fullfile(radiusImDir, 'viewD') ;
radImDirV = fullfile(radiusImDir, 'viewV') ;
radImDirA = fullfile(radiusImDir, 'viewA') ;
radImDirP = fullfile(radiusImDir, 'viewP') ;
radImDirL = fullfile(radiusImDir, 'viewL') ;
radImDirR = fullfile(radiusImDir, 'viewR') ;
radImDirs = fullfile(radiusImDir, 'sphi') ;
radImDiru = fullfile(radiusImDir, 'uphi') ;
dirs2do = {radImDir0, radImDirD, radImDirV, ...
    radImDirA, radImDirP, radImDirL, radImDirR, ...
    radImDirs, radImDiru} ;
for qq = 1:length(dirs2do)
    if ~exist(dirs2do{qq}, 'dir')
        mkdir(dirs2do{qq})
    end
end
% Now compute and plot
for kk=1:length(xp.fileMeta.timePoints)
    t = xp.fileMeta.timePoints(kk) ;
    % Load mcline
    load(sprintf(clineDVhoopBase, t), 'avgpts', 'mcline', 'mss') ;
    mcs_kk = squeeze(mcline_sm(kk, :, :)) ;
    mss_sm = ss_from_xyz(mcs_kk) ;
    inds = pointMatch(avgpts, mcs_kk) ;
    avgpts_projected = mcs_kk(inds, :) ;
    
    % Save the new results with the old
    save(sprintf(clineDVhoopBase, t), 'avgpts', 'mcline', 'mss',...
        'avgpts_projected', 'mcline_sm', 'mss_sm') ;
    
    % Add radii from smoothed mcline to spcutMesh if not already present
    load(sprintf(sphiBase, t), 'spcutMesh') ;
    vrs = ((rot * spcutMesh.v')' + trans) * resolution ;
    radii_in_cutMeshsm = isfield(spcutMesh, 'radii_from_smoothed_mcline') ;
    if ~radii_in_cutMeshsm || overwrite_spcutMesh_smoothradii
        curvesDV = reshape(vrs, [nU, nV, 3]) ;
        radii_from_avgpts_sm = zeros(size(curvesDV, 1), size(curvesDV, 2)) ;
        for jj = 1:size(curvesDV, 1)
            % Consider this hoop
            hoop = squeeze(curvesDV(jj, :, :)) ;
            radii_from_avgpts_sm(jj, :) = vecnorm(hoop - avgpts_projected(jj, :), 2, 2) ;
        end
    else
        radii_from_avgpts_sm = spcutMesh.radii_from_smoothed_mcline ;
    end
    
    % Plot the radii as 3D image
    rad2dfigfn_s = fullfile(radImDirs, sprintf('radius_dvhoop_sphi_%06d.png', t)) ;
    rad2dfigfn_u = fullfile(radImDiru, sprintf('radius_dvhoop_uphi_%06d.png', t)) ;
    rad3dfigfn = fullfile(radImDir0, sprintf('radius_dvhoop_xyz_%06d.png', t)) ;
    if ~exist(rad3dfigfn, 'file') || ~exist(rad2dfigfn_s, 'file') || ...
            ~exist(rad2dfigfn_u, 'file') || overwrite_spcutMesh_smoothradii
        close all
        fig = figure('visible', 'off') ;
        trisurf(spcutMesh.f, vrs(:, 1), vrs(:, 2), ...
            vrs(:, 3), radii_from_avgpts_sm(:), ...
            'Edgecolor', 'none')
        c = colorbar ;
        c.Label.String = 'radius [\mum]' ;
        xlabel('x [\mum]')
        ylabel('y [\mum]')
        zlabel('z [\mum]')
        axis equal
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        title('Radius via DV curves and smoothed centerline')
        saveas(fig, rad3dfigfn)
        title('Radius via DV curves, Dorsal view')
        view(0, 90)
        fn = sprintf('radius_dvhoop_xyD_%06d.png', t) ;
        disp(['Saving ' fn])
        saveas(gcf, fullfile(radImDirD, fn)) 
        title('Radius via DV curves, Ventral view')
        view(0, -90)
        fn = sprintf('radius_dvhoop_xyV_%06d.png', t) ;
        disp(['Saving ' fn])
        saveas(gcf, fullfile(radImDirV, fn)) 
        title('Radius via DV curves, Posterior view')
        view(90, 0)
        fn = sprintf('radius_dvhoop_yzP_%06d.png', t) ;
        disp(['Saving ' fn])
        saveas(gcf, fullfile(radImDirP, fn)) 
        title('Radius via DV curves, Anterior view')
        view(-90, 0)
        fn = sprintf('radius_dvhoop_yzA_%06d.png', t) ;
        disp(['Saving ' fn])
        saveas(gcf, fullfile(radImDirA, fn)) 
        title('Radius via DV curves, Lateral view')
        view(0, 0)
        fn = sprintf('radius_dvhoop_xzL_%06d.png', t) ;
        disp(['Saving ' fn])
        saveas(gcf, fullfile(radImDirL, fn)) 
        title('Radius via DV curves, Lateral view')
        view(0, 180)
        fn = sprintf('radius_dvhoop_xzR_%06d.png', t) ;
        disp(['Saving ' fn])
        saveas(gcf, fullfile(radImDirR, fn)) 

        % Plot the radii as 2D image
        tmp = {spcutMesh.sphi, spcutMesh.uphi} ;
        rad2dfigfns = {rad2dfigfn_s, rad2dfigfn_u} ;
        xlabels = {'AP position, s/L', 'AP position, u/L'} ;
        for qq=1:2
            uu = tmp{qq} ;
            if qq == 1
                uu(:, 1) = uu(:, 1) / max(uu(:, 1)) ;
            else
                uu(:, 1) = 0.5 * uu(:, 1) ;
            end
            close all
            fig = figure('visible', 'off') ;
            trisurf(spcutMesh.f, uu(:, 1), uu(:, 2), ...
                radii_from_avgpts_sm(:), 'EdgeColor', 'none')
            caxis([0, 80]) ;
            % also plot tiled meshes above and below
            hold on;
            trisurf(spcutMesh.f, uu(:, 1), uu(:, 2) + 1, ...
                radii_from_avgpts_sm(:), 'Edgecolor', 'none')
            trisurf(spcutMesh.f, uu(:, 1), uu(:, 2) - 1, ...
                radii_from_avgpts_sm(:), 'Edgecolor', 'none')
            c = colorbar ;
            c.Label.String = 'radius [\mum]' ;
            xlabel(xlabels{qq})
            ylabel('\phi/2\pi')
            ylim([0, 1])
            view(2)
            saveas(fig, rad2dfigfns{qq})
            close all
        end
    end
end
disp('done')
% note: see extract_radius_from_DVhoops for old version of radius plotting



%% PERFORM PIV ============================================================
% Select all frames in PullbackImages_extended_shifted/ 
% Select Sequencing style 1-2, 2-3, ... 
% Image Preprocessing > Select All
% PIV settings: 128 (32 step), 64 (16 step) for two passes
% Post-processing: vector validation: Select velocity limits, 5 stdev
% filter, local median filter: thres=2, eps=0.1
% with default values.
% File > Save > Export all frames as pivresults.mat
disp('Loading pivresults_sphi_pass0.mat ...')
piv = load(fullfile(pivDir, 'pivresults_sphi_pass0.mat')) ;

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
disp('done building dt, meshidx') 
save(fullfile(pivDir, 'timestamps_orbifold.mat'), 'time', 'meshidx')
% 
%% Subtract off the mean flow in y for each frame ========================= 
% meanv = zeros(length(v_filtered), 1) ;
% for i=1:length(v_filtered)
%     tmp = v_filtered{i} ;
%     meanv(i) = mean(tmp(:)) ;
% end
% dy = round(meanv) ;
% shifty = [0; -cumsum(dy)] ;
% disp('computed shifty')
% save(fullfile(pivDir, 'shifty.mat'), 'shifty')

%% Shift each frame by shifty and resave ==================================
% fns = dir(strrep(fullfile([imFolder_sp_e, '/', fileNameBase, '.tif']), '%06d', '*')) ;
% imFolder_es = [imFolder_sp_e '_shifted' filesep] ;
% if ~exist(imFolder_es, 'dir')
%     mkdir(imFolder_es) ;
% end
% 
% % Write the images to disk and get their sizes while being written 
% disp('Resaving shifted images if necessary...')
% imsizes = zeros(length(fns), 2) ;
% for i=1:length(fns)
%     outfn = fullfile(imFolder_sp_es, fns(i).name) ;
%     if ~exist(outfn, 'file')
%         disp(['Reading ' fns(i).name])
%         fileName = split(fns(i).name, '.tif') ;
%         fileName = fileName{1} ;
%         im = imread(fullfile(fns(i).folder, fns(i).name)) ;
%         % Note adaptive histogram equilization
%         im = adapthisteq(circshift(im, shifty(i), 1)) ;
%         imsizes(i, :) = size(im) ;
%         imwrite( im, outfn, 'TIFF' );
%     else
%         im = imread(outfn) ;
%         imsizes(i, :) = size(im) ;
%     end
% end
% disp('done with shifted images')

%% PERFORM PIV ON SHIFTED FRAMES ==========================================
% % Select all frames in PullbackImages_extended_shifted/ 
% % Select Sequencing style 1-2, 2-3, ... 
% % Load settings: piv_set_pass1.mat
% % Image Preprocessing > Select All
% % PIV settings: 128 (32 step), 64 (16 step), 32 (16 step) for three passes
% disp('Loading pivresults_orbifold_pass1.mat ...')
% load(fullfile(pivDir, 'pivresults_orbifold_pass1.mat'))



%% Smooth velocities in time
% Interpolate velocities on grid to smooth them in time
piv3dfn = fullfile(pivDir, 'piv3d.mat') ;
if exist(piv3dfn, 'file')
    load(piv3dfn)
else
    piv3d = cell(length(fns), 1) ;
    for i=1:length(fns) - 1
        % Average 
        x0 = x{i} ;
        y0 = y{i} ;
        uu = u_filtered{i} ;
        vv = v_filtered{i} ; 
        Fx = griddedInterpolant(x0, y0, xvel) ;
        Fy = griddedInterpolant(x0, y0, yvel) ;
        Fz = griddedInterpolant(x0, y0, zvel) ;
        % query the velocity at the advected locations for each PIV gridpt
    end
end


%% MAKE MAP FROM PIXEL TO XYZ =============================================
disp('Making map from pixel to xyz to compute velocities in 3d...')
% Get position of embedding points associated with velocities 
pivOutDir = fullfile(pivDir, 'images') ;
if ~exist(pivOutDir, 'dir')
    mkdir(pivOutDir)
end
use_shifty = false ;

% Compute size of domain (0,1), (0, 1)
if use_shifty
    fns = dir(strrep(fullfile([imFolder_sp_es, '/', fileNameBase, '.tif']), '%06d', '*')) ;
else
    fns = dir(strrep(fullfile([imFolder_sp_e, '/', fileNameBase, '.tif']), '%06d', '*')) ;
end

% For now, assume that all images are the same size
% % Get imsizes
% imsizes = zeros(length(fns), 2) ;
% for i=1:length(fns)
%     im = imread(fullfile(fns(i).folder, fns(i).name)) ;
%     imsizes(i, :) = size(im) ;
% end

im = imread(fullfile(fns(i).folder, fns(i).name)) ;
% size of extended image
esize = size(im) ;
% map extended image size to (0, 1), (-0.5, 1.5)
xesz = esize(2) ;
yesz = esize(1) ;
% map from pixel y to network y (sphi)
Ypix2y = @(ypix) 2.0 * ypix / yesz - 0.5 ;
Xpix2x = @(xpix) 2.0 * xpix / yesz ;
% map from network y to pixel y
% Note that for some reason we need to flip the image (Yscale - stuff)
x2Xpix = @(x, Yscale, xscale) (Yscale * 0.5) * x / xscale ;
dx2dX = @ (y, Yscale, xscale) (Yscale * 0.5) * x / xscale ;
dy2dY = @ (y, Yscale) (Yscale*0.5)*y ;
if use_shifty
    y2Ypix = @(y, h, Yscale) Yscale - (Yscale*0.5)*(y+0.5) + h ;
else
    y2Ypix = @(y, Yscale) Yscale - (Yscale*0.5)*(y+0.5) ;
end

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
        Ysc0 = yesz ;  % imsizes(i, 2) ;
        Ysc1 = yesz ;  % imsizes(i+1, 2) ;

        % Load spcutMesh
        mesh0 = load(sprintf(sphiBase, time(i)), 'spcutMesh') ;
        mesh0 = mesh0.spcutMesh ;
        mesh1 = load(sprintf(sphiBase, time(i + 1)), 'spcutMesh') ;
        mesh1 = mesh1.spcutMesh ;
        assert(time(i) + dt(i) == time(i+ 1))
        xsc0 = max(mesh0.sphi(:, 1)) ;
        xsc1 = max(mesh0.sphi(:, 1)) ;

        % Load the positions of the velocity vectors in pixels
        x0 = x{i} ;
        y0 = y{i} ;
        uu = u_filtered{i} ;
        vv = v_filtered{i} ; 
        % Get position in next timepoint in pixels (in xy plane)
        % Clip the x position to the size of the image, and wrap the y position
        x1 = x0 + uu ;
        eps = 1e-14 ;
        x1 = max(x1, eps) ;
        x1 = min(x1, xesz-eps) ;
        y1 = mod(y0 + vv, Ysc1) ;

        % get embedded vector in R^3 for t0
        % Obtain the equivalent of v2D: ie the 2D vertices: mesh0.u
        mesh0x = x2Xpix(mesh0.sphi(:, 1), Ysc0, xsc0) ;
        if use_shifty
            % The shift in pixels of the current frame = shifty(i)
            mesh0y = y2Ypix(mesh0.sphi(:, 2), shifty(i), Ysc0) ;
        else
            mesh0y = y2Ypix(mesh0.sphi(:, 2), Ysc0) ;
        end
        
        % Create extended mesh (copied above and below), also shifted by shifty
        meshxy = [mesh0x, mesh0y ] ;
        mabove = [mesh0x, mesh0y + dy2dY(1., Ysc0)] ;
        mbelow = [mesh0x, mesh0y - dy2dY(1., Ysc0)] ;
        mabove2 = [mesh0x, mesh0y + dy2dY(2., Ysc0)] ;
        mbelow2 = [mesh0x, mesh0y - dy2dY(2., Ysc0)] ;
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
            hc = triplot(tr0, 'Color', blue) ;
            hold on
            ho = triplot(tr0_orig, 'Color', orange) ;
            axis equal
            title('Triangulation in pullback space')
            ylim([0 yesz])
            xlabel('x [pix]')
            ylabel('y [pix]')
            legend({'tiled', 'original'})
            saveas(fig, fullfile(pivOutDir, ['tri_' timestr '.png']))
            if preview
                set(fig, 'Visible', 'On')
                waitfor(fig)
            end
            close all
        end

        % Find xyz for matching position in t1 xy plane
        % get embedded vector in R^3 for t0
        % The shift in pixels of the current frame = shifty(i)
        mesh1x = x2Xpix(mesh1.sphi(:, 1), Ysc1, xsc1) ;
        if use_shifty
            mesh1y = y2Ypix(mesh1.sphi(:, 2), shifty(i + 1), Ysc1) ;
        else
            mesh1y = y2Ypix(mesh1.sphi(:, 2), Ysc1) ;
        end
        % Create extended mesh (copied above and below), also shifted by shifty
        meshxy = [mesh1x, mesh1y ] ;
        mabove = [mesh1x, mesh1y + dy2dY(1., Ysc1)] ;
        mbelow = [mesh1x, mesh1y - dy2dY(1., Ysc1)] ;
        mabove2 = [mesh1x, mesh1y + dy2dY(2., Ysc1)] ;
        mbelow2 = [mesh1x, mesh1y - dy2dY(2., Ysc1)] ;
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

