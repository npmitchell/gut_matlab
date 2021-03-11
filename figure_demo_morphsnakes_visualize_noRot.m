%% DEPRICATED, SEE QS.VISUALIZEMESHEVOLUTION()

%% Smoothen and plot meshes from morphsnakes demo
% Load all meshes in demo output directory, pass them through meshlab
% filters. Plot them in orthogonal sectioning view in 3D.
% Given plys from morphsnakes in ./growing_plys/, deposit smoothed meshes
% into ./growing_meshlab_plys/, then make figures in ./demofigs

% Add path for read_ply_mod
addpath('/mnt/data/code/imsaneV1.2.3/generalfunctions/mesh/');
% Add path for colors
addpath('/mnt/data/code/gut_matlab/plotting/') ;

% Options for which to do
do_part1 = false ;
do_part2 = true ;
faceon = false ;  % how to visualize the results

%% Colors
colors = define_colors(7) ; 
blue   = colors(1, :) ;
red    = colors(2, :) ;
yellow = colors(3, :) ;
purple = colors(4, :) ;
green  = colors(5, :) ;
sky    = colors(6, :) ;
maroon = colors(6, :) ;
lscolor = sky ;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% PART I: GROWING SEQUENCE ===============================================
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if do_part1
%     %% Options
%     ofn_ply = 'demo_' ;
%     ofn_smoothply = 'meshlabout_' ;
%     rawh5fn = './Time_000110_c1_stab.h5' ;
%     probh5fn = './Time_000110_c1_stab_Probabilities.h5' ;
%     use_pointcloud = false ; % Use the pointcloud from the level set rather
%                              % than the boundary mesh from marching cubes.
%                              % Requires the level sets as h5 files
%     dtype = '.h5' ;
%     ls_outfn = 'ms' ;  % the level set filename to be loaded (if necessary)
%     mlxprogram = '/mnt/data/code/meshlab_codes/surface_rm_resample20k_reconstruct_LS3_1p2pc_ssfactor4.mlx' ;
% 
%     %% Directories
%     demoDir = './' ;
%     mslsDir = fullfile(demoDir, 'growing_plys') ;
%     meshlabOutDir = fullfile(demoDir, 'growing_meshlab_plys') ;
%     if faceon
%         outdir = fullfile(demoDir, 'demofigs') ;
%     else
%         outdir = fullfile(demoDir, 'demofigs_view2') ;
%     end
%     if ~exist(outdir, 'dir')
%         mkdir(outdir)
%     end
% 
%     %% Load the raw data from h5
%     raw = h5read(rawh5fn, '/inputData') ;
%     % For plotting
%     xslice = round(size(raw, 1) * 0.5) ;
%     yslice = round(size(raw, 2) * 0.5) ;
%     zslice = round(size(raw, 3) * 0.5) ;
% 
%     %% Load probabilities
%     % probs = h5read(probh5fn, '/exported_data') ;
% 
%     %% find all ms_...ply files in mslsDir, and smooth them all
%     files_to_smooth = dir(fullfile(mslsDir, [ofn_ply '*.ply'])) ;
%     lsfns_to_smooth = dir(fullfile(mslsDir, [ls_outfn '*' dtype])) ;
% 
%     for i=1:length(files_to_smooth)
%         msls_mesh_outfn = files_to_smooth(i).name ;
%         PCfile = fullfile( mslsDir, msls_mesh_outfn );
%         % Note that LS file is outputLs ;
%         split_fn = strsplit(msls_mesh_outfn, ofn_ply) ;
%         extension_outfn = split_fn{end} ;
%         base_outfn_for_pointcloud = ofn_ply ;
%         mesh_outfn = [ofn_smoothply, extension_outfn] ;
%         outputMesh = fullfile(meshlabOutDir, mesh_outfn);
% 
%         disp(['outputMesh = ', outputMesh])
%         %bad = so_bad
%         if ~exist( outputMesh, 'file')
%             if use_pointcloud
%                 % Use the pointcloud from the level set rather than the
%                 % boundary mesh from marching cubes
%                 %----------------------------------------------------------------------
%                 % Extract the implicit level set as a 3D binary array
%                 %----------------------------------------------------------------------
% 
%                 % The file name of the current time point
%                 ls_outfn_ii = lsfns_to_smooth(i).name ;
%                 % The 3D binay array
%                 bwLS = h5read( ls_outfn_ii, '/implicit_levelset' );
% 
%                 % Extract the (x,y,z)-locations of the level set boundary (in pixel
%                 % space)
%                 bwBdyIDx = bwperim( bwLS );
% 
%                 clear bwBdy
%                 [ bwBdy(:,1), bwBdy(:,2), bwBdy(:,3) ] = ind2sub( size(bwLS), ...
%                     find(bwBdyIDx) );
% 
%                 %----------------------------------------------------------------------
%                 % Create output mesh
%                 %----------------------------------------------------------------------
% 
%                 % Write the points to a .obj file as a point cloud for ouput to Meshlab
%                 clear OBJ
%                 OBJ.vertices = bwBdy;
%                 OBJ.objects(1).type='f';
%                 OBJ.objects(1).data.vertices=[];
% 
%                 pointcloud_fn = [base_outfn_for_pointcloud '.obj'] ;
%                 disp(['Writing point cloud ' pointcloud_fn]);
%                 write_wobj(OBJ, pointcloud_fn );
% 
%                 % Run the meshlab script
%                 system( ['meshlabserver -i ' pointCloudFileName, ...
%                     ' -o ' outputMesh, ...
%                     ' -s ' mlxprogram ' -om vn']);
%             else
%                 % Use the marching cubes mesh surface to smooth
%                 command = ['meshlabserver -i ' PCfile ' -o ' outputMesh, ...
%                     ' -s ' mlxprogram ' -om vn'];
%                 % Either copy the command to the clipboard
%                 clipboard('copy', command);
%                 % or else run it on the system
%                 disp(['running ' command])
%                 system(command)
%             end
%         else
%             disp(['t=', num2str(i) ': smoothed mesh file found...'])
%         end
% 
%     end
% 
%     %% Plot with orthogonal sectioning (note: faceon gives face-on view(2))
%     % brighten the raw data
%     rawclip = raw;
%     thres = 5000 ;
%     clip = 2 * median(raw(raw > thres)) ;
%     rawclip(raw > clip) = clip ;
% 
%     close all
%     if faceon
%         fig = figure('units', 'centimeters', ...
%             'outerposition', [0 0 10 25], 'visible', 'off') ;
%     else
%         fig = figure('units', 'centimeters', ...
%             'outerposition', [0 0 30 30], 'visible', 'off') ;
%     end
%     % X=const slice
%     xraw = double(squeeze(rawclip(xslice, :, :))) ;
%     xraw = 64 * xraw / max(xraw(:)) ;
%     [xpy, xpz] = meshgrid(1:size(xraw, 1), 1:size(xraw, 2)) ;
%     xpx = xslice * ones(size(xpy)) ;
%     % Y=const slice
%     % yraw = double(squeeze(rawclip(:, yslice, :))) ;
%     % yraw = 64 * yraw / max(yraw(:)) ;
%     % [ypx, ypz] = meshgrid(1:size(yraw, 1), 1:size(yraw, 2)) ;
%     % ypy = yslice * ones(size(ypx)) ;
%     % Z=const slice
%     zraw = double(squeeze(rawclip(:, :, zslice))) ;
%     zraw = 64 * zraw / max(zraw(:)) ;
%     [zpx, zpy] = meshgrid(1:size(zraw, 1), 1:size(zraw, 2)) ;
%     zpz = zslice * ones(size(zpx)) ;
%     % Plot them all
%     surface(xpz, xpy, xpx, xraw', 'FaceColor','texturemap', ...
%         'EdgeColor','none','CDataMapping','direct')
%     hold on;
%     % surface(ypz, ypy, ypx, yraw', 'FaceColor','texturemap', ...
%     %     'EdgeColor','none','CDataMapping','direct')
% 
%     % Zslice works well here
%     if ~faceon
%         surface(zpz, zpy, zpx, zraw', 'FaceColor','texturemap', ...
%            'EdgeColor','none','CDataMapping','direct')
%     end
%     % axis equal  % debug: toggle this off?
%     colormap gray
%     % hold on;
% 
%     % Plot each in turn
%     fns = dir(fullfile(meshlabOutDir, [ofn_smoothply '*.ply'])) ;
%     for ii = 1:length(fns)
%         disp(['Reading ply: ' fns(ii).name])
%         mesh = read_ply_mod(fullfile(fns(ii).folder, fns(ii).name)) ;
%         name = strsplit(fns(ii).name, '.ply') ;
%         name = name{1} ;
%         xx = mesh.v(:, 1);
%         yy = mesh.v(:, 2);
%         zz = mesh.v(:, 3);
%         h = trimesh(mesh.f, xx, yy, zz, 'EdgeColor', 'none', 'facecolor', lscolor) ;
%         axis equal
%         axis off
% 
%         lighting gouraud    % preferred method for lighting curved surfaces
%         material dull    % set material to be dull, no specular highlights
% 
%         % Obtain xyzlimits
%         if ii  == 1
%             if faceon
%                 lgt = camlight ;
%             else
%                 view([2,2,2])
%                 lgt = camlight('headlight') ;
%             end
%             xlims = xlim ;
%             ylims = ylim ;
%             zlims = zlim ;
%         end
%         xlim(xlims) ;
%         ylim(ylims) ;
%         if ~faceon
%             zlim(zlims) ;
%         end
% 
%         % xlabel('x')
%         % ylabel('y')
%         % zlabel('z')
%         fn = fullfile(outdir, [name '.png']) ;
%         disp(['Saving ' fn])
%         % saveas(fig, fn) ;
% 
%         patchIm = getframe(gca);
%         imwrite( patchIm.cdata, fn);
% 
%         delete(h)
%     end   
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART II: EVOLVING SEQUENCE =============================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% navigate to the outputdir for figures
ssfactor = 4; 
if do_part2
    if faceon
        outdir = '/mnt/data/analysis/demofigs_view2/' ;
    else
        outdir = '/mnt/data/analysis/demofigs/' ;
    end
    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end
    projectDir = '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/' ;
    meshDir = fullfile(projectDir, 'msls_output/') ;
    rotfn = fullfile(meshDir, 'rotation_APDV.txt') ;
    rot = dlmread(rotfn) ;
    rot = inv(rot) ;
    
    plyfn = 'mesh_apical_stab_' ;
    h5fn = 'Time_%06d_c1_stab.h5' ;
    tiffn = 'Time_%06d_c1_stab.tif' ;
    fns = dir(fullfile(meshDir, [plyfn '*.ply']));
    
    tidx2do = [190] ;
    tidx2do = [tidx2do, setdiff(1:length(fns), tidx2do)] ;
    for ii = tidx2do
        meshfn = fns(ii).name ;
        timestr = strsplit(meshfn, plyfn) ;
        timestr = strsplit(timestr{2}, '.ply') ;
        timestr = timestr{1} ;
        disp(['Considering time ' timestr])
        
        % Plot with h5
        ploth5 = false ;
        if ploth5
            rawh5fn = fullfile(projectDir, sprintf(h5fn, str2num(timestr))) ;
            % Load the raw data from h5
            raw = h5read(rawh5fn, '/inputData') ;
        else
            tiffn_i = fullfile(projectDir, sprintf(tiffn, str2double(timestr))) ;
            % Load the raw data from tiff
            raw = readTiff4D(tiffn_i, 1, 'xyz') ;
        end
        
        % For plotting
        xslice = round(size(raw, 1) * 0.5) ;
        yslice = round(size(raw, 2) * 0.5) ;
        zslice = round(size(raw, 3) * 0.5) ;
        
        %% Plot with orthogonal sectioning (note: faceon gives face-on view(2))
        % brighten the raw data
        rawclip = raw;
        thres = 5000 ; %Lowest intensity to consider a signal
        clip = 2 * median(raw(raw > thres)) ;
        rawclip(raw > clip) = clip ;

        close all
        if faceon
            fig = figure('units', 'centimeters', ...
                'outerposition', [0 0 10 25], 'visible', 'off') ;
        else
            fig = figure('units', 'centimeters', ...
                'outerposition', [0 0 30 30], 'visible', 'off') ;
        end
        
        %% X=const slice
        xraw = double(squeeze(rawclip(xslice, :, :))) ;
        xraw = 64 * xraw / max(xraw(:)) ;
        [xpy, xpz] = meshgrid(1:size(xraw, 1), 1:size(xraw, 2)) ;
        xpx = xslice * ones(size(xpy)) ;
        
        %% Y=const slice
        yraw = double(squeeze(rawclip(:, yslice, :))) ;
        yraw = 64 * yraw / max(yraw(:)) ;
        [ypx, ypz] = meshgrid(1:size(yraw, 1), 1:size(yraw, 2)) ;
        ypy = yslice * ones(size(ypx)) ;
        
        %% Z=const slice
        zraw = double(squeeze(rawclip(:, :, zslice))) ;
        zraw = 64 * zraw / max(zraw(:)) ;
        [zpx, zpy] = meshgrid(1:size(zraw, 1), 1:size(zraw, 2)) ;
        zpz = zslice * ones(size(zpx)) ;
        
        %% Plot them all
        % surface(xpz, xpy, xpx, xraw', 'FaceColor','texturemap', ...
        %     'EdgeColor','none','CDataMapping','direct')
        hold on;
        surface(ypz, ypy, ypx, yraw', 'FaceColor','texturemap', ...
            'EdgeColor','none','CDataMapping','direct')
        hold on;
        
        % Zslice works well here
        if ~faceon
            surface(zpz, zpy, zpx, zraw', 'FaceColor','texturemap', ...
               'EdgeColor','none','CDataMapping','direct')
            hold on;
        end
        % axis equal  % debug: toggle this off?
        colormap gray
        % hold on;

        %% Add the mesh to the figure
        mesh = read_ply_mod(fullfile(fns(ii).folder, fns(ii).name)) ;
        name = strsplit(fns(ii).name, '.ply') ;
        name = name{1} ;
        % when read in as xyz, compared to mesh, we must take
        % x -> z
        % y -> y 
        % z -> -x
        xx = mesh.v(:, 3)  ;
        yy = mesh.v(:, 2)  ;
        zz = mesh.v(:, 1)  ;
        h = trimesh(mesh.f, xx, yy, zz, 'EdgeColor', 'none', 'facecolor', lscolor) ;
        axis equal
        set(gcf, 'visible', 'on')
        xlabel('x') ; ylabel('y'); zlabel('z')
        axis off

        lighting gouraud    % preferred method for lighting curved surfaces
        material dull    % set material to be dull, no specular highlights
        
        % Lighting and view 
        if faceon
            lgt = camlight ;
        else
            view([2,2,2])
            lgt = camlight('headlight') ;
        end
        
        % Obtain xyzlimits
        if ii  == 1
            xlims = xlim ;
            ylims = ylim ;
            zlims = zlim ;
        end
        xlim(xlims) ;
        ylim(ylims) ;
        if ~faceon
            zlim(zlims) ;
        end
        fn = fullfile(outdir, [name '.png']) ;
        disp(['Saving ' fn])
        % saveas(fig, fn) ;
        patchIm = getframe(gca);
        imwrite( patchIm.cdata, fn);

    end
    
end

