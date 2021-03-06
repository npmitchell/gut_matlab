function visualizeMeshEvolution(QS, options)
% Plot meshes as a morphsnakes demo
% Load all meshes in demo output directory, pass them through meshlab
% filters. Plot them in orthogonal sectioning view in 3D.
% Given plys from morphsnakes in ./growing_plys/, deposit smoothed meshes
% into ./growing_meshlab_plys/, then make figures in ./demofigs
% This method is based on figure_demo_morphsnakes_visualize.m script
%
% NPMitchell 2021

%% Unpack options
if nargin < 2
    options = struct() ;
end

% Default options 
plot_growth = false ;       % bool, plot the growth of the mesh to fill the organ for a single (first?) timepoint
plot_evolution = true ;     % bool, plot the evolution of the mesh over time
faceon = false ;            % how to visualize the results
preview = false ;           % preview intermediate results
plotXslice = false ;        % optionally create orthoview of Xslice
brighten = 100 ;            % factor by which to brighten scaled IV data
x0 = 0 ;                    % x value for orthogonal slice
y0 = 20 ;                   % y value for orthogonal slice
z0 = -5 ;                   % z value for orthogonal slice

% texturepatch options
meshFileBase = QS.fullFileBase.mesh ;
normal_shift = QS.normalShift ;
flipy = QS.flipy ;
texture_axis_order = QS.data.axisOrder ;
t0 = QS.t0set() ;

% figure parameters
xwidth = 16 ; % cm
ywidth = 10 ; % cm

if isfield(options, 'plot_growth')
    plot_growth = options.plot_growth ;
end
if isfield(options, 'plot_evolution')
    plot_evolution = options.plot_evolution ;
end
if isfield(options, 'faceon')
    faceon = options.faceon ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'plotXslice')
    plotXslice = options.plotXslice ;
end
if isfield(options, 'adjust_high')
    adjust_high = options.adjust_high ;
end
if isfield(options, 'x0')
    x0 = options.x0 ;
end
if isfield(options, 'y0')
    y0 = options.y0 ;
end
if isfield(options, 'z0')
    z0 = options.z0 ;
end
if isfield(options, 't0')
    t0 = options.t0 ;
end
if isfield(options, 'xwidth')
    xwidth = options.xwidth ;
end
if isfield(options, 'ywidth')
    ywidth = options.ywidth ;
end

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
if plot_growth
    error('this section needs to be edited to be more like the subsequent evolution sequence')
    %% Options
    ofn_ply = 'demo_' ;
    ofn_smoothply = 'meshlabout_' ;
    rawh5fn = './Time_000110_c1_stab.h5' ;
    probh5fn = './Time_000110_c1_stab_Probabilities.h5' ;
    use_pointcloud = false ; % Use the pointcloud from the level set rather
                             % than the boundary mesh from marching cubes.
                             % Requires the level sets as h5 files
    dtype = '.h5' ;
    ls_outfn = 'ms' ;  % the level set filename to be loaded (if necessary)
    mlxprogram = '/mnt/data/code/meshlab_codes/surface_rm_resample20k_reconstruct_LS3_1p2pc_ssfactor4.mlx' ;

    %% Directories
    demoDir = './' ;
    mslsDir = fullfile(demoDir, 'growing_plys') ;
    meshlabOutDir = fullfile(demoDir, 'growing_meshlab_plys') ;
    if faceon
        outdir = fullfile(demoDir, 'demofigs') ;
    else
        outdir = fullfile(demoDir, 'demofigs_view2') ;
    end
    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end

    %% Load the raw data from h5
    raw = h5read(rawh5fn, '/inputData') ;
    % For plotting
    xslice = round(size(raw, 1) * 0.5) ;
    yslice = round(size(raw, 2) * 0.5) ;
    zslice = round(size(raw, 3) * 0.5) ;

    %% find all ms_...ply files in mslsDir, and smooth them all
    files_to_smooth = dir(fullfile(mslsDir, [ofn_ply '*.ply'])) ;
    lsfns_to_smooth = dir(fullfile(mslsDir, [ls_outfn '*' dtype])) ;

    for i=1:length(files_to_smooth)
        msls_mesh_outfn = files_to_smooth(i).name ;
        PCfile = fullfile( mslsDir, msls_mesh_outfn );
        % Note that LS file is outputLs ;
        split_fn = strsplit(msls_mesh_outfn, ofn_ply) ;
        extension_outfn = split_fn{end} ;
        base_outfn_for_pointcloud = ofn_ply ;
        mesh_outfn = [ofn_smoothply, extension_outfn] ;
        outputMesh = fullfile(meshlabOutDir, mesh_outfn);

        disp(['outputMesh = ', outputMesh])
        %bad = so_bad
        if ~exist( outputMesh, 'file')
            if use_pointcloud
                % Use the pointcloud from the level set rather than the
                % boundary mesh from marching cubes
                %----------------------------------------------------------------------
                % Extract the implicit level set as a 3D binary array
                %----------------------------------------------------------------------

                % The file name of the current time point
                ls_outfn_ii = lsfns_to_smooth(i).name ;
                % The 3D binay array
                bwLS = h5read( ls_outfn_ii, '/implicit_levelset' );

                % Extract the (x,y,z)-locations of the level set boundary (in pixel
                % space)
                bwBdyIDx = bwperim( bwLS );

                clear bwBdy
                [ bwBdy(:,1), bwBdy(:,2), bwBdy(:,3) ] = ind2sub( size(bwLS), ...
                    find(bwBdyIDx) );

                %----------------------------------------------------------------------
                % Create output mesh
                %----------------------------------------------------------------------

                % Write the points to a .obj file as a point cloud for ouput to Meshlab
                clear OBJ
                OBJ.vertices = bwBdy;
                OBJ.objects(1).type='f';
                OBJ.objects(1).data.vertices=[];

                pointcloud_fn = [base_outfn_for_pointcloud '.obj'] ;
                disp(['Writing point cloud ' pointcloud_fn]);
                write_wobj(OBJ, pointcloud_fn );

                % Run the meshlab script
                system( ['meshlabserver -i ' pointCloudFileName, ...
                    ' -o ' outputMesh, ...
                    ' -s ' mlxprogram ' -om vn']);
            else
                % Use the marching cubes mesh surface to smooth
                command = ['meshlabserver -i ' PCfile ' -o ' outputMesh, ...
                    ' -s ' mlxprogram ' -om vn'];
                % Either copy the command to the clipboard
                clipboard('copy', command);
                % or else run it on the system
                disp(['running ' command])
                system(command)
            end
        else
            disp(['t=', num2str(i) ': smoothed mesh file found...'])
        end

    end

    %% Plot with orthogonal sectioning (note: faceon gives face-on view(2))
    % brighten the raw data
    rawclip = raw;
    thres = 5000 ;
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
    
    % X=const slice
    xraw = double(squeeze(rawclip(xslice, :, :))) ;
    xraw = 64 * xraw / max(xraw(:)) ;
    [xpy, xpz] = meshgrid(1:size(xraw, 1), 1:size(xraw, 2)) ;
    xpx = xslice * ones(size(xpy)) ;
    
    % Y=const slice
    if ~faceon && plotYslice
        yraw = double(squeeze(rawclip(:, yslice, :))) ;
        yraw = 64 * yraw / max(yraw(:)) ;
        [ypx, ypz] = meshgrid(1:size(yraw, 1), 1:size(yraw, 2)) ;
        ypy = yslice * ones(size(ypx)) ;
    end
    
    % Z=const slice
    if ~faceon
        zraw = double(squeeze(rawclip(:, :, zslice))) ;
        zraw = 64 * zraw / max(zraw(:)) ;
        [zpx, zpy] = meshgrid(1:size(zraw, 1), 1:size(zraw, 2)) ;
        zpz = zslice * ones(size(zpx)) ;
    end
    
    % Plot X slice
    surface(xpz, xpy, xpx, xraw', 'FaceColor','texturemap', ...
        'EdgeColor','none','CDataMapping','direct')
    hold on;

    % Plot Y slice
    if ~faceon && plotYslice
        surface(ypz, ypy, ypx, yraw', 'FaceColor','texturemap', ...
            'EdgeColor','none','CDataMapping','direct')
    end
    
    % Plot Z slice
    if ~faceon
        surface(zpz, zpy, zpx, zraw', 'FaceColor','texturemap', ...
           'EdgeColor','none','CDataMapping','direct')
    end
    % axis equal  % debug: toggle this off?
    colormap gray
    % hold on;

    % Plot each in turn
    fns = dir(fullfile(meshlabOutDir, [ofn_smoothply '*.ply'])) ;
    for ii = 1:length(fns)
        disp(['Reading ply: ' fns(ii).name])
        mesh = read_ply_mod(fullfile(fns(ii).folder, fns(ii).name)) ;
        name = strsplit(fns(ii).name, '.ply') ;
        name = name{1} ;
        xx = mesh.v(:, 1);
        yy = mesh.v(:, 2);
        zz = mesh.v(:, 3);
        h = trimesh(mesh.f, xx, yy, zz, 'EdgeColor', 'none', 'facecolor', lscolor) ;
        axis equal
        axis off

        lighting gouraud    % preferred method for lighting curved surfaces
        material dull    % set material to be dull, no specular highlights

        % Obtain xyzlimits
        if ii  == 1
            if faceon
                lgt = camlight ;
            else
                view([2,2,2])
                lgt = camlight('headlight') ;
            end
            xlims = xlim ;
            ylims = ylim ;
            zlims = zlim ;
        end
        xlim(xlims) ;
        ylim(ylims) ;
        if ~faceon
            zlim(zlims) ;
        end

        % xlabel('x')
        % ylabel('y')
        % zlabel('z')
        fn = fullfile(outdir, [name '.png']) ;
        disp(['Saving ' fn])
        % saveas(fig, fn) ;

        patchIm = getframe(gca);
        imwrite( patchIm.cdata, fn);

        delete(h)
    end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART II: EVOLVING SEQUENCE =============================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% navigate to the outputdir for figures
if plot_evolution
    if faceon
        outdir = fullfile(QS.dir.mesh, 'demo_morphsnakes_figs', 'faceon') ;
    else
        outdir = fullfile(QS.dir.mesh, 'demo_morphsnakes_figs', 'perspective') ;
    end
    if ~exist(fullfile(outdir, 'mesh'), 'dir')
        mkdir(fullfile(outdir, 'mesh'))
        mkdir(fullfile(outdir, 'texture'))
    end
    projectDir = QS.dir.data ;
    meshDir = QS.dir.mesh ;
    [rot, trans] = QS.getRotTrans ;
    pix2um = QS.APDV.resolution ;
        
    % rotfn = fullfile(meshDir, 'rotation_APDV.txt') ;
    % transfn = fullfile(meshDir, 'translation_APDV.txt') ;
    % rot = dlmread(rotfn) ;
    % trans = dlmread(transfn) ;
    % pix2um = 0.2619 ;
    
    plyfn = QS.fullFileBase.mesh ;
    % h5fn = 'Time_%06d_c1_stab.h5' ;
    % tiffn = 'Time_%06d_c1_stab.tif' ;
    % fns = dir(fullfile(meshDir, [plyfn '*.ply']));
    
    timePoints = QS.xp.fileMeta.timePoints ;
    
    first = true ;
    tidx2do = [190] ;
    tidx2do = [tidx2do, setdiff(1:length(timePoints), tidx2do)] ;
    for ii = tidx2do
    
        tp = timePoints(ii) ;
        QS.setTime(tp) ;
        timestr = sprintf('%06d', tp) ;
        disp(['Considering time ' timestr])
        
        fn = fullfile(outdir, 'mesh', [sprintf(QS.fileBase.name, tp) '.png']) ;
        fn_texture = fullfile(outdir, 'texture', [sprintf(QS.fileBase.name, tp) '.png']) ;
        
        if ~exist(fn, 'file')
            % Plot with h5
            ploth5 = false ;
            if ploth5
                rawh5fn = fullfile(projectDir, sprintf(h5fn, str2num(timestr))) ;
                % Load the raw data from h5
                raw = h5read(rawh5fn, '/inputData') ;
            else
                % tiffn_i = fullfile(projectDir, sprintf(tiffn, str2double(timestr))) ;
                % Load the raw data from tiff
                % raw = readTiff4D(tiffn_i, 1, 'xyz') ;
                raw = QS.getCurrentData() ;
                IV = raw{1} ;

            end

            % Raw and APDV aligned Meshes
            rmesh = QS.loadCurrentRawMesh() ;
            amesh = QS.loadCurrentAlignedMesh() ;

            % raw mesh vertices
            xr = rmesh.v(:, 1)  ;
            yr = rmesh.v(:, 2)  ;
            zr = rmesh.v(:, 3)  ;

            % scaled mesh vertices
            xa = amesh.v(:, 1)  ;
            if QS.flipy
                ya = - amesh.v(:, 2)  ;
            else
                ya = amesh.v(:, 2)  ;
            end
            za = amesh.v(:, 3)  ;

            % For plotting
            xslice = round(size(raw, 1) * 0.5) ;
            yslice = round(size(raw, 2) * 0.5) ;
            zslice = round(size(raw, 3) * 0.5) ;

            %% Plot with orthogonal sectioning (note: faceon gives face-on view(2))
            % brighten the raw data
            % rawclip = raw;
            % thres = 5000 ; %Lowest intensity to consider a signal
            % clip = 2 * median(raw(raw > thres)) ;
            % rawclip(raw > clip) = clip ;

            close all
            if faceon
                fig = figure('units', 'centimeters', ...
                    'outerposition', [0 0 10 25], 'visible', 'off') ;
            else
                fig = figure('units', 'centimeters', ...
                    'outerposition', [0 0 30 30], 'visible', 'off') ;
            end

            %% Interpolate the data and evaluate at rotated orthogonal planes     
            active_rotation = false ; % too slow for large data
            if active_rotation
                %% ACTIVE ROTATION: TRANSLATE & ROTATE intensity XYZ points
                disp('trans & rotate...')
                xi = 1:size(rawclip, 1) ;
                yi = 1:size(rawclip, 2) ;
                zi = 1:size(rawclip, 3) ;
                [xr, yr, zr] = meshgrid(xi, yi, zi) ;
                xyzr = (rot * [xr(:), yr(:), zr(:)]')' + trans ;
                xxx = reshape(xyzr(:, 1), size(rawclip)) ;
                yyy = reshape(xyzr(:, 2), size(rawclip)) ;
                zzz = reshape(xyzr(:, 3), size(rawclip)) ;

                %% Interpolate the data and evaluate at rotated orthogonal planes
                disp('defining interpolation...')
                xspace = linspace(min(xxx(:)), max(xxx(:)), 100) ;
                yspace = linspace(min(yyy(:)), max(yyy(:)), 100) ;
                zspace = linspace(min(zzz(:)), max(zzz(:)), 100) ;
                [yq, zq ] = meshgrid(xspace, yspace, zspace) ;
                interpI = scatteredInterpolant(xxx(:), yyy(:), zzz(:), double(rawclip(:))) ;
                newIx = interpI(0*yq, yq, zq) ;
            else
                %% PASSIVE ROTATION            
                % Obtain rotation matrix that undoes the APDV rotation 
                invRot = QS.invertRotation(rot) ;

                % Plane to image
                xMin = -30 ;
                xMax = 300 ;
                yMin = -80 ;
                yMax = 80 ;
                zMin = -80 ;
                zMax = 80 ;
                oversample_factor = 1.1 ;
                nX = round((xMax - xMin) / pix2um * oversample_factor) ;
                nY = round((yMax - yMin) / pix2um * oversample_factor) ; 
                nZ = round((zMax - zMin) / pix2um * oversample_factor) ;
                x1 = x0 * ones(nZ, nY) ;
                y2 = y0 * ones(nZ, nX) ;
                z3 = z0 * ones(nY, nX) ;
                xspace = linspace(xMin, xMax, nX) ;
                yspace = linspace(yMin, yMax, nY) ;
                zspace = linspace(zMin, zMax, nZ) ;
                [y1, z1] = meshgrid(yspace, zspace) ; 
                [x2, z2] = meshgrid(xspace, zspace) ; 
                [x3, y3] = meshgrid(xspace, yspace) ; 
                xyzX = (invRot * (([x1(:), y1(:), z1(:)]/pix2um) - trans)')' ;
                xyzY = (invRot * (([x2(:), y2(:), z2(:)]/pix2um) - trans)')' ;
                xyzZ = (invRot * (([x3(:), y3(:), z3(:)]/pix2um) - trans)')' ;
                xpx = reshape(xyzX(:, 1), [nY, nZ]) ;
                xpy = reshape(xyzX(:, 2), [nY, nZ]) ; 
                xpz = reshape(xyzX(:, 3), [nY, nZ]) ;
                ypx = reshape(xyzY(:, 1), [nX, nZ]) ;
                ypy = reshape(xyzY(:, 2), [nX, nZ]) ; 
                ypz = reshape(xyzY(:, 3), [nX, nZ]) ;
                zpx = reshape(xyzZ(:, 1), [nX, nY]) ;
                zpy = reshape(xyzZ(:, 2), [nX, nY]) ; 
                zpz = reshape(xyzZ(:, 3), [nX, nY]) ;
                assert(all(size(x1) == size(y1)))
                assert(all(size(y1) == size(z1)))
                assert(all(size(x2) == size(y2)))
                assert(all(size(y2) == size(z2)))
                assert(all(size(x3) == size(y3)))
                assert(all(size(y3) == size(z3)))

                % preview with mesh
                if preview
                    clf
                    h = trimesh(rmesh.f, xr, yr, zr, 'EdgeColor', 'none', 'facecolor', lscolor) ;
                    axis equal
                    hold on;
                    scatter3(xpx(:), xpy(:), xpz(:), 10, 'markeredgecolor', 'r')
                    scatter3(ypx(:), ypy(:), ypz(:), 10, 'markeredgecolor', 'g')
                    scatter3(zpx(:), zpy(:), zpz(:), 10, 'markeredgecolor', 'b')
                    waitfor(gcf)
                end
                disp('performing interpolation in Xplane')
                ix = interp3(double(IV), xyzX(:, 2), xyzX(:, 1), xyzX(:, 3)) ;
                ix = reshape(ix, size(x1)) ;

                disp('performing interpolation in Yplane')
                iy = interp3(double(IV), xyzY(:, 2), xyzY(:, 1), xyzY(:, 3)) ;
                iy = reshape(iy, size(x2)) ;

                disp('performing interpolation in Zplane')
                iz = interp3(double(IV), xyzZ(:, 2), xyzZ(:, 1), xyzZ(:, 3)) ;
                iz = reshape(iz, size(x3)) ;

                if preview
                    % Check order of axes wrt IV
                    figure ;
                    subplot(2, 2, 1)
                    imagesc(1:size(IV, 2), 1:size(IV, 3), squeeze(IV(round(mean(xyzX(:, 1))), :, :))')
                    xlabel('size(IV, 2)')
                    ylabel('size(IV, 3)'); axis equal; title('IX (approx, pre-rot)')
                    subplot(2, 2, 2)
                    imagesc(1:size(IV, 1), 1:size(IV, 3), squeeze(IV(:, round(mean(xyzY(:, 2))), :))')
                    xlabel('size(IV, 1)')
                    ylabel('size(IV, 3)'); axis equal; title('IY (approx, pre-rot)')
                    subplot(2, 2, 3)
                    imagesc(1:size(IV, 1), 1:size(IV, 2), squeeze(IV(:, :, round(mean(xyzZ(:, 3)))))')
                    xlabel('size(IV, 1)')
                    ylabel('size(IV, 2)'); axis equal; title('IZ (approx, pre-rot)')           


                    % check interpolation
                    figure; 
                    subplot(2, 2, 1); imagesc(yspace, zspace, ix'); axis equal ;
                    subplot(2, 2, 2); imagesc(xspace, zspace, iy'); axis equal ;
                    subplot(2, 2, 3); imagesc(xspace, yspace, iz'); axis equal ;
                    waitfor(gcf)
                    close all
                end

            end

            ix(isnan(ix)) = 0 ;
            iy(isnan(iy)) = 0 ;
            iz(isnan(iz)) = 0 ;
            ix = brighten * ix / 2^(16) ;
            iy = brighten * iy / 2^(16) ;
            iz = brighten * iz / 2^(16) ;

            %% Plot in APDV
            clf
            
            hold on;
            surface(x2, y2, z2, iy, 'FaceColor','texturemap', ...
                'EdgeColor','none','CDataMapping','direct')
            hold on;
            
            % Plot X slice (optional)
            if ~faceon && plotXslice                
                surface(x1, y1, z1, ix, 'FaceColor','texturemap', ...
                    'EdgeColor','none','CDataMapping','direct')
            end
            
            % Plot Z slice (if not looking laterally)
            if ~faceon
                surface(x3, y3, z3,  iz, 'FaceColor','texturemap', ...
                   'EdgeColor','none','CDataMapping','direct')
                hold on;
            end
            colormap gray
            h = trimesh(amesh.f, xa, ya, za, 'EdgeColor', 'none', ...
                'facecolor', lscolor) ;
            axis equal
            hold on;
            axis off

            lighting gouraud    % preferred method for lighting curved surfaces
            material dull    % set material to be dull, no specular highlights

            % Lighting and view 
            if faceon
                lgt = camlight ;
            else
                view(-20, 20)
                % view([-0.75, -1.0, 0.7]) 
                lgt = camlight('headlight') ;
            end

            % Obtain xyzlimits
            if first
                xlims = xlim ;
                ylims = ylim ;
                zlims = zlim ;
                first = false ;
            end
            xlim(xlims) ;
            ylim(ylims) ;
            if ~faceon
                zlim(zlims) ;
            end
            disp(['Saving ' fn])

            set(fig, 'PaperUnits', 'centimeters');
            set(fig, 'PaperPosition', [0 0 xwidth ywidth]);

            % saveas(fig, fn) ;
            patchIm = getframe(gca);
            imwrite( patchIm.cdata, fn);
            close all
        end
        
        %% Create matching Texturepatch Image
        if ~exist(fn_texture, 'file')
            % Psize is the linear dimension of the grid drawn on each triangular face
            Options.PSize = 5;
            Options.EdgeColor = 'none';
            Options.Rotation = rot ;
            Options.Translation = trans ;
            Options.Dilation = resolution ;
            Options.numLayers = [1, 1];  % [2, 2] marches ~0.5 um in either dir
            Options.layerSpacing = 2 ;

            % Raw and APDV aligned Meshes
            rmesh = QS.loadCurrentRawMesh() ;
            amesh = QS.loadCurrentAlignedMesh() ;
            
            fig = figure('Visible', 'Off') ;
            disp(['creating texture patch ' num2str(tp, '%06d')])

            % Allow for axis swapping
            TV = rmesh.v(:, texture_axis_order) ;
            texture_patch_3d( amesh.f, VV, rmesh.f, TV, IV, Options );

            % format the figure
            disp('formatting figure...')
            axis equal
            xlim([xmin, xmax])
            ylim([ymin, ymax])
            zlim([zmin, zmax])
            colormap bone
            titlestr = ['$t = $' num2str(tp*QS.timeInterval-t0) ' ' QS.timeUnits] ;
            title(titlestr, 'Interpreter', 'Latex', 'Color', 'white') 
            xlabel('AP position [$\mu$m]', 'Interpreter', 'Latex')
            ylabel('lateral position [$\mu$m]', 'Interpreter', 'Latex')
            zlabel('DV position [$\mu$m]', 'Interpreter', 'Latex')
            patchIm = getframe(gca);
            imwrite( patchIm.cdata, fn_texture);
        end
    end
    
end

