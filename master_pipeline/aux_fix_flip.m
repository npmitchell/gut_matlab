
% 
Options.res = 3.0 ;
xwidth = 16 ;                   % width of figure in cm
ywidth = 10 ;                   % height of figure in cm
% Figure options
colors = define_colors ;
blue = colors(1, :) ;
red = colors(2, :) ;
green = colors(5, :) ;

% Unpack parameters into strings
name = [ofn_smoothply '%06d' ] ;
expstr = strrep(num2str(Options.exponent, '%0.1f'), '.', 'p') ;
resstr = strrep(num2str(Options.res, '%0.1f'), '.', 'p') ;
extenstr = ['_exp' expstr '_res' resstr] ;
outdir = fullfile(meshDir, 'centerline') ;
figoutdir = fullfile(outdir, 'images') ;
fig1outdir = fullfile(figoutdir, 'centerline_xy') ;
fig2outdir = fullfile(figoutdir, 'centerline_xz') ;
fig3outdir = fullfile(figoutdir, 'centerline_yz') ;
fig1outname = [fullfile(fig1outdir, name) '_centerline' extenstr '_xy.png'] ;
fig2outname = [fullfile(fig2outdir, name) '_centerline' extenstr '_xz.png'] ;
fig3outname = [fullfile(fig3outdir, name) '_centerline' extenstr '_yz.png'] ;

% Rotate and scale (and mirror if flipy==true)
for tt = timePoints(115:end)
        
    % Load startpoint and endpoint
    tifname = sprintf(fn, tt) ;
    disp(['loading /' tifname '/spt,ept,dpt from ' startendptH5FileName])
    startpt = h5read(startendptH5FileName, ['/' tifname '/spt' ]) ;
    endpt = h5read(startendptH5FileName, ['/' tifname '/ept' ]) ;
    sptrs = h5read(startendptH5FileName, ['/' tifname '/sptrs' ]) ;
    eptrs = h5read(startendptH5FileName, ['/' tifname '/eptrs' ]) ;
    dptrs = h5read(startendptH5FileName, ['/' tifname '/dptrs' ]) ;

    % input filename mesh ending in 'ply'
    name = sprintf(meshFileBase, tt) ;
    if strcmp(name(end-3:end), '.ply')
        name = name(1:end-4) ;
    end

    skel_rs_outfn = [fullfile(outdir, name) '_centerline_scaled' extenstr ] ;
    skeloutfn = [skel_rs_outfn '.txt'] ;
    try
        ssskelrs = importdata(skeloutfn) ;
        redo = true ;
    catch
        disp('WARNING: COULD NOT LOAD FILE') 
        redo = false ;
    end
    
    if redo
        assert(tt > 113)
        assert(tt < 170)
        skelrs = ssskelrs(:, 2:4) ;
        sss = ssskelrs(:, 1) ;

        if flipy
            skelrs(:, 2) = -skelrs(:, 2) ;
        end

        %% Plot and save
        % Save plot of rotated and translated mesh
        disp('Saving rotated & translated figure (xy)...')    
        close all
        % Load rotated mesh
        if useSavedAPDVMeshes
            if strcmp(Options.meshAPDVFileName(end-3:end), '.ply')
                meshAPDVFileName = Options.meshAPDVFileName ;
            else
                meshAPDVFileName = [Options.meshAPDVFileName '.ply'] ;
            end
            disp(['Loading' sprintf(meshAPDVFileName, tt)])
            mesh = read_ply_mod(sprintf(meshAPDVFileName, tt)) ;
            xyzrs = mesh.v ;
            % If we are plotting the mesh, reverse the triangle ordering
            % for ambient occlusion to work properly, but the APDV mesh
            % should ALREADY have the correct y values (mirrored XZ)
            if flipy
                tri = mesh.f(:, [1 3 2]) ;
            end
        else
            xyzrs = ((rot * mesh.v')' + trans) * resolution ;
            if flipy
                xyzrs(:, 2) = -xyzrs(:, 2) ;
                tri = mesh.f(:, [1 3 2]) ;
            end
        end

        % Plot the result
        fig = figure('Visible', 'Off') ;
        tmp = trisurf(tri, xyzrs(:, 1), xyzrs(:,2), xyzrs(:, 3), ...
            'edgecolor', 'none', 'FaceAlpha', 0.1) ;
        % [~,~,~] = apply_ambient_occlusion(tmp, 'SoftLighting', true) ; 
        hold on;
        % plot the skeleton
        for i=1:length(skelrs)
            plot3(skelrs(:,1), skelrs(:,2), skelrs(:,3), ...
                '-','Color',[0,0,0], 'LineWidth', 3);
        end
        % annotate figure with APDV
        plot3(sptrs(1), sptrs(2), sptrs(3), 's', 'color', red)
        plot3(eptrs(1), eptrs(2), eptrs(3), '^', 'color', blue)
        plot3(dptrs(1), dptrs(2), dptrs(3), 'o', 'color', green)
        xlabel('x [$\mu$m]', 'Interpreter', 'Latex');
        ylabel('y [$\mu$m]', 'Interpreter', 'Latex');
        zlabel('z [$\mu$m]', 'Interpreter', 'Latex');
        title(['Centerline using $D^{' num2str(exponent) '}$: ' num2str(tt)], ...
            'Interpreter', 'Latex')
        axis equal
        % xy
        view(2)
        xlim(xyzlim_um(1, :)); 
        ylim(xyzlim_um(2, :)); 
        zlim(xyzlim_um(3, :)) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);
        saveas(fig, sprintf(fig1outname, tt))
        % yz
        disp('Saving rotated & translated figure (yz)...')    
        view(90, 0);
        xlim(xyzlim_um(1, :)); 
        ylim(xyzlim_um(2, :)); 
        zlim(xyzlim_um(3, :)) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=15cm
        saveas(fig, sprintf(fig2outname, tt))
        % xz
        disp('Saving rotated & translated figure (xz)...')  
        view(0, 0)    
        xlim(xyzlim_um(1, :)); 
        ylim(xyzlim_um(2, :)); 
        zlim(xyzlim_um(3, :)) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=15cm
        saveas(fig, sprintf(fig3outname, tt))
        close all


        %% Save the rotated, translated, scaled curve
        disp(['Saving rotated & scaled skeleton to txt: ', skel_rs_outfn, '.txt'])
        skeloutfn = [skel_rs_outfn '.txt'] ;
        fid = fopen(skeloutfn, 'wt');
        % make header
        fprintf(fid, 'Aligned & scaled skeleton: sss [um], skelrs [um]');  
        fclose(fid);
        dlmwrite(skeloutfn, [sss, skelrs])
    end
end
