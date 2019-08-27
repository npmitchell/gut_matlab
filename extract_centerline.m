%% Extract the centerlines from a series of meshes (PLY files)
% Noah Mitchell 2019
% This version relies on Gabriel Peyre's toolbox called
% toolbox_fast_marching/
% Run from the msls_output directory
clear ;

%% First, compile required c code
% mex ./FastMarching_version3b/shortestpath/rk4
close all ;
odir = pwd ;
codepath = '/mnt/data/code/gut_matlab/' ;
if ~exist(codepath, 'dir')
    codepath = [pwd filesep] ;
end
addpath(codepath)
addpath([codepath 'addpath_recurse' filesep]) ;
addpath([codepath 'mesh_handling' filesep]);
addpath([codepath 'inpolyhedron' filesep]);
addpath([codepath 'savgol' filesep])
addpath_recurse('/mnt/data/code/gptoolbox/')
% addpath_recurse([codepath 'gptoolbox' filesep])

toolbox_path = [codepath 'toolbox_fast_marching/toolbox_fast_marching/'];
dtpath = [codepath 'distanceTransform/'] ;
addpath_recurse(toolbox_path)
addpath(dtpath)
% compile_c_files
cd(odir)

%% Parameters
save_figs = false ;
res = 1 ;
resolution = 0.2619 ;
buffer = 5 ;
plot_buffer = 30;
ssfactor = 4; 
weight = 0.1;
normal_step = 0.5 ; 
preview = false ;
eps = 0.01 ;
meshorder = 'zyx' ;
exponent = 1;
% figure parameters
xwidth = 16 ; % cm
ywidth = 10 ; % cm

% Find all meshes to consider
meshdir = pwd ;
cd ../
rootdir = pwd ;
cd(meshdir)
% rootpath = '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/' ;
% rootpath = [rootpath 'Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/'] ;
% if ~exist(rootpath, 'dir')
%     rootpath = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/' ;
%     rootpath = [rootpath 'data/48Ygal4UasCAAXmCherry/201902072000_excellent/'] ;
% end
% meshdir = [rootpath 'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/'];
fns = dir(fullfile(meshdir, 'mesh_apical_stab_0*.ply')) ;
rotname = fullfile(meshdir, 'rotation_APDV') ;
transname = fullfile(meshdir, 'translation_APDV') ;
xyzlimname = fullfile(meshdir, 'xyzlim_APDV') ;

% Name output directory
outdir = [fullfile(fns(1).folder, 'centerline') filesep ];
if ~exist(outdir, 'dir')
    mkdir(outdir) ;
end
figoutdir = [outdir 'images' filesep];
if ~exist(figoutdir, 'dir')
    mkdir(figoutdir) ;
end
% figure 1
fig1outdir = [figoutdir 'centerline_xy' filesep];
if ~exist(fig1outdir, 'dir')
    mkdir(fig1outdir) ;
end
% figure 2
fig2outdir = [figoutdir 'centerline_xz' filesep];
if ~exist(fig2outdir, 'dir')
    mkdir(fig2outdir) ;
end
% figure 3
fig3outdir = [figoutdir 'centerline_yz' filesep];
if ~exist(fig3outdir, 'dir')
    mkdir(fig3outdir) ;
end
radius_vs_s_phi_outdir = [figoutdir 'radius_vs_s_phid' filesep];
if ~exist(radius_vs_s_phi_outdir, 'dir')
    mkdir(radius_vs_s_phi_outdir) ;
end
radius_vs_s_phicd_outdir = [figoutdir 'radius_vs_s_phicd' filesep];
if ~exist(radius_vs_s_phicd_outdir, 'dir')
    mkdir(radius_vs_s_phicd_outdir) ;
end

ii = 1 ;

%% Iterate through each mesh
for ii=1:length(fns)
    %% Name the output centerline
    name_split = strsplit(fns(ii).name, '.ply') ;
    name = name_split{1} ; 
    expstr = strrep(num2str(exponent, '%0.1f'), '.', 'p') ;
    outname = [fullfile(outdir, name) '_centerline_exp' expstr] ;
    polaroutfn = [fullfile(outdir, name) '_polarcoords'] ;
    skel_rs_outfn = [fullfile(outdir, name) '_centerline_scaled_exp' expstr ] ;
    fig1outname = [fullfile(fig1outdir, name) '_centerline_exp' expstr '_xy'] ;
    fig2outname = [fullfile(fig2outdir, name) '_centerline_exp' expstr '_xz'] ;
    fig3outname = [fullfile(fig3outdir, name) '_centerline_exp' expstr '_yz'] ;
    figsmoutname = [fullfile(fig3outdir, name) '_centerline_smoothed_exp' expstr] ;
    tmp = strsplit(name, '_') ;
    timestr = tmp{length(tmp)} ;
    
    %% Read the mesh
    mesh = ply_read(fullfile(fns(ii).folder, fns(ii).name));
    tri = cell2mat(mesh.face.vertex_indices) ;
    
    if strcmp(meshorder, 'zyx')
        xs = mesh.vertex.z / ssfactor ;
        ys = mesh.vertex.y / ssfactor ;
        zs = mesh.vertex.x / ssfactor ; 
    else
        error('Did not code for this order yet')
    end
    
    vertices = [xs, ys, zs] ;
    fv = struct('faces', tri + 1, 'vertices', vertices) ;

    % Must either downsample mesh, compute xyzgrid using ssfactor and
    % pass to options struct.
    % Here, downsampled mesh
    % mesh.vertex.x = xs ;
    % mesh.vertex.y = ys ;
    % mesh.vertex.z = zs ;

    %% Load the AP axis determination
    thres = 0.5 ;
    options.check = false ;
    apfn = fullfile(rootdir, ['Time_' timestr '_c1_stab_Probabilities_apcenterline.h5' ]);
    apdat = h5read(apfn, '/exported_data');
    
    rawfn = fullfile(rootdir, ['Time_' timestr '_c1_stab.h5' ]);
    rawdat = h5read(rawfn, '/inputData');
    adat = squeeze(apdat(1,:,:,:)) ;
    pdat = squeeze(apdat(2,:,:,:)) ;
    % testing
    % adat = 0 * adat ;
    % adat(1:20,1:20,1:20) = 1 ;
    % pdat = 0 * adat ;
    % pdat(1:20,1:60,1:20) = 1 ;
    % ddat = 0 * adat ;
    % ddat(1:20,1:20,1:300) = 1 ;
    % define axis order: 
    % if 1, 2, 3: axes will be yxz
    % if 1, 3, 2: axes will be yzx
    % if 2, 1, 3: axes will be xyz (ie first second third axes, ie --> 
    % so that bright spot at im(1,2,3) gives com=[1,2,3]
    axorder = [2, 1, 3] ;
    adat = permute(adat, axorder) ;
    pdat = permute(pdat, axorder) ;
    [aind, acom] = match_training_to_vertex(adat, thres, vertices, options) ;
    [pind, pcom] = match_training_to_vertex(pdat, thres, vertices, options) ;
    
    % Grab dorsal direction if this is the first timepoint
    if ii == 1    
        dorsal_thres = 0.9 ;
        ddat = permute(squeeze(apdat(4,:,:,:)), axorder) ;
        [dind, dcom] = match_training_to_vertex(ddat,...
            dorsal_thres, vertices, options) ;
        %%%%%%%%%%%%%%%%%%%%%%
        if preview
            % % disp('Showing dorsal segmentation...')
            % clf
            % for slice=1:2:size(ddat, 2)
            %     im = squeeze(ddat(:, slice, :)) ;
            %     % im(im < dorsal_thres) = 0 ;
            %     imshow(im)
            %     xlabel('x')
            %     ylabel('z')
            %     hold on
            %     plot(dcom(:, 1), dcom(:, 3), 'o')
            %     title([num2str(slice) '/' num2str(size(apdat, 3))])
            %     pause(0.001)
            % end
            %%%%%%%%%%%%%%%%%%%%%%
            fig = figure ;
            disp('Displaying mesh in figure ...')
            % iso = isosurface(rawdat, 880) ;
            % patch(iso,'facecolor',[1 0 0],'facealpha',0.1,'edgecolor','none');
            % view(3)
            % camlight
            % hold on;
            tmp = trimesh(fv.faces, ...
                vertices(:, 1), vertices(:,2), vertices(:, 3), ...
                vertices(:, 1)) % , 'edgecolor', 'none', 'FaceAlpha', 0.1) ;
            hold on;
            plot3(acom(1), acom(2), acom(3), 'ro')
            plot3(pcom(1), pcom(2), pcom(3), 'bo')
            plot3(dcom(1), dcom(2), dcom(3), 'ko')
            xlabel('x')
            ylabel('y')
            zlabel('z')
            axis equal
            %%%%%%%%%%%%%%%%%%%%%%
        end
        
        % compute rotation 
        apaxis = pcom - acom ;
        aphat = apaxis / norm(apaxis) ;
        
        % compute rotation matrix using this algorithm: 
        % https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
        xhat = [1, 0, 0] ;
        zhat = [0, 0, 1] ;
        ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0] ;
        RU = @(A,B) eye(3) + ssc(cross(A,B)) + ...
             ssc(cross(A,B))^2*(1-dot(A,B))/(norm(cross(A,B))^2) ;
        % rotz aligns AP to xhat (x axis)
        rotx = RU(aphat, xhat) ;
        
        % Rotate dorsal to the z axis
        % find component of dorsal vector from acom perp to AP
        dvec = rotx * (dcom - acom)' - rotx * (dot(dcom - acom, aphat) * aphat)' ;
        dhat = dvec / norm(dvec) ;
        rotz = RU(dhat, zhat) ;
        rot = rotz * rotx  ;
        
        % Save the rotation matrix
        disp(['Saving rotation matrix to txt: ', rotname, '.txt'])
        dlmwrite([rotname '.txt'], rot)
        
    end
    % acom = acom_ds * ssfactor ;
    % pcom = pcom_ds * ssfactor ;

    % Check if acom is inside mesh. If so, use that as starting point.
    ainside = inpolyhedron(fv, acom(1), acom(2), acom(3)) ;
    pinside = inpolyhedron(fv, pcom(1), pcom(2), pcom(3)) ;

    %% Define start point
    if ainside
        startpt = acom' ;
    else
        % move along the inward normal of the mesh from the matched vertex
        vtx = [vertices(aind, 1), vertices(aind, 2), vertices(aind, 3)]' ;
        normal = [mesh.vertex.nx(aind), ...
                    mesh.vertex.nx(aind), ...
                    mesh.vertex.nx(aind)]' ;
        startpt = vtx + normal;
        if ~inpolyhedron(fv, startpt(1), startpt(2), startpt(3)) 
            % this didn't work, check point in reverse direction
            startpt = vtx - normal * normal_step ;
            if ~inpolyhedron(fv, startpt(1), startpt(2), startpt(3))
                % Can't seem to jitter into the mesh, so use vertex
                disp("Can't seem to jitter into the mesh, so using vertex for startpt")
                startpt = vtx ;
            end
        end
    end 

    % Define end point
    if pinside
        endpt = pcom' ;
    else
        % move along the inward normal of the mesh from the matched vertex
        vtx = [vertices(pind, 1), vertices(pind, 2), vertices(pind, 3)]' ;
        normal = [mesh.vertex.nx(pind), ...
                    mesh.vertex.nx(pind), ...
                    mesh.vertex.nx(pind)]' ;
        endpt = vtx + normal * normal_step;
        if ~inpolyhedron(fv, endpt(1), endpt(2), endpt(3)) 
            % this didn't work, check point in reverse direction
            endpt = vtx - normal * normal_step ;
            if ~inpolyhedron(fv, endpt(1), endpt(2), endpt(3))
                % Can't seem to jitter into the mesh, so use vertex
                disp("Can't seem to jitter into the mesh, so using vertex for endpt")
                endpt = vtx ;
            end
        end
    end 

    % Check out the mesh
    if preview
        hold on
        trimesh(fv.faces, xs, ys, zs)
        % plot3(xs, ys, zs, 'ko')
        scatter3(startpt(1), startpt(2), startpt(3), 'ro')
        scatter3(endpt(1), endpt(2), endpt(3), 'ko')
        xlabel('x')
        ylabel('y')
        zlabel('z')
        hold off
        axis equal
    end

    % permute xy since MATLAB is evil
    % startpt = [startpt(1) startpt(2) startpt(3)]' ;
    % endpt = [endpt(1) endpt(2) endpt(3)]' ;

    %% Get xyz grid for distance transform
    if ii == 1
        xx = 0:res:ceil(max(xs) + buffer) ;
        yy = 0:res:ceil(max(ys) + buffer) ;
        zz = 0:res:ceil(max(zs) + buffer) ;
    end
    disp('Identifying points inside mesh...')
    tic 
    fv.faces = reorient_facets( fv.vertices, fv.faces );
    inside = inpolyhedron(fv, xx, yy, zz) ;
    outside = 1 - inside ;
    toc ; 

    % use the distanceTransform from Yuriy Mishchenko
    disp('Computing distance transform...')
    tic
    Dbw = bwdistsc(outside) ;
    % DD = max(DD(:)) - DD ;
    DD = (Dbw + eps) ./ (max(Dbw(:)) + eps) ;
    % DD = 1 - DD ;
    DD = DD.^(exponent) ; 
    DD(logical(outside)) = eps ;
    toc ; 

    if preview
        % Preview DD
        clf ;
        disp('Previewing the distance transform')
        for ll=1:3
            for kk=1:10:size(DD,1)
                imagesc(squeeze(DD(kk,:,:)))
                title(['DT: z=' num2str(kk)])
                colorbar
                pause(0.001)
            end
        end

        % A better way to plot it
        clf
        p = patch(isosurface(xx,yy,zz,inside,0.5));
        % % isonormals(x,y,z,v,p)
        p.FaceColor = 'red';
        p.EdgeColor = 'none';
        daspect([1 1 1])
        view(3); 
        axis tight
        camlight 
        % lighting gouraud
    end

    % Check points with subsampling
    % ssample = 10 ;
    % xp = X(:); yp=Y(:); zp=Z(:) ; dp=D(:);
    % scatter3(xp(1:ssample:end), yp(1:ssample:end), ...
    %          zp(1:ssample:end), 30, dp(1:ssample:end), 'filled') ;

    %% use Peyre's fast marcher
    disp('Computing centerline from fast marching...')
    tic

    % From example (DD is W, with low values being avoided)
    options.heuristic = weight * DD ;
    startpt_trans = [startpt(2), startpt(1), startpt(3)]' ;
    endpt_trans = [endpt(2), endpt(1), endpt(3)]' ;
    [D2,S] = perform_fast_marching(DD, startpt_trans, options);
    path = compute_geodesic(D2, endpt_trans);
    % plot_fast_marching_3d(D2, S, path, startpt, endpt);

    % Show the intermediate result
    disp('found skel')        
    if preview
        % Preview D2
        clf ;
        for kk=1:10:size(D2,3)
            imshow(squeeze(D2(kk,:,:)))
            title(['D2 for plane z=' num2str(kk)])
            pause(0.001)
        end

        % Preview S
        for kk=1:10:size(S,1)
            imshow(squeeze(S(kk,:,:)))
            title(['S for plane z=' num2str(kk)])
            pause(0.001)
        end
    end

    % options.weight = weight ;
    % [D2,S] = perform_fmstar_3d(DD, endpt, startpt) ;
    % toc ;
    % % Display the skeleton voxels
    % ind = find(skel);
    % [xs, ys, zs] = ind2sub(size(skel), ind);
    % plot3(ys, xs, zs, '.')

    % Convert skeleton's rows to columns and flip start/end
    skel_tmp = fliplr(path)' ;
    % Transpose x<->y back to original
    skel = [ skel_tmp(:,2), skel_tmp(:,1), skel_tmp(:,3) ];
    spt = [startpt(1), startpt(2), startpt(3)] ;
    ept = [endpt(1), endpt(2), endpt(3)] ;
    
    %% Compute the translation to put anterior to origin
    if ii == 1
        % Save translation 
        trans = -(rot * spt')' ;
        disp(['Saving translation vector (post rotation) to txt: ', transname, '.txt'])
        dlmwrite([transname '.txt'], trans)
    end
    
    %% Get axis limits if this is first TP
    if ii == 1
        xyzr = (rot * vertices')' ;
        xmin = min(xyzr(:, 1)) - plot_buffer + trans(1) ;
        ymin = min(xyzr(:, 2)) - plot_buffer + trans(2) ;
        zmin = min(xyzr(:, 3)) - plot_buffer + trans(3) ;
        xmax = max(xyzr(:, 1)) + plot_buffer + trans(1) ;
        ymax = max(xyzr(:, 2)) + plot_buffer + trans(2) ;
        zmax = max(xyzr(:, 3)) + plot_buffer + trans(3) ;
        
        % Save xyzlimits 
        disp('Saving xyzlimits for plotting')
        dlmwrite([xyzlimname '.txt'], [xmin, xmax; ymin, ymax; zmin, zmax])
    end
    
    %% Rotate and translate vertices and endpoints
    xyzr = (rot * vertices')' + trans ;
    skelr = (rot * skel')' + trans ; 
    sptr = (rot * spt')' + trans ; 
    eptr = (rot * ept')' + trans ;
    dptr = (rot * dcom')' + trans ; 
    
    %% Plot and save
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if save_figs
        disp('Saving rotated & translated figure ...')    
        close all
        fig = figure ;
        set(gcf, 'Visible', 'Off')
        tmp = trisurf(tri + 1, xyzr(:, 1), xyzr(:,2), xyzr(:, 3), ...
            xyzr(:, 1), 'edgecolor', 'none', 'FaceAlpha', 0.1) ;
        hold on;
        % plot the skeleton
        for i=1:length(skelr)
            plot3(skelr(:,1), skelr(:,2), skelr(:,3),'-','Color',[0,0,0], 'LineWidth', 3);
        end
        plot3(sptr(1), sptr(2), sptr(3), 'ro')
        plot3(eptr(1), eptr(2), eptr(3), 'bo')
        plot3(dptr(1), dptr(2), dptr(3), 'go')
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title(['Centerline using $D^{' num2str(exponent) '}$: ' timestr], ...
            'Interpreter', 'Latex')
        axis equal
        % xy
        view(2)
        xlim([xmin xmax])
        ylim([ymin ymax])
        zlim([zmin zmax])
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);
        saveas(fig, [fig1outname '.png'])
        % yz
        view(90, 0)
        xlim([xmin xmax])
        ylim([ymin ymax])
        zlim([zmin zmax])
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=15cm
        saveas(fig, [fig2outname '.png'])
        % xz
        view(0, 0)    
        xlim([xmin xmax])
        ylim([ymin ymax])
        zlim([zmin zmax])
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=15cm
        saveas(fig, [fig3outname '.png'])
        close all
    end
    
    % Save centerline as text file
    disp(['Saving centerline to txt: ', outname, '.txt'])
    dlmwrite([outname '.txt'], skel)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display the skeleton
    % disp('Saving figure ...')
    % close all
    % fig = figure ;
    % iso = isosurface(inside, 0.5) ;
    % patch(iso,'facecolor',[1 0 0],'facealpha',0.1,'edgecolor','none');
    % view(3)
    % camlight
    % hold on;
    % % plot the skeleton
    % for i=1:length(skel)
    %     plot3(skel(:,1), skel(:,2), skel(:,3),'-','Color',[0,0,0], 'LineWidth', 10);
    % end
    % plot3(spt(1), spt(2), spt(3), 'ro')
    % plot3(ept(1), ept(2), ept(3), 'bo')
    % xlabel('x')
    % ylabel('y')
    % zlabel('z')
    % axis equal
    % title(['$D^{' num2str(exponent) '}$'], 'Interpreter', 'Latex')
    % view(2)
    % xlim([xmin xmax])
    % ylim([ymin ymax])
    % zlim([zmin zmax])
    % saveas(fig, [fig1outname '.png'])
    % view(10)
    % stopping 
    % saveas(fig, [fig2outname '.png'])
    % saveas(fig, [fig3outname '.png'])
    % close all
    
    %% Associate each vertex with a point on the curve
    % dist2 = (xs - skel(:,1)).^2 + (xs - skel(:,2)).^2 + (xs - skel(:,3)).^2 ;
    [kmatch, dist] = dsearchn(skel, vertices) ;
    
    % A dsearchn() returns closest Euclidean 3D matches. CHeck that the
    % association linesegment between the vertex and the centerline does
    % not leave the mesh
    
    
    % get distance increment
    ds = vecnorm(diff(skel), 2, 2) ;
    % get pathlength at each skeleton point
    ss = [0; cumsum(ds)] ;
    
    % Check the associations
    if preview
        clf ; hold on
        for ijk = 1:500:length(vertices)
            plot3([vertices(ijk, 1) skel(kmatch(ijk), 1)], ...
                [vertices(ijk, 2) skel(kmatch(ijk), 2)], ...
                [vertices(ijk, 3) skel(kmatch(ijk), 3)])
        end
    end
    
    % Compute radius R(s)
    radii = vecnorm(vertices - skel(kmatch), 2, 2) ;
    
    % Compute phi(s), which is just the polar angle in the yz plane - pi/2
    % taken wrt the centerline
    phi_dorsal = mod(atan2(xyzr(:, 3) - sptr(3), xyzr(:, 2) - sptr(2)) - pi * 0.5, 2*pi);
    phi_ctrdorsal = mod(atan2(xyzr(:, 3) - skelr(kmatch, 3), ...
                        xyzr(:, 2) - skelr(kmatch, 2)) - pi * 0.5, 2*pi);

    % Check phi_dorsal
    if preview
        % view global dorsal angle
        fig = figure ;
        set(gcf, 'Visible', 'On')
        tmp = trisurf(tri + 1, xyzr(:, 1), xyzr(:,2), xyzr(:, 3), ...
            phi_dorsal, 'edgecolor', 'none', 'FaceAlpha', 0.1) ;
        hold on;
        
        % view dorsal angle from centerline
        fig = figure ;
        set(gcf, 'Visible', 'On')
        tmp = trisurf(tri + 1, xyzr(:, 1), xyzr(:,2), xyzr(:, 3), ...
            phi_ctrdorsal, 'edgecolor', 'none', 'FaceAlpha', 0.1) ;
        hold on;
    end

    % Save radii and angle wrt dorsal as text file
    disp(['Saving radii to txt: ', polaroutfn, '.txt'])
    dlmwrite([polaroutfn '.txt'], [kmatch, radii, phi_dorsal, phi_ctrdorsal])
    
    % Save the radius data as a plot
    % Color by phi_dorsal
    if save_figs
        close all
        fig = figure;
        set(gcf, 'Visible', 'Off')
        scatter(ss(kmatch) * resolution * ssfactor,...
            radii * resolution * ssfactor, [], ...
            phi_dorsal / pi, 'filled', ...
            'MarkerFaceAlpha', 0.05) ;        
        xlabel('pathlength, $s$ [$\mu$m]', 'Interpreter', 'Latex')
        ylabel('radius, $R$ [$\mu$m]', 'Interpreter', 'Latex')
        cb = colorbar() ;
        ylabel(cb, 'angle w.r.t. dorsal, $\phi / \pi$')
        cb.Label.Interpreter = 'latex';
        cb.Label.FontSize = 12 ;
        title('Midgut radius')
        xlim([0, 525])
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=16cm
        saveas(fig, fullfile(radius_vs_s_phi_outdir, [name '.png']))

        % Color by phi_ctrdorsal
        close all
        fig = figure;
        set(gcf, 'Visible', 'Off')
        scatter(ss(kmatch) * resolution * ssfactor,...
            radii * resolution * ssfactor, [], ...
            phi_dorsal / pi, 'filled', ...
            'MarkerFaceAlpha', 0.05) ;        
        xlabel('pathlength, $s$ [$\mu$m]', 'Interpreter', 'Latex')
        ylabel('radius, $R$ [$\mu$m]', 'Interpreter', 'Latex')
        cb = colorbar() ;
        ylabel(cb, 'angle w.r.t. dorsal, $\phi / \pi$')
        cb.Label.Interpreter = 'latex';
        cb.Label.FontSize = 12 ;
        title('Midgut radius')
        xlim([0, 525])
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=16cm
        saveas(fig, fullfile(radius_vs_s_phicd_outdir, [name '.png']))
        clf
    end
    
    % Get azimuthal angle, phi, from dorsal direction and centerline
    % Save cuts in phi as plots 
    % idx = (phi < eps) | (phi > (2*pi - eps)) ;
    % plot(ss(kmatch(idx)), 
    
    % Save the rotated, translated, scaled curve
    ss_s = ss * resolution * ssfactor ;
    skelr_s = skelr * resolution * ssfactor ;
    disp(['Saving rotated & scaled skeleton to txt: ', skel_rs_outfn, '.txt'])
    dlmwrite([skel_rs_outfn '.txt'], [ss_s, skelr_s])
    
    
end

disp('done')
