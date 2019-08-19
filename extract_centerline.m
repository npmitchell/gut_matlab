%% Extract the centerlines from a series of meshes (PLY files)
% Noah Mitchell 2019
% This version relies on Gabriel Peyre's toolbox called
% toolbox_fast_marching/

%% First, compile required c code
% mex ./FastMarching_version3b/shortestpath/rk4
close all ;
odir = pwd ;
codepath = '/mnt/data/code/gut_matlab/' ;
if ~exist(codepath, 'dir')
    codepath = [pwd filesep] ;
end
addpath(codepath)
addpath([codepath 'addpath_recurse/']) ;
addpath([codepath, 'mesh_handling/']);
addpath([codepath, 'inpolyhedron/']);
% addpath([codepath, 'FastMarching_version3b/']);
% addpath([codepath, 'skeleton3d/']);

toolbox_path = [codepath 'toolbox_fast_marching/toolbox_fast_marching/'];
dtpath = [codepath 'distanceTransform/'] ;
addpath_recurse(toolbox_path)
addpath(dtpath)
% compile_c_files
cd(odir)

%% Parameters
res = 1 ;
buffer = 1 ;
ssfactor = 4; 
weight = 0.05;
normal_step = 2 ; 
preview = false ;
eps = 0.01 ;

% Find all meshes to consider
meshdir = odir ;
cd ../
rootdir = pwd ;
cd(odir) ;
% rootpath = '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/' ;
% rootpath = [rootpath 'Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/'] ;
% if ~exist(rootpath, 'dir')
%     rootpath = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/' ;
%     rootpath = [rootpath 'data/48Ygal4UasCAAXmCherry/201902072000_excellent/'] ;
% end
% meshdir = [rootpath 'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/'];
fns = dir(fullfile(meshdir, 'mesh_apical_0*.ply')) ;
ii = 1 ;
exponent = 1;

% Name output directory
outdir = [fullfile(fns(ii).folder, 'centerline') filesep ];
if ~exist(outdir, 'dir')
    mkdir(outdir) ;
end

% Find the xy limits
buffer = 5 ;
for ii=1:length(fns)
    if mod(ii, 10) == 0
        disp(['finding axis limits: ii = ', num2str(ii)])
        % disp(['xmin = ', num2str(xmin), ' / ymin = ', num2str(ymin)])
    end
    mesh = ply_read(fullfile(fns(ii).folder, fns(ii).name)); 

    if ii > 1
       xmin = min(xmin, min(mesh.vertex.x)) ;
       ymin = min(ymin, min(mesh.vertex.y)) ;
       xmax = max(xmax, max(mesh.vertex.x)) ;
       ymax = max(ymax, max(mesh.vertex.y)) ;
    else
       xmin = min(mesh.vertex.x) ;
       ymin = min(mesh.vertex.y) ;
       xmax = max(mesh.vertex.x) ;
       ymax = max(mesh.vertex.y) ;
    end
end
xmin = xmin / ssfactor - buffer ;
ymin = ymin / ssfactor - buffer ;
xmax = xmax / ssfactor + buffer ;
ymax = ymax / ssfactor + buffer ;

%% Iterate through each mesh
for ii=1:length(fns)
    %% Name the output centerline
    name_split = strsplit(fns(ii).name, '.ply') ;
    name = name_split{1} ; 
    expstr = strrep(num2str(exponent, '%0.1f'), '.', 'p') ;
    outname = [fullfile(outdir, name) '_centerline_exp' expstr] ;
    
    %% Read the mesh
    mesh = ply_read(fullfile(fns(ii).folder, fns(ii).name));
    tri = cell2mat(mesh.face.vertex_indices) ;
    xs = mesh.vertex.x / ssfactor ;
    ys = mesh.vertex.y / ssfactor ;
    zs = mesh.vertex.z / ssfactor ;

    % fv = struct('faces', tri + 1, 'vertices', ...
    %     [mesh.vertex.x, mesh.vertex.y, mesh.vertex.z]) ;    
    fv = struct('faces', tri + 1, 'vertices', [xs, ys, zs]) ;

    % Must either downsample mesh, compute xyzgrid using ssfactor and
    % pass to options struct.
    % Here, downsampled mesh
    % mesh_ds = mesh ;
    % mesh_ds.vertex.x = xs ;
    % mesh_ds.vertex.y = ys ;
    % mesh_ds.vertex.z = zs ;
    mesh.vertex.x = xs ;
    mesh.vertex.y = ys ;
    mesh.vertex.z = zs ;

    %% Load the AP axis determination
    thres = 0.5 ;
    options.check = false ;
    apfn = [rootpath 'Time_000110_c1_stab_Probabilities_apcenterline.h5' ];
    apdat = h5read(apfn, '/exported_data');
    [aind, acom] = match_training_to_vertex(squeeze(apdat(1,:,:,:)), thres, mesh, options) ;
    [pind, pcom] = match_training_to_vertex(squeeze(apdat(2,:,:,:)), thres, mesh, options) ;
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
        vtx = [mesh.vertex.x(aind), mesh.vertex.y(aind), mesh.vertex.z(aind)]' ;
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
        vtx = [mesh.vertex.x(pind), mesh.vertex.y(pind), mesh.vertex.z(pind)]' ;
        normal = [mesh.vertex.nx(pind), ...
                    mesh.vertex.nx(pind), ...
                    mesh.vertex.nx(pind)]' ;
        endpt = vtx + normal * normal_step;
        if ~inpolyhedron(fv, endpt(1), endpt(2), endpt(3)) 
            % this didn't work, check point in reverse direction
            endpt = vtx - normal ;
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
    end

    % permute xy since MATLAB is evil
    startpt = [startpt(2) startpt(1) startpt(3)]' ;
    endpt = [endpt(2) endpt(1) endpt(3)]' ;

    %% Get xyz grid for distance transform
    if ii == 1
        xx = 0:res:ceil(max(xs) + buffer) ;
        yy = 0:res:ceil(max(ys) + buffer) ;
        zz = 0:res:ceil(max(zs) + buffer) ;
        [X, Y, Z] = meshgrid(xx, yy, zz) ;
    end
    disp('Identifying points inside mesh...')
    tic 
    inside = inpolyhedron(fv, xx, yy, zz) ;
    outside = 1 - inside ;
    toc ; 

    % Make trivial example
    % inside = zeros(size(inside)) ;
    % inside(1:150, 50:100,30:60) = 1 ;

    % Old versions
    % skel = Skeleton3D(inside);
    % skel = skeleton(inside) ;

    % New version
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
                clim([0, 1])
                colorbar
                pause(0.001)
            end
        end
    end

    % Check points with subsampling
    % ssample = 10 ;
    % xp = X(:); yp=Y(:); zp=Z(:) ; dp=D(:);
    % scatter3(xp(1:ssample:end), yp(1:ssample:end), ...
    %          zp(1:ssample:end), 30, dp(1:ssample:end), 'filled') ;

    % A better way to plot it
    % p = patch(isosurface(xx,yy,zz,DD,10));
    % % isonormals(x,y,z,v,p)
    % p.FaceColor = 'red';
    % p.EdgeColor = 'none';
    % daspect([1 1 1])
    % view(3); 
    % axis tight
    % camlight 
    % lighting gouraud

    %% use Peyre's fast marcher
    disp('Computing centerline from fast marching...')
    tic

    % From example (DD is W, with low values being avoided)
    options.heuristic = weight * DD ;
    [D2,S] = perform_fast_marching(DD, startpt, options);
    path = compute_geodesic(D2, endpt);
    plot_fast_marching_3d(D2, S, path, startpt, endpt);

    % Show the intermediate result
    disp('found skel')        
    if preview
        % Preview D2
        clf
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

    % Old version 20190809
    % disp('Computing geodesic...')
    % tic
    % options.plot_planes = 0;
    % path = compute_geodesic(D2,endpt);
    % toc ;
    % disp('Plotting result...')
    % clf;
    % plot_fast_marching_3d(D2,S,path, endpt, startpt, options);

    % view(AZ, 30);
    skel_tmp = path' ;
    % Transpose x<->y back to original
    skel = [ skel_tmp(:,2), skel_tmp(:,1), skel_tmp(:,3) ];
    spt = [startpt(2), startpt(1), startpt(3)] ;
    ept = [endpt(2), endpt(1), endpt(3)] ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display the skeleton
    close all
    fig = figure ;
    iso = isosurface(inside, 0.5) ;
    patch(iso,'facecolor',[1 0 0],'facealpha',0.1,'edgecolor','none');
    view(3)
    camlight
    hold on;
    % plot the skeleton
    for i=1:length(skel)
        plot3(skel(:,1), skel(:,2), skel(:,3),'-','Color',[0,0,0], 'LineWidth', 10);
    end
    plot3(spt(1), spt(2), spt(3), 'ro')
    plot3(ept(1), ept(2), ept(3), 'bo')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    title(['$D^{' num2str(exponent) '}$'], 'Interpreter', 'Latex')
    view(2)
    xlim([xmin xmax])
    ylim([ymin ymax])
    saveas(fig, [outname '.png'])
    close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Save centerline as text file
    disp(['Saving centerline to txt: ', outname, '.txt'])
    dlmwrite([outname '.txt'], skel)

end
