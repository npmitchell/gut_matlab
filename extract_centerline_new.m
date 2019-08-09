%% Extract the centerlines from a series of meshes (PLY files)
% Noah Mitchell 2019
% This version relies on Gabriel Peyre's toolbox called
% toolbox_fast_marching/

%% First, compile required c code
% mex ./FastMarching_version3b/shortestpath/rk4
close all ;
odir = pwd ;
codepath = '/mnt/data/code/gut_matlab/' ;
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
res = 10 ;
buffer = 1 ;

% Find all meshes to consider
rootpath = '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/' ;
rootpath = [rootpath 'Time6views_60sec_1.4um_25x_obis1.5_2/data/'] ;
rootpath = [rootpath 'deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/'];
fns = dir(fullfile(rootpath, 'mesh_apical_0*.ply')) ;

%%
weights = [0.05];

%% Iterate through each mesh
for ii=1:1 %length(fns)
    %%
    mesh = ply_read(fullfile(fns(ii).folder, fns(ii).name));
    xs = mesh.vertex.x ;
    ys = mesh.vertex.y ;
    zs = mesh.vertex.z ;
    tri = cell2mat(mesh.face.vertex_indices) ;
    fv = struct('faces', tri + 1, 'vertices', [xs, ys, zs]) ;
    
    % Check out the mesh
    % trimesh(fv.faces, xs, ys, zs)
    % plot3(xs, ys, zs, '.')
    
    % reformat vertex_indices
    if ii == 1
        xx = floor(min(xs) - buffer):res:ceil(max(xs) + buffer) ;
        yy = floor(min(ys) - buffer):res:ceil(max(ys) + buffer) ;
        zz = floor(min(zs) - buffer):res:ceil(max(zs) + buffer) ;
        [X, Y, Z] = meshgrid(xx, yy, zz) ;
    end
    disp('Identifying points inside mesh...')
    tic 
    inside = inpolyhedron(fv, xx, yy, zz) ;
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
    Dbw = bwdistsc(inside) ;
    % D = max(D(:)) - D ;
    D = (Dbw + 1.0) ./ (max(Dbw(:)) + 1.0) ;
    D = sqrt(D) ; 
    %D = 1 - D ;
    toc ; 
    
    % Preview D
    disp('Previewing the distance transform')
    for ll=1:3
        for kk=1:size(D,1)
            imagesc(squeeze(D(kk,:,:)))
            title(['DT: z=' num2str(kk)])
            pause(0.001)
        end
    end
    
    % Check points with subsampling
    % ssample = 10 ;
    % xp = X(:); yp=Y(:); zp=Z(:) ; dp=D(:);
    % scatter3(xp(1:ssample:end), yp(1:ssample:end), ...
    %          zp(1:ssample:end), 30, dp(1:ssample:end), 'filled') ;
    
    % A better way to plot it
    % p = patch(isosurface(xx,yy,zz,D,10));
    % % isonormals(x,y,z,v,p)
    % p.FaceColor = 'red';
    % p.EdgeColor = 'none';
    % daspect([1 1 1])
    % view(3); 
    % axis tight
    % camlight 
    % lighting gouraud

    % Load start and endpoints
    % start_points = (1 / res) * [300, 1045, 238]' ;
    % start_points = [2, 3, 2]' + 0.5 * size(D)';
    % for coil
    % start_points = 0.5 * [90, 300, 50]' ;
    start_points = [25, 35, 25]';
    % start_points = [start_points(2), start_points(1), start_points(3) ]' ;
    % end_points = (1 / res) * [ 364, 239, 347 ]' ;
    % end_points = [11, -3, -2]' + 0.5 * size(D)' ;
    % for coil
    % end_points = 0.5 * [150, 40, 100]' ; 
    end_points = [70, 25, 25]';
    
    % transpose
    % end_points = [end_points(2), end_points(1), end_points(3) ]' ;
    
    % use Peyre's fast marcher
    disp('Computing centerline from fast marching...')
    tic
    
    %% try each weight in turn
    for jj = 1:length(weights)
        options.weight = weights(jj) ;
        [D2,S] = perform_fmstar_3d(D, start_points, end_points) ;
        toc ;

        disp('found skel')        

        % Preview D2
        for kk=1:size(D2,3)
            imshow(squeeze(D2(kk,:,:)))
            title(['D2 for plane z=' num2str(kk)])
            pause(0.001)
        end
        
        % Preview S
        for kk=1:size(S,1)
            imshow(squeeze(S(kk,:,:)))
            title(['S for plane z=' num2str(kk)])
            pause(0.001)
        end
                
       
        % % Display the skeleton voxels
        % ind = find(skel);
        % [xs, ys, zs] = ind2sub(size(skel), ind);
        % plot3(ys, xs, zs, '.')

        disp('Computing geodesic...')
        tic
        options.plot_planes = 0;
        path = compute_geodesic(D2,end_points);
        toc ;
        disp('Plotting result...')
        clf;
        plot_fast_marching_3d(D2,S,path, start_points, end_points, options);

        % view(AZ, 30);
        skel_tmp = path' ;
        % Transpose x<->y
        skel = [ skel_tmp(:,2), skel_tmp(:,1), skel_tmp(:,3) ];
        start_points = [start_points(2), start_points(1), start_points(3)] ;
        end_points = [end_points(2), end_points(1), end_points(3)] ;
        
        % Display the skeleton
        figure,
        iso = isosurface(inside, 0.5) ;
        patch(iso,'facecolor',[1 0 0],'facealpha',0.1,'edgecolor','none');
        view(3)
        camlight
        hold on;
        
        % Display the skeleton
        for i=1:length(skel)
            plot3(skel(:,1), skel(:,2), skel(:,3),'-','Color',rand(1,3));
        end
        plot3(start_points(1), start_points(2), start_points(3), 'ro')
        plot3(end_points(1), end_points(2), end_points(3), 'bo')
        xlabel('x')
        ylabel('y')
        zlabel('z')
        
        % kill script
        bad = good
    end
end

