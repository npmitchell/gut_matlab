%% Extract the centerlines from a series of meshes (PLY files)
% Noah Mitchell 2019
% This version relies on the FastMarching toolbox, which may not compile 

%% First, compile required c code
% mex ./FastMarching_version3b/shortestpath/rk4.c
odir = pwd ;
cd('./FastMarching_version3b')
compile_c_files
cd(odir)

% Parameters
res = 1.

% Find all meshes to consider
rootpath = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/';
codepath = [rootpath, 'codes_matlab/'] ;
addpath(codepath)
addpath([codepath, 'ply_codes/']);
addpath([codepath, 'inpolyhedron/']);
addpath([codepath, 'FastMarching_version3b/']);
addpath([codepath, 'skeleton3d/']);
fns = dir(fullfile(rootpath, 'data/48Ygal4-UAShisRFP/20170329_objFiles/pointCloud_T*_ascii.ply'))

%% Iterate through each mesh
for ii=1:1 %length(fns)
    mesh = ply_read(fullfile(fns(ii).folder, fns(ii).name));
    xs = mesh.vertex.x ;
    ys = mesh.vertex.y ;
    zs = mesh.vertex.z ;
    tri = cell2mat(mesh.face.vertex_indices) ;
    fv = struct('faces', tri + 1, 'vertices', [xs, ys, zs]) 
    % reformat vertex_indices
    if ii == 1
        xx = floor(min(xs)):res:ceil(max(xs)) ;
        yy = floor(min(ys)):res:ceil(max(ys)) ;
        zz = floor(min(zs)):res:ceil(max(zs)) ;
        [X, Y, Z] = meshgrid(xx, yy, zz) ;
    end
    inside = inpolyhedron(fv, xx, yy, zz);
    
    
    skel = Skeleton3D(inside);
    % skel = skeleton(inside) ;
    
    disp('found skel')
    % Display the skeleton
    figure,
    iso = isosurface(inside, 0.5)
    patch(iso,'facecolor',[1 0 0],'facealpha',0.3,'edgecolor','none');
    view(3)
    camlight
    hold on;
    
    % Display the skeleton voxels
    ind = find(skel);
    [xs, ys, zs] = ind2sub(size(skel), ind);
    plot3(ys, xs, zs, '.')
    
    % Display the skeleton
    % for i=1:length(skel)
    %     L=S{i};
    %     plot3(L(:,2),L(:,1),L(:,3),'-','Color',rand(1,3));
    % end

end

