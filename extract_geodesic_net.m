%% Extract the centerlines from a series of meshes (PLY files)
% Noah Mitchell 2019
% This version relies on Gabriel Peyre's toolbox called
% toolbox_fast_marching/
% This code can be run from .../gut_matlab/

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
addpath([codepath, 'mesh_handling' filesep]);
addpath([codepath, 'inpolyhedron' filesep]);
codepath2 = [codepath 'toolbox_fast_marching' filesep ] ; 
codepath2 = [codepath2 'toolbox_fast_marching' filesep ] ; 
addpath([codepath2, 'toolbox' filesep]);
% addpath([codepath2, 'toolbox_graph_data' filesep]);
% addpath([codepath2, 'toolbox_graph_data/off' filesep]);

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
rootpath = '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/' ;
rootpath = [rootpath 'Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/'] ;
if ~exist(rootpath, 'dir')
    rootpath = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/' ;
    rootpath = [rootpath 'data/48Ygal4UasCAAXmCherry/201902072000_excellent/'] ;
end
meshdir = [rootpath 'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/'];
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
ii = 1 ;

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
    mesh.vertex.x = xs ;
    mesh.vertex.y = ys ;
    mesh.vertex.z = zs ;
    vertex = fv.vertices' ;
    faces = fv.faces' ;
    
    % % Check the mesh
    % clf;
    % plot_mesh(vertex, faces,options);
    
    %% Load the AP axis determination
    thres = 0.5 ;
    options.check = false ;
    apfn = [rootpath 'Time_000110_c1_stab_Probabilities_apcenterline.h5' ];
    apdat = h5read(apfn, '/exported_data');
    [aind, acom] = match_training_to_vertex(squeeze(apdat(1,:,:,:)), thres, vertex', options) ;
    [pind, pcom] = match_training_to_vertex(squeeze(apdat(2,:,:,:)), thres, vertex', options) ;
    % acom = acom_ds * ssfactor ;
    % pcom = pcom_ds * ssfactor ;

    %% Define start point & end point
    % start from the matched vertex
    startpt = vertex(:, aind) ;
    % Define end point from the matched vertex
    endpt = vertex(:, pind) ;
    
    %% front propagation on 3D meshes
    % name = 'elephant-50kv';
    % [vertex,faces] = read_mesh(name);
    % [vertex,faces] = read_mesh(fullfile(fns(ii).folder, fns(ii).name));
    % start_points = 24575;
    
    % Rescale vertices
    extents = max(vertex, [], 2) - min(vertex, [], 2) ;
    scale = max(extents) ;
    vertex = vertex / scale ;
    
    options.name = name;
    nverts = max(size(vertex));
    options.end_points = [];
    
    % Check out the mesh
    % scatter3(vertex(1, :), vertex(2, :), vertex(3, :))
    % Alternate mesh check
    plot_mesh(vertex, faces, options);
    hold on 
    scatter3(vertex(1, aind), vertex(2, aind), vertex(3, aind), 1000)
    
    disp('Performing propagation.');
    
    %   D is the distance function to the set of starting points.
    %   S is the final state of the points : -1 for dead (ie the distance
    %       has been computed), 0 for open (ie the distance is only a temporary
    %       value), 1 for far (ie point not already computed). Distance function
    %       for far points is Inf.
    %   Q is the index of the closest point. Q is set to 0 for far points.
    %       Q provide a Voronoi decomposition of the domain. 
    %%
    [D,S,Q] = perform_fast_marching_mesh(vertex, faces, aind, options);
    
    %% compute geodesics
    % npaths = 1;
    % [tmp,I] = sort( D(:) ); 
    % I = I(end:-1:1); 
    % I = I(1:round(nverts*1));
    % end_points = floor( rand(npaths,1)*(length(I)-1) )+1;
    % end_points = I(end_points);
    % % [tmp,I] = sort( D(:) ); end_points(1) = I(end);
    
    options.v2v = compute_vertex_ring(faces);
    options.e2f = compute_edge_face_ring(faces);
    
    disp('Extracting geodesics');
    % options.method = 'discrete';
    options.method = 'continuous';
    paths = compute_geodesic_mesh(D,vertex,faces, pind, options);
    
    options.colorfx = 'equalize';
    plot_fast_marching_mesh(vertex,faces, D, paths, options);
    
    % Save centerline as text file
    disp(['Saving geodesics to txt: ', outname, '.txt'])
    % dlmwrite([outname '.txt'], skel)

end
