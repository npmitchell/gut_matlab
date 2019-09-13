%% Extract the centerlines from a series of meshes (PLY files)
% Noah Mitchell 2019
% This version relies on Gabriel Peyre's toolbox called
% toolbox_fast_marching/
% 
% Prerequisites
% -------------
% extract_centerline.m
% slice_mesh_endcaps.m
%
% Run this script from the meshdir.
% First run slice_mesh_endcaps.m before running this code.

%% First, compile required c code
% mex ./FastMarching_version3b/shortestpath/rk4
clear ;
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
overwrite = false ;  % recompute centerline
save_figs = true ;  % save images of cntrline, etc, along the way
preview = false ;  % display intermediate results
res = 1 ;  % pixels per gridspacing of DT for cntrline extraction
resolution = 0.2619 ;  % um per pixel for full resolution (not subsampled)
dorsal_thres = 0.9 ;  % threshold for extracting Dorsal probability cloud 
buffer = 5 ;  % extra space in meshgrid of centerline extraction, to ensure mesh contained in volume
plot_buffer = 30; 
ssfactor = 4;  % subsampling factor for the h5s used to train for mesh/acom/pcom/dcom
weight = 0.1;  % for speedup of centerline extraction. Larger is less precise
normal_step = 0.5 ;  % how far to move normally from ptmatched vtx if a/pcom is not inside mesh
eps = 0.05 ;  % halfwidth of angular window to use for D,V,or Lateral 
meshorder = 'zyx' ;  % ordering of axes in loaded mesh wrt iLastik output
exponent = 1;  % exponent of DT used for velocity. Good values are ~1-2
axorder = [2, 1, 3] ;  % axis order for APD training output
% figure parameters
xwidth = 16 ; % cm
ywidth = 10 ; % cm

% Find all meshes to consider
meshdir = pwd ;
cd ../
rootdir = pwd ;
cd(meshdir)

% Directories
alignedmeshdir = fullfile(meshdir, ['aligned_meshes' filesep]) ;
if ~exist(alignedmeshdir, 'dir')
    mkdir(alignedmeshdir) ;
end
fns = dir(fullfile(meshdir, 'mesh_apical_stab_0*.ply')) ;
rotname = fullfile(meshdir, 'rotation_APDV') ;
transname = fullfile(meshdir, 'translation_APDV') ;
xyzlimname = fullfile(meshdir, 'xyzlim_APDV') ;
cylmeshdir = fullfile(meshdir, ['cylindercut' filesep]) ;
outapphicd_fn = fullfile(cylmeshdir, 'ap_boundary_phicd_values.h5') ;
boundaryfn = fullfile(cylmeshdir, 'ap_boundary_indices.h5') ;

% Name output directory
cntrlinedir = [fullfile(meshdir, 'centerline') filesep ];
if ~exist(cntrlinedir, 'dir')
    mkdir(cntrlinedir) ;
end
% figure 1
choutdir = fullfile(meshdir, ['geodesic_net' filesep]);
if ~exist(choutdir, 'dir')
    mkdir(choutdir) ;
end
figoutdir = fullfile(choutdir, ['images' filesep]);
if ~exist(figoutdir, 'dir')
    mkdir(figoutdir) ;
end

%% Load transformations
% Load the rotation matrix
rot = importdata([rotname '.txt']) ;
% Load the translation to put anterior to origin
trans = importdata([transname '.txt']) ;
% Load plotting limits
xyzlim = importdata([xyzlimname '.txt']) ;
xmin = xyzlim(1) * resolution ;
ymin = xyzlim(2) * resolution ;
zmin = xyzlim(3) * resolution ;
xmax = xyzlim(4) * resolution ;
ymax = xyzlim(5) * resolution ;
zmax = xyzlim(6) * resolution ;

%% Iterate through each mesh to extract crosssection
ii = 1 ;
F1(length(fns)) = struct('cdata', [], 'colormap', [] );
F2(length(fns)) = struct('cdata', [], 'colormap', [] );
F3(length(fns)) = struct('cdata', [], 'colormap', [] );
F4(length(fns)) = struct('cdata', [], 'colormap', [] );
%%
for ii=1:length(fns)
    %% Name the output centerline
    name_split = strsplit(fns(ii).name, '.ply') ;
    name = name_split{1} ; 
    expstr = strrep(num2str(exponent, '%0.1f'), '.', 'p') ;
    resstr = strrep(num2str(res, '%0.1f'), '.', 'p') ;
    extenstr = ['_exp' expstr '_res' resstr] ;
    outname = [fullfile(cntrlinedir, name) '_centerline' extenstr] ;
    polaroutfn = [fullfile(cntrlinedir, name) '_polarcoords' extenstr] ;
    skel_rs_outfn = [fullfile(cntrlinedir, name) '_centerline_scaled' extenstr ] ;
    tmp = strsplit(name, '_') ;
    timestr = tmp{length(tmp)} ;
    
    % Load the cylinder mesh
    meshfn = fullfile(cylmeshdir, [name '_cylindercut.ply']) ;
    % Load the scaled and rotated mesh
    % meshfn = fullfile(alignedmeshdir, [name '_APDV_um.ply']) ;
    
    %% Read the mesh
    mesh = read_ply_mod(meshfn);
    vtx = mesh.v' ;
    faces = mesh.f ;

    % Load anterior and posterior indices and their phicd values
    ab = h5read(boundaryfn, ['/' name '/anterior_boundary_indices']) ;
    pb = h5read(boundaryfn, ['/' name '/posterior_boundary_indices']) ;
    aphi = h5read(outapphicd_fn, ['/' name '/aphicd']) ;
    pphi = h5read(outapphicd_fn, ['/' name '/pphicd']) ;
    aphi(aphi > pi) = aphi(aphi > pi) - 2 * pi ;
    pphi(pphi > pi) = pphi(pphi > pi) - 2 * pi ;
    
    % Reorder indices to start with phi=0 and go to phi=2*pi
    % find where phicd is smallest
    [~, indx0a] = min(abs(aphi)) ;
    [~, indx0p] = min(abs(pphi)) ;
    aphi = circshift(aphi, -indx0a) ;
    pphi = circshift(pphi, -indx0p) ;
    ab = circshift(ab, -indx0a) ;
    pb = circshift(pb, -indx0p) ;
    
    % % Check the mesh
    % clf;
    % plot_mesh(vertex, faces,options);
    
    %% For each endpoint on posterior determine a suitable startpoint
    % for 
    % index of boundary vertex lists to use for startpt and endpt
    if length(pb) < length(ab)
        pind = pb(1:2:end) ;
        ni = length(pind) ;
        aind = flipud(ab(round(linspace(1, length(ab), ni)))); 
    else
        aind = ab(1:2:end) ;
        ni = length(aind) ;
        pind = flipud(pb(round(linspace(1, length(pb), ni)))) ;
    end
    
    %% Define start point & end point
    % start from the matched vertex
    startpt = vtx(:, aind) ;
    % Define end point from the matched vertex
    endpt = vtx(:, pind) ;
    
    %% front propagation on 3D meshes
    % name = 'elephant-50kv';
    % [vtx, faces] = read_mesh(name);
    % [vtx, faces] = read_mesh(fullfile(fns(ii).folder, fns(ii).name));
    % start_points = 24575;
    
    % Rescale vertices
    extents = max(vtx, [], 2) - min(vtx, [], 2) ;
    scale = max(extents) ;
    vtx = vtx / scale ;
        
    % % Check out the mesh
    % scatter3(vtx(1, :), vtx(2, :), vtx(3, :))
    % % Alternate mesh check
    % plot_mesh(vtx, faces, options);
    % shading interp
    % hold on 
    % scatter3(vtx(1, aind), vtx(2, aind), vtx(3, aind), 100)
    % scatter3(vtx(1, pind), vtx(2, pind), vtx(3, pind), 100)
    
    %% Compute each pair of geodesics
    paths = {} ;
    for jj = 1:length(aind)
        
        % Options for the marching
        startpts = aind(jj) ;
        options.name = name;
        nverts = max(size(vtx));
        options.end_points = [];
        
        msg = ['Performing propagation for pair ' num2str(jj) ] ; 
        if exist('fbar', 'var') 
            if isvalid(fbar)
                waitbar(jj / length(aind), fbar, msg);
            else
                fbar = waitbar(jj / length(aind), msg) ;
            end
        else
            fbar = waitbar(jj / length(aind), msg) ;
        end

        %   D is the distance function to the set of starting points.
        %   S is the final state of the points : -1 for dead (ie the distance
        %       has been computed), 0 for open (ie the distance is only a temporary
        %       value), 1 for far (ie point not already computed). Distance function
        %       for far points is Inf.
        %   Q is the index of the closest point. Q is set to 0 for far points.
        %       Q provide a Voronoi decomposition of the domain. 
        %
        [D,S,Q] = perform_fast_marching_mesh(vtx, faces, startpts, options);

        % compute geodesics
        % npaths = 1;
        % [tmp,I] = sort( D(:) ); 
        % I = I(end:-1:1); 
        % I = I(1:round(nverts*1));
        % end_points = floor( rand(npaths,1)*(length(I)-1) )+1;
        % end_points = I(end_points);
        % % [tmp,I] = sort( D(:) ); end_points(1) = I(end);

        options.v2v = compute_vertex_ring(faces);
        options.e2f = compute_edge_face_ring(faces);

        msg = ['Extracting geodesics for pair ' num2str(jj) ] ; 
        waitbar(jj / length(aind), fbar, msg);
        
        endpts = pind(jj) ;
        % options.method = 'discrete';
        options.method = 'continuous';
        options.verb = true ;
        path = compute_geodesic_mesh(D, vtx, faces, endpts, options);
    
        paths{jj} = path ;
        % Check it
        % for jj = 1:length(paths)
        %     path = paths{jj}' ;
        %     plot3(path(:, 1), path(:, 2), path(:, 3))
        % end
    end
    close(fbar)
    % Save centerline as text file
    disp(['Saving geodesics to txt: ', outname, '.txt'])
    % dlmwrite([outname '.txt'], skel)

    options.colorfx = 'equalize';
    options.edgecolor = 'none' ;
    plot_fast_marching_mesh(vtx, faces, D, paths, options);
    
    saveas(gcf, fullfile(figoutdir, [name '.png']))
    close all
end
