%% Plot crossections from centerlines and a series of meshes (PLY files)
% Noah Mitchell 2019
% Run from the msls_output directory
% This version only works with simple cross sections made from slicing 
% aligned meshes along z=0 and y=0.
%
% First, run extract_centerline.m before running this code.
% Then run extract_geodesic_net.m before running this code.
% In this code, vtx has units of um, as does radii
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
addpath_recurse([codepath 'travelling_salesman' filesep]);
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
eps = 0.5 ;  % halfwidth of angular window to use for D,V,or Lateral 
meshorder = 'zyx' ;  % ordering of axes in loaded mesh wrt iLastik output
exponent = 1;  % exponent of DT used for velocity. Good values are ~1-2
anteriorChannel = 1;  % which channel of APD training is anterior
posteriorChannel = 2;  % which channel of APD training is posterior 
dorsalChannel = 4 ;  % which channel of APD training is dorsal
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

% Name output directory
cntrlinedir = [fullfile(meshdir, 'centerline') filesep ];
if ~exist(cntrlinedir, 'dir')
    mkdir(cntrlinedir) ;
end
% figure 1
choutdir = fullfile(meshdir, ['crosssections' filesep]);
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
    % Load the scaled and rotated mesh
    meshfn = fullfile(alignedmeshdir, [name '_APDV_um.ply']) ;
    tmp = strsplit(name, '_') ;
    timestr = tmp{length(tmp)} ;
    
    % Naming of DV crosssections 
    dv1_figfn = fullfile(figoutdir, ['dvsection_' timestr '.png'] ) ;
    dv2_figfn = fullfile(figoutdir, ['dvsection_fit_' timestr '.png'] ) ;
    lat1_figfn = fullfile(figoutdir, ['latsection_' timestr '.png'] ) ;
    lat2_figfn = fullfile(figoutdir, ['latsection_fit_' timestr '.png'] ) ;
    
    % Report progress
    msg = strrep(['Considering ' fns(ii).name], '_', '\_') ;
    if ii == 1
        fbar = waitbar(ii/length(fns), msg) ;
    else
        waitbar(ii/length(fns), fbar, msg)
    end
    
    % Load centerline 
    disp(['Loading centerline from txt: ', skel_rs_outfn, '.txt'])
    sskelrs = importdata([skel_rs_outfn '.txt']) ;
    ss = sskelrs(:, 1) ;
    ds = diff(ss) ; 
    ds = [ds; ds(length(ds))] ;
    skelrs = sskelrs(:, 2:4) ;
    
    % Load r,phi coordinates
    disp(['Loading radii from txt: ', polaroutfn, '.txt'])
    % Polar outfn has data [kmatch, radii, phi_dorsal, phi_ctrdorsal]
    dat = importdata([polaroutfn '.txt']) ;
    kmatch = dat(:, 1) ;  % indices of skelrs for each vertex
    radii = dat(:, 2) ;  % radii of each vertex
    phi_dorsal = dat(:, 3) ;  
    phi_ctrdorsal = dat(:, 4) ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot crosssections determined by y=0 and z=0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get positions of vertices in space (units = um, APDV aligned)
    mesh = read_ply_mod(meshfn) ;
    tri = mesh.f;
    vtx = mesh.v ;
    dv = (vtx(:, 2) < eps) & (vtx(:, 2) > -eps) ;
    lat = (vtx(:, 3) < eps) & (vtx(:, 3) > -eps) ;
    
    % Find good guess for TSP by indexing with angle from com
    dvcom = mean(vtx(dv, [1 3])) ;
    theta = atan2(vtx(dv, 3) - dvcom(2), vtx(dv, 1) - dvcom(1)) ;
    [~, dvsort] = sort(theta) ;
    
    latcom = mean(vtx(lat, [1 2])) ;
    theta = atan2(vtx(lat, 2) - dvcom(2), vtx(lat, 1) - dvcom(1)) ;
    [~, latsort] = sort(theta) ;
        
    % Find path by solving TSP
    %     USERCONFIG (structure) with zero or more of the following fields:
    %     - XY (float) is an Nx2 matrix of city locations, where N is the number of cities
    %     - DMAT (float) is an NxN matrix of point to point distances/costs
    %     - POPSIZE (scalar integer) is the size of the population (should be divisible by 4)
    %     - NUMITER (scalar integer) is the number of desired iterations for the algorithm to run
    %     - SHOWPROG (scalar logical) shows the GA progress if true
    %     - SHOWRESULT (scalar logical) shows the GA results if true
    %     - SHOWWAITBAR (scalar logical) shows a waitbar if true
    dvxy = vtx(dv, [1 3]) ;
    tspconfig.xy = dvxy(dvsort, :) ;
    tspconfig.numiter = 10000 ;
    tspconfig.showprog = preview ;
    tspconfig.showresult = false ;
    tspconfig.showwaitbar = false ;
    out = tsp_nn(tspconfig) ;
    dvorder = dvsort([out.optRoute, out.optRoute(1)]) ;
    % Now do lateral
    latxy = vtx(lat, [1 2]) ;
    tspconfig.xy = latxy(latsort, :) ;
    out = tsp_nn(tspconfig) ;
    latorder = latsort([out.optRoute, out.optRoute(1)]) ;
        
    %% Plot the realspace positions 
    fig = figure('Visible', 'Off') ;
    plot(dvxy(dvorder, 1), dvxy(dvorder, 2), '.')
    hold on;
    plot(skelrs(:, 1), skelrs(:, 3), 'k-')
    xlim([xmin xmax])
    ylim([zmin zmax])
    axis equal
    xlabel('x [$\mu$m]', 'Interpreter', 'Latex')
    ylabel('z [$\mu$m]', 'Interpreter', 'Latex')
    title(['DV cross section, $t=$',timestr], 'Interpreter', 'Latex')
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperPosition', [0 0 xwidth ywidth]);
    saveas(fig, dv1_figfn)
    % capture the entire figure as a frame
    drawnow
    F1(ii) = getframe(fig);
    cla
        
    % Plot the realspace positions
    % fig = figure('Visible', 'Off') ;
    plot(latxy(latorder, 1), latxy(latorder, 2), '.')
    plot(skelrs(:, 1), skelrs(:, 3), 'k-')
    axis equal
    xlim([xmin xmax])
    ylim([ymin ymax])
    xlabel('x [$\mu$m]', 'Interpreter', 'Latex')
    ylabel('y [$\mu$m]', 'Interpreter', 'Latex')
    title(['Lateral cross section, $t=$',timestr], 'Interpreter', 'Latex')
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperPosition', [0 0 xwidth ywidth]);
    saveas(fig, lat1_figfn)
    % capture the entire figure as a frame
    drawnow
    F2(ii) = getframe(fig);
    cla
    
    %% Plot the realspace positions with fits
    plot(dvxy(dvorder, 1), dvxy(dvorder, 2), '-')
    plot(skelrs(:, 1), skelrs(:, 3), 'k-')
    axis equal
    xlim([xmin xmax])
    ylim([zmin zmax])
    xlabel('x [$\mu$m]', 'Interpreter', 'Latex')
    ylabel('z [$\mu$m]', 'Interpreter', 'Latex')
    title(['DV cross section, $t=$',timestr], 'Interpreter', 'Latex')
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperPosition', [0 0 xwidth ywidth]);
    saveas(fig, dv2_figfn)
    % capture the entire figure as a frame
    drawnow
    F3(ii) = getframe(fig);
    cla
        
    % Plot the realspace positions
    % fig = figure('Visible', 'Off') ;
    plot(latxy(latorder, 1), latxy(latorder, 2), '-')
    plot(skelrs(:, 1), skelrs(:, 3), 'k-')
    axis equal
    xlim([xmin xmax])
    ylim([ymin ymax])
    xlabel('x [$\mu$m]', 'Interpreter', 'Latex')
    ylabel('y [$\mu$m]', 'Interpreter', 'Latex')
    title(['Lateral cross section, $t=$',timestr], 'Interpreter', 'Latex')
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperPosition', [0 0 xwidth ywidth]);
    saveas(fig, lat2_figfn)
    % capture the entire figure as a frame
    drawnow
    F4(ii) = getframe(fig);
    cla

end

close(fig);

%% Write to Video File ====================================================

% Write videos for radii versus pathlength of centerline
v = VideoWriter(fullfile(choutdir, 'dorsalventral_r_vs_s.avi'));
open(v);
for ii = 1:length(fns)
    writeVideo(v, F1(ii));    
    msg = ['Writing DV r(s): ' num2str(ii)] ;
    waitbar(ii/length(fns), fbar, msg)
end
close(v);
v = VideoWriter(fullfile(choutdir, 'lateral_r_vs_s.avi'));
open(v);
for ii = 1:length(fns)
    writeVideo(v, F2(ii));
    msg = ['Writing Lat r(s): ' num2str(ii)] ;
    waitbar(ii/length(fns), fbar, msg)
end
close(v);

% Write videos for raw data crosssections
v = VideoWriter(fullfile(choutdir, 'dorsalventral_xsection.avi'));
open(v);
for ii = 1:length(fns)
    writeVideo(v, F3(ii));
    msg = ['Writing Lat r(s): ' num2str(ii)] ;
    waitbar(ii/length(fns), fbar, msg)
end
close(v);
v = VideoWriter(fullfile(choutdir, 'lateral_xsection.avi'));
open(v);
for ii = 1:length(fns)
    writeVideo(v, F4(ii));
    msg = ['Writing Lat r(s): ' num2str(ii)] ;
    waitbar(ii/length(fns), fbar, msg)
end
close(v);

% Close and exit
close(fbar)
disp('done')
