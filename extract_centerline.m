%% Extract the centerlines from a series of meshes (PLY files) 
% Noah Mitchell 2019
% Also saves scaled meshes with AP along x and DV along y, centered at A=0.
% This version relies on Gabriel Peyre's toolbox called
% toolbox_fast_marching/
%
% Run from the msls_output directory
% Run this code only after training on anterior (A), posterior (P), and 
% dorsal anterior (D) points in different iLastik channels.
% anteriorChannel, posteriorChannel, and dorsalChannel specify the iLastik
% training channel that is used for each specification.
% Name the h5 file output from iLastik as ..._Probabilities_apcenterline.h5
% Train for anterior dorsal (D) only at the first time point, because
% that's the only one that's used.
%
% Note that as of 20191029, the output meshes are mirrored due to the axis
% permutation from iLastik
%
% OUTPUTS
% -------
% xyzlim.txt 
%   xyzlimits of raw meshes in units of full resolution pixels (ie not
%   downsampled)
% xyzlim_APDV.txt 
%   xyzlimits of rotated and translated meshes in units of full resolution 
%   pixels (ie not downsampled)
% xyzlim_APDV_um.txt 
%   xyz limits of rotated and translated meshes in microns
% rotation_APDV.txt
%   rotation matrix to align mesh to APDV frame
% translation_APDV.txt
%   translation vector to align mesh to APDV frame
% xyzlim.txt 
%   raw bounding box in original frame (not rotated), in full res pixels
% xyzlim_APDV.txt
%   bounding box in rotated frame, in full resolution pixels
% xyzlim_APDV_um.txt
%   bounding box in rotated frame, in microns
% apdv_coms_rs.h5
%   Centers of mass for A, P, and D in microns in rotated APDV coord system
% 
% Performs centerline extraction in subsampled units, where ssfactor is the
% same as used for iLastik training to get A, P, and D.
% trans has units of mesh coordinates.
% phi is the angle in the xy plane
% vertices are in units of pixels (at full resolution)
% skel, spt, and ept are all in units of mesh coordinates (at full res)
% startpt, endpt are in subsampled units
% radii are in microns
% To take mesh to rotated + translated mesh in physical units, apply:
%         xs = mesh.vertex.z ;
%         ys = mesh.vertex.y ;
%         zs = mesh.vertex.x ;
%         vtx_rs = (rot * vtx' + trans)' * resolution
%         
% See also
% --------
% To run first:
%    Gut_Pipeline.m
%    align_meshes_APDV.m
% To run after:
%    slice_mesh_endcaps.m
%    extract_chirality_writhe.m
%    Generate_Axisymmetric_Pullbacks_Orbifold.m
%    
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
addpath_recurse([codepath 'mesh_handling' filesep]);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overwrite = false ;  % recompute centerline
overwrite_apdvcoms = false ;  % recompute APDV coms from training
save_figs = true ;  % save images of cntrline, etc, along the way
overwrite_ims = false ;  % overwrite images even if centerlines are not overwritten
preview = false ;  % display intermediate results, for debugging
res = 1 ;  % pixels per gridspacing of DT for cntrline extraction, in units of subsampled pixels
resolution = 0.2619 ;  % um per pixel for full resolution (not subsampled)
dorsal_thres = 0.9 ;  % threshold for extracting Dorsal probability cloud 
buffer = 5 ;  % extra space in meshgrid of centerline extraction, to ensure mesh contained in volume
plot_buffer = 40;  % extra space in plots, in um
ssfactor = 4;  % subsampling factor for the h5s used to train for mesh/acom/pcom/dcom
weight = 0.1;  % for speedup of centerline extraction. Larger is less precise
normal_step = 0.5 ;  % how far to move normally from ptmatched vtx if a/pcom is not inside mesh
eps = 0.01 ;  % value for DT outside of mesh in centerline extraction
meshorder = 'zyx' ;  % ordering of axes in loaded mesh wrt iLastik output
exponent = 1;  % exponent of DT used for velocity. Good values are ~1-2
anteriorChannel = 1;  % which channel of APD training is anterior
posteriorChannel = 2;  % which channel of APD training is posterior 
dorsalChannel = 4 ;  % which channel of APD training is dorsal
axorder = [2, 1, 3] ;  % axis order for APD training output
% figure parameters
xwidth = 16 ; % cm
ywidth = 10 ; % cm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Name output directory
outdir = [fullfile(meshdir, 'centerline') filesep ];
if ~exist(outdir, 'dir')
    mkdir(outdir) ;
end
fns = dir(fullfile(meshdir, 'mesh_apical_stab_0*.ply')) ;
% Ensure that PLY files exist
if isempty(fns)
    error('Found no matching PLY files in ' + meshdir)
end

rotname = fullfile(meshdir, 'rotation_APDV') ;
transname = fullfile(meshdir, 'translation_APDV') ;
xyzlimname_raw = fullfile(meshdir, 'xyzlim') ;
xyzlimname = fullfile(meshdir, 'xyzlim_APDV') ;
xyzlimname_um = fullfile(meshdir, 'xyzlim_APDV_um') ;
outapdvname = fullfile(outdir, 'apdv_coms_rs.h5') ;
dcomname = fullfile(outdir, 'dorsalcom.txt') ;

% Name the directory for outputting figures
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
% Figures for phi_dorsal and phi_cntrdorsal
phi_def_outdir = [figoutdir 'phid_definition' filesep];
if ~exist(phi_def_outdir, 'dir')
    mkdir(phi_def_outdir) ;
end

% No longer saving radius 
% radius_vs_s_phi_outdir = [figoutdir 'radius_vs_s_phid' filesep];
% if ~exist(radius_vs_s_phi_outdir, 'dir')
%     mkdir(radius_vs_s_phi_outdir) ;
% end
% radius_vs_s_phicd_outdir = [figoutdir 'radius_vs_s_phicd' filesep];
% if ~exist(radius_vs_s_phicd_outdir, 'dir')
%     mkdir(radius_vs_s_phicd_outdir) ;
% end

phicd_def_outdir = [figoutdir 'phicd_definition' filesep];
if ~exist(phicd_def_outdir, 'dir')
    mkdir(phicd_def_outdir) ;
end
alignedmeshdir = fullfile(meshdir, ['aligned_meshes' filesep]) ;
if ~exist(alignedmeshdir, 'dir')
    mkdir(alignedmeshdir) ;
end
ii = 1 ;

%% Iterate through each mesh to compute acom(t) and pcom(t). Prepare file.
acoms = zeros(length(fns), 3) ;
pcoms = zeros(length(fns), 3) ;
timepts = zeros(length(fns)) ;
rawapdvname = fullfile(outdir, 'apdv_coms_from_training.h5') ;
load_from_disk = false ;
if exist(rawapdvname, 'file')
    load_from_disk = true ;
    try
        h5create(rawapdvname, ['/' name '/acom_sm'], size(pcom)) ;
        load_from_disk = false ;
    catch
        try
            acom_sm = h5read(rawapdvname, '/acom_sm') ;
            disp('acom_sm already exists')
        catch
            load_from_disk = false;
        end
    end
    try
        h5create(rawapdvname, ['/' name '/pcom_sm'], size(pcom)) ;
        load_from_disk = false ;
    catch
        try
            pcom_sm = h5read(rawapdvname, '/pcom_sm') ;
            disp('pcom_sm already exists')
        catch
            load_from_disk = false;
        end
    end
end
if ~load_from_disk
    disp('acom and pcom not already saved on disk. Compute them')
end


%% Compute acom and pcom if not loaded from disk
if ~load_from_disk || overwrite_apdvcoms
    for ii=1:length(fns)
        %% Get the timestamp string from the name of the mesh
        name_split = strsplit(fns(ii).name, '.ply') ;
        name = name_split{1} ; 
        tmp = strsplit(name, '_') ;
        timestr = tmp{length(tmp)} ;

        %% Load the AP axis determination
        if ~exist('fbar', 'var')
            fbar = waitbar(ii / length(fns), ['Computing acom, pcom for ' timestr ]) ;
        else
            if ~isvalid(fbar)
              fbar = waitbar(ii / length(fns), ['Computing acom, pcom for ' timestr ]) ;
            end
        end
        waitbar(ii / length(fns), fbar, ['Computing acom, pcom for ' timestr ])
        thres = 0.5 ;
        options.check = false ;
        apfn = fullfile(rootdir, ['Time_' timestr '_c1_stab_Probabilities_apcenterline.h5' ]);
        apdat = h5read(apfn, '/exported_data');

        % rawfn = fullfile(rootdir, ['Time_' timestr '_c1_stab.h5' ]);
        % rawdat = h5read(rawfn, '/inputData');
        adat = squeeze(apdat(anteriorChannel,:,:,:)) ;
        pdat = squeeze(apdat(posteriorChannel,:,:,:)) ;
        
        % define axis order: 
        % if 1, 2, 3: axes will be yxz
        % if 1, 3, 2: axes will be yzx
        % if 2, 1, 3: axes will be xyz (ie first second third axes, ie --> 
        % so that bright spot at im(1,2,3) gives com=[1,2,3]
        adat = permute(adat, axorder) ;
        pdat = permute(pdat, axorder) ;
        acom = com_region(adat, thres, options) ;
        pcom = com_region(pdat, thres, options) ;
        % [~, acom] = match_training_to_vertex(adat, thres, vertices, options) ;
        % [~, pcom] = match_training_to_vertex(pdat, thres, vertices, options) ;
        acoms(ii, :) = acom ;
        pcoms(ii, :) = pcom ; 
        timepts(ii) = str2double(timestr) ;
    end
    if isvalid(fbar)
        close(fbar)
    end
    disp('done determining acoms, pcoms')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Smooth the acom and pcom data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Smoothing acom and pcom...')
    timepts = linspace(0, length(fns) - 1, length(fns)) ;
    acom_sm = 0 * acoms ;
    pcom_sm = 0 * acoms ;
    smfrac = 30 / length(timepts) ;  % fraction of data for smoothing window
    acom_sm(:, 1) = smooth(timepts, acoms(:, 1), smfrac, 'rloess');
    pcom_sm(:, 1) = smooth(timepts, pcoms(:, 1), smfrac, 'rloess');
    acom_sm(:, 2) = smooth(timepts, acoms(:, 2), smfrac, 'rloess');
    pcom_sm(:, 2) = smooth(timepts, pcoms(:, 2), smfrac, 'rloess');
    acom_sm(:, 3) = smooth(timepts, acoms(:, 3), smfrac, 'rloess');
    pcom_sm(:, 3) = smooth(timepts, pcoms(:, 3), smfrac, 'rloess');
    
    if preview
        plot(timepts, acoms - mean(acoms,1), '.')
        hold on
        plot(timepts, acom_sm - mean(acoms, 1), '-')
        title('Smoothed COMs for AP')
    end
    clear acom pcom
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save smoothed anterior and posterior centers of mass ===============
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        h5create(rawapdvname, '/acom', size(acoms)) ;
    catch
        disp('acom already exists')
    end
    try
        h5create(rawapdvname, '/pcom', size(pcoms)) ;
    catch
        disp('pcom already exists')
    end
    try
        h5create(rawapdvname, '/acom_sm', size(acom_sm)) ;
    catch
        disp('acom_sm already exists')
    end
    try
        h5create(rawapdvname, '/pcom_sm', size(pcom_sm)) ;
    catch
        disp('pcom_sm already exists')
    end
    h5write(rawapdvname, '/acom', acoms) ;
    h5write(rawapdvname, '/pcom', pcoms) ;
    h5write(rawapdvname, '/acom_sm', acom_sm) ;
    h5write(rawapdvname, '/pcom_sm', pcom_sm) ;
    clear acoms pcoms
else
    disp('Skipping, since already loaded acom_sm and pcom_sm')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get axis limits from looking at all meshes =============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn = [xyzlimname_raw '.txt'] ;
if exist(fn, 'file')
    disp('Loading xyzlimits from disk...')
    xyzlims = dlmread(fn, ',', 1, 0);
    xmin = xyzlims(1);
    ymin = xyzlims(2);
    zmin = xyzlims(3);
    xmax = xyzlims(4);
    ymax = xyzlims(5);
    zmax = xyzlims(6);
else
    disp('Extracting xyzlimits for raw meshes...')
    for ii = 1:length(fns)
        % Get the timestamp string from the name of the mesh
        mesh = read_ply_mod(fns(ii).name) ;

        minx = min(mesh.v) ;
        maxx = max(mesh.v) ;
        if ii == 1
            xmin = minx(1) ;
            ymin = minx(2) ;
            zmin = minx(3) ;
            xmax = maxx(1) ;
            ymax = maxx(2) ;
            zmax = maxx(3) ;
        else
            xmin= min(xmin, minx(1)) ;
            ymin = min(ymin, minx(2)) ;
            zmin = min(zmin, minx(3)) ;
            xmax = max(xmax, maxx(1)) ;
            ymax = max(ymax, maxx(2)) ;
            zmax = max(zmax, maxx(3)) ;
        end 
    end

    % Save xyzlimits 
    disp('Saving raw mesh xyzlimits for plotting')
    header = 'xyzlimits for original meshes in units of full resolution pixels' ; 
    write_txt_with_header(fn, [xmin, xmax; ymin, ymax; zmin, zmax], header) ;
end

% Subsample the xyzlimits
% xmin = xmin / ssfactor; xmax = xmax / ssfactor;
% ymin = ymin / ssfactor; ymax = ymax / ssfactor;
% zmin = zmin / ssfactor; zmax = zmax / ssfactor;
disp('done')

%% With acom and pcom in hand, we compute dorsal and centerlines ==========
xminrs = 0 ; xmaxrs = 0;
yminrs = 0 ; ymaxrs = 0;
zminrs = 0 ; zmaxrs = 0;
for ii=1:length(fns)
    % Pick out the acom and pcom in SUBSAMPLED UNITS from smoothed sequence
    acom = acom_sm(ii, :) ;
    pcom = pcom_sm(ii, :) ; 
    
    %% Name the output centerline
    name_split = strsplit(fns(ii).name, '.ply') ;
    name = name_split{1} ; 
    expstr = strrep(num2str(exponent, '%0.1f'), '.', 'p') ;
    resstr = strrep(num2str(res, '%0.1f'), '.', 'p') ;
    extenstr = ['_exp' expstr '_res' resstr] ;
    outname = [fullfile(outdir, name) '_centerline' extenstr] ;
    polaroutfn = [fullfile(outdir, name) '_polarcoords' extenstr] ;
    skel_rs_outfn = [fullfile(outdir, name) '_centerline_scaled' extenstr ] ;
    fig1outname = [fullfile(fig1outdir, name) '_centerline' extenstr '_xy'] ;
    fig2outname = [fullfile(fig2outdir, name) '_centerline' extenstr '_xz'] ;
    fig3outname = [fullfile(fig3outdir, name) '_centerline' extenstr '_yz'] ;
    figsmoutname = [fullfile(fig3outdir, name) '_centerline_smoothed' extenstr] ;
    tmp = strsplit(name, '_') ;
    timestr = tmp{length(tmp)} ;
    
    %% Read the mesh  
    msg = strrep(['Loading mesh ' fns(ii).name], '_', '\_') ;
    if exist('fbar', 'var')
        if isvalid(fbar)
            waitbar(ii/length(fns), fbar, msg)
        else
            fbar = waitbar(ii/length(fns), msg) ;
        end
    else
        fbar = waitbar(ii/length(fns), msg) ;
    end
    
    mesh = ply_read(fullfile(fns(ii).folder, fns(ii).name));
    tri = cell2mat(mesh.face.vertex_indices) + 1;
    if strcmp(meshorder, 'zyx')
        xs = mesh.vertex.z / ssfactor ;
        ys = mesh.vertex.y / ssfactor ;
        zs = mesh.vertex.x / ssfactor ; 
        vn = [mesh.vertex.nz, mesh.vertex.ny, mesh.vertex.nx] ;
    else
        error('Did not code for this order yet')
    end
    
    vtx_sub = [xs, ys, zs] ;
    fv = struct('faces', tri, 'vertices', vtx_sub, 'normals', vn) ;
    
    % Check normals
    % close all
    % plot3(vtx_sub(1:10:end, 1), vtx_sub(1:10:end, 2), vtx_sub(1:10:end, 3), '.')
    % hold on
    % plot3(vtx_sub(1:10:end, 1) + 10 * vn(1:10:end, 1),...
    %     vtx_sub(1:10:end, 2) + 10 * vn(1:10:end, 2), ...
    %     vtx_sub(1:10:end, 3) + 10 * vn(1:10:end, 3), 'o')
    
    % View the normals a different way
    % close all
    % plot3(vtx_sub(1:10:end, 1), vtx_sub(1:10:end, 2), vtx_sub(1:10:end, 3), '.')
    % for i=1:10:length(vtx_sub)
    %     hold on
    %     plot3([vtx_sub(i, 1), vtx_sub(i, 1) + 10*vn(i, 1)], ... 
    %     [vtx_sub(i, 2), vtx_sub(i, 2) + 10*vn(i, 2)], ...
    %     [vtx_sub(i, 3), vtx_sub(i, 3) + 10*vn(i, 3)], 'r-') 
    % end
    % axis equal
    
    % Must either downsample mesh, compute xyzgrid using ssfactor and
    % pass to options struct.
    % Here, downsampled mesh
    % mesh.vertex.x = xs ;
    % mesh.vertex.y = ys ;
    % mesh.vertex.z = zs ;
    
    % Point match for aind and pind
    msg = strrep(['Point matching mesh ' fns(ii).name], '_', '\_') ;
    waitbar(ii/length(fns), fbar, msg)
    adist2 = sum((vtx_sub - acom) .^ 2, 2);
    %find the smallest distance and use that as an index 
    aind = find(adist2 == min(adist2)) ;
    % Next point match the posterior
    pdist2 = sum((vtx_sub - pcom) .^ 2, 2);
    %find the smallest distance and use that as an index
    pind = find(pdist2 == min(pdist2)) ;
    
    % Check it
    if preview
        trimesh(tri, vtx_sub(:, 1), vtx_sub(:, 2), vtx_sub(:, 3), vtx_sub(:, 1))
        hold on;
        plot3(vtx_sub(aind, 1), vtx_sub(aind, 2), vtx_sub(aind, 3), 'ko')
        plot3(vtx_sub(pind, 1), vtx_sub(pind, 2), vtx_sub(pind, 3), 'ro')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Grab dorsal direction if this is the first timepoint
    if ii == 1    
        disp('Obtaining dorsal direction since this is first TP...')
        if ~exist([rotname '.txt'], 'file') || overwrite_apdvcoms || ~exist(dcomname, 'file')
            apfn = fullfile(rootdir, ['Time_' timestr '_c1_stab_Probabilities_apcenterline.h5' ]);
            apdat = h5read(apfn, '/exported_data');
            ddat = permute(squeeze(apdat(dorsalChannel, :, :, :)), axorder) ;

            options.check = false ;
            dcom = com_region(ddat, dorsal_thres, options) ;
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
                    vtx_sub(:, 1), vtx_sub(:,2), vtx_sub(:, 3), ...
                    vtx_sub(:, 1)) ; % , 'edgecolor', 'none', 'FaceAlpha', 0.1) ;
                hold on;
                plot3(acom(1), acom(2), acom(3), 'ro')
                plot3(pcom(1), pcom(2), pcom(3), 'bo')
                plot4(dcom(1), dcom(2), dcom(3), 'go')
                xlabel('x [subsampled pixels]')
                ylabel('y [subsampled pixels]')
                zlabel('z [subsampled pixels]')
                title('Original mesh in subsampled pixels, with APD marked')
                axis equal
                %%%%%%%%%%%%%%%%%%%%%%
                waitfor(fig)
            end

            % Save dcom to dcomname
            header = 'Dorsal point in xyz subsampled pixel units' ;
            write_txt_with_header(dcomname, dcom, header)
            
            % compute rotation 
            apaxis = pcom - acom ;
            aphat = apaxis / norm(apaxis) ;

            % compute rotation matrix using this procedure: 
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
        else
            disp('Loading rotation from disk since already exists...')
            rot = dlmread([rotname '.txt']) ;
            dcom = dlmread(dcomname) ;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Define start point
    % Check if acom is inside mesh. If so, use that as starting point.
    ainside = inpolyhedron(fv, acom(1), acom(2), acom(3)) ;
    pinside = inpolyhedron(fv, pcom(1), pcom(2), pcom(3)) ;
    
    if ainside
        startpt = acom' ;
    else
        % move along the inward normal of the mesh from the matched vertex
        vtx = [vtx_sub(aind, 1), vtx_sub(aind, 2), vtx_sub(aind, 3)]' ;
        normal = fv.normals(aind, :) ;
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
    % Note: Keep startpt in subsampled units

    % Define end point
    if pinside
        endpt = pcom' ;
    else
        % move along the inward normal of the mesh from the matched vertex
        vtx = [vtx_sub(pind, 1), vtx_sub(pind, 2), vtx_sub(pind, 3)]' ;
        normal = fv.normals(pind, :) ;
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
    % Note: Keep endpt in subsampled units

    % Check out the mesh
    if preview
        hold on
        trimesh(fv.faces, xs, ys, zs)
        % plot3(xs, ys, zs, 'ko')
        scatter3(startpt(1), startpt(2), startpt(3), 'ro')
        scatter3(endpt(1), endpt(2), endpt(3), 'ko')
        xlabel('x [ssampled pixels]')
        ylabel('y [ssampled pixels]')
        zlabel('z [ssampled pixels]')
        hold off
        axis equal
    end

    %% Compute centerline if has not been saved 
    if overwrite || ~exist([outname '.txt'], 'file')
        % Declare in command output
        if exist([outname '.txt'], 'file')
            disp(['OVERWRITING: ' outname '.txt'])
        else
            disp(['Computing for ' outname '.txt'])
        end
        
        %% Get xyz grid for distance transform
        if ~exist('xx', 'var') || ~exist('yy', 'var') || ~exist('zz', 'var')
            xx = 0:res:ceil(max(xs) + buffer) ;
            yy = 0:res:ceil(max(ys) + buffer) ;
            zz = 0:res:ceil(max(zs) + buffer) ;
        end
        msg = strrep(['Identifying pts in mesh ' fns(ii).name], '_', '\_') ;
        waitbar(ii/length(fns), fbar, msg)
        tic 
        fv.faces = reorient_facets( fv.vertices, fv.faces );
        inside = inpolyhedron(fv, xx, yy, zz) ;
        outside = 1 - inside ;
        toc ; 

        % use the distanceTransform from Yuriy Mishchenko
        msg = strrep(['Computing DT for ' fns(ii).name], '_', '\_') ;
        waitbar(ii/length(fns), fbar, msg)
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
            close all ;
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
            % isonormals(x,y,z,v,p)
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
        msg = strrep(['Computing centerline for ' fns(ii).name], '_', '\_') ;
        waitbar(ii/length(fns), fbar, msg)
        tic

        % From example (DD is W, with low values being avoided)
        options.heuristic = weight * DD ;
        % Convert here to the gridspacing of xx,yy,zz
        startpt_transposed = [startpt(2), startpt(1), startpt(3)]' / res ;
        endpt_transposed = [endpt(2), endpt(1), endpt(3)]' / res ;
        [D2,S] = perform_fast_marching(DD, startpt_transposed, options);
        path = compute_geodesic(D2, endpt_transposed);
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

        % Convert skeleton's rows to columns and flip start/end
        skel_tmp = fliplr(path)' ;
        % Transpose x<->y back to original and scale to mesh units
        skel = [ skel_tmp(:,2), skel_tmp(:,1), skel_tmp(:,3) ] * res * ssfactor;
        
        % Save centerline as text file
        fid = fopen([outname '.txt'], 'wt');
        % make header
        fprintf(fid, 'centerline xyz in units of pixels (full resolution, same as mesh)');  
        fclose(fid);
        disp(['Saving centerline to txt: ', outname, '.txt'])
        dlmwrite([outname '.txt'], skel)
    else     
        skel = importdata([outname '.txt']) ;
    end
    
    %% Rescale start point and end point to full resolution
    spt = [startpt(1), startpt(2), startpt(3)] * ssfactor;
    ept = [endpt(1), endpt(2), endpt(3)] * ssfactor;
    
    %% Compute the translation to put anterior to origin
    if ii == 1
        % Save translation in units of mesh coordinates
        trans = -(rot * spt')' ;
        disp(['Saving translation vector (post rotation) to txt: ', transname, '.txt'])
        dlmwrite([transname '.txt'], trans)
    end
    
    %% Rotate and translate vertices and endpoints
    % Note: all in original mesh units (not subsampled)
    xyzr = (rot * vtx_sub')' * ssfactor + trans ;
    skelr = (rot * skel')' + trans ; 
    sptr = (rot * spt')' + trans ; 
    eptr = (rot * ept')' + trans ;
    dptr = (rot * (dcom' * ssfactor))' + trans ; 
    
    % Scale to actual resolution
    xyzrs = xyzr * resolution ;
    % get distance increment
    ds = vecnorm(diff(skel), 2, 2) ;
    % get pathlength at each skeleton point
    ss = [0; cumsum(ds)] ;
    sss = ss * resolution ;
    skelrs = skelr * resolution ;
    sptrs = sptr * resolution ;
    eptrs = eptr * resolution ; 
    dptrs = dptr * resolution ;
    
    % Save the rotated, translated, scaled curve
    disp(['Saving rotated & scaled skeleton to txt: ', skel_rs_outfn, '.txt'])
    fn = [skel_rs_outfn '.txt'] ;
    fid = fopen(fn, 'wt');
    % make header
    fprintf(fid, 'Aligned & scaled skeleton: sss [um], skelrs [um]');  
    fclose(fid);
    dlmwrite(fn, [sss, skelrs])
    
    %% Update our estimate for the true xyzlims
    xminrs = min(xminrs, min(xyzrs(:, 1))) ;
    yminrs = min(yminrs, min(xyzrs(:, 2))) ;
    zminrs = min(zminrs, min(xyzrs(:, 3))) ;
    xmaxrs = max(xmaxrs, max(xyzrs(:, 1))) ;
    ymaxrs = max(ymaxrs, max(xyzrs(:, 2))) ;
    zmaxrs = max(zmaxrs, max(xyzrs(:, 3))) ;
    
    %% Get axis limits if this is first TP
    if ii == 1
        % Check if already saved. If so, load it. Otherwise, guess.
        fntmp = [xyzlimname_um '.txt'] ;
        if exist(fntmp, 'file')
            xyzlims = dlmread(fntmp, ',', 1, 0) ;
            xminrs = xyzlims(1) ;
            yminrs = xyzlims(2) ;
            zminrs = xyzlims(3) ;
            xmaxrs = xyzlims(4) ;
            ymaxrs = xyzlims(5) ;
            zmaxrs = xyzlims(6) ;
            % Note that we can't simply rotate the bounding box, since it will
            % be tilted in the new frame. We must guess xyzlims for plotting
            % and update the actual xyzlims
            % this works for new box: resolution * ((rot * box')' + trans) ;
        end
        % Expand xyzlimits for plots
        xminrs_plot = xminrs - plot_buffer ;
        yminrs_plot = yminrs - plot_buffer ;
        zminrs_plot = zminrs - plot_buffer ;
        xmaxrs_plot = xmaxrs + plot_buffer ;
        ymaxrs_plot = ymaxrs + plot_buffer ;
        zmaxrs_plot = zmaxrs + plot_buffer ;    
    end
    
    %% Check the rotation
    if preview
        close all
        fig = figure ;
        tmp = trisurf(tri, xyzrs(:, 1), xyzrs(:,2), xyzrs(:, 3), ...
                    xyz(:, 1), 'edgecolor', 'none', 'FaceAlpha', 0.1) ;
        hold on;
        xyz = vtx_sub;
        tmp2 = trisurf(tri, xyz(:, 1), xyz(:,2), xyz(:, 3), ...
            xyz(:, 1), 'edgecolor', 'none', 'FaceAlpha', 0.1) ;
        boxx = [xmin, xmin, xmin, xmin, xmax, xmax, xmax, xmax, xmin] ;
        boxy = [ymin, ymax, ymax, ymin, ymin, ymax, ymax, ymin, ymin] ;
        boxz = [zmin, zmin, zmax, zmax, zmax, zmax, zmin, zmin, zmin] ;
        box = [boxx', boxy', boxz'] ;
        box_sub = box / ssfactor ; 
        boxrs = resolution * ((rot * box')' + trans) ;
        plot3(box_sub(:, 1), box_sub(:, 2), box_sub(:, 3), 'k-')
        plot3(boxrs(:, 1), boxrs(:, 2), boxrs(:, 3), 'k-')
        for i=1:3
            plot3([boxrs(i, 1), box_sub(i, 1)], ...
                [boxrs(i, 2), box_sub(i, 2)], ...
                [boxrs(i, 3), box_sub(i, 3)], '--')
        end
          
        % plot the skeleton
        % for i=1:length(skelrs)
        %     plot3(skelrs(:,1), skelrs(:,2), skelrs(:,3),'-','Color',[0,0,0], 'LineWidth', 3);
        % end
        plot3(sptrs(1), sptrs(2), sptrs(3), 'ro')
        plot3(eptrs(1), eptrs(2), eptrs(3), 'bo')
        plot3(dptrs(1), dptrs(2), dptrs(3), 'go')
        plot3(acom(1), acom(2), acom(3), 'rx')
        plot3(pcom(1), pcom(2), pcom(3), 'bx')
        plot3(dcom(1), dcom(2), dcom(3), 'gx')

        xlabel('x [$\mu$m or pix]', 'Interpreter', 'Latex'); 
        ylabel('y [$\mu$m or pix]', 'Interpreter', 'Latex');
        zlabel('z [$\mu$m or pix]', 'Interpreter', 'Latex');
        title('Checking rotation')
        axis equal
        waitfor(fig)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot and save
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Save plot of rotated and translated mesh
    if save_figs && (overwrite || overwrite_ims || ~exist([fig2outname '.png'], 'file'))
        disp('Saving rotated & translated figure (xy)...')    
        close all
        fig = figure ;
        set(gcf, 'Visible', 'Off')
        tmp = trisurf(tri, xyzrs(:, 1), xyzrs(:,2), xyzrs(:, 3), ...
            xyzrs(:, 1), 'edgecolor', 'none', 'FaceAlpha', 0.1) ;
        hold on;
        % plot the skeleton
        for i=1:length(skelrs)
            plot3(skelrs(:,1), skelrs(:,2), skelrs(:,3),'-','Color',[0,0,0], 'LineWidth', 3);
        end
        plot3(sptrs(1), sptrs(2), sptrs(3), 'ro')
        plot3(eptrs(1), eptrs(2), eptrs(3), 'bo')
        plot3(dptrs(1), dptrs(2), dptrs(3), 'go')
        xlabel('x [$\mu$m]', 'Interpreter', 'Latex'); 
        ylabel('y [$\mu$m]', 'Interpreter', 'Latex');
        zlabel('z [$\mu$m]', 'Interpreter', 'Latex');
        title(['Centerline using $D^{' num2str(exponent) '}$: ' timestr], ...
            'Interpreter', 'Latex')
        axis equal
        % xy
        view(2)
        xlim([xminrs_plot xmaxrs_plot]); 
        ylim([yminrs_plot ymaxrs_plot]); 
        zlim([zminrs_plot zmaxrs_plot]) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);
        saveas(fig, [fig1outname '.png'])
        
        % yz
        disp('Saving rotated & translated figure (yz)...')    
        view(90, 0);
        xlim([xminrs_plot xmaxrs_plot]); 
        ylim([yminrs_plot ymaxrs_plot]); 
        zlim([zminrs_plot zmaxrs_plot]) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=15cm
        saveas(fig, [fig2outname '.png'])
        % xz
        disp('Saving rotated & translated figure (xz)...')  
        view(0, 0)    
        xlim([xminrs_plot xmaxrs_plot]); 
        ylim([yminrs_plot ymaxrs_plot]); 
        zlim([zminrs_plot zmaxrs_plot]) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=15cm
        saveas(fig, [fig3outname '.png'])
        close all
    end
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    % xlabel('x'); ylabel('y'); zlabel('z'); axis equal
    % title(['$D^{' num2str(exponent) '}$'], 'Interpreter', 'Latex')
    % view(2)
    % xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax])
    % saveas(fig, [fig1outname '.png'])
    % view(10)
    % stopping 
    % saveas(fig, [fig2outname '.png'])
    % saveas(fig, [fig3outname '.png'])
    % close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Associate each vertex with a point on the curve
    % dist2 = (xs - skel(:,1)).^2 + (xs - skel(:,2)).^2 + (xs - skel(:,3)).^2 ;
    [kmatch, dist] = dsearchn(skel / ssfactor, vtx_sub) ;
    % A dsearchn() returns closest Euclidean 3D matches. Check that the
    % association linesegment between the vertex and the centerline does
    % not leave the mesh
        
    % Check the associations
    if preview
        clf ; hold on
        for ijk = 1:500:length(vtx_sub)
            plot3([vtx_sub(ijk, 1) skel(kmatch(ijk), 1)], ...
                [vtx_sub(ijk, 2) skel(kmatch(ijk), 2)], ...
                [vtx_sub(ijk, 3) skel(kmatch(ijk), 3)])
        end
        close all
    end
    
    % Compute radius R(s) in microns
    % radii = vecnorm(vtx_sub * ssfactor - skel(kmatch), 2, 2) * resolution ;
    
    % Compute phi(s), which is just the polar angle in the yz plane - pi/2
    % taken wrt the centerline
    phi_dorsal = mod(atan2(xyzr(:, 3) - sptr(3), xyzr(:, 2) - sptr(2)) - pi * 0.5, 2*pi);
    phi_ctrdorsal = mod(atan2(xyzr(:, 3) - skelr(kmatch, 3), ...
                        xyzr(:, 2) - skelr(kmatch, 2)) - pi * 0.5, 2*pi);

    % Save radii and angle wrt dorsal as text file
    % disp(['Saving radii to txt: ', polaroutfn, '.txt'])
    % fn = [polaroutfn '.txt'] ;
    % fid = fopen(fn, 'wt');
    % % make header
    % fprintf(fid, 'kmatch, radii [microns], phi_dorsal, phi_ctrdorsal');  
    % fclose(fid);
    % dlmwrite(fn, [kmatch, radii, phi_dorsal, phi_ctrdorsal])
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check phi_dorsal
    fexist1 = exist(fullfile(phi_def_outdir, [name '.png']), 'file') ;
    fexist2 = exist(fullfile(phicd_def_outdir, [name '.png']), 'file') ;
    figs_exist = fexist1 && fexist2 ;
    if save_figs && (overwrite || ~figs_exist)
        % view global dorsal angle
        fig = figure('Visible', 'Off');
        tmp = trisurf(tri, xyzrs(:, 1), xyzrs(:,2), xyzrs(:, 3), ...
            phi_dorsal / pi, 'edgecolor', 'none', 'FaceAlpha', 0.1) ;
        xlabel('x [$\mu$m]', 'Interpreter', 'Latex')
        ylabel('y [$\mu$m]', 'Interpreter', 'Latex')
        zlabel('z [$\mu$m]', 'Interpreter', 'Latex')
        xlim([xminrs_plot xmaxrs_plot])
        ylim([yminrs_plot ymaxrs_plot])
        zlim([zminrs_plot zmaxrs_plot])
        axis equal
        cb = colorbar() ;
        ylabel(cb, 'angle w.r.t. dorsal, $\phi / \pi$')
        cb.Label.Interpreter = 'latex';
        cb.Label.FontSize = 12 ;
        title('$\phi_{\textrm{dorsal}}$', 'Interpreter', 'Latex')
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=16cm
        saveas(fig, fullfile(phi_def_outdir, [name '.png']))
        close all
        
        % view dorsal angle from centerline
        fig = figure('Visible', 'Off');
        tmp = trisurf(tri, xyzrs(:, 1), xyzrs(:,2), xyzrs(:, 3), ...
            phi_ctrdorsal, 'edgecolor', 'none', 'FaceAlpha', 0.1) ;
        xlabel('x [$\mu$m]', 'Interpreter', 'Latex')
        ylabel('y [$\mu$m]', 'Interpreter', 'Latex')
        zlabel('z [$\mu$m]', 'Interpreter', 'Latex')
        xlim([xminrs_plot xmaxrs_plot])
        ylim([yminrs_plot ymaxrs_plot])
        zlim([zminrs_plot zmaxrs_plot])
        axis equal
        cb = colorbar() ;
        ylabel(cb, 'angle w.r.t. dorsal from centerline, $\phi / \pi$')
        cb.Label.Interpreter = 'latex';
        cb.Label.FontSize = 12 ;
        title('$\phi_{\textrm{dorsal}}^c$', 'Interpreter', 'Latex')
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=16cm
        saveas(fig, fullfile(phicd_def_outdir, [name '.png']))
        close all
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save the crude radius data as a plot
    % fexist1 = exist(fullfile(radius_vs_s_phi_outdir, [name '.png']), 'file') ;
    % fexist2 = exist(fullfile(radius_vs_s_phicd_outdir, [name '.png']), 'file') ;
    % figs_exist = fexist1 && fexist2 ;
    % if save_figs && (overwrite || figs_exist) 
    %     % Color by phi_dorsal
    %     close all
    %     fig = figure;
    %     set(gcf, 'Visible', 'Off')
    %     scatter(sss(kmatch), radii, [], ...
    %         phi_dorsal / pi, 'filled', ...
    %         'MarkerFaceAlpha', 0.05) ;        
    %     xlabel('pathlength, $s$ [$\mu$m]', 'Interpreter', 'Latex')
    %     ylabel('radius, $R$ [$\mu$m]', 'Interpreter', 'Latex')
    %     cb = colorbar() ;
    %     ylabel(cb, 'angle w.r.t. dorsal, $\phi / \pi$')
    %     cb.Label.Interpreter = 'latex';
    %     cb.Label.FontSize = 12 ;
    %     title('$\phi_{\textrm{dorsal}}$', 'Interpreter', 'Latex')
    %     xlim([0, 525])
    %     set(gcf, 'PaperUnits', 'centimeters');
    %     set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=16cm
    %     saveas(fig, fullfile(radius_vs_s_phi_outdir, [name '.png']))
    % 
    %     % Color by phi_ctrdorsal
    %     close all
    %     fig = figure;
    %     set(gcf, 'Visible', 'Off')
    %     scatter(ss(kmatch), radii, [], ...
    %         phi_dorsal / pi, 'filled', ...
    %         'MarkerFaceAlpha', 0.05) ;        
    %     xlabel('pathlength, $s$ [$\mu$m]', 'Interpreter', 'Latex')
    %     ylabel('radius, $R$ [$\mu$m]', 'Interpreter', 'Latex')
    %     cb = colorbar() ;
    %     ylabel(cb, 'angle w.r.t. dorsal, $\phi / \pi$')
    %     cb.Label.Interpreter = 'latex';
    %     cb.Label.FontSize = 12 ;
    %     title('Midgut radius')
    %     xlim([0, 525])
    %     set(gcf, 'PaperUnits', 'centimeters');
    %     set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=16cm
    %     saveas(fig, fullfile(radius_vs_s_phicd_outdir, [name '.png']))
    %     clf
    % end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save the rotated, translated, scaled to microns mesh ===============
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Shall we save the aligned mesh?')
    alignedmeshfn = fullfile(alignedmeshdir, [name '_APDV_um.ply']) ;
    if ~exist(alignedmeshfn, 'file') || overwrite
        disp('yes, Saving the aligned mesh...')
        disp([' --> ' alignedmeshfn])
        vtx_rs = (rot * (vtx_sub * ssfactor)' + trans')' * resolution ;
        vn_rs = (rot * fv.normals')' ;
        outfaces = [fv.faces(:, 2), fv.faces(:, 1), fv.faces(:, 3)] ;
        plywrite_with_normals(alignedmeshfn, outfaces, vtx_rs, vn_rs)
    else
        disp('no, already saved')
    end
       
    % Check the normals 
    if preview 
        close all
        plot3(vtx_rs(1:10:end, 1), vtx_rs(1:10:end, 2), vtx_rs(1:10:end, 3), '.')
        for i=1:10:length(vtx_rs)
            hold on
            plot3([vtx_rs(i, 1), vtx_rs(i, 1) + 10*vn_rs(i, 1)], ... 
            [vtx_rs(i, 2), vtx_rs(i, 2) + 10*vn_rs(i, 2)], ...
            [vtx_rs(i, 3), vtx_rs(i, 3) + 10*vn_rs(i, 3)], 'r-') 
        end
        axis equal
    end    
    
    % Save acom, pcom and their aligned counterparts as attributes in an
    % hdf5 file
    acom_rs = ((rot * acom' * ssfactor + trans') * resolution)' ;
    pcom_rs = ((rot * pcom' * ssfactor + trans') * resolution)' ; 
    dcom_rs = ((rot * dcom' * ssfactor + trans') * resolution)' ; 
    try
        h5create(outapdvname, ['/' name '/acom'], size(acom)) ;
    catch
        disp('acom already exists')
    end
    try
        h5create(outapdvname, ['/' name '/pcom'], size(pcom)) ;
    catch
        disp('pcom already exists')
    end
    try 
        h5create(outapdvname, ['/' name '/dcom'], size(dcom)) ;
    catch
        disp('dcom already exists')
    end
    try
        h5create(outapdvname, ['/' name '/acom_rs'], size(acom_rs)) ;
    catch
        disp('acom_rs already exists')
    end
    try
        h5create(outapdvname, ['/' name '/pcom_rs'], size(pcom_rs)) ;
    catch
        disp('pcom_rs already exists')
    end
    try 
        h5create(outapdvname, ['/' name '/dcom_rs'], size(dcom_rs)) ;
    catch
        disp('dcom_rs already exists')
    end
    
    h5write(outapdvname, ['/' name '/acom'], acom) ;
    h5write(outapdvname, ['/' name '/pcom'], pcom) ;
    h5write(outapdvname, ['/' name '/dcom'], dcom) ;
    h5write(outapdvname, ['/' name '/acom_rs'], acom_rs) ;
    h5write(outapdvname, ['/' name '/pcom_rs'], pcom_rs) ;
    h5write(outapdvname, ['/' name '/dcom_rs'], dcom_rs) ;
    % h5disp(outapdvname, ['/' name]);
    
end

% Save xyzlimits 
disp('Saving rot/trans mesh xyzlimits for plotting')
header = 'xyzlimits for rotated translated meshes in units of full resolution pixels' ;
fn = [xyzlimname '.txt'] ;
dat = [xminrs, xmaxrs; yminrs, ymaxrs; zminrs, zmaxrs] / resolution;
write_txt_with_header(fn, dat, header) ;

% Save xyzlimits in um
disp('Saving rot/trans mesh xyzlimits for plotting, in microns')
header = 'xyzlimits for rotated translated meshes in microns' ;
fn = [xyzlimname_um '.txt'] ;
dat = [xminrs, xmaxrs; yminrs, ymaxrs; zminrs, zmaxrs] ;
write_txt_with_header(fn, dat, header) ;

if isvalid(fbar)
    close(fbar)
end
disp('done')
