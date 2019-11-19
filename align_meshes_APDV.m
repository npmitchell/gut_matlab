%% Extract the centerlines from a series of meshes (PLY files) 
% Noah Mitchell 2019
% Saves scaled meshes with AP along x and DV along y, centered at A=0.
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
% Notes
% -----
% vertices are in units of pixels (at full resolution)
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
% To run after:
%    extract_centerline.m
%    slice_mesh_endcaps.m
%    extract_chirality_writhe.m
%    Generate_Axisymmetric_Pullbacks_Orbifold.m
%    
% NPMitchell 2019
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
overwrite = true ;  % recompute APDV rotation, translation
overwrite_apdvcoms = false ;  % recompute APDV coms from training
save_figs = true ;  % save images of cntrline, etc, along the way
overwrite_ims = true ;  % overwrite images even if centerlines are not overwritten
preview = false ;  % display intermediate results, for debugging
resolution = 0.2619 ;  % um per pixel for full resolution (not subsampled)
dorsal_thres = 0.5 ;  % threshold for extracting Dorsal probability cloud 
buffer = 5 ;  % extra space in meshgrid of centerline extraction, to ensure mesh contained in volume
plot_buffer = 40;  % extra space in plots, in um
ssfactor = 4;  % subsampling factor for the h5s used to train for mesh/acom/pcom/dcom
weight = 0.1;  % for speedup of centerline extraction. Larger is less precise
normal_step = 0.5 ;  % how far to move normally from ptmatched vtx if a/pcom is not inside mesh
eps = 0.01 ;  % value for DT outside of mesh in centerline extraction
meshorder = 'zyx' ;  % ordering of axes in loaded mesh wrt iLastik output
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
% Name output directory for apdv info
outdir = [fullfile(meshdir, 'centerline') filesep ];
if ~exist(outdir, 'dir')
    mkdir(outdir) ;
end
% Name the directory for outputting aligned_meshes
alignedmeshdir = fullfile(meshdir, ['aligned_meshes' filesep]) ;
if ~exist(alignedmeshdir, 'dir')
    mkdir(alignedmeshdir) ;
end
fns = dir(fullfile(meshdir, 'mesh_apical_stab_0*.ply')) ;
% Ensure that PLY files exist
if isempty(fns)
    error('Found no matching PLY files in ' + meshdir)
end

% Name the directory for outputting figures
figoutdir = [alignedmeshdir 'images' filesep];
if ~exist(figoutdir, 'dir')
    mkdir(figoutdir) ;
end
% figure 1
fig1outdir = [figoutdir 'aligned_mesh_xy' filesep];
if ~exist(fig1outdir, 'dir')
    mkdir(fig1outdir) ;
end
% figure 2
fig2outdir = [figoutdir 'aligned_mesh_xz' filesep];
if ~exist(fig2outdir, 'dir')
    mkdir(fig2outdir) ;
end
% figure 3
fig3outdir = [figoutdir 'aligned_mesh_yz' filesep];
if ~exist(fig3outdir, 'dir')
    mkdir(fig3outdir) ;
end

rotname = fullfile(meshdir, 'rotation_APDV') ;
transname = fullfile(meshdir, 'translation_APDV') ;
xyzlimname_raw = fullfile(meshdir, 'xyzlim') ;
xyzlimname = fullfile(meshdir, 'xyzlim_APDV') ;
xyzlimname_um = fullfile(meshdir, 'xyzlim_APDV_um') ;
outapdvname = fullfile(outdir, 'apdv_coms_rs.h5') ;

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
        error('breaking to avoid overwrite, debug')
        %% Get the timestamp string from the name of the mesh
        name_split = strsplit(fns(ii).name, '.ply') ;
        name = name_split{1} ; 
        tmp = strsplit(name, '_') ;
        timestr = tmp{length(tmp)} ;
        
        %% Load the AP axis determination
        msg = ['Computing acom, pcom for ' timestr ] ;
        disp(msg)
        if ~exist('fbar', 'var')
            fbar = waitbar(ii / length(fns), msg) ;
        else
            if ~isvalid(fbar)
              fbar = waitbar(ii / length(fns), msg) ;
            end
        end
        waitbar(ii / length(fns), fbar, msg)
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
    if preview
        acom_sm = h5read(rawapdvname, '/acom_sm');
        pcom_sm = h5read(rawapdvname, '/pcom_sm');
        plot3(acom_sm(:, 1), acom_sm(:, 2), acom_sm(:, 3))
        hold on;
        plot3(pcom_sm(:, 1), pcom_sm(:, 2), pcom_sm(:, 3))
        xlabel('x [subsampled pix]')
        ylabel('y [subsampled pix]')
        zlabel('z [subsampled pix]')
        axis equal
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get axis limits from looking at all meshes =============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
fn = [xyzlimname_raw '.txt'] ;
if exist(fn, 'file')
    disp('loading xyzlimits from disk')
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

%% With acom and pcom in hand, we compute dorsal and rot/trans ============
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
    tmp = strsplit(name, '_') ;
    timestr = tmp{length(tmp)} ;
    
    fig1outname = [fullfile(fig1outdir, name) '_xy'] ;
    fig2outname = [fullfile(fig2outdir, name) '_xz'] ;
    fig3outname = [fullfile(fig3outdir, name) '_yz'] ;

    
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
    try
        waitbar(ii/length(fns), fbar, msg)
    catch
        disp(msg)
    end
    adist2 = sum((vtx_sub - acom) .^ 2, 2);
    %find the smallest distance and use that as an index 
    aind = find(adist2 == min(adist2)) ;
    % Next point match the posterior
    pdist2 = sum((vtx_sub - pcom) .^ 2, 2);
    %find the smallest distance and use that as an index
    pind = find(pdist2 == min(pdist2)) ;
    
    % Check it
    if preview
        disp('Previewing mesh in figure window')
        trimesh(tri, vtx_sub(:, 1), vtx_sub(:, 2), vtx_sub(:, 3), vtx_sub(:, 1))
        hold on;
        plot3(vtx_sub(aind, 1), vtx_sub(aind, 2), vtx_sub(aind, 3), 'ko')
        plot3(vtx_sub(pind, 1), vtx_sub(pind, 2), vtx_sub(pind, 3), 'ro')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Grab dorsal direction if this is the first timepoint
    if ii == 1  
        disp('Obtaining dorsal direction since this is first TP')
        if ~exist([rotname '.txt'], 'file') || overwrite
            apfn = fullfile(rootdir, ['Time_' timestr '_c1_stab_Probabilities_apcenterline.h5' ]);
            apdat = h5read(apfn, '/exported_data');
            ddat = permute(squeeze(apdat(dorsalChannel, :, :, :)), axorder) ;

            options.check = preview ;
            options.check_slices = false ;
            search4com = true ;
            % start with a threshold == dorsal_thres, iteratively lower if
            % necessary
            tmp_dorsal_thres = dorsal_thres ;
            while search4com 
                try
                    msg = 'Finding com region of dorsal data: ' ;
                    disp([msg 'thres=' num2str(tmp_dorsal_thres)])
                    dcom = com_region(ddat, tmp_dorsal_thres, options) ;
                    search4com = false ;
                catch
                    disp('no region found, lowering dorsal threshold for prob cloud') ;
                    tmp_dorsal_thres = 0.9 * tmp_dorsal_thres ;
                end
            end
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
                plot3(dcom(1), dcom(2), dcom(3), 'go')
                xlabel('x [subsampled pixels]')
                ylabel('y [subsampled pixels]')
                zlabel('z [subsampled pixels]')
                title('Original mesh in subsampled pixels, with APD marked')
                axis equal
                %%%%%%%%%%%%%%%%%%%%%%
                waitfor(fig)
            end

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
            rot = dlmread([rotname '.txt']) ;
        end
    end
    
    if ii == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Define start point and endpoint for first TP
        disp('Defining start point and endpoint for first TP')
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
        
        %% Rescale start point and end point to full resolution
        spt = [startpt(1), startpt(2), startpt(3)] * ssfactor;
        ept = [endpt(1), endpt(2), endpt(3)] * ssfactor;
    end
    
    
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
    sptr = (rot * spt')' + trans ; 
    eptr = (rot * ept')' + trans ;
    dptr = (rot * (dcom' * ssfactor))' + trans ; 
    
    % Scale to actual resolution
    xyzrs = xyzr * resolution ;
    sptrs = sptr * resolution ;
    eptrs = eptr * resolution ; 
    dptrs = dptr * resolution ;
    
    %% Update our estimate for the true xyzlims
    xminrs = min(xminrs, min(xyzrs(:, 1))) ;
    yminrs = min(yminrs, min(xyzrs(:, 2))) ;
    zminrs = min(zminrs, min(xyzrs(:, 3))) ;
    xmaxrs = max(xmaxrs, max(xyzrs(:, 1))) ;
    ymaxrs = max(ymaxrs, max(xyzrs(:, 2))) ;
    zmaxrs = max(zmaxrs, max(xyzrs(:, 3))) ;
    
    %% Get a guess for the axis limits if this is first TP
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
        plot3(sptrs(1), sptrs(2), sptrs(3), 'ro')
        plot3(eptrs(1), eptrs(2), eptrs(3), 'bo')
        plot3(dptrs(1), dptrs(2), dptrs(3), 'go')
        xlabel('x [$\mu$m]', 'Interpreter', 'Latex'); 
        ylabel('y [$\mu$m]', 'Interpreter', 'Latex');
        zlabel('z [$\mu$m]', 'Interpreter', 'Latex');
        title(['Aligned mesh: ' timestr], 'Interpreter', 'Latex')
        axis equal
        
        % xy
        view(2)
        xlim([xminrs_plot xmaxrs_plot]); 
        ylim([yminrs_plot ymaxrs_plot]); 
        zlim([zminrs_plot zmaxrs_plot]) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);
        disp(['Saving to ' fig1outname '.png'])
        saveas(fig, [fig1outname '.png'])
        
        % yz
        disp('Saving rotated & translated figure (yz)...')    
        view(90, 0);
        xlim([xminrs_plot xmaxrs_plot]); 
        ylim([yminrs_plot ymaxrs_plot]); 
        zlim([zminrs_plot zmaxrs_plot]) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);  % x_width=10cm y_width=15cm
        saveas(fig, [fig2outname '.png'])
        % xz
        disp('Saving rotated & translated figure (xz)...')  
        view(0, 0)    
        xlim([xminrs_plot xmaxrs_plot]); 
        ylim([yminrs_plot ymaxrs_plot]); 
        zlim([zminrs_plot zmaxrs_plot]) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);  % x_width=10cm y_width=15cm
        saveas(fig, [fig3outname '.png'])
        close all
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save the rotated, translated, scaled to microns mesh ===============
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Saving the aligned mesh...')
    alignedmeshfn = fullfile(alignedmeshdir, [name '_APDV_um.ply']) ;
    disp([' --> ' alignedmeshfn])
    vtx_rs = (rot * (vtx_sub * ssfactor)' + trans')' * resolution ;
    vn_rs = (rot * fv.normals')' ;
    outfaces = [fv.faces(:, 2), fv.faces(:, 1), fv.faces(:, 3)] ;
    plywrite_with_normals(alignedmeshfn, outfaces, vtx_rs, vn_rs)
       
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
