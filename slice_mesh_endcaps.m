%% Create cut meshes with endcaps removed 
% Note that a later step involves cutMesh, which slices along AP.
% Noah Mitchell 2019
% This version relies on Gabriel Peyre's toolbox called
% toolbox_fast_marching/
%
% Prerequisites
% -------------
% Gut_Pipeline.m
% align_meshes.m (after training on APD)
% extract_centerline.m (after training on APD)
% 
% Postrequisites (codes to run after)
% -----------------------------------
% extract_chirality_writhe.m
% Generate_Axisymmetric_Pullbacks_Orbifold.m (script)
%
% Run from the msls_output directory
% First run extract_centerline.m before running this code
% Run this code only after training on anterior (A), posterior (P), and 
% dorsal anterior (D) points in different iLastik channels.
% anteriorChannel, posteriorChannel, and dorsalChannel specify the iLastik
% training channel that is used for each specification.
% Name the h5 file output from iLastik as ..._Probabilities_apcenterline.h5
% Train for anterior dorsal (D) only at the first time point, because
% that's the only one that's used.
clear ; close all
% cd /mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/

%% First, compile required c code
% mex ./FastMarching_version3b/shortestpath/rk4
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
cd(odir)

%% Parameters
adist_thres = 40 ;  % distance threshold for cutting off anterior in pix
pdist_thres = 30 ;  % distance threshold for cutting off posterior in pix
overwrite = false ;  % recompute centerline
save_figs = true ;  % save images of cntrline, etc, along the way
preview = false ;  % display intermediate results
res = 1 ;  % pixels per gridspacing of DT for cntrline extraction
resolution = 0.2619 ;  % um per pixel for full resolution (not subsampled)
dorsal_thres = 0.9 ;  % threshold for extracting Dorsal probability cloud 
buffer = 5 ;  % extra space in meshgrid of centerline extraction, to ensure mesh contained in volume
plot_buffer = 15 ;  % extra space around mesh bounding box in 3D plots
ssfactor = 4;  % subsampling factor for the h5s used to train for mesh/acom/pcom/dcom
weight = 0.1;  % for speedup of centerline extraction. Larger is less precise
normal_step = 0.5 ;  % how far to move normally from ptmatched vtx if a/pcom is not inside mesh
eps = 0.01 ;  % value for DT outside of mesh in centerline extraction
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

% Name output directory
outdir = [fullfile(meshdir, 'cylindercut') filesep ];
if ~exist(outdir, 'dir')
    mkdir(outdir) ;
end
figoutdir = [outdir 'images' filesep];
if ~exist(figoutdir, 'dir')
    mkdir(figoutdir) ;
end
fns = dir(fullfile(meshdir, 'mesh_apical_stab_0*.ply')) ;
rotname = fullfile(meshdir, 'rotation_APDV') ;
transname = fullfile(meshdir, 'translation_APDV') ;
xyzlimname = fullfile(meshdir, 'xyzlim_APDV_um') ;
apdvfn = fullfile(meshdir, ['centerline' filesep 'apdv_coms_from_training.h5']) ;
cntrlinedir = fullfile(meshdir, 'centerline') ;
outapd_boundaryfn = fullfile(outdir, 'ap_boundary_dorsalpts.h5') ;
boundaryfn = fullfile(outdir, 'ap_boundary_indices.h5') ;
outapphicd_fn = fullfile(outdir, 'ap_boundary_phicd_values.h5') ;
keepfn = fullfile(outdir, 'cut_endcap_indx.h5') ;

ii = 1 ;

%% Load scaling and limits
% Load the rotation matrix
rot = importdata([rotname '.txt']) ;
% Load the translation to put anterior to origin
trans = importdata([transname '.txt']) ;
% Load plotting limits
xyzlim = dlmread([xyzlimname '.txt'], ',', 1, 0) ;
xmin = xyzlim(1) - plot_buffer;
ymin = xyzlim(2) - plot_buffer;
zmin = xyzlim(3) - plot_buffer;
xmax = xyzlim(4) + plot_buffer;
ymax = xyzlim(5) + plot_buffer;
zmax = xyzlim(6) + plot_buffer;

%% Load AP coms
acom_sm = h5read(apdvfn, '/acom') ;
pcom_sm = h5read(apdvfn, '/pcom') ;

fbar = waitbar(0, 'Initializing...') ;
%% Iterate through each mesh
for ii=1:length(fns)  
    acom = acom_sm(ii, :) ;
    pcom = pcom_sm(ii, :) ;
    
    %% Name the output mesh filename
    name_split = strsplit(fns(ii).name, '.ply') ;
    name = name_split{1} ; 
    expstr = strrep(num2str(exponent, '%0.1f'), '.', 'p') ;
    resstr = strrep(num2str(res, '%0.1f'), '.', 'p') ;
    extenstr = ['_exp' expstr '_res' resstr] ;
    cntrlinefn = [fullfile(cntrlinedir, name) '_centerline' extenstr] ;
    polaroutfn = [fullfile(cntrlinedir, name) '_polarcoords' extenstr] ;
    skel_rs_outfn = [fullfile(cntrlinedir, name) '_centerline_scaled' extenstr ] ;
    tmp = strsplit(name, '_') ;
    timestr = tmp{length(tmp)} ;
    outfn = fullfile(outdir, [name '_cylindercut.ply']);
    
    %% Read the mesh
    msg = strrep(['Loading mesh ' fns(ii).name], '_', '\_') ;
    waitbar(ii/length(fns), fbar, msg) ;
    % if exist('fbar', 'var')
    %     if isvalid(fbar)    
    %         waitbar(ii/length(fns), fbar, msg)
    %     else
    %         fbar = waitbar(ii/length(fns), msg) ;
    %     end
    % else
    %     fbar = waitbar(ii / length(fns), msg) ;
    % end
    
    % Compute the endcaps if not already saved
    if overwrite || ~exist(outfn, 'file') || ~exist(keepfn, 'file')
        mesh = read_ply_mod(fullfile(fns(ii).folder, fns(ii).name));
        tri = mesh.f ;
        if strcmp(meshorder, 'zyx')
            xs = mesh.v(:, 3) / ssfactor ;
            ys = mesh.v(:, 2) / ssfactor ;
            zs = mesh.v(:, 1) / ssfactor ; 
            vn = [mesh.vn(:, 3), mesh.vn(:, 2), mesh.vn(:, 1)] ;
        else
            error('Did not code for this order yet')
        end
        vtx = [xs, ys, zs] ;
        fv = struct('f', tri, 'v', vtx, 'vn', vn) ;
        disp('loaded mesh.')

        %% Remove anterior endcap    
        % Measure distance to the posterior
        adist2 = sum((vtx - acom) .^ 2, 2);

        % Strategy: remove within distance of acom
        pts_to_remove = find(adist2 < adist_thres^2) ;

        % Make sure that we are removing a connected component
        % form a mesh from the piece(s) to be removed
        allpts = linspace(1, length(vtx), length(vtx)) ;
        all_but_acut = uint16(setdiff(allpts', pts_to_remove)) ;
        [ acutfaces, acutvtx, ~] = remove_vertex_from_mesh( fv.f, fv.v, all_but_acut ) ;
        [ acutface, acutvtx, connected_indices, npieces ] = remove_isolated_mesh_components( acutfaces, acutvtx ) ;
        if npieces > 1
            disp('Ensuring that only a single component is removed')
            pts_to_remove = pts_to_remove(connected_indices) ;
        end

        % Remove it
        [faces, vtx, keep_acut] = remove_vertex_from_mesh(fv.f, fv.v, pts_to_remove ) ;
        vn = fv.vn(keep_acut, :) ;
        disp(['Removed ' num2str(length(pts_to_remove)) ' vertices with acut'])
        
        % Inspect
        if preview
            fig = figure;
            trimesh(faces, vtx(:, 1), vtx(:, 2), vtx(:, 3))
            view(30,145)
            waitfor(fig)
        end

        %% Remove the posterior part as well
        % Measure distance to the posterior
        pdist2 = sum((vtx - pcom) .^ 2, 2);

        % STRATEGY 1: within distance of pcom
        pdist_thres_ii = pdist_thres ;
        pcut_done = false ;
        while ~pcut_done
            pts_to_remove = find(pdist2 < pdist_thres^2) ;

            % Make sure that we are removing a connected component
            % form a mesh from the piece(s) to be removed
            allpts = linspace(1, length(vtx), length(vtx)) ;
            all_but_pcut = uint16(setdiff(allpts, pts_to_remove)) ;
            [ pcutfaces, pcutvtx, ~] = remove_vertex_from_mesh( faces, vtx, all_but_pcut ) ;
            [ pcutface, pcutvtx, connected_indices, npieces ] = remove_isolated_mesh_components( pcutfaces, pcutvtx ) ;
            % If there were more than one piece selected, remove the bigger one only
            if npieces > 1
                pts_to_remove = pts_to_remove(connected_indices) ;
            end

            % check it
            if preview
                fig = figure ;
                scatter3(vtx(:, 1), vtx(:, 2), vtx(:, 3))
                hold on;
                scatter3(vtx(pts_to_remove, 1), vtx(pts_to_remove, 2), vtx(pts_to_remove, 3),  'filled')
                plot3(pcom(1), pcom(2), pcom(3), 's')
                waitfor(fig)
            end

            % Remove the posterior cap
            [faces_postpcut, vtx_postpcut, keep_pcut] = remove_vertex_from_mesh(faces, vtx, pts_to_remove ) ;

            disp(['Removed ' num2str(length(pts_to_remove)) ' vertices with pcut'])
            nremain = length(keep_pcut) ;
            nbefore = length(keep_acut) ;
            assert(nbefore - nremain == length(pts_to_remove))
            
            % Repeat if no vertices are removed
            if length(pts_to_remove) > 0
                pcut_done = true ;
            else
                pdist_thres_ii = pdist_thres * 1.1 ;
            end
        end
        faces = faces_postpcut ;
        vtx = vtx_postpcut ;           
        vn = vn(keep_pcut, :) ;
        %% figure out which indices were kept
        keep = keep_acut(keep_pcut) ;
        
        %% Check that the remainder is a single connected component
        % currently have faces, vtx. 
        [ faces, vtx, keep_final_pass, npieces ] = remove_isolated_mesh_components( faces, vtx ) ;
        if any(npieces > 1)
            keep = keep(keep_final_pass) ;
            vn = vn(keep_final_pass, :) ;
        end
        nremain = length(keep_final_pass) ;
        nbefore = length(keep) ;
        disp(['Removed ' num2str(nbefore - nremain) ' vertices with final pass'])

        %% Save the data in units of pixels (same as original mesh)
        faces = [faces(:, 2), faces(:, 1), faces(:, 3) ];
        plywrite_with_normals(outfn, faces, vtx * ssfactor, vn)
        
        % Save the indices to keep when cutting off endcaps
        try
            h5create(keepfn, ['/' name], size(keep)) ;
        catch
            msg = 'keep already exists in keepfn' ;
            waitbar(ii / length(fns), fbar, msg) ;
        end
        h5write(keepfn, ['/' name], keep)
    else
        meshcut = read_ply_mod(outfn) ;
        faces = meshcut.f ;
        vtx = meshcut.v / ssfactor ; 
        vn = meshcut.vn ;
        
        % Load which indices to keep when cutting off endcaps
        keep = h5read(keepfn, ['/' name]) ;
    end
    % % Check the normals 
    % close all
    % plot3(vtx(1:10:end, 1), vtx(1:10:end, 2), vtx(1:10:end, 3), '.')
    % for i=1:10:length(vtx)
    %     hold on
    %     plot3([vtx(i, 1), vtx(i, 1) + 10*vn(i, 1)], ... 
    %     [vtx(i, 2), vtx(i, 2) + 10*vn(i, 2)], ...
    %     [vtx(i, 3), vtx(i, 3) + 10*vn(i, 3)], 'r-') 
    % end
    % axis equal

    %% Load phi
    fn = [polaroutfn '.txt'] ;
    dat = dlmread(fn, ',') ;
    phicd = dat(:, 4) ;
    phicd_kept = phicd(keep) ;
        
    %% Compute dorsal vertex on anterior and posterior free boundaries
    TR = triangulation(faces, vtx) ;
    boundary = freeBoundary(TR) ; 
    nn = length(boundary(:)) ;
    bb = zeros(nn, 1) ;
    bb(1) = boundary(1) ;
    dmyk = 2;
    for kk = 1:length(boundary)
        seg = boundary(kk, :) ;
        if ~any(bb == seg(1))
            bb(dmyk) = seg(1) ;
            dmyk = dmyk + 1;
        end
        if ~any(bb == seg(2))
            bb(dmyk) = seg(2) ;
            dmyk = dmyk + 1 ;
        end
    end
    bb = bb(bb > 0) ;
    
    % Check the boundary segments
    % for kk = 1:length(boundary)
    %     row = boundary(kk, :) ;
    %     plot3(vtx(row, 1), vtx(row, 2), vtx(row, 3), '-')
    %     hold on
    % end
    % Check the slimmed boundary
    % plot3(vtx(bb, 1), vtx(bb, 2), vtx(bb, 3), '-')

    % figure out if boundary is anterior or posterior
    % Note that this assumes an elongated structure
    adist = sum((vtx(bb, :) - acom) .^ 2, 2) ;
    pdist = sum((vtx(bb, :) - pcom) .^ 2, 2) ;
    ab = bb(adist < pdist) ;
    pb = bb(adist > pdist) ;
    % check them
    if preview
        fig = figure;
        trisurf(faces, vtx(:, 1), vtx(:, 2), vtx(:, 3), 'edgecolor', 'none')
        hold on;
        plot3(vtx(bb, 1), vtx(bb, 2), vtx(bb, 3), 'b.-')
        plot3(acom(1), acom(2), acom(3), 'o')
        plot3(pcom(1), pcom(2), pcom(3), 'o')
        waitfor(fig)
    end
    
    % Save the anterior and posterior boundary indices (after endcap cut)
    try
        h5create(boundaryfn, ['/' name '/anterior_boundary_indices'], size(ab)) ;
    catch
        msg = 'ab already exists in boundaryfn' ;
        waitbar(ii / length(fns), fbar, msg) ;
    end
    try
        h5create(boundaryfn, ['/' name '/posterior_boundary_indices'], size(pb)) ;
    catch
        msg = 'pb already exists in boundaryfn' ;
        waitbar(ii / length(fns), fbar, msg) ;
    end
    h5write(boundaryfn, ['/' name '/anterior_boundary_indices'], ab)
    h5write(boundaryfn, ['/' name '/posterior_boundary_indices'], pb)
    
    % choose anterior Dorsal as point with smallest phi
    % (closest to zero or 2pi) if this is first TP. Otherwise, match prev.
    aphicd = phicd_kept(ab) ;
    pphicd = phicd_kept(pb) ;
    if ii == 1
        a_phipi = aphicd ;
        p_phipi = pphicd ;
        a_phipi(a_phipi > pi) = a_phipi(a_phipi > pi) - 2 * pi ;
        p_phipi(p_phipi > pi) = p_phipi(p_phipi > pi) - 2 * pi ;
        [~, adb_tmp] = min(abs(a_phipi)) ;
        [~, pdb_tmp] = min(abs(p_phipi)) ;
        % anterior/posterior dorsal boundary
        adb = ab(adb_tmp) ;
        pdb = pb(pdb_tmp) ;
    else
        % Match previous timepoint
        ka = dsearchn(vtx(ab, :), previous_avtx) ;
        kp = dsearchn(vtx(pb, :), previous_pvtx) ;
        adb = ab(ka) ;
        pdb = pb(kp) ;
    end
    previous_avtx = vtx(adb, :) ;
    previous_pvtx = vtx(pdb, :) ;
    % check 
    % hold on
    % scatter3(vtx(boundary, 1), vtx(boundary, 2), vtx(boundary, 3), [], phicd_kept(boundary), 'filled')
    % plot3(vtx(adb, 1), vtx(adb, 2), vtx(adb, 3), 'o')
    % plot3(vtx(pdb, 1), vtx(pdb, 2), vtx(pdb, 3), 'o')
    
    %% Save the dorsal angles of each endcap
    try
        h5create(outapphicd_fn, ['/' name '/aphicd'], size(aphicd))
    catch 
        msg = 'aphicd already exists' ;
        waitbar(ii / length(fns), fbar, msg) ;
    end
    try
        h5create(outapphicd_fn, ['/' name '/pphicd'], size(pphicd))
    catch 
        msg = 'pphicd already exists' ;
        waitbar(ii / length(fns), fbar, msg) ;
    end
    h5write(outapphicd_fn, ['/' name '/aphicd'], aphicd) ;
    h5write(outapphicd_fn, ['/' name '/pphicd'], pphicd) ;
    
    %% Save the anterior dorsal and posterior dorsal vertex points
    try 
        h5create(outapd_boundaryfn, ['/' name '/adorsal'], size(adb)) ;
    catch
        msg = 'adorsal pt already exists' ;
        waitbar(ii / length(fns), fbar, msg) ;
    end
    try
        h5create(outapd_boundaryfn, ['/' name '/pdorsal'], size(pdb)) ;
    catch
        msg = 'pdorsal pt already exists' ;
        waitbar(ii / length(fns), fbar, msg) ;
    end
    h5write(outapd_boundaryfn, ['/' name '/adorsal'], adb) ;
    h5write(outapd_boundaryfn, ['/' name '/pdorsal'], pdb) ;
    
    %% Plot the result
    figfn = fullfile(figoutdir, [name '.png']) ;
    if preview || (save_figs && ~exist(figfn, 'file'))
        disp(['Saving figure for frame ' num2str(ii)])
        fig = figure('Visible', 'Off');
        vtxrs = (rot * vtx' * ssfactor + trans')' * resolution ;
        % trimesh(faces, vtxr(:, 1), vtxr(:, 2), vtxr(:, 3)) ;
        tmp = trisurf(faces, vtxrs(:, 1), vtxrs(:, 2), vtxrs(:, 3), ...
            phicd_kept, 'edgecolor', 'none', 'FaceAlpha', 0.6) ;
        hold on
        plot3(vtxrs(adb, 1), vtxrs(adb, 2), vtxrs(adb, 3), 'o')
        plot3(vtxrs(pdb, 1), vtxrs(pdb, 2), vtxrs(pdb, 3), 'o')
        axis equal
        xlim([xmin xmax])
        ylim([ymin ymax])
        zlim([zmin zmax])
        xlabel('x [$\mu$m]', 'Interpreter', 'Latex')
        ylabel('y [$\mu$m]', 'Interpreter', 'Latex')
        zlabel('z [$\mu$m]', 'Interpreter', 'Latex')
        title(['Mesh with cylindrical cuts, $t=$' timestr], ...
                'Interpreter', 'Latex')
        %view(50, 145)
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=16cm
        saveas(fig, figfn)

        if preview
            set(fig, 'Visible', 'On') ;
            waitfor(fig) ;
        else
            close(fig)
        end
    end
    
    
end

close(fbar)
disp('done')
