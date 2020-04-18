function sliceMeshEndcaps(QS, opts, methodOpts)
% SLICEMESHENDCAPS()
% Create cut meshes with endcaps removed 
% Note that a later step involves cutMesh, which slices along AP.
% Noah Mitchell 2019
% This version relies on Gabriel Peyre's toolbox called
% toolbox_fast_marching/
%
% Parameters
% ----------
% resolution : float
%   um per pixel for full resolution (not subsampled)
%
% Prerequisites
% -------------
% align_meshes.m or alignMeshesAPDV() after training on APD
% extract_centerline.m or extractCenterlineSeries() after training on APD
% 
% Postrequisites (codes to run after)
% -----------------------------------
% extract_chirality_writhe.m or extractChiralityWrithe
% Generate_Axisymmetric_Pullbacks_Orbifold.m (script) or continue master
%   pipeline script
%
% First run extract_centerline.m or extractCenterlineSeries() before 
% running this code.
% Run this code only after training on anterior (A), posterior (P), and 
% dorsal anterior (D) points in different iLastik channels.
% anteriorChannel, posteriorChannel, and dorsalChannel specify the iLastik
% training channel that is used for each specification.
% Name the h5 file output from iLastik as ..._Probabilities_apcenterline.h5
% Train for anterior dorsal (D) only at the first time point, because
% that's the only one that's used.

%% Parameters from QS
timePoints = QS.xp.fileMeta.timePoints ;

%% Parameters
adist_thres = 20 ;  % distance threshold for cutting off anterior in pix
pdist_thres = 15 ;  % distance threshold for cutting off posterior in pix
overwrite = false ;  % recompute sliced endcaps
save_figs = true ;  % save images of cntrline, etc, along the way
preview = false ;  % display intermediate results
if isfield(opts, 'adist_thres')
    adist_thres = opts.adist_thres ;
end
if isfield(opts, 'pdist_thres')
    pdist_thres = opts.pdist_thres ;
end
if isfield(methodOpts, 'overwrite')
    overwrite = methodOpts.overwrite ;
end
if isfield(methodOpts, 'save_figs')
    save_figs = methodOpts.save_figs ;
end
if isfield(methodOpts, 'preview')
    preview = methodOpts.preview;
end

% figure parameters
xwidth = 16 ; % cm
ywidth = 10 ; % cm

% subsampling factor for the h5s used to train for mesh/acom/pcom/dcom
ssfactor = QS.xp.detector.options.ssfactor ; 
[~, ~, xyzlim_um] = QS.getXYZLims() ;

% Name output directory
outdir = QS.dir.cylinderMesh ;
figoutdir = fullfile(outdir, 'images') ;
if ~exist(figoutdir, 'dir')
    mkdir(figoutdir) ;
end

%% Load AP coms
[acom_sm, pcom_sm] = QS.getAPCOMSm ;

%% Iterate through each mesh
for ii=1:length(timePoints)
    tt = timePoints(ii) ;
    disp(['tt = ' num2str(tt)])
    
    acom = acom_sm(ii, :) ;
    pcom = pcom_sm(ii, :) ;
    
    %% Name the output mesh filename
    name = sprintf(QS.fileBase.mesh, tt) ;
    meshfn = sprintf(QS.fullFileBase.mesh, tt) ;
    outfn = sprintf(QS.fullFileBase.cylinderMesh, tt) ;
    keepfn = sprintf(QS.fullFileBase.cylinderKeep, tt) ; 
    boundaryfn = sprintf(QS.fullFileBase.apBoundary, tt) ;
    outapd_boundaryfn = sprintf(QS.fullFileBase.apBoundaryDorsalPts, tt) ;
    
    
    %% Read the mesh
    disp(['Loading mesh ' name]) ;
    
    % Compute the endcaps if not already saved
    if overwrite || ~exist(outfn, 'file') || ~exist(keepfn, 'file')
        if exist(outfn, 'file')
            disp(['Overwriting cylinder mesh: ' outfn])
        end
        if exist(keepfn, 'file')
            disp(['Overwriting keep indices: ' keepfn])
        end
        disp(['Computing endcaps for ' name])
        error('dfdd')
        mesh = read_ply_mod(sprintf(meshfn, tt));
        % subsample the mesh to match acom, pcom
        vtx = mesh.v / ssfactor ;
        fv = struct('f', mesh.f, 'v', vtx, 'vn', mesh.vn) ;
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
        [ acutface, acutvtx, connected_indices, npieces ] = ...
            remove_isolated_mesh_components( acutfaces, acutvtx ) ;
        if npieces > 1
            disp('Ensuring that only a single component is removed')
            pts_to_remove = pts_to_remove(connected_indices) ;
        end

        % Remove it
        [faces, vtx, keep_acut] = remove_vertex_from_mesh(fv.f, fv.v, pts_to_remove ) ;
        vn = fv.vn(keep_acut, :) ;
        disp(['Removed ' num2str(length(pts_to_remove)) ' vertices with acut'])
        
        % Inspect
        % if preview
        %     fig = figure;
        %     trimesh(faces, vtx(:, 1), vtx(:, 2), vtx(:, 3))
        %     view(30,145)
        %     waitfor(fig)
        % end

        %% Remove the posterior part as well
        % Measure distance to the posterior
        pdist2 = sum((vtx - pcom) .^ 2, 2);

        % STRATEGY 1: within distance of pcom
        pdist_thres_ii = pdist_thres ;
        pcut_done = false ;
        while ~pcut_done
            disp(['finding points within ' num2str(pdist_thres) ' of pcom'])
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
            
            % % check it
            % if preview
            %     fig = figure ;
            %     scatter3(vtx(:, 1), vtx(:, 2), vtx(:, 3))
            %     hold on;
            %     scatter3(vtx(pts_to_remove, 1), vtx(pts_to_remove, 2), vtx(pts_to_remove, 3),  'filled')
            %     plot3(pcom(1), pcom(2), pcom(3), 's')
            %     waitfor(fig)
            % end

            % Remove the posterior cap
            [faces_postpcut, vtx_postpcut, keep_pcut] = ...
                remove_vertex_from_mesh(faces, vtx, pts_to_remove ) ;

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
        nremain = length(keep_final_pass) ;
        nbefore = length(keep) ;
        disp(['Removed ' num2str(nbefore - nremain) ' vertices with final pass'])
        
        if any(npieces > 1)
            % Update keep here to reflect final pass
            keep = keep(keep_final_pass) ;
            vn = vn(keep_final_pass, :) ;
        end

        %% Save the data in units of pixels (same as original mesh)
        disp(['Saving cylinder mesh to ' outfn])
        plywrite_with_normals(outfn, faces, vtx * ssfactor, vn)
        
        % Save the indices to keep when cutting off endcaps
        save(keepfn, 'keep') ;
    else
        meshcut = read_ply_mod(outfn) ;
        faces = meshcut.f ;
        vtx = meshcut.v / ssfactor ; 
        vn = meshcut.vn ;
    end
        
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
    % Note that this assumes an elongated structure so that anterior
    % boundary is closer to acom than pcom and vice versa
    adist = sum((vtx(bb, :) - acom) .^ 2, 2) ;
    pdist = sum((vtx(bb, :) - pcom) .^ 2, 2) ;
    ab = bb(adist < pdist) ;
    pb = bb(adist > pdist) ;
    % check them
    if preview
        fig = figure;
        trisurf(faces, vtx(:, 1), vtx(:, 2), vtx(:, 3),...
            'edgecolor', 'none', 'facealpha', 0.1)
        hold on;
        plot3(vtx(ab, 1), vtx(ab, 2), vtx(ab, 3), 'b.-')
        plot3(vtx(pb, 1), vtx(pb, 2), vtx(pb, 3), 'r.-')
        plot3(acom(1), acom(2), acom(3), 'o')
        plot3(pcom(1), pcom(2), pcom(3), 'o')
        axis equal
        waitfor(fig)
    end
    
    % Save the anterior and posterior boundary indices (after endcap cut)
    save(boundaryfn, 'ab', 'pb')
    
    % choose anterior Dorsal as point with smallest phi
    % (closest to zero or 2pi) if this is first TP. Otherwise, match prev.
    if ii == 1
        % Choose dorsal point based on angle in yz plane. 
        % Since this is the first timepoint, taking angle with wrt x axis
        % is same as wrt AP axis given that AP axis is a straight line
        % presently.        
        % transform to APDV coords
        vrs = QS.xyz2APDV(vtx * ssfactor) ;
        % subtract pi/2 to make dorsal be zero
        a_phipi = mod(atan2(vrs(ab, 3), vrs(ab, 2)) - pi * 0.5, 2*pi) ;
        p_phipi = mod(atan2(vrs(pb, 3), vrs(pb, 2)) - pi * 0.5, 2*pi) ;
        % Now make the range go from -pi to pi with dorsal being zero
        a_phipi(a_phipi > pi) = a_phipi(a_phipi > pi) - 2 * pi ;
        p_phipi(p_phipi > pi) = p_phipi(p_phipi > pi) - 2 * pi ;
        [~, adb_tmp] = min(abs(a_phipi)) ;
        [~, pdb_tmp] = min(abs(p_phipi)) ;
        
        % anterior/posterior dorsal boundary
        adb = ab(adb_tmp) ;
        pdb = pb(pdb_tmp) ;
        
        % check it
        if preview
            clf; set(gcf, 'visible', 'on')
            trisurf(faces, vtx(:, 1), vtx(:, 2), vtx(:, 3), 0*vtx(:, 1),...
                'edgecolor', 'none', 'facealpha', 0.1)
            hold on;
            scatter3(vtx(ab, 1), vtx(ab, 2), vtx(ab, 3), 10, a_phipi)
            scatter3(vtx(pb, 1), vtx(pb, 2), vtx(pb, 3), 10, p_phipi)
            plot3(vtx(adb, 1), vtx(adb, 2), vtx(adb, 3), 'ks')
            plot3(vtx(pdb, 1), vtx(pdb, 2), vtx(pdb, 3), 'k^')
            caxis([-pi, pi])
            colorbar()
            axis equal
            title('Identification of dorsal points through phi')
            waitfor(gcf)
            close all 
        end
    else
        % Match previous timepoint
        ka = dsearchn(vtx(ab, :), previous_avtx) ;
        kp = dsearchn(vtx(pb, :), previous_pvtx) ;
        adb = ab(ka) ;
        pdb = pb(kp) ;
    end
    previous_avtx = vtx(adb, :) ;
    previous_pvtx = vtx(pdb, :) ;
        
    %% Save the anterior dorsal and posterior dorsal vertex points
    try 
        h5create(outapd_boundaryfn, ['/' name '/adorsal'], size(adb)) ;
    catch
        disp('adorsal pt already exists') ;
    end
    try
        h5create(outapd_boundaryfn, ['/' name '/pdorsal'], size(pdb)) ;
    catch
        disp('pdorsal pt already exists') ;
    end
    h5write(outapd_boundaryfn, ['/' name '/adorsal'], adb) ;
    h5write(outapd_boundaryfn, ['/' name '/pdorsal'], pdb) ;
    
    %% Plot the result
    figfn = fullfile(figoutdir, [name '.png']) ;
    if save_figs && (~exist(figfn, 'file') || overwrite)
        disp(['Saving figure for frame ' num2str(ii)])
        fig = figure('Visible', 'Off');
        vrs = QS.xyz2APDV(vtx * ssfactor) ;
        tmp = trisurf(faces, vrs(:, 1), vrs(:, 2), vrs(:, 3), ...
            vrs(:, 3), 'edgecolor', 'none', 'FaceAlpha', 0.1) ;
        hold on
        plot3(vrs(adb, 1), vrs(adb, 2), vrs(adb, 3), 's')
        plot3(vrs(pdb, 1), vrs(pdb, 2), vrs(pdb, 3), '^')
        axis equal
        xlim(xyzlim_um(1, :))
        ylim(xyzlim_um(2, :))
        zlim(xyzlim_um(3, :))
        xlabel('x [$\mu$m]', 'Interpreter', 'Latex')
        ylabel('y [$\mu$m]', 'Interpreter', 'Latex')
        zlabel('z [$\mu$m]', 'Interpreter', 'Latex')
        titlestr = 'Mesh with cylindrical cuts, $t=$' ;
        timestr = sprintf('%03d', tt * QS.timeinterval) ;
        titlestr = [titlestr timestr ' ' QS.timeunits] ;
        title(titlestr, 'Interpreter', 'Latex')
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

disp('done')
