function generateCurrentSPCutMesh(QS, cutMesh, overwrite)
%
%
%
%

% Unpack options
if nargin < 2 || isempty(cutMesh)
    if isempty(QS.currentMesh.cutMesh)
        QS.loadCurrentCutMesh()
    end
    cutMesh = QS.currentMesh.cutMesh ;
    overwrite = false ;
end
if nargin < 3
    overwrite = false ;
end

% Unpack QS
tt = QS.currentTime ;
nU = QS.nU ;
nV = QS.nV ;
spcutMeshfn = sprintf(QS.fullFileBase.spcutMesh, tt) ;
fileNameBase = QS.fileBase.name ; 
phi0fitBase = QS.fullFileBase.phi0fit ;
[rot, trans] = getRotTrans(QS) ;
resolution = QS.APDV.resolution ;
QS.getCleanCntrlines ;
cleanCntrlines = QS.cleanCntrlines ;
cleanCntrline = cleanCntrlines{QS.xp.tIdx(tt)} ;
preview = QS.plotting.preview ;
[~, ~, xyzlim_um] = QS.getXYZLims() ;
clineDVhoopDir = QS.dir.clineDVhoop ;

% Expand dirs for images
clineDVhoopImDir = fullfile(clineDVhoopDir, 'images') ;
clineDVhoopFigBase = fullfile(clineDVhoopImDir, 'clineDVhoop_%06d.png') ;
if ~exist(clineDVhoopImDir, 'dir')
    mkdir(clineDVhoopImDir)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate s,phi coord system for rotated, scaled mesh (rs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Establishing s,phi coord system\n');
if ~exist(spcutMeshfn, 'file') || overwrite
    if overwrite
        disp('Overwriting spcutMesh...')
    else
        disp('spcutMesh not on disk. Generating ...')
    end

    % Transform from u,v coordinates to s, phi coordinates
    % [scoords, phicoords] = generateSPhiFromUV();

    %----------------------------------------------------------------------
    % Generate tiled orbifold triangulation
    %----------------------------------------------------------------------
    tileCount = [1 1];  % how many above, how many below
    cutMeshrs = cutMesh;
    % Rotate and translate TV3D
    cutMeshrs.v = ((rot * cutMesh.v')' + trans) * resolution ;
    cutMeshrs.vn = (rot * cutMesh.vn')' ;
    [ ~, ~, TV3D, TVN3D ] = tileAnnularCutMesh( cutMesh, tileCount );
    [ TF, TV2D, TV3Drs ] = tileAnnularCutMesh( cutMeshrs, tileCount );

    %----------------------------------------------------------------------
    % Calculate abbreviated centerline from cutMesh boundaries
    %----------------------------------------------------------------------
    % Load centerline from cleaned list
    cline = QS.xyz2APDV(cleanCntrline(:, 2:4)) ; 
    ss = cleanCntrline(:, 1) * QS.APDV.resolution ;
    disp('Finding relevant segment of centerline')
    [cseg, acID, pcID, ~, ~] = ...
        centerlineSegmentFromCutMesh(cline, TF, TV2D, TV3Drs) ;

    %----------------------------------------------------------------------
    % Generate surface curves of constant s
    %----------------------------------------------------------------------
    % For lines of constant phi
    disp('Creating crude uv curves with du=const to define uspace by ds(u)')
    % Make grid
    eps = 1e-14 ;
    uspace0 = linspace( eps, cutMesh.umax - eps, nU )' ;
    vspace = linspace( eps, 1-eps, nV )' ;

    disp('Casting crude (equal dU) points into 3D...')
    crude_ringpath_ss = ringpathsGridSampling(uspace0, vspace, TF, TV2D, TV3Drs) ;

    % Resample crude_ringpath_ds made from uspace0 (equal du, not equal ds_3D in u direction)
    [uspace, eq_ringpath_ss] = equidistantSampling1D(linspace(0, 1, nU)', crude_ringpath_ss, nU, 'linear') ;
    % ensure that uspace is nU x 1, not 1 x nU
    uspace = reshape(uspace, [nU, 1]) ; 
    % hedge the first and last point to avoid NaNs
    eps = 1e-13 ;
    uspace(1) = uspace(1) + eps ;
    uspace(end) = uspace(end) - eps ;
    clearvars dsuphi curves3d uspace0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Casting resampled points into 3D (approx equal ds_3D in u dir, but variable ds_3D in v dir)...')
    % NOTE: first dimension indexes u, second indexes v
    curves3d = zeros(nU, nV, 3) ;  % in units of um
    for kk = 1:nU
        if mod(kk, 50) == 0
            disp(['u = ' num2str(kk / nU)])
        end
        uv = [cutMesh.umax * uspace(kk) * ones(size(vspace)), vspace] ;
        curves3d(kk, :, :) = interpolate2Dpts_3Dmesh(TF, TV2D, TV3Drs, uv) ;
    end 

    % Check the 3d curves 
    if preview
        figure ; hold on;
        for kk = 1:nU
            plot3(curves3d(kk, :, 1), curves3d(kk, :, 2), curves3d(kk, :, 3), '.') 
        end
        axis equal
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Compute s(u) and radius(u) for "uniform"--> evenly sample each DV hoop (0,1) so ds_3D=const \n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Resample at evenly spaced dphi in embedding space (rs, in um)
    fprintf('Resampling curves...\n')
    c3d_dsv = zeros(size(curves3d)) ;  % in units of um
    for i=1:nU
        % Note: no need to add the first point to the curve
        % since the endpoints already match exactly in 3d and
        % curvspace gives a curve with points on either
        % endpoint (corresponding to the same 3d location).
        c3d_dsv(i, :, :) = resampleCurvReplaceNaNs(squeeze(curves3d(i, :, :)), nV, true) ;
        if vecnorm(squeeze(c3d_dsv(i, 1, :)) - squeeze(c3d_dsv(i, end, :))) > 1e-7
            error('endpoints do not join! Exiting')
        end

        % Visualization for Troubleshooting:
        % triplot(TF, TV2D(:, 1), TV2D(:, 2))
        % hold on;
        % plot(uv(:, 1), uv(:, 2), '.')
    end

    % Check the 3d curves 
    if preview
        figure ; hold on;
        for kk = 1:nU
            plot3(c3d_dsv(kk, :, 1), c3d_dsv(kk, :, 2), c3d_dsv(kk, :, 3), '.') 
        end
        axis equal
    end

    fprintf('Finding s(u) and r(u) of resampled "uniform" c3ds [uniform ds in V dir]...\n')
    % mcline is the resampled centerline, with mss
    % avgpts is the raw Nx3 averaged hoops, with avgpts_ss
    [mss, mcline, radii_from_mean_uniform_rs, avgpts_ss, avgpts] = srFromDVCurves(c3d_dsv) ;

    % Used to find radius using original centerline
    % [ssv, radii, avgpts, cids] = srFromDVCurvesGivenCenterline(ss, cline, c3ds) ;
    % Could operate just on the centerline segment
    cseg_ss = ss(acID:pcID) ;
    % [ssv, radii, avgpts, cids] = srFromDVCurves(cseg_ss, cseg, c3ds) ;
    % 
    % Adjust the centerline indices to index into the full
    % centerline. Note that cseg_ss already does this for ss.
    % cids = cids + acID ;

    % Plot new centerline
    aux_plot_clineDVhoop(avgpts, avgpts_ss, cseg, cline, cseg_ss, ...
        curves3d, xyzlim_um, clineDVhoopFigBase, tt)

    % Optional: clean curve with polynomial and point match
    % avgpts onto cleaned curve. Skipping for later.

    % Compute ringpath_ss, the mean distance traveled from one
    % line of constant u to the next
    disp('Computing ringpath_ss in "uniform" resampling (equal ds along DV)...')
    % The distance from one hoop to another is the
    % difference in position from (u_i, v_i) to (u_{i+1}, v_i).
    dsuphi = reshape(vecnorm(diff(c3d_dsv), 2, 3), [nU-1, nV]) ;
    ringpath_ds = nanmean(dsuphi, 2) ;
    ringpath_ss = cumsum([0; ringpath_ds]) ;
    clearvars dsuphi ringpath_ds

    % Save new centerline in rotated translated units
    fn = sprintf(clineDVhoopBase, t) ;
    disp(['Saving new centerline to ' fn])
    save(fn, 'mss', 'mcline', 'avgpts', 'avgpts_ss')

    % Note: radii_from_mean_uniform_rs is the radius of 
    % interpolated hoops, not the actual points

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Done making new centerline using uniformly sampled hoops\n') ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Create new3d, the regridded pts at UV, moved to sphi')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    onesUV = ones(nU, nV) ;
    uu = uspace * cutMesh.umax .* onesUV ;
    vv = (vspace .* onesUV')' ;
    uv = [uu(:), vv(:)] ;
    % Note: here interpolate uv in the TV2D coord system, then
    % use uphi as the actual 2D coordinates for these vertices
    % NOTE: unlike curves3d, new3d is NOT rotated/translated/scaled
    new3d = interpolate2Dpts_3Dmesh(TF, TV2D, TV3D, uv) ;

    IVloaded = false ;
    if t == xp.fileMeta.timePoints(1)
        % Store for next timepoint
        phiv = (vspace .* ones(nU, nV))' ;
        phi0s = zeros(size(uspace)) ;
        phi0_fit = phi0s ;
    else
        % Load previous sphi vertices in 3d 
        plotfn = sprintf(phi0fitBase, tt, 0);
        if strcmp(phi_method, '3dcurves')
            % Load the previous spcutMesh and call it prev3d_sphi
            % Also note the previous spcutMesh pullback image's fn
            tmp = load(sprintf(spcutMeshBase, ...
                xp.fileMeta.timePoints(tidx-1)), 'spcutMesh') ;
            prevf = tmp.spcutMesh.f ;
            prev3d_sphi = reshape(tmp.spcutMesh.v, [nU, nV, 3]) ; 
            imfn_sp_prev = sprintf( ...
                fullfile([imFolder_sp, '/', fileNameBase, '.tif']), ...
                xp.fileMeta.timePoints(tidx-1) ) ;

            % fit the shifts in the y direction
            dmyk = 0 ;
            phi0_fit = zeros(size(uspace)) ;
            phi0s = zeros(size(uspace)) ;
            phi0_fit_kk = 1 ; % for first pass                
            phiv_kk = (vspace .* ones(nU, nV))' ;
            ensureDir([sphiDir, '/phi0_correction/'])
            while any(phi0_fit_kk > 0.002) && dmyk < 6
                disp(['Iteration ' num2str(dmyk)])
                plotfn = sprintf(phi0fitBase, tt, dmyk);

                % Will we save check pullbacks to preview the algo?
                if save_phi0patch
                    patchImFn = sprintf( ...
                        fullfile(sphiDir, 'phi0_correction', [fileNameBase, '_prephi0_' num2str(dmyk) '.tif']), ...
                        xp.fileMeta.timePoints(tidx-1) )  ;
                    geomImFn = sprintf( ...
                        fullfile(sphiDir, 'phi0_correction', ['3d' fileNameBase '_prephi0_' num2str(dmyk) '.tif']), ...
                        xp.fileMeta.timePoints(tidx-1) )  ;

                    % Load the intensity data for this timepoint
                    if ~IVloaded
                        % (3D data for coloring mesh pullback)
                        xp.loadTime(t);
                        xp.rescaleStackToUnitAspect();

                        % Raw stack data
                        IV = xp.stack.image.apply();
                        IV = imadjustn(IV{1});         
                        IVloaded = true ;
                    end

                    % Texture patch options
                    Options.PSize = 5;
                    Options.EdgeColor = 'none';
                    % Texture image options
                    Options.imSize = ceil( 1000 .* [ 1 a_fixed ] );
                    Options.yLim = [0 1];

                    % Roll options into a struct
                    patchOpts.patchImFn = patchImFn ;
                    patchOpts.imfn_sp_prev = imfn_sp_prev ;
                    patchOpts.IV = IV ;
                    patchOpts.ringpath_ss = ringpath_ss ;
                    patchOpts.Options = Options ;
                    patchOpts.v3d = new3d ;
                else
                    patchOpts = [] ;
                end

                % Minimize difference in DV hoop positions wrt
                % previous pullback mesh                        
                [phi0_fit_kk, phi0s_kk] = fitPhiOffsetsFromPrevMesh(TF, TV2D, TV3D, ...
                    uspace * cutMesh.umax, phiv_kk, prev3d_sphi, -0.45, 0.45, ...
                    save_ims, plotfn, save_phi0patch, preview, patchOpts) ;

                % Update the result
                dmyk = dmyk + 1;
                phi0_fit = phi0_fit + phi0_fit_kk ;
                phi0s = phi0s + phi0s_kk ;
                phiv_kk = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;


                % plot mesh colored by the phase phi 
                % previous timepoint
                xtmp = prev3d_sphi(:, :, 1) ;
                ytmp = prev3d_sphi(:, :, 2) ;
                ztmp = prev3d_sphi(:, :, 3) ;
                phitmp = (vspace .* ones(nU, nV))' ;
                colormap parula ;
                % cmap = parula ;
                % colors = cmap(max(1, uint8(colortmp(:) * length(parula))), :) ;
                trisurf(prevf, xtmp(:), ytmp(:), ztmp(:), phitmp(:), ...
                    'FaceColor', 'interp',...
                    'EdgeColor', 'none', 'FaceAlpha', 0.25)
                axis equal
                % freezeColors

                % before fitting
                hold on;
                pe0 = find(phitmp(:) < 1e-4 | phitmp(:) > 0.99) ;
                plot3(new3d(pe0, 1), new3d(pe0, 2), ...
                    new3d(pe0, 3), '.')
                % colormap copper
                % trimesh(prevf, new3d(inds, 1), new3d(:, 2), new3d(:, 3),...
                %     phitmp(:), 'FaceColor', 'interp', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
                % freezeColors

                % after fitting
                hold on;                   
                pekk = find(mod(phiv_kk(:), 1) < 1e-4 | mod(phiv_kk(:), 1) > 0.99) ;
                plot3(new3d(pekk, 1), new3d(pekk, 2), ...
                    new3d(pekk, 3), '^')                        
                % colormap summer
                % trimesh(prevf, new3d(:, 1), new3d(:, 2), new3d(:, 3),...
                %     mod(phiv_kk(:), 1), 'FaceColor', 'interp', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
                % freezeColors
                xlabel('x [\mum]')
                ylabel('y [\mum]')
                zlabel('z [\mum]')
                view(2)
                saveas(gcf, geomImFn)
            end

        elseif strcmp(phi_method, 'texture')
            imfn_sp_prev = sprintf( ...
                fullfile([imFolder_sp, '/', fileNameBase, '.tif']), ...
                xp.fileMeta.timePoints(tidx-1) ) ;

            % Load the intensity data            
            % Load 3D data for coloring mesh pullback
            xp.loadTime(t);
            xp.rescaleStackToUnitAspect();

            % Raw stack data
            IV = xp.stack.image.apply();
            IV = imadjustn(IV{1});         
            IVloaded = true ;

            % Texture patch options
            Options.PSize = 5;
            Options.EdgeColor = 'none';
            % Texture image options
            Options.imSize = ceil( 1000 .* [ 1 a_fixed ] );
            Options.yLim = [0 1];

            % fit the shifts in the y direction
            % todo: save uncorrected patchIms,
            % could try tiling twice...
            dmyk = 0 ;
            phi0_fit = zeros(size(uspace)) ;
            phi0s = zeros(size(uspace)) ;
            phi0_fit_kk = 1 ; % for first pass                
            phiv_kk = (vspace .* ones(nU, nV))' ;
            ensureDir([sphiDir, '/phi0_correction/'])
            while any(phi0_fit_kk > 0.002) && dmyk < 6
                disp(['Iteration ' num2str(dmyk)])
                plotfn = sprintf(phi0fitBase, tt, dmyk);
                patchImFn = sprintf( ...
                    fullfile([sphiDir, '/phi0_correction/', fileNameBase, '_prephi0_' num2str(dmyk) '.tif']), ...
                    xp.fileMeta.timePoints(tidx-1) )  ;
                [phi0_fit_kk, phi0s_kk] = fitPhiOffsetsFromPrevPullback(IV, ...
                    new3d, cutMesh.umax, uspace * cutMesh.umax, phiv_kk, ...
                    ringpath_ss, imfn_sp_prev, lowerboundy, upperboundy, ...
                    save_ims, plotfn, Options, ...
                    step_phi0tile, width_phi0tile, potential_sigmay, 'integer', ...
                    patchImFn) ;

                % Update the result
                dmyk = dmyk + 1;
                phi0_fit = phi0_fit + phi0_fit_kk ;
                phi0s = phi0s + phi0s_kk ;
                phiv_kk = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;
            end
        else
            error("Could not recognize phi_method: must be 'texture' or '3dcurves'")
        end
        close all

        % Store to save at this timepoint
        phiv = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;
    end

    % NOTE: We have coordinates u,phiv that we associate with
    % the 3d coordinates already mapped to uv
    uphi = [uu(:), phiv(:)] ;

    % plot(uphi(:, 1), uphi(:, 2), '.')
    % xlabel('u')
    % ylabel('\phi')
    % waitfor(gcf)

    % Recompute radii_from_mean_uniform_rs as radii_from_avgpts 
    % NOTE: all radius calculations done in microns, not pixels
    sphi3d_rs = ((rot * new3d')' + trans) * resolution ;
    radii_from_avgpts = zeros(size(sphi3d_rs, 1), size(sphi3d_rs, 2)) ;
    for jj = 1:nU
        % Consider this hoop
        hoop = squeeze(sphi3d_rs(jj, :, :)) ;
        radii_from_avgpts(jj, :) = vecnorm(hoop - avgpts(jj, :), 2, 2) ;
    end

    % Triangulate the sphigrid and store as its own cutMesh
    % sphiv = zeros(nU, nV, 2) ;
    % sphiv(:, :, 1) = sv ;
    % sphiv(:, :, 2) = phiv ;
    sv = ringpath_ss .* onesUV ;
    % % Triangulate the mesh (topology is already known):
    % tmptri = defineFacesRectilinearGrid(sp, nU, nV) ;
    % % Old version did not assume topology as given:
    % tmptri = delaunay(sv(:), phiv(:)) ;
    % disp('orienting faces of delaunay triangulation (s,phi)')
    % tmptri = bfs_orient( tmptri );

    % Define path pairs for tiling the (s,phi) cut mesh
    spcutP1 = 1:nU;
    spcutP2 = nU*nV - fliplr(0:(nU-1)) ;
    spcutMesh.pathPairs = [ spcutP1', spcutP2' ];

    % Check to see if any members of pathPairs connect to
    % non-Nearest Neighbors. Not necessary now that we assume
    % known gridded mesh topology
    % cleantri = cleanBoundaryPath2D(tmptri, [sv(:), phiv(:)], spcutMesh.pathPairs(:), true) ;

    spcutMesh.f = defineFacesRectilinearGrid(uv, nU, nV) ;
    spcutMesh.nU = nU ;
    spcutMesh.nV = nV ;
    % First resampling
    spcutMesh.v0 = new3d ;
    % spcutMesh.vrs0 = ((rot * new3d')' + trans) * resolution ;
    % Define normals based on the original mesh normals
    spvn03d = interpolate2Dpts_3Dmesh(TF, TV2D, TVN3D, uphi) ;
    spvn03d = spvn03d ./ vecnorm(spvn03d, 2, 2) ;
    spcutMesh.vn0 = spvn03d ;
    spcutMesh.sphi0 = [sv(:), phiv(:)] ;
    spcutMesh.uphi0 = uphi ;
    % Note: uv has no direct relation with cutMesh, just a grid
    % for utility and reference, but it does have unequal 
    % spacing in u in anticipation of building sphi0 as a 
    % near perfect grid.
    spcutMesh.uv = uv ;  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SECOND RESAMPLING
    % Make a new grid
    slin = linspace(0, max(spcutMesh.sphi0(:, 1)), nU) ;
    plin = linspace(0, 1, nV) ;
    [ss, pp] = meshgrid(slin, plin) ;
    % Push the endpoints on each boundary in by epsilon to
    % avoid NaNs
    eps = 1e-14 ;
    ss(:, 1) = eps ;
    ss(:, end) = ss(:, end) - eps ;
    % Transpose so that x increases with increasing index first
    ss = ss' ;
    pp = pp' ;
    sp = [ss(:), pp(:)] ;

    % Tile the spcutMesh
    tileCount = [2, 2] ;
    spcutMesh.u = spcutMesh.sphi0 ;
    spcutMesh.v = spcutMesh.v0 ;
    spcutMesh.vn = spcutMesh.vn0 ;
    [ faces, v2d, v3d, vn3d ] = tileAnnularCutMesh( spcutMesh, tileCount );
    spcutMesh = rmfield(spcutMesh, 'u') ;
    spcutMesh = rmfield(spcutMesh, 'v') ;
    spcutMesh = rmfield(spcutMesh, 'vn') ;
    spv3d = interpolate2Dpts_3Dmesh(faces, v2d, v3d, sp) ;
    % check the pts
    % plot3(spv3d(:, 1), spv3d(:, 2), spv3d(:, 3))  

    % also interpolate the normals
    spvn3d = interpolate2Dpts_3Dmesh(faces, v2d, vn3d, sp) ;
    spvn3d = spvn3d ./ vecnorm(spvn3d, 2, 2) ;

    % Define new faces for second rectilinear resampling
    % NOTE: not necessary since we already defined the topology
    % from the guess [sv(:), phiv(:)] stored as spcutMesh.sphi0
    % spcutMesh.f = defineFacesRectilinearGrid(sp, nU, nV) ;
    spcutMesh.sphi = sp ;
    spcutMesh.v = spv3d ;
    spcutMesh.vn = spvn3d ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    spcutMesh.ringpath_ss = ringpath_ss ;
    spcutMesh.radii_from_mean_uniform_rs = radii_from_mean_uniform_rs ;  % from uniform DV sampling
    spcutMesh.radii_from_avgpts = radii_from_avgpts ;
    spcutMesh.mss = mss ;       % from uniform DV sampling, also stored in centerline
    spcutMesh.mcline = mcline ; % from uniform DV sampling, also stored in centerline
    spcutMesh.avgpts = avgpts ; % from uniform DV sampling, also stored in centerline
    spcutMesh.avgpts_ss = avgpts_ss ; % from uniform sampling, also stored in centerline

    % Define optimal isoareal Affine dilation factor in s
    % tmp = spcutMesh.sphi ;
    % tmp(:, 1) = tmp(:, 1) / max(tmp(:, 1)) ;
    % arsp = minimizeIsoarealAffineEnergy( spcutMesh.f, spcutMesh.v, tmp );
    % clearvars tmp
    spcutMesh.ar = cutMesh.ar ;

    % todo: check that u coords have not shifted upon
    % redefinition of sphi0 -> sphi

    % Save s,phi and their 3D embedding
    spcutMesh.phi0s = phi0s ;
    spcutMesh.phi0_fit = phi0_fit ;
    save(spcutMeshfn, 'spcutMesh') ;
else
    disp('Loading spcutMesh from disk...')
    load(spcutMeshfn, 'spcutMesh') ;
    IVloaded = false ;

    % Load new centerline
    fn = sprintf(clineDVhoopBase, t) ;
    disp(['Loading new centerline from ' fn])
    load(fn, 'mss', 'mcline', 'avgpts', 'avgpts_ss')
end
fprintf('Done with generating S,Phi coords \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
