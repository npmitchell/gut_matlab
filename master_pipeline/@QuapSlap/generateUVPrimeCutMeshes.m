function generateUVPrimeCutMeshes(QS, options)

%% Unpack QS
timePoints = QS.xp.fileMeta.timePoints ;

%% Unpack options
overwrite = true ;
save_ims = false ;
preview = false ;
if isfield(options, 'overwrite')
    overwrite = options.overwrite ; 
end
if isfield(options, 'save_ims')
    save_ims = options.save_ims ; 
end
if isfield(options, 'preview')
    preview = options.preview ; 
end
if save_ims
    imdir = fullfile(QS.dir.uvpcutMesh, 'images') ; 
    mkdir(imdir) ;
end
t0 = QS.t0set() ;


% First load all phi offsets phi0(u=0,1) from spcutMesh
for tidx = 1:length(timePoints)
    tp = timePoints(tidx) ;
    fn = sprintf(QS.fullFileBase.uvpcutMesh, tp) ;
    
    if ~exist(fn, 'file') || overwrite
        QS.setTime(tp) ;    

        % % OPTION 1: Could load spcutMesh.phi0_fit, but this would be before 
        % % smoothing. 
        % % To do so, first load cylMesh for this timepoint to pull back via 
        % conformal map
        % QS.loadCurrentCylinderMeshClean() ;
        % cylMesh = QS.currentMesh.cylinderMeshClean ;
        % % Extract offsets phi0 and phi1 to use for boundaries
        % phi0 = QS.cutMesh.phi0_fit(1) ; phi1 = QS.cutMesh.phi0_fit(end) ;
        % % Could then average phi0s in time by hand.
        % % Prep filter
        % disp('Building tripulse filter equivalent to tripuls(-0.5:0.1:0.5)')
        % tripulse = 0:0.2:1 ;
        % tripulse = [tripulse, fliplr(tripulse(1:end-1))] ;
        % tripulse = tripulse ./ sum(tripulse(:)) ;
        % tripulse = reshape(tripulse, [length(tripulse), 1]) ;

        % % OPTION 2:  Compute phi0 anew
        % uvals_boundary = [1e-9, max(spcutMeshSm.u(:, 1))] ;
        % tileCount = [1,1] ;
        % [ TF, TV2D, TV3D ] = tileAnnularCutMesh( spcutMeshSm, tileCount );
        % phiv_kk = [0 0] ;
        % prev3d_sphi_dsdv = mesh0.v ;
        % save_ims = false ;
        % plotfn = 'none' ;
        % smoothingMethod = 'none' ;
        % phiOpts = struct() ;
        % patchOpts = struct() ;
        % [phi0_fit_kk, phi0s_kk] = fitPhiOffsetsFromPrevMesh(TF,...
        %     TV2D, TV3D, uvals_boundary, phiv_kk, ...
        %     prev3d_sphi_dsdv, -0.45, 0.45, ...
        %     save_ims, plotfn, smoothingMethod, ...
        %     phiOpts, patchOpts) ;

        %% OPTION 3: Load spcutmesh, which already has incorporated phi(u=0,1)
        % Then we only need to conformally map the triangles to a rectilinear
        % domain and voila!
        QS.loadCurrentSPCutMeshSm() ;
        spcutMeshSm = QS.currentMesh.spcutMeshSm ;
        
        % Conformally map to disk
        rawMesh = flattenAnnulus(spcutMeshSm) ;   
        % Relax Affine transformation along x
        rawMesh.ar = minimizeIsoarealAffineEnergy( rawMesh.f, rawMesh.v, rawMesh.u );
        
        % Compute Beltrami coefficient (quasiconformal)
        affine_factor = rawMesh.ar ;
        
        for affine_factor = linspace(rawMesh.ar-0.1, rawMesh.ar+0.1, 10)
            v2d_raw = rawMesh.u ;
            v2d_raw(:, 1) = affine_factor * v2d_raw(:, 1) ;
            rawMesh.mu = bc_metric(rawMesh.f, v2d_raw, rawMesh.v, 3) ;
            disp([num2str(mean(real(rawMesh.mu))) '+/-' std(real(rawMesh.mu))])
        end
         uvpcutMesh.raw = rawMesh ;
        
        % Check against previously saved
        % tmp = load(sprintf(QS.fullFileBase.uvpcutMesh, tp)) ;
        % all(all(tmp.uvpcutMesh.raw.u == rawMesh.u))
        % all(all(tmp.uvpcutMesh.raw.v == rawMesh.v))

        % Make avgpts in pixel space (not RS), from rawMesh
        fprintf('Resampling uvgrid3d curves in pix...\n')
        nU = rawMesh.nU ;
        nV = rawMesh.nV ;
        curves3d_pix = reshape(rawMesh.v, [nU, nV, 3]) ;
        c3d_dsv_pix = zeros(size(curves3d_pix)) ;  % in units of pix
        avgpts_pix = zeros(nU, 3) ;
        radius_pix = zeros(nU, nV) ;
        for i=1:nU
            % Note: no need to add the first point to the curve
            % since the endpoints already match exactly in 3d and
            % curvspace gives a curve with points on either
            % endpoint (corresponding to the same 3d location).
            c3d_dsv_pix(i, :, :) = resampleCurvReplaceNaNs(squeeze(curves3d_pix(i, :, :)), nV, true) ;
            if vecnorm(squeeze(c3d_dsv_pix(i, 1, :)) - squeeze(c3d_dsv_pix(i, end, :))) > 1e-7
                error('endpoints do not join! Exiting')
            end
            avgpts_pix(i, :) = mean(squeeze(c3d_dsv_pix(i, :, :)), 1) ; 
            radius_pix(i, :) = vecnorm(squeeze(curves3d_pix(i, :, :)) - avgpts_pix(i, :), 2, 2) ;
        end
        
        % uvpcutMesh.raw.avgpts_pix = avgpts_pix ;
        % uvpcutMesh.raw.radius_pix = radius_pix ;
        uvpcutMesh.raw.avgpts_um = QS.xyz2APDV(avgpts_pix) ;
        uvpcutMesh.raw.radius_um = radius_pix * QS.APDV.resolution ;

        %% Compute rectified rectilinear grid in (u', v') pullback space   
        % First tile the raw mesh along v'
        tileCount = [1,1] ;
        [ TF, TV2D, TV3D ] = tileAnnularCutMesh( uvpcutMesh.raw, tileCount );

        onesUV = ones(nU, nV) ;
        uspace = linspace(0, 1, nU) ;
        vspace = linspace(0, 1, nV) ;
        uu = (uspace .* onesUV)' ;
        vv = vspace .* onesUV' ;
        % uv is equally spaced in u' and v'
        uv = [uu(:), vv(:)] ;
        uvgrid3d = interpolate2Dpts_3Dmesh(TF, TV2D, TV3D, uv) ;
        resMesh.nU = nU ;
        resMesh.nV = nV ;
        resMesh.f = rawMesh.f ;
        resMesh.u = uv ;
        resMesh.v = uvgrid3d ;
        resMesh.ar = rawMesh.ar ;
        resMesh.pathPairs = rawMesh.pathPairs ;
        v2d_res = resMesh.u ;
        v2d_res(:, 1) = rawMesh.ar * v2d_res(:, 1) ;
        resMesh.mu = bc_metric(resMesh.f, v2d_res, resMesh.v, 3) ;
        uvpcutMesh.resampled = resMesh ;
        assert(all(~any(isnan(resMesh.v))))
        assert(all(~any(isnan(resMesh.u))))

        % Interpolate radius onto resampled grid
        ugrid = reshape(uv(:, 1), [nU, nV]) ;
        radInterp = scatteredInterpolant(rawMesh.u(:, 1), rawMesh.u(:, 2),...
            uvpcutMesh.raw.radius_um(:), 'linear', 'nearest') ;
        resRadius_um = radInterp(resMesh.u(:, 1), resMesh.u(:, 2)) ;
        resRadius_um = reshape(resRadius_um, [nU, nV]) ;
        uvpcutMesh.resampled.radius_pix = resRadius_um ;
        
        % Check orientation
        if save_ims
            set(gcf, 'visible', 'off')
            imagesc(uspace, vspace, resRadius_um') ;
            axis equal
            axis tight
            xlabel('$u$', 'interpreter', 'latex')
            ylabel('$v$', 'interpreter', 'latex')
            cb = colorbar() ;
            ylabel(cb, ['radius [' QS.spaceUnits ']'], 'interpreter', 'latex') 
            imfn = [sprintf(QS.fileBase.uvpcutMesh, tp) '.png'] ;
            saveas(gcf, fullfile(imdir, imfn))
            close all
            
            % Save quasiconformal mesh -- raw
            set(gcf, 'visible', 'off')
            subplot(1, 2, 1)
            options.labels = {'$\Re \mu$', '$\Im \mu$'} ;
            [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
                twoScalarFieldsOnSurface({rawMesh.f, ...
                [v2d_raw(:, 1), v2d_raw(:, 2), 0*v2d_raw(:, 1)]}, ...
                real(rawMesh.mu), imag(rawMesh.mu), options) ;
            sgtitle(['$\mu($embedding, pullback$)$, $t = $', ...
                sprintf('%03d', tp-t0), ' ', QS.timeUnits], ...
                'interpreter', 'latex') ;
            set(gcf,'CurrentAxes', ax1)
            view(2)
            axis off
            set(gcf,'CurrentAxes', ax2)
            view(2)
            axis off
            imfn_raw = ['mu_raw_' sprintf(QS.fileBase.uvpcutMesh, tp) '.png'] ;
            disp(['Saving ' imfn_raw ': ' fullfile(imdir, imfn_raw)])
            saveas(gcf, fullfile(imdir, imfn_raw))
            
            set(gcf, 'visible', 'off')
            subplot(1, 2, 1)
            options.labels = {'$\Re \mu$', '$\Im \mu$'} ;
            [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
                twoScalarFieldsOnSurface({resMesh.f, ...
                [v2d_res(:, 1), v2d_res(:, 2), 0*v2d_res(:, 1)]}, ...
                real(resMesh.mu(tidx, :)), imag(resMesh.mu(tidx, :)), options) ;
            sgtitle(['resampled mesh: $\mu($embedding, pullback$)$, $t = $', ...
                sprintf('%03d', tp-t0), ' ', QS.timeUnits], ...
                'interpreter', 'latex') ;
            set(gcf,'CurrentAxes', ax1)
            view(2)
            set(gcf,'CurrentAxes', ax2)
            view(2)
            imfn_res = ['mu_res_' sprintf(QS.fileBase.uvpcutMesh, tp) '.png'] ;
            saveas(gcf, fullfile(imdir, imfn_res))
        end
        
        if preview
            clf
            scatter(rawMesh.u(:, 1), rawMesh.u(:, 2), 10, radius_um(:)) 
            hold on;
            scatter(resMesh.u(:, 1), resMesh.u(:, 2), 10, resRadius_um(:))
            cb = colorbar() ;
            ylabel(cb, ['radius [' QS.spaceUnits ']'], 'interpreter', 'latex')
            pause(1)
        end
                
        %% Save the result
        disp(['Saving uvpcutMesh t=' num2str(tp) ': ' fn])
        save(fn, 'uvpcutMesh')
    end
end
    
    
% Notes to self:
% uvpcutMesh = flattenAnnulusTiltedBoundaries(cutMesh, phi0, phi1, 'Dirichlet') ;
% ar = minimizeIsoarealAffineEnergy( cutMesh.f, cutMesh.v, cutMesh.u );
