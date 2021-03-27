function plotMetric(QS, options)
%plotMetric(QS, options)
%
% Parameters
% ----------
%
% Returns
% -------
%
% NPMitchell 2020

t0 = QS.t0set() ;
lambda_mesh = QS.smoothing.lambda_mesh ;
overwrite = false ;
makeRawMetricComponentFigures = true ;
coordSys = 'spsm_rs' ;

% Unpack options
if isfield(options, 'coordSys')
    coordSys = options.coordSys ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'makeRawMetricComponentFigures')
    makeRawMetricComponentFigures = options.makeRawMetricComponentFigures ;
end

for tidx = 1:length(QS.xp.fileMeta.timePoints)
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    disp(['t = ', num2str(tp)])
    
    strClims = {'climVariable', 'climUniform'} ;
    outdirs = {sprintf(QS.dir.metric.data, coordSys, lambda_mesh), ...
        fullfile(sprintf(QS.dir.metric.g_images3d, coordSys, lambda_mesh), strClims{1}), ...
        fullfile(sprintf(QS.dir.metric.b_images3d, coordSys, lambda_mesh), strClims{1}), ...
        fullfile(sprintf(QS.dir.metric.g_images3d, coordSys, lambda_mesh), strClims{2}), ...
        fullfile(sprintf(QS.dir.metric.b_images3d, coordSys, lambda_mesh), strClims{2}), ...
        fullfile(sprintf(QS.dir.metric.g_images2d, coordSys, lambda_mesh), strClims{1}), ...
        fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), strClims{1}), ...
        fullfile(sprintf(QS.dir.metric.g_images2d, coordSys, lambda_mesh), strClims{2}), ...
        fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), strClims{2}), ...
        fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), 'director'), ...
        fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), 'Hopf_differential')};
    for qq = 1:length(outdirs)
        if ~exist(outdirs{qq}, 'dir')
            mkdir(outdirs{qq}) ;
        end
    end
    
    outfn = sprintf(QS.fullFileBase.metric, coordSys, lambda_mesh, tp) ;
    
    if ~exist(outfn, 'file') || overwrite
        % Load mesh
        if strcmpi(coordSys, 'spsm')
            mesh = load(sprintf(QS.fullFileBase.spcutMeshSm, tp), 'spcutMeshSm') ;
            mesh = mesh.spcutMeshSm ;
        elseif strcmpi(coordSys, 'spsm_rs')
            mesh = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp), 'spcutMeshSmRS') ;
            mesh = mesh.spcutMeshSmRS ;
        else
            error('handle this coordSys here')
        end
                
        % Smooth the mesh with lambda_mesh
        if lambda_mesh > 0 
            disp('smoothing mesh vertices before computations')
            tri = triangulation(mesh.f, mesh.v) ;
            fbndy = tri.freeBoundary ;
            fbndy = fbndy(:, 1) ;
            mesh.v = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
                lambda_mesh, 'implicit', mesh.v) ;
        end

        % Rescale u dimension for mu relaxation
        uv = mesh.u ;
        uv(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) ;
        mesh.readme = 'mu is computed after rescaling refMesh.u(:, 1) to range from 0 to ar, as determined by mean pure shear in (s,phi)' ;
        mesh.mu = bc_metric(mesh.f, uv, mesh.v, 3) ;

        % Project to APDV coords if not already done
        try 
            mesh.vrs;
        catch
            mesh.vrs = QS.xyz2APDV(mesh.v) ;
        end

        % Rescale and recompute mu
        aspectShear = (1 + mean(real(mesh.mu))) / (1 - mean(real(mesh.mu))) ;
        mesh.u(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) * aspectShear ; 
        mesh.mu = bc_metric(mesh.f, mesh.u, mesh.v, 3) ;

        [gcell, bcell] = constructFundamentalForms(mesh.f, mesh.vrs, mesh.u) ;
        gg = zeros(length(gcell), 4) ;
        bb = zeros(length(gcell), 4) ;
        for qq = 1:length(gcell)
            gg(qq, 1) = gcell{qq}(1, 1) ;
            gg(qq, 2) = gcell{qq}(1, 2) ;
            gg(qq, 3) = gcell{qq}(2, 1) ;
            gg(qq, 4) = gcell{qq}(2, 2) ;
            bb(qq, 1) = bcell{qq}(1, 1) ;
            bb(qq, 2) = bcell{qq}(1, 2) ;
            bb(qq, 3) = bcell{qq}(2, 1) ;
            bb(qq, 4) = bcell{qq}(2, 2) ;
        end

        % Save fundamental forms with mesh
        save(outfn, 'gg', 'bb', 'mesh', 'aspectShear')
    else
        load(outfn, 'gg', 'bb', 'mesh', 'aspectShear')
    end
    
    %% Second fundamental form --> mean curvature and Hopf differential
    fnbb = fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), ...
        'Hopf_differential', sprintf(QS.fileBase.name, tp)) ;        
    fnbb = [fnbb '_b_fundForm_Hopf.png'] ;

    if ~exist(fnbb, 'file') || overwrite || true
        close all
        opts = struct() ;
        
        % glue mesh to topological cylinder
        glueMesh = glueCylinderCutMeshSeam(mesh) ;
        assert(eulerCharacteristic(glueMesh) == 0)
        [~, F2V] = meshAveragingOperators(glueMesh) ;
        
        % Get director as Hopf differential
        LL = bb(:, 1) ;
        MM = bb(:, 2) * 0.5 ;
        NN = bb(:, 4) ;
        QQ = 0.25 * ((LL-NN) - 2*1j*MM ) ;
        
        % Check with bars
        % bc = barycenter(mesh.u, mesh.f) ;
        % bar1 = real(QQ) ;
        % bar1 = cat(2, bar1, imag(QQ)) ;
        % quiver(bc(:, 1), bc(:, 2), bar1(:, 1), bar1(:, 2), 1) ;

        % Check with scatter
        % subplot(2, 1, 1) ;
        % scatter(bc(:, 1), bc(:, 2), 10, real(QQ))
        % subplot(2, 1, 2) ;
        % scatter(bc(:, 1), bc(:, 2), 10, imag(QQ))
        
        % Get director on vertices
        nU = mesh.nU ;
        nV = mesh.nV ;
        QV = zeros(length(mesh.u), 2) ;
        evs = QV ;
        bbV = zeros(length(mesh.u), 4) ;
        bbV = F2V * bb ;
        bbV(nU*(nV-1)+1:nU*nV, :) = bbV(1:nU, :) ;
        
        LV = bbV(:, 1) ;
        MV = bbV(:, 2) * 0.5 ;
        NV = bbV(:, 4) ;
        QV = 0.25 * ((LV-NV) - 2*1j*MV ) ;
        
        % Check against moving QV directly -- checks out!
        % QV2 = F2V * QQ ;
        % QV2(nU*(nV-1)+1:nU*nV, :) = QV(1:nU, :) ;
        
        % Filter a bit
        QVs_r = imgaussfilt(reshape(real(QV), [nU, nV]), 1) ;
        QVs_i = imgaussfilt(reshape(imag(QV), [nU, nV]), 1) ;
        QVs = QVs_r + 1j * QVs_i ;
        
        mags = abs(QVs) ;
        thetas = atan2(imag(QVs), real(QVs)) ;
        uu = reshape(mesh.u(:, 1), [nU, nV]) ;
        vv = reshape(mesh.u(:, 2), [nU, nV]) ;
        QVs = reshape(QVs, [nU, nV]) ;
        
        % Check with polar field
        % opts = struct() ;
        % opts.clim_mag = rms1d(mags(:)) ;
        % plotPolarField(mags', thetas', opts)
        
        % Plot as polar in 2d AND 3d
        m2d = struct() ;
        m2d.f = mesh.f ;
        m2d.v = mesh.u ;
        m2d.v(:, 3) = 0 * mesh.u(:, 1) ;
        opts.clim = rms1d(mags(:, 1)) ;
        opts.axisOff = true ;
        opts.labels = {'', '2nd fundamental form as nematic'} ;
        opts.views = {[0, 0 ], [0, 90 ]};
        opts.cbarlabels = {['Hopf differential, $|Q|$ [' QS.spaceUnits '$^{-2}$]'], ...
                           ['Hopf differential, $|Q|$ [' QS.spaceUnits '$^{-2}$]']} ;
        opts.polarStyle = 'polar';
        close all
        nFieldsOnSurface({mesh, m2d}, ...
            {{mags, thetas}, {mags, thetas}}, opts) ;
        
        sgtitle(['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
                ' ', QS.timeUnits], 'Interpreter', 'latex')
        disp(['Saving ' fnbb])
        saveas(gcf, fnbb)
        close all
        
        
        fnbb2 = fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), ...
            'Hopf_differential', sprintf(QS.fileBase.name, tp)) ;        
        fnbb2 = [fnbb2 '_b_fundForm_Hopf_halfAngle.png'] ;
        opts.polarStyle = 'nematic';
        
        nFieldsOnSurface({mesh, m2d}, ...
            {{mags, thetas + pi*0.5}, {mags, thetas + pi*0.5}}, opts) ;
        sgtitle(['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
                ' ', QS.timeUnits], 'Interpreter', 'latex')
        disp(['Saving ' fnbb2])
        saveas(gcf, fnbb2)
        close all

    end
    
    
    %% SAVE images of second fundamental form director if they dont exist
    fnbb = fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), ...
        'director', sprintf(QS.fileBase.name, tp)) ;        
    fnbb = [fnbb '_b_fundForm_director.png'] ;

    if ~exist(fnbb, 'file') || overwrite
        close all
        opts = struct() ;
        
        % glue mesh to topological cylinder
        glueMesh = glueCylinderCutMeshSeam(mesh) ;
        assert(eulerCharacteristic(glueMesh) == 0)
        [~, F2V] = meshAveragingOperators(glueMesh) ;
        
        % Get director
        QQ = zeros(length(bb), 2) ;
        evs = QQ ;
        for qq = 1:length(bb) 
            b2d = reshape(bb(qq, :), [2, 2]) ;
            [eigvect, eigval] = eig(b2d) ;
            % larger eigvect direction
            [~, id] = max(abs([eigval(1, 1), eigval(2, 2)])) ;
            if sign(eigval(id, id)) > 0
                QQ(qq, :) = eigvect(:, id) ;
            else
                QQ(qq, :) = -eigvect(:, id) ;
            end
            % Note: big, then small
            otherId = setdiff([1,2], id) ;
            evs(qq, :) = [eigval(id, id) eigval(otherId, otherId)] ;
        end
        bc = barycenter(mesh.u, mesh.f) ;
        bar1 = evs(:, 1) .* QQ(:, 1) ;
        bar1 = cat(2, bar1, evs(:, 1) .* QQ(:, 2)) ;
        quiver(bc(:, 1), bc(:, 2), bar1(:, 1), bar1(:, 2), 1)
        
        
        % Get director on vertices
        nU = mesh.nU ;
        nV = mesh.nV ;
        QV = zeros(length(mesh.u), 2) ;
        evs = QV ;
        bbV = zeros(length(mesh.u), 4) ;
        bbV = F2V * bb ;
        bbV(nU*(nV-1)+1:nU*nV, :) = bbV(1:nU, :) ;
        for qq = 1:length(bbV) 
            b2d = reshape(bbV(qq, :), [2, 2]) ;
            [eigvect, eigval] = eig(b2d) ;
            
            % check that the tensor is still symmetric even though it is
            % averaged onto vertices
            assert(b2d(1, 2) == b2d(2, 1))
            
            % larger eigvect direction
            [~, id] = max(abs([eigval(1, 1), eigval(2, 2)])) ;
            if sign(eigval(id, id)) > 0
                QV(qq, :) = eigvect(:, id) ;
            else
                QV(qq, :) = -eigvect(:, id) ;
            end
            % Note: big, then small
            otherId = setdiff([1,2], id) ;
            if sign(eigval(id, id)) > 0
                evs(qq, :) = [eigval(id, id) eigval(otherId, otherId)] ;
            else
                evs(qq, :) = -[eigval(id, id) eigval(otherId, otherId)] ;
            end
        end
        
        %% Bar c2t, s2t 
        thet = atan2(QV(:, 2), QV(:, 1)) ;
        c2t = evs(:, 1) .* cos(2 * thet) ;
        s2t = evs(:, 1) .* sin(2 * thet) ;
        c2t = imgaussfilt(c2t, 1) ;
        s2t = imgaussfilt(s2t, 1) ;
        uu = reshape(mesh.u(:, 1), [nU, nV]) ;
        vv = reshape(mesh.u(:, 2), [nU, nV]) ;
        
        assert(all(evs(:, 1) > 0))
        mags = vecnorm([c2t, s2t], 2, 2) ;
        thetas = mod(atan2(s2t, c2t), pi)  ;
        m2d = struct() ;
        m2d.f = mesh.f ;
        m2d.v = mesh.u ;
        m2d.v(:, 3) = 0 * mesh.u(:, 1) ;
        opts.clim = rms1d(evs(:, 1)) ;
        opts.axisOff = true ;
        opts.labels = {'', '2nd fundamental form as nematic'} ;
        opts.views = {[0, 0 ], [0, 90 ]};
        opts.cbarlabels = {'$b$ eigenvalue, $||b_1||$', ...
                           '$b$ eigenvalue, $||b_1||$'} ;
        nFieldsOnSurface({mesh, m2d}, ...
            {{mags, thetas}, {mags, thetas}}, opts) ;
        
        sgtitle(['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
                ' ', QS.timeUnits], 'Interpreter', 'latex')
        disp(['Saving ' fnbb])
        saveas(gcf, fnbb)
        close all

        %% Bar X,Y
        % bar1V = evs(:, 1) .* QV(:, 1) ;
        % bar1V = cat(2, bar1V, evs(:, 1) .* QV(:, 2)) ;
        % bQX = reshape(bar1V(:, 1), [nU, nV]) ;
        % bQY = reshape(bar1V(:, 2), [nU, nV]) ;
        % % bQXr = imresize(bQX, [0.5 * nU, NaN]) ;
        % % bQYr = imresize(bQY, [0.5 * nU, NaN]) ;
        % % uu = imresize(reshape(mesh.u(:, 1), [nU, nV]), [0.5 * nU, NaN]) ;
        % % vv = imresize(reshape(mesh.u(:, 2), [nU, nV]), [0.5 * nU, NaN]) ;
        % bQXr = imgaussfilt(bQX, 1) ;
        % bQYr = imgaussfilt(bQY, 1) ;
        % uu = reshape(mesh.u(:, 1), [nU, nV]) ;
        % vv = reshape(mesh.u(:, 2), [nU, nV]) ;
        % % quiver(mesh.u(:, 1), mesh.u(:, 2), bar1V(:, 1), bar1V(:, 2), 1)
        % quiver(uu, vv, bQXr, bQYr, 1) ;
        % hold off;
        
        % stream-like plot
        % flipYQ = bar1V(:, 2) < 0 ;
        % bar1V(flipYQ, :) = -bar1V(flipYQ, :) ;
        % streamline(uu', vv',  bQXr, bQYr, uu(1:10:end), vv(1:10:end)) ;
        % 
        % mags = vecnorm([bQXr(:), bQYr(:)], 2, 2) ;
        % thetas = mod(atan2(bQYr(:), bQXr(:)), pi) ;
        % m2d = struct() ;
        % m2d.f = mesh.f ;
        % m2d.v = mesh.u ;
        % m2d.v(:, 3) = 0 * mesh.u(:, 1) ;
        % opts.clim = 0.2 ;
        % opts.view = [0, 90 ];
        % nFieldsOnSurface({mesh, m2d}, ...
        %     {{mags, thetas}, {mags, thetas}}, opts) ;
    end
    
    % SAVE images of first and second fundamental forms if they dont exist
    files_exist = true ;
    for pp = 1:2
        fn1a = fullfile(sprintf(QS.dir.metric.g_images3d, coordSys, lambda_mesh), ...
             strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
        fn1a = [fn1a '_g_fundForm.png'] ;
        
        fn1b = fullfile(sprintf(QS.dir.metric.g_images2d, coordSys, lambda_mesh), ...
             strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
        fn1b = [fn1b '_g_fundForm.png'] ;
        

        fn2a = fullfile(sprintf(QS.dir.metric.b_images3d, coordSys, lambda_mesh), ...
             strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
        fn2a = [fn2a '_b_fundForm.png'] ;
        fn2b = fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), ...
             strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
        fn2b = [fn2b '_b_fundForm.png'] ;
        
        files_exist = files_exist && exist(fn1a, 'file') && ...
            exist(fn1b, 'file') && exist(fn2a, 'file') && ...
            exist(fn2b, 'file') ;
    end
    if ~files_exist && makeRawMetricComponentFigures
        for pp = 1:2
            opts = struct() ;
            if pp == 1
                opts.clims = {max(abs(gg(:, 1))) * [-1, 1], ...
                    max(abs(gg(:, 2))) * [-1, 1], ...
                    max(abs(gg(:, 3))) * [-1, 1], ...
                    max(abs(gg(:, 4))) * [-1, 1]} ;
            else
                opts.clim = max(abs(gg(:))) * [-1, 1] ;
            end

            % 3d case ---------------------------------------------------------
            labels = {'$\mathbf{g}_{\zeta\zeta}$', ...
                '$\mathbf{g}_{\zeta\phi}$', ...
                '$\mathbf{g}_{\phi\zeta}$', ...
                '$\mathbf{g}_{\phi\phi}$'} ;
            m2view = mesh ;
            m2view.v = mesh.vrs ;
            opts.labels = labels ;
            [~, ~, ~, opts.xyzlims] = QS.getXYZLims() ;
            % gg 3d case ------------------------------------------------------
            close all
            [axs, cbs, meshHandles] = ...
                nFieldsOnSurface({m2view, m2view, m2view, m2view}, ...
                {gg(:, 1), gg(:, 2), gg(:, 3), gg(:, 4)}, opts) ;
            for qq = 1:4
                set(gcf,'CurrentAxes', axs{qq})
                view(0, 0)
                axis off
            end
            sgtitle(['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
                ' ', QS.timeUnits], 'Interpreter', 'latex')
            fn = fullfile(sprintf(QS.dir.metric.g_images3d, coordSys, lambda_mesh), ...
                 strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
            fn = [fn '_g_fundForm.png'] ;
            saveas(gcf, fn)

            % gg 2d case ------------------------------------------------------
            close all
            m2d = mesh ;
            m2d.v = [mesh.u(:, 1) / aspectShear, mesh.u(:, 2), mesh.u(:, 1)*0] ;
            opts.labels = labels ;
            opts = rmfield(opts, 'xyzlims') ;
            % gg
            [axs, cbs, meshHandles] = ...
                nFieldsOnSurface({m2d, m2d, m2d, m2d}, ...
                {gg(:, 1), gg(:, 2), gg(:, 3), gg(:, 4)}, opts) ;
            for qq = 1:4
                set(gcf,'CurrentAxes', axs{qq})
                view(2)
                axis off
            end
            sgtitle(['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
                ' ', QS.timeUnits], 'Interpreter', 'latex')
            fn = fullfile(sprintf(QS.dir.metric.g_images2d, coordSys, lambda_mesh), ...
                 strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
            fn = [fn '_g_fundForm.png'] ;
            saveas(gcf, fn)

            % bb 3d case ------------------------------------------------------
            close all
            opts = struct() ;
            if pp == 1
                opts.clims = {max(abs(bb(:, 1))) * [-1, 1], ...
                    max(abs(bb(:, 2))) * [-1, 1], ...
                    max(abs(bb(:, 3))) * [-1, 1], ...
                    max(abs(bb(:, 4))) * [-1, 1]} ;
            else
                opts.clim = max(abs(bb(:))) * [-1, 1] ;
            end
            labels = {'$\mathbf{b}_{\zeta\zeta}$', ...
                '$\mathbf{b}_{\zeta\phi}$', ...
                '$\mathbf{b}_{\phi\zeta}$', ...
                '$\mathbf{b}_{\phi\phi}$'} ;
            opts.labels = labels ;
            [~, ~, ~, opts.xyzlims] = QS.getXYZLims() ;
            [axs, cbs, meshHandles] = ...
                nFieldsOnSurface({m2view, m2view, m2view, m2view}, ...
                {bb(:, 1), bb(:, 2), bb(:, 3), bb(:, 4)}, opts) ;
            % title
            sgtitle(['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
                ' ', QS.timeUnits], 'Interpreter', 'latex')
            fn = fullfile(sprintf(QS.dir.metric.b_images3d, coordSys, lambda_mesh), ...
                 strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
            fn = [fn '_b_fundForm.png'] ;
            for qq = 1:4
                set(gcf,'CurrentAxes', axs{qq})
                view(0, 0)
                axis off
            end
            saveas(gcf, fn)

            % bb 2d case -------------------------------------------------
            close all
            opts = rmfield(opts, 'xyzlims') ;
            [axs, cbs, meshHandles] = ...
                nFieldsOnSurface({m2d, m2d, m2d, m2d}, ...
                {bb(:, 1), bb(:, 2), bb(:, 3), bb(:, 4)}, opts) ;
            for qq = 1:4
                set(gcf,'CurrentAxes', axs{qq})
                view(2)
                axis off
            end
            % title
            sgtitle(['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
                ' ', QS.timeUnits], 'Interpreter', 'latex')
            fn = fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), ...
                 strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
            fn = [fn '_b_fundForm.png'] ;
            saveas(gcf, fn)
        end
    end
    
end
