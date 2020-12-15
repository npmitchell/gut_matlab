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
coordSys = 'spsm' ;

% Unpack options
if isfield(options, 'coordSys')
    coordSys = options.coordSys ;
end

for tidx = 1:length(QS.xp.fileMeta.timePoints)
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    
    strClims = {'climVariable', 'climUniform'} ;
    outdirs = {sprintf(QS.dir.metric.data, coordSys), ...
        fullfile(sprintf(QS.dir.metric.g_images3d, coordSys), strClims{1}), ...
        fullfile(sprintf(QS.dir.metric.b_images3d, coordSys), strClims{1}), ...
        fullfile(sprintf(QS.dir.metric.g_images3d, coordSys), strClims{2}), ...
        fullfile(sprintf(QS.dir.metric.b_images3d, coordSys), strClims{2}), ...
        fullfile(sprintf(QS.dir.metric.g_images2d, coordSys), strClims{1}), ...
        fullfile(sprintf(QS.dir.metric.b_images2d, coordSys), strClims{1}), ...
        fullfile(sprintf(QS.dir.metric.g_images2d, coordSys), strClims{2}), ...
        fullfile(sprintf(QS.dir.metric.b_images2d, coordSys), strClims{2})   };
    for qq = 1:length(outdirs)
        if ~exist(outdirs{qq}, 'dir')
            mkdir(outdirs{qq}) ;
        end
    end
    
    outfn = sprintf(QS.fullFileBase.metric, coordSys, tp) ;
    
    % Load mesh
    if strcmpi(coordSys, 'spsm')
        mesh = load(sprintf(QS.fullFileBase.spcutMeshSm, tp), 'spcutMeshSm') ;
        mesh = mesh.spcutMeshSm ;
    else
        error('handle this coordSys here')
    end
    
    % Rescale u dimension for mu relaxation
    uv = mesh.u ;
    uv(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) ;
    mesh.readme = 'mu is computed after rescaling refMesh.u(:, 1) to range from 0 to 1' ;
    mesh.mu = bc_metric(mesh.f, uv, mesh.v, 3) ;
    
    % Project to APDV coords if not already done
    try 
        mesh.vrs;
    catch
        mesh.vrs = QS.xyz2APDV(mesh.v) ;
    end
    
    % Rescale and recompute mu
    aspectShear = (1 - mean(real(mesh.mu))) / (1 + mean(real(mesh.mu))) ;
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
    save(outfn, 'gg', 'bb', 'mesh')
    
    % SAVE images of first and second fundamental forms
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
        fn = fullfile(sprintf(QS.dir.metric.g_images3d, coordSys), ...
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
        fn = fullfile(sprintf(QS.dir.metric.g_images2d, coordSys), ...
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
        fn = fullfile(sprintf(QS.dir.metric.b_images3d, coordSys), ...
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
        fn = fullfile(sprintf(QS.dir.metric.b_images2d, coordSys), ...
             strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
        fn = [fn '_b_fundForm.png'] ;
        saveas(gcf, fn)
    end
end
