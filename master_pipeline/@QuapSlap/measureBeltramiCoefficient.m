function measureBeltramiCoefficient(QS, options)
%measureBeltramiCoefficient(QS, options)
% 
% Previous version measured instantaneous mu as well, which was the
% Beltrami coefficient to the current (nearly) conformal pullback.
%
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   t0Pathlines : numeric (default = QS.t0) 
%       timepoint used to define Lagrangian/material frame for mu
%       measurements to be made with respect to
%   coordSys : str specifier ('ricci', 'uvprime' or 'sp')
%       coordinate system in which pathlines are computed
%   overwrite : bool, overwrite previous results on disk
%   climit : limit for caxis of scalar fields Re(mu) and Im(mu)
%
% Returns
% -------
%
% NPMitchell 2020

%% Default options
save_ims = true ;
coordSys = 'ricci' ;    % or 'uvprime'
maxIter = 200 ;         % only used if coordSys == 'ricci', #ricci iterations
if isempty(QS.t0)
    t0 = QS.t0set() ;
else
    t0 = QS.t0 ;
end
t0Pathlines = t0 ;
overwrite = false ;
if nargin < 2 
    options = struct() ;
end
climit = 1 ;
[~, ~, ~, xyzlim ] = QS.getXYZLims() ; 

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'save_ims')
    save_ims = options.save_ims ;
end
if isfield(options, 't0Pathlines')
    t0Pathlines = options.t0Pathlines ;
else
    disp('Using default t0 for pathlines')
end
if isfield(options, 'coordSys')
    coordSys = options.coordSys ;
end
if isfield(options, 'climit')
    climit = options.climit ;
end

%% Make sure pathlines are on disk & load them
if contains(coordSys, 'ricci')
    vXYfn = sprintf(QS.fileName.pathlines.vXY, t0Pathlines) ;
    v3dfn = sprintf(QS.fileName.pathlines.v3d, t0Pathlines) ;
    refMeshFn = sprintf(QS.fileName.pathlines.refMesh, t0Pathlines) ;
    if ~exist(vXYfn, 'file') || ~exist(v3dfn, 'file') || ~exist(refMeshFn, 'file')
        QS.measurePullbackPathlines(options)
    end
elseif contains(coordSys, 'uvprime')
    vXYfn = sprintf(QS.fileName.pathlines_uvprime.vXY, t0Pathlines) ;
    v3dfn = sprintf(QS.fileName.pathlines_uvprime.v3d, t0Pathlines) ;
    refMeshFn = sprintf(QS.fileName.pathlines_uvprime.refMesh, t0Pathlines) ;
    if ~exist(vXYfn, 'file') || ~exist(v3dfn, 'file') || ~exist(refMeshFn, 'file')
        QS.measureUVPrimePathlines(options)
    end
else
    vXYfn = sprintf(QS.fileName.pathlines.vXY, t0Pathlines) ;
    v3dfn = sprintf(QS.fileName.pathlines.v3d, t0Pathlines) ;
    refMeshFn = sprintf(QS.fileName.pathlines.refMesh, t0Pathlines) ;
end

% Load pathlines in conformal/(u',v')/other coordinates and their embedding
tmp = load(vXYfn) ;
vP2d = tmp.vertexPathlines ;
tmp = load(v3dfn) ;
vP3d = tmp.v3dPathlines ;
load(refMeshFn, 'refMesh') ;
if contains(coordSys, 'ricci')
    u_material = refMesh.u_ricci ;
    mu0 = 0 + 1j * 0 ;
else
    u_material = [refMesh.u(:, 1), refMesh.u(:, 2)] ;
    mu0 = mean(real(bc_metric(refMesh.f, u_material, refMesh.v, 3))) ;
end
try
    assert(abs(mu0) < 1e-4)
catch
    disp('refMesh is not particularly conformal...')
end

%% 
nTimePoints = length(QS.xp.fileMeta.timePoints) ;

if contains(coordSys, 'ricci')
    fn = sprintf(QS.fileName.pathlines.quasiconformal, t0Pathlines) ;
    imDir = sprintf(QS.dir.pathlines.quasiconformal, t0Pathlines) ;
elseif contains(coordSys, 'uvprime')
    fn = sprintf(QS.fileName.pathlines_uvprime.quasiconformal, t0Pathlines) ;
    imDir = sprintf(QS.dir.pathlines_uvprime.quasiconformal, t0Pathlines) ;
end
if ~exist(fullfile(imDir, '2d'), 'dir')
    mkdir(fullfile(imDir, '2d'))
    mkdir(fullfile(imDir, '3d'))
end

mu_material = zeros(nTimePoints, size(refMesh.f, 1)) ;
mu_material_filtered = zeros(nTimePoints, refMesh.nU * refMesh.nV) ;
% arelax = vP2d.affineRelaxFactors(vP2d.tIdx0) ;

todo1 = [180, 1:10:nTimePoints] ;
tidx2do = [todo1, setdiff(1:nTimePoints, todo1)] ;

first = true ;
if ~exist(fn, 'file') || overwrite
    for tidx = tidx2do
        disp(['tidx = ' num2str(tidx)])
        
        %% Set current time
        tp = QS.xp.fileMeta.timePoints(tidx) ;
        imFn2d_material = sprintf(fullfile(imDir, '2d', 'mu2d_material_%06d.png'), tp);
        imFn3d_material = sprintf(fullfile(imDir, '3d', 'mu3d_material_%06d.png'), tp);
        
        QS.setTime(tp)

        %% Measure Beltrami Coefficient
        v3d = zeros(size(vP3d.vX, 2) * size(vP3d.vX, 3), 3) ;
        v3d(:, 1) = reshape(squeeze(vP3d.vXrs(tidx, :, :)), [], 1) ;
        v3d(:, 2) = reshape(squeeze(vP3d.vYrs(tidx, :, :)), [], 1) ;
        v3d(:, 3) = reshape(squeeze(vP3d.vZrs(tidx, :, :)), [], 1) ;
        
        %% Material Beltrami coefficient
        mu_material(tidx, :) = bc_metric(refMesh.f, u_material, v3d, 3) - mu0 ;

        %% Mode filter mu_material 
        refMesh2glue = refMesh ;
        refMesh2glue.v = v3d ;
        glueMesh = glueCylinderCutMeshSeam(refMesh2glue) ;
        nU = refMesh.nU ;
        nV = refMesh.nV ;
        
        [V2F, F2V] = meshAveragingOperators(glueMesh.f, glueMesh.v) ;
        mu_material_vtx = F2V * squeeze(mu_material(tidx, :))' ;
        mu_material_vtx = reshape(mu_material_vtx, [refMesh.nU, refMesh.nV-1]) ;
        
        %% Check filtered image
        if first 
            nmodes2do = 1:5 ;
            clf
            nrows = nmodes2do(end) + 1 ;
            for nmodes_ii = nmodes2do
                options.widthX = 3 ;
                options.nmodesY = nmodes_ii ;
                muMVfilt_re = modeFilterQuasi1D(real(mu_material_vtx), options) ;
                
                subplot(nrows, 3, 1 + 3*(nmodes_ii-1))
                imagesc(1:nU, 1:nV-1, real(mu_material_vtx)') ;
                axis equal; axis tight; axis off
                caxis([-1,1])
                colormap blueblackred
                if nmodes_ii == 1
                    title('raw $\mu$', 'interpreter', 'latex')
                elseif nmodes_ii == nmodes2do(end)
                    pos = get(gca, 'position') ;
                    cb = colorbar('location', 'southOutside');
                    set(gca, 'position', pos)
                end
                subplot(nrows, 3, 2 + 3*(nmodes_ii-1)) 
                imagesc(1:nU, 1:nV-1, muMVfilt_re')
                caxis([-1,1])
                colormap blueblackred
                axis equal; axis tight; axis off
                if nmodes_ii == 1
                    title('filtered $\mu$', 'interpreter', 'latex')
                elseif nmodes_ii == nmodes2do(end)
                    pos = get(gca, 'position') ;
                    cb = colorbar('location', 'southOutside');
                    set(gca, 'position', pos)
                end
                subplot(nrows, 3, 3 + 3*(nmodes_ii-1)) 
                imagesc(1:nU, 1:nV-1, muMVfilt_re'-real(mu_material_vtx)') ;
                caxis([-1,1])
                colormap blueblackred
                axis equal; axis tight; axis off
                if nmodes_ii == 1
                    title('filtered - raw', 'interpreter', 'latex')
                elseif nmodes_ii == nmodes2do(end)
                    pos = get(gca, 'position') ;
                    cb = colorbar('location', 'southOutside');
                    set(gca, 'position', pos)
                    % set(cb, 'position', get(cb, 'position'));
                end
            end
            sgtitle(['Lowest ', num2str(options.nmodesY), ' modes, ', ...
                '$\sigma= $' num2str(options.widthX)], 'interpreter', 'latex')  
            fnfilter = fullfile(imDir, sprintf('filt_test.png')) ;
            disp(['saving ' fnfilter])
            saveas(gcf, fnfilter)
        end
        
        % Do filtered for real
        filterOptions.widthX = 3 ;
        filterOptions.nmodesY = 5 ;
        filterOptions.preview = false ;
        muMVfilt_re = modeFilterQuasi1D(real(mu_material_vtx), filterOptions) ;
        muMVfilt_im = modeFilterQuasi1D(imag(mu_material_vtx), filterOptions) ;
        % Reshape into [nU, nV] by repeating seam
        muMVfilt_re(:, nV) = muMVfilt_re(:, 1) ;
        muMVfilt_im(:, nV) = muMVfilt_im(:, 1) ;
        
        mu_material_filtered(tidx, :) = muMVfilt_re(:) + 1j * muMVfilt_im(:) ;
        
        if save_ims
        
            close all
            labels = {'$\Re \mu$', '$\Im \mu$'} ;
            
            %% Plot mu_material in 3d
            options.labels = labels ;
            options.clim = climit ;
            [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
                twoScalarFieldsOnSurface({refMesh.f, v3d}, ...
                real(mu_material_filtered(tidx, :)), ...
                imag(mu_material_filtered(tidx, :)), options) ;
            sgtitle(['$\mu($embedding, material frame$)$, $t = $', ...
                sprintf('%03d', tp-t0), ' ', QS.timeUnits], ...
                'interpreter', 'latex') 
            set(gcf,'CurrentAxes', ax1)
            view(0, 0)
            xlim(xyzlim(1, :))
            ylim(xyzlim(2, :))
            zlim(xyzlim(3, :))
            axis off
            set(gcf,'CurrentAxes', ax2)
            view(0, 0)
            xlim(xyzlim(1, :))
            ylim(xyzlim(2, :))
            zlim(xyzlim(3, :))
            axis off
            disp(['saving ' imFn3d_material])
            saveas(gcf, imFn3d_material)
            close all

            %% Save image of mu_material in 2d 
            set(gcf, 'visible', 'off')
            [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
                twoScalarFieldsOnSurface({refMesh.f, ...
                [refMesh.u(:, 1) / max(refMesh.u(:, 1)), ...
                refMesh.u(:, 2), 0 * refMesh.u(:, 2)]}, ...
                real(mu_material_filtered(tidx, :)), ...
                imag(mu_material_filtered(tidx, :)), options) ;
            sgtitle(['$\mu($embedding, material frame$)$, $t = $', ...
                sprintf('%03d', tp-t0), ' ', QS.timeUnits], ...
                'interpreter', 'latex') ;
            set(gcf,'CurrentAxes', ax1)
            view(2)
            set(gcf,'CurrentAxes', ax2)
            view(2)
            disp(['saving ' imFn2d_material])
            saveas(gcf, imFn2d_material)
            close all
        end
        
        % No longer the first pass
        first = false ;
    end

    % Save data as output
    disp(['saving ' fn])
    save(fn, 'mu_material', 'mu_material_filtered', 'filterOptions') 
end
disp('done with measureBeltramiCoefficient')


