function measureBeltramiCoefficient(QS, options)
%measureBeltramiCoefficient(QS, options)
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   t0Pathlines : numeric (default = QS.t0) 
%       timepoint used to define Lagrangian/material frame for mu
%       measurements to be made with respect to
%   overwrite : bool, overwrite previous results on disk
%   climit : limit for caxis of scalar fields Re(mu) and Im(mu)
%
% Returns
% -------
%
% NPMitchell 2020

%% Default options
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

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 't0Pathlines')
    t0Pathlines = options.t0Pathlines ;
else
    disp('Using default t0 for pathlines')
end
if isfield(options, 'climit')
    climit = options.climit ;
end

%% Make sure pathlines are on disk & load them
vXYfn = sprintf(QS.fileName.pathlines_uvprime.vXY, t0Pathlines) ;
v3dfn = sprintf(QS.fileName.pathlines_uvprime.v3d, t0Pathlines) ;
refMeshFn = sprintf(QS.fileName.pathlines_uvprime.refMesh, t0Pathlines) ;
if ~exist(vXYfn, 'file') || ~exist(v3dfn, 'file') || ~exist(refMeshFn, 'file')
    QS.measureUVPrimePathlines(options)
end

% Load pathlines in conformal (u',v') coordinates and their embedding
tmp = load(vXYfn) ;
vP2d = tmp.vertexPathlines ;
tmp = load(v3dfn) ;
vP3d = tmp.v3dPathlines ;
load(refMeshFn, 'refMesh') ;
u_material = [refMesh.u(:, 1), refMesh.u(:, 2)] ;
mu0 = mean(real(bc_metric(refMesh.f, u_material, refMesh.v, 3))) ;

%% 
nTimePoints = length(QS.xp.fileMeta.timePoints) ;
fn = sprintf(QS.fileName.pathlines_uvprime.quasiconformal, t0Pathlines) ;
imDir = sprintf(QS.dir.pathlines_uvprime.quasiconformal, t0Pathlines) ;
if ~exist(fullfile(imDir, '2d'), 'dir')
    mkdir(fullfile(imDir, '2d'))
    mkdir(fullfile(imDir, '3d'))
end

mu_material = zeros(nTimePoints, size(refMesh.f, 1)) ;
mu_instant = zeros(nTimePoints, size(refMesh.f, 1)) ;
arelax = vP2d.affineRelaxFactors(vP2d.tIdx0) ;

todo1 = 1:10:nTimePoints ;
tidx2do = [todo1, setdiff(1:nTimePoints, todo1)] ;

if ~exist(fn, 'file') || overwrite
    for tidx = tidx2do

        %% Set current time
        tp = QS.xp.fileMeta.timePoints(tidx) ;
        imFn2d_material = sprintf(fullfile(imDir, '2d', 'mu2d_material_%06d.png'), tp);
        imFn2d_instant = sprintf(fullfile(imDir, '2d', 'mu2d_instant_%06d.png'), tp);
        imFn3d_material = sprintf(fullfile(imDir, '3d', 'mu3d_material_%06d.png'), tp);
        imFn3d_instant = sprintf(fullfile(imDir, '3d', 'mu3d_instant_%06d.png'), tp);

        QS.setTime(tp)

        %% Measure Beltrami Coefficient
        v3d = zeros(size(vP3d.vX, 2) * size(vP3d.vX, 3), 3) ;
        v3d(:, 1) = reshape(squeeze(vP3d.vXrs(tidx, :, :)), [], 1) ;
        v3d(:, 2) = reshape(squeeze(vP3d.vYrs(tidx, :, :)), [], 1) ;
        v3d(:, 3) = reshape(squeeze(vP3d.vZrs(tidx, :, :)), [], 1) ;
        % 2d mesh
        v2d = zeros(size(vP2d.vX, 2) * size(vP2d.vX, 3), 2) ;
        vU = reshape(vP2d.vU(tidx, :, :), [], 1) ;
        vV = reshape(vP2d.vV(tidx, :, :), [], 1) ;
        v2d(:, 1) = vU ;
        v2d(:, 2) = vV ;
        mu_instant(tidx, :) = bc_metric(refMesh.f, v2d, v3d, 3) ;
        mu_instant(tidx, :) = mu_instant(tidx, :) - nanmean(real(mu_instant(tidx, :))) ;
        mu_material(tidx, :) = bc_metric(refMesh.f, u_material, v3d, 3) - mu0 ;

        %% Plot mu_instant in 3d
        vU = reshape(vP2d.vU(tidx, :, :), [], 1) ;
        close all
        labels = {'$\Re \mu$', '$\Im \mu$'} ;
        options.labels = labels ;
        options.clim = climit ;
        [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
            twoScalarFieldsOnSurface({refMesh.f, v3d}, ...
            real(mu_instant(tidx, :)), imag(mu_instant(tidx, :)), options) ;
        sgtitle(['$\mu($embedding, Lagrangian pullback$)$, $t = $', ...
            sprintf('%03d', tp-t0), ' ', QS.timeUnits], ...
            'interpreter', 'latex') 
        set(gcf,'CurrentAxes', ax1)
        view(0, 0)
        set(gcf,'CurrentAxes', ax2)
        view(0, 0)
        saveas(gcf, imFn3d_instant)
        close all

        %% Plot mu_material in 3d
        options.labels = labels ;
        options.clim = climit ;
        [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
            twoScalarFieldsOnSurface({refMesh.f, v3d}, ...
            real(mu_material(tidx, :)), imag(mu_material(tidx, :)), options) ;
        sgtitle(['$\mu($embedding, material frame$)$, $t = $', ...
            sprintf('%03d', tp-t0), ' ', QS.timeUnits], ...
            'interpreter', 'latex') 
        set(gcf,'CurrentAxes', ax1)
        view(0, 0)
        set(gcf,'CurrentAxes', ax2)
        view(0, 0)
        saveas(gcf, imFn3d_material)
        close all

        %% Save image of mu_instant in 2d
        set(gcf, 'visible', 'off')
        subplot(1, 2, 1)
        [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
            twoScalarFieldsOnSurface({refMesh.f, ...
            [refMesh.u(:, 1), refMesh.u(:, 2), 0 * refMesh.u(:, 2)]}, ...
            real(mu_instant(tidx, :)), imag(mu_instant(tidx, :)), options) ;
        sgtitle(['$\mu($embedding, Lagrangian pullback$)$, $t = $', ...
            sprintf('%03d', tp-t0), ' ', QS.timeUnits], ...
            'interpreter', 'latex') 
        set(gcf,'CurrentAxes', ax1)
        view(2)
        set(gcf,'CurrentAxes', ax2)
        view(2)
        saveas(gcf, imFn2d_instant)
        close all

        %% Save image of mu_material in 2d 
        set(gcf, 'visible', 'off')
        subplot(1, 2, 1)
        [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
            twoScalarFieldsOnSurface({refMesh.f, ...
            [refMesh.u(:, 1), refMesh.u(:, 2), 0 * refMesh.u(:, 2)]}, ...
            real(mu_material(tidx, :)), imag(mu_material(tidx, :)), options) ;
        sgtitle(['$\mu($embedding, material frame$)$, $t = $', ...
            sprintf('%03d', tp-t0), ' ', QS.timeUnits], ...
            'interpreter', 'latex') ;
        set(gcf,'CurrentAxes', ax1)
        view(2)
        set(gcf,'CurrentAxes', ax2)
        view(2)
        saveas(gcf, imFn2d_material)
        close all

    end

    % Save data as output
    save(fn, 'mu_instant', 'mu_material') 
end


