function [ricciMesh, ricciMu] = generateRicciMeshTimePoint(QS, tp, options)
% generateRicciMeshTimePoint(QS, tp, options)
% 
% Parameters
% ----------
% QS : quapSlap class instance
% tp : int (default=QS.t0)
%   timestamp of timepoint for which to generate Ricci flow mesh
% options : optional struct with optional fields
%   maxIter : int (default=100)
%   radiusTolerance : float (default=0.01)
%       maximum allowed fractional deviation of Ricci flow solution inner
%       and outer radius from a true circle with fixed radius (variable
%       inner radius, fixed outer radius of 1)
%   save_ims : bool (default = true)
%       save images of the ricci flow solution
%
% Returns
% -------
%
% Saves to disk
% -------------
%
%
% NPMitchell 2020

%% Default parameters
t0 = QS.t0set() ;
radiusTolerance = 0.01 ;
maxIter = 200 ;
save_ims = true ;

%% Unpack parameters
if nargin < 2
    tp = t0 ;
elseif isempty(tp)
    tp = t0 ;
end
    
%% Load cutmesh to flow
cutMesh = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
cutMesh = cutMesh.spcutMeshSmRS ;
glueMesh = glueCylinderCutMeshSeam(cutMesh) ;
nU = cutMesh.nU ;
nV = cutMesh.nV ;
imDir = fullfile(QS.dir.ricciMesh, 'images', sprintf('%04diter', maxIter)) ;
if ~exist(imDir, 'dir') && save_ims
    mkdir(imDir)
end

%% Generate conformal parameterization in the unit disk
% Load or compute ricci flow solution
ricciFn = sprintf(QS.fullFileBase.ricciSolution, maxIter, tp) ;
try 
    load(ricciFn, 'U')
catch
    [~, U, ~] = DiscreteRicciFlow.EuclideanRicciFlow(glueMesh.f, glueMesh.v, ...
        'BoundaryType', 'Fixed', 'BoundaryShape', 'Circles', ...
        'MaxCircIter', maxIter);
    % [labels, dbonds, topStructTools] = labelRectilinearMeshBonds(cutMesh) ;
    % Let maxIter >= 50
    save(ricciFn, 'U') ;
end

% Plot Ricci result in annulus
outfn = fullfile(imDir, sprintf('%06d_RicciSolution.png', tp)) ;
if ~exist(outfn, 'file') && save_ims
    clf
    triplot(triangulation(glueMesh.f, U), 'color', 'k')
    axis equal
    axis tight
    title('Ricci flow solution', 'interpreter', 'latex')
    xlabel('$\tilde{x}$', 'interpreter', 'latex')
    ylabel('$\tilde{y}$', 'interpreter', 'latex')
    saveas(gcf, outfn)
end


%% Beltrami
mu_annulus = bc_metric(glueMesh.f, U, glueMesh.v, 3) ;

%% plot beltrami for Ricci flow
mesh2d = triangulation(glueMesh.f, [U(:, 1), U(:, 2)]) ; %  ; , 0*U(:, 1)]) ;
options.view = [0, 90] ;
options.axisOff = true ;
options.visible = true ;
options.labels = {'$\Re\mu$', '$\Im\mu$', '$|\mu|$', ...
    '$\Re\mu$', '$\Im\mu$', '$|\mu|$'} ;
outfn = fullfile(imDir, sprintf('%06d_RicciFlowBeltrami.png', tp)) ;
if ~exist(outfn, 'file') && save_ims
    clf
    nFieldsOnSurface({mesh2d, mesh2d, mesh2d, ...
        glueMesh, glueMesh, glueMesh}, ...
        {real(mu_annulus), imag(mu_annulus), abs(mu_annulus), ...
        real(mu_annulus), imag(mu_annulus), abs(mu_annulus)}, ...
        options)
    sgtitle('Beltrami coefficients from Ricci flow', ...
        'interpreter', 'latex')
    saveas(gcf, outfn)
end

% %% Compare to beltrami from UVprime cutmesh
% uvMesh = load('./uvpcutMesh_000160.mat') ;
% uvMesh = uvMesh.uvpcutMesh.resampled ;
% m2duv = struct() ;
% m2duv.f = uvMesh.f ;
% m2duv.v = uvMesh.u ;
% mu_uv = bc_metric(uvMesh.f, uvMesh.u, uvMesh.v, 3) ;
% 
% % plot beltrami
% options.view = [0, 90] ;
% options.axisOff = true ;
% options.visible = true ;
% options.labels = {'$\Re\mu$', '$\Im\mu$', '$|\mu|$', ...
%     '$\Re\mu$', '$\Im\mu$', '$|\mu|$'} ;
% outfn = fullfile(imDir, sprintf('%06d_DirichletBeltrami.png', tp)) ;
% if ~exist(outfn, 'file') && save_ims
%     clf
%     nFieldsOnSurface({m2duv, m2duv, m2duv, ...
%         glueMesh, glueMesh, glueMesh}, ...
%         {real(mu_uv), imag(mu_uv), abs(mu_uv), ...
%         real(mu_uv), imag(mu_uv), abs(mu_uv)}, ...
%         options)
%     sgtitle('Beltrami coefficients from Dirichlet energy minimization', ...
%         'interpreter', 'latex')
%     saveas(gcf, outfn)
% end

%% Plot just the bare triangulation
outfn = fullfile(imDir, sprintf('%06d_RicciFlowSolution.png', tp)) ;
if ~exist(outfn, 'file') && save_ims
    clf
    triplot(triangulation(glueMesh.f, U), 'color', 'k')
    axis equal; axis tight
    hold on
    scatter(0,0,50, 'filled', 'b')
    title('off-center inner boundary', 'interpreter', 'latex')
    xlabel('$\tilde{x}$', 'interpreter', 'latex')
    ylabel('$\tilde{y}$', 'interpreter', 'latex')
    saveas(gcf, outfn)
end


%% Push barycenter of inner hole to origin via Mobius transformation
boundaries = DiscreteRicciFlow.compute_boundaries(glueMesh.f) ;
% which boundary is the inner one? Find which one has shorter circumference
LL = [0, 0] ;
for qq = 1:length(boundaries)
    % make periodic 
    LL(qq) = sum(sqrt(sum((U(boundaries{qq},:)- ...
        circshift(U(boundaries{qq},:), 1, 1)).^2, 2))) ;
end
% which is smaller?
[~, innerID] = min(LL) ;
inner = boundaries{innerID} ;
outer = boundaries{mod(innerID+1, 2)} ;

% barycenter of inner -- THIS IS BIASED
% baryc = mean(U(inner, :), 1) ;
% baryc = complex(baryc(1), baryc(2)) ;

% Take center (not barycenter) of circle
% note: https://www.mathworks.com/matlabcentral/fileexchange/5557-circle-fit
[xc,yc, innerR] = circfit(U(inner, 1), U(inner, 2)); 
baryc = xc + 1j * yc ;

% Covert U to complex
zz = complex(U(:, 1), U(:, 2)) ;

% Mobius transform to place center of annulus at origin
zz = (zz - baryc) ./ (1 - conj(baryc) .* zz) ;
UU = [real(zz), imag(zz) ] ;

% inspect centered mesh
clf
triplot(triangulation(glueMesh.f, UU))
hold on;
plot(UU(boundaries{1}, 1), UU(boundaries{1}, 2), 'co-')
plot(UU(boundaries{2}, 1), UU(boundaries{2}, 2), 'co-')
scatter(0,0, 'r', 'filled')
axis equal
xlim([-2*innerR, 2*innerR])

%% Enforce circularity in inner Boundary
% push inner boundary onto circle with exact radius
phaseInner = atan2(UU(inner, 2), UU(inner, 1)) ; 
% check that this is a minor correction
radii = vecnorm(UU(inner, :), 2, 2) ;
disp(['Correcting radial coordinate by a maximum of ' ...
    num2str(max(abs(radii - innerR) / innerR)*100) '%'])
assert(max(abs(radii - innerR) / innerR) < radiusTolerance)
UU(inner, 1) = innerR * cos(phaseInner) ;
UU(inner, 2) = innerR * sin(phaseInner) ;

% Enforce circularity in outer boundary ==> radius=1
phaseOuter = atan2(UU(outer, 2), UU(outer, 1)) ; 
% check that this is a minor correction
radii = vecnorm(UU(outer, :), 2, 2) ;
disp(['Correcting radial coordinate by a maximum of ' ...
    num2str(max(abs(radii - 1))*100) '%'])
assert(max(abs(radii - 1)) < radiusTolerance)
UU(outer, 1) = cos(phaseOuter) ;
UU(outer, 2) = sin(phaseOuter) ;

% inspect the adjusted mesh
outfn1 = fullfile(imDir, sprintf('%06d_ricci_InnerCorrection.png', tp)) ;
outfn2 = fullfile(imDir, sprintf('%06d_ricci_OuterCorrection.png', tp)) ;
if ~exist(outfn1, 'file')
    triplot(triangulation(glueMesh.f, UU), 'color', 'k')
    hold on;
    plot(UU(boundaries{1}, 1), UU(boundaries{1}, 2), 'b.-')
    plot(UU(boundaries{2}, 1), UU(boundaries{2}, 2), 'b.-')
    scatter(0,0, 'r', 'filled')
    axis equal;
    xlim([-2*innerR, 2*innerR])
    xlabel('$\tilde{x}$', 'interpreter', 'latex')
    ylabel('$\tilde{y}$', 'interpreter', 'latex')
    sgtitle('Ricci flow result with circularity correction', ...
        'interpreter', 'latex')
    saveas(gcf, outfn1)
    ylim([-1,1])
    axis equal; axis tight
    saveas(gcf, outfn2)
end
clf

%% Take log
cutMesh.pathPairs(:, 1)
zz = complex(UU(:, 1), UU(:, 2)) ;
rho = real(log(zz)) ;
phi = imag(log(zz)) ;
triplot(triangulation(glueMesh.f, [rho, phi])) ;
axis equal
rhoM = reshape(rho, [nU, nV-1]) ;
phiM = reshape(phi, [nU, nV-1]) ;

% Check direction
grX = diff(rhoM, 1, 1) ;
grY = diff(phiM, 1, 2) ;
outfn = fullfile(imDir, sprintf('%06d_DrhoDphi_rectified.png', tp)) ;
if ~exist(outfn, 'file') && save_ims
    clf
    subplot(2, 1, 1)
    scatter(reshape(rhoM(1:nU-1, :), [], 1), ...
        reshape(phiM(1:nU-1,:), [], 1), 10, grX(:))
    title('$\Delta \tilde{u}$', 'interpreter', 'latex')
    xlabel('$\tilde{u}$', 'interpreter', 'latex')
    xlabel('$\tilde{v}$', 'interpreter', 'latex')
    caxis([-0.05, 0.05])
    axis equal ; axis tight 
    colormap blueblackred
    cb = colorbar ;
    ylabel(cb, '$\Delta \tilde{u}$', 'interpreter', 'latex')
    subplot(2, 1, 2)
    scatter(reshape(rhoM(:, 1:end-1), [], 1), ...
        reshape(phiM(:,1:end-1), [], 1), 10, grY(:))
    title('$\Delta \tilde{v}$', 'interpreter', 'latex')
    xlabel('$\tilde{u}$', 'interpreter', 'latex')
    xlabel('$\tilde{v}$', 'interpreter', 'latex')
    caxis([-0.05, 0.05])
    axis equal ; axis tight 
    colormap blueblackred
    cb = colorbar() ;
    sgtitle('initial pullback orientation')
    ylabel(cb, '$\Delta \tilde{v}$', 'interpreter', 'latex')
    saveas(gcf, outfn)
end

phi0 = phi ;
% Push y range to approx (0, 2pi), fliping if needed
if all(mean(sign(grY), 2) < 0)
    disp('Flipping angles phi')
    % Push phi to approx (0, 2pi)
    phi = pi - phi ;
elseif all(mean(sign(grY), 2) > 0)
    disp('phi angles properly ordered')
    % Push phi to approx (0, 2pi)
    phi = phi + pi ;
else
    error('Could not determine direction of phi in meshgrid')
end
% Determine whether to flip in x direction
if all(all(sign(grX) < 0))
    rho = -rho ;
end
minrho = min(rho) ;
if abs(minrho) > 0 
    disp('Translating rho to origin')
    rho = rho - minrho ;
end

% Remake grids of rhoM and phiM
rhoM = reshape(rho, [nU, nV-1]) ;
phiM = reshape(phi, [nU, nV-1]) ;

% Check direction
grX = diff(rhoM, 1, 1) ;
grY = diff(phiM, 1, 2) ;

outfn = fullfile(imDir, sprintf('%06d_DrhoDphi_flipped.png', tp)) ;
if ~exist(outfn, 'file') && save_ims
    clf
    subplot(2, 1, 1)
    scatter(reshape(rhoM(1:nU-1, :), [], 1), ...
        reshape(phiM(1:nU-1,:), [], 1), 10, grX(:))
    title('$\Delta \tilde{u}$', 'interpreter', 'latex')
    xlabel('$\tilde{u}$', 'interpreter', 'latex')
    ylabel('$\tilde{v}$', 'interpreter', 'latex')
    caxis([-0.05, 0.05])
    axis equal ; axis tight 
    colormap blueblackred
    cb = colorbar ;
    ylabel(cb, '$\Delta \tilde{u}$', 'interpreter', 'latex')
    subplot(2, 1, 2)
    scatter(reshape(rhoM(:, 1:end-1), [], 1), ...
        reshape(phiM(:,1:end-1), [], 1), 10, grY(:))
    title('$\Delta \tilde{v}$', 'interpreter', 'latex')
    xlabel('$\tilde{u}$', 'interpreter', 'latex')
    ylabel('$\tilde{v}$', 'interpreter', 'latex')
    caxis([-0.05, 0.05])
    axis equal ; axis tight 
    colormap blueblackred
    cb = colorbar() ;
    ylabel(cb, '$\Delta \tilde{v}$', 'interpreter', 'latex')
    sgtitle('reorienting pullback')
    saveas(gcf, outfn)
end

%% plot initial cutpath
outfn = fullfile(imDir, sprintf('%06d_phiOrderInitial.png', tp)) ;
if ~exist(outfn, 'file') && save_ims
    clf
    triplot(triangulation(glueMesh.f, [rho, phi]), 'color', 'k') ;
    hold on;
    scatter(rho(1:nU), phi(1:nU), 'filled', 'c')
    axis equal; axis tight
    title('initial cut path', 'interpreter', 'latex')
    xlabel('$\tilde{u}$', 'interpreter', 'latex')
    ylabel('$\tilde{v}$', 'interpreter', 'latex')
    saveas(gcf, outfn)
end

%% Push all 1:nU vertices to near phi = 0 (cutPath leveling in PB space)
% Make first phi in each row the lowest phi value
clf 
phi_recut = phi ;
for qq = 1:nU
    % Consider this column of the rectilinear pullback structure
    phis = phi(qq:nU:nU*(nV-1)) ;
    rhos = rho(qq:nU:nU*(nV-1)) ;
    phis(phis < phis(1)) = phis(phis < phis(1)) + 2*pi ;

    % Correct all first phis to be in the same register 
    % Note: this means they might not be on the same branch cut!
    if qq == 1
        % For first column, make phi(1) zero
        phi0 = phis(1) ;
        phis = phis - phi0 ;
    else
        % For subsequent columns, find the right branch cut that connects
        % most closely to the adjacent column without allowing
        % a flip in the normal of the mesh. Note that this assumes nothing
        % totally insane is happening to the geometry of the mesh across
        % mesh sampling distance.
        phis = phis - phi0 ;

        % I think it is impossible to be off by more than one branch cut --
        % ie by more than 2pi, so lets find the best match of those
        % possibilities in adjacent branch cuts
        possible = [-2*pi, 0, 2*pi] ;
        [min_dphi, ind_min] = min(abs(phis(1) + possible - prevphi1)) ;
        phis = phis + possible(ind_min) ;
        disp(['Translating phi values by '...
            num2str(round(possible(ind_min)/pi)) 'pi'])

        % todo: perform check on the sign of the cross product of a face
        % including this vertex to check for no mesh intersections
        plot(rhos, phis, '.'); 
        hold on ;
    end

    phi_recut(qq:nU:nU*(nV-1)) = phis ;
    prevphi1 = phis(1) ;
end

%% Minimize offset to all phi values based on previous mesh vertices 
% in pullback space ---> ignore this for example script
tidx = QS.xp.tIdx(tp) ;
tidx0 = QS.xp.tIdx(QS.t0set()) ;
if tidx ~= tidx0 
    if tidx > tidx0
        % load previous timepoint phi values
        prevTP = QS.xp.fileMeta.timePoints(tidx + 1) ;
        phi2match = load(sprintf(QS.fullFileBase.uvtMesh, prevTP)) ;
    else
        assert(tidx < tidx0)
        % load next timepoint phi values
        prevTP = QS.xp.fileMeta.timePoints(tidx - 1) ;
        phi2match = load(sprintf(QS.fullFileBase.uvtMesh, prevTP)) ;
    end
    overall_offset = mean(phi_recut(:) - phi2match(:)) ; 
    phi_recut = phi_recut - overall_offset ;
end

%% plot final cutpath
outfn = fullfile(imDir, sprintf('%06d_phiOrderFinal.png', tp)) ;
if ~exist(outfn, 'file') && save_ims
    clf
    triplot(triangulation(glueMesh.f, [rho, phi_recut]), 'color', 'k') ;
    hold on;
    scatter(rho(1:nU), phi_recut(1:nU), 'filled','c')
    axis equal; axis tight
    title('final cut path', 'interpreter', 'latex')
    xlabel('$\tilde{u}$', 'interpreter', 'latex')
    ylabel('$\tilde{v}$', 'interpreter', 'latex')
    saveas(gcf, outfn)
end

%% Unwrap from annulus into rectangle and save in struct
ricciMesh = struct() ;
ricciMesh.annulus = struct('f', glueMesh.f, 'u', UU, 'v', glueMesh.v, ...
    'nU', nU, 'nV', nV) ;
riccicutMesh = struct() ;
riccicutMesh = ricciMesh.annulus ;
riccicutMesh.u = [rho, phi_recut] ;
opts = struct('vmax', 2 * pi, 'ignoreRectangularConstraint', true) ;
riccicutMesh = cutRectilinearCylMesh(riccicutMesh, opts) ;
ricciMesh.rectangle = struct('f', riccicutMesh.f, 'u', ...
    riccicutMesh.u, 'v', riccicutMesh.v, ...
    'nU', nU, 'nV', nV) ;

% Save ricciMesh
% note: ricciMesh has fields annulus and rectangle
ricciMeshFn = sprintf(QS.fullFileBase.ricciMesh, maxIter, tp) ;
disp(['Saving ricciMesh to ' ricciMeshFn])
save(ricciMeshFn, 'ricciMesh')

%% Preview final ordering
outfn = fullfile(imDir, sprintf('%06d_phiOrderFinal.png', tp)) ;
if ~exist(outfn, 'file') && save_ims
    facecolors = (1:length(ricciMesh.rectangle.f)) ;
    facecolors = mod(facecolors, 2*nV-2) ;
    facecolors(facecolors==0) = 2*nV - 2 ;
    trisurf2d(ricciMesh.rectangle.f, ricciMesh.rectangle.u, 'FaceVertexCData', facecolors(:), 'Facealpha',0.5)

    saveas(gcf, outfn)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recompute mu in rectangular space 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It should be conformal since this is a conformal map of a conformal map
% into the disk
mu_rectangle = bc_metric(riccicutMesh.f, riccicutMesh.u, riccicutMesh.v, 3) ;

%% plot beltrami for Ricci flow
mesh2d = triangulation(riccicutMesh.f, riccicutMesh.u) ; %  ; , 0*U(:, 1)]) ;
options.view = [0, 90] ;
options.axisOff = true ;
options.visible = true ;
options.labels = {'$\Re\mu$', '$\Im\mu$', '$|\mu|$', ...
    '$\Re\mu$', '$\Im\mu$', '$|\mu|$'} ;

outfn = fullfile(imDir, sprintf('%06d_RicciFlowBeltrami_rectangle.png', tp)) ;
if ~exist(outfn, 'file') && save_ims
    clf
    nFieldsOnSurface({mesh2d, mesh2d, mesh2d, ...
        glueMesh, glueMesh, glueMesh}, ...
        {real(mu_rectangle), imag(mu_rectangle), abs(mu_rectangle), ...
        real(mu_rectangle), imag(mu_rectangle), abs(mu_rectangle)}, ...
        options)
    sgtitle('Beltrami coefficients from $\log$ of Ricci flow', ...
        'interpreter', 'latex')
    saveas(gcf, outfn)
end

% %% Think about variation ds, dv wrt conformal coordinates
% [labels, dbonds, topStructTools] = labelRectilinearMeshBonds(cutMesh) ;
% tmp = vecnorm(dbonds.realSpace.v, 2, 2) ;
% trisurf(triangulation(cutMesh.f, cutMesh.v), 'FaceColor', tmp)
% axis equal
% 

%% Histogram |mu| for each case
outfn = fullfile(imDir, sprintf('%06d_BeltramiCoefficients.png', tp)) ;
if ~exist(outfn, 'file') && save_ims
    clf
    maxx = max([max(abs(mu_annulus)), max(abs(mu_rectangle))]) ;
    % subplot(3, 1, 1)
    % histogram(abs(mu_uv))
    % xlim([0, maxx])
    % xlabel('$|\mu|$ for Dirichlet minimization', 'interpreter', 'latex') 
    % ylabel('counts', 'interpreter', 'latex')
    subplot(2, 1, 1)
    histogram(abs(mu_annulus))
    xlim([0, maxx])
    xlabel('$|\mu|$ for Ricci flow to annulus', 'interpreter', 'latex')
    ylabel('counts', 'interpreter', 'latex')
    subplot(2, 1, 2)
    histogram(abs(mu_rectangle))
    xlim([0, maxx])
    xlabel('$|\mu|$ for rectilinear domain from Ricci flow', 'interpreter', 'latex')
    ylabel('counts', 'interpreter', 'latex')
    sgtitle('Conformality test', 'interpreter', 'latex')
    saveas(gcf, outfn)
end


%% Save mu for this #iterations
mufn = sprintf(QS.fullFileBase.ricciMu, maxIter, tp) ;
save(mufn, 'mu_annulus', 'mu_rectangle')

if nargout > 1
    ricciMu = struct();
    ricciMu.mu_annulus = mu_annulus ;
    ricciMu.mu_rectangle = mu_rectangle ;
end
    
