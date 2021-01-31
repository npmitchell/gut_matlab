%% Example for using Ricci flow to conformally map a tube 
% Generates utilde, vtilde pullbacks for a gut mesh
% Example usage demonstration NPMitchell 2021

% Add path if not already done
addpath('../../RicciFlow_MATLAB/')
addpath('../addpath_recurse/')
addpath_recurse('../master_pipeline/')
addpath_recurse('../data_handling/')
addpath_recurse('../mesh_handling/')
addpath_recurse('../geometry/')
addpath_recurse('../plotting/')
disp('done adding paths')

%% Default parameters
radiusTolerance = 0.1 ;
tidx = 151 ;  % which timepoint ID is this in the experiment?
tidx0 = 151 ; % which timepoint ID in the experiment do we register to? (here ignore this)

%% Unpack parameters
% cutMesh = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
cutMesh = load('spcutMeshSmRS_000160.mat') ;
cutMesh = cutMesh.spcutMeshSmRS ;
glueMesh = glueCylinderCutMeshSeam(cutMesh) ;
nU = cutMesh.nU ;
nV = cutMesh.nV ;


maxCircIters = [4,8,10, 16, 25, 50, 100, 200, 400] ;
%%
for ii = 1:length(maxCircIters)
    MaxCircIter = maxCircIters(ii) ;
    outdir = ['example_output_iter' num2str(MaxCircIter)] ;
    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end

    % Check it
    trisurf(triangulation(glueMesh.f, glueMesh.v))
    axis equal

    %% Generate conformal parameterization in the unit disk
    % Load or compute ricci flow solution
    ricciFn = ['ricci_flow_solution_' num2str(MaxCircIter) 'iter.mat'] ;
    try 
        load(ricciFn, 'U')
    catch
        [~, U, ~] = DiscreteRicciFlow.EuclideanRicciFlow(glueMesh.f, glueMesh.v, ...
            'BoundaryType', 'Fixed', 'BoundaryShape', 'Circles', ...
            'MaxCircIter', MaxCircIter);
        % [labels, dbonds, topStructTools] = labelRectilinearMeshBonds(cutMesh) ;
        % Let MaxCircIter >= 50
        save(ricciFn, 'U') ;
    end

    % Plot Ricci result in annulus
    outfn = fullfile(outdir, 'example_00_RicciSolution.png') ;
    if ~exist(outfn, 'file')
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
    mu = bc_metric(glueMesh.f, U, glueMesh.v, 3) ;

    %% plot beltrami for Ricci flow
    mesh2d = triangulation(glueMesh.f, [U(:, 1), U(:, 2)]) ; %  ; , 0*U(:, 1)]) ;
    options.view = [0, 90] ;
    options.axisOff = true ;
    options.visible = true ;
    options.labels = {'$\Re\mu$', '$\Im\mu$', '$|\mu|$', ...
        '$\Re\mu$', '$\Im\mu$', '$|\mu|$'} ;
    outfn = fullfile(outdir, 'example_01_RicciFlowBeltrami_00.png') ;
    if ~exist(outfn, 'file')
        clf
        nFieldsOnSurface({mesh2d, mesh2d, mesh2d, ...
            glueMesh, glueMesh, glueMesh}, ...
            {real(mu), imag(mu), abs(mu), real(mu), imag(mu), abs(mu)}, ...
            options)
        sgtitle('Beltrami coefficients from Ricci flow', ...
            'interpreter', 'latex')
        saveas(gcf, outfn)
    end

    %% Compare to beltrami from UVprime cutmesh
    uvMesh = load('./uvpcutMesh_000160.mat') ;
    uvMesh = uvMesh.uvpcutMesh.resampled ;
    m2duv = struct() ;
    m2duv.f = uvMesh.f ;
    m2duv.v = uvMesh.u ;
    mu_uv = bc_metric(uvMesh.f, uvMesh.u, uvMesh.v, 3) ;

    % plot beltrami
    options.view = [0, 90] ;
    options.axisOff = true ;
    options.visible = true ;
    options.labels = {'$\Re\mu$', '$\Im\mu$', '$|\mu|$', ...
        '$\Re\mu$', '$\Im\mu$', '$|\mu|$'} ;
    outfn = fullfile(outdir, 'example_02_DirichletBeltrami.png') ;
    if ~exist(outfn, 'file')
        clf
        nFieldsOnSurface({m2duv, m2duv, m2duv, ...
            glueMesh, glueMesh, glueMesh}, ...
            {real(mu_uv), imag(mu_uv), abs(mu_uv), ...
            real(mu_uv), imag(mu_uv), abs(mu_uv)}, ...
            options)
        sgtitle('Beltrami coefficients from Dirichlet energy minimization', ...
            'interpreter', 'latex')
        saveas(gcf, outfn)
    end
    
    %% Plot just the bare triangulation
    outfn = fullfile(outdir, 'example_01b_RicciFlowSolution.png') ;
    if ~exist(outfn, 'file')
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
    outfn1 = fullfile(outdir, 'example_04_ricci_InnerCorrection.png') ;
    outfn2 = fullfile(outdir, 'example_05_ricci_OuterCorrection.png') ;
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

    %% Branch cut
    % THIS DOES NOT WORK
    % branch = [0,0, 2, 0] ;
    % E = triangulation(glueMesh.f, UU).edges ; 
    % xy2 = [UU(E(:,1), :), UU(E(:,2), :)] ;
    % intersections = lineSegmentIntersect(branch, xy2) ;
    % % which lines intersect
    % linesThatIntersect = find(intersections.intAdjacencyMatrix) ;
    % % order them by increasing x coordinate along branch 
    % xs = intersections.intMatrixX(linesThatIntersect) ;
    % [xsort, sortIDs] = sort(xs) ;
    % cutPath = E(linesThatIntersect(sortIDs), 1) ;
    % 
    % clf
    % triplot(triangulation(glueMesh.f, UU), 'Color', [0.8,0.8,0.8]) ;
    % hold on;
    % plot(UU(cutPath, 1), UU(cutPath, 2), 'ro-') ;
    % plot([0, 1], [0, 0], 'k--')

    % USE MATLAB graph 

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
    outfn = fullfile(outdir, 'example_06_DrhoDphi_rectified.png') ;
    if ~exist(outfn, 'file')
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
    
    outfn = fullfile(outdir, 'example_07_DrhoDphi_flipped.png') ;
    if ~exist(outfn, 'file')
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
    outfn = fullfile(outdir, 'example_08_phiOrderInitial.png') ;
    if ~exist(outfn, 'file')
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
    outfn = fullfile(outdir, 'example_09_phiOrderFinal.png') ;
    if ~exist(outfn, 'file')
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
    uvtMesh = struct() ;
    uvtMesh.annulus = struct('f', glueMesh.f, 'u', UU, 'v', glueMesh.v, ...
        'nU', nU, 'nV', nV) ;
    uvtcutMesh = struct() ;
    uvtcutMesh = uvtMesh.annulus ;
    uvtcutMesh.u = [rho, phi_recut] ;
    opts = struct('vmax', 2 * pi, 'ignoreRectangularConstraint', true) ;
    uvtcutMesh = cutRectilinearCylMesh(uvtcutMesh, opts) ;
    uvtMesh.rectangle = struct('f', uvtcutMesh.f, 'u', ...
        uvtcutMesh.u, 'v', uvtcutMesh.v, ...
        'nU', nU, 'nV', nV) ;
    
    outfn = fullfile(outdir, 'example_09b_phiOrderFinal.png') ;
    if ~exist(outfn, 'file')
        facecolors = (1:length(uvtMesh.rectangle.f)) ;
        facecolors = mod(facecolors, 2*nV-2) ;
        facecolors(facecolors==0) = 2*nV - 2 ;
        trisurf2d(uvtMesh.rectangle.f, uvtMesh.rectangle.u, 'FaceVertexCData', facecolors(:), 'Facealpha',0.5)

        saveas(gcf, outfn)
    end


    %% Recompute mu in rectangular space 
    % It should be conformal since this is a conformal map of a conformal map
    % into the disk
    muRect = bc_metric(uvtcutMesh.f, uvtcutMesh.u, uvtcutMesh.v, 3) ;

    %% plot beltrami for Ricci flow
    mesh2d = triangulation(uvtcutMesh.f, uvtcutMesh.u) ; %  ; , 0*U(:, 1)]) ;
    options.view = [0, 90] ;
    options.axisOff = true ;
    options.visible = true ;
    options.labels = {'$\Re\mu$', '$\Im\mu$', '$|\mu|$', ...
        '$\Re\mu$', '$\Im\mu$', '$|\mu|$'} ;
    
    outfn = fullfile(outdir, 'example_10_RicciFlowBeltrami_rectangle.png') ;
    if ~exist(outfn, 'file')
        clf
        nFieldsOnSurface({mesh2d, mesh2d, mesh2d, ...
            glueMesh, glueMesh, glueMesh}, ...
            {real(muRect), imag(muRect), abs(muRect), ...
            real(muRect), imag(muRect), abs(muRect)}, ...
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
    outfn = fullfile(outdir, 'example_11_BeltramiCoefficients.png') ;
    if ~exist(outfn, 'file')
        clf
        maxx = max([max(abs(mu_uv)), max(abs(mu)), max(abs(muRect))]) ;
        subplot(3, 1, 1)
        histogram(abs(mu_uv))
        xlim([0, maxx])
        xlabel('$|\mu|$ for Dirichlet minimization', 'interpreter', 'latex') 
        ylabel('counts', 'interpreter', 'latex')
        subplot(3, 1, 2)
        histogram(abs(mu))
        xlim([0, maxx])
        xlabel('$|\mu|$ for Ricci flow to annulus', 'interpreter', 'latex')
        ylabel('counts', 'interpreter', 'latex')
        subplot(3, 1, 3)
        histogram(abs(muRect))
        xlim([0, maxx])
        xlabel('$|\mu|$ for rectilinear domain from Ricci flow', 'interpreter', 'latex')
        ylabel('counts', 'interpreter', 'latex')
        sgtitle('Conformality test', 'interpreter', 'latex')
        saveas(gcf, outfn)
    end
    
    
    %% Save mu for this #iterations
    save(fullfile(outdir, ['mus_iter' num2str(MaxCircIter) '.mat']), 'mu', 'muRect')
end


%% Compare different
for ii = 1:length(maxCircIters)
    MaxCircIter = maxCircIters(ii) ;
    outdir = ['example_output_iter' num2str(MaxCircIter)] ;
    % load and stuff into cells
    load(fullfile(outdir, ['mus_iter' num2str(MaxCircIter) '.mat']), 'mu', 'muRect')
    logmus{ii} = log10(abs(mu)) ;
    logmuRs{ii} = log10(abs(muRect)) ;
    mus{ii} = abs(mu) ;
    muRs{ii} = abs(muRect) ;
    means(ii) = mean(abs(mu)) ;
    meanRs(ii) = mean(abs(muRect)) ;
    medians(ii) = median(abs(mu)) ;
    medianRs(ii) = median(abs(muRect)) ;
    stds(ii) = std(abs(mu)) ;
    stdRs(ii) = std(abs(muRect)) ;
end
colors = define_colors;
color1 = colors(1, :) ;
color2 = colors(2, :) ;
mu2plot = {mus, muRs} ;
logmu2plot = {logmus, logmuRs} ;
mean2plot = {means, meanRs} ;
median2plot = {medians, medianRs} ;
std2plot = {stds, stdRs} ;
titles = {'Conformality convergence for Ricci flow to annular domain', ...
    'Conformality convergence for Ricci flow to rectangle'} ;
fns = {'conformality_scaling_annulus', ...
    'conformality_scaling_rectangle'} ;

for qq = 1:2
    mus = mu2plot{qq} ;
    logmus = logmu2plot{qq} ;
    means = mean2plot{qq} ;
    medians = median2plot{qq} ;
    stds = std2plot{qq} ;
    titleq = titles{qq} ;

    % Plot mus for annuli
    close all
    subplot(2,2, 1)
    violin(mus,'x',log10(maxCircIters), ...
        'facecolor',colors(1:length(maxCircIters), :),'edgecolor','none',...
        'bw',1e-3, 'mc',[],'medc', [])
    hold on;
    plot(log10(maxCircIters), means, '.-', 'color', color1)
    plot(log10(maxCircIters), medians, 'o-', 'color', color2)
    xlabel('$\log_{10}($iterations$)$', 'interpreter', 'latex')
    ylabel('$\mu$', 'interpreter', 'latex')
    ylim([0, 0.1])
    subplot(2,2, 2)
    violin(logmus,'x',log10(maxCircIters), ...
        'facecolor',colors(1:length(maxCircIters), :),'edgecolor','none',...
        'bw',1e-3, 'mc',[],'medc',[])
    hold on;
    plot(log10(maxCircIters), log10(means), '.-', 'color', color1)
    plot(log10(maxCircIters), log10(medians), 'o-', 'color', color2)
    ylim([-3, 0])
    xlabel('$\log_{10}($iterations$)$', 'interpreter', 'latex')
    ylabel('$\log_{10}\mu$', 'interpreter', 'latex')
    ax = subplot(2, 2, 3) ;
    plot(log10(maxCircIters), means, '.-', 'color', color1)
    hold on;
    plot(log10(maxCircIters), medians, 'o-', 'color', color2)
    xlabel('$\log_{10}($iterations$)$', 'interpreter', 'latex')
    ylabel('$\langle\mu\rangle$, median$(\mu)$', 'interpreter', 'latex')
    legend({'mean', 'median'}, 'location', 'best')
    set(ax, 'fontsize', 12)
    ax = subplot(2, 2, 4) ;
    plot(log10(maxCircIters), stds, '.-', 'color', 'k')
    xlabel('$\log_{10}($iterations$)$', 'interpreter', 'latex')
    ylabel('$\sigma_{\mu}$', 'interpreter', 'latex')
    set(ax, 'fontsize', 12) 
    sgtitle(titleq, 'interpreter', 'latex')
    saveas(gcf, [fns{qq} '.pdf'])
    saveas(gcf, [fns{qq} '.png'])
end
