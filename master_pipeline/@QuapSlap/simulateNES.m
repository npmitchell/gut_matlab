function simulateNES(QS, options)

%% Default options
strainstyle = 'hoop' ;
fixVolume = true ;
fixBoundary = false ;
fixCap = false ;
poisson_ratio = 0.5 ;
thickness = 2 ;
nU = QS.nU ;
nV = QS.nV ;
t0Pathlines = QS.t0set() ;
preview = false ;
Dt = 30 ;
Ntotal = 5 ;
plot_faces = true ;
plot_dfaces = true ;
plot_wire = false ;
plot_divcurl = true ;

% Figure Options
figWidth = 16 ; 
figHeight = 10 ;
cmin = 0.95 ;
cmax = 1.05 ;
climit_div = 0.1 ;

%% Unpack options
if isfield(options, 'strainstyle')
    strainstyle = options.strainstyle ;
end
if isfield(options, 'fixVolume')
    fixVolume = options.fixVolume ;
end
if isfield(options, 'fixBoundary')
    fixBoundary = options.fixBoundary ;
end
if isfield(options, 'fixCap')
    fixCap = options.fixCap ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end

%% Add path (todo: make this automatic by putting NES code in gut_matlabl)
NESpath = '/mnt/data/code/NonEuclideanShells/NES/' ;
addpath_recurse(NESpath)

%% Path options
outRoot = fullfile(sprintf(QS.dir.pathlines.data, t0Pathlines), 'simulation') ;
exten = sprintf('_Dt%03d_nu%0.2f_t%0.2f_%03dx%03d', ...
    Dt, poisson_ratio, thickness, nU, nV) ;
exten = strrep(exten, '.', 'p') ;
if fixBoundary
    exten = ['_fixB' exten ] ;
end
if fixVolume
    exten = ['_fixV' exten ] ;
end
if fixCap
    exten = ['_fixC' exten ] ;
end
dirname = [strainstyle exten '_test'] ;
outdir = fullfile(outRoot, dirname) ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

%% 
QS.setTime(t0Pathlines)
QS.getCurrentSPCutMeshSmRS()
cutM = QS.currentMesh.spcutMeshSmRS ;
QS.getCurrentSPCutMeshSmRSC()
mesh = QS.currentMesh.spcutMeshSmRSC ;

% Close the endcaps with a single vertex (may be a problem?)
vtx = reshape(mesh.v, [nU, nV-1, 3]) ;
nvtx = size(mesh.v, 1) ;
endpt1 = mean(squeeze(vtx(1, :, :)), 1) ;
endpt2 = mean(squeeze(vtx(end, :, :)), 1) ;
end1 = 1:nU:nU*(nV-1) ;
endf1 = [];
for qq = 1:length(end1)
    if qq + 1 <= numel(end1) 
        endf1 = [endf1; [end1(qq), end1(qq+1), nvtx + 1]] ;
    else
        endf1 = [endf1; [end1(qq), end1(mod(qq+1, numel(end1))), nvtx + 1]] ;
    end
end
end2 = nU:nU:nU*(nV-1) ;
endf2 = [];
for qq = 1:length(end2)
    if qq + 1 <= numel(end2) 
        endf2 = [endf2; [end2(qq), nvtx+2, end2(qq+1)]] ;
    else
        endf2 = [endf2; [end2(qq), nvtx+2, end2(mod(qq+1, numel(end2)))]] ;
    end
end

% add endcap points to vtx and triangulation
VV = [mesh.v; endpt1; endpt2] ;
FF = [mesh.f; endf1; endf2] ;
tri = triangulation(FF, VV) ;

% Construct Topolgical Structure Tools ===================================
[eIDx, feIDx, bulkEdgeIDx] = topologicalStructureTools(tri) ;

% Check the triangulation
fnormals = faceNormal(tri) ;
trisurf(tri, 'faceVertexCData', fnormals(:, 1), 'edgecolor', 'none')
% Check that closed (chi for sphere = 2) 
chi = size(VV, 1) - size(eIDx, 1) + size(FF, 1) ;
assert(chi == 2)

% Save initial closed mesh
FF = [FF(:, 2), FF(:, 1), FF(:, 3)] ;
save(fullfile(outdir, 'initial_mesh.mat'), 'VV', 'FF', 'eIDx', 'feIDx', 'bulkEdgeIDx')


%% Construct Physical/Target Configurations and Geometries ================

% Construct Initial Configuration -----------------------------------------
% The initial configuration is a weakly buckled spherical cap.  This
% configuration is chosen to break the symmetry of the flat disk

% Directed edge vectors in 3d
eij = VV(eIDx(:,2), :) - VV(eIDx(:,1), :);

% Target edge lengths
eL = sqrt( sum( eij.^2, 2 ) );
eL0 = eL ; 

% Get bond centers in 3d
midx = 0.5 * (VV(eIDx(:, 1), 1) + VV(eIDx(:, 2), 1)) ;
midy = 0.5 * (VV(eIDx(:, 1), 2) + VV(eIDx(:, 2), 2)) ;
midz = 0.5 * (VV(eIDx(:, 1), 3) + VV(eIDx(:, 2), 3)) ;
bc3d = [midx, midy, midz] ;

%% Compute the bond orientation angles in 2d
refMesh = load(sprintf(QS.fileName.pathlines.refMesh, t0Pathlines)) ;
refMesh = refMesh.refMesh ;
V2d = refMesh.u ;
V2d(:, 1) = V2d(:, 1) / max(V2d(:, 1)) ;
% Construct Topolgical Structure Tools ===================================
tri2d = triangulation(refMesh.f, [V2d, V2d(:, 1)*0]) ;
[eIDx2d, feIDx2d, bulkEdgeIDx2d] = topologicalStructureTools(tri2d) ;
% glue it
e2dg = eIDx2d ;
% remove seam bonds
ghostBonds2d = find(all(eIDx2d > nU*(nV-1), 2)) ;
e2dg(ghostBonds2d, :) = [] ;
for qq = 1:nU
    e2dg(e2dg == nU*(nV-1) + qq) = qq ; 
end
% Note that indices must increase from left to right in each row
assert(all(eIDx(:, 2) - eIDx(:, 1) > 0))
toswap = find(e2dg(:, 2) - e2dg(:, 1) < 0) ;
for qq = toswap
    tmp = e2dg(qq, 1) ;
    e2dg(qq, 1) = e2dg(qq, 2) ;
    e2dg(qq, 2) = tmp ;
end
assert(all(e2dg(:, 2) - e2dg(:, 1) > 0))

eij2d = V2d(eIDx2d(:, 2), :) - V2d(eIDx2d(:, 1), :) ;
eL2d = vecnorm(eij2d, 2, 2) ;
dx = abs(eij2d(:, 1)) ;
beta2d = acos(dx ./ eL2d) ;
% make glued version which is missing longitudinal seam bonds
beta2dg = beta2d ;
beta2dg(ghostBonds2d) = [] ;
hoopID = find(abs(beta2d) < 0.01) ;
% Associate the hoop bonds with the 3D triangulation
[intx, i2d, i3d] = intersect(e2dg, eIDx, 'rows', 'stable');
assert(all(i2d' == 1:length(intx)))
% Note that size(i3d, 1) ~= size(beta2d, 1), since we removed longitudinal 
% bonds along seam
bbeta = zeros(size(eIDx, 1), 1) ;
bbeta(i3d) = beta2dg ;

% % Compute the bond orientation angles
% cutTri = triangulation(cutM.f, cutM.u) ;
% [eIDxCut, feIDxCut, bulkEdgeIDxCut] = topologicalStructureTools(cutTri) ;
% drp = rphi(eIDxCut(:, 2)) - rphi(eIDxCut(:, 1)) ;
% dx = ss(eIDxCut(:, 2)) - ss(eIDxCut(:, 1)) ;
% betaCut = atan2(drp, dx) ;
% beta = betaCut(1:(end-nU+1)) ;
% Vc = cutM.u ;
% 
% %%%%%%%

capID = [nU*(nV-1)+1, nU*(nV-1)+2]' ;

%% Make figure of angle wrt axial direction beta as a check
% outfn = fullfile(outdir, 'wire_beta_definition.png') ;
% disp(['Saving ' outfn])
% aux_plot_beta_pattern(betaCut, Vc, eIDxCut, outfn) ;
outfn = fullfile(outdir, 'wire_beta_definition_glued.png') ;
if ~exist(outfn, 'file')
    close all
    cmin_beta = 0; cmax_beta = 1 ;
    cmap = cubehelix(128,1.23,2.98,1.35,1.77,[0.17,0.98],[0.96,0.51]) ; 

    % Draw all bonds colored by angle
    % bondsBelow = find(bc3d(:, 3) < 0) ;
    cID = min(max(1, ...
        sum(abs(sin(bbeta)) > ...
            linspace(cmin_beta, cmax_beta, length(cmap)+1), 2)), ...
        length(cmap)) ;
    ecolors = cmap(cID, :) ;
    figure('visible', 'off')
    plotColoredLinesegs([VV(eIDx(:,1), :), ...
        VV(eIDx(:, 2), :)], ecolors, ...
        'linewidth', 10^3 / length(VV)) ;
    c = colorbar ;
    colormap(cmap)
    caxis([cmin_beta, cmax_beta])
    c.Color = 'w' ;
    c.Label.Interpreter = 'latex' ;
    c.Label.String = '$|\sin(\beta)|$' ;
    % Figure properties
    set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w')
    set(gcf, 'color', 'k')
    title('axial angle $\beta$ definition', 'interpreter', 'latex', 'color', 'w'); 
    axis equal
    axis off
    drawnow
    view(0,0)
    % save figure
    disp(['Saving ' outfn])
    export_fig(outfn, '-r300', '-nocrop') ;
end

%% Initial face areas
a0 = 0.5 * doublearea(VV, FF) ;
V0 = VV ;  % reference vertices
V1 = VV ;  % previous timepoint vertices
    
% Plotting options
cmap = bwr ;

% plot limits
[~,~,~,xyzlim] = QS.getXYZLims() ;

%% save simulation parameters
targetTheta = 0 ; % either 0 or 'quasistatic'
eIDx2dGlued = e2dg ;
save(fullfile(outdir, 'simulation_parameters.mat'), ...
    'targetTheta', 'strainstyle', ...
    'fixBoundary', 'fixVolume', 'fixCap', 'V0', ...
    'thickness', 'poisson_ratio', 'capID', 'xyzlim', 'cutM', ...
    'eIDx2d', 'eIDx2dGlued', 'i3d', 'i2d', 'ghostBonds2d', 'Dt')

%% Assign strain magnitudes based on 2d bond centers
midx = 0.5 * (V2d(eIDx2d(:, 1), 1) + V2d(eIDx2d(:, 2), 1)) ;
midy = 0.5 * (V2d(eIDx2d(:, 1), 2) + V2d(eIDx2d(:, 2), 2)) ;
bxy2d = [midx(:), midy(:)] ;

% glued version of bondcenters
bxy2dg = bxy2d ;
bxy2dg(ghostBonds2d, :) = [] ;

% check it
if preview
    close all; set(gcf, 'visible', 'on')
    triplot(refMesh.f, V2d(:, 1), V2d(:, 2))
    hold on;
    scatter(bxy2d(:, 1), bxy2d(:, 2), 10, 'filled')
    % plot without ghost bonds
    scatter(bxy2dg(:, 1), bxy2dg(:, 2), 50)
    pause(10)
    close all
end

%% Prep directories for images and output
if ~exist(fullfile(outdir, 'vertices'), 'dir')
    mkdir(fullfile(outdir, 'vertices'))
end
if ~exist(fullfile(outdir, 'beltrami'), 'dir')
    mkdir(fullfile(outdir, 'beltrami'));
end
if ~exist(fullfile(outdir, 'divcurl'), 'dir')
    mkdir(fullfile(outdir, 'divcurl'))
end
if ~exist(fullfile(outdir, 'faces'), 'dir')
    mkdir(fullfile(outdir, 'faces'))
end
if ~exist(fullfile(outdir, 'beltrami_images'), 'dir')
    mkdir(fullfile(outdir, 'beltrami_images'));
end
if ~exist(fullfile(outdir, 'wire'), 'dir')
    mkdir(fullfile(outdir, 'wire'))
end
if ~exist(fullfile(outdir, 'dfaces'), 'dir')
    mkdir(fullfile(outdir, 'dfaces'))
end

%% Plot initial beltrami (should be zero, or very nearly so
bfn = fullfile(outdir, 'beltrami', sprintf('beltrami_%03d.mat', 0)) ;
bifn = fullfile(outdir, 'beltrami_images', sprintf('beltrami_%03d.png', 0)) ;
if ~exist(bfn, 'file') || ~exist(bifn, 'file') 
    aux_nes_beltrami(QS, FF, VV, refMesh, capID, nU, nV, 0, bfn, bifn)
end

%% Consider each timestep, which averages Dt timepoints of experiment
for ii = 1:Ntotal
    tp = QS.t0set() + (ii-1)*Dt ;
    % Calculate Target Edge Lengths -------------------------------------------
    % Determine how bonds are strained
    switch strainstyle
        case 'all'
            titlestr = ['simulation with measured $\epsilon$'];
        case 'hoop'
            titlestr = ['simulation with measured $\epsilon_{\phi\phi}$'];  
        case 'axial'
            titlestr = ['simulation with measured $\epsilon_{\zeta\zeta}$'];  
        case 'ring'
            titlestr = ['simulation with measured ring strain $\epsilon_{\phi\phi}$'];  
        case 'line'
            titlestr = ['simulation with measured axial line strain of $\epsilon_{\zeta\zeta}$'];  
    end            
    
    %% Determine hoop strain from experimental dx and/or dy
    for qq = 1:Dt
        disp(['averaging dxs with tp=', num2str(tp + qq -1)])
        tmp = load(sprintf(QS.fileBase.pathlines.dxdyFiltered, t0Pathlines, tp + qq - 1)) ;
        if qq == 1
            dxs = tmp.dxs ;
            dys = tmp.dys ;
        else
            dxs = dxs + tmp.dxs ;
            dys = dys + tmp.dys ;
        end
    end
    dxs = dxs ;
    dys = dys ;

    % hack for now --> go back and repeat with reflected boundaries
    dxs(:, nV) = dxs(:, 1) ;
    dys(:, nV) = dys(:, 1) ;


    % prep for interpolation
    rux = refMesh.u(:, 1) / max(refMesh.u(:, 1)) ;
    ruy = refMesh.u(:, 2) / max(refMesh.u(:, 2)) ;
    rux = reshape(rux, [nU, nV]) ;
    ruy = reshape(ruy, [nU, nV]) ;
    if preview
        close all
        subplot(1, 2, 1)
        imagesc(linspace(0,1,nU), linspace(0,1,nV), dxs')
        xlabel('$\zeta$','interpreter', 'latex')
        ylabel('$\phi$','interpreter', 'latex')
        title('$\delta \zeta$', 'interpreter', 'latex')
        caxis([-max(abs(dxs(:))), max(abs(dxs(:)))])
        colormap bwr
        subplot(1, 2, 2)
        imagesc(linspace(0,1,nU), linspace(0,1,nV), dys')
        xlabel('$\zeta$','interpreter', 'latex')
        ylabel('$\phi$','interpreter', 'latex')
        title('$\delta \phi$', 'interpreter', 'latex')
        caxis([-max(abs(dxs(:))), max(abs(dxs(:)))])
        colormap bwr
    end

    if any(strcmpi(strainstyle, {'all', 'axial', 'line'}))

        rux = refMesh.u(:, 1) / max(refMesh.u(:, 1)) ;
        ruy = refMesh.u(:, 2) / max(refMesh.u(:, 2)) ;
        rux = reshape(rux, [nU, nV]) ;
        ruy = reshape(ruy, [nU, nV]) ;
        dxinterp = griddedInterpolant(rux, ruy, dxs, 'linear') ;
        dxi = dxinterp(bxy2dg(:, 1), bxy2dg(:, 2)) ;
        dxmag = zeros(size(bbeta)) ;
        dxmag(i3d) = dxi ;

        %% Preview 3d assignment
        if preview
            clear all
            scatter3(bc3d(:,1), bc3d(:, 2),bc3d(:,3), 15, dxmag, 'filled')
            caxis([-max(abs(dxs(:))), max(abs(dxs(:)))])
            colormap bwr
            pause(3)
        end
    end
    if any(strcmpi(strainstyle, {'all', 'hoop', 'ring'}))
        dyinterp = griddedInterpolant(rux, ruy, dys, 'linear') ;
        dyi = dyinterp(bxy2dg(:, 1), bxy2dg(:, 2)) ;
        dymag = zeros(size(bbeta)) ;
        dymag(i3d) = dyi ;

        %% Preview 3d assignment
        if preview
            clear all
            midx = 0.5 * (VV(eIDx(:, 1), 1) + VV(eIDx(:, 2), 1)) ;
            midy = 0.5 * (VV(eIDx(:, 1), 2) + VV(eIDx(:, 2), 2)) ;
            midz = 0.5 * (VV(eIDx(:, 1), 3) + VV(eIDx(:, 2), 3)) ;
            scatter3(midx, midy, midz, 15, dymag, 'filled')
            caxis([-max(abs(dys(:))), max(abs(dys(:)))])
            colormap bwr
            pause(3)
        end
    end

    %% check seam
    % % Look at all glued bonds that touch seam
    % clf
    % [rs, cs] = find(e2dg < nU+1) ;
    % rs = unique(rs) ;
    % subplot(2, 1, 1)
    % hold on;
    % for tmpId = 1:length(rs)
    %     rr = rs(tmpId) ;
    %     ids = e2dg(rr, :) ;
    %     plot(V2d(ids(:),1), V2d(ids(:), 2), '.-', ...
    %         'color', [0, 0, 0, 0.25])
    % end
    % subplot(2, 1, 2)
    % hold on;
    % for tmpId = 1:length(rs)
    %     rr = rs(tmpId) ;
    %     ids = e2dg(rr, :) ;
    %     plot3(VV(ids,1), VV(ids, 2), VV(ids, 3), '.-', ...
    %         'color', [0, 0, 0, 0.25])
    % end
    % % Look at all cut bonds that touch seam
    % [rs, cs] = find(eIDx2d < nU+1 | eIDx2d > nU*(nV-1)) ;
    % rs = unique(rs) ;
    % subplot(2, 1, 2)
    % hold on;
    % for tmpId = 1:length(rs)
    %     rr = rs(tmpId) ;
    %     ids = eIDx2d(rr, :) ;
    %     plot(V2d(ids(:),1), V2d(ids(:), 2), 'k.-')
    % end
    % for tmpId = 1:length(rs)
    %     rr = rs(tmpId) ;
    %     ids = e2dg(rr, :) ;
    % scatter(bxy2dg(:, 1), bxy2dg(:, 2), 10, dxi, 'filled')

    % check domain
    % scatter(rux(:), ruy(:), 5, dxs(:))
    % hold on;
    % scatter(bxy2d(:, 1), bxy2d(:, 2), 5, 'filled')

    switch strainstyle
        case 'all'
            scalefactor = 1 + mag ;  
        case 'hoop'
            scalefactor = 1 + dymag .* abs(sin(bbeta)) ;
        case 'axial'
            scalefactor = 1 + mag .* abs(cos(bbeta)) ;
        case 'ring'
            scalefactor = 1 + mag .* abs(sin(bbeta)) ;
            scalefactor(abs(sin(bbeta))<0.99) = 1 ;
        case 'line'
            scalefactor = 1 + mag .* abs(cos(bbeta)) ;
            scalefactor(abs(cos(bbeta))<0.99) = 1 ;
    end

    eL = eL .* scalefactor ;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot initial strain assignment as wire frame
    outfn0 = fullfile(outdir, 'wire0.png') ;
    if ii == 1 && ~exist(outfn0, 'file') && plot_wire
        disp('Saving initial wire deformation fig (may be slow...)')
        % Get indices of the colors to plot edges as
        close all
        figure('visible', 'off')
        colormap bwr
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 figWidth figHeight]);  
        sratio = eL ./ eL0 ;
        cID = max(1, sum(sratio > linspace(cmin, cmax, length(cmap)), 2)) ;
        ecolors = cmap(cID, :) ;
        lsegs = [VV(eIDx(:,1), :), VV(eIDx(:, 2), :)] ;
        plotColoredLinesegs(lsegs, ecolors) ;
        c = colorbar ;
        caxis([cmin, cmax])
        c.Color = 'w' ;
        c.Label.Interpreter = 'latex' ;

        c.Label.String = '$\ell/\ell_0$' ;
        % Figure properties
        set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w')
        set(gcf, 'color', 'k')
        title(titlestr, 'interpreter', 'latex', 'color', 'w'); 
        axis equal
        view(0,0)
        axis off

        % save figure
        export_fig(outfn0, '-nocrop') ;
        disp('saved initial wire deformation')
    end

    % Calculate Target Bending Angles -----------------------------------------
    if strcmpi(targetTheta, 'quasistatic')
        % The unit normal vector for each face
        % see https://www.cs.utexas.edu/users/evouga/uploads/4/5/6/8/45689883/turning.pdf
        fN = faceNormal( triangulation(F, tarV) );

        N1 = fN(edgeFace(:,1), :);
        N2 = fN(edgeFace(:,2), :);

        crossN = cross(N2, N1, 2);
        dotN = dot(N1, N2, 2);

        tarTheta = 2 .* atan2( dot(crossN, eij, 2), 1 + dotN );
    else
        tarTheta = zeros(size(eL, 1), 1) ;
    end
    
    %% Compute initial volume
    % The centroids of each face
    COM = cat( 3, VV(FF(:,1), :), ...
        VV(FF(:,2), :), VV(FF(:,3), :) );
    COM = mean(COM, 3);

    % The area weighted face normal vectors
    ej = VV(FF(:,1), :) - VV(FF(:,3), :);
    ek = VV(FF(:,2), :) - VV(FF(:,1), :);
    ndirpts = cross(ej, ek, 2);

    targetVolume = abs(sum( dot(COM, ndirpts, 2) ) ./ 6 );
    
    %% Run Elastic Relaxation =============================================
    fn = fullfile(outdir, 'vertices', sprintf('vertices_%03d.mat', ii)) ;
    if exist(fn, 'file')
        disp('time step already computed, loading...')
        tmp = load(fn) ;
        VV = tmp.VV ;
        FF = tmp.FF ;
    else
        tic
        if fixBoundary && fixVolume && fixCap
            VV = minimizeElasticEnergy( FF, V1, eL, ...
                'TargetAngles', tarTheta, ...
                'Thickness', thickness, ...
                'Poisson', poisson_ratio, ...
                'MaxIterations', 1000, ...
                'iterDisplay', 1, ...
                'Alpha', 1, ...
                'Beta', 1, ...
                'FixBoundary', ...
                'targetVertices', capID(:), ...
                'targetLocations', V1(capID(:), :), ...
                'FixVolume', 'TargetVolume', targetVolume); 
        elseif fixBoundary && fixCap
            VV = minimizeElasticEnergy( FF, V1, eL, ...
                'TargetAngles', tarTheta, ...
                'Thickness', thickness, ...
                'Poisson', poisson_ratio, ...
                'MaxIterations', 1000, ...
                'iterDisplay', 1, ...
                'Alpha', 1, ...
                'Beta', 1, ...
                'targetVertices', capID(:), ...
                'targetLocations', V1(capID(:), :), ...
                'FixBoundary');    
        elseif fixCap && fixVolume
            VV = minimizeElasticEnergy( FF, V1, eL, ...
                'TargetAngles', tarTheta, ...
                'Thickness', thickness, ...
                'Poisson', poisson_ratio, ...
                'MaxIterations', 1000, ...
                'iterDisplay', 1, ...
                'Alpha', 1, ...
                'Beta', 1, ...
                'targetVertices', capID(:), ...
                'targetLocations', V1(capID(:), :), ...
                'FixVolume', 'TargetVolume', targetVolume);    
        elseif fixBoundary 
            VV = minimizeElasticEnergy( FF, V1, eL, ...
                'TargetAngles', tarTheta, ...
                'Thickness', thickness, ...
                'Poisson', poisson_ratio, ...
                'MaxIterations', 1000, ...
                'iterDisplay', 1, ...
                'Alpha', 1, ...
                'Beta', 1, ...
                'FixBoundary');    
        elseif fixVolume
            VV = minimizeElasticEnergy( FF, V1, eL, ...
                'TargetAngles', tarTheta, ...
                'Thickness', thickness, ...
                'Poisson', poisson_ratio, ...
                'MaxIterations', 1000, ...
                'iterDisplay', 1, ...
                'Alpha', 1, ...
                'Beta', 1, ...
                'FixVolume', ...
                'targetVertices', capID(end-1), ...
                'targetLocations', V1(capID(end-1), :)) ;
        else
            error('not fixing Volume or boundary?')
        end
        toc

        %% Check
        % % Directed edge vectors
        % eijp = VV(eIDx(:,2), :) - VV(eIDx(:,1), :);
        % 
        % % Target edge lengths
        % eLp = sqrt( sum( eijp.^2, 2 ) );

        %% Save vertices
        save(fn, 'VV', 'FF')
    end

    % Area ratio for faces
    areas = 0.5 * doublearea(VV, FF) ;
    a0strain = (areas - a0) ./ a0 ;
    a1 = 0.5 * doublearea(V1, FF) ;
    a1strain = (areas - a1) ./ a1 ;    

    %% View Results ===========================================================   
    % Plot faces
    facesfn = fullfile(outdir, 'faces', sprintf('faces%03d.png', ii)) ;
    if ~exist(facesfn, 'file') && plot_faces
        disp('Creating figure with colored faces by Delta(face area)')
        close all
        figure('visible', 'off')
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 figWidth figHeight]);  

        trisurf(triangulation(FF, VV), a0strain, 'edgecolor', 'none');
        caxis([-0.5, 0.5])
        colormap(bwr)
        c = colorbar ;
        c.Color = 'w' ;
        c.Label.Interpreter = 'latex' ;
        c.Label.String = '$\delta A / A_0$' ;

        % Figure properties
        set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w')
        set(gcf, 'color', 'k')

        title(titlestr, 'interpreter', 'latex', 'color', 'w'); 
        axis equal
        view(0,0)
        axis off

        % save figure
        export_fig(facesfn, '-nocrop') ;
        close all
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot dfaces
    dfacesfn = fullfile(outdir, 'dfaces', sprintf('dfaces%03d.png', ii)) ;
    if ~exist(dfacesfn, 'file') && plot_dfaces
        disp('Creating figure with colored faces by d(facearea)/step')
        close all
        trisurf(triangulation(FF, VV), a1strain, 'edgecolor', 'none');
        caxis([-0.05, 0.05])
        colormap(bwr)
        c = colorbar ;
        c.Color = 'w' ;
        c.Label.Interpreter = 'latex' ;
        c.Label.String = '$(A - A_{prev})/ A_{prev}$' ;

        % Figure properties
        set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w')
        set(gcf, 'color', 'k')

        title(titlestr, 'interpreter', 'latex', 'color', 'w'); 
        axis equal
        view(0,0)
        axis off

        % save figure
        export_fig(dfacesfn, '-nocrop') ;
        close all
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot wire frame
    % Get indices of the colors to plot edges as
    wirefn = fullfile(outdir, 'wire', sprintf('wire%03d.png', ii)) ;
    if ~exist(wirefn, 'file') && plot_wire 
        aux_nes_wire(QS, eIDx, VV, eL, eL0, clims, wirefn, ii*Dt)
    end
        
    %% Plot flow invariants
    divcurlfn = fullfile(outdir, 'divcurl', sprintf('divcurl_%03d.png', ii)) ;
    if ~exist(divcurlfn, 'file') && plot_divcurl
        aux_nes_divcurl(QS, FF, VV, V1, climit_div, divcurlfn, ii*Dt)
    end
    
    %% Compute Beltrami
    bfn = fullfile(outdir, 'beltrami', sprintf('beltrami_%03d.mat', ii)) ;
    bifn = fullfile(outdir, 'beltrami_images', ...
        sprintf('beltrami_%03d.png', ii)) ;
    if ~exist(bfn, 'file') || ~exist(bifn, 'file')
        aux_nes_beltrami(QS, FF, VV, refMesh, capID, nU, nV, ii*Dt, bfn, bifn)
    end
    % Prepare for next timepoint
    V1 = VV ;
end

