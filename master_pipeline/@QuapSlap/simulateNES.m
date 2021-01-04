function simulateNES(QS, options)

%% Default options
strainstyle = 'hoop' ;
fixVolume = true ;
fixBoundary = false ;
fixCap = false ;
poisson_ratio = 0.5 ;
thickness = 0.1 ;
Rad = 50 ;       % radius of the tube
Len = 200 ;      % length of the tube
sigma = 5 ;     % width of the region to shrink
nU = QS.nU ;
nV = QS.nV ;
t0Pathlines = QS.t0set() ;

% Figure Options
figWidth = 16 ; 
figHeight = 10 ;
cmin = 0.95 ;
cmax = 1.05 ;

%% Add path (todo: make this automatic by putting NES code in gut_matlabl)
NESpath = '/mnt/data/code/NonEuclideanShells/NES/' ;
addpath_recurse(NESpath)

%% Path options
outRoot = fullfile(sprintf(QS.dir.pathlines.data, t0Pathlines), 'simulation') ;
exten = sprintf('_L%02dR%02d_sigma%0.3f_nu%0.2f_t%0.2f_strain%0.3f_%03dx%03d', ...
    Len, Rad, sigma, poisson_ratio, thickness, strain, nU, nV) ;
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

% Directed edge vectors in 2d
eij = V2d(eIDx(:,2), :) - V2d(eIDx(:,1), :);

% Target edge lengths
eL = sqrt( sum( eij.^2, 2 ) );
eL0 = eL ; 

% % Compute the bond orientation angles
dx = V2d(eIDx(:, 2)) - V2d(eIDx(:, 1)) ;
beta = acos(abs(eij(:, 1)) ./ eL0) ;

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

%% Make figure of angle wrt axial direction beta as a check
% outfn = fullfile(outdir, 'wire_beta_definition.png') ;
% disp(['Saving ' outfn])
% aux_plot_beta_pattern(betaCut, Vc, eIDxCut, outfn) ;
outfn = fullfile(outdir, 'wire_beta_definition_glued.png') ;
if ~exist(outfn, 'file')
    close all
    cmin_beta = 0; cmax_beta = 1 ;
    cmap = cubehelix(128,1.23,2.98,1.35,1.77,[0.17,0.98],[0.96,0.51]) ; 
    cID = max(1, sum(abs(sin(beta)) > linspace(cmin_beta, cmax_beta, length(cmap)), 2)) ;
    ecolors = cmap(cID, :) ;
    figure('visible', 'off')
    plotColoredLinesegs([VV(eIDx(:,1), :), VV(eIDx(:, 2), :)], ecolors, ...
        'linewidth', 10^3 / length(V)) ;
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
    % save figure
    disp(['Saving ' outfn])
    export_fig(outfn, '-r300', '-nocrop') ;
end

%% Initial face areas
a0 = 0.5 * doublearea(VV, FF) ;
V0 = VV ;
    
% Plotting options
cmap = bwr ;

% plot limits
[~,~,~,xyzlim] = QS.getXYZLims() ;

%% save simulation parameters
targetTheta = 0 ; % either 0 or 'quasistatic'
save(fullfile(outdir, 'simulation_parameters.mat'), ...
    'targetTheta', 'strainstyle', ...
    'fixBoundary', 'fixVolume', 'fixCap', 'V0', ...
    'thickness', 'poisson_ratio', 'capID', 'xyzlim', 'cutM')

%% Assign strain magnitudes based on 2d barycenters
bc2d = barycenter(VV, FF) ;

%% increment strain at rate given by 'strain'
for ii = 1:uint8(1/strain)  
    % Calculate Target Edge Lengths -------------------------------------------
    % Determine how bonds are strained
    switch distribution
        case 'experiment'
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
            titlestr = [titlestr ' $\sigma=$' num2str(sigma)] ;
                    
    end
    
    switch strainstyle
        case 'all'
            scalefactor = 1 + mag ;  
        case 'hoop'
            % Determine hoop strain from experimental dx 
            dxmag = dxinterp()
            scalefactor = 1 + dxmag .* abs(sin(beta)) ;
        case 'axial'
            scalefactor = 1 + mag .* abs(cos(beta)) ;
        case 'ring'
            scalefactor = 1 + mag .* abs(sin(beta)) ;
            scalefactor(abs(sin(beta))<0.99) = 1 ;
        case 'line'
            scalefactor = 1 + mag .* abs(cos(beta)) ;
            scalefactor(abs(cos(beta))<0.99) = 1 ;
    end

    eL = eL .* scalefactor ;

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
    
    %% Run Elastic Relaxation =================================================
    tic
    V0 = VV ;
    if fixBoundary && fixVolume && fixCap
        VV = minimizeElasticEnergy( FF, V0, eL, ...
            'TargetAngles', tarTheta, ...
            'Thickness', thickness, ...
            'Poisson', poisson_ratio, ...
            'MaxIterations', 1000, ...
            'iterDisplay', 1, ...
            'Alpha', 1, ...
            'Beta', 1, ...
            'FixBoundary', ...
            'targetVertices', capID(:), ...
            'targetLocations', V0(capID(:), :), ...
            'FixVolume', 'TargetVolume', targetVolume); 
    elseif fixBoundary && fixCap
        VV = minimizeElasticEnergy( FF, V0, eL, ...
            'TargetAngles', tarTheta, ...
            'Thickness', thickness, ...
            'Poisson', poisson_ratio, ...
            'MaxIterations', 1000, ...
            'iterDisplay', 1, ...
            'Alpha', 1, ...
            'Beta', 1, ...
            'targetVertices', capID(:), ...
            'targetLocations', V0(capID(:), :), ...
            'FixBoundary');    
    elseif fixCap && fixVolume
        VV = minimizeElasticEnergy( FF, V0, eL, ...
            'TargetAngles', tarTheta, ...
            'Thickness', thickness, ...
            'Poisson', poisson_ratio, ...
            'MaxIterations', 1000, ...
            'iterDisplay', 1, ...
            'Alpha', 1, ...
            'Beta', 1, ...
            'targetVertices', capID(:), ...
            'targetLocations', V0(capID(:), :), ...
            'FixVolume', 'TargetVolume', targetVolume);    
    elseif fixBoundary 
        VV = minimizeElasticEnergy( FF, V0, eL, ...
            'TargetAngles', tarTheta, ...
            'Thickness', thickness, ...
            'Poisson', poisson_ratio, ...
            'MaxIterations', 1000, ...
            'iterDisplay', 1, ...
            'Alpha', 1, ...
            'Beta', 1, ...
            'FixBoundary');    
    elseif fixVolume
        VV = minimizeElasticEnergy( FF, V0, eL, ...
            'TargetAngles', tarTheta, ...
            'Thickness', thickness, ...
            'Poisson', poisson_ratio, ...
            'MaxIterations', 1000, ...
            'iterDisplay', 1, ...
            'Alpha', 1, ...
            'Beta', 1, ...
            'FixVolume', ...
            'targetVertices', capLID(1), ...
            'targetLocations', V0(capLID(1), :)) ;
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

    %% View Results ===========================================================
    close all
    figure('visible', 'off')
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 figWidth figHeight]);  
        
    % Area ratio for faces
    areas = 0.5 * doublearea(VV, FF) ;
    astrain = (areas - a0) ./ a0 ;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot faces
    h = trisurf(triangulation(FF, VV), astrain, 'edgecolor', 'none');
    caxis([-1, 1])
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
    view(2);
    axis off
    
    % save figure
    export_fig(fullfile(outdir, sprintf('faces%03d', ii)), '-nocrop') ;
    close all
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot wire frame
    % Get indices of the colors to plot edges as
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
    view(2);
    axis off
    % save figure
    disp('Saving fig')
    export_fig(fullfile(outdir, sprintf('wire%03d', ii)), '-nocrop') ;
     
    disp('done')
    % set(gcf, 'visible', 'on')
    % waitfor(gcf)    
    
    % Save vertices
    save(fullfile(outdir, sprintf('vertices_%03d', ii)), 'VV', 'FF')
    
    %% Plot flow invariants
    close all      
    % Extract face barycenters for this timepoint in 3d
    bc0 = barycenter(V0, FF) ;
    % Extract for prev timpepoint in 3d
    bc = barycenter(VV, FF) ;
    v0 = bc - bc0 ;
    v0_vtx = VV - V0 ;
    % Now resolve the vector field for decomposition
    [v0n, v0t ] = ...
        resolveTangentNormalVelocities(FF, V0, v0, 1:length(FF)) ;

    % Construct DEC class instance
    DEC = DiscreteExteriorCalculus( FF, V0 ) ;
    % Compute divergence and rotational flow
    divv = DEC.divergence(v0t) ;
    rotv = DEC.curl(v0t) ;
    rotv(isnan(rotv)) = 0 ;
    
    % Plot div
    close all
    opts = struct() ;
    opts.clims = {10*strain, 10*strain} ;
    opts.alpha = 1.0 ;
    opts.labels = {'divergence', 'curl'} ;
    opts.cbarlabels = {'$\star$d$\star\left(v_t^\flat\right)$', ...
        '$\star$d$v_t^\flat$'} ;
    opts.xzylims = xyzlim ;
    opts.view = [0,90] ;
    opts.axisOff = true ;
    opts.cmap = bwr ;
    nFieldsOnSurface({FF, 0.5*(V0+VV)}, {divv, rotv}, opts) ;
    saveas(gcf, fullfile(outdir, sprintf('divcurl_%03d.png', ii)))
    close all    
    
    %% Plot Beltrami
    [ff, vv] = remove_vertex_from_mesh(FF, VV, capID) ;
    mesh_ii = struct('f', ff, 'v', vv, 'nU', nU) ;
    cutOpts.ignoreRectilinearConstraint = true ;
    mesh_ii = cutRectilinearCylMesh(mesh_ii, cutOpts) ;
    mu = bc_metric(mesh_ii.f, cutM.u, mesh_ii.v, 3) ;
    
    close all
    opts = struct() ;
    opts.clims = {1, 1} ;
    opts.alpha = 1.0 ;
    opts.labels = {'$\Re\mu$', '$\Im\mu$'} ;
    opts.cbarlabels = {'$\Re\mu$', '$\Im\mu$'} ;
    opts.xzylims = xyzlim ;
    opts.view = [0,90] ;
    opts.axisOff = true ;
    opts.cmap = bwr ;
    nFieldsOnSurface({ff, vv}, {real(mu), imag(mu)}, opts) ;
    saveas(gcf, fullfile(outdir, sprintf('beltrami_%03d.png', ii)))
    close all    
    
    % Prepare for next timepoint
    V0 = VV ;
end

