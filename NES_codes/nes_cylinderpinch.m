%% NonEuclidean_Shell for pinching cylinder ===============================
% This is a script to exhibit the basic functionality of the
% (N)on-(E)uclidean (S)hell (S)imulator (NESS) using a cylindrical geometry
% which is pinched along its midline.

clear; close all; clc;
distribution = 'gaussian' ;
strainstyle = 'ringstrain' ; % 'ringstrain', 'isopulse' 'hoopstrain' 'halfhoop' 'linearhoop'
fixVolume = true ;
fixBoundary = false ;
fixCap = false ;
poisson_ratio = 0.5 ;
thickness = 0.1 ;
Rad = 50 ;       % radius of the tube
Len = 200 ;      % length of the tube
sigma = 5 ;     % width of the region to shrink
strain = 0.005 ; % strain magnitude per increment in time
nU = 50 ;
nV = 50 ;
npts = round((nU / (2*pi))^2 * 4*pi) ;

% Figure Options
figWidth = 16 ; 
figHeight = 10 ;
cmin = 0.95 ;
cmax = 1.05 ;

%% Add paths
if isunix
    % VIP8
    gutDir = '/mnt/data/code/gut_matlab/' ;
    outRoot = '/mnt/data/simulations/NES_cylinder/' ;
    NESpath = '/mnt/data/code/NonEuclideanShells/NES/' ;
elseif ismac
    % KITP MAC OSX Laptop
    NESpath = fullfile(gutDir, '../NES') ; % '~/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/NES/'; 
    gutDir = '../' ;
    outRoot = '~/Desktop/gut/simulations/NES_cylinder/' ;
end

%%
addpath(fullfile(gutDir, 'addpath_recurse')) ;
addpath_recurse(fullfile(gutDir, 'gptoolbox')) ;
addpath_recurse(fullfile(gutDir, 'geometry')) ;
addpath_recurse(fullfile(gutDir, 'mesh_handling')) ;
addpath_recurse(fullfile(gutDir, 'plotting')) ;
addpath_recurse(NESpath)
addpath_recurse(fullfile(gutDir, 'NES_codes')) ;


%% Mex the NES path
thisDir = pwd ;
cd(fullfile(NESpath, 'MATLAB'))
if ~exist('./minimize_elastic_energy.cpp', 'file')
    error('making cpp')
    mex minimize_elastic_energy.cpp
end
cd(thisDir)

%% Path options
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
cd(outdir)

%% --------------------------------------------------------------------------
% Construct 2D mesh corresponding to the planar domain of
% parameterization
%--------------------------------------------------------------------------
% % Load from file
% path = fullfile(NESpath, 'NES_Examples') ;
% mesh = read_ply_mod(fullfile(path, 'tube_simple_h1p00_R0p00_w1p00.ply')) ;
% rmID = [length(mesh.v)-1, length(mesh.v)] ;
% [F, V] = remove_vertex_from_mesh(mesh.f, mesh.v, rmID) ;
% % Center the mesh around the x axis
% midx = 0.5 * (max(mesh.v(:, 1)) + min(mesh.v(:, 1))) ;
% V(:, 1) = V(:, 1) - midx ;

% Build from scratch
[phi, sg] = meshgrid(linspace(0,2*pi,nV), linspace(-Len*0.5, Len*0.5, nU)) ;
ss = sg(:) ;
sp = [ss, Rad * phi(:)] ;
vcut = [ss, Rad * cos(sp(:, 2)), Rad* sin(sp(:, 2))] ;
% Scale the azimuthal direction (phi) to rr
rphi = Rad * phi(:) ;
faces = defineFacesRectilinearGrid(sp, nU, nV) ;
% Now glue it together
cutM.nU = nU; cutM.nV = nV ;
cutM.u = sp ;
cutM.f = faces ;
cutM.v = vcut ;
mesh = glueCylinderCutMeshSeam(cutM) ;
V = mesh.v ;
F = mesh.f ;

% Construct Triangulation
tri = triangulation(F, V) ;

%% Determine boundary
bnd = tri.freeBoundary ;
bnd = bnd(:, 1) ;
% left boundary
left = bnd(find(V(bnd, 1) < 0)) ;
right = bnd(find(V(bnd, 1) > 0)) ;

% Construct Topolgical Structure Tools ===================================
[eIDx, feIDx, bulkEdgeIDx] = topologicalStructureTools(tri) ;

minBondL = min(vecnorm(V(eIDx(:, 1), :) - V(eIDx(:, 2), :), 2, 2)) ;

% Save initial mesh
save(fullfile(outdir, 'initial_mesh_freeB.mat'), 'V', 'F', 'eIDx', 'feIDx', 'bulkEdgeIDx')
% Plot the resulting mesh
lpt = mean(V(left, :)) ;
rpt = mean(V(right, :)) ;
cla
set(gcf, 'visible', 'off')
plot3(V(:, 1), V(:, 2), V(:, 3), 'k.')
hold on; 
plot3(V(bnd, 1), V(bnd, 2), V(bnd, 3), 'co')
plot3(V(left, 1), V(left, 2), V(left, 3), '.')
plot3(V(right, 1), V(right, 2), V(right, 3), '.')
scatter3(lpt(:, 1), lpt(:, 2), lpt(:, 3), 'filled')
scatter3(rpt(:, 1), rpt(:, 2), rpt(:, 3), 'filled')
xlabel('x')
ylabel('y')
zlabel('z')
drawnow
axis equal
saveas(gcf, fullfile(outdir, 'initial_mesh_tube.png'))

%% OPTION 1
% % Now add endcaps vertices (to be made virtual) -- spherical caps
% [ points, sphereTri ] = nSphereMesh(npts, 'ico') ;
% points = points * Rad ;
% [~, oldId, newId_in_old] = remove_duplicate_vertices(points, 0.1) ;
% id2remove = setdiff(1:length(points), oldId) ; 
% [ sphereTri, points] = ...
%     remove_vertex_from_mesh( sphereTri, points, id2remove) ;
% disp(['removed ' num2str(length(id2remove)) ' duplicate vertices']) 
% disp([num2str(length(points)) ' vertices remain'])
% % pts = vogelSphere(npts) ;
% 
% % Divide the points into two hemispheres
% % lcap = points(points(:, 1) < 0, :) - [Len*0.5,0,0] ;
% % rcap = points(points(:, 1) > 0, :) + [Len*0.5,0,0] ;
% 
% % Free boundary of lcap
% [lcapF, lcap] = remove_vertex_from_mesh(sphereTri, points, find(points(:, 1) >= -eps)) ;
% [rcapF, rcap] = remove_vertex_from_mesh(sphereTri, points, find(points(:, 1) <= eps)) ;
% triL = triangulation(lcapF, lcap) ;
% triR = triangulation(rcapF, rcap) ;
% bndL = unique(triL.freeBoundary()) ;
% bndR = unique(triR.freeBoundary()) ;
% 
% % Triangulate the boundary of the cap with the endring via 
% % voronoi in yz plane 
% lring = 1:nU:(nU*(nV-1)) ;
% rring = nU:nU:(nU*(nV-1)) ;
% capL = [V(lring, :) * 2; lcap(bndL, :) ] ;
% capR = [V(rring, :) * 2; rcap(bndR, :) ] ;
% triL = delaunay(capL(:, 2), capL(:, 3)) ;
% triR = delaunay(capR(:, 2), capR(:, 3)) ;
% 
% % check it
% scatter(capL(:, 2), capL(:, 3), 'filled', 'MarkerFaceAlpha', 0.1, ...
%     'Markeredgecolor', 'none')
% triplot(triL, capL) ;
% 
% % Grab bulk endcap triangulation from triL and triR and remove from
% % delaunay version of endcap triangulation
% innerRingIDL = (length(lring)+1):size(capL, 1) ;
% innerRingIDR = (length(lring)+1):size(capR, 1) ;
% 
% keepFacesL = find(~all(ismember(triL, innerRingIDL), 2)) ;
% keepFacesR = find(~all(ismember(triR, innerRingIDR), 2)) ;
% triL = triL(keepFacesL, :) ;
% triR = triR(keepFacesR, :) ;
% error('conntinue here --> add face lists together')
% triL = lcapF

%% OPTION 2: Now add endcaps vertices (to be made virtual) -- spherical caps
closeIt = true ;
if closeIt
    points = vogelSphere(round(npts)) ;
    % points = [-1, 0, 0; 1, 0, 0] ;
    points = points * Rad ;

    % Divide the points into two hemispheres
    epsilon = 0.2 * Rad ;
    lcap = points(points(:, 1) < -epsilon, :) ;
    lcap(:, 1) = lcap(:,1) * 0.25 ;
    lcap = lcap - [Len*0.5,0,0] ;
    rcap = points(points(:, 1) > epsilon, :)  ;
    rcap(:, 1) = rcap(:,1) * 0.25 ;
    rcap = rcap + [Len*0.5,0,0] ;

    % Triangulate the boundary of the cap with the endring via 
    % voronoi in yz plane 
    lring = 1:nU:(nU*(nV-1)) ;
    rring = nU:nU:(nU*(nV-1)) ;
    lRingBig = [V(lring, 1), V(lring, 2)*1.2, V(lring, 3)*1.1] ;
    rRingBig = [V(rring, 1), V(rring, 2)*1.2, V(rring, 3)*1.1] ;
    capL = [lRingBig; lcap ] ;
    capR = [rRingBig; rcap ] ;
    % triL = Delaunay2_5D([capL(:, 2), capL(:, 3), capL(:, 1)]) ;
    % triR = Delaunay2_5D([capR(:, 2), capR(:, 3), capR(:, 1)]) ;
    triL = delaunay(capL(:, 2:3)) ;
    triR = delaunay(capR(:, 2:3)) ;
    
    Ln = faceNormals(triL, capL) ;
    Rn = faceNormals(triR, capR) ;
    if all(Ln(:, 1) > 0)
        disp('flipping Lcap')
        triL = [triL(:, 2), triL(:, 1), triL(:, 3)] ;
    elseif any(Ln(:, 1) > 0)
        error('inconsistent orientation')
    end
    
    if all(Rn(:, 1) < 0)
        disp('flipping Rcap')
        triR = [triR(:, 2), triR(:, 1), triR(:, 3)] ;
    elseif any(Rn(:, 1) < 0)
        error('inconsistent orientation')
    end

    % check it
    clf
    scatter3(capL(:, 1), capL(:, 2), capL(:, 3), 'filled', 'MarkerFaceAlpha', 0.1, ...
        'Markeredgecolor', 'none')
    hold on;
    trisurf(triL, capL(:, 1), capL(:, 2), capL(:, 3)) ;
    axis equal

    % replace tri cap pts with actual vertex IDs
    % first replace the vogel caps
    lupId = find(triL > length(lring)) ;
    rupId = find(triR > length(rring)) ;
    triL(lupId) = triL(lupId) + size(V, 1) - length(lring) ;
    triR(rupId) = triR(rupId) + size(V, 1) + size(lcap, 1) - length(rring) ;
    % Now replace the rings, starting from largest ID to avoid ambiguity
    for qq = fliplr(1:length(lring))
        triL(triL == qq) = lring(qq) ;
    end
    for qq = fliplr(1:length(rring))
        triR(triR == qq) = rring(qq) ;
    end

    % we should have this now closed
    VV = [V; lcap; rcap] ;
    FF = double([F; triL; triR]) ;
    capID = (nU*(nV-1)+1):length(VV) ;
    capLID = (nU*(nV-1)+1) + (0:(length(capL)-1)) ;
    capRID = ((nU*(nV-1)+1) + length(capL)):length(VV) ;

    % Construct Triangulation
    tri = triangulation(double(FF), VV) ;

    % preview the mesh
    clf
    fnormals = faceNormals(FF, VV) ;
    trisurf(tri, fnormals(:, 3))
    axis equal

    % Construct Topolgical Structure Tools ===================================
    [eIDx, feIDx, bulkEdgeIDx] = topologicalStructureTools(tri) ;

    % Check that closed (chi for sphere = 2) 
    chi = size(VV, 1) - size(eIDx, 1) + size(FF, 1) ;
    assert(chi == 2)

    % Save initial closed mesh
    FF = [FF(:, 2), FF(:, 1), FF(:, 3)] ;
    save(fullfile(outdir, 'initial_mesh_sphereCaps.mat'), 'VV', 'FF', 'eIDx', 'feIDx', 'bulkEdgeIDx')
else
    VV = V ;
    FF = [FF(:, 2), FF(:, 1), FF(:, 3)] ;
end
%% Construct Physical/Target Configurations and Geometries ================

% Construct Initial Configuration -----------------------------------------
% The initial configuration is a weakly buckled spherical cap.  This
% configuration is chosen to break the symmetry of the flat disk

% Directed edge vectors
eij = VV(eIDx(:,2), :) - VV(eIDx(:,1), :);

% Target edge lengths
eL = sqrt( sum( eij.^2, 2 ) );
eL0 = eL ; 

% % Compute the bond orientation angles
dx = VV(eIDx(:, 2)) - VV(eIDx(:, 1)) ;
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
% % Make figure of angle wrt axial direction beta as a check
% outfn = fullfile(outdir, 'wire_beta_definition.png') ;
% disp(['Saving ' outfn])
% aux_plot_beta_pattern(betaCut, Vc, eIDxCut, outfn) ;
outfn = fullfile(outdir, 'wire_beta_definition_glued.png') ;
disp(['Saving ' outfn])
close all
aux_plot_beta_pattern(beta, VV, eIDx, outfn) 

% Visualize the strain pattern
strainStruct.sigma = sigma ;
strainStruct.beta = beta ;
outfn = fullfile(outdir, 'wire_strain_pattern.png') ;
disp(['Saving ' outfn])

% Change the target edge lengths for each bond based on its midpt position
midx = 0.5 * (VV(eIDx(:, 1), 1) + VV(eIDx(:, 2), 1)) ;
midy = 0.5 * (VV(eIDx(:, 1), 2) + VV(eIDx(:, 2), 2)) ;
midz = 0.5 * (VV(eIDx(:, 1), 3) + VV(eIDx(:, 2), 3)) ;
midpt = [midx(:), midy(:), midz(:)] ;
clf
aux_plot_strain_pattern(VV, eIDx, eL0, distribution, strainstyle, strainStruct, -1, 1, outfn, Rad)


%% Initial face areas
a0 = 0.5 * doublearea(VV, FF) ;
V0 = VV ;
    
% Plotting options
cmap = bwr ;

% plot limits
xyzlim = [floor(min(V0)) - min(Rad, Len)*0.1; ...
    ceil(max(V0)) + min(Rad, Len)*0.1]' ;

%% save simulation parameters
targetTheta = 0 ; % either 0 or 'quasistatic'
save(fullfile(outdir, 'simulation_parameters.mat'), ...
    'targetTheta', 'distribution', 'strainstyle', 'strain', ...
    'fixBoundary', 'fixVolume', 'fixCap', 'V0', ...
    'thickness', 'poisson_ratio', 'capID', 'xyzlim', 'cutM')

%% increment strain at rate given by 'strain'
for ii = 1:uint8(1/strain)  
    % Calculate Target Edge Lengths -------------------------------------------
    % Determine how bonds are strained
    switch distribution
        case 'gaussian'
            mag = double(ii) * strain * exp(- midx.^2/ (2 * sigma^2)) ;
            switch strainstyle
                case 'isopulse'
                    titlestr = ['isotropic $\epsilon=$' sprintf('%2.2f', -double(ii) * strain) ];
                case 'hoopstrain'
                    titlestr = ['$\epsilon_{\phi\phi}=$' sprintf('%2.2f', -double(ii) * strain) ];
                case 'halfhoop'
                    titlestr = ['half-cylinder $\epsilon_{\phi\phi}=(-y/R)$' sprintf('%2.2f', double(ii) * strain) ];
                case 'linearhoop'
                    titlestr = ['graded cylinder $\epsilon_{\phi\phi}=[(R-y)/2R]$' sprintf('%2.2f', double(ii) * strain) ];
                case 'ringstrain'
                    titlestr = ['filament $\epsilon_{\phi\phi}=$' sprintf('%2.2f', double(ii) * strain) ];
            end            
            titlestr = [titlestr ' $\sigma=$' num2str(sigma)] ;
                    
    end
    
    switch strainstyle
        case 'isopulse'
            scalefactor = 1 - mag ;  
        case 'hoopstrain'
            scalefactor = 1 - mag .* abs(sin(beta)) ;
        case 'halfhoop'
            % target strain is zero at midplane and above, varies to mag at z=-Rad
            scalefactor = 1 + mag .* abs(sin(beta)) .* (midy/Rad);
            scalefactor(scalefactor > 1) = 1 ;
        case 'linearhoop'
            % target strain is zero at dorsal, varies to mag at z=-Rad
            scalefactor = 1 + mag .* abs(sin(beta)) .* (midy-Rad)/(2*Rad) ;
        case 'ringstrain'
            scalefactor = 1 - mag .* abs(sin(beta)) ;
            scalefactor(abs(sin(beta))<0.99) = 1 ;
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
        % debugger
        % [eLp, tarTheta] = calculateEdgeLengthsAndAngles(FF, V0);
        % eL = eLp + 1e-3 * rand(size(eL)) ;
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
    caxis([-0.1, 0.1])
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
    % % Extract face barycenters for prev & current timepoint
    % bc0 = barycenter(V0, FF) ;
    % bc = barycenter(VV, FF) ;
    % v0 = bc - bc0 ;
    % v0_vtx = VV - V0 ;
    % % Now resolve the vector field for decomposition
    % [v0n, v0t ] = ...
    %     resolveTangentNormalVelocities(FF, V0, v0, 1:length(FF)) ;
    % 
    % % Construct DEC class instance
    % DEC = DiscreteExteriorCalculus( FF, V0 ) ;
    % % Compute divergence and rotational flow
    % divv = DEC.divergence(v0t) ;
    % rotv = DEC.curl(v0t) ;
    % rotv(isnan(rotv)) = 0 ;
    
    % quiver scale magnitude
    % qscale = sqrt(Rad*Len) ./ max(vecnorm(v0_vtx, 2, 2)) ;
    % qsub = round(sqrt(nU*nV)*0.1) ;
    %
    % Plot div
    % opts.style = 'diverging' ;
    % opts.sscale = strain ;
    % opts.alpha = 1.0 ;
    % opts.labels = {'$\star$d$\star\left(v_t^\flat\right)$', ...
    %     '$\star$d$v_t^\flat$'} ;
    % opts.titles = {'$\star$d$\star\left(v_t^\flat\right)$', ...
    %     '$\star$d$v_t^\flat$'} ;
    % opts.qsubsample = qsub ;
    % opts.overlay_quiver = true ;
    % opts.qscale = qscale ;
    % opts.outfn = fullfile(outdir, sprintf('divcurl_%03d.png', ii)) ;
    % opts.figWidth = figWidth ;
    % opts.figHeight = figHeight ;
    % scalarVectorFieldsOnSurface(FF, V0, divv, ...
    %    V0(:, 1), V0(:, 2), V0(:, 3), ...
    %    v0_vtx(:, 1), v0_vtx(:, 2), v0_vtx(:, 3), opts)
    close all
end

%% Redo plots on a simulation
load(fullfile(outdir, 'simulation_parameters.mat'))
outdir = pwd ;
for ii = 1:uint8(1/strain) 
    %% Plot flow invariants
    load(fullfile(outdir, sprintf('vertices_%03d.mat', ii)), 'VV') ;

    % Extract face barycenters for this timepoint
    bc0 = barycenter(V0, FF) ;
    % Extract for prev timpepoint
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
