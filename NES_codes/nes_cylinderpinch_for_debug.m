%% NonEuclidean_Shell for pinching cylinder ===============================
% This is a script to exhibit the basic functionality of the
% (N)on-(E)uclidean (S)hell (S)imulator (NESS) using a cylindrical geometry
% which is pinched along its midline.

clear; close all; clc;
distribution = 'gaussian' ;
strainstyle = 'halfhoop' ; % 'isopulse' 'hoopstrain' 'halfhoop'
fixVolume = false ;
fixBoundary = false ;
fixCap = false ;
Rad = 50 ;       % radius of the tube
Len = 20 ;      % length of the tube
sigma = 10 ;     % width of the region to shrink
strain = 0.01 ; % strain magnitude per increment in time
nU = 2 ;
nV = 4 ;
npts = round((nU / (2*pi))^2 * 4*pi) ;

% Figure Options
figWidth = 20 ; 
figHeight = 12 ;
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
exten = sprintf('_L%02dR%02d_sigma%0.3f_strain%0.3f_%03dx%03d', ...
    Len, Rad, sigma, strain, nU, nV) ;
exten = strrep(exten, '.', 'p') ;
if fixBoundary
    exten = ['_fixB' exten ] ;
end
if fixVolume
    exten = ['_fixV' exten ] ;
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
sp = [ss, phi(:)] ;
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
%points = vogelSphere(round(npts)) ;
points = [-1, 0, 0; 1, 0, 0] ;
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

% Construct Triangulation
tri = triangulation(double(FF), VV) ;

% preview the mesh
clf
trisurf(tri)
axis equal

%% Construct Topolgical Structure Tools ===================================
[eIDx, feIDx, bulkEdgeIDx] = topologicalStructureTools(tri) ;

% Check that closed (chi for sphere = 2) 
chi = size(VV, 1) - size(eIDx, 1) + size(FF, 1) ;
assert(chi == 2)

%% Save initial closed mesh
save(fullfile(outdir, 'initial_mesh_sphereCaps.mat'), 'VV', 'FF', 'eIDx', 'feIDx', 'bulkEdgeIDx')

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

% Change the target edge lengths for a ring
midx = 0.5 * (VV(eIDx(:, 1), 1) + VV(eIDx(:, 2), 1)) ;
midy = 0.5 * (VV(eIDx(:, 1), 2) + VV(eIDx(:, 2), 2)) ;
midz = 0.5 * (VV(eIDx(:, 1), 3) + VV(eIDx(:, 2), 3)) ;
midpt = [midx(:), midy(:), midz(:)] ;

clf
aux_plot_strain_pattern(VV, eIDx, eL0, distribution, strainstyle, strainStruct, -1, 1, outfn, Rad)


%% Initial face areas
a0 = 0.5 * doublearea(VV, FF) ;

% Plotting options
cmap = bwr ;

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
    end

    eL = eL .* scalefactor ;

    % Calculate Target Bending Angles -----------------------------------------
    % The unit normal vector for each face
    % see https://www.cs.utexas.edu/users/evouga/uploads/4/5/6/8/45689883/turning.pdf
    % fN = faceNormal( triangulation(F, tarV) );
    % 
    % N1 = fN(edgeFace(:,1), :);
    % N2 = fN(edgeFace(:,2), :);
    % 
    % crossN = cross(N2, N1, 2);
    % dotN = dot(N1, N2, 2);
    % 
    % tarTheta = 2 .* atan2( dot(crossN, eij, 2), 1 + dotN );
    tarTheta = zeros(size(eL, 1), 1) ; 

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
        VV = minimizeElasticEnergy( FF, VV, eL, ...
            'TargetAngles', tarTheta, ...
            'Thickness', 0.1, ...
            'Poisson', 0.3, ...
            'MaxIterations', 1000, ...
            'iterDisplay', 1, ...
            'Alpha', 1, ...
            'Beta', 1, ...
            'FixBoundary', ...
            'targetVertices', capID(:), ...
            'targetLocations', VV(capID(:), :), ...
            'FixVolume', 'TargetVolume', targetVolume); 
    elseif fixBoundary && fixCap
        VV = minimizeElasticEnergy( FF, VV, eL, ...
            'TargetAngles', tarTheta, ...
            'Thickness', 0.1, ...
            'Poisson', 0.3, ...
            'MaxIterations', 1000, ...
            'iterDisplay', 1, ...
            'Alpha', 1, ...
            'Beta', 1, ...
            'targetVertices', capID(:), ...
            'targetLocations', VV(capID(:), :), ...
            'FixBoundary');    
    elseif fixCap && fixVolume
        VV = minimizeElasticEnergy( FF, VV, eL, ...
            'TargetAngles', tarTheta, ...
            'Thickness', 0.1, ...
            'Poisson', 0.3, ...
            'MaxIterations', 1000, ...
            'iterDisplay', 1, ...
            'Alpha', 1, ...
            'Beta', 1, ...
            'targetVertices', capID(:), ...
            'targetLocations', VV(capID(:), :), ...
            'FixVolume', 'TargetVolume', targetVolume);    
    elseif fixBoundary || true
        VV = minimizeElasticEnergy( FF, VV, eL, ...
            'TargetAngles', tarTheta, ...
            'Thickness', 0.1, ...
            'Poisson', 0.3, ...
            'MaxIterations', 1000, ...
            'iterDisplay', 1, ...
            'Alpha', 1, ...
            'Beta', 1, ...
            'FixBoundary');    
    elseif fixVolume
        VV = minimizeElasticEnergy( FF, VV, eL, ...
            'TargetAngles', tarTheta, ...
            'Thickness', 0.1, ...
            'Poisson', 0.3, ...
            'MaxIterations', 1000, ...
            'iterDisplay', 1, ...
            'Alpha', 1, ...
            'Beta', 1, ...
            'FixVolume', 'TargetVolume', targetVolume);    
    else
        error('not fixing Volume or boundary?')
    end
    toc

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
    
    %% Plot flow from previous timepoint
    prev = load(fullfile(outdir, sprintf('vertices_%03d', ii)), 'VV') ;

    % Extract face barycenters for this timepoint
    bcp = barycenter(prev.VV, FF) ;
    % Extract for prev timpepoint
    bc = barycenter(VV, FF) ;

    v0 = bc - bcp ;
    [v0n, v0t ] = ...
        resolveTangentNormalVelocities(FF, VV, v0, 1:length(FF)) ;

    % Plot phase info of 
    opts.style = 'phase' ;
    opts.sscale = strain ;
    opts.alpha = 1.0 ;
    opts.label = '$v$ [$\mu$m / min]' ;
    opts.qsubsample = 10 ;
    opts.overlay_quiver = true ;
    opts.qscale = 10 ;
    opts.outfn = fullfile(outdir, sprintf('vphase_%03d', ii)) ;
    opts.figWidth = figWidth ;
    opts.figHeight = figHeight ;
    scalarVectorFieldsOnSurface(FF, VV, sf, xxv, yyv, zzv, vx,vy,vz, opts)

    % Perform DEC analysis of these fields
    FF = mesh.f ;  % face connectivity list of the mesh
    VV = mesh.v ;  % vertex positions in 3D of the CLOSED mesh

    vel3D = VV - V0 ;
    % Interpolate vel3D onto facebarycenters so we have one vector for each
    % face
    bc = barycenter(VV, FF) ;
    v3Dxi = griddedInterpolant(x0, y0, vel3D(:, 1)) ;
    v3Dyi = griddedInterpolant(x0, y0, vel3D(:, 2)) ;
    v3Dzi = griddedInterpolant(x0, y0, vel3D(:, 3)) ;
    vel3Df = zeros(size(FF, 1), 3) ;
    vel3Df(:, 1) = v3Dxi(bc(:,1), bc(:, 2)) ;
    vel3Df(:, 2) = v3Dyi(bc(:,1), bc(:, 2)) ;
    vel3Df(:, 3) = v3Dzi(bc(:,1), bc(:, 2)) ;

    % Construct DEC class instance
    DEC = DiscreteExteriorCalculus( FF, VV ) ;

    % Now resolve the vector field for decomposition
    [v0n, v0t, v0t2d, jac3d_to_2d, ~, ~, dilation] = ...
        resolveTangentNormalVelocities(FF, VV, vel3Df, 1:size(FF, 1), V2D ) ;

    % Compute divergence and rotational flow
    divv = DEC.divergence(v0t) ;
    rotv = DEC.curl(v0t) ;
    
end
