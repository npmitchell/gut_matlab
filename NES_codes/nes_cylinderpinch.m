%% NonEuclidean_Shell for pinching cylinder ===============================
% This is a script to exhibit the basic functionality of the
% (N)on-(E)uclidean (S)hell (S)imulator (NESS) using a cylindrical geometry
% which is pinched along its midline.

clear; close all; clc;
distribution = 'gaussian' ;
strainstyle = 'halfhoop' ; % 'isopulse' 'hoopstrain' 'halfhoop'
Rad = 5 ;       % radius of the tube
Len = 10 ;      % length of the tube
sigma = 0.5 ;     % width of the region to shrink
strain = 0.01 ; % strain magnitude per increment in time
nU = 50 ;
nV = 50 ;

% Figure Options
figWidth = 20 ; 
figHeight = 12 ;
cmin = 0.95 ;
cmax = 1.05 ;

% Path options
NESpath = '/mnt/data/code/NonEuclideanShells/NES/' ;
addpath_recurse(NESpath)
exten = sprintf('_L%02dR%02d_sigma%0.3f_strain%0.3f_%03dx%03d', ...
    Len, Rad, sigma, strain, nU, nV) ;
exten = strrep(exten, '.', 'p') ;
dirname = [strainstyle exten ] ;
outdir = fullfile('/mnt/data/simulations/NES_cylinder/', dirname) ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Add paths
addpath('/mnt/data/code/gut_matlab/addpath_recurse/')
addpath_recurse('/mnt/data/code/gut_matlab/mesh_handling/')
addpath_recurse('/mnt/data/code/gut_matlab/plotting/')

%--------------------------------------------------------------------------
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

%% Construct Topolgical Structure Tools ===================================
[eIDx, feIDx, bulkEdgeIDx] = topologicalStructureTools(tri) ;


%% Construct Physical/Target Configurations and Geometries ================

% Construct Initial Configuration -----------------------------------------
% The initial configuration is a weakly buckled spherical cap.  This
% configuration is chosen to break the symmetry of the flat disk

% Directed edge vectors
eij = V(eIDx(:,2), :) - V(eIDx(:,1), :);

% Target edge lengths
eL = sqrt( sum( eij.^2, 2 ) );
eL0 = eL ; 

% % Compute the bond orientation angles
dx = ss(eIDx(:, 2)) - ss(eIDx(:, 1)) ;
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
aux_plot_beta_pattern(beta, V, eIDx, outfn) 

% Visualize the strain pattern
strainStruct.sigma = sigma ;
strainStruct.beta = beta ;
outfn = fullfile(outdir, 'wire_strain_pattern.png') ;
disp(['Saving ' outfn])
close all
aux_plot_strain_pattern(V, eIDx, eL0, distribution, strainstyle, strainStruct, -1, 1, outfn)

% boundary
bnd = tri.freeBoundary ;
bnd = bnd(:, 1) ;
% left boundary
left = find(V(bnd, 1) < 0) ;
right = find(V(bnd, 1) > 0) ;

% Initial face areas
a0 = 0.5 * doublearea(V, F) ;

% Change the target edge lengths for a ring
midpts = 0.5 * (V(eIDx(:, 1), 1) + V(eIDx(:, 2), 1)) ;
midy = 0.5 * (V(eIDx(:, 1), 2) + V(eIDx(:, 2), 2)) ;
midz = 0.5 * (V(eIDx(:, 1), 3) + V(eIDx(:, 2), 3)) ;

% Plotting options
cmap = bwr ;

%%
for ii = 1:uint8(1/strain)  
    
    % Calculate Target Edge Lengths -------------------------------------------

    % Determine how bonds are strained
    switch distribution
        case 'gaussian'
            mag = double(ii) * strain * exp(- midpts.^2/ (2 * sigma^2)) ;
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
            % strain is zero at midplane and above, varies to mag at z=-Rad
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


    %% Run Elastic Relaxation =================================================

    tic
    V = minimizeElasticEnergy( F, V, eL, ...
        'TargetAngles', tarTheta, ...
        'Thickness', 0.1, 'Poisson', 0.3, ...
        'MaxIterations', 1000, ...
        'FixBoundary', 'Alpha', 10 );
%     
%     V = minimizeElasticEnergy( F, V, eL, ...
%         'TargetAngles', tarTheta, ...
%         'Thickness', 0.1, 'Poisson', 0.3, ...
%         'MaxIterations', 1000, ...
%         'TargetVertices', left, ...
%         'TargetLocations', V(left, :), ...
%         'Alpha', 10 );
    
    toc

    %% View Results ===========================================================
    close all
    figure('visible', 'off')
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 figWidth figHeight]);  
        
    % Area ratio for faces
    areas = 0.5 * doublearea(V, F) ;
    astrain = (areas - a0) ./ a0 ;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot faces
    h = trisurf(triangulation(F, V), astrain, 'edgecolor', 'none');
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
    lsegs = [V(eIDx(:,1), :), V(eIDx(:, 2), :)] ;
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
    save(fullfile(outdir, sprintf('vertices_%03d', ii)), 'V', 'F')
    
    % Plot flow from previous timepoint
    if ii > 1
        prev = load(fullfile(outdir, sprintf('vertices_%03d', ii)), 'V') ;
        
        % Extract face barycenters for this timepoint
        bcp = barycenter(prev.V, F) ;
        % Extract for prev timpepoint
        bc = barycenter(V, F) ;
        
        v0 = bc - bcp ;
        [v0n, v0t ] = ...
            resolveTangentNormalVelocities(faces, vertices, v0, 1:length(F))
        
        
        
        % Plot phase info of 
        opts.style = 'phase' ;
        opts.sscale = strain ;
        opts.alpha = 1.0 ;
        opts.label = '$v$ [$\mu$m / min]' ;
        opts.qsubsample = 10 ;
        opts.overlay_quiver = true ;
        opts.qscale = 10 ;
        opts.outfn = fullfile(outdir, printf('vphase_%03d', ii)) ;
        opts.figWidth = figWidth ;
        opts.figHeight = figHeight ;
        scalarVectorFieldsOnSurface(F, V, sf, xxv, yyv, zzv, vx,vy,vz, opts)
        error('here')
    end
    
end
