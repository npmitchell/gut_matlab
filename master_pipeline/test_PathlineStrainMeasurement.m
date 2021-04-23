%% Test Pathline Strain Measurement ===============================
% This is a script to exhibit the basic functionality of strain from
% Lagrangian frame. Check quasiconformal results too
%
% TO DO
% -----
% Check against analytic results for gaussian isopulse.
%

clear; close all; clc;
distribution = 'gaussian' ;
strainstyle = 'isopulse' ; % 'twist' 'isopulse' 'hoopstrain' 'halfhoop'
Rad = 5 ;       % radius of the tube
Len = 10 ;      % length of the tube
sigma = 0.5 ;     % width of the region to shrink
strain = 0.01 ; % strain magnitude per increment in time
delta = 0.1 ;
nU = 50 ;
nV = 50 ;

% Figure Options
figWidth = 20 ; 
figHeight = 12 ;
cmin = 0.95 ;
cmax = 1.05 ;
pm256 = phasemap(256) ;
clim_deviatoric = 0.5 ;

% Add paths
gutDir = '/mnt/data/code/gut_matlab/' ;
addpath(fullfile(gutDir, 'addpath_recurse')) ;
addpath_recurse(fullfile(gutDir, 'mesh_handling')) ;
addpath_recurse(fullfile(gutDir, 'plotting')) ;


% Path options
NESpath = '/mnt/data/code/NonEuclideanShells/NES/' ;
addpath_recurse(NESpath)
exten = sprintf('_L%02dR%02d_sigma%0.3f_strain%0.3f_%03dx%03d', ...
    Len, Rad, sigma, strain, nU, nV) ;
exten = strrep(exten, '.', 'p') ;
dirname = [strainstyle exten '_test'] ;
outdir = ['/mnt/data/simulations/strain_tube_' strainstyle '/'] ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

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
refMesh.nU = nU; 
refMesh.nV = nV ;
refMesh.u = sp ;
refMesh.f = faces ;
refMesh.v = vcut ;
refMesh.pathPairs = [ (1:nU)', (nV-1)*nU + (1:nU)' ] ;

gMesh = glueCylinderCutMeshSeam(refMesh) ;
V0 = refMesh.v ;
F0 = refMesh.f ;

% Construct Triangulation
tri = triangulation(F0, V0) ;

%% Compute strain for test case of two meshes with same lagrangian frame
% Load Lagrangian advected vertices as mesh for this timepoint

% Make indentation in cylinder
amps = 0:0.01:0.75 ;
sigma = 1 ;
areas = zeros(length(amps), 1) ;
tres = 0*areas ;
devs = 0*areas ;
for dmyk = 1:length(amps)
    amp = amps(dmyk) ;
    mesh = struct() ;
    mesh.f = refMesh.f ;
    V = V0 ;
    
    if strcmpi(strainstyle, 'isopulse')
        indent = amp * exp(-0.5 * V(:, 1).^2 ./ sigma^2) ;
        V(:, 2) = V(:, 2) .* (1- indent) ;
        V(:, 3) = V(:, 3) .* (1- indent) ;
    elseif strcmpi(strainstyle, 'twist')        
        indent = 0 * V(:, 1).^2 ;
        sp = [ss, phi(:) + amp * ss / Rad ] ;
        V = [ss, Rad * cos(sp(:, 2)), Rad* sin(sp(:, 2))] ;
    end
    mesh.v = V ;
    mesh.u = refMesh.u ;
    mesh.pathPairs = refMesh.pathPairs ;
    mesh.nU = nU ;
    mesh.nV = nV ;

    options = struct() ;
    [strain, tre, dev, theta, outputStruct] = ...
        inducedStrainPeriodicMesh(mesh, refMesh, options) ;

    % for debugging:
    tileCount = [1, 1] ;
    [ TF0, TV2D0, TV3D0] = tileAnnularCutMesh( refMesh, tileCount ) ;
    
    % Face averaging operators
    [V2F, F2V] = meshAveragingOperators(refMesh.f, mesh.v) ;

    % The fields of outputStruct are:
    %   fundForms : struct with fields
    %       gg : first fundamental form on each face for tiled mesh
    %       bb : second fundamental form on each face for tiled mesh
    %   bondDxDy : struct with fields
    %       dx : length of dx in embedding space / length of dx in pullback space
    %       dy : length of dy in embedding space / length of dy in pullback space
    %       dxTiled : length of dx in embedding space / length of dx in
    %           pullback space for tiled mesh
    %       dyTiled : length of dy in embedding space / length of dy in
    %           pullback space for tiled mesh
    %   theta_pb : #tiledFaces x 1 float array
    %       elongation angle, theta, in pullback space (differs from
    %       theta by scaling of the eigenvectors by (dx, dy) returned in
    %       bondDxDy
    %   faceIDs : #faces x 1 int array
    %       indices of tiled faces that return the faces of the
    %       original input mesh, so that strain_orig = strain(faceIDs)
    fundForms = outputStruct.fundForms ;
    bondDxDy = outputStruct.bondDxDy ;
    theta_pb = outputStruct.theta_pb ;
    faceIDs = outputStruct.faceIDs ;
    tre = tre(faceIDs) ;
    dev = dev(faceIDs) ;
    theta = theta(faceIDs) ;
    strain = strain(faceIDs, :) ;

    % qopts : struct with fields
    %   style : str ('diverging' 'phase', default='diverging')
    %       style of overlaid scalar field
    %   sscale : float
    %       scalar field maximum for clims/opacity
    %   alpha : float (optional, used if style ~= 'phase')
    %       opacity of overlaid scalar field
    %   outfn : str
    %       path to save image if given
    %   label : str
    %       colorbar label. Default is '$|v|$ [$\mu$m / min]' 
    %   outfn : str
    %       output filename for figure as png 
    %   figWidth : int (optional, default = 16) 
    %       figure width in cm
    %   figHeight : int (optional, default = 10) 
    %       figure height in cm
    %   cbarPosition : 4 x 1 float array
    %       position of the colorbar
    if dmyk == 1
        area0 = 0.5 * doublearea(mesh.v, mesh.f) ;
    end
    area1 = 0.5 * doublearea(mesh.v, mesh.f) ;
    Darea = double(area1 - area0) ./ double(area0) ;
    
    
    fn = fullfile(outdir, [sprintf('strain_%05d', dmyk) '.png']) ;    
    if ~exist(fn, 'file') || dmyk == 1
        options = struct() ;
        % options.clim_trace = 0.1 ;
        set(gcf, 'visible', 'off')
        options.labels =  {'$\frac{1}{2}\mathrm{Tr} [\bf{g}^{-1}\varepsilon] $', ...
         '$||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\varepsilon\right]\bf{g}||$'} ;
        options.edgecolor = 'black' ;
        [ax1, ax2, cb1, cb2, pbar, meshTr, meshDev] = ...
            plotTraceDeviator(mesh, tre, dev, theta, options) ;
        if dmyk == 1
            set(gcf,'CurrentAxes',ax1)
            xlims = xlim() ;
            ylims = ylim() ;
            zlims = zlim() ;
            axis off
            view(2)
            set(gcf,'CurrentAxes',ax2)
            axis off
            view(2)
        else
            set(gcf,'CurrentAxes',ax1)
            xlim(xlims) ;
            ylim(ylims) ;
            zlim(zlims) ;
            axis off
            view(2)
            set(gcf,'CurrentAxes',ax2)
            xlim(xlims) ;
            ylim(ylims) ;
            zlim(zlims) ;
            axis off
            view(2)
        end
        sgtitle(['$\delta/R = $' sprintf('%0.2f', amp)], 'interpreter', 'latex')

        disp(['Saving ' fn])
        saveas(gcf, fn)
        close all
    end
    
    %% Save image of change in area of triangles
    fn = fullfile(outdir, [sprintf('areas_%05d', dmyk) '.png']) ;    
    if ~exist(fn, 'file')
        options = struct() ;
        % options.clim_trace = 0.1 ;
        set(gcf, 'visible', 'off')
        options.labels =  {'$\Delta A / A_0 $', ...
         '$\frac{1}{2}\mathrm{Tr} [\bf{g}^{-1}\varepsilon] $'} ;
        [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
            twoScalarFieldsOnSurface(mesh, Darea, tre, options) ;
        set(gcf,'CurrentAxes',ax1)
        xlim(xlims) ;
        ylim(ylims) ;
        zlim(zlims) ;
        axis off
        view(2)
        set(gcf,'CurrentAxes',ax2)
        xlim(xlims) ;
        ylim(ylims) ;
        zlim(zlims) ;
        axis off
        view(2)
        sgtitle(['$\delta/R = $' sprintf('%0.2f', amp)], 'interpreter', 'latex')

        disp(['Saving ' fn])
        saveas(gcf, fn)
        close all
    end

    areas(dmyk) = sum(area1 - area0) ./ sum(area0) ;
    tres0(dmyk) = sum(0.5 * tre .* doublearea(refMesh.v, refMesh.f)) / sum(area0) ;
    tres(dmyk) = sum(0.5 * tre .* doublearea(mesh.v, mesh.f)) / sum(area1) ;
    devs(dmyk) = sum(0.5 * dev .* doublearea(mesh.v, mesh.f)) / sum(area1) ;
    theta_cos = nansum(cos(2*theta(:)) .* dev / sum(dev)) ;
    theta_sin = nansum(sin(2*theta(:)) .* dev / sum(dev)) ;
    if theta_cos == 0 && theta_sin == 0
        thetas(dmyk) = NaN ;
    else
        thetas(dmyk) = 0.5 * atan2(theta_sin, theta_cos) ;
    end
        
    %% Test quasiconformal measurement
    % First flatten annulus
    cMesh = flattenAnnulus(mesh) ;      
    % Relax Affine transformation along x
    cMesh.ar = minimizeIsoarealAffineEnergy( cMesh.f, cMesh.v, cMesh.u );
    % Define material frame
    if dmyk == 1
        ar_material = cMesh.ar ;
        refMesh = cMesh ;    
    end
    cMesh.u_affine = cMesh.u ;
    cMesh.u_affine(:, 1) = cMesh.u(:, 1) * cMesh.ar ;
    mu_m = bc_metric(mesh.f, refMesh.u, cMesh.v, 3) ;
    if dmyk == 1
        mu0 = mean(real(mu_m)) ;
    end
    mu_m = mu_m - mu0 ;
    
    %% plot 3d -- affine relaxation map
    fn3d = fullfile(outdir, [sprintf('mu3d_affine_%05d', dmyk) '.png']) ;   
    mu_a = bc_metric(mesh.f, cMesh.u_affine, cMesh.v, 3) ;
    mu_a = mu_a - mean(real(mu_a)) ;
    labels = {'$\Re \mu$', '$\Im \mu$'} ;
    options.labels = labels ;
    [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
        twoScalarFieldsOnSurface({cMesh.f, cMesh.v}, ...
        real(mu_a), imag(mu_a), options) ;
    sgtitle(['conformal pullback: $\delta/R = $' sprintf('%0.2f', amp)], ...
        'interpreter', 'latex')
    disp(['Saving ' fn3d])
    saveas(gcf, fn3d)
    close all
    
    %% plot 3d -- material frame
    fn3d = fullfile(outdir, [sprintf('mu3d_material_%05d', dmyk) '.png']) ;   
    labels = {'$\Re \mu$', '$\Im \mu$'} ;
    options.labels = labels ;
    [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
        twoScalarFieldsOnSurface({cMesh.f, cMesh.v}, ...
        real(mu_m), imag(mu_m), options) ;
    sgtitle(['material pullback: $\delta/R = $' sprintf('%0.2f', amp)], ...
        'interpreter', 'latex')
    disp(['Saving ' fn3d])
    saveas(gcf, fn3d)
    close all
    
    %% plot in 2d -- affine relaxation
    fn2d = fullfile(outdir, [sprintf('mu2d_affine_%05d', dmyk) '.png']) ;  
    labels = {'$\Re \mu$', '$\Im \mu$'} ;
    options.labels = labels ;
    [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
        twoScalarFieldsOnSurface({cMesh.f, [cMesh.u_affine(:, 1), ...
        cMesh.u(:, 2), 0 * cMesh.u(:, 2)]}, ...
        real(mu_a), imag(mu_a), options) ;
    sgtitle(['conformal pullback: $\delta/R = $' sprintf('%0.2f', amp)], ...
        'interpreter', 'latex')
    set(gcf,'CurrentAxes', ax1)
    view(2)
    set(gcf,'CurrentAxes', ax2)
    view(2)
    disp(['Saving ' fn2d])
    saveas(gcf, fn2d)
    close all
        
    %% plot in 2d -- Lagrangian frame
    fn2d = fullfile(outdir, [sprintf('mu2d_material_%05d', dmyk) '.png']) ;   
    labels = {'$\Re \mu$', '$\Im \mu$'} ;
    options.labels = labels ;
    [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
        twoScalarFieldsOnSurface({cMesh.f, [cMesh.u(:, 1), ...
        cMesh.u(:, 2), 0 * cMesh.u(:, 2)]}, ...
        real(mu_m), imag(mu_m), options) ;
    sgtitle(['material pullback: $\delta/R = $' sprintf('%0.2f', amp)], ...
        'interpreter', 'latex')
    set(gcf,'CurrentAxes', ax1)
    view(2)
    set(gcf,'CurrentAxes', ax2)
    view(2)
    disp(['Saving ' fn2d])
    saveas(gcf, fn2d)
    close all
    
    %% Nonlinear mu from isothermal coords 
    % This is like QS.compareBeltramiToConstriction()
    % Grab dzeta and dphi from 
    
    % Check that coords are isothermal
    metricSPhiGridMesh
    
    
    %% Linearized mu approximation to indentation
    fn2d = fullfile(outdir, [sprintf('mu2d_indentation_%05d', dmyk) '.png']) ; 
    % labels = {'$\Re \mu$', '$\rho/2$', '$\frac{\Re\mu-\rho/2}{\langle |\Re \mu| \rangle}$'} ;
    labels = {'$\Re \mu$', '$\rho/2$', '$\Re\mu-\rho/2$'} ;
    options = struct() ;
    options.labels = labels ;
    options.clim = max(max(real(mu_m)), max(indent * 0.5)) ;
    error_mu = real(mu_m) - V2F * indent*0.5 ;
    % if mean(abs(real(mu_m))) > 1e-14
    %     error_mu = error_mu ./ max(abs(real(mu_m)), abs(V2F * indent*0.5)) ;
    % end
    % options.clim3 = 0.1 ;
    options.edgecolor = 'none' ;
    [axs, cbs, meshes] = ...
        threeScalarFieldsOnSurface({cMesh.f, [cMesh.u(:, 1), ...
        cMesh.u(:, 2), 0 * cMesh.u(:, 2)]}, ...
        real(mu_m), indent * 0.5, error_mu, options) ;
    sgtitle(['$\mu=\frac{\rho}{2}$ as linear approximation: $\delta/R = $', ...
        sprintf('%0.2f', amp)], 'interpreter', 'latex')
    set(gcf,'CurrentAxes', axs(1))
    view(2)
    axis off
    set(gcf,'CurrentAxes', axs(2))
    view(2)
    axis off    
    set(gcf,'CurrentAxes', axs(3))
    view(2)
    axis off 
    disp(['Saving ' fn2d])
    saveas(gcf, fn2d)
    close all
    
    %% Smooth a little bit
    
    
%     
%     %% Linearized mu approximation to indentation
%     fn2d = fullfile(outdir, [sprintf('mu2d_6panel_%05d', dmyk) '.png']) ; 
%     labels = {'$\Re \mu$', '$\rho/2$', '$\frac{\Re\mu-\rho/2}{\langle |\Re \mu| \rangle}$', ...
%         '$\Delta A / A_0$', ...
%         '$\frac{1}{2}\mathrm{Tr} [\bf{g}^{-1}\varepsilon] $', ...
%         '$||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\varepsilon\right]\bf{g}||$'} ;
%     options.labels = labels ;
%     options.clim = max(max(real(mu_m)), max(indent * 0.5)) ;
%     error_mu = real(mu_m) - V2F * indent*0.5 ;
%     if mean(abs(real(mu_m))) > 1e-14
%         error_mu = error_mu / mean(abs(real(mu_m))) ;
%     else
%         
%     end
%     c2d = {cMesh.f, [cMesh.u(:, 1), cMesh.u(:, 2), 0 * cMesh.u(:, 2)]} ;
%     c3d = {cMesh.f, [cMesh.v(:, 1), cMesh.v(:, 2), cMesh.v(:, 3)]} ;
%     options.edgecolor = 'none' ;
%     [axs, cbs, meshes] = ...
%         nFieldsOnSurface({c2d,c2d,c2d,c3d,c3d,c3d}, ...
%         {real(mu_m), indent * 0.5, error_mu, ...
%         Darea, tre, {dev,theta}}, options) ;
%     sgtitle(['$\mu=\frac{\rho}{2}$ as linear approximation: $\delta/R = $', ...
%         sprintf('%0.2f', amp)], 'interpreter', 'latex')
%     set(gcf,'CurrentAxes', axs{1})
%     view(2)
%     axis off
%     set(gcf,'CurrentAxes', axs{2})
%     view(2)
%     axis off    
%     set(gcf,'CurrentAxes', axs{3})
%     view(2)
%     axis off 
%     disp(['Saving ' fn2d])
%     saveas(gcf, fn2d)
     

%     %% Lineared equations through derivative of indentation wrt z 
%     % Laplacian(Ru) = d_z rho
%     % d_zeta Iu = - d_phi Ru
%     % Lambda = 1 - 2 d_z Ru
%     
%
%     % Solve the Poisson Equation =============================================
%     % Construct Differential Operators ----------------------------------------
% 
%     % A DEC object for the current mesh
%     DEC = DiscreteExteriorCalculus( F, [ V, zeros(size(V,1), 1) ] );
% 
%     % The full mesh vertex 'mass' operator
%     M = DEC.hd0;
% 
%     % The full mesh (unweighted) Laplacian matrix
%     % ( The weighted Laplacian L = M^(-1) * C )
%     C = DEC.dd1 * DEC.hd1 * DEC.d0;
% 
%     % NOTE: IT IS REQUIRED THAT ALL THE BOUNDARY VERTICES BE AT THE END OF
%     % THE VERTEX COORDINATE LIST (THIS SHOULD BE DONE DURING MESH
%     % CONSTRUCTION)
% 
%     % The 'mass' operator for interior vertices
%     MII = M(1:(min(bdyIDx)-1), 1:(min(bdyIDx)-1));
% 
%     % The (unweighted) Laplacian for interior vertices
%     CII = C(1:(min(bdyIDx)-1), 1:(min(bdyIDx)-1));
% 
%     % The mixed (unweighted) Laplacian
%     CIB = C(1:(min(bdyIDx)-1), min(bdyIDx):end);
% 
%     % Construct the Poisson Kernel --------------------------------------------
% 
%     % The kernel for interior vertices
%     GII = G((1:(min(bdyIDx)-1)).');
% 
%     % Set Dirichlet Boundary Conditions ---------------------------------------
% 
%     % The function values for boundary vertices
%     UB = U(bdyIDx);
% 
%     %--------------------------------------------------------------------------
%     % Solve the Poisson Problem
%     %--------------------------------------------------------------------------
% 
%     % Calculate the velocity components
%     calcU = CII \ (MII * GII - CIB * UB);
% 
%     % Calculate solution residuals (checks if any solution was found - not that
%     % the solution is correct)
%     solveRes = ( CII * calcU - (MII * GII - CIB * UB) );
% 
%     % Include known boundary velocities
%     calcU = [ calcU; UB ];
% 
%     fprintf('Maximum Solution Residual = %0.5e\n', max(abs(solveRes)));
%     

end

% scale = 25 ;
clf
plot(amps, areas)
ylabel('surface area $(A-A_0)/A_0$', 'interpreter', 'latex')
%ylim([min(min(areas), min(scale*tres)), max(max(areas), max(scale*tres))])
yyaxis right
plot(amps, tres)
%ylim([min(min(areas), min(scale*tres)), max(max(areas), max(scale*tres))])
ylabel('$\langle \mathrm{Tr} \varepsilon \rangle$', 'interpreter', 'latex')
xlabel('deformation, $\delta/R$', 'interpreter', 'latex')
fn = fullfile(outdir, ['surfacearea_trace.png']) ;   
saveas(gcf, fn)

% Save deviator
clf
plot(amps, thetas)
ylabel('shear angle, $\langle\theta \rangle$', 'interpreter', 'latex')
yyaxis right
plot(amps, devs)
ylabel('$\langle \mathrm{Dev} \varepsilon \rangle$', 'interpreter', 'latex')
xlabel('deformation, $\delta/R$', 'interpreter', 'latex')
fn = fullfile(outdir, ['deviator.png']) ;   
saveas(gcf, fn)