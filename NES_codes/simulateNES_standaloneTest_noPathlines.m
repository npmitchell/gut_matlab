% DEBUGGING script for function simulateNES(QS, options)
% function simulateNES_standaloneTest_noPathlines()
%
% 1 FOLLOW (s,phi) MESHES --> fix points to new points, fix volume to new vol
% 2 Check that pathline meshes have good triangle quality
% 3 FOLLOW PATHLINE MESHES
% 4 relaxing some contraints 
%
% Try 
% targetEdgeLength = 10 ;
% numIterations = 8 ;
% [newFaces, newVertices, newFaceNormals, newVertexNormals] = ...
%       isotropic_remeshing(faces, vertices, targetEdgeLength, numIterations)
% 
% Try 
% [intx_exist, faces_intx ] = mesh_self_intersection_3d(faces, vertices)
%
% Try
% For pathlines, load smoothed 3d velocity for all vertices/faces, apply to
% current (remeshed) vertices, compute those target lengths/bending angles.
%
% Track triangle quality for pathline meshes
% measure angles for quality inspection
% gptoolbox -->  internalangles( V, F )
% 

clear all

QS = struct() ;
QS.nU = 100 ;
QS.nV = 100 ;
QS.t0 = 123 ;
QS.timeInterval = 1 ;
QS.timeUnits = 'min' ;
options = struct() ;

%% Add path (todo: make this automatic by putting NES code in gut_matlab)

% VIP 8 
codedir = '/mnt/data/code/'; 
NESpath = '/mnt/data/code/NonEuclideanShells/NES/' ; % Debug

% Laptop
% codedir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/'; 
% NESpath = fullfile(codedir, 'NES') ;

% Run once
addpath_recurse(NESpath)
addpath_recurse(fullfile(codedir, 'DEC'))
addpath(fullfile(codedir, 'gut_matlab','NES_codes'))
addpath_recurse(fullfile(codedir, 'gut_matlab', 'mesh_handling'))
addpath_recurse(fullfile(codedir, 'gut_matlab', 'geometry'))
addpath_recurse(fullfile(codedir, 'gut_matlab', 'plotting'))
addpath_recurse(fullfile(codedir, 'gut_matlab', 'NES_codes'))
addpath_recurse(fullfile(codedir, 'gut_matlab', 'master_pipeline', 'nes_functions'))
addpath_recurse(fullfile(codedir, 'gptoolbox'))

% VIP8
dataDir = '/mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/';
pathlineDir = 'msls_output/gridCoords_nU0100_nV0100/piv/pathlines/t0_0123/' ;
dxdyFilteredFn = fullfile(pathlineDir, 'DxDyStrainFiltered/dxdy_strain_%06d.mat') ;
gridCoordDir = fullfile(dataDir, 'msls_output','gridCoords_nU0100_nV0100') ;
pathlineDir = fullfile(gridCoordDir, 'piv','pathlines','t0_0123') ;

% Laptop
% dataDir = '/Users/npmitchell/Desktop/gut/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/' ;
% dxdyFilteredFn = fullfile(pathlineDir, 'DxDyStrainFiltered/dxdy_strain_%06d.mat') ;
ameshFn = fullfile(gridCoordDir, 'sphi_cutMesh_010step','smoothed_rs_closed', 'spcMSmRSC_%06d.mat');

mkdir(fullfile(gridCoordDir, 'sphi_cutMesh_010step', 'smoothed_rs_closed_endcaps')) 
meshFn = fullfile(gridCoordDir, 'sphi_cutMesh_010step', 'smoothed_rs_closed_endcaps', 'spcMSmRSCE_%06d.mat') ;
%%%%%%%%%%%%%%%%%%%%%%%%
%% Make some meshes
% aligned mesh APDV RSC
timePoints = 10:263 ;
for ii = 1:length(timePoints)
    tp = timePoints(ii) ;
    ameshfn_ii = sprintf(ameshFn, tp) ;
    meshfn = sprintf(meshFn, tp) ;
    amesh = load(ameshfn_ii) ;
    closedMesh = amesh.spcutMeshSmRSC ;
    closedMesh.nU = nU ;
    closedMesh.nV = nV ;
    [closedMesh.f, closedMesh.v] = closeMeshEndcaps(closedMesh) ;
    save(meshfn, 'closedMesh')
end
%% end mesh making
%%%%%%%%%%%%%%%%%%%%%%%%

debugFile = './simulationTestFile.mat' ;
tmp = load(debugFile, 'cutM', 'mesh', 'refMesh', 'xyzlim', 'bwr', 'timePoints') ;

% For debugging
ameshfn_ii = sprintf(ameshFn, QS.t0) ;
amesh = load(ameshfn_ii) ;
cutM = amesh.spcutMeshSmRSC ;
cutM.nU = 100 ;
cutM.nV = 100 ;
[cutM.f, cutM.v] = closeMeshEndcaps(cutM) ;
VV = cutM.v ;
FF = cutM.f ;

% cutM = tmp.cutM ;
mesh = tmp.mesh ;
refMesh = tmp.refMesh.refMesh ;
xyzlim = tmp.xyzlim ;
bwr = tmp.bwr ;
timePoints = tmp.timePoints  ;
% save(debugFile, 'cutM', 'mesh', 'refMesh', 'xyzlim', 'bwr', 'timePoints')


%% Default options
strainstyle = 'meshes' ;  % 'meshes' 'total' 'axial' 'hoop' 'hoopCompression' 'ring' 'line'
% 'meshes' : infer metric from a series of meshes and use those as rest
%       lengths
% 'total' : all measured compress/expansion applied to bonds
% 'axial' :  compression/expansion along ap direction on bonds, so
%       epilon_zz applied to bonds with cos(beta)
% 'hoop' : compression/expansion along dv direction for all bonds, so
%       epsilon_phiphi applied to bonds with sin(beta)
% 'hoopCompression' : compression along dv direction
% 'ring' : compression/expansion along dv direction on bonds with 
%       beta==pi/2
% 'line' : compression/expansion along ap direction on bonds with beta==0

targetThetaStyle = 'quasistatic' ; % either 'plate' or 'quasistatic'
Alpha = 0 ;
Beta = 0 ;
fixVolume = false ;
fixBoundary = false ;
fixCap = false ;
restrictGrowth = false ;
poisson_ratio = 0.5 ;
thickness = 1 ;
eL_allow = 0.1 ;
maxIter = 10 ;
nU = QS.nU ;
nV = QS.nV ;
% t0Pathlines = QS.t0set() ; % Debug
t0Pathlines = QS.t0 ;
preview = false ;
Dt = 2 ;
nsteps_per_timepoint = 10 ;
niter_per_Dt = 1 ;
Ntotal = 155 ;
plot_faces = true ;
plot_dfaces = true ;
plot_wire = false ;
plot_divcurl = true ;
% timePoints = QS.xp.fileMeta.timePoints ; % Debug


% Figure Options
figWidth = 16 ; 
figHeight = 10 ;
cmin = 0.95 ;
cmax = 1.05 ;
climit_div = 0.1 ;

%% Unpack options
% Method Options
if isfield(options, 'strainstyle')
    strainstyle = options.strainstyle ;
end
if isfield(options, 'targetThetaStyle')
    targetThetaStyle = options.targetThetaStyle ;
end
if isfield(options, 'restrictGrowth')
    restrictGrowth = options.restrictGrowth ;
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

% Parameters
if isfield(options, 'thickness')
    thickness = options.thickness ;
end
if isfield(options, 'poisson_ratio')
    poisson_ratio = options.poisson_ratio ;
end
if isfield(options, 'Dt')
    Dt = options.Dt ;
end
if isfield(options, 'maxIter')
    maxIter = options.maxIter ;
end

%% Define paths % DEBUG
% pathlineDir = sprintf(QS.dir.pathlines.data, t0Pathlines) ;  % Debug
% dxdyFilteredDir = sprintf(QS.fileBase.pathlines.dxdyFiltered, t0Pathlines) ;  % Debug

%% Path options

outRoot = fullfile(pathlineDir, 'simulation') ; % Debug
exten = sprintf('_Dt%03d_nu%0.2f_t%0.2f_%03dx%03d', ...
    Dt, poisson_ratio, thickness, nU, nV) ;
exten = strrep(exten, '.', 'p') ;
if restrictGrowth
    exten = ['_rGrowth' exten ] ;
end
if fixBoundary
    exten = ['_fixB' exten ] ;
end
if fixVolume
    exten = ['_fixV' exten ] ;
end
if fixCap
    exten = ['_fixC' exten ] ;
end
exten = [ '_' targetThetaStyle exten] ;
dirname = [ strainstyle exten ] ;
outdir = fullfile(outRoot, dirname) ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

%% 
% QS.setTime(t0Pathlines) % Debug
% QS.getCurrentSPCutMeshSmRS() % Debug
% cutM = QS.currentMesh.spcutMeshSmRS ; % Debug
% QS.getCurrentSPCutMeshSmRSC() % Debug
% mesh = QS.currentMesh.spcutMeshSmRSC ; % Debug


%% Close the endcaps with a single vertex (may be a problem?)
% vtx = reshape(mesh.v, [nU, nV-1, 3]) ;
% nvtx = size(mesh.v, 1) ;
% endpt1 = mean(squeeze(vtx(1, :, :)), 1) ;
% endpt2 = mean(squeeze(vtx(end, :, :)), 1) ;
% end1 = 1:nU:nU*(nV-1) ;
% endf1 = [];
% for qq = 1:length(end1)
%     if qq + 1 <= numel(end1) 
%         endf1 = [endf1; [end1(qq), end1(qq+1), nvtx + 1]] ;
%     else
%         endf1 = [endf1; [end1(qq), end1(mod(qq+1, numel(end1))), nvtx + 1]] ;
%     end
% end
% end2 = nU:nU:nU*(nV-1) ;
% endf2 = [];
% for qq = 1:length(end2)
%     if qq + 1 <= numel(end2) 
%         endf2 = [endf2; [end2(qq), nvtx+2, end2(qq+1)]] ;
%     else
%         endf2 = [endf2; [end2(qq), nvtx+2, end2(mod(qq+1, numel(end2)))]] ;
%     end
% end
% 
% % add endcap points to vtx and triangulation
% VV = [mesh.v; endpt1; endpt2] ;
% FF = [mesh.f; endf1; endf2] ;

mesh.nU = nU ;
mesh.nV = nV ;
[FF, VV] = closeMeshEndcaps(mesh) ;
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
[eL, tarTheta] = calculateEdgeLengthsAndAngles(FF, VV);
eL0 = eL ; 

% Get bond centers in 3d
midx = 0.5 * (VV(eIDx(:, 1), 1) + VV(eIDx(:, 2), 1)) ;
midy = 0.5 * (VV(eIDx(:, 1), 2) + VV(eIDx(:, 2), 2)) ;
midz = 0.5 * (VV(eIDx(:, 1), 3) + VV(eIDx(:, 2), 3)) ;
bc3d = [midx, midy, midz] ;

%% Compute the bond orientation angles in 2d
% refMesh = load(sprintf(QS.fileName.pathlines.refMesh, t0Pathlines)) ; % Debug
% refMesh = refMesh.refMesh ; % Debug
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
if ~exist(outfn, 'file') && preview
    disp('Generating wire_beta figure...')
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
% [~,~,~,xyzlim] = QS.getXYZLims() ; % Debug

%% save simulation parameters
eIDx2dGlued = e2dg ;
save(fullfile(outdir, 'simulation_parameters.mat'), ...
    'targetThetaStyle', 'strainstyle', ...
    'fixBoundary', 'fixVolume', 'fixCap', 'V0', ...
    'thickness', 'poisson_ratio', 'capID', 'xyzlim', 'cutM', ...
    'eIDx2d', 'eIDx2dGlued', 'i3d', 'i2d', 'ghostBonds2d', 'Dt', ...
    'niter_per_Dt', 'Alpha', 'Beta', 'eL_allow') ;

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
    aux_nes_beltrami(QS, FF, VV, refMesh, capID, nU, nV, 0, bfn, bifn, xyzlim)
end

%% 
if contains(strainstyle, 'restricted')
    % Calculate the (Initial) Growth Restriction Vector Field ================
    close all; clc;

    [V2F, ~] = meshAveragingOperators(F, V);

    phi = atan2(V(:,2), V(:,1));
    theta = acos(V(:,3));

    growthV0 = V2F * [ -sin(theta) .* sin(phi), ...
        sin(theta) .* cos(phi), zeros(size(phi)) ];

    % Ensure that the restriction field is tangent to the surface
    fN = faceNormal(sphereTri);
    growthV0 = growthV0 - ...
        repmat(dot(growthV0, fN, 2), 1, 3) .* fN;

    growthV0 = normalizerow(growthV0);

    % View results ------------------------------------------------------------
    % COM = barycenter(V, F);
    % ssf = 5;
    % 
    % trisurf(sphereTri, 'FaceColor', [0.8 0.8 0.8]);
    % hold on
    % quiver3( COM(1:ssf:end,1), COM(1:ssf:end,2), COM(1:ssf:end,3), ...
    %     growthV0(1:ssf:end,1), growthV0(1:ssf:end,2), ...
    %     growthV0(1:ssf:end,3), 1, 'LineWidth', 2, 'Color', 'm' );
    %     
    % hold off
    % axis equal

    clear phi theta COM ssf V2F fN

    %% Calculate the Projected Edge Lengths Along the Growth Vector ===========
    % This will serve as the maximum projected length along the azimuthal
    % direction along which mesh edges will not be allowed to expand

    edgeVecX = V(E(:,2), 1) - V(E(:,1), 1);
    edgeVecY = V(E(:,2), 2) - V(E(:,1), 2);
    edgeVecZ = V(E(:,2), 3) - V(E(:,1), 3);

    edgeVecF = cat(3, edgeVecX(feIDx), edgeVecY(feIDx), ...
        edgeVecZ(feIDx) );

    projL0 = repmat(permute(growthV0, [1 3 2]), 1, 3, 1);
    projL0 = abs(dot(projL0, edgeVecF, 3));

    %--------------------------------------------------------------------------
    % NOTE: SETTING THE MAXIMUM EDGE LENGTH EQUAL TO THE INITIAL EDGE LENGTH
    % WILL RESULT IN ABNORMAL TERMINATION OF THE MINIMIZATION. To account for
    % this we set the maximum edge lengths to be slightly larger than the
    % initial edge lengths. We also include a small, but finite threshold for
    % edges that initially lie perpendicular to the growth direction
    %--------------------------------------------------------------------------
    projL0 = 1.1 .* projL0;                             % bonds are allowed to grow a little bit (10%)
    projL0( projL0 < 1e-10 ) = mean(projL0(:)) / 1e2;   % perpendicular bonds do not have zero projection

end



%% Consider each timestep, which averages Dt timepoints of experiment
for ii = 1:Ntotal

    tp = QS.t0 + (ii-1)*Dt ;
    
    assert(tp < max(timePoints) + 1)

    disp(['Considering tp = ' num2str(tp)])
    % Calculate Target Edge Lengths -------------------------------------------
    % Determine how bonds are strained
    switch strainstyle
        case 'all'
            titlestr = ['simulation with measured $\epsilon$'];
        case 'hoop'
            titlestr = ['simulation with measured $\epsilon_{\phi\phi}$'];  
        case 'hoopCompression'
            titlestr = ['compression with measured $\epsilon_{\phi\phi}$'];  
        case 'axial'
            titlestr = ['simulation with measured $\epsilon_{\zeta\zeta}$'];  
        case 'ring'
            titlestr = ['simulation with measured ring strain $\epsilon_{\phi\phi}$'];  
        case 'line'
            titlestr = ['simulation with measured axial line strain of $\epsilon_{\zeta\zeta}$'];  
    end            
    
    %% Determine hoop strain from experimental dx and/or dy
    if strcmpi(strainstyle, 'meshes')
        % aligned mesh APDV RSC
        ameshfn_ii = sprintf(ameshFn, tp) ;
        amesh = load(ameshfn_ii) ;
        tarMesh = amesh.spcutMeshSmRSC ;
        tarMesh.nU = nU ;
        tarMesh.nV = nV ;
        [tarMesh.f, tarMesh.v] = closeMeshEndcaps(tarMesh) ;
        % Calculate target geometry for current time point ----------------
        [eL1, tarTheta] = calculateEdgeLengthsAndAngles(tarMesh.f, tarMesh.v);
        eL1(length(eL1)+1:end) = eL(length(eL1)+1:end) ;
    else
        if strcmpi(strainstyle, 'total')
            error('here')
        else
            for qq = 1:Dt
                disp(['averaging dxs with tp=', num2str(tp + qq -1)])
                tmp = load(sprintf(dxdyFilteredFn, tp + qq - 1)) ;
                if qq == 1
                    dxs = tmp.dxs ;
                    dys = tmp.dys ;
                else
                    dxs = dxs + tmp.dxs ;
                    dys = dys + tmp.dys ;
                end
            end

            % hack for now --> go back and repeat with reflected boundaries
            dxs(:, nV) = dxs(:, 1) ; 
            dys(:, nV) = dys(:, 1) ; 
            % dxs = dxs / Dt ;
            % dys = dys / Dt ;
        end


        %% Roll off the deformation to zero at anterior and posterior ends
        % make something like a tukey window and apply to deformation
        twinAnterior = tukeywin(nU, 0.15) ;
        twinPosterior = tukeywin(nU, 0.25) ;
        twin = twinAnterior ;
        twin(round(nU*0.5):end) = twinPosterior(round(nU*0.5:end)) ;
        dxs = twin .* dxs ;
        dys = twin .* dys ;

        % Smooth dx field with modeFilter
        opt = struct('nModesY', 4) ;
        dxs = modeFilterQuasi1D(dxs, opt) ;
        dys = modeFilterQuasi1D(dys, opt) ;

        % prep for interpolation
        rux = refMesh.u(:, 1) / max(refMesh.u(:, 1)) ;
        ruy = refMesh.u(:, 2) / max(refMesh.u(:, 2)) ;
        rux = reshape(rux, [nU, nV]) ;
        ruy = reshape(ruy, [nU, nV]) ;
        if preview
            tmplim = max(max(dys(:)), max(dxs(:))) ;
            close all
            subplot(1, 2, 1)
            imagesc(linspace(0,1,nU), linspace(0,1,nV), dxs')
            xlabel('$\zeta$','interpreter', 'latex')
            ylabel('$\phi$','interpreter', 'latex')
            title('$\delta \zeta$', 'interpreter', 'latex')
            caxis([-tmplim, tmplim])
            colormap bwr
            colorbar('location', 'southoutside')
            subplot(1, 2, 2)
            imagesc(linspace(0,1,nU), linspace(0,1,nV), dys')
            xlabel('$\zeta$','interpreter', 'latex')
            ylabel('$\phi$','interpreter', 'latex')
            title('$\delta \phi$', 'interpreter', 'latex')
            caxis([-tmplim, tmplim])
            colormap bwr
            colorbar('location', 'southoutside')
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
        if any(contains(lower(strainstyle), {'all', 'hoop', 'ring'}))
            % Interpolation
            dyinterp = griddedInterpolant(rux, ruy, dys, 'linear') ;
            dyi = dyinterp(bxy2dg(:, 1), bxy2dg(:, 2)) ;
            dymag = zeros(size(bbeta)) ;
            dymag(i3d) = dyi ;

            %% Preview 3d assignment
            if preview
                clf
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
            case 'hoopCompression'
                dymag(dymag > 0) = 0 ;
                scalefactor = 1 + dymag .* abs(sin(bbeta)) ;
            case 'axial'
                scalefactor = 1 + dxmag .* abs(cos(bbeta)) ;
            case 'ring'
                scalefactor = 1 + dymag .* abs(sin(bbeta)) ;
                scalefactor(abs(sin(bbeta))<0.99) = 1 ;
            case 'line'
                scalefactor = 1 + dxmag .* abs(cos(bbeta)) ;
                scalefactor(abs(cos(bbeta))<0.99) = 1 ;
        end

        % Calculate target geometry for current time point --------------------
        [eL1, tarTheta] = calculateEdgeLengthsAndAngles(FF, VV);
        eL = eL1 .* scalefactor ;
    end
    
    if restrictGrowth
        % Transform the growth vector into the current time point
        curGrowthV = transformVectorField3Dto3DMesh( growthV0, V, curV, ...
            F, (1:size(F,1)).' );
        curGrowthV = normalizerow(curGrowthV);
    end
    
    %% Check for self-intersections
    [intx_exist, faces_intx ] = mesh_self_intersection_3d(FF, VV) ;
    if intx_exist
        error([num2str(length(find(faces_intx))) 'self intersections exist'])
    end
    
    %% Check for triangle inequality: 
    % sum of edge 2 lengths cannot be less than the third edgelength
    % #faces x 3 array of target edge lengths on faces
    % (i,j)th element gives target lengths of edge opposite j in face i
    feL = eL(feIDx) ;
    % circshift(arr,1,2) shifts by 1 element along dim 2
    try
        assert(all(all(feL - circshift(feL, 1, 2) - circshift(feL, 2, 2) < 0)))
    catch
        dmyk = 0 ;
        bonds_ok = false ;
        while ~bonds_ok
            % which faces violate the triangle inequality?
            fails = find(any(feL - circshift(feL, 1, 2) - circshift(feL, 2, 2) > 0, 2)) ;
            nfails = length(fails) ;
            disp(['triangle inequality is not satisfied in ' num2str(nfails) ' triangles'])

            % Find unique edges to fix
            [badtri, nodes] = find(feL - circshift(feL, 1, 2) - circshift(feL, 2, 2) > 0) ;
            edges2relax = [] ;
            for pq = 1:length(badtri)
                thisTri = badtri(pq,:) ;
                edges2relax = [edges2relax feIDx(thisTri, :)] ;
            end
            edges2relax = unique(edges2relax) ;
            
            % preview the failures
            clf
            colormap cividis
            for axid = 1:4
                subplot(4, 2, axid)
                trisurf(FF, VV(:, 1), VV(:, 2), VV(:, 3), ...
                    sum(sign(feL - circshift(feL, 1, 2) - circshift(feL, 2, 2)), 2), ...
                    'edgecolor', 'none')
                axis equal
                colormap(reds)
                view(0, (axid-1) * 90)
                ylabel('lateral position', 'interpreter', 'latex')
                zlabel('dv position', 'interpreter', 'latex')
                xlabel('ap position', 'interpreter', 'latex')
            end
            subplot(4, 2, 5)
            % Plot eL / eL0 
            histogram(eL ./ eL1)
            title('$e_L / e_L^0$', 'interpreter', 'latex')
            subplot(4, 2, 6)
            histogram(eL(edges2relax) ./ eL1(edges2relax))
            title('$e_L($select$) / e_L^0($select$)$', ...
                'interpreter', 'latex')
            
            sgtitle(['pass ' num2str(dmyk)], 'interpreter', 'latex')
            
            % Also plot all strains post-correction
            
            pause(1)
            clf
            
            % Average with current length
            eL(edges2relax) = eL_allow * eL1(edges2relax) + ...
                 (1-eL_allow) * eL(edges2relax) ;
             
            % If we have been updating this bad bond for too many
            % iterations, just give it the previous bond length
            if dmyk > 50
                eL(edges2relax) = eL1(edges2relax) ;
            end
            
            % Reassign face edge lengths based on new definition of eL
            feL = eL(feIDx) ;
            % Re-check the triangle inequality now with relaxed bonds
            bonds_ok = all(all(feL - circshift(feL, 1, 2) - circshift(feL, 2, 2) < 0)) ;
            dmyk = dmyk + 1 ;
        end
    end
    
    %% check for negative target edgelengths
    assert(all(eL ./ eL0 > 0))
    
    %% Consider gently bounding the minimum/maximum deviation?
    % error('here')

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
    if strcmpi(targetThetaStyle, 'quasistatic')
        % % The unit normal vector for each face
        % % see https://www.cs.utexas.edu/users/evouga/uploads/4/5/6/8/45689883/turning.pdf
        % fN = faceNormal( triangulation(F, tarV) );
        % 
        % N1 = fN(edgeFace(:,1), :);
        % N2 = fN(edgeFace(:,2), :);
        % 
        % crossN = cross(N2, N1, 2);
        % dotN = dot(N1, N2, 2);
        % 
        % tarTheta = 2 .* atan2( dot(crossN, eij, 2), 1 + dotN );
        % [~, tarTheta] = calculateEdgeLengthsAndAngles(F, V) ;
        disp('tarTheta already computed for this timepoint')
    elseif strcmpi(targetThetaStyle, 'plate')
        tarTheta = zeros(size(eL, 1), 1) ;
    else
        erorr(['Did not recognize targetTheta style: ' targetThetaStyle])
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
    
    if fixCap
        fixedIDx = capID(:) ;
    else
        fixedIDx = capID(1) ;
    end
    fixedX = V1(fixedIDx, :) ;
    
    %% Check magnitudes of energies
    [Eb0, Eb0_faces] = calculateBendEnergy(FF, VV, eL, tarTheta, poisson_ratio, thickness) ;
    Efp0 = calculateFixedPointEnergy(FF, VV, fixedIDx, fixedX, Alpha) ;
    Efv0 = calculateFixedVolumeEnergy(FF, VV, targetVolume, Beta) ;
    [Es0, Es0_faces] = calculateStretchEnergy(FF, VV, eL, poisson_ratio) ;

    % Compute per-face energies
    % Eb0_faces = calculateBendEnergyFaces(FF, VV, eL, tarTheta, ...
    %     poisson_ratio, thickness) ; % Debug -- turn on
    % Es0_faces = calculateStretchEnergyFaces(FF, VV, eL, poisson_ratio) ; % Debug -- turn on
    if restrictGrowth
        [Egr0, projL, isValid] = calculateGrowthRestrictionEnergy(F, V, ...
            growthVec, maxProjL, mu) ;
        Egr0_faces = calculateGrowthRestrictionEnergy(F, V, ...
            growthVec, maxProjL, mu) ;
    else
        Egr0 = 0 ;
        Egr0_faces = 0 ;
    end


    %% Run Elastic Relaxation =============================================
    fn = fullfile(outdir, 'vertices', sprintf('vertices_%03d.mat', ii)) ;
    if exist(fn, 'file')
        disp(['time step already computed, loading ' fn])
        tmp = load(fn) ;
        VV = tmp.VV ;
        FF = tmp.FF ;
        Eb = tmp.Eb ;
        Efp = tmp.Efp ;
        Efv = tmp.Efv ;
        Es = tmp.Es ;
        Eb_faces = tmp.Eb_faces ;
        Es_faces = tmp.Es_faces ;
        Eb0_faces = tmp.Eb0_faces ;
        Es0_faces = tmp.Es0_faces ;
        Eb0 = tmp.Eb0 ;
        Efp0 = tmp.Efp0 ;
        Efv0 = tmp.Efv0;
        Es0 = tmp.Es0 ;
    else
        tic
        V4min = V1 ; 
        if fixBoundary && fixVolume 
            error('check params here')
            V4min = minimizeElasticEnergy( FF, V4min, eL, ...
                'TargetAngles', tarTheta, ...
                'Thickness', thickness, ...
                'Poisson', poisson_ratio, ...
                'MaxIterations', maxIter, ...
                'iterDisplay', 1, ...
                'Alpha', Alpha, ...
                'Beta', Beta, ...
                'FixBoundary', ...
                'targetVertices', capID(:), ...
                'targetLocations', V1(capID(:), :), ...
                'FixVolume', 'TargetVolume', targetVolume); 
        elseif fixBoundary 
            error('check params here')
            V4min = minimizeElasticEnergy( FF, V4min, eL, ...
                'TargetAngles', tarTheta, ...
                'Thickness', thickness, ...
                'Poisson', poisson_ratio, ...
                'MaxIterations', maxIter, ...
                'iterDisplay', 1, ...
                'Alpha', Alpha, ...
                'Beta', Beta, ...
                'targetVertices', capID(:), ...
                'targetLocations', V1(capID(:), :), ...
                'FixBoundary');    
        elseif fixVolume
            V4min = minimizeElasticEnergy( FF, V4min, eL, ...
                'TargetAngles', tarTheta, ...
                'Thickness', thickness, ...
                'Poisson', poisson_ratio, ...
                'MaxIterations', maxIter, ...
                'Past', 500, 'Delta', 1e-7, ... % L-BFGS parameters
                'iterDisplay', 10, ...
                'Alpha', Alpha, ...             % fixed vertex coeff
                'Beta', Beta, ...               % fixed volume coeff
                'FixVolume', ...
                'TargetVolume', targetVolume, ...
                'targetVertices', capID(end-1), ...
                'targetLocations', V1(capID(end-1), :)) ;
        else
            disp('WARNING: not fixing Volume or boundary! Fixing one vertex')
            disp(['initial distance from target:' num2str(sum(vecnorm(V4min - V0, 2, 2)))])
            [eL, tarTheta] = calculateEdgeLengthsAndAngles(FF, V4min) ;
            V4min = minimizeElasticEnergy( FF, V4min, eL, ...
                'TargetAngles', tarTheta, ...
                'Thickness', thickness, ...
                'Poisson', poisson_ratio, ...
                'MaxIterations', maxIter, ...
                'Past', 500, 'Delta', 1e-7, ... % L-BFGS parameters
                'iterDisplay', 10, ...
                'Alpha', Alpha, ...             % fixed vertex coeff
                'Beta', Beta, ...               % fixed volume coeff
                'targetVertices', capID(end-1), ...
                'targetLocations', V1(capID(end-1), :)) ;
        end
        VV = V4min ;
        toc
        
        
        % For scale = 1
        % ||dx|| = 14.3817936 ||x|| = 12410.9977 f(x) =   19.0601791136
        % scale = 5
        %  ||dx|| = 776026668.5042272 ||x|| = 60160.2956 f(x) = 1314937352.3637316227
        
        %% Check
        % % Directed edge vectors
        % eijp = VV(eIDx(:,2), :) - VV(eIDx(:,1), :);
        % 
        % % Target edge lengths
        % eLp = sqrt( sum( eijp.^2, 2 ) );

        %% New energies
        [Eb, Eb_faces] = calculateBendEnergy(FF, VV, eL, tarTheta, poisson_ratio, thickness) ;
        Efp = calculateFixedPointEnergy(FF, VV, fixedIDx, fixedX, Alpha) ;
        Efv = calculateFixedVolumeEnergy(FF, VV, targetVolume, Beta) ;
        [Es, Es_faces] = calculateStretchEnergy(FF, VV, eL, poisson_ratio) ;
        
        % Compute per-face energies 
        if restrictGrowth
            [Egr, projL, isValid] = calculateGrowthRestrictionEnergy(F, V, ...
                growthVec, maxProjL, mu) ;
            Egr_faces = calculateGrowthRestrictionEnergy(F, V, ...
                growthVec, maxProjL, mu) ;
        else
            Egr = 0 ;
            Egr_faces = 0 ;
        end
        
        % Check that energies went down?
        total_energy0 = Eb0 + Efp0 + Efv0 + Es0 + Egr0 ;
        total_energy = Eb + Efp + Efv + Es + Egr ;
        
        %% Save vertices
        save(fn, 'VV', 'FF', 'Eb', 'Efp', 'Efv', 'Es', 'Egr',...
            'Eb0', 'Efp0', 'Efv0', 'Es0', 'Egr0', ...
            'Eb0_faces', 'Es0_faces', 'Egr0_faces', ...
            'Eb_faces', 'Es_faces', 'Egr_faces') ;
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
        caxis(max(abs(a0strain(:))) * [-1, 1])
        colormap(bwr)
        c = colorbar ;
        c.Color = 'w' ;
        c.Label.Interpreter = 'latex' ;
        c.Label.String = '$ A / A_0$' ;

        % Figure properties
        set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w')
        set(gcf, 'color', 'k')

        titlestr = ['face dilation, $ (A - A_0) / A_0$: $t=$' num2str(tp-QS.t0)]; 
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
        c.Label.String = '$\delta A/ A_{prev}$' ;

        % Figure properties
        set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w')
        set(gcf, 'color', 'k')

        titlestr = ['face dilation, $\delta A / A_{prev}$: $t=$' num2str(tp-QS.t0)]; 
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
        clims = [] ;
        titlestr = ['$t=$' num2str(tp * QS.timeInterval) ' ' QS.timeUnits] ;
        aux_nes_wire(QS, eIDx, VV, eL, eL0, clims, wirefn, titlestr)
    end
        
    %% Plot flow invariants
    divcurlfn = fullfile(outdir, 'divcurl', sprintf('divcurl_%03d.png', ii)) ;
    if ~exist(divcurlfn, 'file') && plot_divcurl
        aux_nes_divcurl(QS, FF, VV, V1, climit_div, divcurlfn, ii*Dt, xyzlim)
    end
    
    %% Compute Beltrami
    bfn = fullfile(outdir, 'beltrami', sprintf('beltrami_%03d.mat', ii)) ;
    bifn = fullfile(outdir, 'beltrami_images', ...
        sprintf('beltrami_%03d.png', ii)) ;
    if ~exist(bfn, 'file') || ~exist(bifn, 'file')
        aux_nes_beltrami(QS, FF, VV, refMesh, capID, nU, nV, ii*Dt, bfn, bifn, xyzlim)
    end
    
    %% Store energies as we go
    EbV(ii) = Eb ;
    EfpV(ii) = Efp ;
    EfvV(ii) = Efv ;
    EsV(ii) = Es ;
    Eb0V(ii) = Eb0 ;
    Efp0V(ii) = Efp0 ;
    Efv0V(ii) = Efv0;
    Es0V(ii) = Es0 ;
    
    % Plot the energies so far
    clf
    plot(EbV, '.-') ;
    hold on;
    plot(EfpV, '.-') ;
    plot(EfvV, '.-') ;
    plot(EsV, '.-') ;
    legend({'bending', 'point', 'volumetric', 'stretching'})
    pause(0.1)
    %waitfor(gcf)
    
    
    % Prepare for next timepoint
    V1 = VV ;
end

