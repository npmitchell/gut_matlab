%% Compute metric strain for evolving mesh ===============================
% Compute the metric strain for topological equivalent meshes over time
% (same connectivity, different geometry in embedding coordinates)

clear; close all; clc;

% Figure Options
figWidth = 20 ; 
figHeight = 12 ;
cmin = 0.95 ;
cmax = 1.05 ;

% Add paths
gutDir = '/mnt/data/code/gut_matlab/' ;
addpath(fullfile(gutDir, 'addpath_recurse')) ;
addpath_recurse(fullfile(gutDir, 'geometry')) ;
addpath_recurse(fullfile(gutDir, 'mesh_handling')) ;
addpath_recurse(fullfile(gutDir, 'plotting')) ;

% Path options
NESpath = '/mnt/data/code/NonEuclideanShells/NES/' ;
addpath_recurse(NESpath)

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

%% Build from PLY
mslsDir = '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/' ;
mslsDir = [mslsDir 'Time6views_60sec_1.4um_25x_obis1.5_2/data/' ] ;
mslsDir = [mslsDir 'deconvolved_16bit_smallbox/'] ;
mslsDir = [mslsDir 'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1_20200129/'] ;
rscDir = fullfile(mslsDir, 'sphi_cutMesh_010step_nU0100_nV0100/smoothed_rs_closed/') ;
rscPlyFn = fullfile(rscDir, '%06d_spcMSmRSC.ply') ;
rscMatFn = fullfile(rscDir, '%06d_spcMSmRSC.mat') ;
for tp = 110
    % Load timepoint's mesh
    % build PLY from matfile
    tmp = load(sprintf(rscMatFn, tp)) ;
    mesh = tmp.spcutMeshSmRSC ;
    clearvars tmp
    
    % make PLY out of this closed mesh
    % if ~exist(sprintf(rscPlyFn, tp)) ;
    %     disp('Could not load ply, building it here')
    %     % check it
    %     % trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3))
    %     plywrite_with_normals(sprintf(rscPlyFn, tp), ...
    %         mesh.f, mesh.v, mesh.vn)
    % end

    % build PLY from matfile
    tmp = load(sprintf(rscMatFn, tp + 1)) ;
    mesh2 = tmp.spcutMeshSmRSC ;
    clearvars tmp
    
    % Store mesh vertices and faces
    VV = mesh.v ;
    FF = mesh.f ;
    UU = mesh.u ;
    V2 = mesh2.v ;
    F2 = mesh2.f ;
    U2 = mesh2.u ;
    
    % Here manually define nU, nV --> in future, load from QS
    nU = 100 ;
    nV = 100 ;

    % Construct Triangulation
    tri = triangulation(FF, VV) ;
    tri2 = triangulation(FF, VV) ;

    % Plotting options
    cmap = bwr ;

    %% Construct Topolgical Structure Tools ===============================
    % eg = metricStrainSPhiGridMesh(mesh, mesh2) ;
    
    % Metric for mesh
    g0cell = inducedMetric(FF, VV, UU) ;
    % Convert metric tensor to #Faces x 2 x 2
    g0 = zeros(size(g0cell, 1), 2, 2) ;
    for qq = 1:length(size(g0cell, 1)) 
        g0(qq, :, :) = g0cell{qq} ;
    end
    
    % Metric for mesh2
    g1cell = inducedMetric(F2, V2, U2) ;
    % Convert metric tensor to #Faces x 3
    g1 = zeros(size(g1cell, 1), 2, 2) ;
    for qq = 1:length(size(g1cell, 1)) 
        g1(qq, :, :) = g1cell{qq} ;
    end
    
    % Metric strain
    eg = 0.5 * (g1 - g0) ;
    
    % Build deformation tensor over mesh, whose elements are the change in
    % squared length in each direction (u.u,u.v,v.v) 
    dg = 0 * eg ;
    % Grab which bonds are along u and which along v
    mesh.nU = nU ;
    [labels, dbonds, tST] = labelRectilinearMeshBonds(mesh) ;
    du = dbonds.baseSpace.u2D ;
    dv = dbonds.baseSpace.v2D ;
    for qq = 1:length(size(g1cell, 1))
        gq = squeeze(eg(qq, :, :)) ;
        tg(qq, 1, 1) = du(qq, :) * gq * du(qq, :) ;
        tg(qq, 1, 2) = du(qq, :) * gq * dv(qq, :) ;
        tg(qq, 2, 1) = dv(qq, :) * gq * du(qq, :) ;
        tg(qq, 2, 2) = dv(qq, :) * gq * du(qq, :) ;
    end
    
    %% Plot the metric components on trisurf
    denom = sqrt(tg(:, 1) .* tg(:, 3)) ;
    colors = tg(:, 2) ./ denom ;
    trisurf(FF, VV(:, 1), VV(:, 2), VV(:, 3), colors)
    caxis([-1, 1])
    axis equal
    cb = colorbar() ;
    title(['Initial metric, $t=$' sprintf('%03d', tp)], ...
        'Interpreter', 'Latex')
    xlabel('AP position, [$\mu$m]', 'Interpreter', 'Latex')
    ylabel('lateral position, [$\mu$m]', 'Interpreter', 'Latex')
    zlabel('DV position, [$\mu$m]', 'Interpreter', 'Latex')
    ylabel(cb, '$g_{s\phi} / |g_{ss} g_{\phi\phi}|$', 'Interpreter', 'Latex')
    
    error('here')
    
    %% Compute the bond orientation angles
    dx = ss(eIDx(:, 2)) - ss(eIDx(:, 1)) ;
    beta = acos(abs(eij(:, 1)) ./ eL0) ;

    %% Compute new edge lengths

    %% New metric

    %% Difference in metrics

    
end
    
