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
outdir = './output/' ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

%% Default options
Alpha = 0 ;
Beta = 0 ;
fixVolume = false ;
fixBoundary = false ;
fixCap = false ;
restrictGrowth = false ;
poisson_ratio = 0.5 ;
thickness = 1 ;
maxIter = 10 ;
preview = false ;
Dt = 2 ;

%% Construct Physical/Target Configurations and Geometries ================
t0 = 123 ;
timePoints = t0:263 ;
meshFn = './meshes/spcMSmRSCE_%06d.mat' ;
meshfn_ii = sprintf(meshFn, t0) ;
mesh = load(meshfn_ii) ;
FF = mesh.closedMesh.f ;
VV = mesh.closedMesh.v ;
nU = mesh.closedMesh.nU ;
nV = mesh.closedMesh.nV ;
% Target edge lengths
[eLt, tarTheta] = calculateEdgeLengthsAndAngles(FF, VV);
eL0 = eLt ; 

%% Initial face areas
a0 = 0.5 * doublearea(VV, FF) ;
V0 = VV ;  % reference vertices
V1 = VV ;  % previous timepoint vertices (which are reference at the moment)

%% Consider each timestep, which averages Dt timepoints of experiment
for ii = 1:floor(length(timePoints)/Dt)
    tp = t0 + (ii-1)*Dt ;    
    assert(tp < max(timePoints) + 1)
    
    %% Load target mesh from disk
    % aligned mesh APDV RSC
    meshfn_ii = sprintf(meshFn, tp) ;
    tarMesh = load(meshfn_ii) ;
    tarMesh = tarMesh.closedMesh ;
    Ft = tarMesh.f ;
    Vt = tarMesh.v ;
    % Calculate target geometry for current time point ----------------
    [eLt, tarTheta] = calculateEdgeLengthsAndAngles(tarMesh.f, tarMesh.v);
    
    %% Check for self-intersections
    [intx_exist, faces_intx ] = mesh_self_intersection_3d(Ft, Vt) ;
    if intx_exist
        error([num2str(length(find(faces_intx))) 'self intersections exist'])
    end
    
    %% Check for triangle inequality in target Mesh
    % sum of edge 2 lengths cannot be less than the third edgelength
    % #faces x 3 array of target edge lengths on faces
    % Construct Topolgical Structure Tools ===================================
    tri = triangulation(Ft, Vt) ;
    [eIDx, feIDx, bulkEdgeIDx] = topologicalStructureTools(tri) ;    
    % (i,j)th element gives target lengths of edge opposite j in face i
    feLt = eLt(feIDx) ;
    % circshift(arr,1,2) shifts by 1 element along dim 2
    assert(all(all(feLt - circshift(feLt, 1, 2) - circshift(feLt, 2, 2) < 0)))
   
    %% check for negative target edgelengths
    assert(all(eLt ./ eL0 > 0))
    
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
    
    %% Fixed points
    capID = [nU*(nV-1)+1, nU*(nV-1)+2]' ;
    if fixCap
        fixedIDx = capID(:) ;
    else
        fixedIDx = capID(1) ;
    end
    % Can choose to keep vertices fixed or move them to target
    fixedX = Vt(fixedIDx, :) ;
    
    %% Check magnitudes of energies
    [Eb0, Eb0_faces] = calculateBendEnergy(FF, VV, eLt, tarTheta, poisson_ratio, thickness) ;
    Efp0 = calculateFixedPointEnergy(FF, VV, fixedIDx, fixedX, Alpha) ;
    Efv0 = calculateFixedVolumeEnergy(FF, VV, targetVolume, Beta) ;
    [Es0, Es0_faces] = calculateStretchEnergy(FF, VV, eLt, poisson_ratio) ;
    % Compute per-face energies 
    if restrictGrowth
        [Egr0, projL, isValid] = calculateGrowthRestrictionEnergy(FF, VV, ...
            growthVec, maxProjL, mu) ;
        Egr0_faces = calculateGrowthRestrictionEnergy(FF, VV, ...
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
        V1 = tmp.VV ;
        VV = V1 ;
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
        disp('WARNING: not fixing Volume or boundary! Fixing one vertex')
        disp(['initial distance from target:' num2str(sum(vecnorm(VV - V0, 2, 2)))])
        
        [eL_cur, curTheta] = calculateEdgeLengthsAndAngles(FF, VV) ;
        assert(sum(eLt - eL_cur) == 0)
        assert(sum(vecnorm(VV-Vt, 2, 2)) == 0)
        assert(sum(tarTheta - curTheta) == 0)
        
        V1 = minimizeElasticEnergy( FF, VV, eLt, ...
            'TargetAngles', tarTheta, ...
            'Thickness', thickness, ...
            'Poisson', poisson_ratio, ...
            'MaxIterations', maxIter, ...
            'Past', 500, 'Delta', 1e-7, ... % L-BFGS parameters
            'iterDisplay', 10, ...
            'Alpha', Alpha, ...             % fixed vertex coeff
            'Beta', Beta, ...               % fixed volume coeff
            'targetVertices', capID(end-1), ...
            'targetLocations', Vt(capID(end-1), :)) ;

        VV = V1 ;
        toc
        
        %% New energies
        [Eb, Eb_faces] = calculateBendEnergy(FF, VV, eLt, tarTheta, poisson_ratio, thickness) ;
        Efp = calculateFixedPointEnergy(FF, VV, fixedIDx, fixedX, Alpha) ;
        Efv = calculateFixedVolumeEnergy(FF, VV, targetVolume, Beta) ;
        [Es, Es_faces] = calculateStretchEnergy(FF, VV, eLt, poisson_ratio) ;
        
        % Compute per-face energies 
        if restrictGrowth
            [Egr, projL, isValid] = calculateGrowthRestrictionEnergy(FF, VV, ...
                growthVec, maxProjL, mu) ;
            Egr_faces = calculateGrowthRestrictionEnergy(FF, VV, ...
                growthVec, maxProjL, mu) ;
        else
            Egr = 0 ;
            Egr_faces = 0 ;
        end
        
        % Check that energies went down?
        total_energy0 = Eb0 + Efp0 + Efv0 + Es0 + Egr0 ;
        total_energy = Eb + Efp + Efv + Es + Egr ;
        disp(['old total_energy: ' num2str(total_energy0)])
        disp(['new total_energy: ' num2str(total_energy)])
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
    % waitfor(gcf)
    
end

