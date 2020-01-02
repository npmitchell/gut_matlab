%% Load all spcutMesh objects and add a field (ringpath_ss) to each

for t = xp.fileMeta.timePoints
    % Load the spcutMesh for this timepoint
    disp('Loading spcutMesh from disk...')
    load(sprintf(spcutMeshBase, t), 'spcutMesh') ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute ringpath_ss, the mean distance traveled from one
    % line of constant u to the next
    disp('Computing ringpath_ss...')
    uv = spcutMesh.uv ;
    new3d = spcutMesh.v ;
    nU = size(spcutMesh.radii_from_mean, 1) ;
    nV = size(spcutMesh.radii_from_mean, 2) ;

    % The distance from one hoop to another is the
    % difference in position from (u_i, v_i) to (u_{i+1}, v_i).
    dsuphi = reshape([vecnorm(diff(new3d), 2, 2); 0], [nU, nV]) ;
    ringpath_ds = nanmean(dsuphi(1:(end-1), :), 2) ;
    ringpath_ss = cumsum([0; ringpath_ds]) ;
    spcutMesh.ringpath_ss = ringpath_ss ;

    % Resave s,phi and their 3D embedding
    save(sprintf(spcutMeshBase, t), 'spcutMesh') ;
    
    plot(ringpath_ss)
    pause(0.1)
end

%% radii_from_mean -> radii_from_avgpts and add centerline to spcutMesh
% Load all spcutMesh objects and add compute/save radii_from_avgpts.
for t = xp.fileMeta.timePoints
    % Load the spcutMesh for this timepoint
    disp(['Loading spcutMesh from disk... [t = ' num2str(t) ']'])
    load(sprintf(spcutMeshBase, t), 'spcutMesh') ;
    
    % Load new centerline
    fn = sprintf(clineDVhoopBase, t) ;
    disp(['Loading new centerline from ' fn])
    load(fn, 'avgpts')

    % Add nU and nV to spcutMesh
    try
        nU = size(spcutMesh.radii_from_mean, 1) ;
        nV = size(spcutMesh.radii_from_mean, 2) ;
        % Remove radii_from_mean
        spcutMesh.radii_from_mean_uniform = spcutMesh.radii_from_mean ;
        spcutMesh = rmfield(spcutMesh, 'radii_from_mean') ;
    catch
        nU = size(spcutMesh.radii_from_mean_uniform, 1) ;
        nV = size(spcutMesh.radii_from_mean_uniform, 2) ;
    end
    spcutMesh.nU = nU ;
    spcutMesh.nV = nV ;

    % Recompute radii_from_avgpts 
    sphi3d = reshape(spcutMesh.v, [nU, nV, 3]) ;
    radii_from_avgpts = zeros(size(sphi3d, 1), size(sphi3d, 2)) ;
    for jj = 1:nU
        % Consider this hoop
        hoop = squeeze(sphi3d(jj, :, :)) ;
        radii_from_avgpts(jj, :) = vecnorm(hoop - avgpts(jj, :), 2, 2) ;
    end
    spcutMesh.radii_from_avgpts = radii_from_avgpts ;
    
    
    % Resave s,phi and their 3D embedding
    save(sprintf(spcutMeshBase, t), 'spcutMesh') ;
    clearvars spcutMesh radii_from_avgpts radii_from_mean sphi3d hoop
    clearvars avgpts
end


%% mss, mcline, avgpts
% Load all spcutMesh objects and add mss, mcline, avgpts.
for t = xp.fileMeta.timePoints
    % Load the spcutMesh for this timepoint
    disp(['Loading spcutMesh from disk... [t = ' num2str(t) ']'])
    load(sprintf(spcutMeshBase, t), 'spcutMesh') ;
    
    % Load new centerline
    fn = sprintf(clineDVhoopBase, t) ;
    disp(['Loading new centerline from ' fn])
    load(fn, 'avgpts', 'mss', 'mcline')
    
    spcutMesh.mss = mss ;       % from uniform sampling, also stored in centerline
    spcutMesh.mcline = mcline ; % from uniform sampling, also stored in centerline
    spcutMesh.avgpts = avgpts ; % from uniform sampling, also stored in centerline
    
    % Resave s,phi and their 3D embedding
    save(sprintf(spcutMeshBase, t), 'spcutMesh') ;
    clearvars spcutMesh avgpts mcline mss
end

%% avgpts_ss 
% Load all spcutMesh objects and add mss, mcline, avgpts.
for t = xp.fileMeta.timePoints
    % Load the spcutMesh for this timepoint
    disp(['Loading spcutMesh from disk... [t = ' num2str(t) ']'])
    load(sprintf(spcutMeshBase, t), 'spcutMesh') ;
    
    spcutMesh.avgpts_ss = ss_from_xyz(spcutMesh.avgpts) ; 
    
    % Resave s,phi and their 3D embedding
    save(sprintf(spcutMeshBase, t), 'spcutMesh') ;
    clearvars spcutMesh avgpts 
end

%% ringpath_ss 
% Load all spcutMesh objects and add mss, mcline, avgpts.
for t = xp.fileMeta.timePoints
    % Load the spcutMesh for this timepoint
    disp(['Loading spcutMesh from disk... [t = ' num2str(t) ']'])
    load(sprintf(spcutMeshBase, t), 'spcutMesh') ;
    
    spcutMesh.ringpath_ss = spcutMesh.ringpath_ss * resolution ; 
    
    % Resave s,phi and their 3D embedding
    save(sprintf(spcutMeshBase, t), 'spcutMesh') ;
    clearvars spcutMesh avgpts 
end
       

%% spcutMesh.vn
% Load all spcutMesh objects and add mss, mcline, avgpts.
for t = xp.fileMeta.timePoints
    % Load the spcutMesh for this timepoint
    disp(['Loading spcutMesh from disk... [t = ' num2str(t) ']'])
    load(sprintf(spcutMeshBase, t), 'spcutMesh') ;
    uv = spcutMesh.uv ;
    
    % Load the cutMesh for this timepoint
    disp(['Loading cutMesh from disk... [t = ' num2str(t) ']'])
    cutMeshfn = fullfile(cutFolder, [fileNameBase, '_cutMesh.mat']) ;
    cutMeshfn = sprintf(cutMeshfn, t) ;
    load(cutMeshfn) 
    
    % sample the normal vectors in the UV space of cutMesh, apply to the
    % points that live in the same UV space stored in spcutMesh.
    tileCount = [1 1];  % how many above, how many below
    [ TF, TV2D, TV3D, TVN3D ] = tileAnnularCutMesh( cutMesh, tileCount );
    Fnx = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TVN3D(:, 1)) ;
    Fny = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TVN3D(:, 2)) ;
    Fnz = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TVN3D(:, 3)) ;
    spcutMesh.vn = zeros(size(spcutMesh.v)) ;
    spcutMesh.vn(:, 1) = Fnx(uv(:, 1), uv(:, 2)) ;
    spcutMesh.vn(:, 2) = Fny(uv(:, 1), uv(:, 2)) ;
    spcutMesh.vn(:, 3) = Fnz(uv(:, 1), uv(:, 2)) ; 
    
    if isfield(spcutMesh, 'u')
        spcutMesh = rmfield(spcutMesh, 'u') ;
    end
    
    % Resave s,phi and their 3D embedding
    save(sprintf(spcutMeshBase, t), 'spcutMesh') ;
    clearvars spcutMesh avgpts 
end