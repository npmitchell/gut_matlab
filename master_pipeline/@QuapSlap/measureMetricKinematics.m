function [gdot_apM, HH_apM, divv_apM, veln_apM] = ...
    measureMetricKinematics(QS, options)
%[gdot_apM, HH_apM, divv_apM, veln_apM] = measureMetricKinematics(QS, options)
%   Measure degree of incompressibility of the flow on the evolving surface
%   Out-of-plane motion is v_n * 2H, where v_n is normal velocity and H is
%   mean curvature.
%   In-plane motion considered here is div(v_t) where v_t is tangential
%   velocity on the curved surface.
%   The difference div(v_t) - vn*2H = Tr[g^{-1} dot{g}], which is a measure
%   of isotropic metric change over time (dot is the partial derivative wrt
%   time). 

% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%   timePoints : numeric 1D array
%       the timepoints to consider for the measurement. For ex, could
%       choose subset of the QS experiment timePoints
%   alphaVal : float
%       the opacity of the heatmap to overlay
%   invertImage : bool
%       invert the data pullback under the velocity map
%
% Returns 
% -------
%
% Saves to disk
% -------------
% 
% NPMitchell 2020

%% Default options 
overwrite = false ;
plot_Hgdot = true ;
plot_flows = true ;
plot_factors = true ;
lambda = 0.01 ; 
% by default, lambda_mesh = lambda, whether defined here or in options
lambda_err = 0.01 ;
climit = 0.2 ;
climit_err = 0.2 ;
climit_veln = climit * 10 ;
climit_H = climit * 2 ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 

%% Unpack options & assign defaults
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'plot_Hgdot')
    plot_Hgdot = options.plot_Hgdot ;
end
if isfield(options, 'plot_flows')
    plot_flows = options.plot_flows ;
end
if isfield(options, 'plot_factors')
    plot_factors = options.plot_factors ;
end
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_err')
    lambda_err = options.lambda_err ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
else
    % default lambda_mesh is equal to lambda 
    lambda_mesh = lambda ;
end
if isfield(options, 'climit')
    climit = options.climit ;
end
if isfield(options, 'climit_err')
    climit_err = options.climit_err ;
end
if isfield(options, 'climit_veln')
    climit_veln = options.climit_veln ;
end
if isfield(options, 'climit_H')
    climit_H = options.climit_H ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end

%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
    sresStr = '' ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
    sresStr = 'doubleRes_' ;
else 
    error("Could not parse samplingResolution: set to '1x' or '2x'")
end

%% Unpack QS
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;
mKDir = fullfile(QS.dir.metricKinematics, ...
    strrep(sprintf([sresStr 'lambda%0.3f_lerr%0.3f_lmesh%0.3f'], ...
    lambda, lambda_err, lambda_mesh), '.', 'p'));
folds = load(QS.fileName.fold) ;
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;

%% Colormap
bwr256 = bluewhitered(256) ;

%% Load vertex-based velocity measurements
if doubleResolution
    vvsmMfn = fullfile(QS.dir.pivSimAvg2x, 'vvM_simpletimeavg2x.mat')  ;
    tmp = load(vvsmMfn) ;
    vertex_vels = tmp.vvsmM ;
else
    vvsmMfn = fullfile(QS.dir.pivSimAvg, 'vvM_simpletimeavg.mat')  ; 
    tmp = load(vvsmMfn) ;
    vertex_vels = tmp.vvsmM ;   
end
% vfsmMfn = fullfile(QS.dir.pivSimAvg, 'vfM_simpletimeavg.mat') ;
% tmp = load(vfsmMfn) ;
% face_vels = tmp.vfsmM ;

%% Load time offset for first fold, t0
QS.t0set() ;
tfold = QS.t0 ;

%% load from QS
if doubleResolution
    nU = QS.nU * 2 - 1 ;
    nV = QS.nV * 2 - 1 ;
else
    nU = QS.nU ;
    nV = QS.nV ;    
end

%% Test incompressibility of the flow on the evolving surface
% We relate the normal velocities to the divergence / 2 * H.
tps = QS.xp.fileMeta.timePoints(1:end-1) - tfold;

% preallocate for cumulative error
ntps = length(QS.xp.fileMeta.timePoints(1:end-1)) ;
% HH_apM   = zeros(ntps, nU) ;   % dv averaged
% divv_apM = zeros(ntps, nU) ;
% veln_apM = zeros(ntps, nU) ;
% gdot_apM = zeros(ntps, nU) ;
% HH_lM   = zeros(ntps, nU) ;    % left averaged
% divv_lM = zeros(ntps, nU) ;
% veln_lM = zeros(ntps, nU) ;
% gdot_lM = zeros(ntps, nU) ;
% HH_rM   = zeros(ntps, nU) ;    % right averaged
% divv_rM = zeros(ntps, nU) ;
% veln_rM = zeros(ntps, nU) ;
% gdot_rM = zeros(ntps, nU) ;
% HH_dM   = zeros(ntps, nU) ;    % dorsal averaged
% divv_dM = zeros(ntps, nU) ;
% veln_dM = zeros(ntps, nU) ;
% gdot_dM = zeros(ntps, nU) ;
% HH_vM   = zeros(ntps, nU) ;    % ventral averaged
% divv_vM = zeros(ntps, nU) ;
% veln_vM = zeros(ntps, nU) ;
% gdot_vM = zeros(ntps, nU) ;

% Build timepoint list so that we first do every 10, then fill in details
lastIdx = length(QS.xp.fileMeta.timePoints) - 1 ;
coarseIdx = 1:10:lastIdx ;
fineIdx = setdiff(1:lastIdx, coarseIdx) ;
allIdx = [66, coarseIdx, fineIdx ] ;
tp2do = QS.xp.fileMeta.timePoints(allIdx) ;

% Output directory is inside metricKinematics dir
outdir = fullfile(mKDir, 'measurements') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Unit definitions for axis labels
unitstr = [ '[1/' QS.timeunits ']' ];
Hunitstr = [ '[1/' QS.spaceunits ']' ];
vunitstr = [ '[' QS.spaceunits '/' QS.timeunits ']' ];
    
% Compute or load all timepoints
for tp = tp2do
    close all
    disp(['t = ' num2str(tp)])
    tidx = QS.xp.tIdx(tp) ;

    % Check for timepoint measurement on disk
    Hfn = fullfile(outdir, sprintf('HH_series_%06d.mat', tp))   ;
    efn = fullfile(outdir, sprintf('gdot_series_%06d.mat', tp)) ;
    dfn = fullfile(outdir, sprintf('divv_series_%06d.mat', tp)) ;
    nfn = fullfile(outdir, sprintf('veln_series_%06d.mat', tp)) ;
    H2vnfn = fullfile(outdir, sprintf('H2vn_series_%06d.mat', tp)) ;

    redo_comp = overwrite || ~exist(Hfn, 'file') || ~exist(efn, 'file') ;
    redo_comp = redo_comp || ~exist(dfn, 'file') || ~exist(dfn, 'file') ;

    if redo_comp
        tic 
        if doubleResolution
            % Load current mesh
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC2x, tp)) ;
            mesh = tmp.spcutMeshSmRSC2x ;

            % Load cutMesh
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS2x, tp)) ;
            cutMesh = tmp.spcutMeshSmRS2x ;
            clearvars tmp
        else
            % Load current mesh
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp)) ;
            mesh = tmp.spcutMeshSmRSC ;

            % Load cutMesh
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
            cutMesh = tmp.spcutMeshSmRS ;
            clearvars tmp
        end
        
        % Compute mean curvature
        % Smooth the mesh with lambda_mesh
        if lambda_mesh > 0 
            disp('smoothing mesh vertices before computations')
            tri = triangulation(mesh.f, mesh.v) ;
            fbndy = tri.freeBoundary ;
            fbndy = fbndy(:, 1) ;
            mesh.v = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
                lambda_mesh, 'implicit', mesh.v) ;
        end
        DEC = DiscreteExteriorCalculus(mesh.f, mesh.v) ;
        H3d = sum(mesh.vn .* DEC.laplacian(mesh.v), 2) * 0.5 ;
        
        % OPTION1: Could compute divergence anew on smoothed mesh
        % fvel = squeeze(face_vels(tidx, :, :)) ;
        % divv3d_test = DEC.divergence(fvel) ;
        
        % OPTION2: Load divv from disk and convert using smoothed mesh
        % [~, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
        if doubleResolution
            dec_tp = load(sprintf(QS.fullFileBase.dec2x, tp)) ;
        else
            dec_tp = load(sprintf(QS.fullFileBase.dec, tp)) ;
        end
        
        % Smooth divergence(v) [divv]
        divv3d = laplacian_smooth(mesh.v, mesh.f, 'cotan', [],...
                                lambda, 'implicit', dec_tp.divs.raw) ;

        % veln = divv / (2H) ; 
        % veln_pred = DEC.divergence(fvel) ./ (2.0 * H) ;
        % veln is normal velocity field (velocities along vertex normals)
        % mesh.vn are vertex normals

        % Smooth the velocities in space using gptoolbox
        vx = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 1)) ;
        vy = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 2)) ;
        vz = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 3)) ;
        vxs = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
            lambda, 'implicit', vx') ;
        vys = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
            lambda, 'implicit', vy') ;
        vzs = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
            lambda, 'implicit', vz') ;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This should not be necessary, since vertex_vels should already
        % account for APDV framing! Debug this issue earlier
        % if QS.flipy
        %     vys = -vys ;
        % end
        % DEBUG THE FLIP
        % vv = squeeze(vertex_vels(tidx, 1:nU*(nV-1), :)) ;
        % trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
        %     vv(:, 2), 'edgecolor', 'none', 'facealpha', 0.6)
        % pause(1);
        % hold on;
        % quiver3(mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), vv(:, 1), vv(:, 2), vv(:, 3))
        % tmp2 = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp+1)) ;
        % m2 = tmp2.spcutMeshSmRSC ;
        % trisurf(m2.f, m2.v(:, 1), m2.v(:, 2), m2.v(:, 3), ...
        %     vv(:, 2), 'edgecolor', 'k', 'facealpha', 0.6)
        % axis equal
        % figure ; 
        % plot(mesh.v(:, 2) - m2.v(:, 2), vv(:, 2), '.')
        % clearvars tmp

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Actual normal velocity
        veln3d = sum(cat(2, vxs, vys, vzs) .* mesh.vn, 2) ;

        % Predict the divergence
        H2vn3d = 2 * H3d .* veln3d ;

        % Extend to have another row for 2d map
        divv2d = divv3d ;
        divv2d(nU*(nV-1) + 1:nU*nV) = divv3d(1:nU) ;
        H2d = H3d ;
        H2d(nU*(nV-1)+1:(nU*nV)) = H3d(1:nU) ;
        veln2d = veln3d ;
        veln2d(nU*(nV-1)+1:(nU*nV)) = veln3d(1:nU) ;
        H2vn2d = 2 * H2d .* veln2d ;

        % The difference
        if lambda_err > 0 
            gdot3d = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
                lambda_err, 'implicit', divv3d - H2vn3d) ;
            % expand error to 2d map
            gdot2d = gdot3d ;
            gdot2d(nU*(nV-1)+1:nU*nV) = gdot3d(1:nU) ;
        else
            gdot2d = divv2d - H2vn2d ;
            gdot3d = divv3d - H2vn3d ;
        end
    
        %% Store data on disk
        HH = reshape(H2d, [nU,nV]) ;
        gdot = reshape(gdot2d, [nU,nV]) ;
        divv = reshape(divv2d, [nU,nV]) ;
        veln = reshape(veln2d, [nU,nV]) ;
        H2vn = reshape(H2d .* veln2d, [nU, nV]) ;
        
        % Average along DV -- ignore last redudant row at nV
        HH_ap = mean(HH(:, 1:nV-1), 2) ;
        gdot_ap = mean(gdot(:, 1:nV-1), 2) ;
        divv_ap = mean(divv(:, 1:nV-1), 2) ;
        veln_ap = mean(veln(:, 1:nV-1), 2) ;
        H2vn_ap = mean(H2vn(:, 1:nV-1), 2) ;
        
        % quarter bounds
        q0 = round(nV * 0.125) ;
        q1 = round(nV * 0.375) ;
        q2 = round(nV * 0.625) ;
        q3 = round(nV * 0.875) ;
        left = q0:q1 ;
        ventral = q1:q2 ;
        right = q2:q3 ;
        dorsal = [q3:nV, 1:q1] ;
        
        % left quarter
        HH_l = mean(HH(:, left), 2) ;
        gdot_l = mean(gdot(:, left), 2) ;
        divv_l = mean(divv(:, left), 2) ;
        veln_l = mean(veln(:, left), 2) ;
        H2vn_l = mean(H2vn(:, left), 2) ;
        
        % right quarter
        HH_r = mean(HH(:, right), 2) ;
        gdot_r = mean(gdot(:, right), 2) ;
        divv_r = mean(divv(:, right), 2) ;
        veln_r = mean(veln(:, right), 2) ;
        H2vn_r = mean(H2vn(:, right), 2) ;
        
        % dorsal quarter
        HH_d = mean(HH(:, dorsal), 2) ;
        gdot_d = mean(gdot(:, dorsal), 2) ;
        divv_d = mean(divv(:, dorsal), 2) ;
        veln_d = mean(veln(:, dorsal), 2) ;
        H2vn_d = mean(H2vn(:, dorsal), 2) ;
        
        % ventral quarter
        HH_v = mean(HH(:, ventral), 2) ;
        gdot_v = mean(gdot(:, ventral), 2) ;
        divv_v = mean(divv(:, ventral), 2) ;
        veln_v = mean(veln(:, ventral), 2) ;
        H2vn_v = mean(H2vn(:, ventral), 2) ;
        
        %% Save timeseries measurements
        save(Hfn, 'HH', 'HH_ap', 'HH_l', 'HH_r', 'HH_d', 'HH_v')
        save(efn, 'gdot', 'gdot_ap', 'gdot_l', 'gdot_r', 'gdot_d', 'gdot_v')
        save(dfn, 'divv', 'divv_ap', 'divv_l', 'divv_r', 'divv_d', 'divv_v')
        save(nfn, 'veln', 'veln_ap', 'veln_l', 'veln_r', 'veln_d', 'veln_v') 
        save(H2vnfn, 'H2vn', 'H2vn_ap', 'H2vn_l', 'H2vn_r', 'H2vn_d', 'H2vn_v')
        
        save_lambdas = true ;
        toc
        
    else
        % Load timeseries measurements
        load(Hfn, 'HH', 'HH_ap', 'HH_l', 'HH_r', 'HH_d', 'HH_v')
        load(efn, 'gdot', 'gdot_ap', 'gdot_l', 'gdot_r', 'gdot_d', 'gdot_v')
        load(dfn, 'divv', 'divv_ap', 'divv_l', 'divv_r', 'divv_d', 'divv_v')
        load(nfn, 'veln', 'veln_ap', 'veln_l', 'veln_r', 'veln_d', 'veln_v') 
        load(H2vnfn, 'H2vn', 'H2vn_ap', 'H2vn_l', 'H2vn_r', 'H2vn_d', 'H2vn_v') 
        
        H2d = HH ;
        H3d = HH(:, 1:nV-1) ;
        gdot2d = gdot ;
        gdot3d = gdot(:, 1:nV-1);
        divv2d = divv ;
        divv3d = divv(:, 1:nV-1) ;
        veln2d = veln ;
        veln3d = veln(:, 1:nV-1) ;
        H2vn2d = H2vn ;
        H2vn3d = H2vn(:, 1:nV-1) ;
        
        % Ensure that mesh and cutMesh are not assumed to be some other
        % timepoint's meshes
        mesh = []; 
        cutMesh = [] ;
        save_lambdas = false ;
    end    

    %% Save lambdas used to disk
    if save_lambdas    
        lambs = [lambda, lambda_mesh, lambda_err] ;
        header = ['laplacian smoothing parameters: lambda for ', ...
            'velocity & div(v), lambda for mesh, lambda for residual (gdot)'] ;
        filename = fullfile(mKDir, 'lambdas.txt') ;
        write_txt_with_header(filename, lambs, header)
    end
    
    %% Plot results
    % operational plotting options
    pOptions.overwrite = overwrite ;
    pOptions.plot_flows = plot_flows ;
    pOptions.plot_Hgdot = plot_Hgdot ;
    pOptions.plot_factors = plot_factors ;
    % parameter plotting options
    pOptions.doubleResolution = doubleResolution; 
    pOptions.lambda = lambda ;
    pOptions.lambda_err = lambda_err ;
    pOptions.lambda_mesh = lambda_mesh ;
    pOptions.H2vn2d = H2vn2d ;
    pOptions.divv2d = divv2d ;
    pOptions.gdot2d = gdot2d ;
    pOptions.veln2d = veln2d ;
    pOptions.H2d = H2d ;
    pOptions.H2vn3d = H2vn3d ;
    pOptions.divv3d = divv3d ;
    pOptions.gdot3d = gdot3d ;
    pOptions.veln3d = veln3d ;
    pOptions.H3d = H3d ;
    pOptions.cutMesh = cutMesh ;
    pOptions.mesh = mesh ;
    pOptions.tp = tp ;
    pOptions.climit = climit ;
    pOptions.climit_err = climit ;
    pOptions.climit_veln = climit_veln ;
    pOptions.climit_H = climit_H ;
    QS.plotMetricKinematics(pOptions)
    
    % %% Store in matrices
    % % dv averaged
    % HH_apM(tidx, :) = HH_ap ;
    % gdot_apM(tidx, :) = gdot_ap ;
    % divv_apM(tidx, :) = divv_ap ;
    % veln_apM(tidx, :) = veln_ap ;
    % H2vn_apM(tidx, :) = H2vn_ap ;
    % 
    % % left quarter
    % HH_lM(tidx, :) = HH_l ;
    % gdot_lM(tidx, :) = gdot_l ;
    % divv_lM(tidx, :) = divv_l ;
    % veln_lM(tidx, :) = veln_l ;
    % H2vn_lM(tidx, :) = H2vn_l ;
    % 
    % % right quarter
    % HH_rM(tidx, :) = HH_r ;
    % gdot_rM(tidx, :) = gdot_r ;
    % divv_rM(tidx, :) = divv_r ;
    % veln_rM(tidx, :) = veln_r ;
    % H2vn_rM(tidx, :) = H2vn_r ;
    % 
    % % dorsal quarter
    % HH_dM(tidx, :) = HH_d ;
    % gdot_dM(tidx, :) = gdot_d ;
    % divv_dM(tidx, :) = divv_d ;
    % veln_dM(tidx, :) = veln_d ;
    % H2vn_dM(tidx, :) = H2vn_d ;
    % 
    % % ventral quarter
    % HH_vM(tidx, :) = HH_v ;
    % gdot_vM(tidx, :) = gdot_v ;
    % divv_vM(tidx, :) = divv_v ;
    % veln_vM(tidx, :) = veln_v ;
    % H2vn_vM(tidx, :) = H2vn_v ;
end

