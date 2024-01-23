function measurePathlineVorticity(QS, options)
% measurePathlineVorticity(QS, options)
%   Query the covariant curl along lagrangian pathlines.
%   Plot results as kymographs and correlation plots.
%   
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   plot_kymographs : bool
% 
% NPMitchell 2022

%% Default options 
overwrite = false ;

%% Parameter options
lambda = QS.smoothing.lambda ; 
lambda_mesh = QS.smoothing.lambda_mesh ;
nmodes = QS.smoothing.nmodes ;
zwidth = QS.smoothing.zwidth ;
climit = 0.2 ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
averagingStyle = "Lagrangian" ;
QS.t0set() ;
t0Pathline = QS.t0 ;
nU = QS.nU ;
nV = QS.nV ;

%% Unpack options & assign defaults
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
%% parameter options
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
end
if isfield(options, 'nmodes')
    nmodes = options.nmodes ;
end
if isfield(options, 'zwidth')
    zwidth = options.zwidth ;
end
if isfield(options, 'climit')
    climit = options.climit ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'averagingStyle')
    averagingStyle = options.averagingStyle ;
end

%% Operational options
if isfield(options, 'plot_kymographs')
    plot_kymographs = options.plot_kymographs ;
end
if isfield(options, 't0Pathline')
    t0Pathline = options.t0Pathline ;
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
folds = load(QS.fileName.fold) ;
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;

%% Colormap
close all
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = brewermap(256, '*RdBu') ;
bbr256 = blueblackred(256) ;
close all

%% Build timepoint list so that we first do every 10, then fill in details
lastIdx = length(QS.xp.fileMeta.timePoints) - 1 ;
coarseIdx = 1:10:lastIdx ;
fineIdx = setdiff(1:lastIdx, coarseIdx) ;
allIdx = [coarseIdx, fineIdx ] ;
tp2do = QS.xp.fileMeta.timePoints(allIdx) ;

% DONE WITH PREPARATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load pathlines to build Kymographs along pathlines
QS.loadPullbackPathlines(t0Pathline, 'vertexPathlines')
vP = QS.pathlines.vertices ;

% Output directory is inside metricKinematics dir
outdir = sprintf(QS.dir.pathlines.vorticity, t0Pathline) ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Load Lx, Ly by loadingPIV
QS.loadPIV()
Xpiv = QS.piv.raw.x ;
Ypiv = QS.piv.raw.y ;

% Discern if piv measurements are done on a double covering or the meshes
if strcmp(QS.piv.imCoords(end), 'e')
    doubleCovered = true ;
end

% Compute or load all timepoints
for tp = tp2do
    close all
    disp(['t = ' num2str(tp)])
    tidx = QS.xp.tIdx(tp) ;
    QS.setTime(tp) ;
    
    % Check for timepoint measurement on disk
    ofn = fullfile(outdir, sprintf('vorticity_pathline%04d_%06d.mat', t0Pathline, tp))   ;
    files_missing = ~exist(ofn, 'file')  ;
    
    if overwrite || files_missing
        disp('Computing pathline vorticity...')
        % Load timeseries measurements defined on mesh vertices
        % mdatdir = QS.dir.piv.avgDEC.rot_smoothed ;
        % ofnMesh = fullfile(mdatdir, sprintf('vorticity_vertices_%06d.mat', tp))   ;
        %
        % try
        %     load(ofnMesh, 'vorticity_filt')
        %     vorticity = vorticity_filt ;
        % catch
         
        disp('Computing smoothed vorticity from piv/LagrangianAvg/dec/')
        if lambda > 0 || lambda_mesh>0
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
                % tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
                % cutMesh = tmp.spcutMeshSmRS ;
                % clearvars tmp
            end
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
        
        %% OPTION2: Load divv from disk and convert using smoothed mesh
        % [~, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
        if strcmp(averagingStyle, 'Lagrangian')
            if doubleResolution
                dec_tp = load(sprintf(QS.fullFileBase.decAvg2x, tp)) ;
            else
                dec_tp = load(sprintf(QS.fullFileBase.decAvg, tp)) ;
            end
        else
            if doubleResolution
                dec_tp = load(sprintf(QS.fullFileBase.decSimAvg2x, tp)) ;
            else
                dec_tp = load(sprintf(QS.fullFileBase.decSimAvg, tp)) ;
            end
        end

        % Smooth divergence(v) [divv]
        if lambda > 0
            vorticity3d = reshape(dec_tp.rots.rotv, [nU, nV]) ;
            vorticity3d = vorticity3d(:, 1:end-1) ;
            vorticity3d = laplacian_smooth(mesh.v, mesh.f, 'cotan', [],...
                                    lambda, 'implicit', vorticity3d(:)) ;
        else
            vorticity3d = dec_tp.rots.rotv ;
        end

        % Extend to have another row for 2d map
        vorticity2d = vorticity3d ;
        vorticity2d(nU*(nV-1)+1:(nU*nV)) = vorticity3d(1:nU) ;

        %% Store data on disk
        vorticity = reshape(vorticity2d, [nU,nV]) ;

        % Bandpass filter modes
        if nmodes > 0
            filterOptions.nmodes = nmodes ;
            filterOptions.zwidth = zwidth ;
            vorticity_filt = modeFilterQuasi1D(vorticity(:, 1:nV-1), filterOptions) ;   

            % 3d version does not have periodic row
            % vorticity3d_filt = vorticity_filt ;

            % Append periodic row
            vorticity = vorticity_filt ;
            vorticity(:, nV) = vorticity_filt(:, 1) ;

        end
        % end
        
        % Interpolate from vertices onto pathlines
        xx = vP.vX(tidx, :, :) ;
        yy = vP.vY(tidx, :, :) ;
        XY = [xx(:), yy(:)] ;
        Lx = vP.Lx(tidx) ;
        Ly = vP.Ly(tidx) ;
        options.Lx = Lx ;
        options.Ly = Ly ;
        XY = QS.doubleToSingleCover(XY, Ly) ;
        vorticity = QS.interpolateOntoPullbackXY(XY, vorticity, options) ;
        
        % OPTION 1: simply reshape, tracing each XY dot to its t0Pathline
        % grid coordinate
        vorticity = reshape(vorticity, [nU, nV]) ;
        
        %% OPTION 2: the following regrids onto original XY coordinates,
        % rendering the process of following pathlines moot. 
        % Average into AP bins and take mean along 1/4 DV hoop arcs
        % if doubleCovered
        %     vminmax = [0.25 * Ly, 0.75 * Ly] ;
        % else
        %     vminmax = [1, Ly] ;
        % end
        %
        % Note the transposition: to plot as APDV, imshow(m')
        % vorticity = binData2dGrid([XY, vorticity], [1,Lx], vminmax, nU, nV) ;
           
        %% Average along DV -- do not ignore last row at nV since not quite
        % redundant in this version of the algorithm -- is that true?
        vorticity_ap = nanmean(vorticity, 2) ;
        
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
        vorticity_l = nanmean(vorticity(:, left), 2) ;
        
        % right quarter
        vorticity_r = nanmean(vorticity(:, right), 2) ;
        
        % dorsal quarter
        vorticity_d = nanmean(vorticity(:, dorsal), 2) ;
        
        % ventral quarter
        vorticity_v = nanmean(vorticity(:, ventral), 2) ;
        
        % Save results
        save(ofn, 'vorticity', 'vorticity_ap', 'vorticity_l', ...
            'vorticity_r', 'vorticity_d', 'vorticity_v', ...
            'lambda', 'lambda_mesh')
    end
end
disp('done with measuring pathline metric kinematics')

%% Combine DV-averaged profiles into kymographs
apKymoFn = fullfile(outdir, 'apKymographVorticity.mat') ;
lKymoFn = fullfile(outdir, 'leftKymographVorticity.mat') ;
rKymoFn = fullfile(outdir, 'rightKymographVorticity.mat') ;
dKymoFn = fullfile(outdir, 'dorsalKymographVorticity.mat') ;
vKymoFn = fullfile(outdir, 'ventralKymographVorticity.mat') ;
files_exist = exist(apKymoFn, 'file') && ...
    exist(lKymoFn, 'file') && exist(rKymoFn, 'file') && ...
    exist(dKymoFn, 'file') && exist(vKymoFn, 'file') ;
if ~files_exist || overwrite
    for tp = QS.xp.fileMeta.timePoints(1:end-1)
        close all
        disp(['t = ' num2str(tp)])
        tidx = QS.xp.tIdx(tp) ;

        % Check for timepoint measurement on disk
        ofn = fullfile(outdir, sprintf('vorticity_pathline%04d_%06d.mat', t0Pathline, tp))   ;

        % Load timeseries measurements
        load(ofn, 'vorticity', 'vorticity_ap', 'vorticity_l', 'vorticity_r', 'vorticity_d', 'vorticity_v')

        %% Store in matrices
        % dv averaged
        vorticity_apM(tidx, :) = vorticity_ap ;

        % left quarter
        vorticity_lM(tidx, :) = vorticity_l ;

        % right quarter
        vorticity_rM(tidx, :) = vorticity_r ;
        % dorsal quarter
        vorticity_dM(tidx, :) = vorticity_d ;
        % ventral quarter
        vorticity_vM(tidx, :) = vorticity_v ;
    end
    
    disp('Saving DV-averaged kymograph data')
    % Save the DV-averaged kymographs
    save(apKymoFn, 'vorticity_apM')
    save(lKymoFn, 'vorticity_lM')
    save(rKymoFn, 'vorticity_rM')
    save(dKymoFn, 'vorticity_dM')
    save(vKymoFn, 'vorticity_vM')
    
end
disp('done measuring pathline vorticity (rot(vt))')