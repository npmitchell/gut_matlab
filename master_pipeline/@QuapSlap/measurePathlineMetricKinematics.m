function measurePathlineMetricKinematics(QS, options)
% measurePathlineMetricKinematics(QS, options)
%   Query the metric Kinematics along lagrangian pathlines.
%   Plot results as kymographs and correlation plots.
%   Out-of-plane motion is v_n * 2H, where v_n is normal velocity and H is
%   mean curvature.
%   In-plane motion considered here is div(v_t) where v_t is tangential
%   velocity on the curved surface.
%   The difference div(v_t) - vn*2H = Tr[g^{-1} dot{g}], which is a measure
%   of isotropic metric change over time (dot is the partial derivative wrt
%   time). 
% 
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   plot_kymographs : bool
%   plot_kymographs_cumsum : bool
%   plot_correlations : bool
%   plot_gdot_correlations : bool
%   plot_gdot_decomp : bool
% 
% NPMitchell 2020

%% Default options 
overwrite = false ;
plot_kymographs = true ;
plot_kymographs_cumsum = true ;
plot_correlations = true ;
plot_gdot_correlations = false ;
plot_gdot_decomp = true ;

%% Parameter options
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
%% parameter options
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

%% Operational options
if isfield(options, 'plot_kymographs')
    plot_kymographs = options.plot_kymographs ;
end
if isfield(options, 'plot_kymographs_cumsum')
    plot_kymographs_cumsum = options.plot_kymographs_cumsum ;
end
if isfield(options, 'plot_correlations')
    plot_correlations = options.plot_correlations ;
end
if isfield(options, 'plot_gdot_correlations')
    plot_gdot_correlations = options.plot_gdot_correlations ;
end
if isfield(options, 'plot_gdot_decomp')
    plot_gdot_decomp = options.plot_gdot_decomp ;
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

%% Build timepoint list so that we first do every 10, then fill in details
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

    % Load timeseries measurements
    load(Hfn, 'HH', 'HH_ap', 'HH_l', 'HH_r', 'HH_d', 'HH_v')
    load(efn, 'gdot', 'gdot_ap', 'gdot_l', 'gdot_r', 'gdot_d', 'gdot_v')
    load(dfn, 'divv', 'divv_ap', 'divv_l', 'divv_r', 'divv_d', 'divv_v')
    load(nfn, 'veln', 'veln_ap', 'veln_l', 'veln_r', 'veln_d', 'veln_v') 
    load(H2vnfn, 'H2vn', 'H2vn_ap', 'H2vn_l', 'H2vn_r', 'H2vn_d', 'H2vn_v') 
end

    
%% Kymographs of pathlines
t0 = QS.t0set() ;
QS.loadPullbackPathlines('vertexPathlines')
vP = QS.pathlines.vertices ;

% Query each pathline position in relevant mesh and interpolate divv
for 