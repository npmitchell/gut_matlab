function initializeQuapSlap(QS, xp, opts)
%initializeQuapSlap(QS, xp, opts)
%   Hidden method for instantiating QuapSlap class
%
% Parameters
% ----------
% QS : QuapSlap object whose properties to fill in
% xp : Imsane Experiment class instance belonging to QS
% opts : struct with fields
%   xp : ImSAnE class object instance
%   flipy : bool
%   meshDir : str
%   timeUnits : str
%   spaceUnits : str
%   nU : int
%   nV : int
%   lambda : optional float 
%   lambda_mesh : optional float 
%   lambda_err : optional float
%   
%
% NPMitchell 2020

%% PROPERTIES
QS.xp = xp ;
QS.flipy = opts.flipy ;
meshDir = opts.meshDir ;
QS.timeUnits = opts.timeUnits ;
QS.spaceUnits = opts.spaceUnits ;
QS.fileBase.fn = xp.fileMeta.filenameFormat ;
QS.ssfactor = xp.detectOptions(1).ssfactor ;
QS.nU = opts.nU ;
QS.nV = opts.nV ;
if isfield(opts, 'timeinterval')
    QS.timeInterval = opts.timeInterval ;
else
    QS.timeInterval = 1 ;
end

QS.normalShift = opts.normalShift ;
QS.a_fixed = opts.a_fixed ;
if isfield(opts, 'adjustlow')
    QS.data.adjustlow = opts.adjustlow ;
end
if isfield(opts, 'adjusthigh')
    QS.data.adjusthigh = opts.adjusthigh ;
end
if isfield(opts, 'axisOrder')
    QS.data.axisOrder = opts.axisOrder ;
end
if isfield(opts, 'ilastikOutputAxisOrder')
    QS.data.ilastikOutputAxisOrder = opts.ilastikOutputAxisOrder ;
end

% Assign which pullback coordsys is used for velocimetry
if isfield(opts, 'pivPullback')
    QS.pivPullback = pivPullback ;
end

%% NAMING
uvexten = sprintf('_nU%04d_nV%04d', QS.nU, QS.nV) ;

% Naming extension
QS.uvexten = uvexten ;

% APDV coordinate system
QS.APDV.resolution = min(xp.fileMeta.stackResolution) ;
QS.APDV.rot = [] ;
QS.APDV.trans = [] ;

% Plotting settings
QS.plotting.colors = define_colors() ;
QS.plotting.markers = {'o', 's', '^', 'v', '*', '>', '<'} ;

% Directories of the QS object
% Meshes and measurements before gridding into pullback coords
QS.dir.data = xp.fileMeta.dataDir ;
QS.dir.mesh = meshDir ;
QS.dir.maskedData = fullfile(meshDir, 'masked_data') ;
QS.dir.alignedMesh = fullfile(meshDir, 'aligned_meshes') ;
QS.dir.cntrline = fullfile(meshDir, 'centerline') ;
QS.dir.cylinderMesh = fullfile(meshDir, 'cylinder_meshes') ;
QS.dir.cutMesh = fullfile(meshDir, 'cutMesh') ;
QS.dir.cylinderMeshClean = fullfile(QS.dir.cylinderMesh, 'cleaned') ;
QS.dir.texturePatchIm = fullfile(meshDir, 'images_texturepatch') ;

% After gridding into (u,v) / (zeta,phi) pullback coords
uvDir = fullfile(QS.dir.mesh, sprintf('gridCoords_nU%04d_nV%04d', QS.nU, QS.nV)) ;
QS.dir.uvCoord = uvDir ;

% define string for smoothing params
if isfield(opts, 'lambda')
    QS.smoothing.lambda = opts.lambda ;
end
if isfield(opts, 'lambda_mesh')
    QS.smoothing.lambda_mesh = opts.lambda_mesh ;
end
if isfield(opts, 'lambda_err')
    QS.smoothing.lambda_err = opts.lambda_err ;
end
l_lmesh_lerr =  strrep(sprintf('lambda%0.03f_lmesh%0.3f_lerr%0.3f', ...
    QS.smoothing.lambda, QS.smoothing.lambda_mesh, ...
    QS.smoothing.lambda_err), '.', 'p') ;
l_lmesh = strrep(sprintf('lambda%0.03f_lmesh%0.3f', ...
    QS.smoothing.lambda, QS.smoothing.lambda_mesh), '.', 'p') ;

% Metric strain dirs
QS.dir.metricKinematics = struct() ;
QS.dir.metricKinematics.root = fullfile(uvDir, 'metricKinematics') ;
QS.dir.metricKinematics.smoothing = fullfile(uvDir, 'metricKinematics', ...
    l_lmesh_lerr) ;
QS.dir.metricKinematics.measurements = ...
    fullfile(uvDir, 'metricKinematics', l_lmesh_lerr, 'measurements') ;
QS.fullFileBase.metricKinematics = struct() ;
QS.fullFileBase.metricKinematics.divv = ...
    fullfile(QS.dir.metricKinematics.measurements, 'divv_vertices_%06d') ;
QS.fullFileBase.metricKinematics.H2vn = ...
    fullfile(QS.dir.metricKinematics.measurements, 'H2vn_vertices_%06d') ;
QS.fullFileBase.metricKinematics.gdot = ...
    fullfile(QS.dir.metricKinematics.measurements, 'gdot_vertices_%06d') ;

% Metric Kinematics along pathlines
QS.dir.metricKinematics.pathline = struct() ;
QS.dir.metricKinematics.pathline.root = ...
    fullfile(uvDir, 'metricKinematics', l_lmesh_lerr, 'pathline_%04dt0') ;
QS.dir.metricKinematics.pathline.measurements = ...
    fullfile(QS.dir.metricKinematics.pathline.root, 'measurements') ;
QS.fileBase.metricKinematics = struct() ;
QS.fileBase.metricKinematics.pathline.Hfn = 'HH_pathline%04d_%06d.mat' ;
QS.fileBase.metricKinematics.pathline.gdot = 'gdot_pathline%04d_%06d.mat' ;
QS.fileBase.metricKinematics.pathline.divv = 'divv_pathline%04d_%06d.mat' ;
QS.fileBase.metricKinematics.pathline.veln = 'veln_pathline%04d_%06d.mat' ;
QS.fileBase.metricKinematics.pathline.radius = 'radius_pathline%04d_%06d.mat' ;
QS.fileBase.metricKinematics.pathline.kymographs = struct() ;
QS.fileBase.metricKinematics.pathline.kymographs.ap = 'apKymographMetricKinematics.mat' ;
QS.fileBase.metricKinematics.pathline.kymographs.d = 'dKymographMetricKinematics.mat' ;
QS.fileBase.metricKinematics.pathline.kymographs.v = 'vKymographMetricKinematics.mat' ;
QS.fileBase.metricKinematics.pathline.kymographs.l = 'lKymographMetricKinematics.mat' ;
QS.fileBase.metricKinematics.pathline.kymographs.r = 'rKymographMetricKinematics.mat' ;

% Metric deformation ('metric strain' or 'meshMetric)
QS.dir.gstrain = fullfile(uvDir, 'metricStrain') ;
QS.dir.gstrainRate = fullfile(QS.dir.gstrain, 'rateMetric') ;
QS.dir.gstrainRateIm = fullfile(QS.dir.gstrainRate, 'images') ;
QS.dir.gstrainMesh = fullfile(QS.dir.gstrain, 'meshMetric') ;
QS.dir.gstrainMeshIm = fullfile(QS.dir.gstrainMesh, 'images') ;
QS.fileBase.gstrainMesh = 'gstrainMesh_%06d.mat' ;
QS.fileBase.gstrainRate = 'gstrainRate_%06d.mat' ;
QS.fullFileBase.gstrainMesh = fullfile(QS.dir.gstrainMesh, ...
    QS.fileBase.gstrainMesh) ; 
QS.fullFileBase.gstrainRate = fullfile(QS.dir.gstrainRate, ...
    QS.fileBase.gstrainRate) ; 

%% Strain Rate
QS.dir.strainRate = struct() ;
QS.dir.strainRate.root = fullfile(uvDir, 'strainRate') ;
QS.dir.strainRate.smoothing = fullfile(uvDir, 'strainRate', l_lmesh) ;
QS.dir.strainRate.measurements = ...
    fullfile(uvDir, 'strainRate', l_lmesh, 'measurements') ;
QS.fileBase.strainRate = 'strainRate_%06d.mat' ;
QS.fullFileBase.strainRate = fullfile(QS.dir.strainRate.measurements, ...
    QS.fileBase.strainRate) ; 
QS.fileBase.strain = 'strain_%06d.mat' ;
QS.fullFileBase.strain = fullfile(QS.dir.strainRate.measurements, ...
    QS.fileBase.strain) ; 
% Strain rates along pathlines -- to be filled in with smoothing and t0
QS.dir.strainRate.pathline = struct() ;
QS.dir.strainRate.pathline.root = ...
    fullfile(uvDir, 'strainRate', l_lmesh, 'pathline_%04dt0') ;
QS.dir.strainRate.pathline.measurements = ...
    fullfile(QS.dir.strainRate.pathline.root, 'measurements') ;

% shorten variable names for brevity
clineDir = QS.dir.cntrline ;

% fileBases
QS.fileBase.name = xp.fileMeta.filenameFormat(1:end-4) ;
QS.fileBase.mesh = ...
    [xp.detector.options.ofn_smoothply '%06d'] ;
QS.fileBase.alignedMesh = ...
    [QS.fileBase.mesh '_APDV_um'] ;
QS.fileBase.apdProb = [QS.fileBase.name '_Probabilities_apcenterline.h5'] ; ;
QS.fileBase.centerlineXYZ = ...
    [QS.fileBase.mesh '_centerline_exp1p0_res*.txt' ] ;
QS.fileBase.centerlineAPDV = ...
    [QS.fileBase.mesh '_centerline_scaled_exp1p0_res*.txt' ] ;
QS.fileBase.cylinderMesh = ...
    [QS.fileBase.mesh '_cylindercut.ply'] ;
QS.fileBase.apBoundary = 'ap_boundary_indices_%06d.mat';
QS.fileBase.cylinderKeep = 'cylinderMesh_keep_indx_%06d.mat' ;
QS.fileName.apBoundaryDorsalPts = 'ap_boundary_dorsalpts.h5' ;

%% Clean Cylinder Mesh
QS.fileName.aBoundaryDorsalPtsClean = ...
    fullfile(QS.dir.cylinderMeshClean, 'adIDx.h5') ;
QS.fileName.pBoundaryDorsalPtsClean = ...
    fullfile(QS.dir.cylinderMeshClean, 'pdIDx.h5') ;

%% cutMesh
QS.fullFileBase.cutPath = fullfile(QS.dir.cutMesh, 'cutPaths_%06d.txt') ;

%% fileNames
nshift = strrep(sprintf('%03d', QS.normalShift), '-', 'n') ;
shiftstr = ['_' nshift 'step'] ;
QS.fileName.rot = fullfile(meshDir, 'rotation_APDV.txt') ;
QS.fileName.trans = fullfile(meshDir, 'translation_APDV.txt') ;
QS.fileName.xyzlim_raw = fullfile(meshDir, 'xyzlim_raw.txt') ;
QS.fileName.xyzlim_pix = fullfile(meshDir, 'xyzlim_APDV.txt') ;
QS.fileName.xyzlim_um = ...
    fullfile(meshDir, 'xyzlim_APDV_um.txt') ;
QS.fileName.xyzlim_um_buff = ...
    fullfile(meshDir, ['xyzlim_APDV_um' shiftstr '.txt']) ;
% fileNames for APDV and cylinderMesh
QS.fullFileBase.apdProb = fullfile(QS.dir.data, QS.fileBase.apdProb) ;
QS.fileName.apdv = ...
    fullfile(clineDir, 'apdv_coms_from_training.h5') ;
QS.fileName.dcom = fullfile(meshDir, 'dcom_for_rot.txt') ;
QS.fileName.startendPt = fullfile(clineDir, 'startendpt.h5') ;
QS.fileName.cleanCntrlines = ...
    fullfile(clineDir, 'centerlines_anomalies_fixed.mat') ;
QS.fileName.apBoundaryDorsalPts = ...
    fullfile(QS.dir.cylinderMesh, 'ap_boundary_dorsalpts.h5') ;
QS.fileName.endcapOptions = ...
    fullfile(QS.dir.cylinderMesh, 'endcapOptions.mat') ;
QS.fileName.apdBoundary = ...
    fullfile(QS.dir.cylinderMesh, 'ap_boundary_dorsalpts.h5') ;

%% FileNamePatterns
QS.fullFileBase.mesh = ...
    fullfile(QS.dir.mesh, [QS.fileBase.mesh '.ply']) ;
QS.fullFileBase.alignedMesh = ...
    fullfile(QS.dir.alignedMesh, [QS.fileBase.alignedMesh '.ply']) ;
% fileNames for centerlines
QS.fullFileBase.centerlineXYZ = ...
    fullfile(clineDir, QS.fileBase.centerlineXYZ) ;
QS.fullFileBase.centerlineAPDV = ...
    fullfile(clineDir, QS.fileBase.centerlineAPDV) ;
QS.fullFileBase.cylinderMesh = ...
    fullfile(QS.dir.cylinderMesh, QS.fileBase.cylinderMesh) ;
QS.fullFileBase.apBoundary = ...
    fullfile(QS.dir.cylinderMesh, QS.fileBase.apBoundary) ;
% QS.fullFileBase.apBoundaryDorsalPts = QS.fileName.apBoundaryDorsalPts ;
QS.fullFileBase.cylinderKeep = ...
    fullfile(QS.dir.cylinderMesh, QS.fileBase.cylinderKeep) ;
QS.fullFileBase.cylinderMeshClean = ...
    fullfile(QS.dir.cylinderMesh, 'cleaned',...
    [QS.fileBase.mesh '_cylindercut_clean.ply']) ;            

%% Define cutMesh directories
% cutFolder = fullfile(meshDir, 'cutMesh') ;
% cutMeshBase = fullfile(cutFolder, [QS.fileBase.name, '_cutMesh.mat']) ;
imFolderBase = fullfile(uvDir, ['PullbackImages' shiftstr] ) ;
sphiDir = fullfile(uvDir, ['sphi_cutMesh' shiftstr]) ;
sphiSmDir = fullfile(sphiDir, 'smoothed') ;
sphiSmRSDir = fullfile(sphiDir, 'smoothed_rs') ;
% sphiSmRSImDir = fullfile(sphiSmRSDir, 'images') ;
% sphiSmRSPhiImDir = fullfile(sphiSmRSImDir, 'phicolor') ;
sphiSmRSCDir = fullfile(sphiDir, 'smoothed_rs_closed') ;
sphiSmDir2x = fullfile(sphiDir, 'smoothed_doubleResolution') ;
sphiSmRSDir2x = fullfile(sphiDir, 'smoothed_rs_doubleResolution') ;
sphiSmRSCDir2x = fullfile(sphiDir, 'smoothed_rs_closed_doubleResolution') ;

%% Images of pullbacks in smoothed mesh coordinate systems
imFolder_sp = [imFolderBase '_sphi'] ;
imFolder_spe = fullfile(imFolder_sp, 'extended') ;
imFolder_up = [imFolderBase '_uphi'] ;
imFolder_upe = fullfile(imFolder_up, 'extended') ;
% time-averaged meshes
imFolder_spsm = fullfile(imFolder_sp, 'smoothed') ;
imFolder_spsme = fullfile(imFolder_sp, 'smoothed_extended') ;  % raw LUT, no histeq
imFolder_spsmeLUT = fullfile(imFolder_sp, 'smoothed_extended_LUT') ;  % with histeq
imFolder_rsm = fullfile([imFolderBase, '_sphi_relaxed'], 'smoothed');
imFolder_rsme = fullfile([imFolderBase, '_sphi_relaxed'], 'smoothed_extended') ;
imFolder_rsme_stack = fullfile([imFolderBase, '_sphi_relaxed'], ...
    'smoothed_extended_stack') ;  % with histeq?

% Lobe/fold identification paths
lobeDir = fullfile(uvDir, 'lobes') ;
foldHoopImDir = fullfile(lobeDir, 'constriction_hoops') ;
% Folder for curvature measurements
KHSmDir = fullfile(sphiSmRSCDir, 'curvature') ;
KSmDir = fullfile(KHSmDir, 'gauss') ;
HSmDir = fullfile(KHSmDir, 'mean') ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Port into QS
QS.dir.cutFolder = fullfile(meshDir, 'cutMesh') ;
QS.dir.spcutMesh = sphiDir ;
QS.dir.spcutMeshSm = sphiSmDir ;
QS.dir.spcutMeshSmRS = sphiSmRSDir ;
QS.dir.spcutMeshSmRSC = sphiSmRSCDir ;

%% double resolution cutMeshes
QS.dir.spcutMeshSm2x = sphiSmDir2x ;
QS.dir.spcutMeshSmRS2x = sphiSmRSDir2x ;
QS.dir.spcutMeshSmRSC2x = sphiSmRSCDir2x ;

%% centerlines
QS.dir.clineDVhoop = ...
    fullfile(QS.dir.uvCoord, ...
    ['centerline_from_DVhoops' shiftstr]) ;
QS.dir.writhe =  fullfile(QS.dir.clineDVhoop, 'writhe') ;

%% Images
QS.dir.im_uv = [imFolderBase '_uv'] ;
QS.dir.im_uve = [imFolderBase '_uv_extended'] ;
QS.dir.im_r = [imFolderBase '_sphi_relaxed'] ;
QS.dir.im_re = fullfile(QS.dir.im_r, 'extended') ;
QS.dir.im_sp = imFolder_sp ;
QS.dir.im_spe = imFolder_spe ;
QS.dir.im_up = imFolder_up ;
QS.dir.im_upe = imFolder_upe ;
QS.dir.im_sp_sm = imFolder_spsm ;
QS.dir.im_sp_sme = imFolder_spsme ;
QS.dir.im_sp_smeLUT = imFolder_spsmeLUT ;
QS.dir.im_r_sm = imFolder_rsm ;
QS.dir.im_r_sme = imFolder_rsme ;
QS.dir.im_r_sme_stack = imFolder_rsme_stack ;
QS.dir.cellProbabilities = fullfile(QS.dir.im_r_sme, 'cellProbabilities') ;
QS.dir.cellID = fullfile(QS.dir.im_r_sme, 'cellID') ;
QS.dir.lobe = lobeDir ;
QS.dir.foldHoopIm = foldHoopImDir ;
QS.dir.curvatures = KHSmDir ;
QS.dir.meanCurvature = HSmDir ;
QS.dir.gaussCurvature = KSmDir ;
QS.fullFileBase.cutMesh = ...
    fullfile(QS.dir.cutMesh, [QS.fileBase.name, '_cutMesh.mat']) ;
QS.fullFileBase.phi0fit = ...
    fullfile(QS.dir.spcutMesh, 'phi0s_%06d_%02d.png') ; 
QS.fullFileBase.clineDVhoop = ...
    fullfile(QS.dir.clineDVhoop,...
    'centerline_from_DVhoops_%06d.mat');

%% filenames for lobe dynamics
QS.fileName.fold = fullfile(lobeDir, ...
    ['fold_locations_sphi' uvexten '_avgpts.mat']) ;
QS.fileName.lobeDynamics = ...
    fullfile(lobeDir, ['lobe_dynamics' uvexten '.mat']) ;
QS.fileName.writhe = fullfile(QS.dir.writhe, ...
    ['writhe_sphi' uvexten '_avgpts.mat']) ;

%%  spcutMesh and pullbacks
QS.fullFileBase.spcutMesh = ...
    fullfile(sphiDir, 'mesh_apical_stab_%06d_spcutMesh.mat') ;
QS.fileBase.spcutMesh = 'mesh_apical_stab_%06d_spcutMesh' ;
QS.fullFileBase.spcutMeshSm = ...
    fullfile(sphiSmDir, '%06d_spcutMeshSm.mat') ;
QS.fileBase.spcutMeshSm = '%06d_spcutMeshSm' ;
QS.fullFileBase.spcutMeshSmRS = ...
    fullfile(sphiSmRSDir, '%06d_spcutMeshSmRS.mat') ;
QS.fileBase.spcutMeshSmRS = '%06d_spcutMeshSmRS' ;
QS.fullFileBase.spcutMeshSmRSC = ...
    fullfile(sphiSmRSCDir, '%06d_spcMSmRSC.mat') ;
QS.fullFileBase.spcutMeshSmRSCPLY = ...
    fullfile(sphiSmRSCDir, '%06d_spcMSmRSC.ply') ;
QS.fileBase.spcutMeshSmRSC = '%06d_spcMSmRSC' ;

%% double resolution spcutMeshSm's
QS.fullFileBase.spcutMeshSm2x = ...
    fullfile(sphiSmDir2x, '%06d_spcutMeshSm2x.mat') ;
QS.fileBase.spcutMeshSm2x = '%06d_spcutMeshSm2x' ;
QS.fullFileBase.spcutMeshSmRS2x = ...
    fullfile(sphiSmRSDir2x, '%06d_spcutMeshSmRS2x.mat') ;
QS.fileBase.spcutMeshSmRS2x = '%06d_spcutMeshSmRS2x' ;
QS.fullFileBase.spcutMeshSmRSC2x = ...
    fullfile(sphiSmRSCDir2x, '%06d_spcMSmRSC2x.mat') ;
QS.fullFileBase.spcutMeshSmRSCPLY2x = ...
    fullfile(sphiSmRSCDir2x, '%06d_spcMSmRSC2x.ply') ;
QS.fileBase.spcutMeshSmRSC2x = '%06d_spcMSmRSC2x' ;

QS.fileBase.im_uv = [QS.fileBase.name, '_pbuv.tif'] ;
QS.fullFileBase.im_uv = fullfile(QS.dir.im_uv, QS.fileBase.im_uv) ;
QS.fileBase.im_r = [QS.fileBase.name, '_pbr.tif'] ;
QS.fullFileBase.im_r = fullfile(QS.dir.im_r, QS.fileBase.im_r) ;
QS.fileBase.im_re = [QS.fileBase.name, '_pbre.tif'] ;
QS.fullFileBase.im_re = fullfile(QS.dir.im_re, QS.fileBase.im_re) ;
QS.fileBase.im_sp = [QS.fileBase.name, '_pbsp.tif'] ;
QS.fullFileBase.im_sp = fullfile(QS.dir.im_sp, QS.fileBase.im_sp);
QS.fileBase.im_spe = [QS.fileBase.name, '_pbspe.tif'] ;
QS.fullFileBase.im_spe = fullfile(QS.dir.im_spe, QS.fileBase.im_spe);
QS.fileBase.im_speLUT = [QS.fileBase.name, '_pbspe_LUT.tif'] ;
QS.fullFileBase.im_speLUT = fullfile(QS.dir.im_spe, QS.fileBase.im_speLUT);
QS.fileBase.im_up = [QS.fileBase.name, '_pbup.tif'] ;
QS.fullFileBase.im_up = fullfile(QS.dir.im_up, QS.fileBase.im_up) ;

 %% Smoothed pullbacks (pb)
QS.fileBase.im_sp_sm = [QS.fileBase.name, '_pbspsm.tif'] ;
QS.fileBase.im_sp_sme = [QS.fileBase.name, '_pbspsme.tif'] ;
QS.fullFileBase.im_sp_sm = ...
    fullfile(QS.dir.im_sp_sm, QS.fileBase.im_sp_sm);
QS.fullFileBase.im_sp_sme = ...
    fullfile(QS.dir.im_sp_sme, QS.fileBase.im_sp_sme);
QS.fileBase.im_r_sm = [QS.fileBase.name, '_pbrsm.tif'] ;
QS.fileBase.im_r_sme = [QS.fileBase.name, '_pbrsme.tif'] ;
QS.fullFileBase.im_r_sm = ...
    fullfile(QS.dir.im_r_sm, QS.fileBase.im_r_sm);
QS.fullFileBase.im_r_sme = ...
    fullfile(QS.dir.im_r_sme, QS.fileBase.im_r_sme);
QS.fullFileBase.cellProbabilities = ...
     fullfile(QS.dir.cellProbabilities, ...
     [QS.fileBase.name, '_pbrsme_Probabilities.h5']) ;
QS.fileBase.cellID = [QS.fileBase.name, '_pbrsme_cells'] ;
QS.fullFileBase.cellID = fullfile(QS.dir.cellID, ...
     [QS.fileBase.cellID '.mat']) ;

%% PIV
QS.dir.piv = fullfile(uvDir, 'piv') ;
QS.dir.piv3d = fullfile(QS.dir.piv, 'piv3d') ;
QS.dir.pivt2d = fullfile(QS.dir.piv, 'vt2d') ;
QS.dir.pivn2d = fullfile(QS.dir.piv, 'vn2d') ;
QS.dir.pivdilation = fullfile(QS.dir.piv, 'dilation') ;
QS.fileName.pivRaw = struct() ;
QS.fileName.pivRaw.raw = fullfile(QS.dir.piv, 'piv_results.mat') ;
QS.fullFileBase.piv3d = fullfile(QS.dir.piv3d, 'piv3d_%04d.mat') ;
% pathlines data
QS.dir.pathlines = struct() ;
QS.dir.pathlines.data = fullfile(QS.dir.piv, 'pathlines', 't0_%04d') ; 
QS.dir.pathlines.XY = fullfile(QS.dir.pathlines.data, 'images_XY') ;
QS.dir.pathlines.XYZ = fullfile(QS.dir.pathlines.data, 'images_XYZ') ;
QS.dir.pathlines.vXY = fullfile(QS.dir.pathlines.data, 'images_vXY') ;
QS.dir.pathlines.fXY = fullfile(QS.dir.pathlines.data, 'images_fXY') ;
QS.dir.pathlines.v3d = fullfile(QS.dir.pathlines.data, 'images_v3d') ;
QS.dir.pathlines.f3d = fullfile(QS.dir.pathlines.data, 'images_f3d') ;
QS.fileName.pathlines = struct() ;
pdir = QS.dir.pathlines.data ;
QS.fileName.pathlines.featureIDs = fullfile(pdir, 'featureIDs.mat') ;
QS.fileName.pathlines.XY = fullfile(pdir, 'piv_pathlines_XY.mat') ;
QS.fileName.pathlines.XYZ = fullfile(pdir, 'piv_pathlines_XYZ.mat') ;
QS.fileName.pathlines.vXY = fullfile(pdir, 'piv_pathlines_vXY.mat') ;
QS.fileName.pathlines.v3d = fullfile(pdir, 'piv_pathlines_v3d.mat') ;
QS.fileName.pathlines.fXY = fullfile(pdir, 'piv_pathlines_fXY.mat') ;
QS.fileName.pathlines.f3d = fullfile(pdir, 'piv_pathlines_f3d.mat') ;
% pathline velocities
QS.fileName.pathlines.velocities = struct() ;
QS.dir.pathlines.velocities = fullfile(pdir, 'velocities') ; 
pvdir = QS.dir.pathlines.velocities ;
QS.fileName.pathlines.velocities.v3d = fullfile(pvdir, 'vM.mat')  ;
QS.fileName.pathlines.velocities.vn  = fullfile(pvdir, 'vnM.mat') ;
QS.fileName.pathlines.velocities.vv  = fullfile(pvdir, 'vvM.mat') ;
QS.fileName.pathlines.velocities.vf  = fullfile(pvdir, 'vfM.mat') ;
QS.fileName.pathlines.velocities.v2d = fullfile(pvdir, 'v2dM.mat') ;
QS.fileName.pathlines.velocities.v2dum = fullfile(pvdir, 'v2dMum.mat') ;
QS.fileName.pathlines.velocities.v3dsm   = fullfile(pvdir, 'vsmM.mat') ;
QS.fileName.pathlines.velocities.vnsm    = fullfile(pvdir, 'vnsmM.mat') ;
QS.fileName.pathlines.velocities.vvsm    = fullfile(pvdir, 'vvsmM.mat') ;
QS.fileName.pathlines.velocities.vfsm    = fullfile(pvdir, 'vfsmM.mat') ;
QS.fileName.pathlines.velocities.v2dsm   = fullfile(pvdir, 'v2dsmM.mat') ;
QS.fileName.pathlines.velocities.v2dsmum = fullfile(pvdir, 'v2dsmMum.mat') ;

%% Velocities -- Lagrangian averaging
QS.dir.pivAvg = fullfile(QS.dir.piv, 'lagrangianAvg') ;
% QS.dir.pivAvgCurl = fullfile(QS.dir.pivAvg, 'curl') ;
% QS.dir.pivAvgDvg = fullfile(QS.dir.pivAvg, 'dvg') ;
QS.fileName.pivAvg = struct() ;
QS.fileName.pivAvg.v2dum = fullfile(QS.dir.pivAvg, 'v2dMum_avg.mat') ;
QS.fileName.pivAvg.v2d = fullfile(QS.dir.pivAvg, 'v2dM_avg.mat') ;
QS.fileName.pivAvg.vn  = fullfile(QS.dir.pivAvg, 'vnM_avg.mat') ;
QS.fileName.pivAvg.v3d = fullfile(QS.dir.pivAvg, 'vM_avg.mat')  ;
QS.fileName.pivAvg.vv  = fullfile(QS.dir.pivAvg, 'vvM_avg.mat') ;
QS.fileName.pivAvg.vf  = fullfile(QS.dir.pivAvg, 'vfM_avg.mat') ;
% Helmholtz-Hodge and DEC -- Lagrangian averaging
QS.dir.pivAvgDEC = struct() ;
QS.dir.pivAvgDEC.data   = fullfile(QS.dir.pivAvg, 'dec') ;
QS.dir.pivAvgDEC.div2d  = fullfile(QS.dir.pivAvg, 'dec_div2d') ;
QS.dir.pivAvgDEC.div3d  = fullfile(QS.dir.pivAvg, 'dec_div3d') ;
QS.dir.pivAvgDEC.div3dTexture = fullfile(QS.dir.pivAvg, 'dec_div3dTexture') ;
QS.dir.pivAvgDEC.rot2d  = fullfile(QS.dir.pivAvg, 'dec_rot2d') ;
QS.dir.pivAvgDEC.rot3d  = fullfile(QS.dir.pivAvg, 'dec_rot3d') ;
QS.dir.pivAvgDEC.rot3dTexture = fullfile(QS.dir.pivAvg, 'dec_rot3dTexture') ;
QS.dir.pivAvgDEC.harm2d = fullfile(QS.dir.pivAvg, 'dec_harm2d') ;
QS.dir.pivAvgDEC.harm3d = fullfile(QS.dir.pivAvg, 'dec_harm3d') ;
QS.fullFileBase.decAvg = fullfile(QS.dir.pivAvgDEC.data, ...
                              [QS.fileBase.name '_dec.mat'] ) ;
                          
%% Velocities -- no averaging
QS.dir.pivRaw = fullfile(QS.dir.piv, 'noAvg') ;
QS.fileName.pivRaw.v2dum = fullfile(QS.dir.pivRaw, 'v2dMum_simpletimeavg.mat') ;
QS.fileName.pivRaw.v2d = fullfile(QS.dir.pivRaw, 'v2dM_simpletimeavg.mat') ;
QS.fileName.pivRaw.vn  = fullfile(QS.dir.pivRaw, 'vnM_simpletimeavg.mat') ;
QS.fileName.pivRaw.v3d = fullfile(QS.dir.pivRaw, 'vM_simpletimeavg.mat') ;
QS.fileName.pivRaw.vv  = fullfile(QS.dir.pivRaw, 'vvM_simpletimeavg.mat') ;
QS.fileName.pivRaw.vf  = fullfile(QS.dir.pivRaw, 'vfM_simpletimeavg.mat');
% Helmholtz-Hodge and DEC -- simple/surface-Lagrangian averaging
QS.dir.pivRawDEC = struct() ;
QS.dir.pivRawDEC.data   = fullfile(QS.dir.pivRaw, 'dec') ;
QS.dir.pivRawDEC.div2d  = fullfile(QS.dir.pivRaw, 'dec_div2d') ;
QS.dir.pivRawDEC.div3d  = fullfile(QS.dir.pivRaw, 'dec_div3d') ;
QS.dir.pivRawDEC.div3dTexture = fullfile(QS.dir.pivRaw, 'dec_div3dTexture') ;
QS.dir.pivRawDEC.rot2d  = fullfile(QS.dir.pivRaw, 'dec_rot2d') ;
QS.dir.pivRawDEC.rot3d  = fullfile(QS.dir.pivRaw, 'dec_rot3d') ;
QS.dir.pivRawDEC.rot3dTexture = fullfile(QS.dir.pivRaw, 'dec_rot3dTexture') ;
QS.dir.pivRawDEC.harm2d = fullfile(QS.dir.pivRaw, 'dec_harm2d') ;
QS.dir.pivRawDEC.harm3d = fullfile(QS.dir.pivRaw, 'dec_harm3d') ;
QS.fullFileBase.decRaw = fullfile(QS.dir.pivRawDEC.data, ...
                              [QS.fileBase.name '_dec.mat'] ) ;

%% Velocities -- simple/surface-Lagrangian averaging
QS.dir.pivSimAvg = fullfile(QS.dir.piv, 'simpleAvg') ;
% QS.dir.pivSimAvgCurl = fullfile(QS.dir.pivSimAvg, 'rot') ;
% QS.dir.pivSimAvgDvg = fullfile(QS.dir.pivSimAvg, 'dvg') ;
QS.fileName.pivSimAvg = struct() ;
QS.fileName.pivSimAvg.v2dum = fullfile(QS.dir.pivSimAvg, 'v2dMum_simpletimeavg.mat') ;
QS.fileName.pivSimAvg.v2d = fullfile(QS.dir.pivSimAvg, 'v2dM_simpletimeavg.mat') ;
QS.fileName.pivSimAvg.vn  = fullfile(QS.dir.pivSimAvg, 'vnM_simpletimeavg.mat') ;
QS.fileName.pivSimAvg.v3d = fullfile(QS.dir.pivSimAvg, 'vM_simpletimeavg.mat') ;
QS.fileName.pivSimAvg.vv  = fullfile(QS.dir.pivSimAvg, 'vvM_simpletimeavg.mat') ;
QS.fileName.pivSimAvg.vf  = fullfile(QS.dir.pivSimAvg, 'vfM_simpletimeavg.mat');
% Helmholtz-Hodge and DEC -- simple/surface-Lagrangian averaging
QS.dir.pivSimAvgDEC = struct() ;
QS.dir.pivSimAvgDEC.data   = fullfile(QS.dir.pivSimAvg, 'dec') ;
QS.dir.pivSimAvgDEC.div2d  = fullfile(QS.dir.pivSimAvg, 'dec_div2d') ;
QS.dir.pivSimAvgDEC.div3d  = fullfile(QS.dir.pivSimAvg, 'dec_div3d') ;
QS.dir.pivSimAvgDEC.div3dTexture = fullfile(QS.dir.pivSimAvg, 'dec_div3dTexture') ;
QS.dir.pivSimAvgDEC.rot2d  = fullfile(QS.dir.pivSimAvg, 'dec_rot2d') ;
QS.dir.pivSimAvgDEC.rot3d  = fullfile(QS.dir.pivSimAvg, 'dec_rot3d') ;
QS.dir.pivSimAvgDEC.rot3dTexture = fullfile(QS.dir.pivSimAvg, 'dec_rot3dTexture') ;
QS.dir.pivSimAvgDEC.harm2d = fullfile(QS.dir.pivSimAvg, 'dec_harm2d') ;
QS.dir.pivSimAvgDEC.harm3d = fullfile(QS.dir.pivSimAvg, 'dec_harm3d') ;
QS.fullFileBase.decSimAvg = fullfile(QS.dir.pivSimAvgDEC.data, ...
                              [QS.fileBase.name '_dec.mat'] ) ;

%% Double resolution
QS.dir.piv3d2x = fullfile(QS.dir.piv, 'piv3dDoubleRes') ;
QS.dir.pivt2d2x = fullfile(QS.dir.piv, 'vt2dDoubleRes') ;
QS.dir.pivn2d2x = fullfile(QS.dir.piv, 'vn2dDoubleRes') ;
QS.dir.pivdilation2x = fullfile(QS.dir.piv, 'dilationDoubleRes') ;
QS.fullFileBase.piv3d2x = fullfile(QS.dir.piv3d2x, 'piv3dDoubleRes_%04d.mat') ;
%% 2x Velocities -- Lagrangian averaging
QS.dir.pivAvg2x = fullfile(QS.dir.piv, 'lagrangianAvgDoubleRes') ;
QS.fileName.pivAvg2x = struct() ;
QS.fileName.pivAvg2x.v2dum = fullfile(QS.dir.pivAvg2x, 'v2dMum_timeavg2x.mat') ;
QS.fileName.pivAvg2x.v2d = fullfile(QS.dir.pivAvg2x, 'v2dM_timeavg2x.mat') ;
QS.fileName.pivAvg2x.vn  = fullfile(QS.dir.pivAvg2x, 'vnM_timeavg2x.mat') ;
QS.fileName.pivAvg2x.v3d = fullfile(QS.dir.pivAvg2x, 'vM_timeavg2x.mat') ;
QS.fileName.pivAvg2x.vv  = fullfile(QS.dir.pivAvg2x, 'vvM_timeavg2x.mat') ;
QS.fileName.pivAvg2x.vf  = fullfile(QS.dir.pivAvg2x, 'vfM_timeavg2x.mat') ;
QS.dir.pivAvgDEC2x = struct() ;
QS.dir.pivAvgDEC2x.data   = fullfile(QS.dir.pivAvg2x, 'dec') ;
QS.dir.pivAvgDEC2x.div2d  = fullfile(QS.dir.pivAvg2x, 'dec_div2d') ;
QS.dir.pivAvgDEC2x.div3d  = fullfile(QS.dir.pivAvg2x, 'dec_div3d') ;
QS.dir.pivAvgDEC2x.div3dTexture = fullfile(QS.dir.pivAvg2x, 'dec_div3dTexture') ;
QS.dir.pivAvgDEC2x.rot2d  = fullfile(QS.dir.pivAvg2x, 'dec_rot2d') ;
QS.dir.pivAvgDEC2x.rot3d  = fullfile(QS.dir.pivAvg2x, 'dec_rot3d') ;
QS.dir.pivAvgDEC2x.rot3dTexture = fullfile(QS.dir.pivAvg2x, 'dec_rot3dTexture') ;
QS.dir.pivAvgDEC2x.harm2d = fullfile(QS.dir.pivAvg2x, 'dec_harm2d') ;
QS.dir.pivAvgDEC2x.harm3d = fullfile(QS.dir.pivAvg2x, 'dec_harm3d') ;
QS.fullFileBase.decAvg2x = fullfile(QS.dir.pivAvgDEC2x.data, ...
                              [QS.fileBase.name '_dec.mat'] ) ;
%% 2x Velocities -- simple/surface-Lagrangian averaging
QS.dir.pivSimAvg2x = fullfile(QS.dir.piv, 'simpleAvgDoubleRes') ;
QS.fileName.pivSimAvg2x = struct() ;
QS.fileName.pivSimAvg2x.v2dum = ...
    fullfile(QS.dir.pivSimAvg2x, 'v2dMum_simpletimeavg2x.mat') ;
QS.fileName.pivSimAvg2x.v2d = ...
    fullfile(QS.dir.pivSimAvg2x, 'v2dM_simpletimeavg2x.mat') ;
QS.fileName.pivSimAvg2x.vn = ...
    fullfile(QS.dir.pivSimAvg2x, 'vnM_simpletimeavg2x.mat') ;
QS.fileName.pivSimAvg2x.v3d = ...
    fullfile(QS.dir.pivSimAvg2x, 'vM_simpletimeavg2x.mat') ;
QS.fileName.pivSimAvg2x.vv = ...
    fullfile(QS.dir.pivSimAvg2x, 'vvM_simpletimeavg2x.mat') ;
QS.fileName.pivSimAvg2x.vf = ...
    fullfile(QS.dir.pivSimAvg2x, 'vfM_simpletimeavg2x.mat') ;
QS.dir.pivSimAvgDEC2x = struct() ;
QS.dir.pivSimAvgDEC2x.data   = fullfile(QS.dir.pivSimAvg2x, 'dec') ;
QS.dir.pivSimAvgDEC2x.div2d  = fullfile(QS.dir.pivSimAvg2x, 'dec_div2d') ;
QS.dir.pivSimAvgDEC2x.div3d  = fullfile(QS.dir.pivSimAvg2x, 'dec_div3d') ;
QS.dir.pivSimAvgDEC2x.div3dTexture = ...
    fullfile(QS.dir.pivSimAvg2x, 'dec_div3dTexture') ;
QS.dir.pivSimAvgDEC2x.rot2d  = fullfile(QS.dir.pivSimAvg2x, 'dec_rot2d') ;
QS.dir.pivSimAvgDEC2x.rot3d  = fullfile(QS.dir.pivSimAvg2x, 'dec_rot3d') ;
QS.dir.pivSimAvgDEC2x.rot3dTexture = ...
    fullfile(QS.dir.pivSimAvg2x, 'dec_rot3dTexture') ;
QS.dir.pivSimAvgDEC2x.harm2d = fullfile(QS.dir.pivSimAvg2x, 'dec_harm2d') ;
QS.dir.pivSimAvgDEC2x.harm3d = fullfile(QS.dir.pivSimAvg2x, 'dec_harm3d') ;
QS.fullFileBase.decSimAvg2x  = fullfile(QS.dir.pivSimAvgDEC2x.data, ...
                              [QS.fileBase.name '_dec.mat'] ) ;

%% Eulerian kinematics
QS.dir.eulerianKinematics = fullfile(QS.dir.uvCoord, 'eulerianKinematics') ;

%% Ensure directories
dirs2make = struct2cell(QS.dir) ;
for ii=1:length(dirs2make)
    dir2make = dirs2make{ii} ;
    if isa(dir2make, 'struct')
        dirfields = struct2cell(dir2make) ;
        for qq = 1:length(dirfields)
            dir2make = dirfields{qq} ;
            if isa(dir2make, 'struct')
                dirfieldsSub = struct2cell(dir2make) ;
                for pp = 1:length(dirfieldsSub)
                    dir2makeSub = dirfieldsSub{pp} ;
                    if ~exist(dir2makeSub, 'dir') && ~contains(dir2makeSub, '%04d') ...
                            && ~contains(dir2makeSub, '%0.3f')
                        mkdir(dir2makeSub)
                    end
                end
            else
                if ~exist(dir2make, 'dir') && ~contains(dir2make, '%04d') ...
                        && ~contains(dir2make, '%0.3f')
                    mkdir(dir2make)
                end
            end
        end
    else
        if ~exist(dir2make, 'dir')
            mkdir(dir2make)
        end
    end
end