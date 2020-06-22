function initializeQuapSlap(QS, xp, opts)
%
%
%
%
% NPMitchell 2020

%% PROPERTIES
QS.xp = xp ;
QS.flipy = opts.flipy ;
meshDir = opts.meshDir ;
QS.timeinterval = opts.timeinterval ;
QS.timeunits = opts.timeunits ;
QS.fileBase.fn = xp.fileMeta.filenameFormat ;
QS.nV = opts.nV ;
QS.nU = opts.nU ;

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
QS.plotting.markers = {'o', 's', '^', 'v'} ;

% Directories of the QS object
QS.dir.dataDir = xp.fileMeta.dataDir ;
QS.dir.mesh = meshDir ;
QS.dir.alignedMesh = fullfile(meshDir, 'aligned_meshes') ;
QS.dir.cntrline = fullfile(meshDir, 'centerline') ;
QS.dir.cylinderMesh = fullfile(meshDir, 'cylinder_meshes') ;
QS.dir.cutMesh = fullfile(meshDir, 'cutMesh') ;
QS.dir.cylinderMeshClean = fullfile(QS.dir.cylinderMesh, 'cleaned') ;

% Metric strain dirs
QS.dir.gstrain = fullfile(meshDir, 'metric_strain') ;
QS.dir.gstrainRate = fullfile(QS.dir.gstrain, 'rateMetric') ;
QS.dir.gstrainRateIm = fullfile(QS.dir.gstrainRate, 'images') ;
QS.dir.gstrainVel = fullfile(QS.dir.gstrain, 'velMetric') ;
QS.dir.gstrainVelIm = fullfile(QS.dir.gstrainVel, 'images') ;
QS.dir.gstrainMesh = fullfile(QS.dir.gstrain, 'meshMetric') ;
QS.dir.gstrainMeshIm = fullfile(QS.dir.gstrainMesh, 'images') ;
QS.fileBase.gstrainVel = 'gstrainVel_%06d.mat' ;
QS.fileBase.gstrainMesh = 'gstrainMesh_%06d.mat' ;
QS.fileBase.gstrainRate = 'gstrainRate_%06d.mat' ;
QS.fullFileBase.gstrainVel = fullfile(QS.dir.gstrainVel, ...
    QS.fileBase.gstrainVel) ;
QS.fullFileBase.gstrainMesh = fullfile(QS.dir.gstrainMesh, ...
    QS.fileBase.gstrainMesh) ; 
QS.fullFileBase.gstrainRate = fullfile(QS.dir.gstrainRate, ...
    QS.fileBase.gstrainRate) ; 
QS.dir.compressibility = fullfile(QS.dir.mesh, 'compressibility') ;
QS.dir.compressibility2d = ...
    fullfile(QS.dir.compressibility, 'images_2d') ;
QS.dir.compressibility3d = ...
    fullfile(QS.dir.compressibility, 'images_3d') ;
QS.fullFileBase.compressibility2d = ...
    fullfile(QS.dir.compressibility2d, 'compr_2d_%06d.png') ;
QS.fullFileBase.compressibility3d = ...
    fullfile(QS.dir.compressibility3d, 'compr_3d_%06d.png') ;

% shorten variable names for brevity
clineDir = QS.dir.cntrline ;

% fileBases
QS.fileBase.name = xp.fileMeta.filenameFormat(1:end-4) ;
QS.fileBase.mesh = ...
    [xp.detector.options.ofn_smoothply '%06d'] ;
QS.fileBase.alignedMesh = ...
    [QS.fileBase.mesh '_APDV_um'] ;
QS.fileBase.centerlineXYZ = ...
    [QS.fileBase.mesh '_centerline_exp1p0_res*.txt' ] ;
QS.fileBase.centerlineAPDV = ...
    [QS.fileBase.mesh '_centerline_scaled_exp1p0_res*.txt' ] ;
QS.fileBase.cylinderMesh = ...
    [QS.fileBase.mesh '_cylindercut.ply'] ;
QS.fileBase.apBoundary = 'ap_boundary_indices_%06d.mat';
QS.fileBase.cylinderKeep = 'cylinderMesh_keep_indx_%06.mat' ;
QS.fileName.apBoundaryDorsalPts = 'ap_boundary_dorsalpts.h5' ;

% Clean Cylinder Mesh
QS.fileName.aBoundaryDorsalPtsClean = ...
    fullfile(QS.dir.cylinderMeshClean, 'adIDx.h5') ;
QS.fileName.pBoundaryDorsalPtsClean = ...
    fullfile(QS.dir.cylinderMeshClean, 'pdIDx.h5') ;

% cutMesh
QS.fullFileBase.cutPath = fullfile(QS.dir.cutMesh, 'cutPaths_%06d.txt') ;

% fileNames
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
QS.fileName.apdv = ...
    fullfile(clineDir, 'apdv_coms_from_training.h5') ;
QS.fileName.startendPt = fullfile(clineDir, 'startendpt.h5') ;
QS.fileName.cleanCntrlines = ...
    fullfile(clineDir, 'centerlines_anomalies_fixed.mat') ;
QS.fileName.apBoundaryDorsalPts = ...
    fullfile(QS.dir.cylinderMesh, 'ap_boundary_dorsalpts.h5') ;
QS.fileName.endcapOptions = ...
    fullfile(QS.dir.cylinderMesh, 'endcapOptions.mat') ;
QS.fileName.apdBoundary = ...
    fullfile(QS.dir.cylinderMesh, 'ap_boundary_dorsalpts.h5') ;

% FileNamePatterns
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
QS.fullFileBase.apBoundaryDorsalPts = ...
    fullfile(QS.dir.cylinderMesh, QS.fileName.apBoundaryDorsalPts) ;
QS.fullFileBase.cylinderKeep = ...
    fullfile(QS.dir.cylinderMesh, QS.fileBase.cylinderKeep) ;
QS.fullFileBase.cylinderMeshClean = ...
    fullfile(QS.dir.cylinderMesh, 'cleaned',...
    [QS.fileBase.mesh '_cylindercut_clean.ply']) ;            

% Define cutMesh directories
% cutFolder = fullfile(meshDir, 'cutMesh') ;
% cutMeshBase = fullfile(cutFolder, [QS.fileBase.name, '_cutMesh.mat']) ;
imFolderBase = fullfile(meshDir, ['PullbackImages' shiftstr uvexten] ) ;
sphiDir = fullfile(meshDir, ['sphi_cutMesh' shiftstr uvexten]) ;
sphiSmDir = fullfile(sphiDir, 'smoothed') ;
sphiSmRSDir = fullfile(sphiDir, 'smoothed_rs') ;
% sphiSmRSImDir = fullfile(sphiSmRSDir, 'images') ;
% sphiSmRSPhiImDir = fullfile(sphiSmRSImDir, 'phicolor') ;
sphiSmRSCDir = fullfile(sphiDir, 'smoothed_rs_closed') ;
imFolder_sp = [imFolderBase '_sphi'] ;
imFolder_spe = fullfile(imFolder_sp, 'extended') ;
imFolder_up = [imFolderBase '_uphi'] ;
imFolder_upe = fullfile(imFolder_up, 'extended') ;
% time-averaged meshes
imFolder_spsm = fullfile(imFolder_sp, 'smoothed') ;
imFolder_spsme = fullfile(imFolder_sp, 'smoothed_extended') ;  % raw LUT, no histeq
imFolder_spsme2 = fullfile(imFolder_sp, 'smoothed_extended_LUT') ;  % with histeq?
imFolder_rsm = fullfile([imFolderBase, '_sphi_relaxed'], 'smoothed');
imFolder_rsme = fullfile([imFolderBase, '_sphi_relaxed'], 'smoothed_extended') ;
imFolder_rsme_stack = fullfile([imFolderBase, '_sphi_relaxed'], 'smoothed_extended_stack') ;  % with histeq?

% Lobe/fold identification paths
lobeDir = fullfile(meshDir, 'lobes') ;
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
QS.dir.clineDVhoop = ...
    fullfile(QS.dir.cntrline, ...
    ['centerline_from_DVhoops' shiftstr uvexten]) ;
QS.dir.writhe =  fullfile(QS.dir.clineDVhoop, 'writhe') ;
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
QS.dir.im_sp_sme2 = imFolder_spsme2 ;
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
% filenames for lobe dynamics
QS.fileName.fold = fullfile(lobeDir, ...
    ['fold_locations_sphi' uvexten '_avgpts.mat']) ;
QS.fileName.lobeDynamics = ...
    fullfile(lobeDir, ['lobe_dynamics' uvexten '.mat']) ;

%  spcutMesh and pullbacks
QS.fullFileBase.spcutMesh = ...
    fullfile(sphiDir, 'mesh_apical_stab_%06d_spcutMesh.mat') ;
QS.fileBase.spcutMesh = 'mesh_apical_stab_%06d_spcutMesh' ;
QS.fullFileBase.spcutMeshSm = ...
    fullfile(sphiSmDir, '%06d_spcutMeshSm.mat') ;
QS.fileBase.spcutMeshSm = '%06d_spcMeshSm' ;
QS.fullFileBase.spcutMeshSmRS = ...
    fullfile(sphiSmRSDir, '%06d_spcutMeshSmRS.mat') ;
QS.fileBase.spcutMeshSmRS = '%06d_spcMeshSmRS' ;
QS.fullFileBase.spcutMeshSmRSC = ...
    fullfile(sphiSmRSCDir, '%06d_spcMSmRSC.mat') ;
QS.fullFileBase.spcutMeshSmRSCPLY = ...
    fullfile(sphiSmRSCDir, '%06d_spcMSmRSC.ply') ;
QS.fileBase.spcutMeshSmRSC = '%06d_spcMSmRSC' ;
QS.fileBase.im_uv = [QS.fileBase.name, '_pbuv.tif'] ;
QS.fullFileBase.im_uv = ...
    fullfile(QS.dir.im_uv, QS.fileBase.im_uv) ;
QS.fileBase.im_r = [QS.fileBase.name, '_pbr.tif'] ;
QS.fullFileBase.im_r = ...
    fullfile(QS.dir.im_r, QS.fileBase.im_r) ;
QS.fileBase.im_re = [QS.fileBase.name, '_pbre.tif'] ;
QS.fullFileBase.im_re =  ...
    fullfile(QS.dir.im_re, QS.fileBase.im_re) ;
QS.fileBase.im_sp = [QS.fileBase.name, '_pbsp.tif'] ;
QS.fullFileBase.im_sp = ...
    fullfile(QS.dir.im_sp, QS.fileBase.im_sp);
QS.fileBase.im_up = [QS.fileBase.name, '_pbup.tif'] ;
QS.fullFileBase.im_up = ...
     fullfile(QS.dir.im_up, QS.fileBase.im_up) ;
% Smoothed pullbacks (pb)
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

% PIV
QS.dir.piv = fullfile(meshDir, 'piv') ;
QS.dir.piv3d = fullfile(QS.dir.piv, 'piv3d') ;
QS.dir.pivt2d = fullfile(QS.dir.piv, 'vt2d') ;
QS.dir.pivn2d = fullfile(QS.dir.piv, 'vn2d') ;
QS.dir.pivSimAvg = fullfile(QS.dir.piv, 'simpleAvg') ;
QS.dir.pivSimAvgCurl = fullfile(QS.dir.pivSimAvg, 'curl') ;
QS.dir.pivSimAvgDvg = fullfile(QS.dir.pivSimAvg, 'dvg') ;
QS.fileName.pivRaw = fullfile(QS.dir.piv, 'piv_results.mat') ;
QS.fullFileBase.piv3d = fullfile(QS.dir.piv3d, 'piv3d_%04d.mat') ;
QS.fileName.pivSimAvg = struct() ;
QS.fileName.pivSimAvg.v2dMum = fullfile(QS.dir.pivSimAvg, 'v2dMum_simpletimeavg.mat') ;
QS.fileName.pivSimAvg.v2dM = fullfile(QS.dir.pivSimAvg, 'v2dM_simpletimeavg.mat') ;
QS.fileName.pivSimAvg.vnM = fullfile(QS.dir.pivSimAvg, 'vnM_simpletimeavg.mat') ;
QS.fileName.pivSimAvg.vM = fullfile(QS.dir.pivSimAvg, 'vM_simpletimeavg.mat') ;
QS.fileName.pivSimAvg.vvM = fullfile(QS.dir.pivSimAvg, 'vvM_simpletimeavg.mat') ;
QS.fileName.pivSimAvg.vfM = fullfile(QS.dir.pivSimAvg, 'vfM_simpletimeavg.mat') ;

% Ensure directories
dirs2make = struct2cell(QS.dir) ;
for ii=1:length(dirs2make)
    dir2make = dirs2make{ii} ;
    if ~exist(dir2make, 'dir')
        mkdir(dir2make)
    end
end