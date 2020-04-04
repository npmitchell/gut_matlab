classdef QuapSlap < handle
    % Quasi-Axisymmetric Pipeline for Surface Lagrangian Pullbacks class
    %
    % flipy         % APDC coord system is mirrored XZ wrt raw data
    % xyzlim        % mesh limits in full resolution pixels, in data space
	% xyzlim_um     % mesh limits in lab APDV frame in microns
    % resolution    % resolution of pixels in um
    % rot           % APDV rotation matrix
    % trans         % APDV translation 
    properties
        xp
        timeinterval
        timeunits
        dir
        fileName
        fileBase
        fullFileBase
        APDV = struct('resolution', [], ...
            'rot', [], ...
            'trans', [])
        flipy 
        nV 
        nU
        normalShift
        a_fixed
        endcapOptions
        plotting = struct('preview', false, ...
            'xyzlim', [], ...
            'xyzlim_raw', [], ...
            'xyzlim_um', [] )
        apdvCOM = struct('acom', [], ...
            'pcom', [], ... 
            'acom_sm', [], ...
            'pcom_sm', [], ... 
            'dcom', [], ... 
            'acom_rs', [], ... 
            'pcom_rs', [], ... 
            'dcom_rs', [])
        apdvCOMOptions
        currentTime
        currentMesh = struct('cylinderMesh', [], ...
            'cutMesh', [], ...
            'cutPath', [], ...
            'spcutMesh', []) 
        cleanCntrlines
        
    end
    
    methods
        function QS = QuapSlap(xp, opts)
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
            dvexten = sprintf('_nU%04d_nV%04d', QS.nU, QS.nV) ;
            
            % APDV coordinate system
            QS.APDV.resolution = min(xp.fileMeta.stackResolution) ;
            QS.APDV.rot = [] ;
            QS.APDV.trans = [] ;
            
            % dirs
            QS.dir.mesh = meshDir ;
            QS.dir.alignedMesh = fullfile(meshDir, 'aligned_meshes') ;
            QS.dir.cntrline = fullfile(meshDir, 'centerline') ;
            QS.dir.cylinderMesh = fullfile(meshDir, 'cylinder_meshes') ;
            QS.dir.cutMesh = fullfile(meshDir, 'cutMesh') ;
            QS.dir.cylinderMeshClean = fullfile(QS.dir.cylinderMesh, 'cleaned') ;
            
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
            QS.fileBase.apBoundaryDorsalPts = 'ap_boundary_dorsalpts.h5' ;
            
            % Clean Cylinder Mesh
            QS.fileName.aBoundaryDorsalPtsClean = ...
                fullfile(QS.dir.cylinderMeshClean, 'adIDx.h5') ;
            QS.fileName.pBoundaryDorsalPtsClean = ...
                fullfile(QS.dir.cylinderMeshClean, 'pdIDx.h5') ;
            
            % cutMesh
            QS.fullFileBase.cutPath = fullfile(QS.dir.cutMesh, 'cutPaths_%06d.txt') ;
            
            % fileNames
            QS.fileName.rot = fullfile(meshDir, 'rotation_APDV.txt') ;
            QS.fileName.trans = fullfile(meshDir, 'translation_APDV.txt') ;
            QS.fileName.xyzlim_raw = fullfile(meshDir, 'xyzlim.txt') ;
            QS.fileName.xyzlim = fullfile(meshDir, 'xyzlim.txt') ;
            QS.fileName.xyzlim_um = fullfile(meshDir, 'xyzlim_APDV_um.txt') ;
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
                fullfile(QS.dir.cylinderMesh, QS.fileBase.apBoundaryDorsalPts) ;
            QS.fullFileBase.cylinderKeep = ...
                fullfile(QS.dir.cylinderMesh, QS.fileBase.cylinderKeep) ;
            QS.fullFileBase.cylinderMeshClean = ...
                fullfile(QS.dir.cylinderMesh, 'cleaned',...
                [QS.fileBase.mesh 'cylindercut_clean.ply']) ;
            
            % Define cutMesh directories
            nshift = strrep(sprintf('%03d', QS.normalShift), '-', 'n') ;
            shiftstr = ['_' nshift 'step'] ;
            % cutFolder = fullfile(meshDir, 'cutMesh') ;
            % cutMeshBase = fullfile(cutFolder, [QS.fileBase.name, '_cutMesh.mat']) ;
            imFolder = fullfile(meshDir, ['PullbackImages' shiftstr] ) ;
            sphiDir = fullfile(meshDir, ['sphi_cutMesh' shiftstr dvexten]) ;
            spcutMeshBase = fullfile(sphiDir, 'mesh_apical_stab_%06d_spcutMesh.mat') ;
            sphiSmDir = fullfile(sphiDir, 'smoothed') ;
            sphiSmRSDir = fullfile(sphiDir, 'smoothed_rs') ;
            % sphiSmRSImDir = fullfile(sphiSmRSDir, 'images') ;
            % sphiSmRSPhiImDir = fullfile(sphiSmRSImDir, 'phicolor') ;
            sphiSmRSCDir = fullfile(sphiDir, 'smoothed_rs_closed') ;
            spcutMeshSmBase = fullfile(sphiSmDir, '%06d_spcutMeshSm.mat') ;
            spcutMeshSmRSBase = fullfile(sphiSmRSDir, '%06d_spcutMeshSmRS.mat') ;
            spcutMeshSmRSCBase = fullfile(sphiSmRSCDir, '%06d_spcMSmRSC.mat') ;
            imFolder_sp = [imFolder '_sphi' dvexten] ;
            imFolder_sp_e = [imFolder '_sphi' dvexten '_extended'] ;
            imFolder_up = [imFolder '_uphi' dvexten] ;
            imFolder_up_e = [imFolder '_uphi' dvexten '_extended'] ;
            % time-averaged meshes
            imFolder_spsm = fullfile(imFolder_sp, 'smoothed') ;
            imFolder_spsm_e = fullfile(imFolder_sp, '_extended_smoothed') ;
            imFolder_spsm_e2 = fullfile(imFolder_sp, '_extended_LUT_smoothed') ;  % raw LUT, no histeq
            imFolder_rsm = fullfile(imFolder, '_relaxed_smoothed') ;
            imFolder_rsm_e = fullfile(imFolder, '_relaxed_extended_smoothed') ;
            % Lobe identification paths
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
            QS.dir.spcutMeshSm = sphiSmRSDir ;
            QS.dir.spcutMeshSmRS = sphiSmRSCDir ;
            QS.dir.im = imFolder ;
            QS.dir.im_e = [imFolder '_extended'] ;
            QS.dir.im_r = [imFolder '_relaxed'] ;
            QS.dir.im_re = [imFolder '_relaxed_extended'] ;
            QS.dir.piv = fullfile(meshDir, 'piv') ;
            QS.dir.lobe = lobeDir ;
            QS.fullFileBase.cutMesh = ...
                fullfile(QS.dir.cutMesh, [QS.fileBase.name, '_cutMesh.mat']) ;
            QS.fullFileBase.spcutMesh = spcutMeshBase ;
            QS.fullFileBase.spcutMesh = spcutMeshBase ;
            QS.fullFileBase.spcutMeshSm = spcutMeshSmBase ;
            QS.fullFileBase.spcutMeshSmRS = spcutMeshSmRSBase ;
            QS.fullFileBase.spcutMeshSmRSC = spcutMeshSmRSCBase ;
            QS.fullFileBase.phi0fit = ...
                fullfile(QS.dir.spcutMesh, 'phi0s_%06d_%02d.png') ; 
            QS.dir.clineDVhoop = ...
                fullfile(QS.dir.cntrline, ...
                ['centerline_from_DVhoops' shiftstr dvexten]) ;
            
            % Ensure directories
            dirs2make = struct2cell(QS.dir) ;
            for ii=1:length(dirs2make)
                dir2make = dirs2make{ii} ;
                if ~exist(dir2make, 'dir')
                    mkdir(dir2make)
                end
            end
        end
        
        function setTime(QS, tt)
            if tt ~= QS.currentTime
                QS.currentMesh.cylinderMesh = [] ;
                QS.currentMesh.cutMesh = [] ;
                QS.currentMesh.cutPath = [] ;
                QS.currentMesh.spcutMesh = [] ;
                QS.currentMesh.cutMesh = [] ;
                QS.currentMesh.spcutMesh = [] ;
            end
            QS.currentTime = tt ;
            QS.xp.setTime(tt) ;
        end
        
        function [acom_sm, pcom_sm] = getAPCOMSm(QS) 
            acom_sm = h5read(QS.fileName.apdv, '/acom') ;
            pcom_sm = h5read(QS.fileName.apdv, '/pcom') ;
        end
        
        function [rot, trans] = getRotTrans(QS)
            % Load the translation to put anterior to origin
            if ~isempty(QS.APDV.trans)
                % return from self
                trans = QS.APDV.trans ;
            else
                % load from disk
                trans = importdata(QS.fileName.trans) ;
                QS.APDV.trans = trans ;
            end
            % Load the rotation from XYZ to APDV coordinates
            if ~isempty(QS.APDV.rot)
                rot = QS.APDV.rot ;
            else
                % Load the rotation matrix
                rot = importdata(QS.fileName.rot) ;
                QS.APDV.rot = rot ;
            end
        end
        
        function [xyzlim_raw, xyzlim, xyzlim_um] = getXYZLims(QS)
            % Grab each xyzlim from self, otherwise load from disk
            if ~isempty(QS.plotting.xyzlim_raw)
                xyzlim_raw = QS.plotting.xyzlim_raw ;
            else
                xyzlim_raw = dlmread(QS.fileName.xyzlim_raw, ',', 1, 0) ; 
            end
            if ~isempty(QS.plotting.xyzlim)
                xyzlim = QS.plotting.xyzlim ;
            else
                xyzlim = dlmread(QS.fileName.xyzlim, ',', 1, 0) ;
            end
            if ~isempty(QS.plotting.xyzlim_um)
                xyzlim_um = QS.plotting.xyzlim_um ;
            else
                xyzlim_um = dlmread(QS.fileName.xyzlim_um, ',', 1, 0) ;
            end
        end
        
        function getAPDCOMs(QS, apdvCOMOptions)
            computeAPDCOMs(QS, apdvCOMOptions)
        end
        
        function ars = xyz2APDV(QS, a)
            % transform 3d coords from XYZ data space to APDV coord sys
            [ro, tr] = QS.getRotTrans() ;
            ars = ((ro * a')' + tr) * QS.APDV.resolution ;
            if QS.flipy
                ars(:, 2) = - ars(:, 2) ;
            end
        end
        
        % APDV methods
        function setAPDVCOMOptions(QS, apdvCOMOpts)
            QS.apdvCOMOptions = apdvCOMOpts ;
        end        
        function apdvCOMOptions = loadAPDVCOMOptions(QS)
            load(QS.fileName.apdvCOMOptions, 'apdvCOMOptions')
            QS.apdvCOMOptions = apdvCOMOptions ;
        end        
        function apdvCOMOptions = saveAPDVCOMOptions(QS)
            apdvCOMOptions = QS.APDVCOMs.apdvCOMOptions ;
            save(QS.fileName.apdvCOMOptions, 'apdvCOMOptions')
        end
        
        % Endcap methods
        function setEndcapOptions(QS, endcapOpts)
            QS.endcapOptions = endcapOpts ;
        end        
        function loadEndcapOptions(QS)
            load(QS.fileName.endcapOptions, 'endcapOptions')
        end        
        function saveEndcapOptions(QS)
            endcapOptions = QS.endcapOptions ;
            save(QS.fileName.endcapOptions, 'endcapOptions')
        end
        
        generateCleanCntrlines(QS, idOptions)
        
        function getCleanCntrlines(QS)
            if isempty(QS.cleanCntrlines)
                try
                    tmp = load(QS.fileName.cleanCntrlines, 'cntrlines') ;
                    QS.cleanCntrlines = tmp.cntrlines ;
                    disp('Loaded clean centerlines from disk')
                catch
                    disp('No clean centerlines on disk, generating...')
                    QS.cleanCntrlines = generateCleanCntrlines(QS, idOptions) ;
                end
            end
        end
            
        function loadCurrentCylinderMesh(QS)
            cylmeshfn = ...
                sprintf( QS.fullFileBase.cylinderMesh, QS.currentTime ) ;
            QS.currentMesh.cylinderMesh = read_ply_mod( cylmeshfn );
        end
        
        generateCurrentCutMesh(QS)
        
        plotCutPath(QS, cutMesh, cutPath)
        
        function loadCurrentCutMesh(QS)
            cutMeshfn = sprintf(QS.fullFileBase.cutMesh, QS.currentTime) ;
            cutPfn = sprintf(QS.fullFileBase.cutPath, QS.currentTime) ;
            cutMesh = load(cutMeshfn, 'cutMesh') ;
            cutMesh.v = cutMesh.v + cutMesh.vn * QS.normalShift ;
            QS.currentMesh.cutMesh = cutMesh ;
            QS.currentMesh.cutPath = dlmread(cutPfn, ',', 1, 0) ;
        end
        
        generateCurrentSPCutMesh(QS, cutMesh, overwrite)
        
        function loadCurrentSPCutMesh(QS)
            spcutMeshfn = sprintf(QS.fullFileBase.spcutMesh, QS.currentTime) ;
            QS.currentMesh.spcutMesh = load(spcutMeshfn, 'cutMesh') ;
        end
    end
   
    methods (Static)
    end
end