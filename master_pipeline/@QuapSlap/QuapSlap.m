lassdef QuapSlap < handle
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
        axisOrderIV2Mesh
        a_fixed
        phiMethod = '3dcurves'  % must be '3dcurves' or 'texture'
        endcapOptions
        plotting = struct('preview', false, ...
            'save_ims', true, ...
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
            'cylinderMeshClean', [], ...
            'cutMesh', [], ...
            'cutPath', [], ...
            'spcutMesh', []) 
        data = struct('adjustlow', 0, ...
            'adjusthigh', 0, ...
            'axisOrder', [1 2 3]) 
        currentData = struct('IV', [], ...
            'adjustlow', 0, ...
            'adjusthigh', 0 ) 
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
            if isfield(opts, 'adjustlow')
                QS.data.adjustlow = opts.adjustlow ;
            end
            if isfield(opts, 'adjusthigh')
                QS.data.adjusthigh = opts.adjusthigh ;
            end
            if isfield(opts, 'axisOrder')
                QS.data.axisOrder = opts.axisOrder ;
            end
            dvexten = sprintf('_nU%04d_nV%04d', QS.nU, QS.nV) ;
            
            % APDV coordinate system
            QS.APDV.resolution = min(xp.fileMeta.stackResolution) ;
            QS.APDV.rot = [] ;
            QS.APDV.trans = [] ;
            
            % dirs
            QS.dir.dataDir = xp.fileMeta.dataDir ;
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
                [QS.fileBase.mesh '_cylindercut_clean.ply']) ;
            
            % Define cutMesh directories
            nshift = strrep(sprintf('%03d', QS.normalShift), '-', 'n') ;
            shiftstr = ['_' nshift 'step'] ;
            % cutFolder = fullfile(meshDir, 'cutMesh') ;
            % cutMeshBase = fullfile(cutFolder, [QS.fileBase.name, '_cutMesh.mat']) ;
            imFolder = fullfile(meshDir, ['PullbackImages' shiftstr] ) ;
            sphiDir = fullfile(meshDir, ['sphi_cutMesh' shiftstr dvexten]) ;
            sphiSmDir = fullfile(sphiDir, 'smoothed') ;
            sphiSmRSDir = fullfile(sphiDir, 'smoothed_rs') ;
            % sphiSmRSImDir = fullfile(sphiSmRSDir, 'images') ;
            % sphiSmRSPhiImDir = fullfile(sphiSmRSImDir, 'phicolor') ;
            sphiSmRSCDir = fullfile(sphiDir, 'smoothed_rs_closed') ;
            imFolder_sp = [imFolder '_sphi' dvexten] ;
            imFolder_sp_e = fullfile(imFolder_sp, '_extended') ;
            imFolder_up = [imFolder '_uphi' dvexten] ;
            imFolder_up_e = fullfile(imFolder_up, '_extended') ;
            % time-averaged meshes
            imFolder_spsm = fullfile(imFolder_sp, '_smoothed') ;
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
            QS.dir.clineDVhoop = ...
                fullfile(QS.dir.cntrline, ...
                ['centerline_from_DVhoops' shiftstr dvexten]) ;
            QS.dir.im = imFolder ;
            QS.dir.im_e = [imFolder '_extended'] ;
            QS.dir.im_r = [imFolder '_relaxed'] ;
            QS.dir.im_re = [imFolder '_relaxed_extended'] ;
            QS.dir.im_sp = imFolder_sp ;
            QS.dir.im_sp_e = imFolder_sp_e ;
            QS.dir.im_up = imFolder_up ;
            QS.dir.im_up = imFolder_up_e ;
            QS.dir.im_spsm = imFolder_spsm ;
            QS.dir.im_spsm_e = imFolder_spsm_e ;
            QS.dir.im_spsm_e2 = imFolder_spsm_e2 ;
            QS.dir.im_rsm = imFolder_rsm ;
            QS.dir.im_rsm_e = imFolder_rsm_e ;
            QS.dir.piv = fullfile(meshDir, 'piv') ;
            QS.dir.lobe = lobeDir ;
            QS.fullFileBase.cutMesh = ...
                fullfile(QS.dir.cutMesh, [QS.fileBase.name, '_cutMesh.mat']) ;
            QS.fullFileBase.phi0fit = ...
                fullfile(QS.dir.spcutMesh, 'phi0s_%06d_%02d.png') ; 
            QS.fullFileBase.clineDVhoop = ...
                fullfile(QS.dir.clineDVhoop,...
                'centerline_from_DVhoops_%06d.mat');
            
            %  spcutMesh and pullbacks
            QS.fullFileBase.spcutMesh = ...
                fullfile(sphiDir, 'mesh_apical_stab_%06d_spcutMesh.mat') ;
            QS.fullFileBase.spcutMeshSm = ...
                fullfile(sphiSmDir, '%06d_spcutMeshSm.mat') ;
            QS.fullFileBase.spcutMeshSmRS = ...
                fullfile(sphiSmRSDir, '%06d_spcutMeshSmRS.mat') ;
            QS.fullFileBase.spcutMeshSmRSC = ...
                fullfile(sphiSmRSCDir, '%06d_spcMSmRSC.mat') ;
            QS.fullFileBase.im = ...
                fullfile([QS.dir.im, '/', QS.fileBase.name, '_pb.tif']) ;
            QS.fullFileBase.im_r = ...
                fullfile([QS.dir.im_r, '/', QS.fileBase.name, '_pr.tif']) ;
            QS.fullFileBase.im_re =  ...
                fullfile([QS.dir.im_re, '/', QS.fileBase.name, '_pre.tif']) ;
            QS.fullFileBase.im_sp = ...
                fullfile([QS.dir.im_sp, '/', QS.fileBase.name, '_pbsp.tif']);
            QS.fullFileBase.im_up = ...
                 fullfile([QS.dir.im_up, '/', QS.fileBase.name, '_pbup.tif']) ;
            
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
                QS.currentMesh.cylinderMeshClean = [] ;
                QS.currentMesh.cutMesh = [] ;
                QS.currentMesh.cutPath = [] ;
                QS.currentMesh.spcutMesh = [] ;
                QS.currentMesh.cutMesh = [] ;
                QS.currentMesh.spcutMesh = [] ;
                QS.currentData.IV = [] ;
                QS.currentData.adjustlow = 0 ;
                QS.currentData.adjusthigh = 0 ;
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
        
        function data = loadBioFormats(QS, fullFileName)
            r = bfGetReader(fullFileName);
            r.setSeries(QS.xp.fileMeta.series-1);
            nChannelsUsed = numel(QS.xp.expMeta.channelsUsed);
            if QS.xp.fileMeta.swapZT == 0
                stackSize = [r.getSizeX(), r.getSizeY(), r.getSizeZ(), r.getSizeT()];
            else
                stackSize = [r.getSizeX(), r.getSizeY(), r.getSizeT(), r.getSizeZ()];
            end
            debugMsg(2, ['stack size (xyzt) ' num2str(stackSize) '\n']);

            xSize = stackSize(1);
            ySize = stackSize(2);
            zSize = stackSize(3);
            
            % number of channels
            nChannels = r.getSizeC();
            nTimePts = stackSize(4);
            
            data = zeros([ySize xSize zSize nChannelsUsed], 'uint16');
            for i = 1:r.getImageCount()

                ZCTidx = r.getZCTCoords(i-1) + 1;
                
                % in the fused embryo data coming out of the python script,
                % Z and T are swaped. In general this isn't the case, thus
                % introduce a file metaField swapZT
                if QS.xp.fileMeta.swapZT == 0
                    zidx = ZCTidx(1);
                    tidx = ZCTidx(3);
                else 
                    zidx = ZCTidx(3);
                    tidx = ZCTidx(1);
                end
                cidx = ZCTidx(2);

                % see above: if there is only one timepoint all the planes
                % should be read, if there are multiple timepoints, only
                % the correct time should be read
                if nTimePts == 1 || (nTimePts > 1 && this.currentTime == tidx-1)
                    
                    debugMsg(1,'.');
                    if rem(i,80) == 0
                        debugMsg(1,'\n');
                    end

                    dataCidx = find(QS.xp.expMeta.channelsUsed == cidx);
                    if ~isempty(dataCidx)
                        data(:,:, zidx, dataCidx) = bfGetPlane(r, i);
                    else
                        disp('skipping channel and z plane')
                    end
                end
            end
        end
        
        function getCurrentData(QS)
            if isempty(QS.currentTime)
                error('No currentTime set. Use QuapSlap.setTime()')
            end
            if isempty(QS.currentData.IV)
                % Load 3D data for coloring mesh pullback
                QS.xp.loadTime(QS.currentTime);
                QS.xp.rescaleStackToUnitAspect();
                IV = QS.xp.stack.image.apply() ;
                adjustlow = QS.data.adjustlow ;
                adjusthigh = QS.data.adjusthigh ;
                QS.currentData.IV = QS.adjustIV(IV, adjustlow, adjusthigh) ;
            end
        end
        
        function IV = adjustIV(QS, IV, adjustlow, adjusthigh)
            if nargin > 2 
                adjustlow = QS.data.adjustlow ;
                adjusthigh = QS.data.adjusthigh ;
            end
            if nargin < 2 
                if ~isempty(QS.currentData.IV) 
                    IV = QS.currentData.IV ;
                end
            end
            % custom image intensity adjustment
            if adjustlow == 0 && adjusthigh == 0
                disp('Using default limits for imadjustn')
                for ii = 1:length(IV)
                    IV{ii} = imadjustn(IV{ii});
                end
            else
                disp('Taking custom limits for imadjustn')
                for ii = 1:length(IV)
                    IVii = IV{ii} ;
                    vlo = double(prctile( IVii(:) , adjustlow )) / double(max(IVii(:))) ;
                    vhi = double(prctile( IVii(:) , adjusthigh)) / double(max(IVii(:))) ;
                    disp(['--> ' num2str(vlo) ', ' num2str(vhi)])
                    IV{ii} = imadjustn(IVii, [double(vlo); double(vhi)]) ;
                end
            end
            if nargout > 0
                disp('Attributing to self.currentData.IV')
                QS.currentData.IV = IV ;
                QS.currentData.adjustlow = adjustlow ;
                QS.currentData.adjustlow = adjusthigh ;
            else
                disp('WARNING: returning IV instead of attributing to self')
            end
        end
        
        % APDV methods
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
        
        % Centerlines & cylinderMesh
        extractCenterlineSeries(QS, cntrlineOpts)
        function setEndcapOptions(QS, endcapOpts)
            QS.endcapOptions = endcapOpts ;
        end        
        function loadEndcapOptions(QS)
            QS.endcapOptions = ...
                load(QS.fileName.endcapOptions, 'endcapOptions');
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
        function loadCurrentCylinderMeshlean(QS)
            cylmeshfn = ...
                sprintf( QS.fullFileBase.cylinderMesh, QS.currentTime ) ;
            QS.currentMesh.cylinderMesh = read_ply_mod( cylmeshfn );
        end
        function loadCurrentCylinderMeshClean(QS)
            cylmeshfn = ...
                sprintf( QS.fullFileBase.cylinderMeshClean, QS.currentTime ) ;
            disp(['Loading cylinderMeshClean ' cylmeshfn])
            QS.currentMesh.cylinderMeshClean = read_ply_mod( cylmeshfn );
        end
        
        % cutMesh
        generateCurrentCutMesh(QS)
        plotCutPath(QS, cutMesh, cutPath)
        function loadCurrentCutMesh(QS)
            if isempty(QS.currentTime)
                error('No currentTime set. Use QuapSlap.setTime()')
            end
            cutMeshfn = sprintf(QS.fullFileBase.cutMesh, QS.currentTime) ;
            cutPfn = sprintf(QS.fullFileBase.cutPath, QS.currentTime) ;
            tmp = load(cutMeshfn, 'cutMesh') ;
            tmp.cutMesh.v = tmp.cutMesh.v + tmp.cutMesh.vn * QS.normalShift ;
            QS.currentMesh.cutMesh = tmp.cutMesh ;
            QS.currentMesh.cutPath = dlmread(cutPfn, ',', 1, 0) ;
        end
        
        % spcutMesh
        generateCurrentSPCutMesh(QS, cutMesh, overwrite)
        function loadCurrentSPCutMesh(QS)
            spcutMeshfn = sprintf(QS.fullFileBase.spcutMesh, QS.currentTime) ;
            tmp = load(spcutMeshfn, 'spcutMesh') ;
            QS.currentMesh.spcutMesh = tmp.spcutMesh ;
        end
        
        % spcutMeshSmStack
        generateSPCutMeshSmStack(QS, spcutMeshSmStackOptions)
        measureThickness(QS, thicknessOptions)

    end
    methods (Static)
    end
end