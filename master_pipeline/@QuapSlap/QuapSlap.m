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
        dirBase
        fileName
        fileBase
        fullFileBase
        APDV = struct('resolution', [], ...
            'rot', [], ...
            'trans', [])
        flipy 
        nV 
        nU
        uvexten                 % naming extension with nU and nV like '_nU0100_nV0100'
        t0                      % reference time in the experiment
        normalShift
        features = struct('folds', [], ...  % #timepoints x #folds int, indices of nU sampling of folds
            'fold_onset', [], ...           % #folds x 1 float, timestamps (not indices) of fold onset
            'ssmax', [], ...                % #timepoints x 1 float, maximum length of the centerline at each timepoint
            'ssfold', [], ...               % #timepoints x #folds float, positional pathlength along centerline of folds
            'rssmax', [], ...               % #timepoints x 1 float, maximum proper length of the surface over time
            'rssfold', []) ;                % #timepoints x #folds float, positional proper length along surface of folds
        a_fixed
        phiMethod = '3dcurves'  % must be '3dcurves' or 'texture'
        endcapOptions
        plotting = struct('preview', false, ... % display intermediate results
            'save_ims', true, ...               % save images
            'xyzlim_um_buff', [], ...           % xyzlimits in um in RS coord sys with buffer
            'xyzlim_raw', [], ...               % xyzlimits in pixels
            'xyzlim_pix', [], ...               % xyzlimits in pixels RS
            'xyzlim_um', [], ...                % xyzlimits in um in RS coord sys
            'colors', [])                       % color cycle for QS
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
            'spcutMesh', [], ...
            'spcutMeshSm', []) 
        data = struct('adjustlow', 0, ...
            'adjusthigh', 0, ...
            'axisOrder', [1 2 3]) 
        currentData = struct('IV', [], ...
            'adjustlow', 0, ...
            'adjusthigh', 0 ) 
        currentVelocity = struct('piv3d', struct()) 
        velocitySimpleAverage = struct() 
        cleanCntrlines
        
    end
    
    % Some methods are hidden from public view. These are used internally
    % to the class.
    methods (Hidden)
        function QS = QuapSlap(xp, opts)
            QS.initializeQuapSlap(xp, opts)
        end
        initializeQuapSlap(QS, xp, opts)
        plotSPCutMeshSmSeriesUtility(QS, coordsys, options)
    end
    
    % Public methods, accessible from outside the class and reliant on 
    % properties of the class instance
    methods
        function setTime(QS, tt)
            if tt ~= QS.currentTime
                QS.currentMesh.cylinderMesh = [] ;
                QS.currentMesh.cylinderMeshClean = [] ;
                QS.currentMesh.cutMesh = [] ;
                QS.currentMesh.cutPath = [] ;
                QS.currentMesh.spcutMesh = [] ;
                QS.currentMesh.cutMesh = [] ;
                QS.currentMesh.spcutMesh = [] ;
                QS.currentMesh.spcutMeshSm = [] ;
                QS.currentData.IV = [] ;
                QS.currentData.adjustlow = 0 ;
                QS.currentData.adjusthigh = 0 ;
                QS.currentVelocity.piv3d = struct() ;
            end
            QS.currentTime = tt ;
            QS.xp.setTime(tt) ;
        end
        
        function t0 = t0set(QS, t0)
            % t0set(QS, t0) Set time offset to 1st fold onset or manually 
            if nargin < 2
                try
                    % Note that fold_onset is in units of timepoints, not 
                    % indices into timepoints
                    load(QS.fileName.fold, 'fold_onset') ;
                    QS.t0 = min(fold_onset) ;
                catch
                    error('No folding times saved to disk')
                end
            else
                QS.t0 = t0 ;
            end
            t0 = QS.t0 ;
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
        
        function [xyzlim_raw, xyzlim_pix, xyzlim_um, xyzlim_um_buff] = ...
                getXYZLims(QS)
            % Grab each xyzlim from self, otherwise load from disk
            % full resolution pix
            if ~isempty(QS.plotting.xyzlim_raw)
                xyzlim_raw = QS.plotting.xyzlim_raw ;
            else
                xyzlim_raw = dlmread(QS.fileName.xyzlim_raw, ',', 1, 0) ; 
                QS.plotting.xyzlim_raw = xyzlim_raw ;
            end
            % rotated scaled in full resolution pix
            if ~isempty(QS.plotting.xyzlim_pix)
                xyzlim_pix = QS.plotting.xyzlim_pix ;
            else
                xyzlim_pix = dlmread(QS.fileName.xyzlim_pix, ',', 1, 0) ; 
                QS.plotting.xyzlim_pix = xyzlim_pix ;
            end
            % rotated scaled APDV in micron
            if ~isempty(QS.plotting.xyzlim_um)
                xyzlim_um = QS.plotting.xyzlim_um ;
            else
                xyzlim_um = dlmread(QS.fileName.xyzlim_um, ',', 1, 0) ;
                QS.plotting.xyzlim_um = xyzlim_um ;
            end
            % rotated scaled APDV in micron, with padding
            if ~isempty(QS.plotting.xyzlim_um_buff)
                xyzlim_um_buff = QS.plotting.xyzlim_um_buff ;
            else
                xyzlim_um_buff = dlmread(QS.fileName.xyzlim_um_buff, ',', 1, 0) ;
                QS.plotting.xyzlim_um = xyzlim_um_buff ;
            end
        end
        
        function getFeatures(QS, varargin)
            %GETFEATURES(QS, varargin)
            %
            if nargin > 1
                for qq=1:length(varargin)
                    if isempty(eval(['QS.features.' varargin{qq}]))
                        disp(['Loading feature: ' varargin{qq}])
                        QS.loadFeatures(varargin{qq})
                    end
                end
            else
                QS.loadFeatures() ;
            end
        end
        function loadFeatures(QS, varargin)
            % Load all features stored in QS.features
            % 
            % Parameters
            % ----------
            % varargin : optional string list/cell
            %   which specific features to load
            %
            if nargin > 1
                if any(strcmp(varargin, {'folds', 'fold_onset', ...
                    'ssmax', 'ssfold', 'rssmax', 'rssfold'}))
                    
                    % Load all features relating to folds
                    disp('Loading folding features')
                    load(QS.fileName.fold, 'folds', 'fold_onset', ...
                        'ssmax', 'ssfold', 'rssmax', 'rssfold') ;
                    QS.features.folds = folds ;
                    QS.features.fold_onset = fold_onset ; 
                    QS.features.ssmax = ssmax ; 
                    QS.features.ssfold = ssfold ;
                    QS.features.rssmax = rssmax ;
                    QS.features.rssfold = rssfold ;
                else
                    error('Feature not recognized')
                end
            else
                % Load all features
                load(QS.fileName.fold, 'folds', 'fold_onset', ...
                    'ssmax', 'ssfold', 'rssmax', 'rssfold') ;
                QS.features.folds = folds ;
                QS.features.fold_onset = fold_onset ; 
                QS.features.ssmax = ssmax ; 
                QS.features.ssfold = ssfold ;
                QS.features.rssmax = rssmax ;
                QS.features.rssfold = rssfold ;
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
        
        % Get velocity
        function getCurrentVelocity(QS, varargin)
            if isempty(QS.currentTime)
                error('No currentTime set. Use QuapSlap.setTime()')
            end
            if isempty(varargin) 
                do_all = true ;
            else
                do_all = false ;
            end
            
            no_piv3d = isempty(fieldnames(QS.currentVelocity.piv3d)) ;
            if (do_all || contains(varargin, 'piv3d')) && no_piv3d
                % Load 3D data for piv results
                piv3dfn = QS.fullFileBase.piv3d ;
                load(sprintf(piv3dfn, QS.currentTime), 'piv3dstruct') ;
                QS.currentVelocity.piv3d = piv3dstruct ;
            end
        end
        
        % APDV methods
        function getAPDCOMs(QS, apdvCOMOptions)
            computeAPDCOMs(QS, apdvCOMOptions)
        end
        
        function ars = xyz2APDV(QS, a)
            %ars = xyz2APDV(QS, a)
            %   Transform 3d coords from XYZ data space to APDV coord sys
            [rot, trans] = QS.getRotTrans() ;
            ars = ((rot * a')' + trans) * QS.APDV.resolution ;
            if QS.flipy
                ars(:, 2) = - ars(:, 2) ;
            end
        end
        
        function dars = dx2APDV(QS, da)
            %dars = dx2APDV(QS, da)
            %   Transform 3d difference vector from XYZ data space to APDV 
            %   coord sys
            [rot, ~] = QS.getRotTrans() ;
            dars = ((rot * da')') * QS.APDV.resolution ;
            if QS.flipy
                dars(:, 2) = - dars(:, 2) ;
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
        
        % Surface Area and Volume over time
        measureSurfaceAreaVolume(QS, options)
        
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
        sliceMeshEndcaps(QS, endcapOpts, methodOpts)
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
        function loadCurrentSPCutMeshSm(QS)
            spcutMeshfn = sprintf(QS.fullFileBase.spcutMeshSm, QS.currentTime) ;
            tmp = load(spcutMeshfn, 'spcutMeshSm') ;
            QS.currentMesh.spcutMeshSm = tmp.spcutMeshSm ;
        end
        
        % Pullbacks
        generateCurrentPullbacks(QS, cutMesh, spcutMesh, spcutMeshSm, pbOptions)
        function doubleCoverPullbackImages(QS, options)
            % options : struct with fields
            %   coordsys : ('sp', 'uv', 'up')
            %       coordinate system to make double cover 
            %   overwrite : bool, default=false
            %       whether to overwrite current images on disk
            %   histeq : bool, default=true
            %       perform histogram equilization during pullback
            %       extension
            %   ntiles : int, default=50 
            %       The number of bins in each dimension for histogram equilization for a
            %       square original image. That is, the extended image will have (a_fixed *
            %       ntiles, 2 * ntiles) bins in (x,y).
            %   a_fixed : float, default=QS.a_fixed
            %       The aspect ratio of the pullback image: Lx / Ly       
            
            if nargin > 1
                % unpack options
                if isfield(options, 'coordsys')
                    coordsys = options.coordsys ;
                    options = rmfield(options, 'coordsys') ;
                    if strcmp(coordsys, 'sp')
                        imDir = QS.dir.im_sp ;
                        imDir_e = QS.dir.im_spe ;
                        fn0 = QS.fileBase.im_sp ;
                        ofn = QS.fileBase.im_spe ;
                    elseif strcmp(coordsys, 'spsm') || strcmp(coordsys, 'sp_sm')
                        imDir = QS.dir.im_sp_sm ;
                        imDir_e = QS.dir.im_sp_sme ;
                        fn0 = QS.fileBase.im_sp_sm ;
                        ofn = QS.fileBase.im_sp_sme ;
                    elseif strcmp(coordsys, 'rsm') || strcmp(coordsys, 'r_sm')
                        imDir = QS.dir.im_r_sm ;
                        imDir_e = QS.dir.im_r_sme ;
                        fn0 = QS.fileBase.im_r_sm ;
                        ofn = QS.fileBase.im_r_sme ;
                    elseif strcmp(coordsys, 'uv')
                        imDir = QS.dir.im_uv ;
                        imDir_e = QS.dir.im_uve ;
                        fn0 = QS.fileBase.im_uv ;
                        ofn = QS.fileBase.im_uv_sme ;
                    elseif strcmp(coordsys, 'up')
                        imDir = QS.dir.im_up ;
                        imDir_e = QS.dir.im_upe ;
                        fn0 = QS.fileBase.im_up ;
                        ofn = QS.fileBase.im_up_sme ;
                    end
                else
                    % Default value of coordsys = 'sp' ;
                    imDir = QS.dir.im_sp ;
                    imDir_e = QS.dir.im_spe ;
                    fn0 = QS.fileBase.im_sp ;
                    ofn = QS.fileBase.im_spe ;
                end
                
                % pack options if missing fields
                if ~isfield(options, 'histeq')
                    % equalize the histogram in patches of the image
                    options.histeq = true ;
                end
                if ~isfield(options, 'a_fixed')
                    % Assign the aspect ratio for histogram equilization
                    options.a_fixed = QS.a_fixed ;
                end
                if ~isfield(options, 'ntiles')
                    % Number of tiles per circumference and per unit ap
                    % length, so that if aspect ratio is two, there will be
                    % 2*ntiles samplings for histogram equilization along
                    % the ap axis
                    options.ntiles = 50 ;
                end
                if ~isfield(options, 'overwrite')
                    options.overwrite = false ;
                end
            else
                % Default options
                % Default value of coordsys = 'sp' ;
                options = struct() ;
                imDir = QS.dir.im_sp ;
                imDir_e = QS.dir.im_spe ;
                fn0 = QS.fileBase.im_sp ;
                ofn = QS.fileBase.im_spe ;
                options.histeq = true ;
                options.a_fixed = QS.a_fixed ;
                options.ntiles = ntiles ;
            end
            options.outFnBase = ofn ;
            extendImages(imDir, imDir_e, fn0, QS.xp.fileMeta.timePoints, options)
            disp(['done ensuring extended tiffs for ' imDir ' in ' imDir_e])
        end
        
        % measure writhe
        measureWrithe(QS, options)
        
        % folds & lobes
        identifyFolds(QS, options)
        [lengths, areas, volumes] = measureLobeDynamics(QS, options)
        plotLobes(QS, options) 
        function plotConstrictionDynamics(QS, overwrite)
            % Plot the location of the constrictions over time along with
            % centerlines over time
            QS.getXYZLims() ;
            QS.getRotTrans() ;
            QS.t0set() ;
            QS.loadFeatures() ;
            % Plot motion of avgpts at folds in yz plane over time
            aux_plot_avgptcline_lobes(QS.features.folds, ...
                QS.features.fold_onset, QS.dir.lobe, ...
                QS.uvexten, QS.plotting.save_ims, ...
                overwrite, QS.xp.fileMeta.timePoints - QS.t0,...
                QS.xp.fileMeta.timePoints, ...
                QS.fullFileBase.spcutMesh, QS.fullFileBase.clineDVhoop)
            
            % Plot motion of DVhoop at folds in yz plane over time
            aux_plot_constriction_DVhoops(QS.features.folds, ...
                QS.features.fold_onset, QS.dir.foldHoopIm,...
                QS.uvexten, QS.plotting.save_ims, ...
                overwrite, QS.xp.fileMeta.timePoints - QS.t0,...
                QS.xp.fileMeta.timePoints, QS.fullFileBase.spcutMesh, ...
                QS.fullFileBase.alignedMesh, ...
                QS.normalShift, QS.APDV.rot, QS.APDV.trans, QS.APDV.resolution, ...
                QS.plotting.colors, QS.plotting.xyzlim_um, QS.flipy)
        end
        
        % Smooth meshes in time
        [v3dsmM, nsmM] = smoothDynamicSPhiMeshes(QS, options) ;
        function plotSPCutMeshSm(QS, options) 
            plotSPCutMeshSmSeriesUtility(QS, 'spcutMeshSm', options)
        end
        function plotSPCutMeshSmRS(QS, options) 
            plotSPCutMeshSmSeriesUtility(QS, 'spcutMeshSmRS', options)
        end
        function plotSPCutMeshSmRSC(QS, options) 
            plotSPCutMeshSmSeriesUtility(QS, 'spcutMeshSmRSC', options)
        end
        function [v3dsmM, nsmM] = loadSPCutMeshSm(QS) 
            timePoints = QS.xp.fileMeta.timePoints ;
            v3dsmM = zeros(length(timePoints), QS.nU*QS.nV, 3);
            nsmM = zeros(length(timePoints), QS.nU*QS.nV, 3) ;
            % Load each mesh into v3dsmM and nsmM    
            for qq = 1:length(timePoints)
                load(sprintf(QS.fullFileBase.spcutMeshSm, ...
                    timePoints(qq)), 'spcutMeshSm') ;
                v3dsmM(qq, :, :) = spcutMeshSm.v ;
                nsmM(qq, :, :) = spcutMeshSm.vn ;
            end
        end
        
        % Mean & Gaussian curvature videos
        measureCurvatures(QS, options)
        
        % density of cells -- nuclei or membrane based
        measureCellDensity(QS, nuclei_or_membrane, options)
        function loadCurrentCellDensity(QS)
            if QS.currentData
                disp('Loading from self')
            else
                disp('Loading from disk')
            end
        end
        plotCellDensity(QS, options)
        plotCellDensityKymograph(QS, options)
        
        % spcutMeshSmStack
        generateSPCutMeshSmStack(QS, spcutMeshSmStackOptions)
        measureThickness(QS, thicknessOptions)
        phi0_fit = fitPhiOffsetsViaTexture(QS, uspace_ds_umax, vspace,...
            phi0_init, phi0TextureOpts)
       
        % spcutMeshSm coordinate system demo
        coordinateSystemDemo(QS)
        
        % flow measurements
        measurePIV3D(QS, options)
        timeAverageVelocitiesSimple(QS, options)
        plotTimeAvgVelSimple(QS, options)
        
        % compressible/incompressible flow on evolving surface
        [cumerr, HHs, divvs, velns] = measureCompressibility(QS, lambda, lambda_err)
        
    end
    
    methods (Static)
        function uv = XY2uv(im, XY, doubleCovered, umax, vmax)
            %XY2uv(im, XY, doubleCovered, umax, vmax)
            %   Map pixel positions (1, sizeImX) and (1, sizeImY) to
            %   (0, umax) and (0, vmax) of pullback space if singleCover,
            %   or y coords are mapped to (-0.5, 1.5)*vmax if doubleCover
            % 
            % NOTE THAT INVERSE MAP IS
            % x--> (xy(:, 1) * (Xsz-1)) / (1*umax) + 1 , ...
            % y--> (xy(:, 2) * (Ysz-1)) / (2*vmax) + 1 + (Ysz-1)*0.25 ;
            %
            % Parameters
            % ----------
            % im : NxM numeric array
            %   image in which pixel coordinates are defined
            % XY : Qx2 numeric array
            %   pixel coordinates to convert to pullback space
            % doubleCovered : bool
            %   the image is a double cover of the pullback (extended/tiled
            %   so that the "top" half repeats below the bottom and the
            %   "bottom" half repeats above the top
            % umax : float
            %   maximum extent of the pullback x coordinate
            % vmax : float
            %   maximum extent of the pullback y coordinate (BEFORE double 
            %   covering/tiling)
            
            if nargin < 3
                doubleCovered = true ;
            end
            if nargin < 4
                umax = 1.0 ;
            end
            if nargin < 5
                vmax = 1.0 ;
            end
            % size of extended image
            esize = size(im') ;
            % map extended image size to (0, 1), (-0.5, 1.5) if double
            % covered. 
            % subtract 1 since pixel positions range from (1, sizeIm)
            xesz = double(esize(1) - 1) ;
            yesz = double(esize(2) - 1) ;
            % map from pixel y to network y (sphi)
            uv = zeros(size(XY)) ;
            % convert x axis
            uv(:, 1) = umax * (XY(:, 1) - 1) / xesz ;
            % convert y axis
            % subtract 1 since pixel positions range from (1, sizeIm)
            if doubleCovered
                uv(:, 2) = vmax * 2.0 * (XY(:, 2) - 1) / yesz - 0.5 ;
            else
                uv(:, 2) = vmax * (XY(:, 2) - 1) / yesz ;
            end
        end
        
        function XY = uv2XY(im, uv, doubleCovered, umax, vmax) 
            % XY = uv2XY(im, uv, doubleCovered, umax, vmax) 
            %   Map from pullback uv u=(0,1), v=(0,1) to pixel XY
            % x--> (xy(:, 1) * (size(im, 2)-1)) / (1*umax) + 1 , ...
            % y--> (xy(:, 2) * (size(im, 1)-1)) / (2*vmax) + 0.75 + (size(im,1)-1)*0.25
            %
            if nargin < 3
                doubleCovered = true ;
            end
            if nargin < 4
                umax = 1.0 ;
            end
            if nargin < 5
                vmax = 1.0 ;
            end
            
            Xsz = size(im, 2) ;
            Ysz = size(im, 1) ;
            XY = 0*uv ;
            XY(:, 1) = uv(:, 1) * (Xsz-1) / (1*umax) + 1 ;
            if doubleCovered
                % image is double cover of physical object (periodic
                % cylinder)
                XY(:, 2) = uv(:, 2) * (Ysz-1) / (2*vmax) + 1 + (Ysz-1)*0.25 ;
            else
                % singleCover image of physical cylindrical object
                XY(:, 2) = uv(:, 2) * (Ysz-1) / (1*vmax) + 1  ;
            end        
        end
        
        function uv2pix_old(im, aspect)
            % map from network xy to pixel xy
            % Note that we need to flip the image (Yscale - stuff) since saved ims had
            % natural ydirection.
            % Assume here a coord system xy ranging from (0, xscale) and (0, 1) 
            % maps to a coord system XY ranging from (0, Yscale * 0.5) and (0, Yscale)
            x2Xpix = @(x, Yscale, xscale) (Yscale * 0.5) * aspect * x / xscale ;
            % y2Ypix = @(y, h, Yscale) Yscale - (Yscale*0.5)*(y+0.5) + h ;
            y2Ypix = @(y, h, Yscale) (Yscale*0.5)*(y+0.5) + h ;
            
            dx2dX = @ (y, Yscale, xscale) (Yscale * 0.5) * aspect * x / xscale ;
            dy2dY = @ (y, Yscale) (Yscale*0.5)*y ;
            
            % Now map the coornates
        end
    end
    
end