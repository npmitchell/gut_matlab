classdef QuapSlap < handle
    % Quasi-Axisymmetric Pipeline for Surface Lagrangian Pullbacks class
    %
    % xyzlim        : 3x2 float, mesh limits in full resolution pixels, in data space
	% xyzlim_um     : 3x2 float, mesh limits in lab APDV frame in microns
    % resolution    : float, resolution of pixels in um
    % rot           : 3x3 float, APDV rotation matrix
    % trans         : 3x1 float, APDV translation 
    % a_fixed       : float, aspect ratio for fixed-width pullbacks
    % phiMethod     : str, '3dcurves' or 'texture'
    % flipy         : bool, APDV coord system is mirrored XZ wrt raw data
    % nV            : int, sampling number along circumferential axis
    % nU            : int, sampling number along longitudinal axis
    % uvexten       : str, naming extension with nU and nV like '_nU0100_nV0100'
    % t0            : int, reference timePoint in the experiment
    % features : struct with fields
    %   folds : #timepoints x #folds int, 
    %       indices of nU sampling of folds
    %   fold_onset : #folds x 1 float
    %       timestamps (not indices) of fold onset
    %   ssmax : #timepoints x 1 float
    %       maximum length of the centerline at each timepoint
    %   ssfold : #timepoints x #folds float
    %       positional pathlength along centerline of folds
    %   rssmax : #timepoints x 1 float
    %       maximum proper length of the surface over time
    %   rssfold : #timepoints x #folds float
    %       positional proper length along surface of folds
    % velocityAverage : struct with fields
    %   vsmM : (#timePoints-1) x (nX*nY) x 3 float array
    %       3d velocities at PIV evaluation coordinates in um/dt rs
    %   vfsmM : (#timePoints-1) x (2*nU*(nV-1)) x 3 float array
    %       3d velocities at face barycenters in um/dt rs
    %   vnsmM : (#timePoints-1) x (nX*nY) float array
    %       normal velocity at PIV evaluation coordinates in um/dt rs
    %   vvsmM : (#timePoints-1) x (nU*nV) x 3 float array
    %       3d velocities at (1x resolution) mesh vertices in um/min rs
    %   v2dsmM : (#timePoints-1) x (nX*nY) x 2 float array
    %       2d velocities at PIV evaluation coordinates in pixels/ min
    %   v2dsmMum : (#timePoints-1) x (nX*nY) x 2 float array
    %       2d velocities at PIV evaluation coordinates in scaled pix/min, but 
    %       proportional to um/min (scaled by dilation of map)
    properties
        xp
        timeInterval
        timeUnits
        spaceUnits
        dir
        dirBase
        fileName
        fileBase
        fullFileBase
        ssfactor                % subsampling factor for probabilities 
        APDV = struct('resolution', [], ...
            'rot', [], ...
            'trans', [])
        flipy                   % whether data is mirror image of lab frame coordinates
        nV                      % sampling number along circumferential axis
        nU                      % sampling number along longitudinal axis
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
            'axisOrder', [1 2 3]) % options for scaling and transposing image intensity data
        currentData = struct('IV', [], ...
            'adjustlow', 0, ...
            'adjusthigh', 0 )    % image intensity data in 3d and scaling
        currentVelocity = struct('piv3d', struct(), ...
            'piv3d2x', struct()) ;     
        piv = struct( ...
            'imCoords', 'sp_sme', ... % image coord system for measuring PIV / optical flow) ;
            'Lx', [], ...
            'Ly', [], ...
            'raw', struct()) ;  
        velocityAverage = struct('v3d', [], ... 
            'v2d', [], ...
            'v2dum', [], ...
            'vn', [], ...
            'vf', [], ...       
            'vv', []) ;          % velocity field after Lagrangian avg
        velocityAverage2x = struct('v3d', [], ...
            'v2d', [], ...
            'v2dum', [], ...
            'vn', [], ...
            'vf', [], ...       
            'vv', []) ;          % velocity field after Lagrangian avg
        velocitySimpleAverage = struct('v3d', [], ...
            'v2d', [], ...
            'v2dum', [], ...
            'vn', [], ...
            'vf', [], ...       
            'vv', []) ;          % velocity field after in-place (uv) avg
        velocitySimpleAverage2x = struct('v3d', [], ...
            'v2d', [], ...
            'v2dum', [], ...
            'vn', [], ...
            'vf', [], ...       
            'vv', []) ;          % velocity field after in-place (uv) avg
        cleanCntrlines
        pivPullback = 'sp_sme' ; % coordinate system used for velocimetry
        pathlines = struct('t0', [], ...   % timestamp (not an index) at which pathlines form regular grid in space
            'piv', [], ...                 % Lagrangian pathlines from piv coords
            'vertices', [], ...            % Lagrangian pathlines from mesh vertices
            'faces', []) ;                 % Lagrangian pathlines from mesh face barycenters
    end
    
    % Some methods are hidden from public view. These are used internally
    % to the class.
    methods (Hidden)
        function QS = QuapSlap(xp, opts)
            QS.initializeQuapSlap(xp, opts)
        end
        initializeQuapSlap(QS, xp, opts)
        plotSPCutMeshSmSeriesUtility(QS, coordsys, options)
        plotMetricKinematicsTimePoint(QS, tp, options)
        [XX, YY] = pullbackPathlines(QS, x0, y0, t0, options) 
        plotAverageVelocitiesTimePoint(QS, tp, options)
        plotPathlineVelocitiesTimePoint(QS, tp, options)
        plotStrainRateTimePoint(QS, tp, options) 
    end
    
    % Public methods, accessible from outside the class and reliant on 
    % properties of the class instance
    methods
        function setTime(QS, tt)
            % Set the current time of the dataset and clear current data
            % which was associated with the previously considered time
            %
            % Parameters
            % ----------
            % tt : int or float
            %   timePoint to set to be current, from available times in
            %   QS.xp.fileMeta.timePoints
            %
            if tt ~= QS.currentTime
                QS.clearTime() ;
            end
            QS.currentTime = tt ;
            QS.xp.setTime(tt) ;
        end
        
        function clearTime(QS)
            % clear current timepoint's data for QS instance
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
            QS.currentVelocity.piv3d2x = struct() ;
        end
        
        function t0 = t0set(QS, t0)
            % t0set(QS, t0) Set time offset to 1st fold onset or manually 
            if nargin < 2
                if exist(QS.fileName.fold, 'file')
                    % Note that fold_onset is in units of timepoints, not 
                    % indices into timepoints
                    load(QS.fileName.fold, 'fold_onset') ;
                    QS.t0 = min(fold_onset) ;
                else
                    error('No folding times saved to disk')
                end
            else
                QS.t0 = t0 ;
            end
            t0 = QS.t0 ;
        end
        
        function [acom_sm, pcom_sm] = getAPCOMSm(QS) 
            % Load the anterior and posterior 'centers of mass' ie the
            % endpoints of the object's centerline
            try
                acom_sm = h5read(QS.fileName.apdv, '/acom_sm') ;
                pcom_sm = h5read(QS.fileName.apdv, '/pcom_sm') ;
            catch
                [acom_sm, pcom_sm] = QS.computeAPDCOMs() ;
            end
        end
        
        function [rot, trans] = getRotTrans(QS)
            % Load the translation to put anterior to origin and AP axis
            % along x axis 
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
            %[raw, pix, um, um_buff] = GETXYZLIMS(QS)
            % Grab each xyzlim from self, otherwise load from disk
            % full resolution pix
            if ~isempty(QS.plotting.xyzlim_raw)
                xyzlim_raw = QS.plotting.xyzlim_raw ;
            else
                try
                    xyzlim_raw = dlmread(QS.fileName.xyzlim_raw, ',', 1, 0) ; 
                    QS.plotting.xyzlim_raw = xyzlim_raw ;
                catch
                    [QS.plotting.xyzlim_raw, QS.plotting.xyzlim_pix, ...
                        QS.plotting.xyzlim_um, ...
                        QS.plotting.xyzlim_um_buff] = ...
                        QS.measureXYZLims() ;
                    xyzlim_raw = QS.plotting.xyzlim_raw ;
                end
            end
            % rotated scaled in full resolution pix
            if ~isempty(QS.plotting.xyzlim_pix)
                xyzlim_pix = QS.plotting.xyzlim_pix ;
            else
                try
                    xyzlim_pix = dlmread(QS.fileName.xyzlim_pix, ',', 1, 0) ; 
                    QS.plotting.xyzlim_pix = xyzlim_pix ;
                catch
                    [~, QS.plotting.xyzlim_pix, ...
                        QS.plotting.xyzlim_um, ...
                        QS.plotting.xyzlim_um_buff] = ...
                        QS.measureXYZLims() ;
                    xyzlim_pix = QS.plotting.xyzlim_pix ;
                end
            end
            % rotated scaled APDV in micron
            if ~isempty(QS.plotting.xyzlim_um)
                xyzlim_um = QS.plotting.xyzlim_um ;
            else
                try
                    xyzlim_um = dlmread(QS.fileName.xyzlim_um, ',', 1, 0) ;
                    QS.plotting.xyzlim_um = xyzlim_um ;
                catch
                    [~, ~, QS.plotting.xyzlim_um, ...
                        QS.plotting.xyzlim_um_buff] = ...
                        QS.measureXYZLims() ;
                    xyzlim_um = QS.plotting.xyzlim_um ;
                end
            end
            % rotated scaled APDV in micron, with padding
            if ~isempty(QS.plotting.xyzlim_um_buff)
                xyzlim_um_buff = QS.plotting.xyzlim_um_buff ;
            else
                try
                    xyzlim_um_buff = dlmread(QS.fileName.xyzlim_um_buff, ',', 1, 0) ;
                    QS.plotting.xyzlim_um_buff = xyzlim_um_buff ;
                catch
                    [~, ~, ~, QS.plotting.xyzlim_um_buff] = ...
                        QS.measureXYZLims() ;
                    xyzlim_um_buff = QS.plotting.xyzlim_um_buff ;
                end
            end
        end
        
        function getFeatures(QS, varargin)
            %GETFEATURES(QS, varargin)
            %   Load features of the QS object (those specied, or all of 
            %   them). Features include {'folds', 'fold_onset', 'ssmax', 
            %   'ssfold', 'rssmax', 'rssfold'}. 
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
        
        function setDataLimits(QS, tp, adjustlow_pctile, adjusthigh_pctile)
            % Use timepoint (tp) to obtain hard values for intensity limits
            % so that data is rescaled to fixed limits instead of
            % percentile. This is useful to avoid flickering of overall
            % intensity in data in which a few voxels vary a lot in
            % intensity.
            QS.xp.loadTime(tp);
            QS.xp.rescaleStackToUnitAspect();
            IV = QS.xp.stack.image.apply() ;
            try
                assert(adjusthigh_pctile > 0 && adjustlow_pctile < 100)
                assert(adjusthigh_pctile > adjustlow_pctile)
            catch
                error('adjustment values must be 0<=val<=100 and increasing')
            end
            adjustlow = prctile(IV{1}(:), adjustlow_pctile) ;
            adjusthigh = prctile(IV{1}(:), adjusthigh_pctile) ;
            QS.data.adjustlow = adjustlow ;
            QS.data.adjusthigh = adjusthigh ;
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
            elseif adjustlow < 100 && adjusthigh < 100
                disp('Taking custom limits for imadjustn as prctile')
                for ii = 1:length(IV)
                    IVii = IV{ii} ;
                    vlo = double(prctile( IVii(:) , adjustlow )) / double(max(IVii(:))) ;
                    vhi = double(prctile( IVii(:) , adjusthigh)) / double(max(IVii(:))) ;
                    disp(['--> ', num2str(vlo), ', ', num2str(vhi), ...
                        ' for ', num2str(adjustlow), '/', num2str(adjusthigh)])
                    IV{ii} = imadjustn(IVii, [double(vlo); double(vhi)]) ;
                end
            else
                % adjusthigh is > 100, so interpret as an intensity value
                disp('Taking custom limits for imadjustn as direct intensity limit values')
                for ii = 1:length(IV)
                    IVii = IV{ii} ;
                    vlo = double(adjustlow) ;
                    vhi = double(adjusthigh) ;
                    disp(['--> ', num2str(vlo), ', ', num2str(vhi), ...
                        ' for ', num2str(adjustlow), '/', num2str(adjusthigh)])
                    tmp = (double(IVii) - vlo) / (vhi - vlo) ;
                    tmp(tmp > (vhi - vlo)) = 1.0 ;
                    IV{ii} = uint16(2^16 * tmp) ;
                    % cast(tmp, class(IVii)) ;  
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
            no_piv3d2x = isempty(fieldnames(QS.currentVelocity.piv3d2x)) ;
            if (do_all || contains(varargin, 'piv3d2x')) && no_piv3d2x
                % Load 3D data for piv results
                piv3d2xfn = QS.fullFileBase.piv3d2x ;
                load(sprintf(piv3d2xfn, QS.currentTime), 'piv3dstruct') ;
                QS.currentVelocity.piv3d2x = piv3dstruct ;
            end
        end
        
        % APDV methods
        [acom_sm, pcom_sm] = computeAPDCOMs(QS, opts)
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
        [rot, trans, xyzlim_raw, xyzlim, xyzlim_um, xyzlim_um_buff] = ...
            alignMeshesAPDV(QS, alignAPDVOpts) 
        
        % Masked Data
        generateMaskedData(QS)
        alignMaskedDataAPDV(QS)
        plotSeriesOnSurfaceTexturePatch(QS, overwrite, metadat, ...
                                        TexturePatchOptions)
        
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
        generateCurrentCutMesh(QS, options)
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
        function mesh = getCurrentSPCutMeshSm(QS)
            if isempty(QS.currentMesh.spcutMesh)
                QS.loadCurrentSPCutMeshSm() ;
            end
            mesh = QS.currentMesh.spcutMeshSm ;
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
        
        % spcutMeshSm at DoubleRes (2x resolution)
        function generateSPCutMeshSm2x(QS, overwrite)
            % Double the resolution of spcutMeshSm meshes
            %
            if nargin < 2
                overwrite = false ;
            end
            for tp = QS.xp.fileMeta.timePoints 
                sp2xfn = sprintf(QS.fullFileBase.spcutMeshSm2x, tp) ;
                if overwrite || ~exist(sp2xfn, 'file')
                    mesh1x = load(sprintf(QS.fullFileBase.spcutMeshSm, tp),...
                        'spcutMeshSm') ;
                    mesh1x = mesh1x.spcutMeshSm ;
                    spcutMeshSm2x = QS.doubleResolution(mesh1x) ;
                    disp(['saving ' sp2xfn])
                    save(sp2xfn, 'spcutMeshSm2x')
                end
                sp2xfn = sprintf(QS.fullFileBase.spcutMeshSmRS2x, tp) ;
                spC2fn = sprintf(QS.fullFileBase.spcutMeshSmRSC2x, tp) ;
                if overwrite || ~exist(sp2xfn, 'file') || ...
                        ~exist(spC2fn, 'file')
                    mesh1x = load(...
                        sprintf(QS.fullFileBase.spcutMeshSmRS, tp), ...
                        'spcutMeshSmRS') ;
                    mesh1x = mesh1x.spcutMeshSmRS ;
                    [spcutMeshSmRS2x, spcutMeshSmRSC2x] = ...
                        QS.doubleResolution(mesh1x) ;
                    disp(['saving ' sp2xfn])
                    save(sp2xfn, 'spcutMeshSmRS2x')
                    disp(['saving ' spC2fn])
                    save(spC2fn, 'spcutMeshSmRSC2x')
                end
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
        function getPIV(QS, options)
            % Load PIV results and store in QS.piv if not already loaded
            if isempty(fieldnames(QS.piv.raw)) || isempty(QS.piv.Lx) ...
                    || isempty(QS.piv.Lx) 
                % Load raw PIV results
                if nargin > 1
                    QS.loadPIV(options)
                else
                    QS.loadPIV() 
                end
            end
        end
        function loadPIV(QS, options)
            % Load PIV results from disk and store in QS.piv
            if ~isempty(fieldnames(QS.piv.raw)) && isempty(QS.piv.Lx) ...
                    && isempty(QS.piv.Lx) 
                disp("WARNING: Overwriting QS.piv with piv from disk")
            end
            QS.piv.raw = load(QS.fileName.pivRaw) ;  
            timePoints = QS.xp.fileMeta.timePoints ;
            if strcmp(QS.piv.imCoords, 'sp_sme')
                im0 = imread(sprintf(QS.fullFileBase.im_sp_sme, ...
                    timePoints(1))) ;
                % for now assume all images are the same size
                QS.piv.Lx = size(im0, 1) * ones(length(timePoints), 1) ;
                QS.piv.Ly = size(im0, 2) * ones(length(timePoints), 1) ;
            end
        end
        measurePIV3d(QS, options)
        measurePIV3dDoubleResolution(QS, options)
        % Note: To timeAverage Velocities at Double resolution, pass
        % options.doubleResolution == true
        
        %% Pathlines
        measurePullbackPathlines(QS, options)
        function getPullbackPathlines(QS, t0, varargin)
            % Discern if we must load pathlines or if already loaded
            if nargin > 1 
                if QS.pathlines.t0 ~= t0
                    % The timestamp at which pathlines form grid that is 
                    % requested is different than the one that is loaded,
                    % if any are indeed already loaded. Therefore, we load
                    % anew
                    if nargin > 2
                        % pass varargin along to load method
                        QS.loadPullbackPathlines(t0, varargin)
                    else
                        QS.loadPullbackPathlines(t0)
                    end
                end
            else
                % No t0 supplied, assume t0 is the same as what is stored
                % in QS.pathlines.t0, if any is already stored (ie if any
                % pathlines are already loaded)
                if isempty(QS.pathlines.t0)
                    % no pathlines loaded. Load here
                    if nargin > 2
                        QS.loadPullbackPathlines(t0, varargin)
                    elseif narargin > 1
                        QS.loadPullbackPathlines(t0)
                    else
                        QS.loadPullbackPathlines()
                    end
                else
                    % There are pathlines loaded already. Which are
                    % requested here in varargin? First check if varargin 
                    % is empty or not  
                    if nargin > 2            
                        if any(contains(varargin, 'pivPathlines'))
                            if isempty(QS.pathlines.piv)
                                QS.loadPullbackPathlines(t0, 'pivPathlines')
                            end          
                        end
                        if any(contains(varargin, 'vertexPathlines'))
                            if isempty(QS.pathlines.vertices)
                                QS.loadPullbackPathlines(t0, 'vertexPathlines')
                            end            
                        end            
                        if any(contains(varargin, 'facePathlines'))
                            if isempty(QS.pathlines.faces)
                                QS.loadPullbackPathlines(t0, 'facePathlines')
                            end
                        end
                    else
                        % varargin is not supplied, so load all three if
                        % not already loaded
                        % First grab t0
                        if nargin < 2
                            t0 = QS.t0set() ;
                        end
                        if isempty(QS.pathlines.piv)
                            QS.loadPullbackPathlines(t0, 'pivPathlines')
                        end            
                        if isempty(QS.pathlines.vertices)
                            QS.loadPullbackPathlines(t0, 'vertexPathlines')
                        end            
                        if isempty(QS.pathlines.faces)
                            QS.loadPullbackPathlines(t0, 'facePathlines')
                        end
                    end
                end
                        
                    
            end
        end
        function loadPullbackPathlines(QS, t0, varargin)
            if nargin < 2
                t0 = QS.t0set() ;
            elseif isempty(t0)
                t0 = QS.t0set() ;
            else
                try
                    assert(isnumeric(t0))
                catch
                    error('t0 supplied must be numeric')
                end
            end
            % assign t0 as the pathline t0
            QS.pathlines.t0 = t0 ;
            if nargin < 3
                varargin = {'pivPathlines', 'vertexPathlines', ...
                    'facePathlines'} ;
            end
            if any(contains(varargin, 'pivPathlines'))
                load(sprintf(QS.fileName.pathlines.XY, t0), 'pivPathlines')
                QS.pathlines.piv = pivPathlines ;
            end
            if any(contains(varargin, 'vertexPathlines'))
                load(sprintf(QS.fileName.pathlines.vXY, t0), 'vertexPathlines')
                QS.pathlines.vertices = vertexPathlines ;
            end
            if any(contains(varargin, 'facePathlines'))
                load(sprintf(QS.fileName.pathlines.fXY, t0), 'facePathlines')
                QS.pathlines.faces = facePathlines ;
            end
        end
                
        %% Velocities -- Lagrangian Averaging
        timeAverageVelocities(QS, samplingResolution, options)
        function loadVelocityAverage(QS, varargin)
            % Load and pack into struct
            if isempty(varargin)
                varargin = {'v3d', 'v2dum', 'v2d', 'vn', 'vf', 'vv'};
            end
            if any(strcmp(varargin, 'v3d'))
                load(QS.fileName.pivAvg.v3d, 'vsmM') ;
                QS.velocityAverage.v3d = vsmM ;
            end
            if any(strcmp(varargin, 'v2dum'))
                load(QS.fileName.pivAvg.v2dum, 'v2dsmMum') ;
                QS.velocityAverage.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'vn'))
                load(QS.fileName.pivAvg.vn, 'vnsmM') ;
                QS.velocityAverage.vn = vnsmM ;
            end
            if any(strcmp(varargin, 'vf'))
                load(QS.fileName.pivAvg.vf, 'vfsmM') ;
                QS.velocityAverage.vf = vfsmM ;
            end
            if any(strcmp(varargin, 'vv'))
                load(QS.fileName.pivAvg.vv, 'vvsmM') ;
                QS.velocityAverage.vv = vvsmM ;
            end
        end
        function getVelocityAverage(QS, varargin)
            % todo: check if all varargin are already loaded
            loadVelocityAverage(QS, varargin{:})
        end
        function loadVelocityAverage2x(QS, varargin)
            % Load and pack into struct
            if isempty(varargin)
                varargin = {'v3d', 'v2dum', 'v2d', 'vn', 'vf', 'vv'};
            end
            if any(strcmp(varargin, 'v3d'))
                load(QS.fileName.pivSimAvg2x.v3d, 'vsmM') ;
                QS.velocityAverage2x.v3d = vsmM ;
            end
            if any(strcmp(varargin, 'v2dum'))
                load(QS.fileName.pivSimAvg2x.v2dum, 'v2dsmMum') ;
                QS.velocityAverage2x.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'v2d'))
                load(QS.fileName.pivSimAvg2x.v2dum, 'v2dsmMum') ;
                QS.velocityAverage2x.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'vn'))
                load(QS.fileName.pivSimAvg2x.vn, 'vnsmM') ;
                QS.velocityAverage2x.vn = vnsmM ;
            end
            if any(strcmp(varargin, 'vf'))
                load(QS.fileName.pivSimAvg2x.vf, 'vfsmM') ;
                QS.velocityAverage2x.vf = vfsmM ;
            end
            if any(strcmp(varargin, 'vv'))
                load(QS.fileName.pivSimAvg2x.vv, 'vvsmM') ;
                QS.velocityAverage2x.vv = vvsmM ;
            end
        end
        function getVelocityAverage2x(QS, varargin)
            if isempty(QS.velocityAverage2x.v3d)
                loadVelocityAverage2x(QS, varargin{:})
            end
        end
        plotTimeAvgVelocities(QS, options)
        helmholtzHodge(QS, options)
        measurePathlineVelocities(QS, options)
        plotPathlineVelocities(QS, options)
        
        %% Velocities -- simple/surface-Lagrangian averaging
        timeAverageVelocitiesSimple(QS, samplingResolution, options)
        function loadVelocitySimpleAverage(QS, varargin)
            % Load and pack into struct
            if any(strcmp(varargin, 'v3d'))
                load(QS.fileName.pivSimAvg.v3d, 'vsmM') ;
                QS.velocitySimpleAverage.v3d = vsmM ;
            end
            if any(strcmp(varargin, 'v2dum'))
                load(QS.fileName.pivSimAvg.v2dum, 'v2dsmMum') ;
                QS.velocitySimpleAverage.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'vn'))
                load(QS.fileName.pivSimAvg.vn, 'vnsmM') ;
                QS.velocitySimpleAverage.vn = vnsmM ;
            end
            if any(strcmp(varargin, 'vf'))
                load(QS.fileName.pivSimAvg.vf, 'vfsmM') ;
                QS.velocitySimpleAverage.vf = vfsmM ;
            end
            if any(strcmp(varargin, 'v2v'))
                load(QS.fileName.pivSimAvg.vf, 'vvsmM') ;
                QS.velocitySimpleAverage.vv = vvsmM ;
            end
        end
        function getVelocitySimpleAverage(QS, varargin)
            % todo: check if all varargin are already loaded
            loadVelocitySimpleAverage(QS, varargin{:})
        end
        function loadVelocitySimpleAverage2x(QS, varargin)
            % Load and pack into struct
            if isempty(varargin)
                varargin = {'v3d', 'v2dum', 'v2d', 'vn', 'vf', 'vv'};
            end
            if any(strcmp(varargin, 'v3d'))
                load(QS.fileName.pivSimAvg2x.v3d, 'vsmM') ;
                QS.velocitySimpleAverage2x.v3d = vsmM ;
            end
            if any(strcmp(varargin, 'v2dum'))
                load(QS.fileName.pivSimAvg2x.v2dum, 'v2dsmMum') ;
                QS.velocitySimpleAverage2x.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'v2d'))
                load(QS.fileName.pivSimAvg2x.v2dum, 'v2dsmMum') ;
                QS.velocitySimpleAverage2x.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'vn'))
                load(QS.fileName.pivSimAvg2x.vn, 'vnsmM') ;
                QS.velocitySimpleAverage2x.vn = vnsmM ;
            end
            if any(strcmp(varargin, 'vf'))
                load(QS.fileName.pivSimAvg2x.vf, 'vfsmM') ;
                QS.velocitySimpleAverage2x.vf = vfsmM ;
            end
            if any(strcmp(varargin, 'vv'))
                load(QS.fileName.pivSimAvg2x.vv, 'vvsmM') ;
                QS.velocitySimpleAverage2x.vv = vvsmM ;
            end
        end
        function getVelocitySimpleAverage2x(QS, varargin)
            if isempty(QS.velocitySimpleAverage2x.v3d)
                loadVelocitySimpleAverage2x(QS, varargin{:})
            end
        end
        % NOTE: the following have a simple option for averagingStyle
        % plotTimeAvgVelocities(QS, options)
        % helmholtzHodge(QS, options)
        
        %% compressible/incompressible flow on evolving surface
        measureMetricKinematics(QS, options)
        plotMetricKinematics(QS, options)
        measurePathlineMetricKinematics(QS, options)
        plotPathlineMetricKinematics(QS, options)
        
        %% 
        measureMetricStrainRate(QS, options)
        measureStrainRate(QS, options)
        plotStrainRate(QS, options)
        measurePathlineStrainRate(QS, options)
        plotPathlineStrainRate(QS, options)
        
        %% timepoint-specific coordinate transformations
        sf = interpolateOntoPullbackXY(QS, XY, scalar_field, options)
        
    end
    
    methods (Static)
        function uv = XY2uv(im, XY, doubleCovered, umax, vmax)
            %XY2uv(im, XY, doubleCovered, umax, vmax)
            %   Map pixel positions (1, sizeImX) and (1, sizeImY) to
            %   (0, umax) and (0, vmax) of pullback space if singleCover,
            %   or y coords are mapped to (-0.5, 1.5)*vmax if doubleCover
            % 
            % NOTE THAT MAP IS
            % [xesz, yesz] = [size(im, 1), size(im, 2)]
            % uv(:, 1) = umax * (XY(:, 1) - 1) / xesz ;
            % uv(:, 2) = vmax * 2.0 * (XY(:, 2) - 1) / yesz - 0.5 ;
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
            %   "bottom" half repeats above the top. That is, 
            %   consider im to be a double cover in Y (periodic in Y and 
            %   covers pullback space twice (-0.5 * Ly, 1.5 * Ly)
            % umax : float
            %   extent of pullback mesh coordinates in u direction (X)
            % vmax : float
            %   extent of pullback mesh coordinates in v direction (Y)
            %   before double covering/tiling
            
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
            % Parameters
            % ----------
            % im : NxM numeric array or length(2) int array
            %   2D image into whose pixel space to map or size(im)
            % uv : Q*2 numeric array
            %   mesh coordinates to convert to pullback pixel space (XY)
            % doubleCovered: bool
            %   the image is a double cover of the pullback (extended/tiled
            %   so that the "top" half repeats below the bottom and the
            %   "bottom" half repeats above the top. That is, 
            %   consider im to be a double cover in Y (periodic in Y and 
            %   covers pullback space twice (-0.5 * Ly, 1.5 * Ly)
            % umax : float
            %   extent of pullback mesh coordinates in u direction (X)
            % vmax : float 
            %   extent of pullback mesh coordinates in v direction (Y)
            %   before double covering/tiling
            %
            % Returns
            % -------
            % XY : N x 2 float array
            %   positions of uv coordinates in pullback pixel space
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
            
            if any(size(im) > 2) 
                Xsz = size(im, 2) ;
                Ysz = size(im, 1) ;
            else
                % Interpret im as imsize
                Xsz = im(1) ;
                Ysz = im(2) ;
            end
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
        
        function [xx, yy] = clipXY(xx, yy, Lx, Ly)
            % Clip x at (1, Lx) and clip Y as periodic (1=Ly, Ly=1), for
            % image that is periodic in Y. Consider Y in [1, Ly]. If the
            % pullback is a doubleCover, Ly = 2*mesh width in pixels
            %
            % Parameters
            % ----------
            % xx : 
            % yy : 
            % Lx : int
            %   number of pixels along x dimension
            % Ly : int
            %   number of pixels along y dimension
            % 
            % Returns
            % -------
            % [xx, yy] : x and y coordinates clipped to [1, Lx] and [1, Ly]
            xx(xx > Lx) = Lx ;
            xx(xx < 1 ) = 1 ;
            yy(yy > Ly) = yy(yy > Ly) - Ly + 1;
            yy(yy < 1) = yy(yy < 1) + Ly ;
        end
        
        function XY = doubleToSingleCover(XY, Ly)
            % detect if XY is passed as a pair of grids
            if length(size(XY)) > 2 && size(XY, 2) > 2
                % XY is a pair of position grids each as 2D arrays. Clip Y
                tmp = XY(:, :, 2) ;
                tmp(tmp < Ly * .25) = tmp(tmp < Ly * .25) + Ly * 0.5 ;
                tmp(tmp > Ly * .75) = tmp(tmp > Ly * .75) - Ly * 0.5 ;
                XY(:, :, 2) = tmp ;
            elseif size(XY, 2) == 2
                % XY is input as Nx2 array
                XY(XY(:, 2) < Ly * .25, 2) = XY(XY(:, 2) < Ly * .25, 2) + Ly * 0.5 ;
                XY(XY(:, 2) > Ly * .75, 2) = XY(XY(:, 2) > Ly * .75, 2) - Ly * 0.5 ;
            else
                error('XY must be passed as Nx2 or QxRx2 array')
            end
        end
        
        % function uv2pix_old(im, aspect)
        %     % map from network xy to pixel xy
        %     % Note that we need to flip the image (Yscale - stuff) since saved ims had
        %     % natural ydirection.
        %     % Assume here a coord system xy ranging from (0, xscale) and (0, 1) 
        %     % maps to a coord system XY ranging from (0, Yscale * 0.5) and (0, Yscale)
        %     x2Xpix = @(x, Yscale, xscale) (Yscale * 0.5) * aspect * x / xscale ;
        %     % y2Ypix = @(y, h, Yscale) Yscale - (Yscale*0.5)*(y+0.5) + h ;
        %     y2Ypix = @(y, h, Yscale) (Yscale*0.5)*(y+0.5) + h ;
        % 
        %     dx2dX = @ (y, Yscale, xscale) (Yscale * 0.5) * aspect * x / xscale ;
        %     dy2dY = @ (y, Yscale) (Yscale*0.5)*y ;
        % 
        %     % Now map the coornates
        % end
        
        [cutMesh, cutMeshC] = doubleResolution(cutMesh, preview)
        
        function [mag_ap, theta_ap] = dvAverageNematic(magnitude, theta)
            % 
            ap_x = magnitude .* cos(2*theta) ;
            ap_y = magnitude .* sin(2*theta) ;
            ap_xy = [mean(ap_x, 2) , mean(ap_y, 2)] ;
            mag_ap = vecnorm(ap_xy, 2, 2) ;
            theta_averages = atan2(ap_xy(:, 2), ap_xy(:, 1)) ;
            theta_ap = 0.5 * mod(theta_averages, 2*pi) ;
        end
       
    end
    
end