function helmholtzHodge(QS, options) 
%helmholtzHodge(QS, options) 
%   Take divergence and "curl" on 2d evolving surface in 3d
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%   averagingStyle : str ('Lagrangian', 'simple')
%       style in which velocities were averaged over time
%   alphaVal : float
%       the opacity of the heatmap to overlay
%   invertImage : bool
%       invert the data pullback under the velocity map
%
%
% NPMitchell 2020

%% Default options
overwrite = false ;
overwrite_images = false ;
qsubU = 5 ; 
qsubV = 10 ;
niter_smoothing = [1000, 1000, 1000] ;
plot_dec_pullback = true ;
plot_dec_texturepatch = false ;
preview = false ;
pivimCoords = QS.piv.imCoords ;
lambda_smooth = 0.01 ;
lambda_mesh = 0.001 ;
samplingResolution = '1x' ;
averagingStyle = 'Lagrangian' ;

%% Unpack options
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'qsubU')
    qsubU = options.qsubU ;
end
if isfield(options, 'qsubV')
    qsubV = options.qsubV ;
end
if isfield(options, niter_smoothing)
    niter_smoothing = options.niter_smoothing ;
end
if isfield(options, plot_dec_pullback)
    plot_dec_pullback = options.plot_dec_pullback ;
end
if isfield(options, plot_dec_texturepatch)
    plot_dec_texturepatch = options.plot_dec_texturepatch ;
end
if isfield(options, 'pivimCoords')
    pivimCoords = options.pivimCoords ;
end
if isfield(options, 'lambda')
    lambda_smooth = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
end

% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
else 
    error("Could not parse samplingResolution: set to '1x' or '2x'")
end

%% Unpack QS
timePoints = QS.xp.fileMeta.timePoints ;
fname = QS.fileBase.name ;
[rot, trans] = QS.getRotTrans() ;
resolution = QS.APDV.resolution ;
if doubleResolution
    nU = 2 * QS.nU - 1 ;
    nV = 2 * QS.nV - 1 ;
    piv3dFileBase = QS.fullFileBase.piv3d2x ;
else
    nU = QS.nU ;
    nV = QS.nV ;
    piv3dFileBase = QS.fullFileBase.piv3d ;
end
if strcmp(averagingStyle, 'Lagrangian')
    if doubleResolution
        QS.getVelocityAverage2x('vf', 'v2dum') ;
        vf = QS.velocityAverage2x.vf ;
        v2dum = QS.velocityAverage2x.v2dum ;
        decDirRoot = QS.dir.pivAvgDEC2x ;
        decFnBase = QS.fullFileBase.decAvg2x ;
    else
        QS.getVelocityAverage('vf', 'v2dum') ;
        vf = QS.velocityAverage.vf ;
        v2dum = QS.velocityAverage.v2dum ;
        decDirRoot = QS.dir.pivAvgDEC ;
        decFnBase = QS.fullFileBase.decAvg ;    
    end
elseif strcmp(averagingStyle, 'simple') 
    if doubleResolution
        QS.getVelocitySimpleAverage2x('vf', 'v2dum') ;
        vf = QS.velocitySimpleAverage2x.vf ;
        v2dum = QS.velocitySimpleAverage2x.v2dum ;
        decDirRoot = QS.dir.pivSimAvgDEC2x ;
        nU = 2 * QS.nU - 1 ;
        nV = 2 * QS.nV - 1 ;
        piv3dFileBase = QS.fullFileBase.piv3d2x ;
        decFnBase = QS.fullFileBase.decSimAvg2x ;
    else
        QS.getVelocitySimpleAverage('vf', 'v2dum') ;
        vf = QS.velocitySimpleAverage.vf ;
        v2dum = QS.velocitySimpleAverage.v2dum ;
        decDirRoot = QS.dir.pivSimAvgDEC ;
        nU = QS.nU ;
        nV = QS.nV ;
        piv3dFileBase = QS.fullFileBase.piv3d ;
        decFnBase = QS.fullFileBase.decSimAvg ;
    end
else
    error('averagingStyle not recognized. Use Lagrangian or simple')
end
 
t0 = QS.t0set() ;
[~, ~, ~, xyzlim] = QS.getXYZLims() ;

%% define the output dirs
decDir = decDirRoot.data ;
dDir2d = decDirRoot.div2d ;
dDir3d = decDirRoot.div3d ;
dDir3dt = decDirRoot.div3dTexture ;
rDir2d = decDirRoot.rot2d ;
rDir3d = decDirRoot.rot3d ;
rDir3dt = decDirRoot.rot3dTexture ;
hDir2d = decDirRoot.harm2d ;
hDir3d = decDirRoot.harm3d ;
% create the output dirs
dirs2do = {decDir, dDir2d, dDir3d, dDir3dt, ...
    rDir2d, rDir3d, rDir3dt, ...
    hDir2d, hDir3d} ;
for i = 1:length(dirs2do)
    ensureDir(dirs2do{i}) ;
end

%% Perform Decomposition
% Consider each timepoint and plot the div and curl
tidxAll = 1:length(timePoints(1:end-1)) ;
tidx2do = 1:10:length(timePoints(1:end-1)) ;
tidxOther = setdiff(tidxAll, tidx2do) ;
tidx2do = [tidx2do, tidxOther] ;
for tidx = tidx2do
    tp = timePoints(tidx) ;
    disp(['t = ', num2str(tp)])
    
    % DEC data filename to save
    decDataFn = sprintf(decFnBase, tp) ;
    
    % Prepare filenames
    div2dfn = fullfile(dDir2d, [sprintf(fname, tp) '_div2d.png']) ; 
    div3dfn_lateral = fullfile(dDir3d, ['lateral_' sprintf(fname, tp) '_div3d.png']) ;
    div3dfn_ventral = fullfile(dDir3d, ['ventral_' sprintf(fname, tp) '_div3d.png']) ;
    div3dtfn = fullfile(dDir3dt, [sprintf(fname, tp) '_divt3d.png']) ;
    rot2dfn = fullfile(rDir2d, [sprintf(fname, tp) '_rot2d.png']) ;
    rot3dfn_lateral = fullfile(rDir3d, ['lateral_' sprintf(fname, tp) '_rot3d.png']) ;
    rot3dfn_ventral = fullfile(rDir3d, ['ventral_' sprintf(fname, tp) '_rot3d.png']) ;
    rot3dtfn = fullfile(rDir3dt, [sprintf(fname, tp) '_rot3d.png']) ;
    harm2dfn = fullfile(hDir2d, [sprintf(fname, tp) '_harm2d.png']) ; 
    harm3dfn = fullfile(hDir3d, [sprintf(fname, tp) '_harm3d.png']) ;

    % Discern whether to plot results in 2d and 3d or in texturepatch
    plot_dec_pullback_tidx = plot_dec_pullback && ...
        (~exist(div2dfn, 'file') || ...
        ~exist(div3dfn_lateral, 'file') || ...
        ~exist(div3dfn_ventral, 'file') || overwrite || overwrite_images) ;
    plot_dec_texturepatch_tidx = plot_dec_texturepatch && ...
        (overwrite || ~exist(div2dfn, 'file') || ...
        ~exist(div3dfn_lateral, 'file') || ...
        ~exist(div3dfn_ventral, 'file') || overwrite || overwrite_images) ;
    
    if ~exist(decDataFn, 'file') || overwrite
        if ~exist(decDataFn, 'file')
            disp(['Creating DEC: ' decDataFn])
        else
            disp(['Overwriting DEC: ', decDataFn])
        end
        
        % Load piv3d
        piv3d = load(sprintf(piv3dFileBase, tp)) ;
        piv3d = piv3d.piv3dstruct ;
    
        % Obtain smoothed velocities on all faces
        vfsm = squeeze(vf(tidx, :, :)) ;
        v2dsmum_ii = squeeze(v2dum(tidx, :, :)) ;
        
        % Use current time's tiled smoothed mesh
        % Note: vfsmM is in um/min rs
        FF = piv3d.m0f ;   % #facesx3 float: mesh connectivity list
        V2D = piv3d.m0XY ; % Px2 float: 2d mesh vertices in pullback image pixel space
        v3drs = QS.xyz2APDV(piv3d.m0v3d) ;
        
        % Check it
        % figure;
        % set(gcf, 'visible', 'on')
        % plot(v3drs)
        % pause(1)
        % if tp > timePoints(1)
        %     plot(v3drs - oldv)
        %     disp('here')
        % end
        % oldv = v3drs ;
        % close all
        
        % Grab cutMesh from piv3d. Could grab from disk instead...
        cutM.f = FF ;
        cutM.u = V2D ;
        cutM.v = v3drs ;
        cutM.nU = nU ;
        cutM.nV = nV ;
        
        %% Compute divs and rots. Note that smoothing occurs inside func
        disp('Decomposing flow into div/rot...')
        Options = struct() ;
        Options.lambda = lambda_smooth ;   
        Options.lambda_mesh = lambda_mesh ;        
        [divs, rots, harms, glueMesh] = ...
            helmHodgeDECRectGridPullback(cutM, vfsm, Options,...
            'niterSmoothing', niter_smoothing, ...
            'clipDiv', [-5, 5], 'clipRot', [-0.5, 0.5], ...
            'preview', preview, 'method', 'both') ;
        
        %% save divs, rots, and harms as structs in .mat file
        disp(['Saving DEC t=' num2str(tp) ': ' decDataFn])
        save(decDataFn, 'divs', 'rots', 'harms')
    elseif plot_dec_pullback_tidx || plot_dec_texturepatch_tidx
        disp(['Loading DEC t=' num2str(tp) ': ' decDataFn])
        load(decDataFn, 'divs', 'rots', 'harms')
    
        % Load piv3d
        piv3d = load(sprintf(piv3dFileBase, tp)) ;
        piv3d = piv3d.piv3dstruct ;

        % Obtain smoothed velocities on all faces
        vfsm = squeeze(vf(tidx, :, :)) ;
        v2dsmum_ii = squeeze(v2dum(tidx, :, :)) ;

        % Use current time's tiled smoothed mesh
        % Note: vfsmM is in um/min rs
        FF = piv3d.m0f ;   % #facesx3 float: mesh connectivity list
        V2D = piv3d.m0XY ; % Px2 float: 2d mesh vertices in pullback image pixel space
        v3drs = QS.xyz2APDV(piv3d.m0v3d) ;

        % Grab cutMesh from piv3d. Could grab from disk instead...
        cutM.f = FF ;
        cutM.u = V2D ;
        cutM.v = v3drs ;
        cutM.nU = nU ;
        cutM.nV = nV ;
    end

    %% Plot results
    disp('Plotting div/rot...')
    if plot_dec_pullback_tidx
        close all
        set(gcf, 'visible', 'off')
        if strcmp(pivimCoords, 'sp_sme')
            im = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
            ylims = [0.25 * size(im, 1), 0.75 * size(im, 1)] ;
        else
            error(['Have not coded for this pivimCoords option. Do so here: ' pivimCoords])
        end
        im = cat(3, im, im, im) ;  % convert to rgb for no cmap change
        addTitleStr = [': $t=$', num2str((tp - t0)*QS.timeInterval), ...
                       ' ', QS.timeUnits] ;
        Options = struct() ;
        Options.addTitleStr = addTitleStr ;
        Options.div2dfn = div2dfn ;
        Options.div3dfn_lateral = div3dfn_lateral ;
        Options.div3dfn_ventral = div3dfn_ventral ;
        Options.rot2dfn = rot2dfn ;
        Options.rot3dfn_lateral = rot3dfn_lateral ;
        Options.rot3dfn_ventral = rot3dfn_ventral ;
        Options.qsubU = qsubU ; 
        Options.qsubV = qsubV ;
        Options.sscaleDiv = 0.1 ;
        Options.sscaleRot = 0.1 ;
        Options.qscaleDiv3d = 0 ;
        Options.qscaleRot3d = 0 ;
        Options.qscaleDiv2d = 0 ;
        Options.qscaleRot2d = 0 ;
        Options.xyzlim = xyzlim ;
        opts2d.xlim = [0, size(im, 1)] ;
        opts2d.ylim = [0.25 * size(im, 2), 0.75 * size(im, 2) ] ;
        xy = {piv3d.x0, piv3d.y0} ;
        plotHelmHodgeDECPullback(im, cutM, vfsm, xy, v2dsmum_ii, ...
            divs, rots, Options, opts2d) ;
    end
    if plot_dec_texturepatch_tidx
        % load current timepoint
        % (3D data for coloring mesh pullback)
        QS.setTime(tp) ;
        QS.getCurrentData() ;     
        IV = QS.currentData.IV ;
        IV = imcomplement(IV) ; % does this work? debug 2020
        % IV = max(IV(:)) - IV ; % used to do this.

        addTitleStr = ['$t=$', num2str(tp - t0)] ;
        Options.addTitleStr = addTitleStr ;
        Options.div3dfn = div3dtfn ;
        Options.rot3dfn = rot3dtfn ;
        Options.sscaleDiv = 0.5 ;
        Options.sscaleRot = 0.2 ;
        Options.xyzlim = xyzlim ;

        % Load the cutmesh vertices and normals
        cutM.v = piv3d.m0v3d ;
        cutM.v3drs = ((rot * cutM.v')' + trans) * resolution ;
        plotHelmHodgeDECTexture3d(IV, cutM, divs, rots, rot, trans, Options) ;
    end
end
