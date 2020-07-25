function timeAverageVelocities(QS, options)
%measurePullbackStreamlines(QS, options)
%   Use pathlines of optical flow in pullback space to query velocities
%   and average along pathlines.
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%
% NPMitchell 2020

%% Default options
overwrite = false ;
timePoints = QS.xp.fileMeta.timePoints ;
samplingResolution = '1x' ;
imethod = 'linear' ;

%% Unpack options
% Default values for options are to use sphi smoothed extended coords
% as PIV reference coord sys
pivimCoords = QS.piv.imCoords ;

if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if strcmp(pivimCoords(end), 'e')
    doubleCovered = true ;
else
    doubleCovered = false ;
end

%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
else 
    error("Could not parse samplingResolution: set to '1x' or '2x'")
end

%% Unpack QS
pivDir = QS.dir.piv ;
piv3dfn = QS.fullFileBase.piv3d ;
ntps = length(timePoints) ;
% [rot, ~] = QS.getRotTrans() ;
% resolution = QS.APDV.resolution ; 
[~, ~, ~, xyzlim_APDV] = QS.getXYZLims() ;
axis_order = QS.data.axisOrder ;
blue = QS.plotting.colors(1, :) ;
red = QS.plotting.colors(2, :) ;
green = QS.plotting.colors(4, :) ;
t0 = QS.t0set() ;
timePoints = QS.xp.fileMeta.timePoints ;

%% Perform/Load simple averaging
disp('Performing/Loading simple averaging')
% Create directories
fileNames = QS.fileName.pathlines ;
% Apply t0 to fileNames
fieldnames = fields(fileNames) ;
for qq = 1:length(fieldnames)
    fileNames.(fieldnames{qq}) = sprintf(fileNames.(fieldnames{qq}), t0) ;
end

% Check if the time smoothed velocities exist already
% 2d velocities (pulled back), scaled by dilation of metric are v2dum 
% 2D velocities (pulled back) are v2d
% normal velocities on fieldfaces are vn
% 3d velocities on fieldfaces are v3d
% vertex-based velocities are vv
% face-based velocities are vf

QS.clearTime() ;

%% Build grids for averaging
if ~exist(fileNames.v2dum, 'file') || ~exist(fileNames.v2d, 'file') || ...
        ~exist(fileNames.vn, 'file') || ~exist(fileNames.v3d, 'file') || ...
        ~exist(fileNames.vf, 'file') || ~exist(fileNames.vv, 'file') || ...
        ~exist(fileNames.v2dsmum, 'file') || ~exist(fileNames.v2dsm, 'file') || ...
        ~exist(fileNames.vnsm, 'file') || ~exist(fileNames.v3dsm, 'file') || ...
        ~exist(fileNames.vfsm, 'file') || ~exist(fileNames.vvsm, 'file') || ...
        overwrite
    
    disp('Could not find time-smoothed velocities on disk')
    disp('Computing them...')
    
    % Load lagrangian pathlines
    disp('Loading pathlines: XY')
    load(sprintf(QS.fileName.pathlines.XY, t0), 'pivPathlines')
    load(sprintf(QS.fileName.pathlines.vXY, t0), 'vertexPathlines')
    load(sprintf(QS.fileName.pathlines.fXY, t0), 'facePathlines')
    
    ntps = length(timePoints)-1;
    
    %% Now load 3d piv and smooth in Lagrangian coords (along streamlines)
    first = true ; 
    for tidx = 1:ntps
        tp = timePoints(tidx) ;
        % dt = timePoints(i + 1) - tp ;
        
        disp(['Filling in velocity matrices, t=' num2str(tp)])
        QS.setTime(tp) ;
        if doubleResolution
            QS.getCurrentVelocity('piv3d2x') ;
            piv3d = QS.currentVelocity.piv3d2x ;
        else
            QS.getCurrentVelocity('piv3d') ;
            piv3d = QS.currentVelocity.piv3d ;
        end
        
        % Allocate memory if this is the first timestep. Assume all grids
        % are equally sized.
        if first
            vM = zeros(ntps, size(piv3d.v0_rs, 1), size(piv3d.v0_rs, 2));
            vfM = zeros(ntps, size(piv3d.v3dfaces, 1), size(piv3d.v3dfaces, 2)); 
            vvM = zeros(ntps, ...
                size(piv3d.v3dvertices, 1), ...
                size(piv3d.v3dvertices, 2)); 
            vnM = zeros(ntps, size(piv3d.v0n_rs, 1), size(piv3d.v0n_rs, 2));
            v2dM = zeros(ntps, size(piv3d.v0t2d, 1), size(piv3d.v0t2d, 2));
            v2dMum = zeros(ntps, size(piv3d.v0t2d, 1), size(piv3d.v0t2d, 2));
            first = false ;
        end
        
        % Assert no NaNs    
        try
            assert(~any(isnan(piv3d.v0_rs(:))))
            assert(~any(isnan(piv3d.v3dfaces_rs(:))))
            assert(~any(isnan(piv3d.v0n_rs(:))))
            assert(~any(isnan(piv3d.v0t2d(:))))
        catch
           % disp('inpainting NaNs in pt0 & pt1')
           error(['There are NaNs in the velocity data. Could use ', ...
               'inpaint_nans, but why/how are they there?'])
           % pt0 = inpaint_nans(pt0) ;
           % pt1 = inpaint_nans(pt1) ;
           close all
           figure ;
           scatter(piv3d.x0(:), piv3d.y0(:), 10, piv3d.v0_rs(:, 1))
           bad = find(isnan(piv3d.v0_rs(:, 1))) ;
           hold on; 
           xx1 = piv3d.x0(:) ;
           yy1 = piv3d.y0(:) ;
           scatter(xx1(bad), yy1(bad), 30, 'k')
        end
        
        %% For all pathlines, query velocity at this current time
        % Query velocities at eval coords, mesh vertices, mesh faces via 
        % interpolation. Also query vn, v2d (ie vt), and v2dum (vt/|g|). 
               
        % Load streamline positions at this timePoint
        XX = pivPathlines.XX(tidx, :, :) ;
        YY = pivPathlines.YY(tidx, :, :) ;
        fX = facePathlines.fX(tidx, :, :) ;
        fY = facePathlines.fY(tidx, :, :) ;
        vX = vertexPathlines.vX(tidx, :, :) ;
        vY = vertexPathlines.vY(tidx, :, :) ;
        
        % Store x0,y0 PIV evaluation positions in pixels
        x0 = piv3d.x0 ;
        y0 = piv3d.y0 ;
        
        % 1. Interpolate v0_rs
        v0_rsX = reshape(piv3d.v0_rs(:, 1), size(x0)) ;
        v0_rsY = reshape(piv3d.v0_rs(:, 2), size(x0)) ;
        v0_rsZ = reshape(piv3d.v0_rs(:, 3), size(x0)) ;
        Fx = griddedInterpolant(x0', y0', v0_rsX', imethod, 'nearest') ;
        Fy = griddedInterpolant(x0', y0', v0_rsY', imethod, 'nearest') ;
        Fz = griddedInterpolant(x0', y0', v0_rsZ', imethod, 'nearest') ;
        % Query velocities
        v0_rs = [Fx(XX(:), YY(:)), Fy(XX(:), YY(:)), Fz(XX(:), YY(:))] ;
        
        % 2. Interpolate velocities onto face barycenter pathlines 
        % Query velocities
        v3dfaces_rs = [Fx(fX(:), fY(:)), Fy(fX(:), fY(:)), Fz(fX(:), fY(:))] ;

        % 3. Interpolate v0n_rs
        v0nrsN = reshape(piv3d.v0n_rs, size(x0)) ;
        Fn = griddedInterpolant(x0', y0', v0nrsN', imethod, 'nearest') ;
        % Query velocities
        v0n_rs = Fx(XX(:), YY(:)) ;

        % 4. Interpolate onto mesh vertices (v3dvertices)
        % Query velocities
        v3dvertices = [Fx(vX(:), vY(:)), Fy(vX(:), vY(:)), Fz(vX(:), vY(:))] ;

        % 5. Interpolate v0t2d
        v0t2dX = reshape(piv3d.v0t2d(:, 1), size(x0)) ;
        v0t2dY = reshape(piv3d.v0t2d(:, 2), size(x0)) ;
        Fx = griddedInterpolant(x0', y0', v0t2dX', imethod, 'nearest') ;
        Fy = griddedInterpolant(x0', y0', v0t2dY', imethod, 'nearest') ;
        % Query velocities
        v0t2d = [Fx(XX(:), YY(:)), Fy(XX(:), YY(:))] ;

        % 6. Interpolate v0t2dum
        v0t2dx = reshape(piv3d.v0t2d(:, 1) ./ piv3d.dilation, size(x0)) ;
        v0t2dy = reshape(piv3d.v0t2d(:, 2) ./ piv3d.dilation, size(x0)) ;
        Fx = griddedInterpolant(x0', y0', v0t2dx', imethod, 'nearest') ;
        Fy = griddedInterpolant(x0', y0', v0t2dy', imethod, 'nearest') ;
        % Query velocities
        v0t2dum = [Fx(XX(:), YY(:)), Fy(XX(:), YY(:))] ;

        %% BUILD ARRAYS
        vM(tidx, :, :) = v0_rs ;             % in um/dt rs at PIV evaluation points
        vfM(tidx, :, :) = v3dfaces_rs ;      % in um/dt rs at face barycenters
        vnM(tidx, :, :) = v0n_rs ;           % in um/dt rs at 
        vvM(tidx, :, :) = v3dvertices ;      % in um/min rs
        v2dM(tidx, :, :) = v0t2d ;           % in pixels/ min
        v2dMum(tidx, :, :) = v0t2dum ;       % in scaled pix/min, but proportional to um/min
        
    end
    
    clearvars first 
    disp('built v0 matrix')
    
    %% Save raw matrices
    pdir  = QS.dir.pathlines.data ;
    fvM   = fullfile(pdir, sprintf('vM_%04dt0.mat', t0)) ;
    fvfM  = fullfile(pdir, sprintf('vfM_%04dt0.mat', t0)) ;
    fvnM  = fullfile(pdir, sprintf('vnM_%04dt0.mat', t0)) ;
    fvvM  = fullfile(pdir, sprintf('vvM_%04dt0.mat', t0)) ;
    fv2dM = fullfile(pdir, sprintf('v2dM_%04dt0.mat', t0)) ;
    fv2dMum = fullfile(pdir, sprintf('v2dMum_%04dt0.mat', t0)) ;
    save(fvM, 'vM') 
    save(fvfM, 'vfM') 
    save(fvnM, 'vnM') 
    save(fvvM, 'vvM') 
    save(fv2dM, 'v2dM') 
    save(fv2dMum, 'v2dMum') 
    
    %% Filter in time axis -- light time smoothing tripulse of length 5
    % linfilt = 0.1 * ones(10, 1, 1) ;
    % ellipsoid = fspecial3('ellipsoid', [5, 1, 1]) ;
    disp('Building tripulse filter equivalent to tripuls()')
    tripulse3 = [ 0.3333; 0.6666; 1; 0.6666; 0.3333];
    tripulse3 = tripulse3 ./ sum(tripulse3(:)) ;
    tripulse3 = reshape(tripulse3, [length(tripulse3), 1]) ;
    
    % Check that no NaNs
    assert(~any(isnan(vM(:))))
    assert(~any(isnan(vvM(:))))
    assert(~any(isnan(vnM(:))))
    assert(~any(isnan(v2dM(:))))
    assert(~any(isnan(v2dMum(:))))
    assert(~any(isnan(vfM(:))))
    
    disp('filtering velocity matrices')
    vsmM = imfilter(vM, tripulse3, 'replicate');         % in um/min, rs
    vvsmM = imfilter(vvM, tripulse3, 'replicate');       % in um/min, rs
    vnsmM = imfilter(vnM, tripulse3, 'replicate');       % in um/min
    v2dsmM = imfilter(v2dM, tripulse3, 'replicate');     % in pix/min
    v2dsmMum = imfilter(v2dMum, tripulse3, 'replicate'); % in scaled pix/min, proportional to um/min  
    vfsmM = imfilter(vfM, tripulse3, 'replicate') ;      % in um/min, rs

    % Save the simpleminded averaging
    disp('Saving the time-smoothed velocities to disk')
    fvsmM   = QS.fileName.pathlines.v3dsm ;
    fvfsmM  = QS.fileName.pathlines.vfsm ;
    fvnsmM  = QS.fileName.pathlines.vnsm ;
    fvvsmM  = QS.fileName.pathlines.vvsm  ;
    fv2dsmM = QS.fileName.pathlines.v2dsm  ;
    fv2dsmMum = QS.fileName.pathlines.v2dsmum  ;
    save(fv2dsmMum, 'v2dsmMum') ;  % in scaled pix/min, proportional to um/min 
    save(fv2dsmM, 'v2dsmM') ;      % in pix/min
    save(fvnsmM, 'vnsmM') ;        % in um/min
    save(fvsmM, 'vsmM') ;          % in um/min
    save(fvvsmM, 'vvsmM') ;        % in um/min, rs
    save(fvfsmM, 'vfsmM') ;        % in um/min, rs
    
    disp('done')
else
    disp('Lagrangian-pathline-averaged velocities already found on disk.')
end

