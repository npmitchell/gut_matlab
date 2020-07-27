function timeAverageVelocities(QS, options)
%measurePullbackStreamlines(QS, options)
%   Use pathlines of optical flow in pullback space to query velocities
%   and average along Lagrangian pathlines. 
%   Default is to weight velocities contributions with a time width 5 
%   weighted by tripulse filter.
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
pivimCoords = QS.piv.imCoords ;
samplingResolution = '1x' ;
imethod = 'linear' ;  % interpolation method for velocities onto pathlines
twidth = 2 ;          % average over (t-twidth, t+twidth) timepoints

%% Unpack options
% Default values for options are to use sphi smoothed extended coords
% as PIV reference coord sys
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
if doubleResolution
    fileNames = QS.fileName.pivAvg2x ;
else
    fileNames = QS.fileName.pivAvg ;
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
        overwrite
    
    disp('Could not find time-smoothed velocities on disk')
    disp('Computing them...')
                
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
            vsmM = zeros(ntps, size(piv3d.v0_rs, 1), size(piv3d.v0_rs, 2));
            vfsmM = zeros(ntps, size(piv3d.v3dfaces, 1), size(piv3d.v3dfaces, 2)); 
            vvsmM = zeros(ntps, ...
                size(piv3d.v3dvertices, 1), ...
                size(piv3d.v3dvertices, 2)); 
            vnsmM = zeros(ntps, size(piv3d.v0n_rs, 1), size(piv3d.v0n_rs, 2));
            v2dsmM = zeros(ntps, size(piv3d.v0t2d, 1), size(piv3d.v0t2d, 2));
            v2dsmMum = zeros(ntps, size(piv3d.v0t2d, 1), size(piv3d.v0t2d, 2));
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
        
        % Compute short lagrangian pathlines (t-twidth, t+twidth)
        % Note that we clip the timepoints to compute at 1 and ntps.
        disp(['Computing pathlines around t=', num2str(tp)])
        tp2do = (tp - twidth):(tp + twidth) ;
        tp2do = min(max(timePoints(1), tp2do), timePoints(end-1)) ;
        popts.timePoints = tp2do ;
        % PIV pathlines
        [XXpath, YYpath] = QS.pullbackPathlines(piv3d.x0, piv3d.y0, tp, popts) ;
        % face pathlines
        if strcmp(QS.piv.imCoords, 'sp_sme')
            im0 = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
            mesh0 = load(sprintf(QS.fullFileBase.spcutMeshSm, tp), 'spcutMeshSm') ;
            mesh0 = mesh0.spcutMeshSm ;
            umax = max(mesh0.u(:, 1)) ;
            vmax = max(mesh0.u(:, 2)) ;
            mXY = QS.uv2XY(im0, mesh0.u, doubleCovered, umax, vmax) ;
        end
        bc = barycenter(mXY, mesh0.f) ; % this is in gptoolbox
        [fXpath, fYpath] = QS.pullbackPathlines(bc(:, 1), bc(:, 2), tp, popts) ;
        % vertex pathlines
        [vXpath, vYpath] = QS.pullbackPathlines(mXY(:, 1), mXY(:, 2), tp, popts) ;
        
        % Check eval points
        clf
        plot(mXY(:, 1), mXY(:, 2), '.')
            
        %% Run through all timepoints in the window and compute a grid of
        % velocities to average over -- Grab PIV for each timepoint and 
        % interpolate onto pathline segment for all timepoints in tp2do]
        
        %% Filter in time: multiply each contribution by weight of tripulse
        % linfilt = 0.1 * ones(10, 1, 1) ;
        % ellipsoid = fspecial3('ellipsoid', [5, 1, 1]) ;
        disp('Building tripulse filter equivalent to tripuls()')
        if twidth == 2
            tripulse = [ 0.3333; 0.6666; 1; 0.6666; 0.3333];
            tripulse = tripulse ./ sum(tripulse(:)) ;
            % tripulse = reshape(tripulse, [length(tripulse), 1]) ;
        else
            error(['build tripulse of twidth ' num2str(twidth) ' here'])
        end
        
        % Preallocate
        vMtp  = zeros(size(vsmM, 2), size(vsmM, 3)) ;    % in um/dt rs at PIV evaluation points
        vfMtp = zeros(size(vfsmM, 2), size(vfsmM, 3)) ;  % in um/dt rs at face barycenters
        vnMtp = zeros(size(vnsmM, 2), size(vnsmM, 3)) ;  % in um/dt rs at 
        vvMtp = zeros(size(vvsmM, 2), size(vvsmM, 3)) ;  % in um/min rs
        v2dMtp   = zeros(size(v2dsmM, 2), size(v2dsmM, 3)) ;  % in pixels/ min
        v2dMumtp = zeros(size(v2dsmMum, 2), size(v2dsmMum, 3)) ;  % in scaled pix/min, but proportional to um/min
        for qq = 1:length(tp2do)
            disp(['timeAverageVelocities: Interpolating t0=', ...
                num2str(tp) ' coords onto pathline at t=', ...
                num2str(tp2do(qq))])
            % Load streamline positions at this timePoint
            XX = XXpath(qq, :, :) ;
            YY = YYpath(qq, :, :) ;
            fX = fXpath(qq, :, :) ;
            fY = fYpath(qq, :, :) ;
            vX = vXpath(qq, :, :) ;
            vY = vYpath(qq, :, :) ;
            
            % Grab PIV for this timepoint and interpolate onto pathline
            % segment for all timepoints in tp2do          
            QS.setTime(tp2do(qq)) ;
            
            % Load piv if this is not the last timePoint
            if doubleResolution 
                QS.getCurrentVelocity('piv3d2x') ;
                piv3d = QS.currentVelocity.piv3d2x ;
            else
                QS.getCurrentVelocity('piv3d') ;
                piv3d = QS.currentVelocity.piv3d ;
            end

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

            %% BUILD ARRAYS -- AVERAGE ARRAY over time dimension w/ weight
            vMtp = vMtp + tripulse(qq) * v0_rs ;           % in um/dt rs at PIV evaluation points
            vfMtp = vfMtp + tripulse(qq) * v3dfaces_rs ;   % in um/dt rs at face barycenters
            vnMtp = vnMtp + tripulse(qq) * v0n_rs ;        % in um/dt rs at 
            vvMtp = vvMtp + tripulse(qq) * v3dvertices ;   % in um/min rs
            v2dMtp = v2dMtp + tripulse(qq) * v0t2d ;       % in pixels/ min
            v2dMumtp = v2dMumtp + tripulse(qq) * v0t2dum ; % in scaled pix/min, but proportional to um/min
        end
        
        %% BUILD ARRAYS by collating timepoints
        vsmM(tidx, :, :) = v0_rs ;             % in um/dt rs at PIV evaluation points
        vfsmM(tidx, :, :) = v3dfaces_rs ;      % in um/dt rs at face barycenters
        vnsmM(tidx, :, :) = v0n_rs ;           % in um/dt rs at 
        vvsmM(tidx, :, :) = v3dvertices ;      % in um/min rs
        v2dsmM(tidx, :, :) = v0t2d ;           % in pixels/ min
        v2dsmMum(tidx, :, :) = v0t2dum ;       % in scaled pix/min, but proportional to um/min
        
    end
    
    clearvars first 
    disp('built v0 matrix')
    
    % Check that no NaNs
    assert(~any(isnan(vsmM(:))))
    assert(~any(isnan(vvsmM(:))))
    assert(~any(isnan(vnsmM(:))))
    assert(~any(isnan(v2dsmM(:))))
    assert(~any(isnan(v2dsmMum(:))))
    assert(~any(isnan(vfsmM(:))))
    
    %% Save raw matrices
    save(fileNames.v3d, 'vsmM') 
    save(fileNames.vf, 'vfsmM') 
    save(fileNames.vn, 'vnsmM') 
    save(fileNames.vv, 'vvsmM') 
    save(fileNames.v2d, 'v2dsmM') 
    save(fileNames.v2dum, 'v2dsmMum') 
    
    disp('done')

end