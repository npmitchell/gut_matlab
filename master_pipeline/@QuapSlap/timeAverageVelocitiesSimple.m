function timeAverageVelocitiesSimple(QS, options)
%timeAverageVelocities(QS, options)
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%   timePoints : numeric 1D array
%       the timepoints to consider for the measurement. For ex, could
%       choose subset of the QS experiment timePoints
%
%
% Returns
% -------
%
%
% NPMitchell 2020

%% Default options
% Declare plotting options for limits
vtscale = 5 ;      % um / min
vnscale = 2 ;       % um / min
vscale = 2 ;        % um / min
alphaVal = 0.7 ;    % alpha for normal velocity heatmap
qsubsample = 5 ;   % quiver subsampling in pullback space 

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
else
    overwrite = false ;
end
if isfield(options, 'preview')
    preview = options.preview ;
else
    preview = false ;
end
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
else
    timePoints = QS.xp.fileMeta.timePoints ;
end
if isfield(options, 'pivimCoords')
    pivimCoords = options.pivimCoords ;
    if strcmp(pivimCoords(-1), 'e')
        doubleCovered = true ;
    else
        doubleCovered = false ;
    end
else
    pivimCoords = 'sp_sme' ;
    doubleCovered = true;
end
if isfield(options, 'vtscale')
    vtscale = options.vtscale ;
end
if isfield(options, 'vnscale')
    vnscale = options.vnscale ;
end
if isfield(options, 'vscale')
    vscale = options.vscale ;
end
if isfield(options, 'plot_vxyz')
    plot_vxyz = false ;
end

%% Unpack QS
pivDir = QS.dir.piv ;
piv3dfn = QS.fullFileBase.piv3d ;
ntps = length(timePoints) ;
[rot, ~] = QS.getRotTrans() ;
resolution = QS.APDV.resolution ; 
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
pivSimAvgDir = QS.dir.pivSimAvg ;

% Check if the time smoothed velocities exist already
v2dsmMumfn = QS.fileName.pivSimAvg.v2dMum ;
v2dsmMfn = QS.fileName.pivSimAvg.v2dM ;
vnsmMfn = QS.fileName.pivSimAvg.vnM ;
vsmMfn = QS.fileName.pivSimAvg.vM ;
% vertex-based velocities
vvsmMfn = QS.fileName.pivSimAvg.vvM ;
% face-based velocities
vfsmMfn = QS.fileName.pivSimAvg.vfM ;
if ~exist(v2dsmMumfn, 'file') && exist(v2dsmMfn, 'file') && ...
        exist(vnsmMfn, 'file') && exist(vsmMfn, 'file') && ...
        exist(vfsmMfn, 'file') && exist(vvsmMfn, 'file') && ...
        ~overwrite
    
    disp('Could not find time-smoothed velocities on disk')
    disp('Computing them...')
    first = true ;
    ntps = length(timePoints)-1;
    for i = 1:ntps
        tp = timePoints(i) ;
        dt = timePoints(i + 1) - tp ;
        
        disp(['Filling in velocity matrices, t=' num2str(tp)])
        QS.setTime(tp) ;
        QS.getCurrentVelocity('piv3d') ;
        piv3d = QS.currentVelocity.piv3d ;
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
        vM(i, :, :) = piv3d.v0_rs ;          % in um/min rs
        % try
        vfM(i, :, :) = piv3d.v3dfaces_rs ;   % in um/min rs
        % catch
        %     v3dfaces = piv3d.v3dfaces ;
        %     vfM(i, :, :) = (rot * v3dfaces')' * resolution / dt ; 
        % end
        vnM(i, :, :) = piv3d.v0n_rs ;        % in rs coords, unit length
        vvM(i, :, :) = (rot * piv3d.v3dvertices')' * resolution ;
        v2dM(i, :, :) = piv3d.v0t2d ;        % in pixels/ min
        v2dMum(i, :, 1) = piv3d.v0t2d(:, 1) ./ piv3d.dilation ; % in pix/min, scaled as um/min
        v2dMum(i, :, 2) = piv3d.v0t2d(:, 2) ./ piv3d.dilation ;
        v2dMum(i, :, 2) = piv3d.v0t2d(:, 2) ./ piv3d.dilation ;
    end
    clearvars first 
    disp('built v0 matrix')
    % Filter in time axis
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
    save(fullfile(pivSimAvgDir, 'v2dMum_simpletimeavg.mat'), 'v2dsmMum') ;  % in scaled pix/min, proportional to um/min 
    save(fullfile(pivSimAvgDir, 'v2dM_simpletimeavg.mat'), 'v2dsmM') ;      % in pix/min
    save(fullfile(pivSimAvgDir, 'vnM_simpletimeavg.mat'), 'vnsmM') ;        % in um/min
    save(fullfile(pivSimAvgDir, 'vM_simpletimeavg.mat'), 'vsmM') ;          % in um/min
    save(fullfile(pivSimAvgDir, 'vvM_simpletimeavg.mat'), 'vvsmM') ;          % in um/min, rs
    save(fullfile(pivSimAvgDir, 'vfM_simpletimeavg.mat'), 'vfsmM') ;        % in um/min, rs

end
disp('done with timeAverageVelocitiesSimple()')