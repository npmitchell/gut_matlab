function measurePullbackPathlines(QS, options)
%measurePullbackPathlines(QS, options)
%   Measure pathlines of optical flow in pullback space. 
%   These pathlines can then be used to query velocities or other
%   properties that are parameterized by pullback coordinates. 
%   For example, we average velocites along Lagrangian pathlines.
%   For another example, to build a Lagrangian-frame measure of divergence 
%   of tangential velocity, div(v_t), we can query div(v_t) at the coords
%   of the PullbackPathlines (via interpolation of div(v_t) defined on
%   pullback mesh vertices).
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
% Saves to disk
% -------------
% sprintf(QS.fileName.pathlines.XY, t0) ;
% sprintf(QS.fileName.pathlines.vXY, t0) ;
% sprintf(QS.fileName.pathlines.fXY, t0) ;
% sprintf(QS.fileName.pathlines.XYZ, t0) ;
%
% See also
% --------
% timeAverageVelocities(QS, options)
% pullbackPathlines(QS, options)
%
% NPMitchell 2020

%% Default options
overwrite = false ;
preview = false ;
debug = false ;
pivimCoords = QS.piv.imCoords ;
timePoints = QS.xp.fileMeta.timePoints ;
samplingResolution = '1x' ;
nY2plot = 30 ;
scatterSz = 2 ;
movieScatterSz = 2 ;
if nargin < 2
    options = struct() ;
end

%% Unpack options
% Default values for options are to use sphi smoothed extended coords
% as PIV reference coord sys
if isfield(options, 't0Pathlines')
    t0 = options.t0Pathlines ;
    disp(['Setting t0 for Pathlines to be: ' num2str(t0)])
elseif isfield(options, 't0')
    t0 = options.t0 ;
    disp(['Setting t0 for Pathlines to be: ' num2str(t0)])
else
    disp('Using default t0 for pathlines')
    t0 = QS.t0set() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'debug')
    debug = options.debug ;
end
if isfield(options, 'nY2plot')
    nY2plot = options.nY2plot ;
end
if strcmp(pivimCoords(end), 'e')
    doubleCovered = true ;
    Yoffset = 0.25 ;
else
    doubleCovered = false ;
    Yoffset = 0.0 ;
end

%% Unpack QS
% [rot, ~] = QS.getRotTrans() ;
% resolution = QS.APDV.resolution ; 
[~, ~, ~, xyzlim] = QS.getXYZLims() ;
axis_order = QS.data.axisOrder ;
blue = QS.plotting.colors(1, :) ;
red = QS.plotting.colors(2, :) ;
green = QS.plotting.colors(4, :) ;

% Define t0 at which pathlines form grid
tIdx0 = QS.xp.tIdx(t0) ;

%% Create directory for pathlines
pathlineDir = QS.dir.pathlines ;
pdir = sprintf(pathlineDir.data, t0) ;
XYdir = sprintf(pathlineDir.XY, t0) ;
XYZdir = sprintf(pathlineDir.XYZ, t0) ;
vXYdir = sprintf(pathlineDir.vXY, t0) ;
fXYdir = sprintf(pathlineDir.fXY, t0) ;
v3ddir = sprintf(pathlineDir.v3d, t0) ;
f3ddir = sprintf(pathlineDir.f3d, t0) ;
dirs2make = {pdir, XYdir, XYZdir, vXYdir, fXYdir, v3ddir, f3ddir} ;
for qq = 1:length(dirs2make)
    dir2make = dirs2make{qq} ;
    if ~exist(dir2make, 'dir')
        disp(['Making dir: ' dir2make])
        mkdir(dir2make)
    end
end

%% Perform/Load pathline calculations

% Check if the time smoothed velocities exist already
% 2d velocities (pulled back), scaled by dilation of metric are v2dum 
% 2D velocities (pulled back) are v2d
% normal velocities on fieldfaces are vn
% 3d velocities on fieldfaces are v3d
% vertex-based velocities are vv
% face-based velocities are vf

QS.clearTime() ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIV evaluation coordinates in XY pixel space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plineXY = sprintf(QS.fileName.pathlines.XY, t0) ;
if ~exist(plineXY, 'file') || overwrite
    disp('Computing piv Pathlines for XY(t0)')
    % Load 'initial' positions for pathlines to intersect at t=t0
    QS.getPIV() 
    piv = QS.piv.raw ;
    x0 = piv.x{QS.xp.tIdx(t0)} ;
    y0 = piv.y{QS.xp.tIdx(t0)} ;
    
    % Create pathlines emanating from (x0,y0) at t=t0
    [XX, YY] = QS.pullbackPathlines(x0, y0, t0, options) ;
    Lx = QS.piv.Lx ;
    Ly = QS.piv.Ly ;
    
    % Save pathlines of PIV evaluation points XY
    pivPathlines = struct() ;
    pivPathlines.XX = XX ;
    pivPathlines.YY = YY ;
    pivPathlines.Lx = Lx ;
    pivPathlines.Ly = Ly ;
    pivPathlines.t0 = t0 ;
    pivPathlines.tIdx0 =  tIdx0 ;
    save(plineXY, 'pivPathlines')
    
    computed_XY = true ;
else
    disp(['PIV XY pathlines already on disk: ' plineXY])
    computed_XY = false ;
end

% Save image of result
plineFig = [plineXY(1:end-4) '.png'] ;
if ~exist(plineFig, 'file') || overwrite
    disp(['Saving PIV XY pathline image to disk: ' plineFig])
    if ~computed_XY 
        load(plineXY, 'pivPathlines')
        XX = pivPathlines.XX ;
        YY = pivPathlines.YY ;
        t0 = pivPathlines.t0 ;
        tIdx0 = pivPathlines.tIdx0 ;
        Lx = pivPathlines.Lx ;
        Ly = pivPathlines.Ly ;
        computed_XY = true ;
    end
    
    % Plot the pathlines -- NOTE that we flip YY coordinates (MATLAB)
    qsubX = round(size(XX, 2) / nY2plot * QS.a_fixed) ;
    if doubleCovered
        qsubY = round(size(XX, 3) / nY2plot * 2);
    else
        qsubY = round(size(XX, 3) / nY2plot) ;
    end
    colors = parula(length(timePoints)) ;
    clf
    for qq = fliplr(1:length(timePoints))
        xx = XX(qq, 1:qsubX:end, 1:qsubY:end) ;
        yy = YY(qq, 1:qsubX:end, 1:qsubY:end) ;
        scatter(xx(:) / double(Lx(qq)), ...
            (1 - (yy(:) / double(Ly(qq))) - Yoffset) / (1 - 2*Yoffset), ...
            scatterSz, ...
            'markerfacecolor', colors(qq, :), ...
            'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.3) ;
        hold on;
    end
    axis equal
    ylim([0, 1])
    xlim([0, 1])
    daspect([1 QS.a_fixed 1])
    title('pullback pathlines: PIV coordinates', 'Interpreter', 'Latex')
    ylabel('circumferential position, $\phi / 2\pi$', ...
        'Interpreter', 'Latex')
    xlabel('ap position, $\zeta$', 'Interpreter', 'Latex')
    cb = colorbar() ;
    caxis([(timePoints(1) - t0) * QS.timeInterval, ...
           (timePoints(end) - t0) * QS.timeInterval]) ;
    ylabel(cb, ['time [' QS.timeUnits ']'], ...
        'Interpreter', 'Latex')
    saveas(gcf, plineFig)
else
    disp(['PIV XY pathline image already on disk: ' plineFig])
end

% Make movie
for tidx = 1:length(timePoints)    
    if mod(tidx, 10) == 0
        disp(['t = ', num2str(timePoints(tidx))])
    end
    fn = fullfile(sprintf(pathlineDir.XY, t0), ...
        [sprintf(QS.fileBase.name, timePoints(tidx)) '.png']) ;
    
    if ~exist(fn, 'file') || overwrite
        % Load pathlines if not already in RAM
        if ~computed_XY 
            load(plineXY, 'pivPathlines')
            XX = pivPathlines.XX ;
            YY = pivPathlines.YY ;
            t0 = pivPathlines.t0 ;
            tIdx0 = pivPathlines.tIdx0 ;
            Lx = pivPathlines.Lx ;
            Ly = pivPathlines.Ly ;
            computed_XY = true ;
        end
    
        clf ;
        movieScatterSz = 3 ;
        scatter(XX(tidx, :) / double(Lx(tidx)), ...
            (1 - (YY(tidx, :) / double(Ly(tidx))) - Yoffset) / (1 - 2*Yoffset), ...
            movieScatterSz, 'markerfacecolor', QS.plotting.colors(1, :), ...
            'markeredgecolor', 'none', 'MarkerFaceAlpha', 1.0) ;
        axis equal
        ylim([0, 1])
        xlim([0, 1])
        daspect([1 QS.a_fixed 1])
        time_in_units = (timePoints(tidx) - t0)* QS.timeInterval ;
        tstr = ['$t = $' num2str(time_in_units) ' ' QS.timeUnits] ;
        title(['pullback pathlines, ' tstr], 'Interpreter', 'Latex')
        ylabel('circumferential position, $\phi / 2\pi$', ...
            'Interpreter', 'Latex')
        xlabel('ap position, $\zeta$', 'Interpreter', 'Latex')
        saveas(gcf, fn)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh vertex coordinates in XY pixel space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plinevXY = sprintf(QS.fileName.pathlines.vXY, t0) ;
if ~exist(plinevXY, 'file') || overwrite 
    % Create pathlines emanating from vertex positions (vX,vY) at t=t0
    % Load pullback mesh vertex positions to get Lx and Ly (will be done
    % later in pullbackPathlines() if not now, so no extra cost to do it
    % now instead). 
    QS.getPIV() 
    Lx = QS.piv.Lx ;
    Ly = QS.piv.Ly ;
    
    if strcmp(QS.piv.imCoords, 'sp_sme')
        mesh0 = load(sprintf(QS.fullFileBase.spcutMeshSm, t0), 'spcutMeshSm') ;
    else
        error('handle this coord sys here')
    end
    mesh0 = mesh0.spcutMeshSm ;
    umax0 = max(mesh0.u(:, 1)) ;
    vmax0 = max(mesh0.u(:, 2)) ;
    m0XY = QS.uv2XY([Lx(tIdx0), Ly(tIdx0)], mesh0.u, doubleCovered, ...
        umax0, vmax0) ;
    m0X = m0XY(:, 1) ;
    m0Y = m0XY(:, 2) ;
    
    % Create pathlines emanating from vertex positions (vX,vY) at t=t0
    [vX, vY] = QS.pullbackPathlines(m0X, m0Y, t0, options) ;
    vX = reshape(vX, [length(timePoints), QS.nU, QS.nV]) ;
    vY = reshape(vY, [length(timePoints), QS.nU, QS.nV]) ;
    
    % Save pathlines of mesh locations u in XY coords, vXY
    vertexPathlines = struct() ;
    vertexPathlines.vX = vX ;
    vertexPathlines.vY = vY ;
    vertexPathlines.t0 = t0 ;
    vertexPathlines.tIdx0 = tIdx0 ;
    vertexPathlines.Lx = Lx ;
    vertexPathlines.Ly = Ly ;
    save(plinevXY, 'vertexPathlines')
    
    computed_vXY = true ;
else
    disp(['Mesh vertex XY pathlines already on disk: ' plinevXY])
    computed_vXY = false ;
end

% Save image of result
plineFig = [plinevXY(1:end-4) '.png'] ;
if ~exist([plinevXY(1:end-4) '.png'], 'file') || overwrite
    disp(['Saving mesh vertex XY pathline image to disk: ' plineFig])
    if ~computed_vXY 
        load(QS.fileName.pathlines.vXY, 'vertexPathlines')
        vX = vertexPathlines.vX ;
        vY = vertexPathlines.vY ;
        t0 = vertexPathlines.t0 ;
        tIdx0 = vertexPathlines.tIdx ;
        Lx = vertexPathlines.Lx ;
        Ly = vertexPathlines.Ly ;
    end
    
    % Plot the pathlines -- NOTE that we flip YY coordinates (MATLAB)
    qsubX = round(size(vX, 2) / nY2plot) ;
    qsubY = round(size(vX, 3) / nY2plot * QS.a_fixed) ;
    colors = parula(length(timePoints)) ;
    clf
    for qq = fliplr(1:length(timePoints))
        xx = vX(qq, 1:qsubX:end, 1:qsubY:end) ;
        yy = vY(qq, 1:qsubX:end, 1:qsubY:end) ;
        scatter(xx(:) / double(Lx(qq)), ...
            (1 - (yy(:) / double(Ly(qq))) - Yoffset) / (1 - 2*Yoffset), ...
            scatterSz, ...
            'markerfacecolor', colors(qq, :), ...
            'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.3) ;
        hold on;
        pause(0.1) ;
    end
    axis equal
    ylim([0, 1])
    xlim([0, 1])
    daspect([1 QS.a_fixed 1])
    title('pullback pathlines: mesh vertices', 'Interpreter', 'Latex')
    ylabel('circumferential position, $\phi / 2\pi$', ...
        'Interpreter', 'Latex')
    xlabel('ap position, $\zeta$', 'Interpreter', 'Latex')
    cb = colorbar() ;
    caxis([(timePoints(1) - t0) * QS.timeInterval, ...
           (timePoints(end) - t0) * QS.timeInterval]) ;
    ylabel(cb, ['time [' QS.timeUnits ']'], ...
        'Interpreter', 'Latex')
    saveas(gcf, plineFig)
else
    disp(['Mesh vertex XY pathline image already on disk: ' plineFig])
end

% Make movie
for tidx = 1:length(timePoints)    
    if mod(tidx, 10) == 0
        disp(['t = ', num2str(timePoints(tidx))])
    end
    fn = fullfile(sprintf(pathlineDir.vXY, t0), ...
        [sprintf(QS.fileBase.name, timePoints(tidx)) '.png']) ;
    
    if ~exist(fn, 'file') || overwrite
        % Load pathlines if not already in RAM
        if ~computed_vXY 
            load(plineXY, 'pivPathlines')
            vX = pivPathlines.vX ;
            vY = pivPathlines.vY ;
            t0 = pivPathlines.t0 ;
            tIdx0 = pivPathlines.tIdx0 ;
            Lx = pivPathlines.Lx ;
            Ly = pivPathlines.Ly ;
            computed_vXY = true ;
        end
    
        clf ;
        movieScatterSz = 3 ;
        scatter(vX(tidx, :) / double(Lx(tidx)), ...
            (1 - (vY(tidx, :) / double(Ly(tidx))) - Yoffset) / (1 - 2*Yoffset), ...
            movieScatterSz, 'markerfacecolor', QS.plotting.colors(1, :), ...
            'markeredgecolor', 'none', 'MarkerFaceAlpha', 1.0) ;
        axis equal
        ylim([0, 1])
        xlim([0, 1])
        daspect([1 QS.a_fixed 1])
        time_in_units = (timePoints(tidx) - t0)* QS.timeInterval ;
        tstr = ['$t = $' num2str(time_in_units) ' ' QS.timeUnits] ;
        title(['pullback pathlines, ' tstr], 'Interpreter', 'Latex')
        ylabel('circumferential position, $\phi / 2\pi$', ...
            'Interpreter', 'Latex')
        xlabel('ap position, $\zeta$', 'Interpreter', 'Latex')
        saveas(gcf, fn)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh face barycenters in XY pixel space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plineFXY = sprintf(QS.fileName.pathlines.fXY, t0) ;
if ~exist(plineFXY, 'file') || overwrite
    % Create pathlines emanating from barycenters (bcX,bcY) at t=t0
    % Load pullback mesh vertex positions to get Lx and Ly (will be done
    % later in pullbackPathlines() if not now, so no extra cost to do it
    % now instead).
    QS.loadPIV() 
    Lx = QS.piv.Lx ;
    Ly = QS.piv.Ly ;
    
    if strcmp(QS.piv.imCoords, 'sp_sme')
        mesh0 = load(sprintf(QS.fullFileBase.spcutMeshSm, t0), 'spcutMeshSm') ;
    else
        error('handle coord sys here')
    end
    mesh0 = mesh0.spcutMeshSm ;
    umax0 = max(mesh0.u(:, 1)) ;
    vmax0 = max(mesh0.u(:, 2)) ;
    m0XY = QS.uv2XY([Lx(tIdx0), Ly(tIdx0)], mesh0.u, doubleCovered, umax0, vmax0) ;
    
    % Obtain barycenters of faces
    f0XY = barycenter(m0XY, mesh0.f) ;  % this is in gptoolbox
    f0X = f0XY(:, 1) ;
    f0Y = f0XY(:, 2) ;
    
    % Create pathlines emanating from barycenters (bcX,bcY) at t=t0
    [fX, fY] = QS.pullbackPathlines(f0X, f0Y, t0, options) ;
    
    % Save pathlines of mesh face barycenters f in XY coords, fXY
    facePathlines = struct() ;
    facePathlines.fX = fX ;
    facePathlines.fY = fY ;
    facePathlines.t0 = t0 ;
    facePathlines.tIdx0 = tIdx0 ;
    facePathlines.Lx = Lx ;
    facePathlines.Ly = Ly ;
    save(plineFXY, 'facePathlines')
    
    computed_fXY = true ;
else
    computed_fXY = false ;
end

% Save image of result
plineFig = [plineFXY(1:end-4) '.png'] ;
if ~exist(plineFig, 'file') || overwrite
    if ~computed_fXY 
        load(QS.fileName.pathlines.fXY, 'facePathlines')
        fX = facePathlines.fX ;
        fY = facePathlines.fY ;
        t0 = facePathlines.t0 ;
        tIdx0 = facePathlines.tIdx0  ;
    end
    
    % Plot the pathlines -- NOTE that we flip YY coordinates (MATLAB)
    qsubX = round(sqrt(size(fX, 2)) / nY2plot * 8) ;
    colors = parula(length(timePoints)) ;
    clf
    for qq = fliplr(1:length(timePoints))
        xx = fX(qq, 1:qsubX:end) ;
        yy = fY(qq, 1:qsubX:end) ;
        scatter(xx(:) / double(Lx(qq)), ...
            (1 - (yy(:) / double(Ly(qq))) - Yoffset) / (1 - 2*Yoffset), ...
            scatterSz, ...
            'markerfacecolor', colors(qq, :), ...
            'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.3) ;
        hold on;
    end
    axis equal
    ylim([0, 1])
    xlim([0, 1])
    daspect([1 QS.a_fixed 1])
    title('pullback pathlines: mesh vertices', 'Interpreter', 'Latex')
    ylabel('circumferential position, $\phi / 2\pi$', ...
        'Interpreter', 'Latex')
    xlabel('ap position, $\zeta$', 'Interpreter', 'Latex')
    cb = colorbar() ;
    caxis([(timePoints(1) - t0) * QS.timeInterval, ...
           (timePoints(end) - t0) * QS.timeInterval]) ;
    ylabel(cb, ['time [' QS.timeUnits ']'], ...
        'Interpreter', 'Latex')
    saveas(gcf, plineFig)
else
    disp(['Mesh face barycenter XY pathline image already on disk: ' plineFig])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Push forward from pullback to embedding coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plineXYZ = sprintf(QS.fileName.pathlines.XYZ, t0) ;
if ~exist(plineXYZ, 'file') || overwrite
    disp('Compute 3d pathlines in embedding space for PIV eval coords')
    % Build 3d pathlines in embedding space for all vertex locations
    % Load 2d pathlines in pullback space if not already in RAM
    if ~computed_XY 
        load(plineXY, 'pivPathlines')
        XX = pivPathlines.XX ;
        YY = pivPathlines.YY ;
        t0 = pivPathlines.t0 ;
        tIdx0 = pivPathlines.tIdx0 ;
        Lx = pivPathlines.Lx ;
        Ly = pivPathlines.Ly ;
        computed_XY = true ;
    end
    
    % Preallocate v3d = (XX, YY, ZZ)
    X3 = zeros(size(XX)) ;
    Y3 = zeros(size(XX)) ;
    Z3 = zeros(size(XX)) ;
    fieldfacesCell = cell(length(timePoints), 1) ;
    
    % For each timepoint, push forward into embedding space
    for tidx = 1:length(timePoints)
        tp = timePoints(tidx) ;
        if mod(tidx, 10) == 0
            disp(['t = ', num2str(tp)])
        end
        
        % Load this timepoint's cutMesh
        mesh0 = load(sprintf(QS.fullFileBase.spcutMeshSm, tp), ...
                             'spcutMeshSm') ;
        mesh0 = mesh0.spcutMeshSm ;
        umax0 = max(mesh0.u(:, 1)) ;
        vmax0 = max(mesh0.u(:, 2)) ;
        tileCount = [2 2] ;
        [tm0f, tm0v2d, tm0v3d, ~] = tileAnnularCutMesh(mesh0, tileCount);
        tm0XY = QS.uv2XY([Lx(tidx) Ly(tidx)], tm0v2d, ...
                         doubleCovered, umax0, vmax0) ;
        Xeval = XX(tidx, :, :) ;
        Yeval = YY(tidx, :, :) ;
        [pt0, fieldfaces] = interpolate2Dpts_3Dmesh(tm0f, tm0XY,...
                                        tm0v3d, [Xeval(:), Yeval(:)]) ;
        % Store in grid 
        X3(tidx, :, :) = reshape(pt0(:, 1), [size(X3, 2), size(X3, 3)]) ;
        Y3(tidx, :, :) = reshape(pt0(:, 2), [size(Y3, 2), size(Y3, 3)]) ;
        Z3(tidx, :, :) = reshape(pt0(:, 3), [size(Z3, 2), size(Z3, 3)]) ;
        fieldfacesCell{tidx} = fieldfaces ;
    end
    
    % Save pathlines of PIV evaluation coords in embedding coords, XYZ
    piv3dPathlines = struct() ;
    piv3dPathlines.XX = X3 ;
    piv3dPathlines.YY = Y3 ;
    piv3dPathlines.ZZ = Z3 ;
    piv3dPathlines.t0 = t0 ;
    piv3dPathlines.tIdx0 = tIdx0 ;
    save(plineXYZ, 'piv3dPathlines')
end                             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Push forward vertex pathlines from pullback to embedding coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plinev3d = sprintf(QS.fileName.pathlines.v3d, t0) ;
if ~exist(plinev3d, 'file') || overwrite || true
    disp('Compute 3d pathlines in embedding space for Lagrangian/advected spcutMesh vertex coords')
    % Build 3d pathlines in embedding space for all vertex locations
    % Load 2d pathlines in pullback space if not already in RAM
    if ~computed_vXY 
        load(sprintf(QS.fileName.pathlines.vXY, t0), 'vertexPathlines')
        vX = vertexPathlines.vX ;
        vY = vertexPathlines.vY ;
        t0 = vertexPathlines.t0 ;
        tIdx0 = vertexPathlines.tIdx0 ;
        Lx = vertexPathlines.Lx ;
        Ly = vertexPathlines.Ly ;
        computed_vXY = true ;
    end
    
    % Preallocate v3d = (XX, YY, ZZ)
    vX3 = zeros(size(vX)) ;
    vY3 = zeros(size(vX)) ;
    vZ3 = zeros(size(vX)) ;
    vX3rs = zeros(size(vX)) ;
    vY3rs = zeros(size(vX)) ;
    vZ3rs = zeros(size(vX)) ;
    fieldfacesCell = cell(length(timePoints), 1) ;
    
    % For each timepoint, push advected vertices forward into embedding space
    for tidx = 1:length(timePoints)
        tp = timePoints(tidx) ;
        if mod(tidx, 10) == 0
            disp(['t = ', num2str(tp)])
        end
        
        % Load this timepoint's cutMesh
        mesh0 = load(sprintf(QS.fullFileBase.spcutMeshSm, tp), ...
                             'spcutMeshSm') ;
        mesh0 = mesh0.spcutMeshSm ;
        umax0 = max(mesh0.u(:, 1)) ;
        vmax0 = max(mesh0.u(:, 2)) ;
        tileCount = [2 2] ;
        [tm0f, tm0v2d, tm0v3d, ~] = tileAnnularCutMesh(mesh0, tileCount);
        tm0XY = QS.uv2XY([Lx(tidx) Ly(tidx)], tm0v2d, ...
                         doubleCovered, umax0, vmax0) ;
                     
        % Evaluate at advected Lagrangian-frame vertex positions
        Xeval = vX(tidx, :, :) ;
        Yeval = vY(tidx, :, :) ;
        
        % Adjust the XY positions so that none are out of bounds
        epsilon = 1e-7 ;
        minXallowed = 1 + epsilon ;
        Xeval = Xeval(:) ;
        Yeval = Yeval(:) ;
        Xeval(Xeval(:, 1) < minXallowed) = minXallowed ;
        Xeval(Xeval(:, 1) > Lx(tidx) - epsilon) = Lx(tidx) - epsilon ;
        
        [pt0, fieldfaces] = interpolate2Dpts_3Dmesh(tm0f, tm0XY,...
                                        tm0v3d, [Xeval(:), Yeval(:)]) ;
        assert(all(isfinite(pt0(:)))) 
        
        % Check that points are inside
        % plot(tm0XY(:, 1), tm0XY(:, 2), '-')
        % hold on ;
        % plot(Xeval(:), Yeval(:), '.')
        
        % Rotate and scale    
        vXYZrs = QS.xyz2APDV(pt0) ;
        
        % Store in grid 
        vX3(tidx, :, :) = reshape(pt0(:, 1), [size(vX3, 2), size(vX3, 3)]) ;
        vY3(tidx, :, :) = reshape(pt0(:, 2), [size(vY3, 2), size(vY3, 3)]) ;
        vZ3(tidx, :, :) = reshape(pt0(:, 3), [size(vZ3, 2), size(vZ3, 3)]) ;
        vX3rs(tidx, :, :) = reshape(vXYZrs(:, 1), [size(vX3, 2), size(vX3, 3)]) ;
        vY3rs(tidx, :, :) = reshape(vXYZrs(:, 2), [size(vY3, 2), size(vY3, 3)]) ;
        vZ3rs(tidx, :, :) = reshape(vXYZrs(:, 3), [size(vZ3, 2), size(vZ3, 3)]) ;
        fieldfacesCell{tidx} = fieldfaces ;
    end
    
    % Save pathlines of PIV evaluation coords in embedding coords, XYZ
    v3dPathlines = struct() ;
    v3dPathlines.vX = vX3 ;
    v3dPathlines.vY = vY3 ;
    v3dPathlines.vZ = vZ3 ;
    v3dPathlines.vXrs = vX3rs ;
    v3dPathlines.vYrs = vY3rs ;
    v3dPathlines.vZrs = vZ3rs ;
    v3dPathlines.t0 = t0 ;
    v3dPathlines.tIdx0 = tIdx0 ;
    save(plinev3d, 'v3dPathlines')
    computed_v3d = true ;
else
    computed_v3d = false ;
end

%% Plot pathlines in 3d
% Save image of result
[~,~,~,xyzlims] = QS.getXYZLims() ;
nTimePoints = 45 ;
first = true ;
slab = 20 ;
black_figs = true ;

if ~computed_v3d 
    load(sprintf(QS.fileName.pathlines.v3d, t0), 'v3dPathlines')
    vX3rs = v3dPathlines.vXrs ;
    vY3rs = v3dPathlines.vYrs ;
    vZ3rs = v3dPathlines.vZrs ;
    t0 = v3dPathlines.t0 ;
    tIdx0 = v3dPathlines.tIdx0 ;
    computed_v3d = false ;
end

% Get subsampling
qsubX = round(size(vX3rs, 2) / nY2plot) - 1 ;
qsub = 1 ;

% Obtain indices for subsampling
idmat = false([size(vX3rs, 2), size(vX3rs, 3)]) ;
% This would make gridlike sampling
% ap2do = 1:qsubX:(size(vX3rs, 2)-1)

% Saggital sampling
halfV = floor(QS.nV*0.5) ;
ap2do = [1, halfV] ;  % 2, QS.nV-1,halfV-1, halfV, halfV + 1] ;
for qq = ap2do 
    idmat(1:qsubX:end, qq) = true; 
end
inds = find(idmat) ;
clearvars idmat
xx = vX3rs(tIdx0, inds) ;
yy = vY3rs(tIdx0, inds) ;
zz = vZ3rs(tIdx0, inds) ;

% % Define the lagrangian points that lie in saggital plane at t=0
% yslab = find(abs(yy)<slab) ;
% plot3(xx(yslab), yy(yslab), zz(yslab), '.')

for tidx = 1:length(QS.xp.fileMeta.timePoints)
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    
    plineFig = fullfile(sprintf(pathlineDir.v3d, t0), ...
        [sprintf(QS.fileBase.name, tp) '.png']) ;
    
    if ~exist(plineFig, 'file') || overwrite
        % Plot the pathlines -- NOTE that we flip YY coordinates (MATLAB)
        close all
        set(gcf, 'visible', 'off')
        tp2plot = max(1, tidx-nTimePoints):tidx ;
        
        if black_figs
            colors = bone(length(tp2plot)) ;
            colormap(bone)
        else
            colors = flipud(bone(length(tp2plot))) ;
            colormap(flipud(bone))
        end
        if size(colors, 1) == 1
            colors = bone(2) ;
            colors = colors(2, :) ;
        end
        
        alphas = linspace(0, 1, length(tp2plot)) ;
        for qind = 1:length(tp2plot)
            qq = tp2plot(qind) ;
            xx = vX3rs(qq, inds) ;
            yy = vY3rs(qq, inds) ;
            zz = vZ3rs(qq, inds) ;
            % xx = xx(yslab) ;
            % yy = yy(yslab) ;
            % zz = zz(yslab) ;
            scatter3(xx(:), yy(:), zz(:), ...
                scatterSz * 2, ...
                'markerfacecolor', colors(qind, :), ...
                'markeredgecolor', 'none', 'markerfacealpha', alphas(qind)) ;
            % s.MarkerFaceAlpha = (max(1.0 - 0.1 * yy(:), 0)) ;
            hold on;
        end
        view(0,0)
        axis equal
        scatter([], [])
        xlim(xyzlims(1, :))
        ylim(xyzlims(2, :))
        zlim(xyzlims(3, :))
        daspect([1 QS.a_fixed 1])
        % title('pathlines: mesh vertices', 'Interpreter', 'Latex')
        xlabel(['ap position [' QS.spaceUnits ']'], 'Interpreter', 'Latex')
        ylabel(['lateral position [' QS.spaceUnits ']'], 'Interpreter', 'Latex')
        zlabel(['dv position [' QS.spaceUnits ']'], 'Interpreter', 'Latex')
        
        axis off
        
        if black_figs
            set(gcf, 'color', 'k')
            set(gca, 'color', 'k', 'xcol', 'k', 'ycol', 'k', 'zcol', 'k')
        else
            set(gca, 'color', 'w')
        end

        if black_figs
            cb = colorbar('eastOutside', 'color', 'w') ;
        else
            cb = colorbar('eastOutside') ;
        end
        
        % Adjust size of colorbar
        cbheight = cb.Position(4) ;
        cb.Position(1) = cb.Position(1) + 0.02 ;
        cb.Position(4) = cb.Position(4) * 0.6 ;
        cb.Position(2) = cb.Position(2) + 0.15 * cbheight ;
        ax = gca ;
        ax.Position(1) = ax.Position(1) - 0.05 ; 
        
        caxis([(QS.xp.fileMeta.timePoints(max(tp2plot))-t0-nTimePoints) * QS.timeInterval, ...
           (QS.xp.fileMeta.timePoints(max(tp2plot))-t0) * QS.timeInterval]) ;
        ylabel(cb, ['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')

        % saveas(gcf, plineFig)
        disp(['Saving ' plineFig])
        export_fig(plineFig, '-nocrop', '-r150')
        close all
    else
        disp(['Mesh advected vertex XYZ pathline image already on disk: ' plineFig])
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Streamline version has issues but included here for reference:
% % Build pathlines starting at t=t0 
% % NOTE: If timePoints does not have a uniform dt, there is no problem:
% % we define the brick of worldlines with uniform spacing and stepsize=1
% % so that paths follow optical flow between frames and vz is 1: ie, the
% % velocities connect adjacent frames, whatever the dt between them.
% xG = repmat(x0, [1, 1, ntps]) ;
% yG = repmat(y0, [1, 1, ntps]) ;
% % make z grid increasing along third (z) dimension
% zG = ones(size(xG)) ;
% for qq = 1:ntps
%     zG(:, :, qq) = qq ;
% end
% % check each xyzgrid as image
% if preview
%     for qq = 1:ntps
%         imagesc(squeeze(xG(:, :, qq)))
%         title(['t=' num2str(qq)])
%         pause(0.1)
%     end
%     for qq = 1:ntps
%         imagesc(squeeze(yG(:, :, qq)))
%         title(['t=' num2str(qq)])
%         pause(0.1)
%     end
%     for qq = 1:ntps
%         imagesc(squeeze(zG(:, :, qq)))
%         title(['t=' num2str(qq)])
%         caxis([0, max(max(zG(:)), max(xG(:)))])
%         pause(0.1)
%     end
% end
%
% vzdat = ones(size(xG)) ;
% sLineOptions = [1] ;  % stepsize is dt 
% % z0 is the starting time array in the worldline loaf
% z0 = QS.xp.tIdx(t0) * ones(size(x0)) ; 
% 
% strLineM = streamline(xG, yG, zG, ...
%     squeeze(vPIV(:, :, :, 1)), ...
%     squeeze(vPIV(:, :, :, 2)), vzdat, ...
%     x0(:), y0(:), z0(:), sLineOptions) ;
% % reshape pathlines into grid format
% strLineM = reshape(strLineM, [ntps, size(x0, 1), size(x0, 2)], 2) ;

