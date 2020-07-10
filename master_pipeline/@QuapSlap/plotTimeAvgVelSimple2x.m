function plotTimeAvgVelSimple2x(QS, samplingResolution, options) 
%plotTimeAvgVelSimple(QS, options) 
%   Save images of the velocity field over time that has been "simply"
%   averaged over time in-place in the surface-Lagrangian pullback. That
%   is, the velocity field at location (u,v) has been averaged in time with
%   previous and later timepoints at the same pullback location (u,v).
%   
%
% Parameters
% ----------
% QS : QuapSlap class instance
% samplingResolution : str specifier ('1x' or '2x', 'single' or 'double')
%   the resolution in which to sample the Pullback for velocities
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%   timePoints : numeric 1D array
%       the timepoints to consider for the measurement. For ex, could
%       choose subset of the QS experiment timePoints
%   alphaVal : float
%       the opacity of the heatmap to overlay
%   invertImage : bool
%       invert the data pullback under the velocity map
%
%
%
% NPMitchell 2020

%% Default options
% Declare plotting options for limits
vtscale = 5 ;      % um / min
vnscale = 2 ;      % um / min
vscale = 2 ;       % um / min
alphaVal = 0.7 ;   % alpha for normal velocity heatmap
washout2d = 0.5 ;  % lightening factor for data
qsubsample = 5 ;   % quiver subsampling in pullback space 
plot_vxyz = false ;
overwrite = false ;
preview = false ;
pivimCoords = 'sp_sme' ;
doubleCovered = true;

%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
else 
    error("Could not parse samplingResolution: set to '1x' or '2x'")
end

%% Unpack options
if isfield(options, 'plot_vxyz')
    plot_vxyz = options.plot_vxyz ;
end
if isfield(options, 'pivimCoords')
    pivimCoords = options.pivimCoords ;
    if strcmp(pivimCoords(-1), 'e')
        doubleCovered = true ;
    else
        doubleCovered = false ;
    end
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
    
%% Unpack QS
if doubleResolution
    pivSADir = QS.dir.pivSimAvg2x ;
else
    pivSADir = QS.dir.pivSimAvg ;
end
timePoints = QS.xp.fileMeta.timePoints ;
QS.getFeatures('ssfold')
ssfold_frac = QS.features.ssfold / QS.nU ;
t0 = QS.t0set() ;
tunit = [' ' QS.timeunits] ;

% Auxiliary function for plotting smoothed meshes

%% Define directories
pivSAImXDir = fullfile(pivSADir, 'vx') ;
pivSAImYDir = fullfile(pivSADir, 'vy') ;
pivSAImZDir = fullfile(pivSADir, 'vz') ;
pivSAImTDir = fullfile(pivSADir, 'vtH') ;  % Heatmap
pivSAImGDir = fullfile(pivSADir, 'vtG ') ;  % Gaussian smoothed in space
pivSAImSDir = fullfile(pivSADir, 'vmag ') ;  % speed |v_3D|
pivSAImNDir = fullfile(pivSADir, 'vn') ;
% pivSimAvgImQDir = fullfile(pivSimAvgDir, 'vtQ') ;  % Quiverplot
% pivSimAvgImDvgDir = fullfile(pivSimAvgDir, 'dvg') ;
% pivSimAvgImCurlDir = fullfile(pivSimAvgDir, 'curl') ;
% pivSimAvgImShearDir = fullfile(pivSimAvgDir, 'shear_dvphi_ds') ;

pivDir = QS.dir.piv ;
dilDir = fullfile(pivDir, 'dilation') ;
vxyorigDir = fullfile(pivDir, 'vxyorig') ;
if plot_vxyz
    ensureDir(pivSAImXDir)
    ensureDir(pivSAImYDir)
    ensureDir(pivSAImZDir)
end
dirs2make = {pivSAImTDir, pivSAImGDir,...
    pivSAImNDir, dilDir, vxyorigDir, pivSAImSDir} ;
for pp = 1:length(dirs2make)
    ensureDir(dirs2make{pp}) ;
end

%% Load piv results
disp('Loading raw piv to get field size')
piv = load(QS.fileName.pivRaw) ;
% get size of images to make
gridsz = size(piv.x{1}) ;

%% Load simple average piv results
if doubleResolution
    QS.getVelocitySimpleAverage2x()
    velstruct = QS.velocitySimpleAverage2x ;
else
    QS.getVelocitySimpleAverage()
    velstruct = QS.velocitySimpleAverage ;
end
vsmM = velstruct.v3d ;
v2dsmMum = velstruct.v2dum ;
vnsmM = velstruct.vn ;

%% Make plots
% Display the velocities
close all
fig = figure('visible', 'off') ;
for i = 1:size(vsmM, 1)
    tp = timePoints(i) ;
    
    % Check if normal velocity plot exists
    vnfn = fullfile(pivSAImNDir, [sprintf(QS.fileBase.name, tp) '.png']) ;
    vthfn = fullfile(pivSAImTDir, [sprintf(QS.fileBase.name, tp) '.png']) ;
    vtgfn = fullfile(pivSAImGDir, [sprintf(QS.fileBase.name, tp) '.png']) ;
    vxyorigfn = fullfile(vxyorigDir, [sprintf(QS.fileBase.name, tp) '.png']) ;
    speedfn = fullfile(pivSAImSDir, [sprintf(QS.fileBase.name, tp) '.png']) ;
    
    % grab the tangential velocity for this timestep
    vsm_ii = squeeze(vsmM(i, :, :)) ;
    v2dsmum_ii = squeeze(v2dsmMum(i, :, :)) ;
    vnsm_ii = squeeze(vnsmM(i, :, :)) ;

    % Load the image to put flow on top
    if strcmp(pivimCoords, 'sp_sme')
        im = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
        ylims = [0.25 * size(im, 1), 0.75 * size(im, 1)] ;
    else
        error(['Have not coded for this pivimCoords option. Do so here: ' pivimCoords])
    end
    im = cat(3, im, im, im) ;  % convert to rgb for no cmap change

    % Define Nx1 and Mx1 float arrays for xspace and yspace
    xx = piv.x{i}(1, :) ;
    yy = piv.y{i}(:, 1) ;

    % Define the proper coordinates (approx)
    % if ~exist(dvgfn, 'file') || ~exist(curlfn, 'file') || overwrite || true
    %     dilation_allfaces = zeros(length(jac), 1) ;
    %     for f = 1:length(jac)
    %         qg = jac{f} * jac{f}' ;
    %         dilation_allfaces(f) = sqrt(det(qg)) ;
    %     end
    % end
    %     % load the spcutMesh
    %     load(sprintf(spcutMeshBase, tp), 'spcutMesh')
    %     ss = spcutMesh.sphi(:, 1) ;
    %     pp = spcutMesh.sphi(:, 1) ;
    %     [ss, pp] = gridDistancesInterpMesh3D() ;
    % end

    if plot_vxyz
        aux_plot_vxyz_simpleavg(im, vsm_ii, xx, yy, tp, vscale, ...
            pivSAImXDir, pivSAImYDir, pivSAImZDir) ;
    end

    % Plot the magnitude of the velocity 
    if ~exist(speedfn, 'file') || overwrite
        disp(['Saving ' speedfn])
        close all
        fig = figure('units', 'normalized', ...
                'outerposition', [0 0 1 1], 'visible', 'off') ;
        colormap parula ;
        scalarFieldOnImage(im, xx, yy, reshape(vecnorm(vsm_ii, 2, 2), gridsz),...
            alphaVal, vtscale, '$|v|$ [$\mu$m/min]', 'Style', 'Positive') ;
        ylim([0.25 * size(im, 1), 0.75 * size(im, 1)])
        title(['speed, $|v|$: $t=$' num2str(tp - t0) tunit], ...
            'Interpreter', 'Latex')
        saveas(fig, speedfn) ;
        close all
    end

    % % Check normals
    % quiver3(piv3d{i}.pt0(:, 1), piv3d{i}.pt0(:, 2), piv3d{i}.pt0(:, 3), ...
    %     piv3d{i}.normals(:, 1), piv3d{i}.normals(:, 2), piv3d{i}.normals(:, 3)) ;
    % 
    % % Check normals rotated and scaled
    % pt0_rs = ((rot * piv3d{i}.pt0')' + trans) * resolution  ;
    % quiver3(pt0_rs(:, 1), pt0_rs(:, 2), pt0_rs(:, 3), ...
    %     piv3d{i}.normals_rs(:, 1), piv3d{i}.normals_rs(:, 2), piv3d{i}.normals_rs(:, 3)) ;
    %

    % Look at smoothed 2d velocity fields
    vn = reshape(vnsm_ii, gridsz) ;
    vx = reshape(v2dsmum_ii(:, 1), gridsz) ;
    vy = reshape(v2dsmum_ii(:, 2), gridsz) ;

    % % Check normal velocity and rotations
    % quiver3(piv.x{i}, piv.y{i}, 0*piv.x{i}, vx, vy, vn, 0)

    % Get lobes for this timepoint
    foldx = ssfold_frac(i, :) * size(im, 2) ;

    % Plot the normal velocity on top
    if ~exist(vnfn, 'file') || overwrite
        disp(['Saving ' vnfn])
        close all
        fig = figure('units', 'normalized', ...
                'outerposition', [0 0 1 1], 'visible', 'off') ;
        scalarFieldOnImage(im, xx, yy, vn, alphaVal, vnscale, ...
            '$v_n$ [$\mu$m/min]') ;
        ylim(ylims)
        title(['normal velocity, $v_n$: $t=$' num2str(tp - t0) tunit], ...
            'Interpreter', 'Latex')
        saveas(fig, vnfn) ;
        close all
    end

    % Plot the tangential velocity as heatmap on top of the image
    if ~exist(vthfn, 'file') || overwrite
        disp(['Saving ' vthfn])
        imw = im * washout2d + max(im(:)) * (1-washout2d) ;
        qopts.overlay_quiver = false ;
        qopts.qsubsample = qsubsample ;
        qopts.overlay_quiver = true ;
        qopts.qscale = 10 ;
        qopts.label = '$v_t$ [$\mu$m/min]' ;
        qopts.title = ['tangential velocity, $v_t$: $t=$' num2str(tp - t0) tunit] ;
        qopts.outfn = vthfn ;
        qopts.ylim = ylims ;
        vectorFieldHeatPhaseOnImage(imw, xx, yy, vx, vy, vtscale, qopts) ;
        clear qopts 
    end

    % Gaussian smooth the velocities
    if ~exist(vtgfn, 'file') || overwrite
        disp(['Saving ' vtgfn])
        vxb = imgaussfilt(vx, 4) ;
        vyb = imgaussfilt(vy, 4) ;
        imw = im * washout2d + max(im(:)) * (1-washout2d) ;
        qopts.qsubsample = qsubsample ;
        qopts.overlay_quiver = true ;
        qopts.qscale = 10 ;
        qopts.label = '$v_t$ [$\mu$m/min]' ;
        qopts.title = ['Smoothed velocity, $v_t$: $t=$', ...
            num2str(tp - t0), tunit] ;
        qopts.outfn = vtgfn ;
        qopts.ylim = ylims ;
        % Plot the coarse-grained tang velocity as heatmap on top of the image
        vectorFieldHeatPhaseOnImage(imw, xx, yy, vxb, vyb, vtscale, qopts) ;    
        clearvars qopts
    end

    % Find hyperbolic fixed points
    %

    % Plot original velocity -- note no smoothing done here ;)
    if ~exist(vxyorigfn, 'file') || overwrite
        imw = im * washout2d + max(im(:)) * (1-washout2d) ;
        opts.label = '$\tilde{v}$ [pix/min]' ;
        opts.outfn = vxyorigfn ;
        opts.qscale = 15 ;
        qopts.ylim = ylims ;
        qopts.title = ['piv results, $t=$', num2str(tp - t0), tunit] ;
        vectorFieldHeatPhaseOnImage(imw, xx, yy, ...
            piv.u_filtered{i}, piv.v_filtered{i}, 15, opts) ;
        clearvars opts
    end
        
end




% RETIRED CODE: DIVERGENCE, CURL, and SHEAR (these are crude, replaced
% by DEC) 
%
%
% shearfn = fullfile(pivSimAvgImShearDir, [sprintf('%04d', tp) '.png']) ;
%
% % Plot divergence
% if ~exist(dvgfn, 'file') || overwrite 
% 
%     vxb = imgaussfilt(vx, 10) ;
%     vyb = imgaussfilt(vy, 10) ;
% 
%     % Interpolate dilation onto locations where curl is defined
%     % Di = scatteredInterpolant(tm0X, tm0Y, dilation_allfaces) ;
%     % tr0 = triangulation(tm0f, [tm0X, tm0Y]) ;
%     % [subfieldfaces, ~] = pointLocation(tr0, [xx(:), yy(:)]) ;
%     % dilv = dilation_allfaces(subfieldfaces) ;
% 
%     dilum = reshape(piv3d{i}.dilation / resolution, size(vxb)) ;
%     dvg = divergence(piv.x{i}, piv.y{i}, vxb, vyb) .* dilum;
%     opts.label = '$\nabla \cdot v_t$ [min$^{-1}$]' ;
%     opts.title = [ '$t=$' num2str(tp) ' min' ] ;
%     opts.outfn = dvgfn ;
%     opts.qscale = 10 ;
%     opts.sscale = 1.0 ;
%     opts.alpha = 0.8 ;
%     opts.ylim = [size(im, 2) * 0.25, size(im, 2) * 0.75] ;
%     scalarVectorFieldsOnImage(im, xx, yy, ...
%         dvg, xx, yy, vxb, vyb, opts) ;
%     clearvars opts
% end
% 
% % Plot curl
% if ~exist(curlfn, 'file') || overwrite || true
%     % xd = xx(1:10:end) ;
%     % yd = yy(1:10:end) ;
%     % xdgrid = xd .* ones(length(xd), length(yd)) ; 
%     % ydgrid = (xd .* ones(length(yd), length(xd)))' ; 
%     % check it
%     % scatter(xdgrid(:), ydgrid(:), 10, xdgrid(:))
% 
%     % Interpolate dilation onto locations where curl is defined
%     % Di = scatteredInterpolant(tm0X, tm0Y, dilation_allfaces) ;
%     % tr0 = triangulation(tm0f, [tm0X, tm0Y]) ;
%     % [subfieldfaces, ~] = pointLocation(tr0, [xx(:), yy(:)]) ;
%     % dilv = dilation_allfaces(subfieldfaces) ;
% 
%     dilum = reshape(piv3d{i}.dilation / resolution, size(vxb)) ;
%     curlv = curl(piv.x{i}, piv.y{i}, vxb, vyb) ;
%     curlv = curlv .* dilum ;
%     opts.title = [ '$t=$' num2str(tp) ' min' ] ;
%     opts.label = '$\nabla \times v_t$ [min$^{-1}$]' ;
%     opts.outfn = curlfn ;
%     opts.qscale = 10 ;
%     opts.sscale = 0.5 ;
%     opts.alpha = 0.8 ;
%     opts.ylim = [size(im, 2) * 0.25, size(im, 2) * 0.75] ;
%     scalarVectorFieldsOnImage(im, xx, yy, ...
%         curlv, xx, yy, vxb, vyb, opts) ;
%     clearvars opts
% end
% 
% % Plot strain dv_phi/ds on image
% if ~exist(shearfn, 'file') || overwrite || true
%     % xd = xx(1:10:end) ;
%     % yd = yy(1:10:end) ;
%     % xdgrid = xd .* ones(length(xd), length(yd)) ; 
%     % ydgrid = (xd .* ones(length(yd), length(xd)))' ; 
%     % check it
%     % scatter(xdgrid(:), ydgrid(:), 10, xdgrid(:))
% 
%     % Interpolate dilation onto locations where curl is defined
%     % Di = scatteredInterpolant(tm0X, tm0Y, dilation_allfaces) ;
%     % tr0 = triangulation(tm0f, [tm0X, tm0Y]) ;
%     % [subfieldfaces, ~] = pointLocation(tr0, [xx(:), yy(:)]) ;
%     % dilv = dilation_allfaces(subfieldfaces) ;
% 
%     dilum = reshape(piv3d{i}.dilation / resolution, size(vxb)) ;
%     dvphidX = gradient(vyb, xx(2) - xx(1), yy(2) - yy(1)) ;
%     dvphids = dvphidX .* dilum ;
%     opts.title = [ '$t=$' num2str(tp) ' min' ] ;
%     opts.label = '$\nabla_s v_{\phi}$ [min$^{-1}$]' ;
%     opts.outfn = shearfn ;
%     opts.qscale = 10 ;
%     opts.sscale = 0.5 ;
%     opts.alpha = 0.8 ;
%     opts.ylim = [size(im, 2) * 0.25, size(im, 2) * 0.75] ;
%     scalarVectorFieldsOnImage(im, xx, yy, ...
%         dvphids, xx, yy, vxb, vyb, opts) ;
%     clearvars opts
% end