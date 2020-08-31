function plotStrainRate(QS, options)
%plotStrainRateTimePoint(QS, tp, options)
%   Plot the traceful and traceless components of the strain rate tensor
%   defined on each face, both for individual timepoints, and 
%   also over time as kymographs
%
% Parameters
% ----------
% QS : QuapSlap class instance
% tp : int 
%   timepoint in units of (1/QS.timeInterval) * QS.timeUnits
% options: struct with fields
%   
% Returns
% -------
% 
% 
% NPMitchell 2020

%% Unpack required params
lambda = options.lambda ;
lambda_mesh = options.lambda_mesh ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
debug = false ;

%% Parameters
overwrite = false ;
clim_trace = 0.05 ;
clim_deviatoric = 0.05 ;
averagingStyle = 'Lagrangian' ;
skipTimePoint = false ;
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'mesh')
    mesh = options.mesh ;
end
if isfield(options, 'cutMesh')
    cutMesh = options.cutMesh ;
end
if isfield(options, 'clim_trace')
    clim_trace = options.clim_trace ;
end
if isfield(options, 'clim_deviatoric')
    clim_deviatoric = options.clim_deviatoric ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'averagingStyle')
    averagingStyle = options.averagingStyle ;
end
if isfield(options, 'debug')
    debug = options.debug ;
end
if isfield(options, 'skipTimePoint')
    skipTimePoint = options.skipTimePoint ;
end

%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
    sresStr = '' ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
    sresStr = 'doubleRes_' ;
else 
    error("Could not parse samplingResolution: set to '1x' or '2x'")
end

%% Unpack QS
t0 = QS.t0set() ;
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
srlambdaDir = fullfile(QS.dir.strainRate.root, ...
    strrep(sprintf('lambda%0.3f_lmesh%0.3f', ...
    lambda, lambda_mesh), '.', 'p')) ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;
nU = QS.nU ;
nV = QS.nV ;
folds = load(QS.fileName.fold) ;
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;

%% Prepare directories for images
dirs2make = { srlambdaDir, ...
    fullfile(srlambdaDir, 'strainRate3d'), ...
    fullfile(srlambdaDir, 'strainRate2d') } ;
for ii = 1:length(dirs2make)
    dir2make = dirs2make{ii} ;
    if ~exist(dir2make, 'dir')
        mkdir(dir2make)
    end
end

%% Colormap
close all
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;
clf
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([0, 1])
pos256 = bluewhitered(256) ;
close all
pm256 = phasemap(256) ;

%% Collate data from each quarter of the gut
if strcmp(averagingStyle, 'simple')
    sKDir = fullfile(QS.dir.strainRateSimple, ...
        strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f'], ...
        lambda, lambda_mesh), '.', 'p'));
else
    sKDir = fullfile(QS.dir.strainRate.root, ...
        strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f'], ...
        lambda, lambda_mesh), '.', 'p'));
end

for tp = QS.xp.fileMeta.timePoints(1:end-1)
    tidx = QS.xp.tIdx(tp) ;
    
    %% Define metric strain filename        
    estrainFn = fullfile(srlambdaDir, 'measurements', ...
        sprintf(QS.fileBase.strainRate, tp)) ;
    disp(['t=' num2str(tp) ': Loading strainrate results from disk: ' estrainFn])
    load(estrainFn, 'strainrate', 'treps', 'dvtre', 'theta', ...
        'dvtre_ap', 'dvtre_l', 'dvtre_r', 'dvtre_d', 'dvtre_v', ...
        'theta_ap', 'theta_l', 'theta_r', 'theta_d', 'theta_v', ...
        'lambda', 'lambda_mesh')
    
    %% Plot strain rate for this timepoint
    if ~skipTimePoint
        tpOpts.lambda = lambda ;
        tpOpts.lambda_mesh = lambda_mesh ;
        tpOpts.overwrite = overwrite ;
        QS.plotStrainRateTimePoint(tp, tpOpts)
    end
    
    %% Store in matrices    
    % dv averaged
    dvtre_apM(tidx, :) = dvtre_ap ;
    theta_apM(tidx, :) = theta_ap ;
    
    % left quarter
    dvtre_lM(tidx, :) = dvtre_l ;
    theta_lM(tidx, :) = theta_l ;
    
    % right quarter
    dvtre_rM(tidx, :) = dvtre_r ;
    theta_rM(tidx, :) = theta_r ;

    % dorsal quarter
    dvtre_dM(tidx, :) = dvtre_d ;
    theta_dM(tidx, :) = theta_d ;

    % ventral quarter
    dvtre_vM(tidx, :) = dvtre_v ;
    theta_vM(tidx, :) = theta_v ;
    
end

%% Store kymograph data in cell arrays
dvtresK = {dvtre_apM, dvtre_lM, dvtre_rM, dvtre_dM, dvtre_vM} ;
thetasK = {theta_apM, theta_lM, theta_rM, theta_dM, theta_vM} ;

%% Now plot different measured quantities as kymographs
% Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
dvDir = fullfile(sKDir, 'avgDV') ;
lDir = fullfile(sKDir, 'avgLeft') ;
rDir = fullfile(sKDir, 'avgRight') ;
dDir = fullfile(sKDir, 'avgDorsal') ;
vDir = fullfile(sKDir, 'avgVentral') ;
outdirs = {dvDir, lDir, rDir, dDir, vDir} ;
titleadd = {': circumferentially averaged', ...
    ': left side', ': right side', ': dorsal side', ': ventral side'} ;
tps = QS.xp.fileMeta.timePoints(1:end-1) ;
titlestr = '$||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\varepsilon\right]\bf{g}||$' ;
for qq = 1:length(outdirs)
    % Prep the output directory for this averaging
    odir = outdirs{qq} ;
    if ~exist(odir, 'dir')
        mkdir(odir)
    end

    % denom = sqrt(tg(:, 1, 1) .* tg(:, 2, 2)) ;
    % NOTE: \varepsilon --> ${\boldmath${\varepsilon}$}$
    label = '$||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\varepsilon\right]\bf{g}||$' ;
    name = 'dvtre' ;

    %% Plot DV-averaged/quarter-avgeraged kymograph            
    % Check if images already exist on disk
    fn = fullfile(odir, [ name '.png']) ;
    fn_zoom = fullfile(odir, [name '_zoom_early.png']) ;
    zoomstrs = {'', '_zoom'} ;
    
    if ~exist(fn, 'file') || ~exist(fn_zoom, 'file') || overwrite
        for pp = 1:2
            zoomstr = zoomstrs{pp} ;        
            if pp == 1
                clim = clim_deviatoric ;
            else
                clim = clim_deviatoric * 0.3 ;
            end
            
            % Unpack what to plot (averaged kymographs, vary averaging region)
            dvtreK = dvtresK{qq} ;
            thetaK = thetasK{qq} ;

            % Map intensity from dvtre and color from the theta
            indx = max(1, round(mod(2*thetaK(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
            colors = pm256(indx, :) ;
            dvtreKclipped = min(dvtreK / clim, 1) ;
            colorsM = dvtreKclipped(:) .* colors ;
            colorsM = reshape(colorsM, [size(dvtreK, 1), size(dvtreK, 2), 3]) ;

            % Plot the kymograph
            close all
            set(gcf, 'visible', 'off')
            imagesc((1:nU)/nU, tps, colorsM)
            caxis([-clim, clim])
            colormap(bwr256)
            % Add folds to plot
            hold on;
            fons1 = max(1, fons(1)) ;
            fons2 = max(1, fons(2)) ;
            fons3 = max(1, fons(3)) ;
            plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
            plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
            plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

            % Colorbar and phasewheel
            colormap(gca, phasemap)
            phasebar('colormap', phasemap, ...
                'location', [0.82, 0.7, 0.1, 0.135], 'style', 'nematic')
            ax = gca ;
            get(gca, 'position')
            cb = colorbar('location', 'eastOutside') ;
            drawnow
            axpos = get(ax, 'position') ;
            cbpos = get(cb, 'position') ;
            set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3), cbpos(4)*0.6])
            set(ax, 'position', axpos) 
            hold on;
            caxis([0, clim])
            colormap(gca, gray)

            % title and save
            title([titlestr, titleadd{qq}], 'Interpreter', 'Latex')
            ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
            xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
            ylabel(cb, label, 'Interpreter', 'Latex')  
            fn = fullfile(odir, [ name zoomstr '.png']) ;
            disp(['saving ', fn])
            export_fig(fn, '-png', '-nocrop', '-r200')   

            if pp == 2
                % Zoom in on early times
                ylim([min(tps), max(fons) + 10])
                fn = fullfile(odir, [name zoomstr '_early.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
            end
        end
    end
end
disp('done with plotting strain rate')


