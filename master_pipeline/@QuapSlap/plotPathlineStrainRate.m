function plotPathlineStrainRate(QS, options)
% plotPathlineStrainRate(QS, options)
%   Plot the strain rate along pathlines as kymographs and 
%   correlation plots.
% 
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   plot_kymographs         : bool
%   plot_kymographs_cumsum  : bool
%   plot_correlations       : bool 
%   plot_fold_kinematics    : bool
%   plot_lobe_kinematics    : bool
% 
% Returns
% -------
% <none>
%
% NPMitchell 2020

%% Default options 
overwrite = false ;
plot_kymographs = true ;
plot_kymographs_strain = true ;
plot_fold_kinematics = true ;
plot_lobe_kinematics = true ;
t0 = QS.t0set() ;

%% Parameter options
lambda = 0.02 ;
lambda_mesh = 0.002 ;
climit = 0.2 ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 

%% Unpack options & assign defaults
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
%% parameter options
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_err')
    lambda_err = options.lambda_err ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
end
if isfield(options, 't0')
    t0 = options.t0 ;
end
if isfield(options, 'climit')
    climit = options.climit ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end

%% Operational options
if isfield(options, 'plot_kymographs')
    plot_kymographs = options.plot_kymographs ;
end
if isfield(options, 'plot_kymographs_cumprod')
    plot_kymographs_strain = options.plot_kymographs_cumprod ;
end
if isfield(options, 'plot_fold_kinematics')
    plot_fold_kinematics = options.plot_fold_kinematics ;
end
if isfield(options, 'plot_lobe_kinematics')
    plot_lobe_kinematics = options.plot_lobe_kinematics ;
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
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;
sKDir = fullfile(QS.dir.strainRate.root, ...
    strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f'], ...
    lambda, lambda_mesh), '.', 'p'));
folds = load(QS.fileName.fold) ;
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;

%% Colormap
close all
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;
% clf
% set(gcf, 'visible', 'off')
% imagesc([-1, 0, 1; -1, 0, 1])
% caxis([0, 1])
% pos256 = bluewhitered(256) ;
close all
pm256 = phasemap(256) ;


%% Load time offset for first fold, t0
QS.t0set() ;
tfold = QS.t0 ;

%% load from QS
if doubleResolution
    nU = QS.nU * 2 - 1 ;
    nV = QS.nV * 2 - 1 ;
else
    nU = QS.nU ;
    nV = QS.nV ;    
end

%% Test incompressibility of the flow on the evolving surface
% We relate the normal velocities to the divergence / 2 * H.
tps = QS.xp.fileMeta.timePoints(1:end-1) - tfold;

% preallocate for cumulative error
ntps = length(QS.xp.fileMeta.timePoints(1:end-1)) ;
dv_apM   = zeros(ntps, nU) ;   % dv averaged
tr_apM   = zeros(ntps, nU) ;   
th_apM   = zeros(ntps, nU) ;   
dv_lM   = zeros(ntps, nU) ;    % left averaged
tr_lM = zeros(ntps, nU) ;
th_lM = zeros(ntps, nU) ;
dv_rM   = zeros(ntps, nU) ;    % right averaged
tr_rM = zeros(ntps, nU) ;
th_rM = zeros(ntps, nU) ;
dv_dM   = zeros(ntps, nU) ;    % dorsal averaged
tr_dM = zeros(ntps, nU) ;
th_dM = zeros(ntps, nU) ;
dv_vM   = zeros(ntps, nU) ;    % ventral averaged
tr_vM = zeros(ntps, nU) ;
th_vM = zeros(ntps, nU) ;

% Output directory is inside metricKinematics dir
mKPDir = fullfile(sKDir, sprintf('pathline_%04dt0', t0)) ;
datdir = fullfile(mKPDir, 'measurements') ;
% Data for kinematics on meshes (defined on vertices) [not needed here]
% mdatdir = fullfile(mKDir, 'measurements') ;

% Unit definitions for axis labels
unitstr = [ '[1/' QS.timeUnits ']' ];

%% Compute or load all timepoints
apKymoFn = fullfile(datdir, 'apKymographPathlineStrainRate.mat') ;
lKymoFn = fullfile(datdir, 'leftKymographPathlineStrainRate.mat') ;
rKymoFn = fullfile(datdir, 'rightKymographPathlineStrainRate.mat') ;
dKymoFn = fullfile(datdir, 'dorsalKymographPathlineStrainRate.mat') ;
vKymoFn = fullfile(datdir, 'ventralKymographPathlineStrainRate.mat') ;
files_exist = exist(apKymoFn, 'file') && ...
    exist(lKymoFn, 'file') && exist(rKymoFn, 'file') && ...
    exist(dKymoFn, 'file') && exist(vKymoFn, 'file') ;
if files_exist
    load(apKymoFn, 'tr_apM', 'dv_apM', 'th_apM', ...
        'str_apM', 'sdv_apM', 'sth_apM')
    load(lKymoFn, 'tr_lM', 'dv_lM', 'th_lM', ...
        'str_lM', 'sdv_lM', 'sth_lM')
    load(rKymoFn, 'tr_rM', 'dv_rM', 'th_rM', ...
        'str_rM', 'sdv_rM', 'sth_rM')
    load(dKymoFn, 'tr_dM', 'dv_dM', 'th_dM', ...
        'str_dM', 'sdv_dM', 'sth_dM')
    load(vKymoFn, 'tr_vM', 'dv_vM', 'th_vM', ...
        'str_vM', 'sdv_vM', 'sth_vM')
else
    for tp = QS.xp.fileMeta.timePoints(1:end-1)
        close all
        disp(['t = ' num2str(tp)])
        tidx = QS.xp.tIdx(tp) ;

        % Check for timepoint measurement on disk
        srfn = fullfile(datdir, sprintf('strainRate_series_%06d.mat', tp))   ;

        % Load timeseries measurements
        load(srfn, 'treps_ap', 'treps_l', 'treps_r', 'treps_d', 'treps_v', ...
            'dvtre_ap', 'dvtre_l', 'dvtre_r', 'dvtre_d', 'dvtre_v', ...
            'theta_ap', 'theta_l', 'theta_r', 'theta_d', 'theta_v', ...
            'strain_tr_ap', 'strain_tr_l', 'strain_tr_r', ...
            'strain_tr_d', 'strain_tr_v', ...
            'strain_theta_ap', 'strain_theta_l', 'strain_theta_r', ....
            'strain_theta_d', 'strain_theta_v', ...
            'strain_dv_ap', 'strain_dv_l', 'strain_dv_r', ...
            'strain_dv_d', 'strain_dv_v') ;

        %% Store in matrices
        % dv averaged
        tr_apM(tidx, :) = treps_ap ;
        dv_apM(tidx, :) = dvtre_ap ;
        th_apM(tidx, :) = theta_ap ;

        % left quarter
        tr_lM(tidx, :) = treps_l ;
        dv_lM(tidx, :) = dvtre_l ;
        th_lM(tidx, :) = theta_l ;

        % right quarter
        tr_rM(tidx, :) = treps_r ;
        dv_rM(tidx, :) = dvtre_r ;
        th_rM(tidx, :) = theta_r ;

        % dorsal quarter
        tr_dM(tidx, :) = treps_d ;
        dv_dM(tidx, :) = dvtre_d ;
        th_dM(tidx, :) = theta_d ;

        % ventral quarter
        tr_vM(tidx, :) = treps_v ;
        dv_vM(tidx, :) = dvtre_v ;
        th_vM(tidx, :) = theta_v ;

        %% Store accumulated strain in matrices
        % dv averaged
        str_apM(tidx, :) = strain_tr_ap ;
        sdv_apM(tidx, :) = strain_dv_ap ;
        sth_apM(tidx, :) = strain_theta_ap ;

        % left quarter
        str_lM(tidx, :) = strain_tr_l ;
        sdv_lM(tidx, :) = strain_dv_l ;
        sth_lM(tidx, :) = strain_theta_l ;

        % right quarter
        str_rM(tidx, :) = strain_tr_r ;
        sdv_rM(tidx, :) = strain_dv_r ;
        sth_rM(tidx, :) = strain_theta_r ;

        % dorsal quarter
        str_dM(tidx, :) = strain_tr_d ;
        sdv_dM(tidx, :) = strain_dv_d ;
        sth_dM(tidx, :) = strain_theta_d ;

        % ventral quarter
        str_vM(tidx, :) = strain_tr_v ;
        sdv_vM(tidx, :) = strain_dv_v ;
        sth_vM(tidx, :) = strain_theta_v ;
    end
    
    %% Save kymographs
    save(apKymoFn, 'tr_apM', 'dv_apM', 'th_apM', ...
        'str_apM', 'sdv_apM', 'sth_apM')
    save(lKymoFn, 'tr_lM', 'dv_lM', 'th_lM', ...
        'str_lM', 'sdv_lM', 'sth_lM')
    save(rKymoFn, 'tr_rM', 'dv_rM', 'th_rM', ...
        'str_rM', 'sdv_rM', 'sth_rM')
    save(dKymoFn, 'tr_dM', 'dv_dM', 'th_dM', ...
        'str_dM', 'sdv_dM', 'sth_dM')
    save(vKymoFn, 'tr_vM', 'dv_vM', 'th_vM', ...
        'str_vM', 'sdv_vM', 'sth_vM')
end

%% Store kymograph data in cell arrays
trsK = {tr_apM, tr_lM, tr_rM, tr_dM, tr_vM} ;
dvsK = {dv_apM, dv_lM, dv_rM, dv_dM, dv_vM} ;
thsK = {th_apM, th_lM, th_rM, th_dM, th_vM} ;

%% Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
dvDir = fullfile(mKPDir, 'avgDV') ;
lDir = fullfile(mKPDir, 'avgLeft') ;
rDir = fullfile(mKPDir, 'avgRight') ;
dDir = fullfile(mKPDir, 'avgDorsal') ;
vDir = fullfile(mKPDir, 'avgVentral') ;
outdirs = {dvDir, lDir, rDir, dDir, vDir} ;

%% To grab fold locationsrobustly, find minima of div(v) from ap average
% and grab folds and lobes indexed in Lagrangian coords
try
    metricKDir = QS.dir.metricKinematics.root ;
    dirs = dir(metricKDir) ; 
    for qq = 1:length(dirs)
        if contains(dirs(qq).name, strrep(sprintf('lambda%0.3f_', lambda), '.', 'p'))
            if contains(dirs(qq).name, strrep(sprintf('lmesh%0.3f', lambda_mesh), '.', 'p'))
                metricKDir = fullfile(metricKDir, dirs(qq).name) ;
            end
        end
    end
    loadDir = fullfile(metricKDir, sprintf('pathline_%04dt0', t0), ...
        'measurements') ;
    apKymoFn = fullfile(loadDir, 'apKymographsMetricKinematics.mat') ;
    load(apKymoFn, 'divv_apM')
catch
    error('Run QS.plotPathlineMetricKinematics() before QS.plotPathlineStrainRate()')
end
div1d = mean(divv_apM(tps > 20 & tps < 50, :), 1) ;
div1dsm = savgol(div1d, 2, 11) ;
[~, valleys] = maxk(-islocalmin(div1dsm) .* div1dsm, 3) ;
valleys = sort(valleys) ;

%% Now plot different measured quantities as kymographs
clim_zoom = climit * 0.3 ;
if plot_kymographs
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    for qq = 1:length(outdirs)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end

        % Unpack what to plot (averaged kymographs, vary averaging region)
        trK = trsK{qq} ;
        dvK = dvsK{qq} ;
        thK = thsK{qq} ;
        
        titles = {'dilation rate, $\textrm{Tr}[g^{-1}\dot{\varepsilon}]$',...
            'shear rate, $||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\dot{\varepsilon}\right]\bf{g}||$'} ;
        labels = {['$\textrm{Tr}[g^{-1}\dot{\varepsilon}]$ ' unitstr], ...
            ['$||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\dot{\varepsilon}\right]\bf{g}||$' unitstr]} ;
        names = {'dilation', 'deviator'} ;
        
        %% Plot trace/deviator DV-averaged kymograph
        for pp = 1:2
            
            % Check if images already exist on disk
            fn = fullfile(odir, [ names{pp} '.png']) ;
            fn_zoom = fullfile(odir, [names{pp} '_zoom_early.png']) ;
            
            if ~exist(fn, 'file') || ~exist(fn_zoom, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                if pp == 1
                    imagesc((1:nU)/nU, tps, trK)
                    caxis([-climit, climit])
                else
                    % Map intensity from dvtre and color from the theta
                    indx = max(1, round(mod(2*thK(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
                    colors = pm256(indx, :) ;
                    dvtreKclipped = min(dvK / climit, 1) ;
                    colorsM = dvtreKclipped(:) .* colors ;
                    colorsM = reshape(colorsM, [size(dvK, 1), size(dvK, 2), 3]) ;
                    imagesc((1:nU)/nU, tps, colorsM)
                    caxis([0, climit])
                end
                colormap(bwr256)

                % Plot fold identifications
                hold on;
                fons1 = max(1, fons(1)) ;
                fons2 = max(1, fons(2)) ;
                fons3 = max(1, fons(3)) ;
                t1ones = ones(size(tps(fons1:end))) ;
                t2ones = ones(size(tps(fons2:end))) ;
                t3ones = ones(size(tps(fons3:end))) ;
                tidx0 = QS.xp.tIdx(t0) ;

                % OPTION 1: use QS.folds
                % plot(folds.folds(tidx0, 1) * t1ones / nU, tps(fons1:end))
                % plot(folds.folds(tidx0, 2) * t2ones / nU, tps(fons2:end))
                % plot(folds.folds(tidx0, 3) * t3ones / nU, tps(fons3:end))

                % OPTION 1: use identified div(v) < 0
                plot(valleys(1) * t1ones / nU, tps(fons1:end))
                plot(valleys(2) * t2ones / nU, tps(fons2:end))
                plot(valleys(3) * t3ones / nU, tps(fons3:end))
                
                % % Add folds to plot
                % hold on;
                % tidx0 = QS.xp.tIdx(t0) ;
                % % Which is the first fold (at t0)?
                % [~, minID] = min(fons) ;
                % thisFoldTimes = tps(fons(minID):end) ;
                % t_ones = ones(size(thisFoldTimes)) ;
                % plot(folds.folds(tidx0, minID) * t_ones / nU, thisFoldTimes)
                % % Now plot the later folds at their Lagrangian locations
                % % determined by div(v) being very negative.
                % laterID = setdiff(1:length(fons), minID) ;
                % for ii = laterID
                %     thisFoldTimes = tps(fons(ii):end) ;
                %     t_ones = ones(size(thisFoldTimes)) ;
                %     plot(valleys(ii) * t_ones / nU, thisFoldTimes)
                % end
                
                % Titles 
                title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')

                if pp == 1
                    cb = colorbar() ;
                elseif pp == 2
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
                    caxis([0, climit])
                    colormap(gca, gray)
                end
                
                % title and save
                ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
                fn = fullfile(odir, [ names{pp} '.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   

                if pp == 1
                    % Zoom in on small values
                    caxis([-clim_zoom, clim_zoom])
                    colormap(bwr256)
                    fn = fullfile(odir, [names{pp} '_zoom.png']) ;
                    disp(['saving ', fn])
                    export_fig(fn, '-png', '-nocrop', '-r200')   
                    % Zoom in on early times
                    ylim([min(tps), max(fons) + 10])
                    caxis([-clim_zoom, clim_zoom])
                    colormap(bwr256)
                    fn = fullfile(odir, [names{pp} '_zoom_early.png']) ;
                    disp(['saving ', fn])
                    export_fig(fn, '-png', '-nocrop', '-r200')   
                elseif pp == 2
                    % Zoom in on color limits by changing intensity clip
                    close all
                    set(gcf, 'visible', 'off')
                    
                    % Map intensity from dvtre and color from the theta
                    indx = max(1, round(mod(2*thK(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
                    colors = pm256(indx, :) ;
                    dvtreKclipped = min(dvK / clim_zoom, 1) ;
                    colorsM = dvtreKclipped(:) .* colors ;
                    colorsM = reshape(colorsM, [size(dvK, 1), size(dvK, 2), 3]) ;
                    imagesc((1:nU)/nU, tps, colorsM)
                    caxis([0, clim_zoom])
                    
                    % Add folds to plot
                    hold on;
                    
                    % Plot fold identifications
                    hold on;
                    fons1 = max(1, fons(1)) ;
                    fons2 = max(1, fons(2)) ;
                    fons3 = max(1, fons(3)) ;
                    t1ones = ones(size(tps(fons1:end))) ;
                    t2ones = ones(size(tps(fons2:end))) ;
                    t3ones = ones(size(tps(fons3:end))) ;
                    tidx0 = QS.xp.tIdx(t0) ;
                    
                    % OPTION 1: use QS.folds
                    % plot(folds.folds(tidx0, 1) * t1ones / nU, tps(fons1:end))
                    % plot(folds.folds(tidx0, 2) * t2ones / nU, tps(fons2:end))
                    % plot(folds.folds(tidx0, 3) * t3ones / nU, tps(fons3:end))
                    
                    % OPTION 1: use identified div(v) < 0
                    plot(valleys(1) * t1ones / nU, tps(fons1:end))
                    plot(valleys(2) * t2ones / nU, tps(fons2:end))
                    plot(valleys(3) * t3ones / nU, tps(fons3:end))
                
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
                    caxis([0, clim_zoom])
                    colormap(gca, gray)

                    % title and save
                    title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
                    ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                    xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                    ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
                    
                    % Zoom in on small values
                    fn = fullfile(odir, [names{pp} '_zoom.png']) ;
                    disp(['saving ', fn])
                    export_fig(fn, '-png', '-nocrop', '-r200')   
                    % Zoom in on early times
                    ylim([min(tps), max(fons) + 10])
                    fn = fullfile(odir, [names{pp} '_zoom_early.png']) ;
                    disp(['saving ', fn])
                    export_fig(fn, '-png', '-nocrop', '-r200')   
                end
            end
        end
    end
end

%% Kymographs of cumulative products along pathlines 
strK = {str_apM, str_lM, str_rM, str_dM, str_vM} ;
sdvK = {sdv_apM, sdv_lM, sdv_rM, sdv_dM, sdv_vM} ;
sthK = {sth_apM, sth_lM, sth_rM, sth_dM, sth_vM} ;

if plot_kymographs_strain
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;
    
    
    for qq = 1:length(outdirs)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end

        % Unpack what to plot (averaged kymographs, vary averaging region)
        trK_pos = cumprod(1 + strK{qq}(tps>eps, :), 1) ;
        trK_neg = flipud(cumprod(flipud(1 ./ (1 + strK{qq}(tps<eps, :))), 1)) ;
        trK = cat(1, trK_neg, trK_pos) ;
        
        labels = {['$\Pi$d$t \, \left[1+\textrm{Tr}[g^{-1}\dot{\varepsilon}]\right)$ ' unitstr], ...
            ['$\Pi \,$d$t \, \left( 1 + ||\varepsilon-\frac{1}{2}\mathrm{Tr}\left[\mathbf{g}^{-1}\dot{\varepsilon}\right)\bf{g}||\right)$' unitstr]} ;
        
        titles = {'dilation, $\Pi$d$t \, \left[1+\textrm{Tr}[g^{-1}\dot{\varepsilon}]\right)$ ',...
            'shear, $\Pi \,$d$t \, \left( 1 + ||\varepsilon-\frac{1}{2}\mathrm{Tr}\left[\mathbf{g}^{-1}\dot{\varepsilon}\right)\bf{g}||\right)$'} ;
        names = {'Idilation_t0', 'Ideviator_t0'} ;
        climitWide = climit * 3; 
        
        %% Plot dilation DV-averaged Lagrangian pathline kymograph 
        % Check if images already exist on disk
        fn = fullfile(odir, [ names{1} '.png']) ;
        fn_zoom = fullfile(odir, [names{1} '_zoom_early.png']) ;
        if ~exist(fn, 'file') || ~exist(fn_zoom, 'file') || overwrite
            close all
            set(gcf, 'visible', 'off')
            imagesc((1:nU)/nU, tps, trK)
            caxis([1-climitWide, 1+climitWide])
            colormap(bwr256)

            % Plot fold identifications
            hold on;
            fons1 = max(1, fons(1)) ;
            fons2 = max(1, fons(2)) ;
            fons3 = max(1, fons(3)) ;
            t1ones = ones(size(tps(fons1:end))) ;
            t2ones = ones(size(tps(fons2:end))) ;
            t3ones = ones(size(tps(fons3:end))) ;
            tidx0 = QS.xp.tIdx(t0) ;

            % OPTION 1: use QS.folds
            % plot(folds.folds(tidx0, 1) * t1ones / nU, tps(fons1:end))
            % plot(folds.folds(tidx0, 2) * t2ones / nU, tps(fons2:end))
            % plot(folds.folds(tidx0, 3) * t3ones / nU, tps(fons3:end))

            % OPTION 1: use identified div(v) < 0
            plot(valleys(1) * t1ones / nU, tps(fons1:end))
            plot(valleys(2) * t2ones / nU, tps(fons2:end))
            plot(valleys(3) * t3ones / nU, tps(fons3:end))

            % title and save
            title([titles{1}, titleadd{1}], 'Interpreter', 'Latex')
            ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
            xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
            cb = colorbar() ;
            ylabel(cb, labels{1}, 'Interpreter', 'Latex')  
            disp(['saving ', fn])
            export_fig(fn, '-png', '-nocrop', '-r200')   

            % Zoom in on small values
            caxis([1-climitWide, 1+climitWide])
            colormap(bwr256)
            fn = fullfile(odir, [names{1} '_zoom.png']) ;
            tmp = strsplit(fn, filesep) ;
            disp(['saving ', tmp{end}])
            export_fig(fn, '-png', '-nocrop', '-r200')   
            % Zoom in on early times
            ylim([min(tps), max(fons) + 10])
            caxis([1-climitWide, 1+climitWide])
            colormap(bwr256)
            tmp = strsplit(fn_zoom, filesep) ;
            disp(['saving ', tmp{end}])
            export_fig(fn_zoom, '-png', '-nocrop', '-r200')   
        end
        
        %% Plot shear DV-averaged Lagrangian pathline kymograph 
        % Check if images already exist on disk
        fn = fullfile(odir, [ names{2} '.png']) ;
        fn_zoom = fullfile(odir, [names{2} '_zoom_early.png']) ;
        if ~exist(fn, 'file') || ~exist(fn_zoom, 'file') || overwrite
            close all
            set(gcf, 'visible', 'off')
            % Map intensity from dvtre and color from the theta
            indx = max(1, round(mod(2*sthK{qq}(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
            colors = pm256(indx, :) ;
            dvtreKclipped = min(sdvK{qq} / climitWide, 1) ;
            colorsM = dvtreKclipped(:) .* colors ;
            colorsM = reshape(colorsM, [size(sdvK{qq}, 1), size(sdvK{qq}, 2), 3]) ;
            imagesc((1:nU)/nU, tps, colorsM)
            caxis([1-climitWide, 1+climitWide])

            % Plot fold identifications
            hold on;
            fons1 = max(1, fons(1)) ;
            fons2 = max(1, fons(2)) ;
            fons3 = max(1, fons(3)) ;
            t1ones = ones(size(tps(fons1:end))) ;
            t2ones = ones(size(tps(fons2:end))) ;
            t3ones = ones(size(tps(fons3:end))) ;
            tidx0 = QS.xp.tIdx(t0) ;

            % OPTION 1: use QS.folds
            % plot(folds.folds(tidx0, 1) * t1ones / nU, tps(fons1:end))
            % plot(folds.folds(tidx0, 2) * t2ones / nU, tps(fons2:end))
            % plot(folds.folds(tidx0, 3) * t3ones / nU, tps(fons3:end))

            % OPTION 1: use identified div(v) < 0
            plot(valleys(1) * t1ones / nU, tps(fons1:end))
            plot(valleys(2) * t2ones / nU, tps(fons2:end))
            plot(valleys(3) * t3ones / nU, tps(fons3:end))

            % Titles 
            title([titles{2}, titleadd{2}], 'Interpreter', 'Latex')
            ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
            xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')

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
            caxis([0, climitWide])
            colormap(gca, gray)

            % title and save
            ylabel(cb, labels{2}, 'Interpreter', 'Latex')  
            fn = fullfile(odir, [ names{2} '.png']) ;
            disp(['saving ', fn])
            export_fig(fn, '-png', '-nocrop', '-r200')   
            
            % Zoom in on small values
            fn = fullfile(odir, [names{2} '_zoom.png']) ;
            disp(['saving ', fn])
            export_fig(fn, '-png', '-nocrop', '-r200')   
            % Zoom in on early times
            ylim([min(tps), max(fons) + 10])
            fn = fullfile(odir, [names{2} '_early.png']) ;
            disp(['saving ', fn])
            export_fig(fn, '-png', '-nocrop', '-r200')   
        end
        
    end
end

%% Plot 1D curves for region around each fold
% Sample divv near each fold
foldw = 0.05 ;
endw = 0.10 ;
cut = [round(endw*nU), valleys(1)-round(foldw*nU), ...
        valleys(1)+round(foldw*nU), valleys(2)-round(foldw*nU), ...
        valleys(2)+round(foldw*nU), valleys(3)-round(foldw*nU), ... 
        valleys(3)+round(foldw*nU), round((1-endw) * nU)] ;
lobes = { cut(1):cut(2), cut(3):cut(4), cut(5):cut(6), cut(7):cut(8) } ;
avgStrings = {'dv-averaged', 'left side', 'right side', ...
    'dorsal side', 'ventral side'} ;
avgLabel = {'dv', 'left', 'right', 'dorsal', 'ventral'} ;
titleFoldBase = 'Lagrangian metric kinematics along folds, ' ;
titleLobeBase = 'Lagrangian metric kinematics along lobes, ' ;

foldYlabels = {'anterior fold', 'middle fold', 'posterior fold'} ;
lobeYlabels = {'lobe 1', 'lobe 2', 'lobe 3', 'lobe 4'} ;
for qq = 1:5 
    divv = divvsK{qq} ;
    % Explore a range of widths for plotting
    for width = 1:round(0.05 * nU)
        fn = fullfile(outdirs{qq}, ...
            [sprintf('fold_kinematics_w%03d_', 2*width+1), ...
            avgLabel{qq}, '.png']) ;
        if ( ~exist(fn, 'file') || overwrite  ) && plot_fold_kinematics
            close all
            ymin = 0 ;
            ymax = 0 ;
            % Each fold is valley+/- width
            for jj = 1:length(valleys)
                axisColl{jj} = subplot(length(valleys), 1, jj) ;
                valley = (valleys(jj)-width):(valleys(jj)+width) ;
                dvj = mean(divv(:, valley), 2) ;
                Hvj = mean(H2vn(:, valley), 2) ;
                plot(tps, dvj, '.-', 'Color', QS.plotting.colors(1, :))
                hold on;
                plot(tps, Hvj, '.-', 'Color', QS.plotting.colors(2, :))
                if jj == 1
                    % Title and labels
                    sgtitle([titleFoldBase, avgStrings{qq}, ', ', ...
                        '$w_{\textrm{fold}}=', ...
                        num2str(100*(2*width + 1)/ nU), '$\%$\, L_\zeta$'], ...
                        'Interpreter', 'Latex')
                    legend({'$\nabla\cdot\mathbf{v}_\parallel$', ...
                        '$v_n 2H$'}, 'Interpreter', 'Latex', ...
                        'location', 'eastOutside')  
                    drawnow
                    pos = get(gca, 'position') ;
                elseif jj == length(valleys)
                    xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                end
                ylabel(foldYlabels{jj}, 'Interpreter', 'Latex')
                ylims = ylim() ;
                ymin = min(ymin, ylims(1)) ;
                ymax = max(ymax, ylims(2)) ;
            end

            for jj = 1:length(valleys)
                axes(axisColl{jj})
                pos2 = get(gca, 'position') ;
                set(gca, 'position', [pos2(1) pos2(2) pos(3) pos2(4)])
                ylim([ymin, ymax])

                % Mark wherever the divv<0 and also vn2H<0
                valley = (valleys(jj)-width):(valleys(jj)+width) ;
                dvj = mean(divv(:, valley), 2) ;
                Hvj = mean(H2vn(:, valley), 2) ;
                dvpos = dvj > 0 ;
                Hvpos = Hvj > 0 ;
                scatter(tps(dvpos), ...
                    (ymin-0.05*(ymax-ymin)) * ones(size(tps(dvpos))), 5, ...
                    'markeredgecolor', 'none', 'markerFaceAlpha', 0.6, ...
                    'markerFaceColor', QS.plotting.colors(1, :), ...
                    'HandleVisibility', 'off')
                hold on;
                scatter(tps(Hvpos), ...
                    ymin * ones(size(tps(Hvpos))), 5,  ...
                    'markeredgecolor', 'none', 'markerFaceAlpha', 0.6, ...
                    'markerFaceColor', QS.plotting.colors(2, :), ...
                    'HandleVisibility', 'off')
                ylim([ymin-0.1*(ymax-ymin), ymax])
            end

            % Save figure
            disp(['Saving figure: ', fn])
            saveas(gcf, fn)
        end
    end
    
    %% Fold kinematics -- cumprod gdot
    for width = 1:round(0.05 * nU)
        fn = fullfile(outdirs{qq}, ...
            [sprintf('fold_kinematics_gdot_cumprod_w%03d_', 2*width+1), ...
            avgLabel{qq}, '.png']) ;
        if ( ~exist(fn, 'file') || overwrite ) && plot_fold_kinematics
            close all
            ymin = 0 ;
            ymax = 0 ;
            % Each fold is valley+/- width
            for jj = 1:length(valleys)
                axisColl{jj} = subplot(length(valleys), 1, jj) ;
                valley = (valleys(jj)-width):(valleys(jj)+width) ;
                
                % div(v), H*2*vn, gdot
                ddj = mean(divv(:, valley), 2) ;
                Hdj = mean(H2vn(:, valley), 2) ;
                gdj = mean(divv(:, valley) - H2vn(:, valley), 2) ;
                
                % Take cumulative product marching forward from t0
                gpj_pos = cumprod(1 + gdj(tps > eps)) ;
                gpj_neg = flipud(cumprod(flipud(1 ./ (1 + gdj(tps < eps))))) ;
                gpj = cat(1, gpj_neg, gpj_pos) ;               
                % Take cumulative product marching forward from t0
                dpj_pos = cumprod(1 + ddj(tps > eps)) ;
                dpj_neg = flipud(cumprod(flipud(1 ./ (1 + ddj(tps < eps))))) ;
                dpj = cat(1, dpj_neg, dpj_pos) ;
                % Take cumulative product marching forward from t0
                Hpj_pos = cumprod(1 + Hdj(tps > eps)) ;
                Hpj_neg = flipud(cumprod(flipud(1 ./ (1 + Hdj(tps < eps))))) ;
                Hpj = cat(1, Hpj_neg, Hpj_pos) ;
                
                % Plot all three
                plot(tps, dpj, '.-', 'Color', QS.plotting.colors(1, :))
                hold on;
                plot(tps, Hpj, '.-', 'Color', QS.plotting.colors(2, :))
                plot(tps, gpj, '.-', 'Color', QS.plotting.colors(3, :))
                
                if jj == 1
                    % Title and labels
                    sgtitle([titleFoldBase, avgStrings{qq}, ', ', ...
                        '$w_{\textrm{fold}}=', ...
                        num2str(100*(2*width + 1)/ nU), '$\%$\, L_\zeta$'], ...
                        'Interpreter', 'Latex')
                    legend({'$\Pi(1+\nabla\cdot\mathbf{v}_\parallel)$', ...
                        '$\Pi(1+v_n 2H)$', ...
                        '$\Pi(1+\frac{1}{2}\mathrm{Tr}\left[g^{-1} \dot{g} \right])$'}, ...
                        'Interpreter', 'Latex', 'location', 'eastOutside')  
                    drawnow
                    pos = get(gca, 'position') ;
                elseif jj == length(valleys)
                    xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                end
                ylabel(foldYlabels{jj}, 'Interpreter', 'Latex')
                ylims = ylim() ;
                ymin = min(ymin, ylims(1)) ;
                ymax = max(ymax, ylims(2)) ;
            end

            for jj = 1:length(valleys)
                axes(axisColl{jj})
                pos2 = get(gca, 'position') ;
                set(gca, 'position', [pos2(1) pos2(2) pos(3) pos2(4)])
                ylim([ymin, ymax])

                % Mark wherever the divv<0 and also vn2H<0
                valley = (valleys(jj)-width):(valleys(jj)+width) ;
                dvj = mean(divv(:, valley), 2) ;
                Hvj = mean(H2vn(:, valley), 2) ;
                dvpos = dvj > 0 ;
                Hvpos = Hvj > 0 ;
                scatter(tps(dvpos), ...
                    (ymin-0.05*(ymax-ymin)) * ones(size(tps(dvpos))), 5, ...
                    'markeredgecolor', 'none', 'markerFaceAlpha', 0.6, ...
                    'markerFaceColor', QS.plotting.colors(1, :), ...
                    'HandleVisibility', 'off')
                hold on;
                scatter(tps(Hvpos), ...
                    ymin * ones(size(tps(Hvpos))), 5,  ...
                    'markeredgecolor', 'none', 'markerFaceAlpha', 0.6, ...
                    'markerFaceColor', QS.plotting.colors(2, :), ...
                    'HandleVisibility', 'off')
                ylim([ymin-0.1*(ymax-ymin), ymax])
            end

            % Save figure
            disp(['Saving figure: ', fn])
            saveas(gcf, fn)
        end
        
        % Plot all gdots on one axis
        fn = fullfile(outdirs{qq}, ...
            [sprintf('fold_kinematics_gdot_cumprod_compare_w%03d_', 2*width+1), ...
            avgLabel{qq}, '.png']) ;
        if ~exist(fn, 'file') || overwrite 
            close all
            % Each fold is valley+/- width
            for jj = 1:length(valleys)
                valley = (valleys(jj)-width):(valleys(jj)+width) ;
                
                % div(v), H*2*vn, gdot
                gdj = mean(divv(:, valley) - H2vn(:, valley), 2) ;
                
                % Take cumulative product marching forward from t0
                gpj_pos = cumprod(1 + gdj(tps > eps)) ;
                gpj_neg = flipud(cumprod(flipud(1 ./ (1 + gdj(tps < eps))))) ;
                gpj = cat(1, gpj_neg, gpj_pos) ;               

                % Plot this fold
                plot(tps, gpj, '.-', 'Color', QS.plotting.colors(jj+3, :))
                hold on;
            end    
            
            % Title and labels
            sgtitle(['Tissue dilation in folds, ', avgStrings{qq}, ', ', ...
                '$w_{\textrm{fold}}=', ...
                num2str(100*(2*width + 1)/ nU), '$\%$\, L_\zeta$'], ...
                'Interpreter', 'Latex')
            legend(foldYlabels, 'Interpreter', 'Latex', 'location', 'eastOutside')  
            drawnow
            xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
            ylabel('$\frac{1}{2}\mathrm{Tr} \left[g^{-1} \dot{g}\right]$', ...
                'Interpreter', 'Latex')
        end
        
        % Save figure
        disp(['Saving figure: ', fn])
        saveas(gcf, fn)
    end
   
    %% Plot lobe Kinematics
    fn = fullfile(outdirs{qq}, ...
        ['lobe_kinematics_' avgLabel{qq} '.png']) ;
    if ( ~exist(fn, 'file') || overwrite ) && plot_lobe_kinematics
        close all
        ymin = 0 ;
        ymax = 0 ;
        % Each fold is valley+/- width
        for jj = 1:length(lobes)
            axisColl{jj} = subplot(length(lobes), 1, jj) ;
            dvj = mean(divv(:, lobes{jj}), 2) ;
            Hvj = mean(H2vn(:, lobes{jj}), 2) ;
            plot(tps, dvj, '.-', 'Color', QS.plotting.colors(1, :))
            hold on;
            plot(tps, Hvj, '.-', 'Color', QS.plotting.colors(2, :))
            if jj == 1
                % Title and labels
                sgtitle([titleLobeBase, avgStrings{qq}], ...
                    'Interpreter', 'Latex')
                legend({'$\nabla\cdot\mathbf{v}_\parallel$', ...
                    '$v_n 2H$'}, 'Interpreter', 'Latex', ...
                    'location', 'eastOutside')  
                drawnow
                pos = get(gca, 'position') ;
            elseif jj == length(valleys)
                xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
            end
            ylabel(lobeYlabels{jj}, 'Interpreter', 'Latex')
            ylims = ylim() ;
            ymin = min(ymin, ylims(1)) ;
            ymax = max(ymax, ylims(2)) ;
        end

        for jj = 1:length(lobes)
            axes(axisColl{jj})
            pos2 = get(gca, 'position') ;
            set(gca, 'position', [pos2(1) pos2(2) pos(3) pos2(4)])
            ylim([ymin, ymax])

            % Mark wherever the divv<0 and also vn2H<0
            dvj = mean(divv(:, lobes{jj}), 2) ;
            Hvj = mean(H2vn(:, lobes{jj}), 2) ;
            dvpos = dvj > 0 ;
            Hvpos = Hvj > 0 ;
            scatter(tps(dvpos), ...
                (ymin-0.05*(ymax-ymin)) * ones(size(tps(dvpos))), 5, ...
                'markeredgecolor', 'none', 'markerFaceAlpha', 0.6, ...
                'markerFaceColor', QS.plotting.colors(1, :), ...
                'HandleVisibility', 'off')
            hold on;
            scatter(tps(Hvpos), ...
                ymin * ones(size(tps(Hvpos))), 5,  ...
                'markeredgecolor', 'none', 'markerFaceAlpha', 0.6, ...
                'markerFaceColor', QS.plotting.colors(2, :), ...
                'HandleVisibility', 'off')
            ylim([ymin-0.1*(ymax-ymin), ymax])
        end

        % Save figure
        disp(['Saving figure: ', fn])
        saveas(gcf, fn)
    end
    
    %% Lobe kinematics -- cumprod gdot
    fn = fullfile(outdirs{qq}, ...
        sprintf('lobe_kinematics_gdot_cumprod.png')) ;
    if ( ~exist(fn, 'file') || overwrite ) && plot_lobe_kinematics
        close all
        ymin = 0 ;
        ymax = 0 ;
        for jj = 1:length(lobes)
            axisColl{jj} = subplot(length(lobes), 1, jj) ;
            
            % div(v), H*2*vn, gdot
            ddj = mean(divv(:, lobes{jj}), 2) ;
            Hdj = mean(H2vn(:, lobes{jj}), 2) ;
            gdj = mean(divv(:, lobes{jj}) - H2vn(:, lobes{jj}), 2) ;

            % Take cumulative product marching forward from t0
            gpj_pos = cumprod(1 + gdj(tps > eps)) ;
            gpj_neg = flipud(cumprod(flipud(1 ./ (1 + gdj(tps < eps))))) ;
            gpj = cat(1, gpj_neg, gpj_pos) ;               
            % Take cumulative product marching forward from t0
            dpj_pos = cumprod(1 + ddj(tps > eps)) ;
            dpj_neg = flipud(cumprod(flipud(1 ./ (1 + ddj(tps < eps))))) ;
            dpj = cat(1, dpj_neg, dpj_pos) ;
            % Take cumulative product marching forward from t0
            Hpj_pos = cumprod(1 + Hdj(tps > eps)) ;
            Hpj_neg = flipud(cumprod(flipud(1 ./ (1 + Hdj(tps < eps))))) ;
            Hpj = cat(1, Hpj_neg, Hpj_pos) ;

            % Plot all three
            plot(tps, dpj, '.-', 'Color', QS.plotting.colors(1, :))
            hold on;
            plot(tps, Hpj, '.-', 'Color', QS.plotting.colors(2, :))
            plot(tps, gpj, '.-', 'Color', QS.plotting.colors(3, :))

            if jj == 1
                % Title and labels
                sgtitle([titleLobeBase, avgStrings{qq}], ...
                    'Interpreter', 'Latex')
                legend({'$\Pi(1+\nabla\cdot\mathbf{v}_\parallel)$', ...
                    '$\Pi(1+v_n 2H)$', ...
                    '$\Pi(1+\frac{1}{2}\mathrm{Tr}\left[g^{-1} \dot{g} \right])$'}, ...
                    'Interpreter', 'Latex', 'location', 'eastOutside')  
                drawnow
                pos = get(gca, 'position') ;
            elseif jj == length(valleys)
                xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
            end
            ylabel(lobeYlabels{jj}, 'Interpreter', 'Latex')
            ylims = ylim() ;
            ymin = min(ymin, ylims(1)) ;
            ymax = max(ymax, ylims(2)) ;
        end

        for jj = 1:length(lobes)
            axes(axisColl{jj})
            pos2 = get(gca, 'position') ;
            set(gca, 'position', [pos2(1) pos2(2) pos(3) pos2(4)])
            ylim([ymin, ymax])

            % Mark wherever the divv<0 and also vn2H<0
            dvj = mean(divv(:, lobes{jj}), 2) ;
            Hvj = mean(H2vn(:, lobes{jj}), 2) ;
            dvpos = dvj > 0 ;
            Hvpos = Hvj > 0 ;
            scatter(tps(dvpos), ...
                (ymin-0.05*(ymax-ymin)) * ones(size(tps(dvpos))), 5, ...
                'markeredgecolor', 'none', 'markerFaceAlpha', 0.6, ...
                'markerFaceColor', QS.plotting.colors(1, :), ...
                'HandleVisibility', 'off')
            hold on;
            scatter(tps(Hvpos), ...
                ymin * ones(size(tps(Hvpos))), 5,  ...
                'markeredgecolor', 'none', 'markerFaceAlpha', 0.6, ...
                'markerFaceColor', QS.plotting.colors(2, :), ...
                'HandleVisibility', 'off')
            ylim([ymin-0.1*(ymax-ymin), ymax])
        end

        % Save figure
        disp(['Saving figure: ', fn])
        saveas(gcf, fn)
    end

    % Plot all gdots on one axis
    fn = fullfile(outdirs{qq}, ...
        sprintf('lobe_kinematics_gdot_cumprod_compare.png')) ;
    if ( ~exist(fn, 'file') || overwrite ) && plot_lobe_kinematics
        close all
        % Each fold is valley+/- width
        for jj = 1:length(lobes)
            % div(v), H*2*vn, gdot
            gdj = mean(divv(:, lobes{jj}) - H2vn(:, lobes{jj}), 2) ;

            % Take cumulative product marching forward from t0
            gpj_pos = cumprod(1 + gdj(tps > eps)) ;
            gpj_neg = flipud(cumprod(flipud(1 ./ (1 + gdj(tps < eps))))) ;
            gpj = cat(1, gpj_neg, gpj_pos) ;               

            % Plot this fold
            plot(tps, gpj, '.-', 'Color', QS.plotting.colors(jj+3, :))
            hold on;
        end    

        % Title and labels
        sgtitle(['Tissue dilation in lobes, ', avgStrings{qq}, ...
            ', $\frac{1}{2}\mathrm{Tr} \left[g^{-1} \dot{g}\right]$'], ...
            'Interpreter', 'Latex')
        legend(lobeYlabels, 'Interpreter', 'Latex', 'location', 'eastOutside')  
        drawnow
        xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
        ylabel('$\frac{1}{2}\mathrm{Tr} \left[g^{-1} \dot{g}\right]$', ...
            'Interpreter', 'Latex')
        
        % Save figure
        disp(['Saving figure: ', fn])
        saveas(gcf, fn)
    end

end

disp('done')



