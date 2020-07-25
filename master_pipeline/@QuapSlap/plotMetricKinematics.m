function plotMetricKinematics(QS, options)
% plotMetricKinematics(QS, options)
%   Plot the metric Kinematics as kymographs and correlation plots
%   Out-of-plane motion is v_n * 2H, where v_n is normal velocity and H is
%   mean curvature.
%   In-plane motion considered here is div(v_t) where v_t is tangential
%   velocity on the curved surface.
%   The difference div(v_t) - vn*2H = Tr[g^{-1} dot{g}], which is a measure
%   of isotropic metric change over time (dot is the partial derivative wrt
%   time). 
% 
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   plot_kymographs : bool
%   plot_kymographs_cumsum : bool
%   plot_correlations : bool
%   plot_gdot_correlations : bool
%   plot_gdot_decomp : bool
% 
% NPMitchell 2020

%% Default options 
overwrite = false ;
plot_kymographs = true ;
plot_kymographs_cumsum = true ;
plot_correlations = true ;
plot_gdot_correlations = false ;
plot_gdot_decomp = true ;

%% Parameter options
lambda = 0.01 ; 
% by default, lambda_mesh = lambda, whether defined here or in options
lambda_err = 0.01 ;
climit = 0.2 ;
climit_err = 0.2 ;
climit_veln = climit * 10 ;
climit_H = climit * 2 ;
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
else
    % default lambda_mesh is equal to lambda 
    lambda_mesh = lambda ;
end
if isfield(options, 'climit')
    climit = options.climit ;
end
if isfield(options, 'climit_err')
    climit_err = options.climit_err ;
end
if isfield(options, 'climit_veln')
    climit_veln = options.climit_veln ;
end
if isfield(options, 'climit_H')
    climit_H = options.climit_H ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end

%% Operational options
if isfield(options, 'plot_kymographs')
    plot_kymographs = options.plot_kymographs ;
end
if isfield(options, 'plot_kymographs_cumsum')
    plot_kymographs_cumsum = options.plot_kymographs_cumsum ;
end
if isfield(options, 'plot_correlations')
    plot_correlations = options.plot_correlations ;
end
if isfield(options, 'plot_gdot_correlations')
    plot_gdot_correlations = options.plot_gdot_correlations ;
end
if isfield(options, 'plot_gdot_decomp')
    plot_gdot_decomp = options.plot_gdot_decomp ;
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
mKDir = fullfile(QS.dir.metricKinematics, ...
    strrep(sprintf([sresStr 'lambda%0.3f_lerr%0.3f_lmesh%0.3f'], ...
    lambda, lambda_err, lambda_mesh), '.', 'p'));
folds = load(QS.fileName.fold) ;
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;

%% Colormap
bwr256 = bluewhitered(256) ;

%% Load vertex-based velocity measurements
if doubleResolution
    vvsmMfn = fullfile(QS.dir.pivSimAvg2x, 'vvM_simpletimeavg2x.mat')  ;
    tmp = load(vvsmMfn) ;
    vertex_vels = tmp.vvsmM ;
else
    vvsmMfn = fullfile(QS.dir.pivSimAvg, 'vvM_simpletimeavg.mat')  ; 
    tmp = load(vvsmMfn) ;
    vertex_vels = tmp.vvsmM ;   
end
% vfsmMfn = fullfile(QS.dir.pivSimAvg, 'vfM_simpletimeavg.mat') ;
% tmp = load(vfsmMfn) ;
% face_vels = tmp.vfsmM ;

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
HH_apM   = zeros(ntps, nU) ;   % dv averaged
divv_apM = zeros(ntps, nU) ;
veln_apM = zeros(ntps, nU) ;
gdot_apM = zeros(ntps, nU) ;
HH_lM   = zeros(ntps, nU) ;    % left averaged
divv_lM = zeros(ntps, nU) ;
veln_lM = zeros(ntps, nU) ;
gdot_lM = zeros(ntps, nU) ;
HH_rM   = zeros(ntps, nU) ;    % right averaged
divv_rM = zeros(ntps, nU) ;
veln_rM = zeros(ntps, nU) ;
gdot_rM = zeros(ntps, nU) ;
HH_dM   = zeros(ntps, nU) ;    % dorsal averaged
divv_dM = zeros(ntps, nU) ;
veln_dM = zeros(ntps, nU) ;
gdot_dM = zeros(ntps, nU) ;
HH_vM   = zeros(ntps, nU) ;    % ventral averaged
divv_vM = zeros(ntps, nU) ;
veln_vM = zeros(ntps, nU) ;
gdot_vM = zeros(ntps, nU) ;

% Build timepoint list so that we first do every 10, then fill in details
lastIdx = length(QS.xp.fileMeta.timePoints) - 1 ;
coarseIdx = 1:10:lastIdx ;
fineIdx = setdiff(1:lastIdx, coarseIdx) ;
allIdx = [66, coarseIdx, fineIdx ] ;
tp2do = QS.xp.fileMeta.timePoints(allIdx) ;

% Output directory is inside metricKinematics dir
outdir = fullfile(mKDir, 'measurements') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Unit definitions for axis labels
unitstr = [ '[1/' QS.timeunits ']' ];
Hunitstr = [ '[1/' QS.spaceunits ']' ];
vunitstr = [ '[' QS.spaceunits '/' QS.timeunits ']' ];
    
% Compute or load all timepoints
for tp = tp2do
    close all
    disp(['t = ' num2str(tp)])
    tidx = QS.xp.tIdx(tp) ;

    % Check for timepoint measurement on disk
    Hfn = fullfile(outdir, sprintf('HH_series_%06d.mat', tp))   ;
    efn = fullfile(outdir, sprintf('gdot_series_%06d.mat', tp)) ;
    dfn = fullfile(outdir, sprintf('divv_series_%06d.mat', tp)) ;
    nfn = fullfile(outdir, sprintf('veln_series_%06d.mat', tp)) ;
    H2vnfn = fullfile(outdir, sprintf('H2vn_series_%06d.mat', tp)) ;

    % Load timeseries measurements
    load(Hfn, 'HH', 'HH_ap', 'HH_l', 'HH_r', 'HH_d', 'HH_v')
    load(efn, 'gdot', 'gdot_ap', 'gdot_l', 'gdot_r', 'gdot_d', 'gdot_v')
    load(dfn, 'divv', 'divv_ap', 'divv_l', 'divv_r', 'divv_d', 'divv_v')
    load(nfn, 'veln', 'veln_ap', 'veln_l', 'veln_r', 'veln_d', 'veln_v') 
    load(H2vnfn, 'H2vn', 'H2vn_ap', 'H2vn_l', 'H2vn_r', 'H2vn_d', 'H2vn_v') 

    % separate 2d/3d data
    H2d = HH ;
    H3d = HH(:, 1:nV-1) ;
    gdot2d = gdot ;
    gdot3d = gdot(:, 1:nV-1);
    divv2d = divv ;
    divv3d = divv(:, 1:nV-1) ;
    veln2d = veln ;
    veln3d = veln(:, 1:nV-1) ;
    H2vn2d = H2vn ;
    H2vn3d = H2vn(:, 1:nV-1) ;
    
    %% Plot results
    % operational plotting options
    pOptions.overwrite = overwrite ;
    pOptions.plot_flows = plot_flows ;
    pOptions.plot_Hgdot = plot_Hgdot ;
    pOptions.plot_factors = plot_factors ;
    % parameter plotting options
    pOptions.doubleResolution = doubleResolution; 
    pOptions.lambda = lambda ;
    pOptions.lambda_err = lambda_err ;
    pOptions.lambda_mesh = lambda_mesh ;
    pOptions.H2vn2d = H2vn2d ;
    pOptions.divv2d = divv2d ;
    pOptions.gdot2d = gdot2d ;
    pOptions.veln2d = veln2d ;
    pOptions.H2d = H2d ;
    pOptions.H2vn3d = H2vn3d ;
    pOptions.divv3d = divv3d ;
    pOptions.gdot3d = gdot3d ;
    pOptions.veln3d = veln3d ;
    pOptions.H3d = H3d ;
    pOptions.cutMesh = cutMesh ;
    pOptions.mesh = mesh ;
    pOptions.tp = tp ;
    pOptions.climit = climit ;
    pOptions.climit_err = climit ;
    pOptions.climit_veln = climit_veln ;
    pOptions.climit_H = climit_H ;
    QS.plotMetricKinematicsTimePoint(pOptions)
    
    %% Store in matrices
    % dv averaged
    HH_apM(tidx, :) = HH_ap ;
    gdot_apM(tidx, :) = gdot_ap ;
    divv_apM(tidx, :) = divv_ap ;
    veln_apM(tidx, :) = veln_ap ;
    H2vn_apM(tidx, :) = H2vn_ap ;

    % left quarter
    HH_lM(tidx, :) = HH_l ;
    gdot_lM(tidx, :) = gdot_l ;
    divv_lM(tidx, :) = divv_l ;
    veln_lM(tidx, :) = veln_l ;
    H2vn_lM(tidx, :) = H2vn_l ;

    % right quarter
    HH_rM(tidx, :) = HH_r ;
    gdot_rM(tidx, :) = gdot_r ;
    divv_rM(tidx, :) = divv_r ;
    veln_rM(tidx, :) = veln_r ;
    H2vn_rM(tidx, :) = H2vn_r ;

    % dorsal quarter
    HH_dM(tidx, :) = HH_d ;
    gdot_dM(tidx, :) = gdot_d ;
    divv_dM(tidx, :) = divv_d ;
    veln_dM(tidx, :) = veln_d ;
    H2vn_dM(tidx, :) = H2vn_d ;

    % ventral quarter
    HH_vM(tidx, :) = HH_v ;
    gdot_vM(tidx, :) = gdot_v ;
    divv_vM(tidx, :) = divv_v ;
    veln_vM(tidx, :) = veln_v ;
    H2vn_vM(tidx, :) = H2vn_v ;
end

%% Store kymograph data in cell arrays
HHsK = {HH_apM, HH_lM, HH_rM, HH_dM, HH_vM} ;
gdotsK = {gdot_apM, gdot_lM, gdot_rM, gdot_dM, gdot_vM} ;
divvsK = {divv_apM, divv_lM, divv_rM, divv_dM, divv_vM} ;
velnsK = {veln_apM, veln_lM, veln_rM, veln_dM, veln_vM} ;
H2vnsK = {H2vn_apM, H2vn_lM, H2vn_rM, H2vn_dM, H2vn_vM} ;

%% Now plot different measured quantities as kymographs
if plot_kymographs
    % Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
    dvDir = fullfile(mKDir, 'avgDV') ;
    lDir = fullfile(mKDir, 'avgLeft') ;
    rDir = fullfile(mKDir, 'avgRight') ;
    dDir = fullfile(mKDir, 'avgDorsal') ;
    vDir = fullfile(mKDir, 'avgVentral') ;
    outdirs = {dvDir, lDir, rDir, dDir, vDir} ;
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    for qq = 1:length(outdirs)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end

        % Unpack what to plot (averaged kymographs, vary averaging region)
        HHK = HHsK{qq} ;
        gdotK = gdotsK{qq} ;
        divvK = divvsK{qq} ;
        velnK = velnsK{qq} ;
        H2vnK = H2vnsK{qq} ;
        m2plot = {gdotK, HHK, divvK, velnK, H2vnK} ;
        titles = {'$\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$',...
            'mean curvature, $H$', ...
            'divergence of flow, $\nabla \cdot \mathbf{v}$', ...
            'normal velocity, $v_n$', ...
            'normal motion, $v_n 2 H$'} ;
        labels = {['$\textrm{Tr}[g^{-1}\dot{g}]$ ' unitstr], ...
            ['mean curvature, $H$ ' Hunitstr], ...
            ['$\nabla \cdot \mathbf{v}$ ' unitstr], ...
            ['normal velocity, $v_n$ ' vunitstr] , ...
            ['normal motion, $v_n 2 H $ ' unitstr]} ;
        names = {'gdot', 'HH', 'divv', 'veln', 'H2vn'} ;
        climits = [climit, climit_H, climit, climit_veln, climit_err] ;

        %% Plot gdot/HH/divv/veln/H2vn DV-averaged kymograph
        for pp = 1:length(m2plot)
            
            % Check if images already exist on disk
            fn = fullfile(odir, [ names{pp} '.png']) ;
            fn_zoom = fullfile(odir, [names{pp} '_zoom_early.png']) ;
            
            if ~exist(fn, 'file') || ~exist(fn_zoom, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                colormap bwr
                imagesc((1:nU)/nU, tps, m2plot{pp})
                caxis([-climits(pp), climits(pp)])
                % Add folds to plot
                hold on;
                fons1 = max(1, fons(1)) ;
                fons2 = max(1, fons(2)) ;
                fons3 = max(1, fons(3)) ;
                plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
                plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
                plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

                % title and save
                title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeunits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
                fn = fullfile(odir, [ names{pp} '.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   

                % Zoom in on small values
                caxis([-climits(pp)/3, climits(pp)/3])
                fn = fullfile(odir, [names{pp} '_zoom.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
                % Zoom in on early times
                ylim([min(tps), max(fons) + 10])
                caxis([-climits(pp)/3, climits(pp)/3])
                fn = fullfile(odir, [names{pp} '_zoom_early.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
            end
        end
    end
end


%% Kymographs of cumulative sums along pathlines
if plot_kymographs_cumsum
    % Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
    dvDir = fullfile(mKDir, 'avgDV') ;
    lDir = fullfile(mKDir, 'avgLeft') ;
    rDir = fullfile(mKDir, 'avgRight') ;
    dDir = fullfile(mKDir, 'avgDorsal') ;
    vDir = fullfile(mKDir, 'avgVentral') ;
    outdirs = {dvDir, lDir, rDir, dDir, vDir} ;
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    for qq = 1:length(outdirs)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end

        % Unpack what to plot (averaged kymographs, vary averaging region)
        HHK = HHsK{qq} ;
        gdotK = cumsum(gdotsK{qq}, 1) ;
        divvK = cumsum(divvsK{qq}, 1) ;
        velnK = cumsum(velnsK{qq}, 1) ;
        H2vnK = cumsum(H2vnsK{qq}, 1) ;
        m2plot = {gdotK, divvK, velnK, H2vnK} ;
        titles = {'$\int_0^t\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$',...
            'divergence of flow, $\nabla \cdot \mathbf{v}$', ...
            'normal velocity, $\int_0^tv_n$', ...
            'normal motion, $\int_0^t v_n 2 H$'} ;
        labels = {['$\int_0^t\textrm{Tr}[g^{-1}\dot{g}]$ ' unitstr], ...
            ['$\int_0^t \nabla \cdot \mathbf{v}$ ' unitstr], ...
            ['normal velocity, $\int_0^tv_n$ ' vunitstr] , ...
            ['normal motion, $\int_0^t v_n 2 H $ ' unitstr]} ;
        names = {'Igdot', 'Idivv', 'Iveln', 'IH2vn'} ;
        climits = [climit, climit, climit_veln, climit] ;
        climits = climits * 3; 
        
        %% Plot gdot/HH/divv/veln/H2vn DV-averaged kymograph
        for pp = 1:length(m2plot)
            % Check if images already exist on disk
            fn = fullfile(odir, [ names{pp} '.png']) ;
            fn_zoom = fullfile(odir, [names{pp} '_zoom_early.png']) ;
            if ~exist(fn, 'file') || ~exist(fn_zoom, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                colormap bwr
                imagesc((1:nU)/nU, tps, m2plot{pp})
                caxis([-climits(pp), climits(pp)])
                % Add folds to plot
                hold on;
                fons1 = max(1, fons(1)) ;
                fons2 = max(1, fons(2)) ;
                fons3 = max(1, fons(3)) ;
                plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
                plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
                plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

                % title and save
                title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeunits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   

                % Zoom in on small values
                caxis([-climits(pp)/3, climits(pp)/3])
                fn = fullfile(odir, [names{pp} '_zoom.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
                % Zoom in on early times
                ylim([min(tps), max(fons) + 10])
                caxis([-climits(pp)/3, climits(pp)/3])
                disp(['saving ', fn_zoom])
                export_fig(fn_zoom, '-png', '-nocrop', '-r200')   
            end
        end
    end
end

%% Metric Kinematic Correlations
% Plot both all time and select times
timeSpans = {tps, tps(tps < max(fons) + 11)} ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot correlation between terms div(v) and 2Hvn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
padN = round(0.1 * nU) ;
cols = padN:nU - padN ;
if plot_correlations
    corrDir = fullfile(mKDir, 'correlations') ;
    if ~exist(corrDir, 'dir')
        mkdir(corrDir)
    end
    
    for sigma = 0:3
        outputFileNames = {fullfile(corrDir, ...
            sprintf('correlation_sigma%02d_alltime_div_2Hvn', sigma)), ...
            fullfile(corrDir, ...
            sprintf('correlation_sigma%02d_earlytimes_div_2Hvn', sigma))} ;
        alphaVal = 0.6 ;
        sz = 10 ;
        cmap = parula ;
        close all
        set(gcf, 'visible', 'off')
        % Consider each timespan (early or entire series)
        for tspanIdx = 1:2
            fnout = outputFileNames{tspanIdx} ;
            timeSpan_i = timeSpans{tspanIdx} ;
            ntspan = length(timeSpan_i) ;
            titles = {'left lateral', 'right lateral', 'dorsal', 'ventral'} ;
            markers = QS.plotting.markers ;
            colors = mapValueToColor(1:ntspan, [1, ntspan], cmap) ;
            close all
            sphCollection = cell(4, 1) ;
            sposCollection = cell(4, 1) ;
            sphCollection2 = cell(4, 1) ;
            sposCollection2 = cell(4, 1) ;
            sphCollection3 = cell(4, 1) ;
            sposCollection3 = cell(4, 1) ;
            for qq = 1:4  % consider left, right, dorsal, ventral
                disp(['qq = ', num2str(qq), ': ', titles{qq}])
                divv = divvsK{qq + 1}(:, cols) ;    
                H2vn = H2vnsK{qq + 1}(:, cols) ;

                % Optional: smooth here
                if sigma > 0
                    divv = imgaussfilt(divv, sigma);            
                    H2vn = imgaussfilt(H2vn, sigma);  
                end
                
                % Check the smoothing on kymographs
                figure(2) ;
                sphCollection2{qq} = subplot(2, 2, qq) ;
                imagesc(cols/nU, tps, divv); 
                caxis([-climit, climit])
                colormap(bwr256)
                figure(3) ;
                sphCollection3{qq} = subplot(2, 2, qq) ;
                imagesc(cols/nU, tps, H2vn); 
                caxis([-climit, climit])
                colormap(bwr256)

                figure(1) ;
                sphCollection{qq} = subplot(2, 2, qq) ;
                for row = 1:ntspan
                    disp(['row = ', num2str(row)])
                    scatter(divv(row, :), H2vn(row, :), sz, ...
                        markers{qq}, 'MarkerFaceColor', 'none', ...
                        'MarkerEdgeColor', colors(row, :), ...
                        'MarkerEdgeAlpha', alphaVal) ;
                    hold on ;
                end

                % Label the x axis if on the bottom row
                if qq > 2
                    figure(1)
                    xlabel(['$\nabla \cdot \bf{v}_\parallel$ ' unitstr], ...
                            'Interpreter', 'Latex') ;
                    figure(2)
                    xlabel(['ap position, $\zeta/L$'], ...
                        'Interpreter', 'Latex') ;
                    figure(3)
                    xlabel(['ap position, $\zeta/L$'], ...
                        'Interpreter', 'Latex') ;
                end
                figure(1)
                axis equal
                % Add dashed y=x line
                xlims = get(gca, 'xlim') ;
                ylims = get(gca, 'ylim') ;
                xlim([max(-2 * climit, xlims(1)), min(2 * climit, xlims(2))])
                ylim([max(-2 * climit, ylims(1)), min(2 * climit, ylims(2))])
                xlims = get(gca, 'xlim') ;
                ylims = get(gca, 'ylim') ;
                leftdot = max(xlims(1), ylims(1)) ;
                rightdot = min(xlims(2), ylims(2)) ;
                plot([leftdot, rightdot], [leftdot, rightdot], 'k--')

                % Label the y axis if on the left column
                ylabel(['$2Hv_n$ ' unitstr], 'Interpreter', 'Latex') ;
                title(titles{qq}, 'Interpreter', 'Latex')
                
                figure(2)
                ylabel(['time [' QS.timeunits, ']'], 'Interpreter', 'Latex') ;
                title(titles{qq}, 'Interpreter', 'Latex')
                figure(3)
                ylabel(['time [' QS.timeunits, ']'], 'Interpreter', 'Latex') ;
                title(titles{qq}, 'Interpreter', 'Latex')

                % Grab axis position
                sposCollection{qq} = get(sphCollection{qq}, 'Position');
                sposCollection2{qq} = get(sphCollection2{qq}, 'Position');
                sposCollection3{qq} = get(sphCollection3{qq}, 'Position');
            end

            % Move subplots left a bit for colorbar space
            for qq = 1:length(sphCollection)
                spos = sposCollection{qq} ;
                wh = min(spos(3)-0.05, spos(4)) ;
                if mod(qq, 2) == 1
                    set(sphCollection{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                    set(sphCollection2{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                    set(sphCollection3{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                else
                    set(sphCollection{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                    set(sphCollection2{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                    set(sphCollection3{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                end
            end

            % master titles (suptitles)
            figure(1) ;
            sgtitle(['$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$, ',...
                '$\sigma=$', num2str(sigma), ' ', QS.timeunits], ...
                'Interpreter', 'Latex') ;
            
            figure(2) ;
            sgtitle(['$\nabla \cdot \bf{v}_\parallel,$ $\sigma=$', ...
                num2str(sigma), ' ', QS.timeunits], ...
                    'Interpreter', 'Latex') ;
            figure(3) ;
            sgtitle(['$2Hv_n$ $\sigma=$', ...
                num2str(sigma), ' ', QS.timeunits], 'Interpreter', 'Latex') ;
                
            % Add colorbar
            figure(1) ;
            c = colorbar('Position',[.9 .333 .02 .333]) ;
            % Make colorbar share the alpha of the image
            % Manually flush the event queue and force MATLAB to render the colorbar
            % necessary on some versions
            drawnow
            % Get the color data of the object that correponds to the colorbar
            cdata = c.Face.Texture.CData;
            % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
            cdata(end,:) = uint8(alphaVal * cdata(end,:));
            % Ensure that the display respects the alpha channel
            c.Face.Texture.ColorType = 'truecoloralpha';
            % Update the color data with the new transparency information
            c.Face.Texture.CData = cdata;
            c.Label.Interpreter = 'Latex' ;
            c.Label.String = ['time [' QS.timeunits ']'] ;
            c.Ticks = [0, 1] ;
            c.TickLabels = [tps(1), max(timeSpan_i)] ;
            
            figure(2) ;
            c = colorbar('Position',[.9 .333 .02 .333]) ;
            figure(3) ;
            c = colorbar('Position',[.9 .333 .02 .333]) ;
            
            % Save figure
            figure(1)
            saveas(gcf, [fnout '.png']) ;
            figure(2)
            saveas(gcf, [fnout '_kymo_divv.png']) ;
            figure(3)
            saveas(gcf, [fnout '_kymo_H2vn.png']) ;
            close all
            set(gcf, 'visible', 'off')
        end
        disp('done with correlation plots betweeen divv and H2vn')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot correlations between gdot and each term in the sum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_gdot_correlations
    outputFileNames = cell(2, 1) ;
    outputFileNames{1} = {fullfile(mKDir, 'correlation_alltime_div_gdot'), ...
                       fullfile(mKDir, 'correlation_earlytimes_div_gdot')} ;
    outputFileNames{2} = {fullfile(mKDir, 'correlation_alltime_2Hvn_gdot'), ...
                       fullfile(mKDir, 'correlation_earlytimes_2Hvn_gdot')} ;
    alphaVal = 0.6 ;
    sz = 10 ;
    cmap = parula ;
    close all
    set(gcf, 'visible', 'off')
    for pairIdx = 1:2
        for tspanIdx = 1:2
            fnout = outputFileNames{pairIdx}{tspanIdx} ;
            timeSpan_i = timeSpans{tspanIdx} ;
            ntspan = length(timeSpan_i) ;
            titles = {'left lateral', 'right lateral', 'dorsal', 'ventral'} ;
            markers = QS.plotting.markers ;
            colors = mapValueToColor(1:ntspan, [1, ntspan], cmap) ;
            close all
            cols = round(nV * [0.2, 0.85]) ;
            sphCollection = cell(4, 1) ;
            sposCollection = cell(4, 1) ;
            for qq = 1:4  % consider left, right, dorsal, ventral
                disp(['qq = ', num2str(qq), ': ', titles{qq}])
                if pairIdx == 1
                    divv = divvsK{qq + 1} ;
                    gdot = gdotsK{qq + 1} ;
                    sphCollection{qq} = subplot(2, 2, qq) ;
                    for row = 1:ntspan
                        disp(['row = ', num2str(row)])
                        scatter(divv(row, cols), gdot(row, cols), sz, ...
                            markers{qq}, 'MarkerFaceColor', 'none', ...
                            'MarkerEdgeColor', colors(row, :), ...
                            'MarkerEdgeAlpha', alphaVal) ;
                        hold on ;
                    end

                    % Label the x axis if on the bottom row
                    if qq > 2
                        xlabel(['$\nabla \cdot \bf{v}_\parallel$ ' unitstr], ...
                                'Interpreter', 'Latex') ;
                    end
                    axis equal
                    % Add dashed y=x line
                    xlims = get(gca, 'xlim') ;
                    ylims = get(gca, 'ylim') ;
                    leftdot = max(xlims(1), ylims(1)) ;
                    rightdot = min(xlims(2), ylims(2)) ;
                    plot([leftdot, rightdot], [leftdot, rightdot], 'k--')
                else
                    H2vn = H2vnsK{qq + 1} ;
                    gdot = gdotsK{qq + 1} ;
                    sphCollection{qq} = subplot(2, 2, qq) ;
                    for row = 1:ntspan
                        disp(['row = ', num2str(row)])
                        scatter(H2vn(row, cols), gdot(row, cols), sz, ...
                            markers{qq}, 'MarkerFaceColor', 'none', ...
                            'MarkerEdgeColor', colors(row, :), ...
                            'MarkerEdgeAlpha', alphaVal) ;
                        hold on ;
                    end

                    % Label the x axis if on the bottom row
                    if qq > 2
                        xlabel(['$2Hv_n$ ' unitstr], 'Interpreter', 'Latex') ;
                    end
                    axis equal
                    % Add dashed y=x line
                    xlims = get(gca, 'xlim') ;
                    ylims = get(gca, 'ylim') ;
                    leftdot = max(xlims(1), -ylims(2)) ;
                    rightdot = min(xlims(2), -ylims(1)) ;
                    plot([leftdot, rightdot], [-leftdot, -rightdot], 'k--')
                end

                % Label the y axis if on the left column
                if qq == 1 || qq == 3
                    ylabel(['$\textrm{Tr}[g^{-1} \dot{g}]$ ' unitstr], ...
                            'Interpreter', 'Latex')
                end
                title(titles{qq}, 'Interpreter', 'Latex')

                % Grab axis position
                sposCollection{qq} = get(sphCollection{qq}, 'Position');
            end

            % Move subplots left a bit for colorbar space
            for qq = 1:length(sphCollection)
                spos = sposCollection{qq} ;
                wh = min(spos(3)-0.05, spos(4)) ;
                if mod(qq, 2) == 1
                    set(sphCollection{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                else
                    set(sphCollection{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                end
            end

            % Add colorbar
            c = colorbar('Position',[.9 .333 .02 .333]) ;
            % Make colorbar share the alpha of the image
            % Manually flush the event queue and force MATLAB to render the colorbar
            % necessary on some versions
            drawnow
            % Get the color data of the object that correponds to the colorbar
            cdata = c.Face.Texture.CData;
            % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
            cdata(end,:) = uint8(alphaVal * cdata(end,:));
            % Ensure that the display respects the alpha channel
            c.Face.Texture.ColorType = 'truecoloralpha';
            % Update the color data with the new transparency information
            c.Face.Texture.CData = cdata;
            c.Label.Interpreter = 'Latex' ;
            c.Label.String = ['time [' QS.timeunits ']'] ;
            c.Ticks = [0, 1] ;
            c.TickLabels = [tps(1), max(timeSpan_i)] ;

            % Save figure
            saveas(gcf, [fnout '.png']) ;
            saveas(gcf, [fnout '.pdf']) ;
            close all
            set(gcf, 'visible', 'off')
        end
    end
    disp('done')
end


%% Metric Kinematics -- decompose into isotropic and other component
if plot_gdot_decomp
    % Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
    dvDir = fullfile(mKDir, 'avgDV') ;
    lDir = fullfile(mKDir, 'avgLeft') ;
    rDir = fullfile(mKDir, 'avgRight') ;
    dDir = fullfile(mKDir, 'avgDorsal') ;
    vDir = fullfile(mKDir, 'avgVentral') ;
    outdirs = {dvDir, lDir, rDir, dDir, vDir} ;
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    gdotK0 = gdotsK{1} ;
    isogrowth = sum(gdotK0, 2) ;
    
    %% Plot as 1d Curve
    plot(tps, isogrowth); 
    xlabel(['time [' QS.timeunits ']'], 'Interpreter', 'Latex')
    ylabel('Isotropic component of growth') 
    saveas(gcf, fullfile(odir, 'average_growth.png'))  
    
    %% Plot quadrant contributions
    for qq = 1:length(outdirs)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end
        % Unpack what to plot (averaged kymographs, vary averaging region)
        HHK = HHsK{qq} ;
        gdotK = gdotsK{qq} ;
        divvK = divvsK{qq} ;
        velnK = velnsK{qq} ;
        H2vnK = H2vnsK{qq} ;
        m2plot = {gdotK, HHK, divvK, velnK, H2vnK} ;
        titles = {'$\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$',...
            'mean curvature, $H$', ...
            'divergence of flow, $\nabla \cdot \mathbf{v}$', ...
            'normal velocity, $v_n$', ...
            'normal motion, $v_n 2 H$'} ;
        labels = {['$\textrm{Tr}[g^{-1}\dot{g}]$ ' unitstr], ...
            ['mean curvature, $H$ ' Hunitstr], ...
            ['$\nabla \cdot \mathbf{v}$ ' unitstr], ...
            ['normal velocity, $v_n$ ' vunitstr] , ...
            ['normal motion, $v_n 2 H $ ' unitstr]} ;
        names = {'gdot', 'HH', 'divv', 'veln', 'H2vn'} ;
        climits = [climit, climit_H, climit, climit_veln, climit_err] ;

        %% Plot gdot/HH/divv/veln/H2vn DV-averaged kymograph
        for pp = 1:length(m2plot)            
            fn = fullfile(odir, [ names{pp} '.png']) ;
            if ~exist(fn, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                colormap bwr
                imagesc((1:nU)/nU, tps, m2plot{pp})
                caxis([-climits(pp), climits(pp)])
                % Add folds to plot
                hold on;
                fons1 = max(1, fons(1)) ;
                fons2 = max(1, fons(2)) ;
                fons3 = max(1, fons(3)) ;
                plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
                plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
                plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

                % title and save
                title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeunits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
            end
        end
    end
end
