%% figure_script_dynamicAtlas_gut_correlations
%% dynamicAtlas paper

% Plotting settings

clearvars ax cut timePoints timestamps cmp dmyk overwrite 
clearvars v3dcorr v3ddotprod
cut = 10 ;
timePoints = QS.xp.fileMeta.timePoints ;
timestamps = timePoints(1:end-1) - QS.t0set() ;
cmp = brewermap(256, '*RdBu') ;
cmp = bwr ;
dmyk = 2;
N = 2 ;
M = 3 ;
overwrite = false ; 
close all
ax{1} = subplot(N, M, 1) ;

%% Load 3d velocity correlation
vfn = fullfile(QS.dir.uvCoord, 'autocorrelation_pathline_velocities_vvsm.mat') ;
if ~exist(vfn, 'file') || overwrite 
    tmp = load(sprintf(QS.fileName.pathlines.velocities.vvsm, QS.t0set())) ;
    vvsmM = tmp.vvsmM ;
    v3dcorr = zeros(size(vvsmM, 1), size(vvsmM, 1)) ;
    assert(size(vvsmM, 1) == length(timestamps)) ;
    for tidx = 1:size(vvsmM, 1)
        disp(['tidx = ' num2str(tidx) '/' num2str(size(vvsmM, 1))])
        vi = reshape(squeeze(vvsmM(tidx, :, :)), [QS.nU, QS.nV, 3]) ;
        vi = vi(cut:QS.nU-cut, 1:end-1, :) ;
        vi = reshape(vi, [size(vi, 1)*size(vi, 2), 3]) ;

        % load mesh
        tp = QS.xp.fileMeta.timePoints(tidx) ;
        QS.setTime(tp) ;
        mesh = QS.getCurrentSPCutMeshSmRSC() ;
        faceAreas = doublearea(mesh.v,mesh.f)*0.5;
        vertexAreas = full(sparse( mesh.f, 1, repmat(faceAreas, 1, 3), size(mesh.v,1), 1))/3;
        vertexAreas = reshape(vertexAreas, [QS.nU, QS.nV - 1]) ;
        vertexAreas = vertexAreas(cut:QS.nU-cut, :) ;
        totArea = sum(vertexAreas(:)) ;
        weights = vertexAreas(:) ./ totArea ;
        % weights = ones(size(vi, 1), 1) ;

        for jj = 1:size(vvsmM, 1)
            vj = reshape(squeeze(vvsmM(jj, :, :)), [QS.nU, QS.nV, 3]) ;
            vj = vj(cut:QS.nU-cut, 1:end-1, :) ;
            vj = reshape(vj, [size(vj, 1)*size(vj, 2), 3]) ;
            v3dcorr(tidx, jj) = corr(vi(:), vj(:)) ;
            v3ddotprod(tidx, jj) = sum(weights .* dot(vi, vj, 2) ./ sqrt(dot(vi,vi,2).*dot(vj,vj,2)))/ sum(weights) ;
        end
    end
    
    % Get falloff with time
    v3dcorrROI = v3dcorr(timestamps>=0 & timestamps<=90, timestamps>=0 & timestamps<=90) ;
    v3ddotprodROI = v3ddotprod(timestamps>=0 & timestamps<=90, timestamps>=0 & timestamps<=90) ;
    v3dcorr_falloff = nan(size(v3dcorrROI)) ;
    v3ddprod_falloff = nan(size(v3dcorrROI)) ;
    counts = zeros(size(v3dcorrROI)) ;
    for ii = 1:size(v3dcorrROI, 1)
        strip = v3dcorrROI(ii, ii:end) ;
        nn = length(strip) ;
        v3dcorr_falloff(ii, 1:nn) = strip(:) ;
        strip = v3ddotprodROI(ii, ii:end) ;
        v3ddprod_falloff(ii, 1:nn) = strip(:) ;
        counts(ii, 1:nn) = counts(1:nn) + 1 ;
    end
    v3dcorr_falloff_mean = mean(v3dcorr_falloff, 'omitnan') ;
    v3ddprod_falloff_mean = mean(v3ddprod_falloff, 'omitnan') ;
    v3dcorr_falloff_std = std(v3dcorr_falloff, 'omitnan') ;
    v3ddprod_falloff_std = std(v3ddprod_falloff, 'omitnan') ;
    
    disp('saving the correlations in v3d')
    save(vfn, 'v3dcorr', 'vvsmM', 'v3ddotprod', 'timestamps', ...
        'v3dcorr_falloff', 'v3ddprod_falloff', ...
        'v3dcorr_falloff_mean', 'v3ddprod_falloff_mean', ...
        'v3dcorr_falloff_std', 'v3ddprod_falloff_std', 'counts')
else
    disp('loading v3d correlations from disk')
    load(vfn, 'v3dcorr', 'vvsmM', 'v3ddotprod', 'timestamps', ...
        'v3dcorr_falloff', 'v3ddprod_falloff', ...
        'v3dcorr_falloff_mean', 'v3ddprod_falloff_mean', ...
        'v3dcorr_falloff_std', 'v3ddprod_falloff_std', 'counts')
end
    
% plot it
ax{dmyk} = subplot(N, M, dmyk) ;
imagesc(timestamps, timestamps, v3dcorr)
colormap(cmp)
caxis([-1,1])
colorbar()
axis equal; axis tight
xlabel(['time [' QS.timeUnits, ']'])
ylabel(['time [' QS.timeUnits, ']'])
title('v3dcorr')


set(gcf, 'CurrentAxes', ax{1})
lineProps = {'-','color', colors(dmyk, :)} ;
h1 = shadedErrorBar(1:length(v3dcorr_falloff_mean), v3dcorr_falloff_mean, v3dcorr_falloff_std, 'lineProps', lineProps) ;
xlabel(['time difference, \deltat [' QS.timeUnits, ']'])
ylabel('correlation')
title('v3dcorr')
dmyk = dmyk + 1 ;


ax{dmyk} = subplot(N, M, dmyk) ;
imagesc(timestamps, timestamps, v3ddotprod)
colormap(cmp)
caxis([-1,1])
colorbar()
axis equal; axis tight
xlabel(['time [' QS.timeUnits, ']'])
ylabel(['time [' QS.timeUnits, ']'])
title('v3ddotprod')

set(gcf, 'CurrentAxes', ax{1})
lineProps = {'-','color', colors(dmyk, :)} ;
h2 = shadedErrorBar(1:length(v3ddprod_falloff_mean), v3ddprod_falloff_mean, v3ddprod_falloff_std, 'lineProps', lineProps) ;
xlabel(['time difference, \deltat [' QS.timeUnits, ']'])
ylabel('correlation')
title('v3ddotprod')
dmyk = dmyk + 1 ;

%% v2d smoothed in um per minute
% v2fn = fullfile(QS.dir.uvCoord, 'autocorrelation_pathline_velocities_v2dsmum.mat') ;
% if ~exist(v2fn, 'file')
%     tmp = load(sprintf(QS.fileName.pathlines.velocities.vv2dsmum, QS.t0set())) ;
%     v2dsmumM = tmp.vv2dsmMum ;
%     v2dsmumcorr = zeros(size(v2dsmumM, 1), size(v2dsmumM, 1)) ;
%     assert(size(v2dsmumM_faces, 1) == length(timestamps)) ;
%     for tidx = 1:size(v2dsmumM_faces, 1)
%         disp(['tidx = ' num2str(tidx) '/' num2str(size(v2dsmumM_faces, 1))])
%         tp = QS.xp.fileMeta.timePoints(tidx) ;
%         % load mesh
%         % QS.setTime(tp) ;
%         % mesh = QS.getCurrentSPCutMeshSmRSC() ;
%         % [V2F, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
%         % v2dsmumM(tidx, :, :) = F2V * squeeze(v2dsmumM_faces(tidx, :, :));
%         
%         vi = reshape(squeeze(v2dsmumM(tidx, :, :)), [QS.nU, QS.nV, 3]) ;
%         vi = vi(cut:QS.nU-cut, 1:end-1, :) ;
%         for jj = 1:size(v2dsmumM, 1)
%             vj = reshape(squeeze(v2dsmumM(jj, :, :)), [QS.nU, QS.nV, 3]) ;
%             vj = vj(cut:QS.nU-cut, 1:end-1, :) ;
%             v2dsmumcorr(tidx, jj) = corr(vi(:), vj(:)) ;
%         end
%     end
%     save(v2fn, 'v2dsmumcorr', 'v2dsmumM', 'timestamps')
% else
%     disp('loading v2dsmum correlations from disk')
%     load(v2fn, 'v2dsmumcorr', 'v2dsmumM', 'timestamps')
% end
%     
% % plot it
% ax{dmyk} = subplot(N, M, dmyk) ;
% imagesc(timestamps, timestamps, v2dsmumcorr)
% colormap(cmp)
% caxis([-1,1])
% colorbar()
% axis equal; axis tight
% xlabel(['time [' QS.timeUnits, ']'])
% ylabel(['time [' QS.timeUnits, ']'])
% dmyk = dmyk + 1 ;

%% correlation of pathline out-of-plane deformation
hfn = fullfile(QS.dir.uvCoord, 'autocorrelation_pathline_H2vn.mat') ;
if ~exist(hfn, 'file') || overwrite
    fn = fullfile(sprintf(QS.dir.metricKinematics.pathline.measurements, QS.t0set()), ...
        'H2vn_pathline0123_%06d.mat') ;
    h2vnM = zeros(length(timestamps), QS.nU-2*cut+1, QS.nV-1) ;
    for tidx = 1:length(timestamps)
        tp = QS.xp.fileMeta.timePoints(tidx) ;
        h2vn = load(sprintf(fn, tp)) ;
        h2vnM(tidx, :, :) = h2vn.H2vn(cut:end-cut, 1:QS.nV-1) ;
    end

    % build correlation matrix
    h2vncorr = zeros(length(timestamps), length(timestamps)) ;
    for tidx = 1:length(timestamps)
        disp(['tidx = ' num2str(tidx) '/' num2str(length(timestamps))])
        h2vni = squeeze(h2vnM(tidx, :, :)) ;
        for jj = 1:length(timestamps)    
            h2vnj = squeeze(h2vnM(jj, :, :)) ;
            h2vncorr(tidx, jj) = corr(h2vni(:), h2vnj(:)) ;
        end
    end
    
    
    
    % Get falloff with time
    h2vncorrROI = h2vncorr(timestamps>=0 & timestamps<=90, timestamps>=0 & timestamps<=90) ;
    h2vncorr_falloff = nan(size(h2vncorrROI)) ;
    counts = zeros(size(h2vncorrROI)) ;
    for ii = 1:size(h2vncorrROI, 1)
        strip = h2vncorrROI(ii, ii:end) ;
        nn = length(strip) ;
        h2vncorr_falloff(ii, 1:nn) = strip(:) ;
        counts(ii, 1:nn) = counts(1:nn) + 1 ;
    end
    h2vncorr_falloff_mean = mean(h2vncorr_falloff, 'omitnan') ;
    h2vncorr_falloff_std = std(h2vncorr_falloff, 'omitnan') ;
    
    disp('saving h2vn correlation')
    save(hfn, 'h2vncorr', 'h2vnM', 'timestamps',...
        'h2vncorr_falloff', 'h2vncorr_falloff_mean', 'h2vncorr_falloff_std', 'counts')
else
    disp('loading H2vn correlations from disk')
    load(hfn, 'h2vncorr', 'h2vnM', 'timestamps', ...
        'h2vncorr_falloff', 'h2vncorr_falloff_mean', 'h2vncorr_falloff_std', 'counts')
end

% plot it
ax{dmyk} = subplot(N, M, dmyk) ;
imagesc(timestamps, timestamps, h2vncorr)
try
    colormap(vik)
catch
    colormap(bwr)
    % colormap(brewermap(256, '*RdBl'))
end
caxis([-1,1])
colorbar()
axis equal; axis tight
xlabel(['time [' QS.timeUnits, ']'])
ylabel(['time [' QS.timeUnits, ']'])
title('h2vn')

set(gcf, 'CurrentAxes', ax{1})
lineProps = {'-','color', colors(dmyk, :)} ;
h3 = shadedErrorBar(1:length(h2vncorr_falloff_mean), h2vncorr_falloff_mean, h2vncorr_falloff_std, 'lineProps', lineProps) ;
xlabel(['time difference, \deltat [' QS.timeUnits, ']'])
ylabel('correlation')
title('h2vn')

dmyk = dmyk + 1 ;

%% Correlation of divergence
dfn = fullfile(QS.dir.uvCoord, 'autocorrelation_pathline_divv.mat') ;
if ~exist(dfn, 'file') || overwrite
    fn = fullfile(sprintf(QS.dir.metricKinematics.pathline.measurements, QS.t0set()), ...
        'divv_pathline0123_%06d.mat') ;
    divvM = zeros(length(timestamps), QS.nU-2*cut+1, QS.nV-1) ;
    for tidx = 1:length(timestamps)
        tp = QS.xp.fileMeta.timePoints(tidx) ;
        divv = load(sprintf(fn, tp)) ;
        divvM(tidx, :, :) = divv.divv(cut:end-cut, 1:QS.nV-1) ;
    end

    % build correlation matrix
    divvcorr = zeros(length(timestamps), length(timestamps)) ;
    for tidx = 1:length(timestamps)
        disp(['tidx = ' num2str(tidx) '/' num2str(length(timestamps))])
        divvi = squeeze(divvM(tidx, :, :)) ;
        for jj = 1:length(timestamps)    
            divvj = squeeze(divvM(jj, :, :)) ;
            divvcorr(tidx, jj) = corr(divvi(:), divvj(:)) ;
        end
    end
    
    % Get falloff with time
    divvcorrROI = divvcorr(timestamps>=0 & timestamps<=90, timestamps>=0 & timestamps<=90) ;
    divvcorr_falloff = nan(size(divvcorrROI)) ;
    counts = zeros(size(divvcorrROI)) ;
    for ii = 1:size(divvcorrROI, 1)
        strip = divvcorrROI(ii, ii:end) ;
        nn = length(strip) ;
        divvcorr_falloff(ii, 1:nn) = strip(:) ;
        counts(ii, 1:nn) = counts(1:nn) + 1 ;
    end
    divvcorr_falloff_mean = mean(divvcorr_falloff, 'omitnan') ;
    divvcorr_falloff_std = std(divvcorr_falloff, 'omitnan') ;
    
    disp('saving divv')
    save(dfn, 'divvcorr', 'divvM', 'timestamps', ...
        'divvcorr_falloff', 'divvcorr_falloff_mean', 'divvcorr_falloff_std', 'counts')
else
    disp('Loading divv from disk')
    load(dfn, 'divvcorr', 'divvM', 'timestamps')
    
end

% plot it
ax{dmyk} = subplot(N, M, dmyk) ;
imagesc(timestamps, timestamps, divvcorr)
colormap(cmp)
caxis([-1,1])
colorbar()
axis equal; axis tight
xlabel(['time [' QS.timeUnits, ']'])
ylabel(['time [' QS.timeUnits, ']'])
title('divv')

set(gcf, 'CurrentAxes', ax{1})
lineProps = {'-','color', colors(dmyk, :)} ;
h4 = shadedErrorBar(1:length(divvcorr_falloff_mean), divvcorr_falloff_mean, divvcorr_falloff_std, 'lineProps', lineProps) ;
xlabel(['time difference, \deltat [' QS.timeUnits, ']'])
ylabel('correlation')
title('divv')

dmyk = dmyk + 1;

%% correlation of pathline curl
% options = struct('lambda', 0.002, 'lambda_mesh', 0) ;
% options.overwrite = true ;
% QS.measurePathlineVorticity(options) ;
% QS.plotPathlineVorticity() ;

ofn = fullfile(QS.dir.uvCoord, 'autocorrelation_pathline_vorticity.mat') ;
if ~exist(ofn, 'file') || overwrite
    fn = fullfile(sprintf(QS.dir.pathlines.vorticity, QS.t0set()), ...
        'vorticity_pathline0123_%06d.mat') ;
    omegaM = zeros(length(timestamps), QS.nU-2*cut+1, QS.nV-1) ;
    for tidx = 1:length(timestamps)
        tp = QS.xp.fileMeta.timePoints(tidx) ;
        omi = load(sprintf(fn, tp)) ;
        omegaM(tidx, :, :) = omi.vorticity(cut:end-cut, 1:QS.nV-1) ;
    end
    
    % Despeckle by running median filter
    % first augment the matrix
    nn = 10 ;
    omegaMbig = zeros(size(omegaM, 1), size(omegaM, 2), size(omegaM, 3)+nn*2) ;
    omegaMbig(:, :, nn+1:nn+nV-1) = omegaM ;
    omegaMbig(:, :, nn+nV:end) = omegaM(:, :, 1:nn) ;
    omegaMbig(:, :, 1:nn) = omegaM(:, :, end-nn+1:end) ;
    omegaMbig = medfilt3(omegaMbig, [3,3,3]) ;
    omegaMbig = medfilt3(omegaMbig, [3,3,3]) ;
    omegaM = omegaMbig(:, :, nn+1:nn+nV-1) ;

    % Check it
    % for tidx = 1:length(timestamps)
    %     imagesc(squeeze(omegaM(tidx,:,:))')
    %     axis equal
    %     caxis([-0.2, 0.2])
    %     colormap(bwr)
    %     pause(0.01)
    % end

    % build correlation matrix
    omegacorr = zeros(length(timestamps), length(timestamps)) ;
    for tidx = 1:length(timestamps)
        disp(['tidx = ' num2str(tidx) '/' num2str(length(timestamps))])
        omegai = squeeze(omegaM(tidx, :, :)) ;
        for jj = 1:length(timestamps)    
            omegaj = squeeze(omegaM(jj, :, :)) ;
            omegacorr(tidx, jj) = corr(omegai(:), omegaj(:)) ;
        end
    end
    
    % Get falloff with time
    omegacorrROI = omegacorr(timestamps>=0 & timestamps<=90, timestamps>=0 & timestamps<=90) ;
    omegacorr_falloff = nan(size(omegacorrROI)) ;
    counts = zeros(size(omegacorrROI)) ;
    for ii = 1:size(omegacorrROI, 1)
        strip = omegacorrROI(ii, ii:end) ;
        nn = length(strip) ;
        omegacorr_falloff(ii, 1:nn) = strip(:) ;
        counts(ii, 1:nn) = counts(1:nn) + 1 ;
    end
    omegacorr_falloff_mean = mean(omegacorr_falloff, 'omitnan') ;
    omegacorr_falloff_std = std(omegacorr_falloff, 'omitnan') ;
    
    disp('saving omega....')
    save(ofn, 'omegacorr', 'omegaM', 'timestamps', ...
        'omegacorr_falloff', 'omegacorr_falloff_mean', 'omegacorr_falloff_std', 'counts')
else
    disp('Loading omega from disk')
    load(ofn, 'omegacorr', 'omegaM', 'timestamps', ...
        'omegacorr_falloff', 'omegacorr_falloff_mean', 'omegacorr_falloff_std', 'counts')
end

% plot it
ax{dmyk} = subplot(N, M, dmyk) ;
imagesc(timestamps, timestamps, omegacorr)
colormap(cmp)
caxis([-1,1])
colorbar()
axis equal; axis tight
xlabel(['time [' QS.timeUnits, ']'])
ylabel(['time [' QS.timeUnits, ']'])
title('omega')

set(gcf, 'CurrentAxes', ax{1})
lineProps = {'-','color', colors(dmyk, :)} ;
h5 = shadedErrorBar(1:length(omegacorr_falloff_mean), omegacorr_falloff_mean, omegacorr_falloff_std, 'lineProps', lineProps) ;
xlabel(['time difference, \deltat [' QS.timeUnits, ']'])
ylabel('correlation')
title('omega')

dmyk = dmyk + 1 ;


%% legend for ax1
legend({'v3dcorr', 'v3ddotprod', '2hvn', 'divv', 'omega'})


%% Save figure
saveas(gcf, fullfile(QS.dir.uvCoord, 'autocorrelations_bwr.pdf')) 
for dd = 2:6
    set(gcf, 'CurrentAxes', ax{dd})
    xlim([-5, 139.5])
    
    ylim([-5, 139.5])
    
end
saveas(gcf, fullfile(QS.dir.uvCoord, 'autocorrelations_zoom_bwr.pdf')) 

for dd = 2:6
    set(gcf, 'CurrentAxes', ax{dd})
    colormap(brewermap(256, '*RdBu'))    
end
saveas(gcf, fullfile(QS.dir.uvCoord, 'autocorrelations_zoom_brewermap.pdf')) 


%% Correlation times
tau_v3d_mean = min(find(v3dcorr_falloff_mean < 1/exp(1))) ;
tau_v3d_upper = min(find(v3dcorr_falloff_mean + v3dcorr_falloff_std < 1/exp(1))) ;
tau_v3d_lower = min(find(v3dcorr_falloff_mean - v3dcorr_falloff_std < 1/exp(1))) ;
tau_v3ddprod_mean = min(find(v3ddprod_falloff_mean < 1/exp(1))) ;
tau_v3ddprod_upper = min(find(v3ddprod_falloff_mean + v3ddprod_falloff_std < 1/exp(1))) ;
tau_v3ddprod_lower = min(find(v3ddprod_falloff_mean - v3ddprod_falloff_std < 1/exp(1))) ;
tau_h2vn_mean = min(find(h2vncorr_falloff_mean < 1/exp(1))) ;
tau_h2vn_upper = min(find(h2vncorr_falloff_mean + h2vncorr_falloff_std < 1/exp(1))) ;
tau_h2vn_lower = min(find(h2vncorr_falloff_mean - h2vncorr_falloff_std < 1/exp(1))) ;
tau_divv_mean = min(find(divvcorr_falloff_mean < 1/exp(1))) ;
tau_divv_upper = min(find(divvcorr_falloff_mean + divvcorr_falloff_std < 1/exp(1))) ;
tau_divv_lower = min(find(divvcorr_falloff_mean - divvcorr_falloff_std < 1/exp(1))) ;
tau_omega_mean = min(find(omegacorr_falloff_mean < 1/exp(1))) ;
tau_omega_upper = min(find(omegacorr_falloff_mean + omegacorr_falloff_std < 1/exp(1))) ;
tau_omega_lower = min(find(omegacorr_falloff_mean - omegacorr_falloff_std < 1/exp(1))) ;

close all
figure;
errorbar([1,2,3,4,5],...
    [tau_v3d_mean, tau_v3ddprod_mean, tau_h2vn_mean, tau_divv_mean, tau_omega_mean], ...
    [tau_v3d_lower, tau_v3ddprod_lower, tau_h2vn_lower, tau_divv_lower, tau_omega_lower], ...
    [tau_v3d_upper, tau_v3ddprod_upper, tau_h2vn_upper, tau_divv_upper, tau_omega_upper], ...
    'marker', 'o', 'linestyle', 'none')
xlim([0, 6])
xticks([1,2,3,4,5])
xticklabels({'r_v', 'v_i \cdot v_j', 'r_{2Hv_n}', ...
    'r_{\nabla \cdot v}', 'r_{\nabla \times v}'})
ylabel('mean correlation time \tau')
xtickangle(gca, 45)
saveas(gcf, fullfile(QS.dir.uvCoord, 'autocorrelations_timescales.pdf'))
save(fullfile(QS.dir.uvCoord, 'autocorrelations_timescales.mat'), ...
    'tau_v3d_mean', 'tau_v3ddprod_mean', 'tau_h2vn_mean', 'tau_divv_mean', 'tau_omega_mean', ...
    'tau_v3d_lower', 'tau_v3ddprod_lower', 'tau_h2vn_lower', 'tau_divv_lower', 'tau_omega_lower', ...
    'tau_v3d_upper', 'tau_v3ddprod_upper', 'tau_h2vn_upper', 'tau_divv_upper', 'tau_omega_upper')



