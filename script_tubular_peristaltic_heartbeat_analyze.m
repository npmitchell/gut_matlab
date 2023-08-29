% Look at peristaltic motion of zebrafish heart
lambda = 0 ;
lambda_mesh = 0 ;
lambda_err = 0; 
nmodes = 5 ;
zwidth = 2 ; 
period = 11 ;

sresStr = '' ;
mKDir = fullfile(tubi.dir.metricKinematics.root, ...
    strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d'], ...
    lambda, lambda_mesh, lambda_err, nmodes, zwidth), '.', 'p'));
outdir = fullfile(mKDir, 'measurements') ;

% Compute or load all timepoints
ntps = length(tubi.xp.fileMeta.timePoints(1:end-1)) ;
divv_apM = zeros(ntps, tubi.nU) ;
H2vn_apM = zeros(ntps, tubi.nU) ;
for tp = tubi.xp.fileMeta.timePoints(1:end-1)
    close all
    disp(['t = ' num2str(tp)])
    tidx = tubi.xp.tIdx(tp) ;

    % Check for timepoint measurement on disk
    Hfn = fullfile(outdir, sprintf('HH_vertices_%06d.mat', tp))   ;
    efn = fullfile(outdir, sprintf('gdot_vertices_%06d.mat', tp)) ;
    dfn = fullfile(outdir, sprintf('divv_vertices_%06d.mat', tp)) ;
    nfn = fullfile(outdir, sprintf('veln_vertices_%06d.mat', tp)) ;
    H2vnfn = fullfile(outdir, sprintf('H2vn_vertices_%06d.mat', tp)) ;
    rfn = fullfile(outdir, sprintf('radius_vertices_%06d.mat', tp)) ;

    % Load timeseries measurements
    load(dfn, 'divv_filt', 'divv_ap', 'divv_l', 'divv_r', 'divv_d', 'divv_v')
    load(H2vnfn, 'H2vn_filt', 'H2vn_ap', 'H2vn_l', 'H2vn_r', 'H2vn_d', 'H2vn_v') 
    
    %% Store in matrices
    % dv averaged
    divv_apM(tidx, :) = divv_ap ;
    H2vn_apM(tidx, :) = H2vn_ap ;
end

% Convert to period and cut off endcap values
divv_apM = divv_apM(:, 10:tubi.nV - 10) * period;
H2vn_apM = H2vn_apM(:, 10:tubi.nV - 10) * period ;
ddiv = diff(divv_apM, 1,  1) ;
dH2vn = diff(H2vn_apM, 1,  1) ;

divv_cut = divv_apM(1:end-1, :) ;
H2vn_cut = H2vn_apM(1:end-1, :) ;

%% Density plots
gridx1 = -0.8:.05:0.8;
gridx2 = -0.8:.05:0.8;
[x1g,x2g] = meshgrid(gridx1, gridx2);
x1 = x1g(:);
x2 = x2g(:);
xi = [x1 x2];
[densA, xi, h3] = ksdensity([divv_cut(:), dH2vn(:)], xi);
[densB, xi, h3] = ksdensity([H2vn_cut(:), ddiv(:)], xi);
xx = reshape(xi(:, 1), size(x1g)) ;
yy = reshape(xi(:, 2), size(x1g)) ;
densA = reshape(densA, size(x1g)) ;
densB = reshape(densB, size(x1g)) ;

% colormap 
gamma = 1.5 ;
tmp = cubehelix(256,0.,1,1,gamma,[0,1],[0,1]) ;

h1 = subplot(1, 2, 1)
imagesc(gridx1, gridx2, densA)
axis equal ;
axis tight
xlabel('$\nabla \cdot \textbf{v}_{\parallel}$', 'interpreter', 'latex')
ylabel('$\partial_t \left[ 2Hv_n \right]$', 'interpreter', 'latex')
% caxis([0, 2.5])
h2 = subplot(1, 2, 2)
imagesc(gridx1, gridx2, densB)
colormap(flipud(tmp))
axis equal
axis tight
xlabel('$2Hv_n$', 'interpreter', 'latex')
ylabel('$\partial_t \left[ \nabla \cdot \textbf{v}_{\parallel}\right]$', 'interpreter', 'latex')
% caxis([0, 2.5])

[fitA, covA ] = fit(divv_cut(:), dH2vn(:), 'poly1') ;
[fitB, covB ] = fit(H2vn_cut(:), ddiv(:), 'poly1') ;

confsA = confint(fitA) ;
confsB = confint(fitB) ;
uncA = (fitA.p1 - confsA(1, 1)) / 1.96 ;
uncB = (fitB.p1 - confsB(1, 1)) / 1.96  ;

subplot(1, 2, 1) 
set(h1, 'Ydir', 'normal')
title(['$m=$', sprintf('%0.2f', fitA.p1), '$\pm$' sprintf('%0.2f', uncA)], ...
    'interpreter', 'latex')

subplot(1, 2, 2) 
set(h2, 'Ydir', 'normal')
title(['$m=$', sprintf('%0.2f', fitB.p1), '$\pm$' sprintf('%0.2f', uncB)], ...
    'interpreter', 'latex')

saveas(gcf, fullfile(mKDir, 'metricKinematics_local_peristaltic.png'))
saveas(gcf, fullfile(mKDir, 'metricKinematics_local_peristaltic.pdf'))
