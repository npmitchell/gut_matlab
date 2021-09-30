% aux_plotCellSegmentation3DStats

%% Define some colors
colors = define_colors() ;
bluecol = colors(1, :) ;
redcol = colors(2, :) ;
yelcol = colors(3, :) ;
CScolors = colors([9,5], :) ;
figW = 18 ;
figH = 9 ;

%% Plot mean +/- pctile over time
imfn = fullfile(QS.dir.segmentation, segSubDir, 'cell_anisotropy') ;
timestamps = timePoints - t0 ;
if contains(lower(QS.timeUnits), 'min')
    timestamps = timestamps / 60 ;
    timeunits = 'hr' ;
else
    timeunits = QS.timeUnits ;
end
xlims = [min(timestamps) - 0.25*min(diff(timestamps)), ...
        max(timestamps) + 0.25*min(diff(timestamps))] ;

close all 
fig = figure('units', 'centimeters', 'position', [0, 0, figW, figH]) ;
stdStyles = {'_quartiles', '_std'} ;
for stdStyle = 1:2
    % shade(timePoints - t0, bndlow, timePoints, bndhigh)
    x2 = [timestamps, fliplr(timestamps)] ;
    subplot(1, 2, 1)
    if stdStyle == 1
        fill(x2, [mratio_low25, fliplr(mratio_high75)], bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
    else
        fill(x2, [mean_mratio - mratio_std, ...
            fliplr(mean_mratio + mratio_std)], ...
            bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
    end
    hold on;
    errorbar(timestamps, mean_mratio, mratio_ste, ...
        '-', 'color', bluecol)
    xlim(xlims) 
    xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
    ylabel('cell aspect ratio, $\sqrt{I_{1}/I_{2}}$',   'interpreter', 'latex')
    title('endoderm anisotropy over time', 'interpreter', 'latex')
    
    % Cos/Sin(2theta)
    set(fig,'defaultAxesColorOrder', CScolors);
    subplot(1, 2, 2)
    yyaxis left
    if stdStyle == 1
        fill(x2, [c2t_low25, fliplr(c2t_high75)], CScolors(1, :), 'facealpha', 0.3, 'edgecolor', 'none');
    else
        fill(x2, [min(1, max(-1, mc2t - c2t_std)), ...
            fliplr(min(1, max(-1,mc2t + c2t_std)))],...
            CScolors(1, :), 'facealpha', 0.3, 'edgecolor', 'none');
    end
    hold on;
    % shadedErrorBar(timePoints - t0, mean(y,1),std(y),'lineProps','g');
    errorbar(timestamps, mc2t, ...
        min(c2t_ste, mc2t+1), ...  % bound from below at -1      
        min(c2t_ste, 1-mc2t), ... % bound from above at 1
        '-', 'color', CScolors(1, :))
    ylim([-1,1])
    yyaxis right ;
    if stdStyle == 1
        fill(x2, [s2t_low25, fliplr(s2t_high75)], CScolors(2, :), 'facealpha', 0.3, 'edgecolor', 'none');
    else
        fill(x2, [min(1, max(-1, ms2t - s2t_std)), ...
            fliplr(min(1, max(-1,ms2t + s2t_std)))], ...
            CScolors(2, :), 'facealpha', 0.3, 'edgecolor', 'none');
    end
    errorbar(timestamps, ms2t, ...
        min(s2t_ste, ms2t+1), ... % bound from below at -1
        min(s2t_ste, 1-ms2t), ... % bound from above at 1
        '-', 'color', CScolors(2, :))
    ylim([-1,1])

    yyaxis left
    ylabel('cell orientation, $\cos 2\theta$',   'interpreter', 'latex')
    yyaxis right
    ylabel('cell orientation, $\sin 2\theta$',   'interpreter', 'latex')
    
    xlim(xlims)
    xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
    title('endoderm orientation over time', 'interpreter', 'latex')

    saveas(gcf, [imfn, stdStyles{stdStyle} '.png'])
    saveas(gcf, [imfn, stdStyles{stdStyle} '.pdf'])
    clf
end

%% Plot nematic strength and direction for each lobe 
imfn = fullfile(QS.dir.segmentation, segSubDir, ...
    'cell_anisotropy_lobes') ;
close all
fig = figure('units','centimeters','position',[0,0,figH, figH], 'visible', 'off') ;
for lobe = 1:nLobes
    subplot(ceil(nLobes * 0.5), 2, lobe)
    midline = squeeze(meanQLobeAspects(lobe, :)) ;
    uncs = squeeze(meanQLobeAspectStds(lobe, :)) ;
    timestamps = timePoints - t0 ;
    if contains(QS.timeUnits, 'min')
        timestamps = timestamps / 60 ;
    end
    x2 = [timestamps, fliplr(timestamps)] ;
    fill(x2,[midline-uncs, fliplr(midline+uncs)], ...
        bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
    hold on;
    plot(timestamps, midline, '.-', 'color', bluecol)
    ylim([1, Inf])
    
    yyaxis right
    % fill(x2, [c2t_low, fliplr(c2t_high)], redcol, 'facealpha', 0.3, 'edgecolor', 'none');
    hold on;
    plot(timestamps, mod(meanQLobeThetas(lobe, :), pi)/pi, '.-')
    % plot(timestamps, sin(2*meanQLobeThetas(lobe, :)), '.-')
    % 'color', redcol)
    ylim([0, 1])
    
    if mod(lobe, 2) == 1 
        yyaxis left
        ylabel(['lobe ' num2str(lobe) ' aspect ratio'],   'interpreter', 'latex')
    else
        yyaxis right
        ylabel('nematic orientation $\theta/\pi$',   'interpreter', 'latex')
    end
    
    % Time label
    if lobe > nLobes - 2
        if contains(QS.timeUnits, 'min')
            xlabel('time [hr]', 'interpreter', 'latex')  
        else
            xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')  
        end
    end
end
sgtitle('endoderm orientation over time', 'interpreter', 'latex')
saveas(gcf, [imfn '.png'])
saveas(gcf, [imfn '.pdf'])

%% Plot nematic strength and direction for each lobe as shaded errorplot
imfn = fullfile(QS.dir.segmentation, segSubDir, ...
    'cell_anisotropy_lobes_signed') ;
close all
fig = figure('units','centimeters','position',[0,0,figH,figH]) ;
for lobe = fliplr(1:nLobes)
    % subplot(ceil(nLobes * 0.5), 2, lobe)
    c2t = cos(2*meanQLobeThetas(lobe, :))  ;
    midline = 0.5 *c2t .* (squeeze(meanQLobeAspects(lobe, :)) - 1) ;
    uncs = 0.5 *c2t .* (squeeze(meanQLobeAspectStds(lobe, :)) - 1) ;
    timestamps = timePoints - t0 ;
    if contains(QS.timeUnits, 'min')
        timestamps = timestamps / 60 ;
    end
    x2 = [timestamps, fliplr(timestamps)] ;
    fill(x2,[midline-abs(uncs), fliplr(midline+abs(uncs))], ...
        colors(lobe, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
        'HandleVisibility', 'off');
    hold on;
    hs{lobe} = plot(timestamps, midline, '.-', 'color', colors(lobe, :)) ;
    
    legendentries{lobe} = ['chamber ' num2str(lobe)] ;
end
% ylims = ylim() ;
% ylim([-max(abs(ylims)), max(abs(ylims))])

% Mark zero line
plot(timestamps, 0*timestamps, 'k--', 'HandleVisibility','off')
% Labels
legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
ylabel('elongation, $Q_{xx}= \left(\sqrt{I_2/I_1}-1\right) \cos 2\theta$', ...
    'interpreter', 'latex')
if contains(QS.timeUnits, 'min')
    xlabel('time [hr]', 'interpreter', 'latex')  
else
    xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')  
end
sgtitle('endoderm orientation over time', 'interpreter', 'latex')

saveas(gcf, [imfn, '.png'])
saveas(gcf, [imfn, '.pdf'])


%% Plot histograms
imfn = fullfile(QS.dir.segmentation, segSubDir, 'cell_anisotropy_hist') ;

close all
fig = figure('units','centimeters','position',[0,0,7, 7]) ;
colormap(inferno)
subplot(2, 2, 1)
imagesc(timePoints - t0, edges, cos2thetaM)
set(gca,'YDir','normal')
caxis([0, 0.05])
xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
ylabel('cell orientation, $\cos2\theta$',   'interpreter', 'latex')
subplot(2, 2, 2)
imagesc(timePoints - t0, edges, sin2thetaM)
set(gca,'YDir','normal')
caxis([0, 0.05])
xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
ylabel('cell orientation, $\sin2\theta$',   'interpreter', 'latex')
subplot(2, 2, 3)
imagesc(timePoints - t0, edgesAR, aspectM)
set(gca,'YDir','normal')
caxis([0, 0.05])
xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
ylabel('cell aspect ratio',   'interpreter', 'latex')
subplot(4, 2, 6)
cb = colorbar('location', 'south') ;
ylabel(cb, 'probability', 'interpreter', 'latex')
caxis([0, 0.05])
axis off
saveas(gcf, [imfn '.png'])
saveas(gcf, [imfn '.pdf']) 

%% Aspect ratio distributions, variation of mean shape aspect
imfn = fullfile(QS.dir.segmentation, segSubDir, 'cell_anisotropy_mratio') ;
clf
colors = define_colors() ;
bluecol = colors(1, :) ;
% shade(timePoints - t0, bndlow, timePoints, bndhigh)
x2 = [timePoints - t0, fliplr(timePoints - t0)] ;
fill(x2, [mratio_low25, fliplr(mratio_high75)], ...
    bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
hold on;
% shadedErrorBar(timePoints - t0, mean(y,1),std(y),'lineProps','g');
plot(timePoints - t0, median_moiratio, '.-')
plot(timePoints - t0, mean_moiratio, '.-')
plot(timePoints - t0, mean_mratio, '.-')
legend({'25-75\%', 'median $\sqrt{I_1/I_2}$', ...
    '$\sqrt{\lambda_1^{\langle I \rangle}/\lambda_2^{\langle I \rangle}}$', ...
    '$2||\langle Q\rangle|| + 1$'}, 'interpreter', 'latex')

xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
ylabel('aspect ratio',   'interpreter', 'latex')
title('endoderm orientation over time', 'interpreter', 'latex')

saveas(gcf, [imfn, '.png'])
saveas(gcf, [imfn, '.pdf'])



%% Plot as a function of AP position and time (kymograph)
imfn = fullfile(QS.dir.segmentation, segSubDir, 'ap_kymographs_qc2t_qs2t') ;
clf
if ~exist(imfn, 'file') || overwriteImages
    subplot(1, 2, 1)
    % imagesc(mid_ap, timePoints-t0, medfilt2(mean_qc2ts, [3, 1])) ;
    imagesc(mid_ap, timestamps, mean_qc2ts) ;
    caxis([-2.5, 2.5])
    colormap(blueblackred)
    % gray out NaNs
    whiteout = double(isnan(mean_qc2ts)) ;
    hold on; 
    wh = imagesc(mid_ap,timestamps, cat(3, whiteout, whiteout, whiteout)) ;
    set(wh, 'alphaData', whiteout)
    xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
    ylabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
    cb = colorbar('location', 'southOutside') ;
    ylabel(cb, '$|Q|\cos 2\theta$', ...
        'interpreter', 'latex')
    subplot(1, 2, 2)
    % imagesc(mid_ap, timestamps, medfilt2(mean_qs2ts, [3, 1])) ;
    imagesc(mid_ap, timestamps, mean_qs2ts) ;
    caxis([-2.5, 2.5])
    colormap(blueblackred)
    hold on; 
    wh = imagesc(mid_ap,timestamps, cat(3, whiteout, whiteout, whiteout)) ;
    set(wh, 'alphaData', whiteout)
    xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
    ylabel(['time [' timeunits ']'], 'interpreter', 'latex')
    cb = colorbar('location', 'southOutside') ;
    ylabel(cb, '$|Q|\sin 2\theta$', ...
        'interpreter', 'latex')
    sgtitle('cell anisotropy kymographs', 'interpreter', 'latex')
    saveas(gcf, [imfn, '.png'])
    saveas(gcf, [imfn, '.pdf'])
end

%% Time derivative of filtered image AP position
imfn = fullfile(QS.dir.segmentation, segSubDir, 'ap_kymographs_dc2t_ds2t') ;
clf
clim = 0.15 ;
if ~exist(imfn, 'file') || overwriteImages
    % cfiltered = medfilt2(mean_qc2ts, [3, 1]) ;
    [~, dc2t] = gradient(mean_qc2ts, diff(timeunits)) ;
    % sfiltered = medfilt2(mean_qs2ts, [3, 1]) ;
    [~, ds2t] = gradient(mean_qs2ts, diff(timeunits)) ;
    subplot(1, 2, 1)
    imagesc(mid_ap, timestamps, dc2t) ;
    % gray out NaNs
    whiteout = double(isnan(ds2t)) ;
    hold on; 
    wh = imagesc(mid_ap,timestamps, cat(3, whiteout, whiteout, whiteout)) ;
    set(wh, 'alphaData', whiteout)
    caxis([-clim, clim])
    xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
    ylabel(['time [' timeunits ']'], 'interpreter', 'latex')
    colormap(blueblackred)
    cb = colorbar('location', 'southOutside') ;
    ylabel(cb, ...
        ['$\partial_t\left(|Q| \cos 2\theta\right)$ [' timeunits '$^{-1}$]'], ...
        'interpreter', 'latex')
    subplot(1, 2, 2)
    imagesc(mid_ap, timestamps, ds2t) ;
    % gray out NaNs
    hold on; 
    wh = imagesc(mid_ap,timestamps, cat(3, whiteout, whiteout, whiteout)) ;
    set(wh, 'alphaData', whiteout)
    caxis([-clim clim])
    xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
    ylabel(['time [' timeunits ']'], 'interpreter', 'latex')
    cb = colorbar('location', 'southOutside') ;
    ylabel(cb,  ...
        ['$\partial_t\left(|Q| \sin 2\theta\right)$ [' timeunits '$^{-1}$]'], ...
        'interpreter', 'latex')
    sgtitle('cell anisotropy kymographs', 'interpreter', 'latex')
    saveas(gcf, [imfn '.png'])
    saveas(gcf, [imfn '.pdf'])
end






% OLD CODE FROM generateCellSegmentationPathlines3D.m for plotting
% imfn = fullfile(QS.dir.segmentation, 'pathlines', [segSubDir 'cell_anisotropy.png']) ;
% timestamps = timePoints - t0 ;
% if contains(lower(QS.timeUnits), 'min')
%     timestamps = timestamps / 60 ;
%     timeunits = 'hr' ;
% else
%     timeunits = QS.timeUnits ;
% end
% xlims = [min(timestamps) - 0.25*min(diff(timestamps)), ...
%         max(timestamps) + 0.25*min(diff(timestamps))] ;
% 
% close all 
% fig = figure('units', 'centimeters', 'position', [0, 0, figW, figH]) ;
% stdStyles = {'_quartiles', '_std'} ;
% for stdStyle = 1:2
%     % shade(timePoints - t0, bndlow, timePoints, bndhigh)
%     x2 = [timestamps, fliplr(timestamps)] ;
%     subplot(1, 2, 1)
%     if stdStyle == 1
%         fill(x2, [mratio_low25, fliplr(mratio_high75)], bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
%     else
%         fill(x2, [mean_mratio - mratio_std, ...
%             fliplr(mean_mratio + mratio_std)], ...
%             bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
%     end
%     hold on;
%     errorbar(timestamps, mean_mratio, mratio_ste, ...
%         '-', 'color', bluecol)
%     xlim(xlims) 
%     xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
%     ylabel('cell aspect ratio, $\sqrt{I_{1}/I_{2}}$',   'interpreter', 'latex')
%     title('advected endoderm anisotropy over time', 'interpreter', 'latex')
%     
%     % Cos/Sin(2theta)
%     set(fig,'defaultAxesColorOrder', CScolors);
%     subplot(1, 2, 2)
%     yyaxis left
%     if stdStyle == 1
%         fill(x2, [c2t_low25, fliplr(c2t_high75)], CScolors(1, :), 'facealpha', 0.3, 'edgecolor', 'none');
%     else
%         fill(x2, [min(1, max(-1, mc2t - c2t_std)), ...
%             fliplr(min(1, max(-1,mc2t + c2t_std)))],...
%             CScolors(1, :), 'facealpha', 0.3, 'edgecolor', 'none');
%     end
%     hold on;
%     % shadedErrorBar(timePoints - t0, mean(y,1),std(y),'lineProps','g');
%     errorbar(timestamps, mc2t, ...
%         min(c2t_ste, mc2t+1), ...  % bound from below at -1      
%         min(c2t_ste, 1-mc2t), ... % bound from above at 1
%         '-', 'color', CScolors(1, :))
%     ylim([-1,1])
%     yyaxis right ;
%     if stdStyle == 1
%         fill(x2, [s2t_low25, fliplr(s2t_high75)], CScolors(2, :), 'facealpha', 0.3, 'edgecolor', 'none');
%     else
%         fill(x2, [min(1, max(-1, ms2t - s2t_std)), ...
%             fliplr(min(1, max(-1,ms2t + s2t_std)))], ...
%             CScolors(2, :), 'facealpha', 0.3, 'edgecolor', 'none');
%     end
%     errorbar(timestamps, ms2t, ...
%         min(s2t_ste, ms2t+1), ... % bound from below at -1
%         min(s2t_ste, 1-ms2t), ... % bound from above at 1
%         '-', 'color', CScolors(2, :))
%     ylim([-1,1])
% 
%     yyaxis left
%     ylabel('cell orientation, $\cos 2\theta$',   'interpreter', 'latex')
%     yyaxis right
%     ylabel('cell orientation, $\sin 2\theta$',   'interpreter', 'latex')
%     
%     xlim(xlims)
%     xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
%     title('advected endoderm orientation over time', 'interpreter', 'latex')
% 
%     saveas(gcf, [imfn, stdStyles{stdStyle} '.png'])
%     saveas(gcf, [imfn, stdStyles{stdStyle} '.pdf'])
%     clf
% end
% 
% 
% %% Plot nematic strength and direction for each lobe 
% imfn = fullfile(QS.dir.segmentation, 'pathlines', 'cell_anisotropy_lobes') ;
% close all
% fig = figure('units','centimeters','position',[0,0,figH, figH], 'visible', 'off') ;
% 
% for lobe = 1:nLobes
%     subplot(ceil(nLobes * 0.5), 2, lobe)
%     midline = squeeze(meanQLobeAspects(lobe, :)) ;
%     uncs = squeeze(meanQLobeAspectStds(lobe, :)) ;
%     timestamps = timePoints - t0 ;
%     if contains(QS.timeUnits, 'min')
%         timestamps = timestamps / 60 ;
%     end
%     x2 = [timestamps, fliplr(timestamps)] ;
%     fill(x2,[midline-uncs, fliplr(midline+uncs)], ...
%         bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
%     hold on;
%     plot(timestamps, midline, '.-', 'color', bluecol)
%     ylim([1, Inf])
%     
%     yyaxis right
%     % fill(x2, [c2t_low, fliplr(c2t_high)], redcol, 'facealpha', 0.3, 'edgecolor', 'none');
%     hold on;
%     plot(timestamps, mod(meanQLobeThetas(lobe, :), pi)/pi, '.-')
%     % plot(timestamps, sin(2*meanQLobeThetas(lobe, :)), '.-')
%     % 'color', redcol)
%     ylim([0, 1])
%     
%     if mod(lobe, 2) == 1 
%         yyaxis left
%         ylabel(['lobe ' num2str(lobe) ' aspect ratio'],   'interpreter', 'latex')
%     else
%         yyaxis right
%         ylabel('nematic orientation $\theta/\pi$',   'interpreter', 'latex')
%     end
%     
%     % Time label
%     if lobe > nLobes - 2
%         if contains(QS.timeUnits, 'min')
%             xlabel('time [hr]', 'interpreter', 'latex')  
%         else
%             xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')  
%         end
%     end
% end
% sgtitle('advected endoderm orientation over time', 'interpreter', 'latex')
% saveas(gcf, [imfn '.png'])
% saveas(gcf, [imfn '.pdf'])
% 
% %% Plot nematic strength and direction for each lobe 
% imfn = fullfile(QS.dir.segmentation, 'pathlines', ...
%     'cell_anisotropy_lobes_signed') ;
% close all
% fig = figure('units','centimeters','position',[0,0,figH, figH]) ;
% 
% for lobe = 1:nLobes
%     % subplot(ceil(nLobes * 0.5), 2, lobe)
%     c2t = cos(2*meanQLobeThetas(lobe, :))  ;
%     midline = 0.5 *c2t .* (squeeze(meanQLobeAspects(lobe, :)) - 1) ;
%     uncs = 0.5 *c2t .* (squeeze(meanQLobeAspectStds(lobe, :)) - 1) ;
%     timestamps = timePoints - t0 ;
%     if contains(QS.timeUnits, 'min')
%         timestamps = timestamps / 60 ;
%     end
%     x2 = [timestamps, fliplr(timestamps)] ;
%     fill(x2,[midline-abs(uncs), fliplr(midline+abs(uncs))], ...
%         colors(lobe, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
%         'HandleVisibility', 'off');
%     hold on;
%     hs{lobe} = plot(timestamps, midline, '.-', 'color', colors(lobe, :)) ;
%     
%     legendentries{lobe} = ['chamber ' num2str(lobe)] ;
% end
% % ylims = ylim() ;
% % ylim([-max(abs(ylims)), max(abs(ylims))])
% 
% % Mark zero line
% plot(timestamps, 0*timestamps, 'k--', 'HandleVisibility','off')
% % Labels
% legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
% ylabel('elongation, $2||Q|| \cos 2\theta$',   'interpreter', 'latex')
% if contains(QS.timeUnits, 'min')
%     xlabel('time [hr]', 'interpreter', 'latex')  
% else
%     xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')  
% end
% sgtitle('endoderm orientation over time', 'interpreter', 'latex')
% 
% saveas(gcf, [imfn, '.png'])
% saveas(gcf, [imfn, '.pdf'])
% 
% %% Plot histograms
% imfn = fullfile(QS.dir.segmentation, 'pathlines', 'cell_anisotropy_hist.png') ;
% clf
% colormap(cividis)
% subplot(2, 2, 1)
% imagesc(timePoints - t0, edges, cos2thetaM)
% set(gca,'YDir','normal')
% caxis([0, 0.05])
% xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
% ylabel('nematic orientation $\cos2\theta$',   'interpreter', 'latex')
% subplot(2, 2, 2)
% imagesc(timePoints - t0, edges, sin2thetaM)
% set(gca,'YDir','normal')
% caxis([0, 0.05])
% xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
% ylabel('nematic orientation $\sin2\theta$',   'interpreter', 'latex')
% subplot(2, 2, 3)
% imagesc(timePoints - t0, edgesAR, aspectM)
% set(gca,'YDir','normal')
% caxis([0, 0.05])
% xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
% ylabel('aspect ratio $\sqrt{I_{1}/I_{2}}$',   'interpreter', 'latex')
% subplot(4, 2, 6)
% cb = colorbar('location', 'south') ;
% ylabel(cb, 'probability', 'interpreter', 'latex')
% caxis([0, 0.05])
% axis off
% saveas(gcf, imfn) 
% 
% %% Aspect ratio distributions, variation of mean shape aspect
% imfn = fullfile(QS.dir.segmentation, 'pathlines', 'cell_anisotropy_mratio.png') ;
% clf
% colors = define_colors() ;
% bluecol = colors(1, :) ;
% % shade(timePoints - t0, bndlow, timePoints, bndhigh)
% x2 = [timePoints - t0, fliplr(timePoints - t0)] ;
% fill(x2, [mratio_low, fliplr(mratio_high)], ...
%     bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
% hold on;
% % shadedErrorBar(timePoints - t0, mean(y,1),std(y),'lineProps','g');
% plot(timePoints - t0, median_moiratio, '.-')
% plot(timePoints - t0, mean_moiratio, '.-')
% plot(timePoints - t0, mean_mratio, '.-')
% legend({'25-75\%', 'median $\sqrt{I_1/I_2}$', ...
%     '$\sqrt{\lambda_1^{\langle I \rangle}/\lambda_2^{\langle I \rangle}}$', ...
%     '$2||\langle Q\rangle|| + 1$'}, 'interpreter', 'latex')
% 
% xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
% ylabel('aspect ratio',   'interpreter', 'latex')
% title('endoderm orientation over time', 'interpreter', 'latex')
% saveas(gcf, imfn)
% 
% 
% 
% %% Plot as a function of AP position and time (kymograph)
% imfn = fullfile(QS.dir.segmentation, 'pathlines', 'ap_kymographs_c2t_s2t.png') ;
% clf
% subplot(1, 2, 1)
% imagesc(mid_ap, timePoints-t0, medfilt2(mean_c2ts, [3, 1])) ;
% caxis([-10, 10])
% colormap(blueblackred)
% xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
% ylabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
% cb = colorbar('location', 'southOutside') ;
% ylabel(cb, '$\sqrt{\lambda_1/\lambda_2}\cos 2\theta$', ...
%     'interpreter', 'latex')
% subplot(1, 2, 2)
% imagesc(mid_ap, timePoints-t0, medfilt2(mean_s2ts, [3, 1])) ;
% caxis([-10, 10])
% colormap(blueblackred)
% xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
% ylabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
% cb = colorbar('location', 'southOutside') ;
% ylabel(cb, '$\sqrt{\lambda_1/\lambda_2}\sin 2\theta$', ...
%     'interpreter', 'latex')
% sgtitle('cell anisotropy kymographs', 'interpreter', 'latex')
% saveas(gcf, imfn)
% 
% %% Time derivative of filtered image AP position
% imfn = fullfile(QS.dir.segmentation, 'pathlines', 'ap_kymographs_dc2t_ds2t.png') ;
% clf
% cfiltered = medfilt2(mean_c2ts, [3, 1]) ;
% [~, dc2t] = gradient(cfiltered) ;
% sfiltered = medfilt2(mean_s2ts, [3, 1]) ;
% [~, ds2t] = gradient(sfiltered) ;
% subplot(1, 2, 1)
% imagesc(mid_ap, timePoints-t0, imgaussfilt(medfilt2(dc2t, [3, 1]), 0.5)) ;
% caxis([-3, 3])
% xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
% ylabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
% colormap(blueblackred)
% cb = colorbar('location', 'southOutside') ;
% ylabel(cb, '$\partial_t\left(\sqrt{\lambda_1/\lambda_2} \cos 2\theta\right)$', ...
%     'interpreter', 'latex')
% subplot(1, 2, 2)
% imagesc(mid_ap, timePoints-t0, imgaussfilt(medfilt2(ds2t, [3, 1]), 0.5)) ;
% caxis([-3, 3])
% xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
% ylabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
% cb = colorbar('location', 'southOutside') ;
% ylabel(cb, '$\partial_t\left(\sqrt{\lambda_1/\lambda_2} \sin 2\theta\right)$', ...
%     'interpreter', 'latex')
% sgtitle('cell anisotropy kymographs', 'interpreter', 'latex')
% saveas(gcf, imfn)