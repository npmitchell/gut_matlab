function aux_plotPathlineMetricKinematicsLobes_subpanels(QS, m2plot, ...
    fn, fn_withH, lobes, tps, divv, H2vn, HH, lobeYlabels, avgString,...
    titleLobeBase, divvcolor, H2vncolor, ...
    Hposcolor, Hnegcolor, Hsz, overwrite)
% Plot instantaneous kinematics
%
%
% NPMitchell 2020

close all
set(gcf, 'visible', 'off')
ymin = 0 ;
ymax = 0 ;
for jj = 1:length(lobes)
    axisColl{jj} = subplot(length(lobes), 1, jj) ;
    dvj = mean(divv(:, lobes{jj}), 2) ;
    Hvj = mean(H2vn(:, lobes{jj}), 2) ;
    plot(tps, dvj, '.-', 'Color', divvcolor)
    hold on;
    plot(tps, Hvj, '.-', 'Color', H2vncolor)
    if jj == 1
        % Title and labels
        sgtitle([titleLobeBase, avgString], ...
            'Interpreter', 'Latex')
        legend({'$\nabla\cdot\mathbf{v}_\parallel$', ...
            '$v_n 2H$'}, 'Interpreter', 'Latex', ...
            'location', 'eastOutside')  
        drawnow
        pos = get(gca, 'position') ;
    elseif jj == length(lobes)
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
        'markerFaceColor', divvcolor, ...
        'HandleVisibility', 'off')
    hold on;
    scatter(tps(Hvpos), ...
        ymin * ones(size(tps(Hvpos))), 5,  ...
        'markeredgecolor', 'none', 'markerFaceAlpha', 0.6, ...
        'markerFaceColor', H2vncolor, ...
        'HandleVisibility', 'off')
    ylim([ymin-0.1*(ymax-ymin), ymax])
end

% Save figure
disp(['Saving figure: ', fn])
saveas(gcf, fn)

%% Add mean curvature to each plot -- lobes, instantaneous
for jj = 1:length(lobes)
    axes(axisColl{jj})

    % Mark the instantaneous mean curvature
    Hj = mean(HH(:, lobes{jj}), 2) ;
    scatter(tps(Hj>0), Hj(Hj>0), Hsz, 's', 'filled', ...
        'markeredgecolor', Hposcolor, ...
        'markerfacecolor', Hposcolor)
    scatter(tps(Hj<0), Hj(Hj<0), Hsz, 's', ...
        'markeredgecolor', Hnegcolor)
    ymin2 = min(min(Hj)-0.1*(ymax-ymin), ymin-0.1*(ymax-ymin)) ;
    ymax2 = max(ymax, max(Hj)+0.1*(ymax-ymin)) ;
    ylim([ymin2, ymax2])
end
axes(axisColl{1})
legend({'$\nabla\cdot\mathbf{v}_\parallel$', ...
    '$v_n 2H$', '$H > 0$', '$H<0$'}, 'Interpreter', 'Latex', ...
    'location', 'eastOutside')  

% Save figure with mean curvature
disp(['Saving figure: ', fn_withH])
saveas(gcf, fn_withH)

