function aux_plotPathlineMetricKinematicsFolds_subpanels(QS, fn, fn_withH, ...
            foldIds, width, nU, tps, divv, H2vn, HH, titleFoldBase, ...
            foldYlabels, avgString, divvcolor, H2vncolor, ...
            Hposcolor, Hnegcolor, Hsz, overwrite)
%
%
%
%
% NPMitchell

if ~exist(fn, 'file') || ~exist(fn_withH, 'file') || overwrite 
    close all
    set(gcf, 'visible', 'off')
    ymin = 0 ;
    ymax = 0 ;
    % Each fold is valley+/- width
    for jj = 1:length(foldIds)
        axisColl{jj} = subplot(length(foldIds), 1, jj) ;
        valley = (foldIds(jj)-width):(foldIds(jj)+width) ;
        dvj = mean(divv(:, valley), 2) ;
        Hvj = mean(H2vn(:, valley), 2) ;
        plot(tps, dvj, '.-', 'Color', divvcolor)
        hold on;
        plot(tps, Hvj, '.-', 'Color', H2vncolor)
        if jj == 1
            % Title and labels
            sgtitle([titleFoldBase, avgString, ', ', ...
                '$w_{\textrm{fold}}=', ...
                num2str(100*(2*width + 1)/ nU), '$\%$\, L_\zeta$'], ...
                'Interpreter', 'Latex')
            legend({'$\nabla\cdot\mathbf{v}_\parallel$', ...
                '$v_n 2H$'}, 'Interpreter', 'Latex', ...
                'location', 'eastOutside')  
            drawnow
            pos = get(gca, 'position') ;
        elseif jj == length(foldIds)
            xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
        end
        ylabel(foldYlabels{jj}, 'Interpreter', 'Latex')
        ylims = ylim() ;
        ymin = min(ymin, ylims(1)) ;
        ymax = max(ymax, ylims(2)) ;
    end

    for jj = 1:length(foldIds)
        axes(axisColl{jj})
        pos2 = get(gca, 'position') ;
        set(gca, 'position', [pos2(1) pos2(2) pos(3) pos2(4)])
        ylim([ymin, ymax])

        % Mark wherever the divv<0 and also vn2H<0
        valley = (foldIds(jj)-width):(foldIds(jj)+width) ;
        dvj = mean(divv(:, valley), 2) ;
        Hvj = mean(H2vn(:, valley), 2) ;
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

    %% Add mean curvature to each plot
    for jj = 1:length(foldIds)
        valley = (foldIds(jj)-width):(foldIds(jj)+width) ;
        axes(axisColl{jj})
        % yyaxis right

        % Mark the instantaneous mean curvature
        Hj = mean(HH(:, valley), 2) ;
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
end