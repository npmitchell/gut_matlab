function aux_plotPathlineMetricKinematicsLobes_integrated_subplots(QS, ...
    fn, fn_withH, lobes, tps, divv, H2vn, lobeYlabels, ...
    avgString, titleLobeBase, overwrite)
% Plot each lobe's integrated dynamics on a separate axis.
%
%
%
% NPMitchell 2020

if ~exist(fn, 'file') || ~exist(fn_withH, 'file') || overwrite
    close all
    set(gcf, 'visible', 'off')
    ymin = 0 ;
    ymax = 0 ;
    for jj = 1:length(lobes)
        axisColl{jj} = subplot(length(lobes), 1, jj) ;

        % div(v), H*2*vn, gdot means across the integrated ring section
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
            sgtitle([titleLobeBase, avgString], ...
                'Interpreter', 'Latex')
            legend({'$\Pi(1+\nabla\cdot\mathbf{v}_\parallel)$', ...
                '$\Pi(1+v_n 2H)$', ...
                '$\Pi(1+\frac{1}{2}\mathrm{Tr}\left[g^{-1} \dot{g} \right])$'}, ...
                'Interpreter', 'Latex', 'location', 'eastOutside')  
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

    % Add formatting to each
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