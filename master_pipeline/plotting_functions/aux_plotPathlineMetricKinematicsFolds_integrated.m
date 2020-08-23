function aux_plotMetricKinematicsFolds_integrated(...
    fn, fn_withH, foldIds, width, divv, H2vn, HH, foldYlabels, ...
    avgString, m2plot, overwrite)
%aux_plotMetricKinematicsFolds_2panel_withH
%(fn, divv, H2vn, HH, width, foldYlabels)
%
% Parameters
% ----------
% fn : str
%   path to output figure filename
% divv : 
%
% NPMitchell 2020


 if ~exist(fn, 'file') || overwrite 
    close all
    % Each fold is valley+/- width
    for jj = 1:length(foldIds)
        valley = (foldIds(jj)-width):(foldIds(jj)+width) ;

        % div(v), H*2*vn, gdot
        switch lower(m2plot)
            case 'gdot'
                gdj = mean(divv(:, valley) - H2vn(:, valley), 2) ;
            case 'divv'
                gdj = mean(divv(:, valley), 2) ;
            case 'h2vn'
                gdj = mean(H2vn(:, valley), 2) ;
        end

        % Take cumulative product marching forward from t0
        gpj_pos = cumprod(1 + gdj(tps > eps)) ;
        gpj_neg = flipud(cumprod(flipud(1 ./ (1 + gdj(tps < eps))))) ;
        gpj = cat(1, gpj_neg, gpj_pos) ;               

        % Plot this fold
        plot(tps, gpj, '.-', 'Color', QS.plotting.colors(jj, :))
        hold on;
    end    

    % Title and labels
    sgtitle(['Tissue dilation in folds, ', avgString, ', ', ...
        '$w_{\textrm{fold}}=', ...
        num2str(100*(2*width + 1)/ nU), '$\%$\, L_\zeta$'], ...
        'Interpreter', 'Latex')
    legend(foldYlabels, 'Interpreter', 'Latex', 'location', 'eastOutside')  
    drawnow
    xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
    ylabel('$\Pi\big(1+$d$t\, \frac{1}{2}\mathrm{Tr} \left[g^{-1} \dot{g}\right]\big)$', ...
        'Interpreter', 'Latex')

    % Save figure
    disp(['Saving figure: ', fn])
    saveas(gcf, fn)
end


%% Same plot but add mean curvature as second panel
if ~exist(fn_withH, 'file') || overwrite 
    close all
    subplot(2, 1, 1)
    % Each fold is valley+/- width
    for jj = 1:length(foldIds)
        valley = (foldIds(jj)-width):(foldIds(jj)+width) ;

        % div(v), H*2*vn, gdot
        switch lower(m2plot)
            case 'gdot'
                gdj = mean(divv(:, valley) - H2vn(:, valley), 2) ;
            case 'divv'
                gdj = mean(divv(:, valley), 2) ;
            case 'h2vn'
                gdj = mean(H2vn(:, valley), 2) ;
        end
        % Take cumulative product marching forward from t0
        gpj_pos = cumprod(1 + gdj(tps > eps)) ;
        gpj_neg = flipud(cumprod(flipud(1 ./ (1 + gdj(tps < eps))))) ;
        gpj = cat(1, gpj_neg, gpj_pos) ;               

        % Plot this fold
        plot(tps, gpj, '.-', 'Color', QS.plotting.colors(jj, :))
        hold on;
    end    

    % Title and labels
    sgtitle(['Tissue deformation in folds, ', avgString, ', ', ...
        '$w_{\textrm{fold}}=', ...
        num2str(100*(2*width + 1)/ nU), '$\%$\, L_\zeta$'], ...
        'Interpreter', 'Latex')
    legend(foldYlabels, 'Interpreter', 'Latex', 'location', 'eastOutside')  
    drawnow
    xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
    
    switch lower(m2plot)
        case 'gdot'
            ylabel('$\Pi\big(1+$d$t\, \frac{1}{2}\mathrm{Tr} \left[g^{-1} \dot{g}\right]\big)$', ...
                    'Interpreter', 'Latex')
        case 'divv'
            ylabel('$\Pi\big(1+$d$t\,  \nabla \cdot \mathbf{v}_{\parallel} \big)$', ...
                    'Interpreter', 'Latex')
        case 'h2vn'
            ylabel('$\Pi\big(1+$d$t\, 2H v_n \big)$', ...
                    'Interpreter', 'Latex')
    end
    pos = get(gca, 'pos') ;

    % SECOND PANEL -- mean curvature
    subplot(2, 1, 2)
    for jj = 1:length(foldIds)
        valley = (foldIds(jj)-width):(foldIds(jj)+width) ;

        % Mark the instantaneous mean curvature
        Hj = mean(HH(:, valley), 2) ;
        scatter(tps(Hj>0), Hj(Hj>0), Hsz, 'filled', 's', ...
            'markeredgecolor', QS.plotting.colors(jj, :), ...
            'markerfacecolor', QS.plotting.colors(jj, :))
        hold on
        scatter(tps(Hj<0), Hj(Hj<0), Hsz, 's', ...
            'markeredgecolor', QS.plotting.colors(jj, :))

        foldHlabels{jj*2-1} = ['$H>0$, ' foldYlabels{jj} ] ;
        foldHlabels{jj*2} = ['$H<0$, ' foldYlabels{jj} ] ;

        % Identify time(s) at which H changes sign
        crossH{jj} = find(diff(sign(Hj)) < 0) ;
    end
    % Add indicators for changing sign
    for jj = 1:length(foldIds)
        for kk = 1:length(crossH{jj})
            subplot(2, 1, 1)
            ylims = get(gca, 'ylim') ;
            plot(tps(crossH{jj}(kk))*[1,1], ylims, '--', ...
                'color', QS.plotting.colors(jj, :), ...
                'HandleVisibility','off') ;
            ylim(ylims)
            subplot(2, 1, 2)
            Hlims = get(gca, 'ylim') ;
            plot(tps(crossH{jj}(kk))*[1,1], Hlims, '--', ...
                'color', QS.plotting.colors(jj, :), ...
                'HandleVisibility','off') ;
            ylim(Hlims)
        end
    end
    plot(tps, 0*tps, 'k--', 'HandleVisibility','off') ;
    legend(foldHlabels, 'Interpreter', 'Latex', 'location', 'eastOutside')  
    xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
    ylabel('fold curvature, $H$', 'Interpreter', 'Latex')
    pos2 = get(gca, 'pos') ;
    set(gca, 'pos', [pos(1) pos2(2) pos(3) pos(4)]) ;

    % Save figure with mean curvature
    disp(['Saving figure: ', fn_withH])
    saveas(gcf, fn_withH)
end