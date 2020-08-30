function aux_plotPathlineMetricKinematicsFolds_integrated_subpanels(QS, ...
    fn, foldIds, width, nU, tps, divv, H2vn, titleFoldBase, foldYlabels, ...
    avgString, divvcolor, H2vncolor, gdotcolor, overwrite, cumsum_cumprod)
% aux_plotMetricKinematicsFold_integrated_subpanels(fn, ...
%     fn_withH, foldIds, width, divv, H2vn, titleFoldBase, avgString)
%
% Parameters
% ----------
% fn : str
% fn_withH : str
% cumsum_cumprod : str ('cumsum' or 'cumprod')
%   Take cumulative sum (Euler integration) or cumulative product of
%   (1 + dt * Tr[epsilon_dot]) for integration
%
%
% NPMitchell 2020

close all
set(gcf, 'visible', 'off')
ymin = 0 ;
ymax = 0 ;
% Each fold is valley+/- width
if ~exist(fn, 'file') || overwrite
    for jj = 1:length(foldIds)
        axisColl{jj} = subplot(length(foldIds), 1, jj) ;
        valley = (foldIds(jj)-width):(foldIds(jj)+width) ;

        % div(v), H*2*vn, gdot
        ddj = mean(divv(:, valley), 2) ;
        Hdj = mean(H2vn(:, valley), 2) ;
        gdj = mean(divv(:, valley) - H2vn(:, valley), 2) ;

        if strcmpi(cumsum_cumprod, 'cumsum')
            % Take cumulative product marching forward from t0
            gpj_pos = cumsum(gdj(tps > eps)) ;
            gpj_neg = flipud(cumsum(flipud(-gdj(tps < eps)))) ;
            gpj = cat(1, gpj_neg, gpj_pos) ;               
            % Take cumulative product marching forward from t0
            dpj_pos = cumsum(ddj(tps > eps)) ;
            dpj_neg = flipud(cumsum(flipud(-ddj(tps < eps)))) ;
            dpj = cat(1, dpj_neg, dpj_pos) ;
            % Take cumulative product marching forward from t0
            Hpj_pos = cumsum(Hdj(tps > eps)) ;
            Hpj_neg = flipud(cumsum(flipud(-Hdj(tps < eps)))) ;
            Hpj = cat(1, Hpj_neg, Hpj_pos) ;
        else
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
        end
        
        % Plot all three
        plot(tps, dpj, '.-', 'Color', divvcolor)
        hold on;
        plot(tps, Hpj, '.-', 'Color', H2vncolor)
        plot(tps, gpj, '.-', 'Color', gdotcolor)

        if jj == 1
            % Title and labels
            sgtitle([titleFoldBase, avgString, ', ', ...
                '$w_{\textrm{fold}}=', ...
                num2str(100*(2*width + 1)/ nU), '$\%$\, L_\zeta$'], ...
                'Interpreter', 'Latex')
            
            if strcmp(cumsum_cumprod, 'cumsum')
                legend({'$\int $d$t \, \nabla\cdot\mathbf{v}_\parallel$', ...
                    '$\int $d$t \, v_n 2H$', ...
                    '$\int $d$t \, \frac{1}{2}\mathrm{Tr}\left[g^{-1} \dot{g} \right]$'}, ...
                    'Interpreter', 'Latex', 'location', 'eastOutside')  
            else
                legend({'$\Pi(1+\nabla\cdot\mathbf{v}_\parallel)$', ...
                    '$\Pi(1+v_n 2H)$', ...
                    '$\Pi(1+\frac{1}{2}\mathrm{Tr}\left[g^{-1} \dot{g} \right])$'}, ...
                    'Interpreter', 'Latex', 'location', 'eastOutside')                  
            end
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

        % % Mark wherever the divv<0 and also vn2H<0
        % valley = (foldIds(jj)-width):(foldIds(jj)+width) ;
        % dvj = mean(divv(:, valley), 2) ;
        % Hvj = mean(H2vn(:, valley), 2) ;
        % dvpos = dvj > 0 ;
        % Hvpos = Hvj > 0 ;
        % scatter(tps(dvpos), ...
        %     (ymin-0.05*(ymax-ymin)) * ones(size(tps(dvpos))), 5, ...
        %     'markeredgecolor', 'none', 'markerFaceAlpha', 0.6, ...
        %     'markerFaceColor', divvcolor, ...
        %     'HandleVisibility', 'off')
        % hold on;
        % scatter(tps(Hvpos), ...
        %     ymin * ones(size(tps(Hvpos))), 5,  ...
        %     'markeredgecolor', 'none', 'markerFaceAlpha', 0.6, ...
        %     'markerFaceColor', H2vncolor, ...
        %     'HandleVisibility', 'off')
        % ylim([ymin-0.1*(ymax-ymin), ymax])
    end

    % Save figure
    disp(['Saving figure: ', fn])
    saveas(gcf, fn)
end