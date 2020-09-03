function aux_plotPathlineStrainFeatures_subpanels(QS, fn, fn_withH, ...
            featureIDs, width, nU, tps, trK, dvK, thK, HHK, ...
            avgString, titleFoldBase, foldYlabels, ...
            trace_label, deviator_label, ...
            trecolor, Hposcolor, Hnegcolor, Hsz, phasecmap, options)
%
% Parameters
% ----------
% 
%
% Returns
% -------
%
%
% NPMitchell 2020

% Default options 
overwrite = false ;
H_on_yyaxis = true ;
Halpha = 0.5 ;

% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'H_on_yyaxis')
    H_on_yyaxis = options.H_on_yyaxis ;
end
if isfield(options, 'Halpha')
    Halpha = options.Halpha ;
end

if ( ~exist(fn, 'file') || ~exist(fn_withH, 'file') || overwrite ) 
    close all
    ymin = 0 ;
    ymax = 0 ;
    % Each fold is valley+/- width
    for jj = 1:length(featureIDs)
        axisColl{jj} = subplot(length(featureIDs), 1, jj) ;
        valley = (featureIDs(jj)-width):(featureIDs(jj)+width) ;

        % trace, deviator, theta_deviator
        trj = mean(trK(:, valley), 2) ;
        [dvj, thj] = ...
            QS.dvAverageNematic(dvK(:, valley), thK(:, valley)) ;

        % Color deviator by theta
        plot(tps, trj, '.-', 'Color', trecolor)
        hold on;
        scatter(tps, dvj, 8, thj, 'filled')
        colormap(phasecmap)
        caxis([0, pi])

        if jj == 1
            % Title and labels
            sgtitle([titleFoldBase, avgString, ', ', ...
                '$w_{\textrm{fold}}=', ...
                num2str(100*(2*width + 1)/ nU), '$\%$\, L_\zeta$'], ...
                'Interpreter', 'Latex')
            legend({ trace_label, deviator_label}, ...
                'Interpreter', 'Latex', 'location', 'eastOutside')  
            drawnow
            pos = get(gca, 'position') ;    
        elseif jj == length(featureIDs)
            xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
        end
        ylabel(foldYlabels{jj}, 'Interpreter', 'Latex')
        ylims = ylim() ;
        ymin = min(ymin, ylims(1)) ;
        ymax = max(ymax, ylims(2)) ;
    end

    % Adjust axis positions 
    for jj = 1:length(featureIDs)
        axes(axisColl{jj}) ;
        pos2 = get(gca, 'position') ;
        set(gca, 'position', [pos2(1) pos2(2) pos(3) pos2(4)])
        ylim([ymin, ymax])
    end
    
    % Add phase colorbar
    axes(axisColl{1})
    phasebar('colormap', phasecmap, ...
            'location', [0.82, 0.5, 0.1, 0.135], 'style', 'nematic') ;
    
    % Save figure
    disp(['Saving figure: ', fn])
    saveas(gcf, fn)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Add mean curvature to each plot
    Hscale = 0 ;
    for jj = 1:length(featureIDs)
        valley = (featureIDs(jj)-width):(featureIDs(jj)+width) ;
        axes(axisColl{jj})

        % Put H on the same or a dual axis
        if H_on_yyaxis
            yyaxis right
        end
        
        % Mark the instantaneous mean curvature
        Hj = mean(HHK(:, valley), 2) ;
        scatter(tps(Hj>0), Hj(Hj>0), Hsz, 's', ...
            'markeredgecolor', Hposcolor, ...
            'MarkerEdgeAlpha', Halpha)
        scatter(tps(Hj<0), Hj(Hj<0), Hsz, 's', ...
            'markeredgecolor', Hnegcolor, ...
            'MarkerEdgeAlpha', Halpha)
        if ~H_on_yyaxis
            % Set ylimits either from mean curvature or from either H or strain
            ymin2 = min(min(Hj)-0.1*(ymax-ymin), ymin-0.1*(ymax-ymin)) ;
            ymax2 = max(ymax, max(Hj)+0.1*(ymax-ymin)) ;
            ylim([ymin2, ymax2])
        else
            Hscale = max(Hscale, max(abs(Hj))) ;
            yyaxis left
            ylims = ylim ;
            ylim([-max(abs(ylims)), max(abs(ylims))])
        end
    end
    
    % Make all mean curvature plots (yyaxes) same scale
    if H_on_yyaxis
        for jj = 1:length(featureIDs)
            axes(axisColl{jj})
            yyaxis right
            ylim([-Hscale, Hscale])
            ylabel('curvature, $H$', 'interpreter', 'latex')
            axisColl{jj}.YAxis(1).Color = 'k';
            axisColl{jj}.YAxis(2).Color = 'k';
        end
    end
    
    axes(axisColl{1})
    legend({trace_label, deviator_label, ...
        ['$H>0$ [' QS.spaceUnits '$^{-1}$]'], ...
        ['$H<0$ [' QS.spaceUnits '$^{-1}$]']}, 'Interpreter', 'Latex', ...
        'location', 'eastOutside')  
    
    % Save figure with mean curvature
    disp(['Saving figure: ', fn_withH])
    saveas(gcf, fn_withH)

end
