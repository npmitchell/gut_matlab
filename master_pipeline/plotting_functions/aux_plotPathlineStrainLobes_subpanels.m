function aux_plotPathlineStrainLobes_subpanels(QS, fn, fn_withH, ...
            lobes, tps, trK, dvK, thK, HHK, ...
            avgString, titleBase, lobeYlabels, ...
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

% Create plots
if ( ~exist(fn, 'file') || ~exist(fn_withH, 'file') || overwrite ) 
    close all
    ymin = 0 ;
    ymax = 0 ;
    % Each fold is valley+/- width
    for jj = 1:length(lobes)
        axisColl{jj} = subplot(length(lobes), 1, jj) ;
        
        % trace, deviator, theta_deviator
        trj = mean(trK(:, lobes{jj}), 2) ;
        [dvj, thj] = ...
            QS.dvAverageNematic(dvK(:, lobes{jj}), thK(:, lobes{jj})) ;

        % Color deviator by theta
        plot(tps, trj, '.-', 'Color', trecolor)
        hold on;
        scatter(tps, dvj, 8, thj, 'filled')
        colormap(phasecmap)
        caxis([0, pi])

        if jj == 1
            % Title and labels
            sgtitle([titleBase, avgString], ...
                'Interpreter', 'Latex')
            legend({ trace_label, deviator_label}, ...
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

    % Adjust axis positions 
    for jj = 1:length(lobes)
        axes(axisColl{jj})
        pos2 = get(gca, 'position') ;
        set(gca, 'position', [pos2(1) pos2(2) pos(3) pos2(4)])
        ylim([ymin, ymax])
    end
    
    axes(axisColl{1})
    phasebar('colormap', phasecmap, ...
             'location', [0.82, 0.5, 0.1, 0.135], 'style', 'nematic')

    % Save figure
    disp(['Saving figure: ', fn])
    saveas(gcf, fn)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Add mean curvature to each plot
    Hscale = 0 ;
    if length(fn_withH) > 0 || ~strcmpi(fn_withH, 'none')
        for jj = 1:length(lobes)
            axes(axisColl{jj})
            if H_on_yyaxis
                yyaxis right
            end
            
            % Mark the instantaneous mean curvature
            Hj = mean(HHK(:, lobes{jj}), 2) ;
            scatter(tps(Hj>0), Hj(Hj>0), Hsz, 's', ...
                'markeredgecolor', Hposcolor, ...
                'MarkerEdgeAlpha', Halpha)
            scatter(tps(Hj<0), Hj(Hj<0), Hsz, 's', ...
                'markeredgecolor', Hnegcolor, ...
                'MarkerEdgeAlpha', Halpha)
            
            % Set ylims so each axis is proportional
            if H_on_yyaxis
                Hscale = max(Hscale, max(abs(Hj))) ;
            else                
                ymin2 = min(min(Hj)-0.1*(ymax-ymin), ymin-0.1*(ymax-ymin)) ;
                ymax2 = max(ymax, max(Hj)+0.1*(ymax-ymin)) ;
                ylim([ymin2, ymax2])
            end
        end
        
        % Make all mean curvature plots (yyaxes) same scale
        if H_on_yyaxis
            for jj = 1:length(lobes)
                axes(axisColl{jj})
                yyaxis right
                ylims = ylim() ;
                ylim([-ylims(1), Hscale])
                set_proportional_yylims(0, 'right')
                yyaxis right
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
end