function aux_plotPathlineStrainRegions_subpanels(QS, fns, ...
            regions, data, labels, plotOptions)
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
tp_select = [] ;
trecolor = 'k' ;
Hposcolor = [0.4660    0.6740    0.1880] ; % green
Hnegcolor = [0.4940    0.1840    0.5560] ; % purple
Hsz = 3 ;  % size of scatterplot dots for mean curvature in yyaxis right
phasecmap = phasemap(256) ;  % colormap for nematic angle

% Unpack fns
fn = fns.fn ;
fn_withH = fns.fn_withH ;
fn_selectTimes = fns.fn_selectTimes ; 
% unpack data
tps = data.timepoints ;
trK = data.tr ;  % note this is 1/2 * tr(dot(epsilon))
dvK = data.dv ;
thK = data.th ;
HHK = data.HH ;
tp_select = data.selectTimes ;
% unpack labels
avgString = labels.avgString ;
titleBase = labels.titleBase ;
Ylabels = labels.ylabels  ;
trace_label = labels.trace ;
deviator_label = labels.deviator ;

% Unpack plotOptions
if isfield(plotOptions, 'overwrite')
    overwrite = plotOptions.overwrite ;
end
if isfield(plotOptions, 'H_on_yyaxis')
    H_on_yyaxis = plotOptions.H_on_yyaxis ;
end
if isfield(plotOptions, 'Halpha')
    Halpha = plotOptions.Halpha ;
end
if isfield(plotOptions, 'selectTimes')
    tp_select = plotOptions.selectTimes ;
end
if isfield(plotOptions, 'trecolor')
    trecolor = plotOptions.trecolor ;
end
if isfield(plotOptions, 'Hposcolor')
    Hposcolor = plotOptions.Hposcolor ;
end
if isfield(plotOptions, 'Hnegcolor')
    Hnegcolor = plotOptions.Hnegcolor ;
end
if isfield(plotOptions, 'Hsz')
    Hsz = plotOptions.Hsz ;
end
if isfield(plotOptions, 'phasecmap')
    phasecmap = plotOptions.phasecmap ;
end

% Create plots
if ( ~exist(fn, 'file') || ~exist(fn_withH, 'file') || overwrite ) 
    close all
    set(gcf, 'visible', 'off')
    ymin = 0 ;
    ymax = 0 ;
    % Each fold is valley+/- width
    for jj = 1:length(regions)
        axisColl{jj} = subplot(length(regions), 1, jj) ;
        
        % trace, deviator, theta_deviator
        trj = mean(trK(:, regions{jj}), 2) ;
        [dvj, thj] = ...
            QS.dvAverageNematic(dvK(:, regions{jj}), thK(:, regions{jj})) ;

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
        elseif jj == length(regions)
            xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
        end
        ylabel(Ylabels{jj}, 'Interpreter', 'Latex')
        ylims = ylim() ;
        ymin = min(ymin, ylims(1)) ;
        ymax = max(ymax, ylims(2)) ;
    end

    % Adjust axis positions 
    for jj = 1:length(regions)
        set(gcf, 'currentAxes', axisColl{jj})
        pos2 = get(gca, 'position') ;
        set(gca, 'position', [pos2(1) pos2(2) pos(3) pos2(4)])
        ylim([ymin, ymax])
    end
    
    set(gcf, 'currentAxes', axisColl{1})
    phasebar('colormap', phasecmap, ...
             'location', [0.82, 0.5, 0.1, 0.135], 'style', 'nematic')

    % Save figure
    disp(['Saving figure: ', fn])
    saveas(gcf, fn)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Add mean curvature to each plot
    Hscale = 0 ;
    if ~isempty(fn_withH) || ~strcmpi(fn_withH, 'none')
        for jj = 1:length(regions)
            set(gcf, 'currentAxes', axisColl{jj})
            if H_on_yyaxis
                yyaxis right
            end
            
            % Mark the instantaneous mean curvature
            Hj = mean(HHK(:, regions{jj}), 2) ;
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
            for jj = 1:length(regions)
                set(gcf, 'currentAxes', axisColl{jj})
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

        set(gcf, 'currentAxes', axisColl{1})
        legend({trace_label, deviator_label, ...
            ['$H>0$ [' QS.spaceUnits '$^{-1}$]'], ...
            ['$H<0$ [' QS.spaceUnits '$^{-1}$]']}, 'Interpreter', 'Latex', ...
            'location', 'eastOutside')  

        % Save figure with mean curvature
        disp(['Saving figure: ', fn_withH])
        saveas(gcf, fn_withH)
    end

    %% Save figure with select times only (typically early times)
    if ~isempty(tp_select)
        for jj = 1:length(regions)
            set(gcf, 'currentAxes', axisColl{jj}) ;
            xlim([min(tp_select), max(tp_select)])
        end
        disp(['Saving figure (select times): ', fn_selectTimes])
        saveas(gcf, fn_selectTimes)
    end
end