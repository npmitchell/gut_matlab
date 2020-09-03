function aux_plotPathlineStrainRegions(QS, ...
            fns, measurements, data, labels, options)
%AUX_PLOTPATHLINESTRAINFEATURES(QS, fns, measurements, data, labels,
%   options)
% Plot each feature's (fold's) strain on one axis.  
%
% Parameters
% ----------
% fn : str
%   path to output figure filename
% fons : #features x 1 numeric array
%   timepoints for the onset of features
% divv : 
% cumsum_cumprod : str ('cumsum' or 'cumprod')
%   Take cumulative sum (Euler integration) or cumulative product of
%   (1 + dt * Tr[epsilon_dot]) for integration
%
% NPMitchell 2020

% Default options
overwrite = false ;

% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'phasecmap')
    phasecmap = options.phasecmap ;
else
    phasecmap = phasemap(256) ;
end
% unpack fns
fn = fns.fn ;
fn_early = fns.early ;
fn_norms = fns.norms ;
fn_norms_early = fns.norms_early ;
fn_withH = fns.withH ;
% unpack measurements
fons = measurements.fons ;
regions = measurements.regions ;
tps = data.timepoints ;
trK = data.tr ;
dvK = data.dv ;
thK = data.th ;
HHK = data.HH ;
avgString = labels.avgString ;
titleBase = labels.titleBase ;
Ylabel = labels.ylabel ;
ratioLabel = labels.ylabel_ratios ;
foldYlabels = labels.legend  ;
foldYlabelsRatio = labels.legend_ratios ;
trace_label = labels.trace ;
deviator_label = labels.deviator ;

% Plot the data
if ~exist(fn, 'file') || overwrite 
    close all
    % Each fold is regions{jj}
    for jj = 1:length(regions)

        % trace, deviator
        trj = mean(trK(:, regions{jj}), 2) ;
        [dvj, thj] = ...
            QS.dvAverageNematic(dvK(:, regions{jj}), thK(:, regions{jj})) ;

        % Plot this fold
        plot(tps, trj, QS.plotting.markers{jj}, ...
            'Color', QS.plotting.colors(jj+3, :))
        hold on;
        scatter(tps, dvj, 10, thj, QS.plotting.markers{jj}, 'filled')
        colormap(phasecmap)
        caxis([0, pi])
    end    

    % Title and labels
    sgtitle([titleBase, avgString], 'Interpreter', 'Latex')
    legend(foldYlabels, 'Interpreter', 'Latex', 'location', 'eastOutside')  
    drawnow
    xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
    ylabel(Ylabel, 'Interpreter', 'Latex')
    
    % Adjust axis positions 
    phasebar('colormap', phasemap, ...
             'location', [0.82, 0.2, 0.1, 0.135], 'style', 'nematic')

    % Save figure
    disp(['Saving figure: ', fn])
    saveas(gcf, fn)
    xlim([min(tps), min(max(tps), max(fons) + 10)])
    disp(['Saving figure: ', fn_early])
    saveas(gcf, fn_early)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the magnitudes 
if ~exist(fn, 'file') || overwrite 
    close all
    maxr = 0 ;
    % Each fold is regions{jj}
    for jj = 1:length(regions)
        
        % trace, deviator
        trj = mean(trK(:, regions{jj}), 2) ;
        [dvj, ~] = ...
            QS.dvAverageNematic(dvK(:, regions{jj}), thK(:, regions{jj})) ;

        % Plot this fold
        plot(tps, abs(dvj) ./ abs(trj), [QS.plotting.markers{jj} '-'],...
            'Color', QS.plotting.colors(jj+3, :))
        maxr = max(maxr, max(abs(dvj) ./ abs(trj))) ;
        hold on;
    end    

    % Title and labels
    sgtitle([titleBase, avgString], 'Interpreter', 'Latex')
    legend(foldYlabelsRatio, 'Interpreter', 'Latex', 'location', 'eastOutside')  
    drawnow
    xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
    ylabel(ratioLabel, 'Interpreter', 'Latex')
    % Save figure
    disp(['Saving figure: ', fn_norms])
    ylim([0, min(10, maxr)])
    xlim([min(tps), max(tps)])
    saveas(gcf, fn_norms)
    xlim([min(tps), min(max(tps), max(fons) + 10)])
    disp(['Saving figure: ', fn_norms_early])
    saveas(gcf, fn_norms_early)
end