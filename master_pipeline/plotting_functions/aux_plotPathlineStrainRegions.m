function aux_plotPathlineStrainRegions(QS, ...
            fns, measurements, data, labels, options)
%AUX_PLOTPATHLINESTRAINFEATURES(QS, fns, measurements, data, labels,
%   options)
% Plot each feature's (fold's) strain on one axis.  
% NOTE: trace is passed to this function as 1/2*tr(epsilon), not tr(eps)
%
% Parameters
% ----------
% fns : struct with fields
%   fn : str
%       path to output figure filename
%   fn_early : str
%       path to output figure filename (zoom in x axis)
%   norms_early : str
%       path to output figure filename (zoom in x axis)
%   fn_withH : str
%       path to output figure filename with mean curvature on yyaxis right
% data : struct with fields
%   timepoints : Nx1 numeric array
%       timepoints indexing the dynamics being plotted
%   tr : Nx1 float array
%       1/2 * trace of epsilon (half the trace of the strain rate), ie mean
%       on-diagonal strain rate
%   dv : Nx1 float array
%   th : Nx1 float array
%   HH : Nx1 float array
% fons : #features x 1 numeric array
%   timepoints for the onset of features
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
% unpack data
tps = data.timepoints ;
trK = data.tr ;  % note this is 1/2 * tr(dot(epsilon))
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
tp_early = [max(-30,min(tps)), min(max(tps), max(120))] ;

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
        plot(tps, trj, QS.plotting.markers{jj}, 'markersize', 4, 'color', 'k')
        %    'Color', QS.plotting.colors(jj+3, :))
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
    xlim(tp_early)
    disp(['Saving figure: ', fn_early])
    saveas(gcf, fn_early)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the magnitudes 
if ~exist(fn, 'file') || overwrite 
    close all
    set(gcf, 'visible', 'off')
    maxr = 0 ;
    minr = 0 ;
    % Each fold is regions{jj}
    for jj = 1:length(regions)
        
        % trace, deviator
        trj = mean(trK(:, regions{jj}), 2) ;
        [dvj, ~] = ...
            QS.dvAverageNematic(dvK(:, regions{jj}), thK(:, regions{jj})) ;
        trj = movmean(trj, 3) ;
        dvj = movmean(dvj, 3) ;

        % Plot this fold -- note that trK is supplied 1/2*tr(epsilondot)
        plot(tps, medfilt1m(trj ./ abs(dvj),5), [QS.plotting.markers{jj} '-'],...
            'Color', QS.plotting.colors(jj+3, :))
        maxr = max(maxr, max(trj ./ abs(dvj))) ;
        minr = min(minr, min(trj ./ abs(dvj))) ;
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
    ylim([max(-2, minr), min(2, maxr)])
    plot([min(tps), max(tps)], [1,1], 'k--')
    plot([min(tps), max(tps)], [-1,-1], 'k--')
    xlim([min(tps), max(tps)])
    saveas(gcf, fn_norms)
    xlim(tp_early)
    disp(['Saving figure: ', fn_norms_early])
    saveas(gcf, fn_norms_early)
end