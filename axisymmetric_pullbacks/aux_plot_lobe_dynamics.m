function aux_plot_lobe_dynamics(length_lobes, area_lobes, volume_lobes, ...
    timePoints, fold_onset, colors, lobe_dynamics_figfn, ...
    lobedyn_figfn_scaled, t0, nlobes)
%AUX_PLOT_LOBE_DYNAMICS(length_lobes, area_lobes, volume_lobes, ...
%    timePoints, fold_onset, colors, lobe_dynamics_figfn)
% Auxiliary function for plotting the geometric dynamics of each lobe. 
% Also plots a scaled version where each geometric feature is divided by
% its initial value at t=min(fold_onset)
%
% Parameters
% ----------
% length_lobes : 
% area_lobes : 
% volume_lobes : 
% timePoints : #timepoints x 1 increasing float or int array
%   the timepoints of the dynamic dataset, in minutes
% fold_onset : 3x1 float array
%   the timepoints in minutes at which folds appear in the sample
% colors : 
% lobe_dynamics_figfn : 
% t0 : int 
%   timepoint to set to zero time
% flipy : bool
%   data is reflected wrt lab frame
% nlobes : int
%
% Returns
% -------
% <none>
%
% NPMitchell 2020 

close all;
fig = figure('visible', 'off');
fh = cell(1, nlobes) ;
tp = timePoints - t0 ; 
for lobe = 1:nlobes
    % Length
    subplot(3,2,1)
    plot(tp, length_lobes(:, lobe), '.', 'Color', colors(lobe, :)) ;
    hold on

    % Area
    subplot(3,2,3)
    fh{lobe} = plot(tp, area_lobes(:, lobe), '.', 'Color', colors(lobe, :)) ;
    hold on

    % Volume
    subplot(3,2,5)
    plot(tp, volume_lobes(:, lobe), '.', 'Color', colors(lobe, :)) ;
    hold on
end

subplot(3, 2, 1)
xlim([min(tp), max(tp)])
ylabel('Length [\mum]')

subplot(3, 2, 3)
xlim([min(tp), max(tp)])
ylabel('Area [\mum^2]')

if nlobes == 4 
    legend({'\color[rgb]{ 0,0.4470,0.7410} lobe 1', ...
        '\color[rgb]{0.8500,0.3250,0.0980} lobe 2', ...
        '\color[rgb]{0.9290,0.6940,0.1250} lobe 3', ...
        '\color[rgb]{0.4940,0.1840,0.5560} lobe 4'}, 'Position', [0.55 0.4 0.1 0.2])
elseif nlobes == 3
    legend({'\color[rgb]{ 0,0.4470,0.7410} lobe 1', ...
        '\color[rgb]{0.8500,0.3250,0.0980} lobe 2', ...
        '\color[rgb]{0.9290,0.6940,0.1250} lobe 3'}, 'Position', [0.55 0.4 0.1 0.2])    
else
    error('have not coded for this number of lobes')
end

subplot(3, 2, 5)
xlim([min(tp), max(tp)])
ylabel('Volume [\mum^3]')

xlabel('time [min]')
disp(['Saving summary to ' lobe_dynamics_figfn])
saveas(fig, lobe_dynamics_figfn)
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plotting lobe dynamics scaled...')
close all;
fig = figure('visible', 'off');
fh = cell(1, 4) ;
tp = timePoints - min(fold_onset) ;
for lobe = 1:nlobes
    % Length
    subplot(3,2,1)
    y = length_lobes(:, lobe) / length_lobes(tp==0, lobe) ;
    plot(tp, y, '.', 'Color', colors(lobe, :)) ;
    hold on

    % Area
    subplot(3,2,3)
    y = area_lobes(:, lobe)/ area_lobes(tp==0, lobe) ;
    fh{lobe} = plot(tp, y, '.', 'Color', colors(lobe, :)) ;
    hold on

    % Volume
    subplot(3,2,5)
    y = volume_lobes(:, lobe) / volume_lobes(tp==0, lobe) ;
    plot(tp, y, '.', 'Color', colors(lobe, :)) ;
    hold on
end

subplot(3, 2, 1)
xlim([min(tp), max(tp)])
ylabel('Length / L_0')
ylim([0, 3])

subplot(3, 2, 3)
xlim([min(tp), max(tp)])
ylabel('Area / A_0')
legend({'\color[rgb]{ 0,0.4470,0.7410} lobe 1', ...
    '\color[rgb]{0.8500,0.3250,0.0980} lobe 2', ...
    '\color[rgb]{0.9290,0.6940,0.1250} lobe 3', ...
    '\color[rgb]{0.4940,0.1840,0.5560} lobe 4'}, 'Position', [0.55 0.4 0.1 0.2])
% ylim([0, 2])

subplot(3, 2, 5)
xlim([min(tp), max(tp)])
ylabel('Volume / V_0')
ylim([0, 2])

xlabel('time [min]')
disp(['Saving summary to ' lobedyn_figfn_scaled])
saveas(fig, lobedyn_figfn_scaled)