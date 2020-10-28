%% Compare Flow Fields from All Datasets
% Isaac B Breinyn 2020
%
% This is a script that reads in all the data from comparing flow fields
% of various datasets and combines them in a single workspace. All data for
% this script should already exist. If some data is missing from a
% particular dataset, run the compare_flow_fields_____.m script that
% corresponds to that dataset. The only manual input in this script is time
% offset between datasets (calculated by qualitatively extracting fold
% depth using ImageJ)
%
%% Clean up Matlab and define Datasets

clear; close all; clc ;

% We now have to collect all of the relevant information and data necessary
% to creating a single workspace. We will do this by organizing the
% following information into a structure with fields: data name, data set
% location, binned velocities, relative velocities, the full flow data, the
% phases, and the normed dot products. All of these are stored in the
% datasets' file directories.

% Where should master plots go?
data.masterDir = 'L:\\Streichan\\data\\master_plots\\' ;

% First, the base name of each dataset (will pop up in naming conventions)
data.dset{1} = '1644_folding' ;
data.dset{2} = '2038_e5' ;
data.dset{3} = 'midfold_e4' ;
data.dset{4} = '1006_e2' ;
data.dset{5} = '2101_e1' ;

% Next, the directories of each dataset
data.dir{1} = 'L:\\Streichan\\data\\202003111630_mef2gal4klarUASCAAXmChHiFP_wo_great\\1644_folding_30s_2um_la8_4p5zoom_63x\\' ;
data.dir{2} = 'L:\\Streichan\\data\\202003111830_mef2gal4klarUASCAAXmChHiFP_wo_63x_e4e5\\2038_e5\\' ;
data.dir{3} = 'L:\\Streichan\\data\\202003111830_mef2gal4klarUASCAAXmChHiFP_wo_63x_e4e5\\midfold_e4\\' ;
data.dir{4} = 'L:\\Streichan\\data\\202007142100_mef2GAL4klarH2AGFPCAAXmCh_e1e2\\1006_e2\\' ;
data.dir{5} = 'L:\\Streichan\\data\\202007142100_mef2GAL4klarH2AGFPCAAXmCh_e1e2\\2101_e1\\' ;

% Just to save the information, let's store the (t,x,y,z) resolution
data.res{1} = [30 .0802 .0802 1.9994] ;
data.res{2} = [120 .0401 .0401 1.4996] ;
data.res{3} = [60 .0833 .0833 .9999] ;
data.res{4} = [60 .1137 .1137 1.4996] ;
data.res{5} = [60 .1264 .1264 1.4996] ;

% Save the last time point of each dataset
data.ltp{1} = 32 ;
data.ltp{2} = 24 ;
data.ltp{3} = 33 ;
data.ltp{4} = 60 ;
data.ltp{5} = 58 ;

for dd = 1:length(data.dset)
    try
        data.binvels{dd} = importdata([data.dir{dd} 'binvels.m']) ; % The binned velocities
        data.binvels2{dd} = importdata([data.dir{dd} 'binvels2.m']) ; % The two lobe binned velocities
        data.relvels{dd} = importdata([data.dir{dd} 'relvels.m']) ; % The relative velocities
        data.relvels2{dd} = importdata([data.dir{dd} 'relvels2.m']) ; % The two lobe relative velocities
        data.fullflow{dd} = importdata([data.dir{dd} 'fullflow.m']) ; % The full flow data (as outputted by PIVlab)
        data.phases{dd} = importdata([data.dir{dd} 'phases.m']) ; % The phases
        data.phases2{dd} = importdata([data.dir{dd} 'phases2.m']) ; % The two lobe phases
        data.surfaceArray{dd} = importdata([data.dir{dd} 'surfaceArray.m']) ;% The surface array (used to find time offsets)
        data.normdots{dd} = importdata([data.dir{dd} 'normdots.m']) ; % The normed dot products
        data.normdots2{dd} = importdata([data.dir{dd} 'normdots2.m']) ; % The two lobe normed dot products
    catch
        disp(['Data missing from dataset ' data.dset{dd} ', please rerun compare_flow_MASTER.m on it now']) ;
    end
end

%% Load data into RAM to plot

% phases
for dd = 1:length(data.dset)
    if dd == 1
        maxfsize = size(data.phases{dd}, 1) ;
    else
        if size(data.phases{dd}, 1) > maxfsize
            maxfsize = size(data.phases{dd}, 1) ;
        end
    end
end

for dd = 1:length(data.dset)
    data.phases{dd} = padarray(data.phases{dd}, [maxfsize - size(data.phases{dd}, 1) 0], NaN, 'post') ;
end

phases = cat(2, data.phases{:}) ;
clear maxfsize

% 2 lobe phases
for dd = 1:length(data.dset)
    if dd == 1
        maxfsize = size(data.phases2{dd}, 2) ;
    else
        if size(data.phases2{dd}, 2) > maxfsize
            maxfsize = size(data.phases2{dd}, 2) ;
        end
    end
end

for dd = 1:length(data.dset)
    data.phases2{dd} = padarray(data.phases2{dd}, [0, maxfsize - size(data.phases2{dd}, 2), 0], NaN, 'post') ;
end

phases2 = cat(3, data.phases2{:}) ;
clear maxfsize

% normed dot products
for dd = 1:length(data.dset)
    if dd == 1
        maxfsize = size(data.normdots{dd}, 1) ;
    else
        if size(data.normdots{dd}, 1) > maxfsize
            maxfsize = size(data.normdots{dd}, 1) ;
        end
    end
end

for dd = 1:length(data.dset)
    data.normdots{dd} = padarray(data.normdots{dd}, [maxfsize - size(data.normdots{dd}, 1), 0], NaN, 'post') ;
end

normdots = cat(2, data.normdots{:}) ;
clear maxfsize

% two lobe normed dot products
for dd = 1:length(data.dset)
    if dd == 1
        maxfsize = size(data.normdots2{dd}, 2) ;
    else
        if size(data.normdots2{dd}, 2) > maxfsize
            maxfsize = size(data.normdots2{dd}, 2) ;
        end
    end
end

for dd = 1:length(data.dset)
    data.normdots2{dd} = padarray(data.normdots2{dd}, [0, maxfsize - size(data.normdots2{dd}, 2), 0], NaN, 'post') ;
end

normdots2 = cat(3, data.normdots2{:}) ;
clear maxfsize

% binned velocities
for dd = 1:length(data.dset)
    if dd == 1
        maxfsize = size(data.binvels{dd}, 2) ;
    else
        if size(data.binvels{dd}, 2) > maxfsize
            maxfsize = size(data.binvels{dd}, 2) ;
        end
    end
end

for dd = 1:length(data.dset)
    data.binvels{dd} = padarray(data.binvels{dd}, [0, maxfsize - size(data.binvels{dd}, 2), 0], NaN, 'post') ;
end

binvels = cat(3, data.binvels{:}) ;
clear maxfsize

% two lobe binned velocities
for dd = 1:length(data.dset)
    if dd == 1
        maxfsize = size(data.binvels2{dd}, 3) ;
    else
        if size(data.binvels2{dd}, 3) > maxfsize
            maxfsize = size(data.binvels2{dd}, 3) ;
        end
    end
end

for dd = 1:length(data.dset)
    data.binvels2{dd} = padarray(data.binvels2{dd}, [0, 0, maxfsize - size(data.binvels2{dd}, 3), 0], NaN, 'post') ;
end

binvels2 = cat(4, data.binvels2{:}) ;
clear maxfsize

% relative velocities
for dd = 1:length(data.dset)
    if dd == 1
        maxfsize = size(data.relvels{dd}, 2) ;
    else
        if size(data.relvels{dd}, 2) > maxfsize
            maxfsize = size(data.relvels{dd}, 2) ;
        end
    end
end

for dd = 1:length(data.dset)
    data.relvels{dd} = padarray(data.relvels{dd}, [0, maxfsize - size(data.relvels{dd}, 2), 0], NaN, 'post') ;
end

relvels = cat(3, data.relvels{:}) ;
clear maxfsize

% two lobe relative velocities
for dd = 1:length(data.dset)
    if dd == 1
        maxfsize = size(data.relvels2{dd}, 3) ;
    else
        if size(data.relvels2{dd}, 3) > maxfsize
            maxfsize = size(data.relvels2{dd}, 3) ;
        end
    end
end

for dd = 1:length(data.dset)
    data.relvels2{dd} = padarray(data.relvels2{dd}, [0, 0, maxfsize - size(data.relvels2{dd}, 3), 0], NaN, 'post') ;
end

relvels2 = cat(4, data.relvels2{:}) ;
clear maxfsize
%% Create Master Plots (overall and considering lobes seperately)
% The master plot will include a 95% confidence interval ellipse.
% Code for that was taken from Vision Dummy.

close all
lin = linspace(-1, 1, 3) ;
colors = zeros(data.ltp{dd}-1, 3) ;

for dd = 1:length(data.dset)
    if dd == 1
        maxltp = data.ltp{dd} ;
    elseif data.ltp{dd} > maxltp
        maxltp = data.ltp{dd} ;
    end
end

phases(phases > pi)  = phases(phases > pi) - 2*pi ;
phases2(phases2 > pi)  = phases2(phases2 > pi) - 2*pi ;

figure('units','normalized','outerposition', [0 0 1 1]) ;
for tt = 1:maxltp-1
    for dd = 1:length(data.dset)
        if tt <= size(data.binvels{dd}, 2)
            cmap = parula(maxltp) ;
            color = cmap(tt, :) ;
            colors(tt,:) = color ;
            subplot(2, 3, 1)
            title('u-comp velocities') ;
            ylabel('V_u for Membrane data [um/min]') ;
            xlabel('V_u for Nuclear data [um/min]') ;
            hold on
            uu1 = squeeze(data.binvels{dd}(1, tt, :)) ;
            uu2 = squeeze(data.binvels{dd}(3, tt,:)) ;
            uu1 = filloutliers(uu1, 'nearest', 'percentiles', [15 85]) ;
            uu2 = filloutliers(uu2, 'nearest', 'percentiles', [15 85]) ;
            plot(lin, lin, '--')
            scatter(uu2, uu1, 5, color, 'filled', 'markeredgecolor', 'none') ;
            xlim([-1 1])
            ylim([-1 1])
            axis equal
            
            subplot(2, 3, 2)
            title('v-comp velocities') ;
            ylabel('V_v for Membrane data [um/min]') ;
            xlabel('V_v for Nuclear data [um/min]') ;
            hold on
            vv1 = squeeze(data.binvels{dd}(2, tt, :)) ;
            vv2 = squeeze(data.binvels{dd}(4, tt,:)) ;
            vv1 = filloutliers(vv1, 'nearest', 'percentiles', [15 85]);
            vv2 = filloutliers(vv2, 'nearest', 'percentiles', [15 85]) ;
            plot(lin, lin, '--')
            scatter(vv2, vv1, 5, color, 'filled', 'markeredgecolor', 'none') ;
            xlim([-1 1])
            ylim([-1 1])
            axis equal
            
            subplot(2, 3, 3)
            title('Relative motion') ;
            ylabel('V_v_,_r_e_l') ;
            xlabel('V_u_,_r_e_l') ;
            hold on
            scatter(squeeze(data.relvels{dd}(1, tt, :)), squeeze(data.relvels{dd}(2, tt, :)), 5, color, 'filled', 'markeredgecolor', 'none') ;
            xlim([-1 1])
            ylim([-1 1])
        end
    end
end

subplot(2, 3, 4)
title('Relative Phase') ;
ylabel('Count') ;
xlabel('Phase') ;
xlim([-pi pi]);
hold on
N21 = zeros(maxltp-1, 9) ;
for tt = 1:maxltp-1
    phase = phases(tt, phases(tt, :) ~= 0) ;
    [~,edges] = histcounts(phase, linspace(-pi, pi, 10));
    N21(tt, :) = histcounts(phase ,edges); % Bin using the same edges
end
ctrs1 = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
b1 = bar(ctrs1, N21, 'stacked', 'edgecolor', 'none') ;
for kk = 1:size(N21,1)
    cmap = parula(maxltp) ;
    color = cmap(kk, :) ;
    b1(kk).FaceColor = color ;
end

subplot(2, 3, 5)
title('Relative Magnitude |V_m_e_m|/|V_n_u_c|') ;
ylabel('Count') ;
xlabel('Relative Magnitude') ;
relmags = zeros(maxltp-1, size(binvels, 3)) ;
N22 = zeros(maxltp-1, 9) ;
for tt = 1:maxltp-1
    relmags(tt, :) = sqrt((squeeze(binvels(1, tt, :))).^2 + (squeeze(binvels(2, tt, :))).^2)./sqrt((squeeze(binvels(3, tt, :))).^2 + (squeeze(binvels(4, tt, :))).^2) ;
    hold on
    relmag = relmags(tt, relmags(tt,:) ~=0) ;
    [~,edges] = histcounts(relmag, linspace(0, 2, 10));
    N22(tt, :) = histcounts(relmag ,edges); % Bin using the same edges
end
ctrs2 = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
b2 = bar(ctrs2, N22, 'stacked', 'FaceColor', 'flat', 'edgecolor', 'none') ;
for kk = 1:size(N22, 1)
    cmap = parula(maxltp) ;
    color = cmap(kk, :) ;
    b2(kk).CData = color ;
end

subplot(2, 3, 6)
title('Normed Dot Products') ;
ylabel('Count') ;
xlabel('Normed Dot Product Values');
N23 = zeros(maxltp-1, 9) ;
for tt = 1:maxltp-1
    normdot = normdots(tt, normdots(tt,:) ~=0) ;
    [~,edges] = histcounts(normdot, linspace(-1, 1, 10));
    N23(tt, :) = histcounts(normdot ,edges); % Bin using the same edges
end
ctrs3 = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
b3 = bar(ctrs3, N23, 'stacked', 'FaceColor', 'flat', 'edgecolor', 'none') ;
for kk = 1:maxltp-1
    cmap = parula(maxltp) ;
    color = cmap(kk, :) ;
    b3(kk).CData = color ;
end

sgtitle('Master Correlation Plot') ;

% relvelucomp = squeeze(data.relvels{dd}(1,:,:)) ;
% relvelvcomp = squeeze(data.relvels{dd}(2,:,:)) ;
% data = horzcat(relvelucomp(:),relvelvcomp(:)) ;
% [r_ellipse, X0, Y0] = error_ellipse(data) ;
% subplot(2, 3, 3)
% p = plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'color', 'r','linestyle', '-') ;
% legend(p, '95% Confidence Interval') ;

saveas(gcf, fullfile(data.masterDir, 'master_corr.png')) ;
%%
for ll = 1:2
    
    close all
    colors = zeros(data.ltp{dd}-1, 3) ;
    
    figure('units','normalized','outerposition', [0 0 1 1]) ;
    for tt = 1:maxltp-1
        for dd = 1:length(data.dset)
            if tt <= size(data.binvels2{dd}, 3)
                cmap = parula(maxltp) ;
                color = cmap(tt, :) ;
                colors(tt,:) = color ;
                subplot(2, 3, 1)
                title('u-comp velocities') ;
                ylabel('V_u for Membrane data [um/min]') ;
                xlabel('V_u for Nuclear data [um/min]') ;
                hold on
                uu1 = squeeze(data.binvels2{dd}(ll, 1, tt, :)) ;
                uu2 = squeeze(data.binvels2{dd}(ll, 3, tt,:)) ;
                uu1 = filloutliers(uu1, 'nearest', 'percentiles', [15 85]) ;
                uu2 = filloutliers(uu2, 'nearest', 'percentiles', [15 85]) ;
                uu1(uu1 == 0) = NaN ;
                uu2(uu2 == 0) = NaN ;
                plot(lin, lin, '--')
                scatter(uu2, uu1, 5, color, 'filled', 'markeredgecolor', 'none') ;
                xlim([-1 1])
                ylim([-1 1])
                axis equal
                
                subplot(2, 3, 2)
                title('v-comp velocities') ;
                ylabel('V_v for Membrane data [um/min]') ;
                xlabel('V_v for Nuclear data [um/min]') ;
                hold on
                vv1 = squeeze(data.binvels2{dd}(ll, 2, tt, :)) ;
                vv2 = squeeze(data.binvels2{dd}(ll, 4, tt,:)) ;
                vv1 = filloutliers(vv1, 'nearest', 'percentiles', [15 85]);
                vv2 = filloutliers(vv2, 'nearest', 'percentiles', [15 85]) ;
                vv1(vv1 == 0) = NaN ;
                vv2(vv2 == 0) = NaN ;
                plot(lin, lin, '--')
                scatter(vv2, vv1, 5, color, 'filled', 'markeredgecolor', 'none') ;
                xlim([-1 1])
                ylim([-1 1])
                axis equal
                
                subplot(2, 3, 3)
                title('Relative motion') ;
                ylabel('V_v_,_r_e_l') ;
                xlabel('V_u_,_r_e_l') ;
                hold on
                scatter(squeeze(data.relvels2{dd}(ll, 1, tt, :)), squeeze(data.relvels2{dd}(ll, 2, tt, :)), 5, color, 'filled', 'markeredgecolor', 'none') ;
                xlim([-1 1])
                ylim([-1 1])
            end
        end
    end
    
    subplot(2, 3, 4)
    title('Relative Phase') ;
    ylabel('Count') ;
    xlabel('Phase') ;
    xlim([-pi pi]);
    hold on
    N21 = zeros(maxltp-1, 9) ;
    for tt = 1:maxltp-1
        phase = phases2(ll, tt, phases2(ll, tt, :) ~= 0) ;
        [~,edges] = histcounts(phase, linspace(-pi, pi, 10));
        N21(tt, :) = histcounts(phase ,edges); % Bin using the same edges
    end
    ctrs1 = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
    b1 = bar(ctrs1, N21, 'stacked', 'edgecolor', 'none') ;
    for kk = 1:size(N21,1)
        cmap = parula(maxltp) ;
        color = cmap(kk, :) ;
        b1(kk).FaceColor = color ;
    end
    
    subplot(2, 3, 5)
    title('Relative Magnitude |V_m_e_m|/|V_n_u_c|') ;
    ylabel('Count') ;
    xlabel('Relative Magnitude') ;
    relmags = zeros(maxltp-1, size(binvels2, 4)) ;
    N22 = zeros(maxltp-1, 9) ;
    for tt = 1:maxltp-1
        relmags(tt, :) = sqrt((squeeze(binvels2(ll, 1, tt, :))).^2 + (squeeze(binvels2(ll, 2, tt, :))).^2)./sqrt((squeeze(binvels2(ll, 3, tt, :))).^2 + (squeeze(binvels2(ll, 4, tt, :))).^2) ;
        hold on
        relmag = relmags(tt, relmags(tt,:) ~=0) ;
        [~,edges] = histcounts(relmag, linspace(0, 2, 10));
        N22(tt, :) = histcounts(relmag ,edges); % Bin using the same edges
    end
    ctrs2 = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
    b2 = bar(ctrs2, N22, 'stacked', 'FaceColor', 'flat', 'edgecolor', 'none') ;
    for kk = 1:size(N22, 1)
        cmap = parula(maxltp) ;
        color = cmap(kk, :) ;
        b2(kk).CData = color ;
    end
    
    subplot(2, 3, 6)
    title('Normed Dot Products') ;
    ylabel('Count') ;
    xlabel('Normed Dot Product Values');
    N23 = zeros(maxltp-1, 9) ;
    for tt = 1:maxltp-1
        normdot = normdots2(ll, tt, normdots2(ll, tt,:) ~=0) ;
        [~,edges] = histcounts(normdot, linspace(-1, 1, 10));
        N23(tt, :) = histcounts(normdot ,edges); % Bin using the same edges
    end
    ctrs3 = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
    b3 = bar(ctrs3, N23, 'stacked', 'FaceColor', 'flat', 'edgecolor', 'none') ;
    for kk = 1:maxltp-1
        cmap = parula(maxltp) ;
        color = cmap(kk, :) ;
        b3(kk).CData = color ;
    end
    
    % relvelucomp = squeeze(data.relvels{dd}(1,:,:)) ;
    % relvelvcomp = squeeze(data.relvels{dd}(2,:,:)) ;
    % data = horzcat(relvelucomp(:),relvelvcomp(:)) ;
    % [r_ellipse, X0, Y0] = error_ellipse(data) ;
    % subplot(2, 3, 3)
    % p = plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'color', 'r','linestyle', '-') ;
    % legend(p, '95% Confidence Interval') ;
    
    sgtitle(['Master Correlation Plot Lobe ' num2str(ll)]) ;
    saveas(gcf, fullfile(data.masterDir, ['master_corr_lobe' num2str(ll) '.png'])) ;
end