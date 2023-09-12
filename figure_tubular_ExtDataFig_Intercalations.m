%% Create Extended Data Figure with Convergent Extension 
% Decompose convergent extension into cell shape change and intercalations



%% Clear workspace ========================================================
% We start by clearing the memory and closing all figures
clear; close all; clc;
cd /mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/
zwidth =1 ;

dataDir = cd ;

%% ADD PATHS TO THIS ENVIRONMENT ==========================================
origpath = matlab.desktop.editor.getActiveFilename;
cd('/mnt/data/code/tubular/example')
addpath(genpath('../'))
addpath(genpath('../utility'))
addpath(genpath('/mnt/data/code/gptoolbox'))
addpath('../TexturePatch')
addpath('../DECLab')
addpath(fullfile('../utility','plotting'))
addpath(fullfile('../utility','plotting'))
% go back to the data
cd(dataDir)

disp('loading xp struct from disk')
load(fullfile(dataDir, 'xp.mat'), 'xp', 'opts')

%% TubULAR class instance
disp('defining TubULAR class instance (tubi= tubular instance)')
tubi = TubULAR(xp, opts) ;
disp('done defining TubULAR instance')


%% Now generate the figure panel c with tissue tectonics
options = struct() ;
options.timePoints = [96:10:206] ;
options.overwrite = false ;
tubi.generateCellSegmentationPathlines3D(options)

%% Get number of cells used to make this measurement

% number in the advected pathlines
cellVertexPathlineFn = fullfile(tubi.dir.segmentation, 'pathlines', ...
    sprintf('cellVertexPathlines_%06dt0.mat', t0)) ;
load(cellVertexPathlineFn, 'segVertexPathlines2D', ...
            'segVertexPathlines3D', 'cellIDs')

% number in the segmented timepoints



%% Below is the key code that is run by the above method
% Compare to true segmentation GLOBALLY -- nematic strength and direction 
timePoints = 96:10:206 ;
t0 = tubi.t0set() ;
figW = 9 ;
figH = 4.5 ;
useCorrected = true ;



timeList = {timePoints, timePoints(timePoints < t0 + 75 & timePoints > -30)} ;
timeStr = {'', '_tlimit'} ;
stdsteStr = {'_std', '_ste'};
for std_ste = 1:2
    for pp = 1:2
        close all 
        fig = figure('units', 'centimeters', 'position', [0,0,figH,figH]) ;
        time2do = timeList{pp} ;
        pInds = ismember(timePoints, time2do) ;

        imfn = fullfile(tubi.dir.segmentation, 'pathlines', ...
            ['cell_anisotropy_global_signed_COMPARE', timeStr{pp}]) ;
        % Collate results
        kk = 1;
        meanQAspectsTrue = zeros(length(time2do), 1) ;
        meanQAspectStdsTrue = zeros(length(time2do), 1) ;
        meanQAspectStesTrue = zeros(length(time2do), 1) ;
        meanQThetasTrue = meanQAspectsTrue ;
        for tp = time2do
            tubi.setTime(tp) ;
            if useCorrected     
                seg3d = tubi.getCurrentSegmentation3DCorrected() ; 
            else
                seg3d = tubi.getCurrentSegmentation3D() ; 
            end
            meanQAspectsTrue(kk) = seg3d.seg3d.statistics.meanQ.aspectWeighted.meanQAspect ;
            meanQAspectStdsTrue(kk) = seg3d.seg3d.statistics.meanQ.aspectWeighted.meanQAspectStd ;
            meanQAspectStesTrue(kk) = seg3d.seg3d.statistics.meanQ.aspectWeighted.meanQAspectSte ;
            meanQThetasTrue(kk) = seg3d.seg3d.statistics.meanQ.aspectWeighted.meanQTheta ;
            meanQThetaStdsTrue(kk) = seg3d.seg3d.statistics.meanQ.aspectWeighted.meanQThetaStd ;
            meanQThetaStesTrue(kk) = seg3d.seg3d.statistics.meanQ.aspectWeighted.meanQThetaSte ;
            kk = kk + 1;
        end


        c2tTrue = cos(2*meanQThetasTrue(pInds))  ;
        s2tTrue = sin(2*meanQThetasTrue(pInds))  ;
        trueAn = meanQAspectsTrue(pInds) - 1 ;
        trueAnStd = meanQAspectStdsTrue(pInds)  ;
        trueAnSte = meanQAspectStesTrue(pInds)  ;
        trueThetaStd = meanQThetaStdsTrue(pInds) ;
        trueThetaSte = meanQThetaStesTrue(pInds) ;
        trueLine = c2tTrue .* trueAn ;
        trueStds = sqrt( (c2tTrue .* trueAnStd).^2 + ...
            (trueAn .* s2tTrue * 2 .* trueThetaStd(:)).^2) ;
        trueStes = sqrt( (c2tTrue .* trueAnSte).^2 + ...
            (trueAn .* s2tTrue * 2 .* trueThetaSte(:)).^2) ;
        % ratios
        trueAR = meanQAspectsTrue(pInds) ;
        trueLineR = c2tTrue .* trueAR ;
        trueStdsR = sqrt( (c2tTrue .* trueAnStd).^2 + ...
            (trueAR .* s2tTrue * 2 .* trueThetaStd(:)).^2) ;
        trueStesR = sqrt( (c2tTrue .* trueAnSte).^2 + ...
            (trueAR .* s2tTrue * 2 .* trueThetaSte(:)).^2) ;

        c2t = cos(2*mQAtheta(pInds))  ;
        s2t = sin(2*mQAtheta(pInds))  ;
        An = squeeze(mQAar(pInds)) - 1 ;
        AnStd = squeeze(mQAarStd(pInds)) ;
        AnSte = squeeze(mQAarSte(pInds)) ;
        ThetaStd = squeeze(mQAthetaStd(pInds)) ;
        ThetaSte = squeeze(mQAthetaSte(pInds)) ;
        midline = c2t .* An ;
        stds = sqrt( (c2t .* AnStd).^2 + (An .* s2t *2 .* ThetaStd).^2) ;
        stes = sqrt( (c2t .* AnSte).^2 + (An .* s2t *2 .* ThetaSte).^2) ;
        % ratios
        AR = squeeze(mQAar(pInds)) ;
        midlineR = c2t .* AR ;
        stdsR = sqrt( (c2t .* AnStd).^2 + (AR .* s2t *2 .* ThetaStd).^2) ;
        stesR = sqrt( (c2t .* AnSte).^2 + (AR .* s2t *2 .* ThetaSte).^2) ;


        timestamps = timePoints - t0 ;
        if contains(tubi.timeUnits, 'min')
            timestamps = timestamps / 60 ;
            timeunits = 'hr';
        else
            timeunits = tubi.timeUnits ;
        end

        % Transformation (flip sign or subtract initial value)
        midline0 = midline ; 
        trueLine0 = trueLine ;
        midline = midline - midline(1) ;
        trueLine = trueLine - trueLine(1) ;


        if std_ste == 1
            hs = errorbar(midline, trueLine, trueStds, trueStds, stds, stds) ;
        else
            [midline0, inds] = sort(midline) ;
            trueLine0 = trueLine(inds) ;
            trueStds0 = trueStds(inds) ;
            trueStes0 = trueStes(inds) ;
            bounds = [trueLine0-abs(trueStds0), trueLine0+abs(trueStds0)] ;
            lowerB = bounds(:, 1)' ;
            upperB = bounds(:, 2)' ;
            fill([midline0, fliplr(midline0)], ...
                [lowerB, fliplr(upperB)], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs = errorbar(midline0, trueLine0, ...
                trueStes0, trueStes0, 'color', colors(1, :)) ;
        end
        hold on;
        ylims = [ylim() xlim()] ;

        % Mark zero line
        plot([0,2], [0,2], 'k--', 'HandleVisibility','off')
        % Labels
        % legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        ylabel('cell shear, $\Delta(a/b)$',   'interpreter', 'latex')
        xlabel('tissue shear, $\Delta(a/b)$',   'interpreter', 'latex')
        axis equal
        % ylim([-max(abs(ylims)), max(abs(ylims))])
        % xlim([-max(abs(ylims)), max(abs(ylims))])
        
        ylim([-Inf, max(abs(ylims))])
        xlim([-Inf, max(abs(ylims))])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '.pdf'])
        


        %% Global Difference over time
        clf
        x2 = [timestamps(pInds), fliplr(timestamps(pInds))] ;
        if std_ste == 1
            hs2 = errorbar(timestamps(pInds), trueLine(:),trueStds(:)) ;
            hold on;
            hs = errorbar(timestamps(pInds), midline(:), stds(:)) ;
        else
            fill(x2, [midline - stds, fliplr(midline + stds)], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            fill(x2, [trueLine - trueStds; flipud(trueLine + trueStds)], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs2 = errorbar(timestamps(pInds), trueLine(:),trueStes(:), 'color', colors(1, :)) ;
            hold on;
            hs = errorbar(timestamps(pInds), midline(:), stes(:), 'color', colors(2, :)) ;
        end
        % Labels
        % legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        xlabel(['time [' timeunits ']'],   'interpreter', 'latex')
        ylabel('tissue shear, cell shear',   'interpreter', 'latex')
        legend([hs, hs2], {'tissue shear', 'cell shape change'})
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff.pdf'])
        
        
        
        %% Global curves and difference over time -- shaded ste or shaded std
        clf
        x2 = [timestamps(pInds), fliplr(timestamps(pInds))] ;
        diffLine = midline(:) - trueLine(:) ;
        diffStds = sqrt(stds(:).^2 + trueStds(:).^2) ;
        diffStes = sqrt(stes(:).^2 + trueStes(:).^2) ;
        if std_ste == 1
            % hs2 = errorbar(timestamps(pInds), trueLine(:),trueStds(:)) ;
            % hold on;
            % hs = errorbar(timestamps(pInds), midline(:), stds(:)) ;
            hs = fill(x2, [midline - stds, fliplr(midline + stds)], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs2 = fill(x2, [trueLine - trueStds; flipud(trueLine + trueStds)], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [diffLine - diffStds; flipud(diffLine + diffStds)], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            errorbar(timestamps(pInds), trueLine(:),trueStes(:), 'color', colors(1, :)) ;
            hold on;
            errorbar(timestamps(pInds), midline(:), stes(:), 'color', colors(2, :)) ;
            errorbar(timestamps(pInds), diffLine(:), diffStes(:), 'color', colors(3, :)) ;
        else
            hs = fill(x2, [midline - stes, fliplr(midline + stes)], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs2 = fill(x2, [trueLine - trueStes; flipud(trueLine + trueStes)], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [diffLine - diffStes; flipud(diffLine + diffStes)], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
        end
        % Labels
        % legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        xlabel(['time [' timeunits ']'],   'interpreter', 'latex')
        ylabel('tissue shear, cell shear, intercalations',   'interpreter', 'latex')
         legend([hs, hs2, hs3], {'tissue shear', 'cell shape change', 'inferred intercalations'})
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2.pdf'])
        set(gca, 'YScale', 'log')
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_semilog.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_semilog.pdf'])
        
        
        %% Global curves and difference over time, relative to start of plot, moving mean -- shaded ste or shaded std
        clf
        x2 = [timestamps(pInds), fliplr(timestamps(pInds))] ;
        diffLine = (midline(:) - trueLine(:)) ;
        diffStds = sqrt(stds(:).^2 + trueStds(:).^2) ;
        diffStes = sqrt(stes(:).^2 + trueStes(:).^2) ;
        if std_ste == 1
            % hs2 = errorbar(timestamps(pInds), trueLine(:),trueStds(:)) ;
            % hold on;
            % hs = errorbar(timestamps(pInds), midline(:), stds(:)) ;
            hs = fill(x2, [movmean(midline(:) - stds(:), 3); ...
                flipud(movmean(midline(:) + stds(:), 3))], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs2 = fill(x2, [movmean(trueLine(:) - trueStds(:), 3); ...
                flipud(movmean(trueLine(:) + trueStds(:), 3))], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [movmean(diffLine - diffStds, 3); ...
                flipud(movmean(diffLine + diffStds, 3))], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            errorbar(timestamps(pInds), movmean(trueLine, 3),...
                movmean(trueStes(:), 3), 'color', colors(1, :)) ;
            hold on;
            errorbar(timestamps(pInds), movmean(midline, 3), ...
                movmean(stes(:), 3), 'color', colors(2, :)) ;
            errorbar(timestamps(pInds), movmean(diffLine(:), 3), ...
                movmean(diffStes(:), 3), 'color', colors(3, :)) ;
        else
            hs = fill(x2, [movmean(midline - stes, 3), ...
                fliplr(movmean(midline + stes, 3))], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs2 = fill(x2, [movmean(trueLine - trueStes, 3); ...
                flipud(movmean(trueLine + trueStes, 3))], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [movmean(diffLine - diffStes, 3); ...
                flipud(movmean(diffLine + diffStes, 3))], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
        end
        % Labels
        % legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        xlabel(['time [' timeunits ']'],   'interpreter', 'latex')
        ylabel('tissue shear, cell shear, intercalations',   'interpreter', 'latex')
        try
         legend([hs, hs2, hs3], {'tissue shear', 'cell shape change', 'inferred intercalations'})
        catch
         legend([hs, hs2(1), hs3], {'tissue shear', 'cell shape change', 'inferred intercalations'})
        end
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean.pdf'])
        set(gca, 'YScale', 'log')
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog.pdf'])
        xlim([-0.25, Inf])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog_tx.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog_tx.pdf'])
        save([imfn, stdsteStr{std_ste}, '_diff2_movmean_DATA.mat'], 'dataPlotted')
        
        %% Global curves and difference over time, relative to t0, with moving mean -- shaded ste or shaded std
        clf
        idx = find(timestamps < 0, 1, 'last') ;
        x2 = [timestamps(pInds), fliplr(timestamps(pInds))] ;
        midlineT0 = midline0(:) - midline0(idx) ;
        trueLineT0 = trueLine0(:)-trueLine0(idx) ;
        diffLine = (midlineT0 - trueLineT0) ;
        diffStds = sqrt(stds(:).^2 + trueStds(:).^2) ;
        diffStes = sqrt(stes(:).^2 + trueStes(:).^2) ;
        dataPlotted = struct('tissueShear', movmean(midlineT0(:), 3) , ...
            'cellShear', movmean(trueLineT0(:), 3), ...
            'intercalations', movmean(diffLine(:), 3), ...
            'std_tissueShear', movmean(stds(:), 3) , ...
            'std_cellShear', movmean(trueStds(:), 3), ...
            'std_intercalations', movmean(diffStds(:), 3), ...
            'ste_tissueShear', movmean(stes(:), 3) , ...
            'ste_cellShear', movmean(trueStes(:), 3), ...
            'ste_intercalations', movmean(diffStes(:), 3)) ;
        if std_ste == 1
            % hs2 = errorbar(timestamps(pInds), trueLine(:),trueStds(:)) ;
            % hold on;
            % hs = errorbar(timestamps(pInds), midline(:), stds(:)) ;
            hs2 = fill(x2, [movmean(midlineT0(:) - stds(:), 3); ...
                flipud(movmean(midlineT0(:) + stds(:), 3))], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs = fill(x2, [movmean(trueLineT0(:) - trueStds(:), 3); ...
                flipud(movmean(trueLineT0(:) + trueStds(:), 3))], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [movmean(diffLine - diffStds, 3); ...
                flipud(movmean(diffLine + diffStds, 3))], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            errorbar(timestamps(pInds), movmean(trueLineT0, 3),...
                movmean(trueStes(:), 3), 'color', colors(1, :)) ;
            hold on;
            errorbar(timestamps(pInds), movmean(midlineT0, 3), ...
                movmean(stes(:), 3), 'color', colors(2, :)) ;
            errorbar(timestamps(pInds), movmean(diffLine(:), 3), ...
                movmean(diffStes(:), 3), 'color', colors(3, :)) ;
        else
            hs2 = fill(x2, [movmean(midlineT0 - stes(:), 3), ...
                flipud(movmean(midlineT0 + stes(:), 3))], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs = fill(x2, [movmean(trueLineT0 - trueStes, 3); ...
                flipud(movmean(trueLineT0 + trueStes, 3))], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [movmean(diffLine - diffStes, 3); ...
                flipud(movmean(diffLine + diffStes, 3))], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
        end
        % Labels
        % legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        xlabel(['time [' timeunits ']'],   'interpreter', 'latex')
        ylabel('tissue shear, cell shear, intercalations',   'interpreter', 'latex')
        try
         legend([hs2, hs, hs3], {'tissue shear', 'cell shape change', 'inferred intercalations'})
        catch
         legend([hs2(1), hs, hs3], {'tissue shear', 'cell shape change', 'inferred intercalations'})
        end
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_t0.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_t0.pdf'])
        set(gca, 'YScale', 'log')
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog_t0.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog_t0.pdf'])
        xlim([0, Inf])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog_t0x.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog_t0x.pdf'])
        save([imfn, stdsteStr{std_ste}, '_diff2_movmean_t0_DATA.mat'], 'dataPlotted')
        
        
        %% Global curves and difference over time divided by first entry RATIO -- shaded ste or shaded std
        clf
        x2 = [timestamps(pInds), fliplr(timestamps(pInds))] ;
        diffLine1 = midline(:) - trueLine(:) + 1;
        diffStds = sqrt(stds(:).^2 + trueStds(:).^2) ;
        diffStes = sqrt(stes(:).^2 + trueStes(:).^2) ;
        if std_ste == 1
            % hs2 = errorbar(timestamps(pInds), trueLine(:),trueStds(:)) ;
            % hold on;
            % hs = errorbar(timestamps(pInds), midline(:), stds(:)) ;
            hs = fill(x2, [midlineR - stdsR, fliplr(midlineR + stdsR)], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs2 = fill(x2, [trueLineR - trueStdsR; flipud(trueLineR + trueStdsR)], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [diffLine1 - diffStds; flipud(diffLine1 + diffStds)], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            errorbar(timestamps(pInds), trueLineR(:),trueStesR(:), 'color', colors(1, :)) ;
            hold on;
            errorbar(timestamps(pInds), midlineR(:), stesR(:), 'color', colors(2, :)) ;
            errorbar(timestamps(pInds), diffLine1(:), diffStes(:), 'color', colors(3, :)) ;
        else
            hs = fill(x2, [midlineR - stesR, fliplr(midlineR + stesR)], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs2 = fill(x2, [trueLineR - trueStesR; flipud(trueLineR + trueStesR)], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [diffLine1 - diffStes; flipud(diffLine1 + diffStes)], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
        end
        % Labels
        % legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        xlabel(['time [' timeunits ']'],   'interpreter', 'latex')
        ylabel('tissue shear, cell shear, intercalations',   'interpreter', 'latex')
         legend([hs, hs2, hs3], {'tissue stretch ratio', 'cell stretch ratio', 'inferred intercalations'})
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diffRatio.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diffRatio.pdf'])
        set(gca, 'YScale', 'log')
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diffRatio_semilog.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diffRatio_semilog.pdf'])
        
    end
end


