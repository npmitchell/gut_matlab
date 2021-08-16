%% Compare GCaMP signals in Heterozygotes to mutants 
% NPMitchell 2021
clear all

%% Global plotting options
preview = false ;
% Fixed xlimits in microns
ylimFix = [] ;
ylimFixAll = [-15, 90] ;
cLimits = {[0, 5.5], [0, 2.75]} ;
xlabels = {'ap position from proventriculus [$\mu$m]', ...
    'ap position from anterior fold [$\mu$m]'} ;
pausetimeCurve = 0.1 ;
pausetimeKymo = 0.25 ;
shoulder_um = 5 ;



%% Options
datdir = '/mnt/data/confocal_data/gut/48YGAL4GCaMP6sIIAntpNSRC3/analysis_AntpNSRC3_48YGAL4GCaMP6sII/' ;

resDirControlFn = 'anteriorMutantHeterozygotes' ;
resDirFn = 'anteriorMutantResults' ;
compareDirFn = 'anteriorMutantControlComparison' ;
if ~exist(fullfile(datdir, compareDirFn), 'dir')
    mkdir(fullfile(datdir, compareDirFn)) ;
end

% tOffset is the offset time between morphological stage of posterior
% Folding and anterior fold onset in heterozygotes (controls)
load(fullfile(datdir, 'anteriorMutantHeterozygotes', 'spacetime_offsets_controls.mat'), ...
    'xoffs', 'xOffset', 'xOffset_std', 'xOffset_unc', ...
    'toffs', 'tOffset', 'tOffset_std', 'tOffset_unc') ;
controlRes = load(fullfile(datdir, resDirControlFn, 'AntpNSRC3HeterozSettings.mat'), ...
    'pcPower', 'dz', 'anteriorFoldXs', 'anteriorFoldTs', ...
    'posteriorTs', 'dates',...
    'expts', 'dts', 'xfixed_antFoldRef', 'xfixed_nonFoldRef', 'pix2um') ;
controlDates = controlRes.dates ; 
controlExpts = controlRes.expts ;
controlDts = controlRes.dts ;
controlAnteriorFoldTs = controlRes.anteriorFoldTs ;
controlPosteriorTs = controlRes.posteriorTs ;
controlPix2um = controlRes.pix2um ;
nControlExpts = length(controlExpts) ;
controlXFixed_antFoldRef = controlRes.xfixed_antFoldRef ;
controlXFixed_nonFoldRef = controlRes.xfixed_nonFoldRef ;
mutRes = load(fullfile(datdir, resDirFn, 'AntpNSRC3HeterozSettings.mat'), ...
    'pcPower', 'dz', 'dates', 'dts', 'expts', 'xfixed', 'pix2um') ;
mutDates = mutRes.dates ; 
expts = mutRes.expts ;
dts = mutRes.dts ;
pix2um = mutRes.pix2um ;

% First get mean shoulders from control as bg
% Get mean peak for scale

% NOTE: There are two clusters of bgVals: 1-6, and 7 onward.
% These correspond to two different imaging conditions used.
% 202107172106 onward used a different condition than before!
expts2include = 1:length(expts) ;
expts2includeControl = 1:length(controlExpts);
for refID = [2,1]
    if refID == 1
        % onset of POSTERIOR folding in minutes
        t0 = controlDts .* controlPosteriorTs ;
        refStr = '_nonFoldRef' ;
        xfixedGlobal = controlXFixed_nonFoldRef ;
        fixXticks = [0, 25, 50, 75, 100] ;
        xlimFix = [0, 80] ;
        xOff = 0 ;
    else
        % onset of ANTERIOR folding in minutes
        t0 = controlDts .* controlAnteriorFoldTs ;
        refStr = '_antFoldRef' ;
        xfixedGlobal = controlXFixed_antFoldRef ;
        fixXticks = [-40, -20, 0, 20, 40] ;
        xlimFix = [-40,40] ;
        xOff = xOffset ;
    end

    for clipyPairIdx = [2,1]
        colors = define_colors(nControlExpts) ;
        fixTimeStamps = -20:1.5:90.5 ;  % minutes
        kymoM = zeros(length(xfixedGlobal), length(fixTimeStamps)) ;
        ksz = size(kymoM) ;
        nsamples = kymoM ;
        kymoMControl = kymoM ;
        nsamplesControl = nsamples ;

        close all
        
        % For different averaging
        avgMin = [15, 20, 25, 30, 45, 60] ;
        allres = [] ;
        for avgID = 1:length(avgMin)
            edmy = 1 ;
            disp(['Performing ' num2str(avgMin(avgID)) ' min average'])
            close all
            
            if refID == 1
                outdir = fullfile(datdir, compareDirFn, ...
                    'nonFoldRef', sprintf('YClip%0d', clipyPairIdx)) ;
                % , ...  sprintf('avg%0dMin', avgMin(avgID))) ;
            else
                outdir = fullfile(datdir, compareDirFn, ...
                    'antFoldRef', sprintf('YClip%0d', clipyPairIdx)) ;
                % , ...  sprintf('avg%0dMin', avgMin(avgID))) ;
            end
            if ~exist(outdir, 'dir')
                mkdir(outdir)
            end
            
            % Get NORMALIZATION as function of bgVal
            bgVals = zeros(1, length(expts2includeControl)) ;
            clf
            for ee = expts2includeControl
                disp(['loading ee = ' num2str(ee)])
                outfn = sprintf([controlExpts{ee} '_results_Yrange%d' refStr '.mat'], clipyPairIdx) ;
                resfn = fullfile(datdir, resDirControlFn, outfn) ;
                load(resfn,  'activity15', ...
                    'activity20', 'activity25', 'activity30', 'activity45', 'activity60',...
                    'timeGrid', 'dinterp')
                timestamps = unique(timeGrid) ;
                acts = {activity15, activity20, activity25, activity30, ...
                    activity45, activity60} ;

                w1um = round(1 / controlPix2um(ee)) ;
                activity = movmedian(acts{avgID}, w1um) ;
                padd = round(5/controlPix2um(ee)) ;

                % no need for this anymore -- all x are xfixed already
                % keep = xfixed > xlimFix(1) & xfixed < xlimFix(2) ;
                % idx2keep = find(keep) ;
                % normIdx = (idx2keep(1) : idx2keep(1)+padd) ;
                % normIdx = [normIdx, idx2keep(end)-padd:idx2keep] ;

                normVal = nanmedian(activity) ;
                bgVal = nanmean(mink([activity(1:padd); activity(end-padd:end)], 10)) ;
                % assert(~isnan(bgVal))
                maxVal = maxk(activity(1+padd:end-padd), 10) ;
                maxVal = mean(maxVal) ;
                anorm = (activity - bgVal) / (maxVal - bgVal); 
                bgVals(ee) = bgVal ; 
                maxVals(ee) = maxVal ; 
                
                if preview
                    figure(2)
                    plot(activity); hold on;
                    pause(pausetimeCurve)
                end
            end
            
            % Fit a model for the normalization
            pp = polyfit(bgVals(~isnan(bgVals)), maxVals(~isnan(bgVals)), 1) ;
            fig = figure( 'units','centimeters','position',[0,0,9,9]);
            clf
            scatter(bgVals, maxVals, 25, 'filled'); hold on;
            plot([min(bgVals), max(bgVals)], polyval(pp, [min(bgVals) max(bgVals)]), 'k--')
            xlabel('background $\delta I$ [a.u.]', ...
                'interpreter', 'latex')
            ylabel('brightest DV-averaged signal, $\delta I_{\mathrm{max}}$  [a.u.]', ...
                'interpreter', 'latex')
            title('Normalization model for $\delta I\rightarrow (\delta I-\delta I_{\mathrm{bg}})/(\delta I_{\mathrm{max}} - \delta I_{\mathrm{bg}})$', ...
                'interpreter', 'latex')
            ylims = ylim ;
            ylim([0, ylims(2)])
            saveas(gcf, fullfile(outdir, sprintf('noise_model_avgMin%d.pdf', avgMin(avgID)))) ;
            
            % GET Control after normalization
            statAllControl = zeros(length(xfixedGlobal), length(expts2includeControl)) ;
            clf
            edmy = 1 ;
            areaCurvControl = zeros(1, length(expts2includeControl)) ;
            midValsControl = zeros(1, length(expts2includeControl)) ;
            for ec = expts2includeControl
                disp(['loading ee = ' num2str(ec)])
                outfn = sprintf([controlExpts{ec} '_results_Yrange%d' refStr '.mat'], clipyPairIdx) ;
                resfn = fullfile(datdir, resDirControlFn, outfn) ;
                load(resfn, 'xfixed', 'activity15', ...
                    'activity20', 'activity25', 'activity30', ...
                    'activity45', 'activity60', ...
                    'timeGrid', 'dinterp')
                timestamps = unique(timeGrid) ;
                acts = {activity15, activity20, activity25, ...
                    activity30, activity45, activity60} ;

                w1um = round(1 / controlPix2um(ec)) ;
                activity = movmedian(acts{avgID}, w1um) ;
                padd = round(5/ controlPix2um(ec)) ;

                % no need for this anymore -- all x are xfixed already
                % keep = xfixed > xlimFix(1) & xfixed < xlimFix(2) ;
                % idx2keep = find(keep) ;
                % normIdx = (idx2keep(1) : idx2keep(1)+padd) ;
                % normIdx = [normIdx, idx2keep(end)-padd:idx2keep] ;

                normVal = nanmedian(activity) ;
                bgVal = mean(mink([activity(1:padd); activity(end-padd:end)], 10)) ;
                % maxVal = maxk(activity(1+padd:end-padd), 10) ;
                maxVal = polyval(pp, bgVal) ;
                
                anorm = (activity - bgVal) / (maxVal - bgVal); 
                
                % Look at center for t-test
                keep = (abs(xfixedGlobal)< shoulder_um) ;
                midValsControl(ec) = mean(anorm(keep)) ;
                
                % Measure area under the curve
                dxs = diff(xfixed) ;
                dx = dxs(1) ;
                assert(all(abs(dxs - dx) < 1e-5 * dx))
                areaCurvControl(ec) = sum(anorm .* dx) ;
                
                % Check for consistent sizing across stored samples
                assert(all(size(kymoMControl) == ksz))
                assert(all(size(xfixed) == size(xfixedGlobal)))

                % AVERAGE KYMO
                if avgID == 1
                    dnorm = (dinterp - normVal) / (maxVal - normVal) ;
                    intrp = griddedInterpolant({timestamps, xfixed'}, dnorm, 'nearest', 'none') ;
                    [tt, xx] = ndgrid(fixTimeStamps, xfixedGlobal) ;
                    newKymo = intrp(tt, xx) ;
                    nsamplesControl = nsamplesControl + ~isnan(newKymo)' ;
                    newKymo(isnan(newKymo)) = 0 ;
                    kymoMControl = kymoMControl + newKymo' ;
                    
                    if preview
                        figure(3)
                        clf
                        subplot(1, 2, 1)
                        imagesc(kymoMControl')
                        subplot(1, 2, 2)
                        imagesc(nsamplesControl')
                        pause(pausetimeKymo)
                    end
                end
            
                clf
                statAllControl(:, edmy) = anorm ;
                
                if preview
                    figure(2)
                    subplot(1, 2, 1)
                    hold on;
                    plot(xfixed, anorm, 'color', colors(edmy, :)); 
                    subplot(1, 2, 2)
                    hold on;
                    plot(xfixed, (activity - bgVal)' ./ ( maxk(activity(1+padd:end-padd), 10) - bgVal), ...
                        'color', colors(edmy, :)); 
                    pause(pausetimeCurve)
                end

                if ee == expts2include(1)
                    xall = xfixed ;
                end

                edmy = edmy + 1;
            end
            
            %% Now load MUTANT, shifted in space and time, apply normalization
            statAll = zeros(length(xfixed), length(expts2include)) ;
            midValsMutant = zeros(1, length(expts2include)) ;
            clf
            edmy = 1; 
            areaCurv = zeros(1, length(expts2include)) ;
            for ee = expts2include
                disp(['loading ee = ' num2str(ee)])
                outfn = sprintf([expts{ee} '_results_Yrange%d.mat'], clipyPairIdx) ;
                resfn = fullfile(datdir, resDirFn, outfn) ;
                load(resfn, 'xfixed', 'timeGrid', 'dinterp')
                timestamps = unique(timeGrid) ;
                
                % RECOMPUTE Activities with offset time
                [~, tidx0] = min(abs(timestamps+tOffset)) ;
                end5 = find(timestamps+tOffset - 5 > 0, 1);
                end10 = find(timestamps+tOffset - 10 > 0, 1);
                end15 = find(timestamps+tOffset - 15 > 0, 1);
                end20 = find(timestamps+tOffset - 20 > 0, 1);
                end25 = find(timestamps+tOffset - 25 > 0, 1);
                end30 = find(timestamps+tOffset - 30 > 0, 1);
                end45 = find(timestamps+tOffset - 45 > 0, 1);
                end60 = find(timestamps+tOffset - 60 > 0, 1);
                activity00 = mean(dinterp(1:tidx0, :))' ;
                activity05 = mean(dinterp(tidx0:end5, :))' ;
                activity10 = mean(dinterp(tidx0:end10, :))' ;
                activity15 = mean(dinterp(tidx0:end15, :))' ;
                activity20 = mean(dinterp(tidx0:end20, :))' ;
                activity25 = mean(dinterp(tidx0:end25, :))' ;
                activity30 = mean(dinterp(tidx0:end30, :))' ;
                activity45 = mean(dinterp(tidx0:end45, :))' ;
                activity60 = mean(dinterp(tidx0:end60, :))' ;                
                % Cat activities into cell
                acts = {activity15, activity20, activity25, activity30, ...
                    activity45, activity60} ;

                w1um = round(1 / pix2um(ee)) ;
                activity = movmedian(acts{avgID}, w1um) ;
                padd = round(5/pix2um(ee)) ;

                % no need for this anymore -- all x are xfixed already
                % keep = xfixed > xlimFix(1) & xfixed < xlimFix(2) ;
                % idx2keep = find(keep) ;
                % normIdx = (idx2keep(1) : idx2keep(1)+padd) ;
                % normIdx = [normIdx, idx2keep(end)-padd:idx2keep] ;

                normVal = nanmedian(activity) ;
                bgVal = mean(mink([activity(1:padd); activity(end-padd:end)], 10)) ;
                % bgVal = mean(bgVals) ;
                % maxVal = maxk(activity(1+padd:end-padd), 10) ;
                % maxVal = mean(maxVal) ;
                maxVal = polyval(pp, bgVal) ;
                anorm = (activity - bgVal) / (maxVal - bgVal); 
                
                % Look at center for t-test                
                keep = (abs(xfixedGlobal-xOff)< shoulder_um) ;
                midValsMutant(ee) = mean(anorm(keep)) ;
                
                statAll(:, edmy) = anorm ;
                
                if preview || true
                    figure(2) ;
                    hold on;
                    plot(xfixed, anorm, 'color', colors(edmy, :)); 
                end
                
                % measure area under the curve
                dxs = diff(xfixed) ;
                dx = dxs(1) ;
                assert(all(abs(dxs - dx) < 1e-5 * dx))
                areaCurv(ee) = sum(anorm .* dx) ;
                  
                if ee == expts2include(1)
                    xall = xfixed ;
                end

                % AVERAGE KYMO
                if avgID == 1
                    dnorm = (dinterp - normVal) / (maxVal - normVal) ;
                    intrp = griddedInterpolant({timestamps-tOffset, (xfixed-xOff)'}, dnorm, 'nearest', 'none') ;
                    [tt, xx] = ndgrid(fixTimeStamps, xfixedGlobal) ;
                    newKymo = intrp(tt, xx) ;
                    nsamples = nsamples + ~isnan(newKymo)' ;
                    newKymo(isnan(newKymo)) = 0 ;
                    kymoM = kymoM + newKymo' ;
                    assert(all(size(kymoM) == ksz))
                    assert(all(size(xfixed) == size(xfixedGlobal)))
                    
                    if preview
                        figure(3) ;
                        subplot(1, 2, 1)
                        imagesc(xfixedGlobal, fixTimeStamps, kymoM');
                        subplot(1, 2, 2) ;
                        imagesc(nsamples') ;
                        pause(pausetimeKymo) ;
                    end
                end
                
                edmy = edmy + 1;
            end
            
            %% T-test 
            % tobs = (x-mu) / (std / sqrt(n)) 
            % z score = (xmean1 - xmean2) / sqrt(std1^2 + std2^2)
            % assert(sum(~isnan(areaCurvControl)) == nControlExpts)
            
            num = mean(areaCurv) - nanmean(areaCurvControl) ;
            denom = sqrt(nanvar(areaCurvControl) / sum(~isnan(areaCurvControl)) + ...
                nanvar(areaCurv) / sum(~isnan(areaCurv))) ;
            zscore = num / denom ;
            pvalue = normcdf(zscore);
            
            
            %% ALternative T-test
            num = mean(midValsMutant) - nanmean(midValsControl) ;
            denom = sqrt(nanvar(midValsControl) / sum(~isnan(midValsControl)) + ...
                nanvar(midValsMutant) / sum(~isnan(midValsMutant))) ;
            zscoreMidVals = num / denom ;
            pvalueMidVals = normcdf(zscoreMidVals);
            
            %% Plot stats to compare each
            figure(2)
            close
            fig = figure('units','centimeters','position',[0,0,6,15]);

            % kymo panel CONTROL
            subplot(3, 1, 1)
            kymoMControl(~isfinite(kymoMControl)) = NaN;
            kymoMeanControl = (kymoMControl - min(kymoMControl(:))) ./ nsamplesControl  ;
            imagesc(xfixedGlobal, fixTimeStamps, kymoMeanControl' ); 
            caxis(cLimits{refID})
            % xlabel('ap position from anterior fold [$\mu$m]', 'interpreter', 'latex')
            ylabel('time from fold onset [min]', 'interpreter', 'latex')

            hold on;
            redcolor = colors(2, :) ;
            plot(xlimFix, 0*xlimFix + 5, '--', 'color', redcolor)
            plot(xlimFix, 0*xlimFix - 5, '--', 'color', redcolor)
            plot(xlimFix, 0*xlimFix, '-', 'color', redcolor)
            xlim(xlimFix)
            xticks(fixXticks)
            %set(gca, 'xticklabels', [])
            ylim(ylimFixAll)
            colormap(viridis)
            
            % kymo panel MUTANT
            subplot(3, 1, 2)
            kymoM(~isfinite(kymoM)) = NaN;
            kymoMean = (kymoM - min(kymoM(:))) ./ nsamples  ;
            imagesc(xfixedGlobal, fixTimeStamps, kymoMean' ); 
            caxis(cLimits{refID})
            % xlabel('ap position from anterior fold [$\mu$m]', 'interpreter', 'latex')
            ylabel('morphological time [min]', 'interpreter', 'latex')

            hold on;
            redcolor = colors(2, :) ;
            plot(xlimFix, 0*xlimFix + 5, '--', 'color', redcolor)
            plot(xlimFix, 0*xlimFix - 5, '--', 'color', redcolor)
            plot(xlimFix, 0*xlimFix, '-', 'color', redcolor)
            xlim(xlimFix)
            xticks(fixXticks)
            %set(gca, 'xticklabels', [])
            ylim(ylimFixAll)
            colormap(viridis)

            % curve panel -- Control
            subplot(3, 1, 3)
            cla
            avgact = mean(statAllControl, 2, 'omitnan') ;
            stdact = std(statAllControl, 0, 2, 'omitnan') ;
            lineProps = {'-','color', colors(1, :)} ;
            finalNorm = max(avgact(:)) ;
            factor = 1.0 / finalNorm ;
            avgact = movmean(avgact, 3) ;
            stdact = movmedian(stdact, 15) ;
            h1=shadedErrorBar(xfixedGlobal, avgact*factor, stdact*factor, 'lineProps', lineProps) ;
            
            % curve panel -- Mutant
            avgact = mean(statAll, 2, 'omitnan') ;
            stdact = std(statAll, 0, 2, 'omitnan') ;
            lineProps = {'-','color', colors(2, :)} ;
            factor = 1.0 / finalNorm ;
            avgact = movmean(avgact, 3) ;
            stdact = movmedian(stdact, 15) ;
            h1=shadedErrorBar(xfixedGlobal, avgact*factor, stdact*factor, 'lineProps', lineProps) ;
            % ylim([min(avgact*factor - stdact*factor), max(avgact*factor + stdact*factor)])
            ylim([-0.25, 1.5])
            xlim(xlimFix)
            xticks(fixXticks)
            xlabel(xlabels{refID}, 'interpreter', 'latex')
            ylabel('GCaMP activity [a.u.]', ...
                'interpreter', 'latex')
            title([sprintf('%d min average:', avgMin(avgID)), ...
                ' $p=$' num2str(pvalue, '%.3e')], ...
                'interpreter', 'latex')
            
            resfn = fullfile(outdir, 'gcamp_mean_results') ;
            saveas(gcf, [resfn sprintf(['_clipY%d_avg%d' refStr '.png'], clipyPairIdx, avgMin(avgID))])
            saveas(gcf, [resfn sprintf(['_clipY%d_avg%d' refStr '.pdf'], clipyPairIdx, avgMin(avgID))])

            
            % Alternative t-test plot
            title([sprintf('%d min average:', avgMin(avgID)), ...
               ' $p_{|x|<', sprintf('%01d', shoulder_um), ' \mu \mathrm{m}} =$', ...
                num2str(pvalueMidVals, '%.3e')], ...
                'interpreter', 'latex')
            resfn = fullfile(outdir, 'gcamp_mean_results_midValPValue') ;
            saveas(gcf, [resfn sprintf(['_clipY%d_avg%d' refStr '.png'], clipyPairIdx, avgMin(avgID))])
            saveas(gcf, [resfn sprintf(['_clipY%d_avg%d' refStr '.pdf'], clipyPairIdx, avgMin(avgID))])

            
        end
    end
   
    
    disp('done with this refID')
end

