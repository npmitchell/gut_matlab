% Analyze the GCAMP pulses for width of anterior fold
% NPMitchell 2021

datdir = '/Users/npmitchell/Desktop/gut/48YGAL4klarGCAMPGFP/analysis' ;
expts = {'20210206_e1', '20210206_e2' } ;

clipY0s = {[100, 350], [165, 415]} ;
clipYs = {{clipY0s{1}, clipY0s{1} + [50,-50], clipY0s{1} + [100, -100]}, ...
    {clipY0s{2}, clipY0s{2} + [50,-50], clipY0s{2} + [100, -100]}} ;
clipXs = {[1, 410], [1, 410]} ;
dts = [54/60, 54/60] ;
pix2um = [0.37895, 0.37895] ;
foldXs = [181, 160] ; % fold X position at frame foldT, in pixels
foldTs = [20, 34] ; % fold frame # at which position is measured
t0 = [1.5, 13.5] ;


xlimFix = [-50, 50] ;
ylimFix = [] ;


% Antp domain is about 40 microns, broadens on dorsal side.

% If drift in X, define it here as [frame#, offsetX] 
offsets = {[], [1, 27], } ;
for ee = 1:length(expts)
    fns1 = dir(fullfile(datdir, expts{ee}, 'images', '*c001.png')) ;
    fns2 = dir(fullfile(datdir, expts{ee}, 'images', '*c002.png')) ;
    fns3 = dir(fullfile(datdir, expts{ee}, 'images', '*c003.png')) ;
    
    % For each yrange
    for clipyPairIdx = 1:length(clipYs{ee})
        % this y range is defined
        clipY = clipYs{ee}{clipyPairIdx} ;
        clipX = clipXs{ee} ;
        foldX = foldXs(ee) ;
        foldT = foldTs(ee) ;
        
        % Look at differences along y for each timepoint
        ntps = length(fns1) ;
        cmap = viridis(ntps) ;
        timestamps = dts(ee) * (0:ntps-1) ;
        timestamps = timestamps - t0(ee) ;
        if clipyPairIdx == 1
            deltaX = zeros(ntps, 1) ;
        end
        clf
        for ii = 1:ntps
            im1 = imread(fullfile(fns1(ii).folder, fns1(ii).name)) ;
            im2 = imread(fullfile(fns3(ii).folder, fns2(ii).name)) ;
            im3 = imread(fullfile(fns3(ii).folder, fns3(ii).name)) ;

            c1 = im1(clipY(1):clipY(2), clipX(1):clipX(2)) ;
            c2 = im2(clipY(1):clipY(2), clipX(1):clipX(2)) ;
            c3 = im3(clipY(1):clipY(2), clipX(1):clipX(2)) ;
            % mean image
            curr = (c1 + c2 + c3) / 3 ; 

            if ii == 1 
                wX = clipX(2)-clipX(1) + 1 ;
                dd = ones(length(fns1), wX) ;
                ss = ones(length(fns1), wX) ;
                xshifted = zeros(length(fns1), wX) ;
                xx = pix2um(ee) * (0:wX-1) ;
                xfixed = xx ;
                dinterp = dd ; 
                sinterp = dd ;
            elseif clipyPairIdx == 1 
                % Find offset in x wrt previous image
                [deltaX(ii), ~] = ExtPhaseCorrelation(curr, prev) ;
            end    

            s1 = sum(c1, 1) ;
            s2 = sum(c2, 1) ;
            s3 = sum(c3, 1) ;

            d12 = abs(s1 - s2); 
            d13 = abs(s1 - s3) ;
            d23 = abs(s2 - s3) ;

            % Save transient signal in matrix
            ddii = d12 + d13 + d23 ;
            dd(ii, :) = ddii ;
            sii = s1 + s2 + s3 ;
            ss(ii, :) = sii ;
            
            % Save shifted ap positions in matrix
            xshiftii = xx+sum(deltaX(1:ii)) ;
            xshifted(ii, :) = xshiftii ;
            
            % Interpolate shifted results
            dinterp(ii, :) = interp1(xshiftii, ddii, xfixed) ;
            sinterp(ii, :) = interp1(xshiftii, sii, xfixed) ;
            
            % Factor in translation of the embryo over the timecourse
            plot(xx+sum(deltaX(1:ii)), dd(ii, :), 'color', cmap(ii, :)) ; hold on;

            % Save image for phase correlation with next image
            prev = curr ;
        end
        
        % Register the position so that 0um is anterior fold
        xshiftPix = xshifted(foldT, foldX) ;
        xshifted = xshifted - xshiftPix ;
        xfixed = xfixed - xshiftPix ;
        
        % First plot -- raw curves
        xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
        ylabel('GCAMP transient signal [a.u.]', 'interpreter', 'latex')    
        title('Transient GCAMP activity')
        cb = colorbar ;
        caxis([min(timestamps), max(timestamps)])
        ylabel(cb, 'time [min]')
        outfigfn = sprintf(fullfile(datdir, expts{ee}, '00_difference_Yrange%d.png'), clipyPairIdx) ;
        saveas(gcf, outfigfn) 
        

        % Kymograph -- pre shift
        clf
        imagesc(xx, timestamps, dd)
        xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
        ylabel('time [min]', 'interpreter', 'latex')
        title('Transient GCAMP activity')
        cb = colorbar ;
        ylabel(cb, '$ \delta I $ [a.u.]', 'interpreter', 'latex')
        outfigfn = sprintf(fullfile(datdir, expts{ee}, '01_difference_kymo_Yrange%d.png'), clipyPairIdx) ;
        saveas(gcf, outfigfn)
        
        % Shifted kymo
        clf
        for row = 1:ntps
            imagesc(xshifted(row, :), timestamps(row), dd(row, :));
            hold on;
        end
        xlabel('stabilized ap position [$\mu$m]', 'interpreter', 'latex')
        ylabel('time [min]', 'interpreter', 'latex')
        title('Transient GCAMP signal', ...
            'interpreter', 'latex')
        xlim(xlimFix)
        if ~isempty(ylimFix)
            ylim(ylimFix)
        else
            ylim([min(timestamps), max(timestamps)])
        end
        cb = colorbar ;
        ylabel(cb, '$ \delta I $ [a.u.]', 'interpreter', 'latex')
        outfigfn = sprintf(fullfile(datdir, expts{ee}, '02_difference_kymoShift_Yrange%d.png'), clipyPairIdx) ;
        saveas(gcf, outfigfn)
        
        % Plot interpolated results
        clf
        timeGrid = timestamps' .* ones(size(dd)) ;
        apGrid = xfixed .* ones(size(dd)) ;
        imagesc(xfixed, timestamps, dinterp) ;
        xlim(xlimFix)
        if ~isempty(ylimFix)
            ylim(ylimFix)
        else
            ylim([min(timestamps), max(timestamps)])
        end
        xlabel('stabilized, re-gridded ap position [$\mu$m]', 'interpreter', 'latex')
        ylabel('time [min]', 'interpreter', 'latex')
        title('Transient GCAMP signal [registered frame]', ...
            'interpreter', 'latex')
        cb = colorbar ;
        ylabel(cb, '$ \delta I $ [a.u.]', 'interpreter', 'latex')
        outfigfn = sprintf(fullfile(datdir, expts{ee}, ...
            '03_difference_kymoShiftInterp_Yrange%d.png'), clipyPairIdx) ;
        saveas(gcf, outfigfn)
        
        % Filter out pulses
        close all
        background = movmedian(sqrt(sinterp),  ceil(15 / dts(ee))) ;
        imagesc(xfixed, timestamps, background)
        xlabel('stabilized, re-gridded ap position [$\mu$m]', 'interpreter', 'latex')
        ylabel('time [min]', 'interpreter', 'latex')
        xlim(xlimFix)
        if ~isempty(ylimFix)
            ylim(ylimFix)
        else
            ylim([min(timestamps), max(timestamps)])
        end
        title('Square root of median background signal [registered frame]', ...
            'interpreter', 'latex')
        cb = colorbar() ;
        ylabel(cb, '$\sqrt{\langle I \rangle_{15\mathrm{ min}}}$ [a.u.]', ...
            'interpreter', 'latex')
        outfigfn = sprintf(fullfile(datdir, expts{ee}, ...
            '04_mean_filtered15min_Yrange%d.png'), clipyPairIdx) ;
        saveas(gcf, outfigfn)
        
        % Plot interpolated results - background
        close all
        imagesc(xfixed, timestamps, dinterp-background) ;
        xlim(xlimFix)
        if ~isempty(ylimFix)
            ylim(ylimFix)
        else
            ylim([min(timestamps), max(timestamps)])
        end
        xlabel('stabilized, re-gridded ap position [$\mu$m]', 'interpreter', 'latex')
        ylabel('time [min]', 'interpreter', 'latex')
        title('Transient GCAMP activity [registered frame]')
        cb = colorbar ;
        ylabel(cb, '$ \delta I $ [a.u.]', 'interpreter', 'latex')
        outfigfn = sprintf(fullfile(datdir, expts{ee}, ...
            '05_difference_kymoShiftInterpBgsub_Yrange%d.png'), clipyPairIdx) ;
        saveas(gcf, outfigfn)
        
        % Average over time
        dbgsug = dinterp-background ;
        
        % Get t=0 idx
        [~, tidx0] = min(abs(timestamps)) ;
        
        signal00 = mean(dinterp(1:tidx0, :))' ;
        signal05 = mean(dinterp(1:tidx0+5, :))' ;
        signal10 = mean(dinterp(1:tidx0+10, :))' ;
        signal15 = mean(dinterp(1:tidx0+15, :))' ;
        signal20 = mean(dinterp(1:tidx0+20, :))' ;
        signal25 = mean(dinterp(1:tidx0+25, :))' ;
        signal30 = mean(dinterp(1:tidx0+30, :))' ;
        activity00 = mean(dbgsug (1:tidx0, :))' ;
        activity05 = mean(dbgsug (1:tidx0+5, :))' ;
        activity10 = mean(dbgsug (1:tidx0+10, :))' ;
        activity15 = mean(dbgsug (1:tidx0+15, :))' ;
        activity20 = mean(dbgsug (1:tidx0+20, :))' ;
        activity25 = mean(dbgsug (1:tidx0+25, :))' ;
        activity30 = mean(dbgsug (1:tidx0+30, :))' ;
        legends = {'$\langle t<$0 minutes$\rangle$', ...
            '$\langle t<$5 minutes$\rangle$', ...
            '$\langle t<$10 minutes$\rangle$', ...
            '$\langle t<$15 minutes$\rangle$', ...
            '$\langle t<$20 minutes$\rangle$', ...
            '$\langle t<$25 minutes$\rangle$', ...
            '$\langle t<$30 minutes$\rangle$'} ;
        colors = flipud(viridis(7)) ;
        
        clf
        plot(xfixed, signal00, 'color', colors(1, :)); hold on;
        plot(xfixed, signal05, 'color', colors(2, :)) 
        plot(xfixed, signal10, 'color', colors(3, :)) 
        plot(xfixed, signal15, 'color', colors(4, :))
        plot(xfixed, signal20, 'color', colors(5, :))
        plot(xfixed, signal25, 'color', colors(6, :))
        plot(xfixed, signal30, 'color', colors(7, :))
        legend(legends, 'interpreter', 'latex')
        xlabel('stabilized ap position [$\mu$m]', 'interpreter', 'latex')
        ylabel('GCAMP transient signal [a.u.]', 'interpreter', 'latex')    
        title('Transient GCAMP signal', ...
            'interpreter', 'latex')
        xlim(xlimFix)
        outfigfn = sprintf(fullfile(datdir, expts{ee}, '06_difference_mean_Yrange%d.png'), clipyPairIdx) ;
        saveas(gcf, outfigfn) 
        
        clf
        plot(xfixed, activity00, 'color', colors(1, :)); hold on;
        plot(xfixed, activity05, 'color', colors(2, :))
        plot(xfixed, activity10, 'color', colors(3, :))
        plot(xfixed, activity15, 'color', colors(4, :))
        plot(xfixed, activity20, 'color', colors(5, :))
        plot(xfixed, activity25, 'color', colors(6, :))
        plot(xfixed, activity30, 'color', colors(7, :))
        legend(legends, 'interpreter', 'latex')
        xlabel('stabilized ap position [$\mu$m]', 'interpreter', 'latex')
        ylabel('GCAMP transient signal [a.u.]', 'interpreter', 'latex')    
        title('Transient GCAMP activity', ...
            'interpreter', 'latex')
        xlim(xlimFix)
        outfigfn = sprintf(fullfile(datdir, expts{ee}, '07_difference_meanBgSub_Yrange%d.png'), clipyPairIdx) ;
        saveas(gcf, outfigfn) 
        
        % Filter in space
        close all
        w1um = round(1 / pix2um(ee)) ;
        plot(xfixed, movmedian(activity00, w1um), 'color', colors(1, :)); 
        hold on;
        plot(xfixed, movmedian(activity05, w1um), 'color', colors(2, :)); 
        plot(xfixed, movmedian(activity10, w1um), 'color', colors(3, :)); 
        plot(xfixed, movmedian(activity15, w1um), 'color', colors(4, :)); 
        plot(xfixed, movmedian(activity20, w1um), 'color', colors(5, :)); 
        plot(xfixed, movmedian(activity25, w1um), 'color', colors(6, :)); 
        plot(xfixed, movmedian(activity30, w1um), 'color', colors(7, :)); 
        legend(legends, 'interpreter', 'latex')
        xlabel('stabilized ap position [$\mu$m]', 'interpreter', 'latex')
        ylabel('GCAMP transient signal [a.u.]', 'interpreter', 'latex')    
        title('Transient GCAMP activity, filtered in space 1 $\mu$m', ...
            'interpreter', 'latex')
        xlim(xlimFix)
        outfigfn = sprintf(fullfile(datdir, expts{ee}, '08_difference_meanBgSubSm_Yrange%d.png'), clipyPairIdx) ;
        saveas(gcf, outfigfn) 
        
        
        % Save result
        outfn = sprintf([expts{ee} '_results_Yrange%d.mat'], clipyPairIdx) ;
        save(fullfile(datdir, outfn), 'xx', 'xshifted', 'dd', 'deltaX', ...
            'xfixed', 'dinterp', 'timeGrid', 'apGrid', ...
            'signal00', 'signal05', 'signal10', ...
            'signal15', 'signal20', 'signal25', 'signal30', ...
            'activity00', 'activity05', 'activity10', ...
            'activity15', 'activity20', 'activity25', 'activity30')
    end
    
end


%% Average all experiments together
clipyPairIdx = 1; 
colors = define_colors ;
close all
for ee = 1:length(expts)
    outfn = sprintf([expts{ee} '_results_Yrange%d.mat'], clipyPairIdx) ;
    resfn = fullfile(datdir, outfn) ;
    load(resfn, 'xfixed', 'activity30')
    
    w1um = round(1 / pix2um(ee)) ;
    keep = xfixed > xlimFix(1) & xfixed < xlimFix(2) ;
    activity = movmedian(activity30, w1um) ;
    idx2keep = find(keep) ;
    normIdx = [idx2keep(1): idx2keep(1)+round(5/pix2um(ee))] ;
    normIdx = [normIdx, idx2keep(end)-round(5/pix2um(ee)):idx2keep] ;
    normVal = mean(activity(normIdx)) ;
    anorm = activity / normVal ; 
    
    plot(xfixed(keep), anorm(keep), 'color', colors(ee, :)); 
    hold on;
end
xlabel('ap position from anterior fold [$\mu$m]', 'interpreter', 'latex')
ylabel('normalized transient GCAMP activity [a.u.]', ...
    'interpreter', 'latex')
title('Transient calcium activity', 'interpreter', 'latex')
resfn = fullfile(datdir, 'gcamp_activity_results') ;
saveas(gcf, [resfn '.png'])
saveas(gcf, [resfn '.pdf'])
