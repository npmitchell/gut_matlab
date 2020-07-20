function plotMetricKinematics(QS, options)
% plotMetricKinematics(QS, options)
%   Plot the metric Kinematics as kymographs and correlation plots
% 
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   plot_kymographs : bool
%   plot_kymographs_cumsum : bool
%   plot_correlations : bool
%   plot_gdot_correlations : bool
%   plot_gdot_decomp : bool
% 
% NPMitchell 2020

plot_kymographs = true ;
plot_kymographs_cumsum = true ;
plot_correlations = true ;
plot_gdot_correlations = false ;
plot_gdot_decomp = true ;


if isfield(options, 'plot_kymographs')
    plot_kymographs = options.plot_kymographs ;
end
if isfield(options, 'plot_kymographs_cumsum')
    plot_kymographs_cumsum = options.plot_kymographs_cumsum ;
end
if isfield(options, 'plot_correlations')
    plot_correlations = options.plot_correlations ;
end
if isfield(options, 'plot_gdot_correlations')
    plot_gdot_correlations = options.plot_gdot_correlations ;
end
if isfield(options, 'plot_gdot_decomp')
    plot_gdot_decomp = options.plot_gdot_decomp ;
end

%% Store kymograph data in cell array
HHsK = {HH_apM, HH_lM, HH_rM, HH_dM, HH_vM} ;
gdotsK = {gdot_apM, gdot_lM, gdot_rM, gdot_dM, gdot_vM} ;
divvsK = {divv_apM, divv_lM, divv_rM, divv_dM, divv_vM} ;
velnsK = {veln_apM, veln_lM, veln_rM, veln_dM, veln_vM} ;
H2vnsK = {H2vn_apM, H2vn_lM, H2vn_rM, H2vn_dM, H2vn_vM} ;

%% Now plot different measured quantities as kymographs
if plot_kymographs
    % Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
    dvDir = fullfile(mKDir, 'avgDV') ;
    lDir = fullfile(mKDir, 'avgLeft') ;
    rDir = fullfile(mKDir, 'avgRight') ;
    dDir = fullfile(mKDir, 'avgDorsal') ;
    vDir = fullfile(mKDir, 'avgVentral') ;
    outdirs = {dvDir, lDir, rDir, dDir, vDir} ;
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    for qq = 1:length(outdirs)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end

        % Unpack what to plot (averaged kymographs, vary averaging region)
        HHK = HHsK{qq} ;
        gdotK = gdotsK{qq} ;
        divvK = divvsK{qq} ;
        velnK = velnsK{qq} ;
        H2vnK = H2vnsK{qq} ;
        m2plot = {gdotK, HHK, divvK, velnK, H2vnK} ;
        titles = {'$\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$',...
            'mean curvature, $H$', ...
            'divergence of flow, $\nabla \cdot \mathbf{v}$', ...
            'normal velocity, $v_n$', ...
            'normal motion, $v_n 2 H$'} ;
        labels = {['$\textrm{Tr}[g^{-1}\dot{g}]$ ' unitstr], ...
            ['mean curvature, $H$ ' Hunitstr], ...
            ['$\nabla \cdot \mathbf{v}$ ' unitstr], ...
            ['normal velocity, $v_n$ ' vunitstr] , ...
            ['normal motion, $v_n 2 H $ ' unitstr]} ;
        names = {'gdot', 'HH', 'divv', 'veln', 'H2vn'} ;
        climits = [climit, climit_H, climit, climit_veln, climit_err] ;

        %% Plot gdot/HH/divv/veln/H2vn DV-averaged kymograph
        for pp = 1:length(m2plot)
            
            % Check if images already exist on disk
            fn = fullfile(odir, [ names{pp} '.png']) ;
            fn_zoom = fullfile(odir, [names{pp} '_zoom_early.png']) ;
            
            if ~exist(fn, 'file') || ~exist(fn_zoom, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                colormap bwr
                imagesc((1:nU)/nU, tps, m2plot{pp})
                caxis([-climits(pp), climits(pp)])
                % Add folds to plot
                hold on;
                fons1 = max(1, fons(1)) ;
                fons2 = max(1, fons(2)) ;
                fons3 = max(1, fons(3)) ;
                plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
                plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
                plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

                % title and save
                title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeunits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
                fn = fullfile(odir, [ names{pp} '.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   

                % Zoom in on small values
                caxis([-climits(pp)/3, climits(pp)/3])
                fn = fullfile(odir, [names{pp} '_zoom.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
                % Zoom in on early times
                ylim([min(tps), max(fons) + 10])
                caxis([-climits(pp)/3, climits(pp)/3])
                fn = fullfile(odir, [names{pp} '_zoom_early.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
            end
        end
    end
end

%% kymographs of cumulative sums
if plot_kymographs_cumsum
    % Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
    dvDir = fullfile(mKDir, 'avgDV') ;
    lDir = fullfile(mKDir, 'avgLeft') ;
    rDir = fullfile(mKDir, 'avgRight') ;
    dDir = fullfile(mKDir, 'avgDorsal') ;
    vDir = fullfile(mKDir, 'avgVentral') ;
    outdirs = {dvDir, lDir, rDir, dDir, vDir} ;
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    for qq = 1:length(outdirs)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end

        % Unpack what to plot (averaged kymographs, vary averaging region)
        HHK = HHsK{qq} ;
        gdotK = cumsum(gdotsK{qq}, 1) ;
        divvK = cumsum(divvsK{qq}, 1) ;
        velnK = cumsum(velnsK{qq}, 1) ;
        H2vnK = cumsum(H2vnsK{qq}, 1) ;
        m2plot = {gdotK, divvK, velnK, H2vnK} ;
        titles = {'$\int_0^t\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$',...
            'divergence of flow, $\nabla \cdot \mathbf{v}$', ...
            'normal velocity, $\int_0^tv_n$', ...
            'normal motion, $\int_0^t v_n 2 H$'} ;
        labels = {['$\int_0^t\textrm{Tr}[g^{-1}\dot{g}]$ ' unitstr], ...
            ['$\int_0^t \nabla \cdot \mathbf{v}$ ' unitstr], ...
            ['normal velocity, $\int_0^tv_n$ ' vunitstr] , ...
            ['normal motion, $\int_0^t v_n 2 H $ ' unitstr]} ;
        names = {'Igdot', 'Idivv', 'Iveln', 'IH2vn'} ;
        climits = [climit, climit, climit_veln, climit] ;
        climits = climits * 3; 
        
        %% Plot gdot/HH/divv/veln/H2vn DV-averaged kymograph
        for pp = 1:length(m2plot)
            % Check if images already exist on disk
            fn = fullfile(odir, [ names{pp} '.png']) ;
            fn_zoom = fullfile(odir, [names{pp} '_zoom_early.png']) ;
            if ~exist(fn, 'file') || ~exist(fn_zoom, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                colormap bwr
                imagesc((1:nU)/nU, tps, m2plot{pp})
                caxis([-climits(pp), climits(pp)])
                % Add folds to plot
                hold on;
                fons1 = max(1, fons(1)) ;
                fons2 = max(1, fons(2)) ;
                fons3 = max(1, fons(3)) ;
                plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
                plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
                plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

                % title and save
                title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeunits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   

                % Zoom in on small values
                caxis([-climits(pp)/3, climits(pp)/3])
                fn = fullfile(odir, [names{pp} '_zoom.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
                % Zoom in on early times
                ylim([min(tps), max(fons) + 10])
                caxis([-climits(pp)/3, climits(pp)/3])
                disp(['saving ', fn_zoom])
                export_fig(fn_zoom, '-png', '-nocrop', '-r200')   
            end
        end
    end
end

%% Metric Kinematic Correlations
% Plot both all time and select times
timeSpans = {tps, tps(tps < max(fons) + 11)} ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot correlation between terms div(v) and 2Hvn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
padN = round(0.1 * nU) ;
cols = padN:nU - padN ;
if plot_correlations
    corrDir = fullfile(mKDir, 'correlations') ;
    if ~exist(corrDir, 'dir')
        mkdir(corrDir)
    end
    
    for sigma = 0:3
        outputFileNames = {fullfile(corrDir, ...
            sprintf('correlation_sigma%02d_alltime_div_2Hvn', sigma)), ...
            fullfile(corrDir, ...
            sprintf('correlation_sigma%02d_earlytimes_div_2Hvn', sigma))} ;
        alphaVal = 0.6 ;
        sz = 10 ;
        cmap = parula ;
        close all
        set(gcf, 'visible', 'off')
        % Consider each timespan (early or entire series)
        for tspanIdx = 1:2
            fnout = outputFileNames{tspanIdx} ;
            timeSpan_i = timeSpans{tspanIdx} ;
            ntspan = length(timeSpan_i) ;
            titles = {'left lateral', 'right lateral', 'dorsal', 'ventral'} ;
            markers = QS.plotting.markers ;
            colors = mapValueToColor(1:ntspan, [1, ntspan], cmap) ;
            close all
            sphCollection = cell(4, 1) ;
            sposCollection = cell(4, 1) ;
            sphCollection2 = cell(4, 1) ;
            sposCollection2 = cell(4, 1) ;
            sphCollection3 = cell(4, 1) ;
            sposCollection3 = cell(4, 1) ;
            for qq = 1:4  % consider left, right, dorsal, ventral
                disp(['qq = ', num2str(qq), ': ', titles{qq}])
                divv = divvsK{qq + 1}(:, cols) ;    
                H2vn = H2vnsK{qq + 1}(:, cols) ;

                % Optional: smooth here
                if sigma > 0
                    divv = imgaussfilt(divv, sigma);            
                    H2vn = imgaussfilt(H2vn, sigma);  
                end
                
                % Check the smoothing on kymographs
                figure(2) ;
                sphCollection2{qq} = subplot(2, 2, qq) ;
                imagesc(cols/nU, tps, divv); 
                caxis([-climit, climit])
                colormap(bwr256)
                figure(3) ;
                sphCollection3{qq} = subplot(2, 2, qq) ;
                imagesc(cols/nU, tps, H2vn); 
                caxis([-climit, climit])
                colormap(bwr256)

                figure(1) ;
                sphCollection{qq} = subplot(2, 2, qq) ;
                for row = 1:ntspan
                    disp(['row = ', num2str(row)])
                    scatter(divv(row, :), H2vn(row, :), sz, ...
                        markers{qq}, 'MarkerFaceColor', 'none', ...
                        'MarkerEdgeColor', colors(row, :), ...
                        'MarkerEdgeAlpha', alphaVal) ;
                    hold on ;
                end

                % Label the x axis if on the bottom row
                if qq > 2
                    figure(1)
                    xlabel(['$\nabla \cdot \bf{v}_\parallel$ ' unitstr], ...
                            'Interpreter', 'Latex') ;
                    figure(2)
                    xlabel(['ap position, $\zeta/L$'], ...
                        'Interpreter', 'Latex') ;
                    figure(3)
                    xlabel(['ap position, $\zeta/L$'], ...
                        'Interpreter', 'Latex') ;
                end
                figure(1)
                axis equal
                % Add dashed y=x line
                xlims = get(gca, 'xlim') ;
                ylims = get(gca, 'ylim') ;
                xlim([max(-2 * climit, xlims(1)), min(2 * climit, xlims(2))])
                ylim([max(-2 * climit, ylims(1)), min(2 * climit, ylims(2))])
                xlims = get(gca, 'xlim') ;
                ylims = get(gca, 'ylim') ;
                leftdot = max(xlims(1), ylims(1)) ;
                rightdot = min(xlims(2), ylims(2)) ;
                plot([leftdot, rightdot], [leftdot, rightdot], 'k--')

                % Label the y axis if on the left column
                ylabel(['$2Hv_n$ ' unitstr], 'Interpreter', 'Latex') ;
                title(titles{qq}, 'Interpreter', 'Latex')
                
                figure(2)
                ylabel(['time [' QS.timeunits, ']'], 'Interpreter', 'Latex') ;
                title(titles{qq}, 'Interpreter', 'Latex')
                figure(3)
                ylabel(['time [' QS.timeunits, ']'], 'Interpreter', 'Latex') ;
                title(titles{qq}, 'Interpreter', 'Latex')

                % Grab axis position
                sposCollection{qq} = get(sphCollection{qq}, 'Position');
                sposCollection2{qq} = get(sphCollection2{qq}, 'Position');
                sposCollection3{qq} = get(sphCollection3{qq}, 'Position');
            end

            % Move subplots left a bit for colorbar space
            for qq = 1:length(sphCollection)
                spos = sposCollection{qq} ;
                wh = min(spos(3)-0.05, spos(4)) ;
                if mod(qq, 2) == 1
                    set(sphCollection{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                    set(sphCollection2{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                    set(sphCollection3{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                else
                    set(sphCollection{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                    set(sphCollection2{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                    set(sphCollection3{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                end
            end

            % master titles (suptitles)
            figure(1) ;
            sgtitle(['$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$, ',...
                '$\sigma=$', num2str(sigma), ' ', QS.timeunits], ...
                'Interpreter', 'Latex') ;
            
            figure(2) ;
            sgtitle(['$\nabla \cdot \bf{v}_\parallel,$ $\sigma=$', ...
                num2str(sigma), ' ', QS.timeunits], ...
                    'Interpreter', 'Latex') ;
            figure(3) ;
            sgtitle(['$2Hv_n$ $\sigma=$', ...
                num2str(sigma), ' ', QS.timeunits], 'Interpreter', 'Latex') ;
                
            % Add colorbar
            figure(1) ;
            c = colorbar('Position',[.9 .333 .02 .333]) ;
            % Make colorbar share the alpha of the image
            % Manually flush the event queue and force MATLAB to render the colorbar
            % necessary on some versions
            drawnow
            % Get the color data of the object that correponds to the colorbar
            cdata = c.Face.Texture.CData;
            % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
            cdata(end,:) = uint8(alphaVal * cdata(end,:));
            % Ensure that the display respects the alpha channel
            c.Face.Texture.ColorType = 'truecoloralpha';
            % Update the color data with the new transparency information
            c.Face.Texture.CData = cdata;
            c.Label.Interpreter = 'Latex' ;
            c.Label.String = ['time [' QS.timeunits ']'] ;
            c.Ticks = [0, 1] ;
            c.TickLabels = [tps(1), max(timeSpan_i)] ;
            
            figure(2) ;
            c = colorbar('Position',[.9 .333 .02 .333]) ;
            figure(3) ;
            c = colorbar('Position',[.9 .333 .02 .333]) ;
            
            % Save figure
            figure(1)
            saveas(gcf, [fnout '.png']) ;
            figure(2)
            saveas(gcf, [fnout '_kymo_divv.png']) ;
            figure(3)
            saveas(gcf, [fnout '_kymo_H2vn.png']) ;
            close all
            set(gcf, 'visible', 'off')
        end
        disp('done with correlation plots betweeen divv and H2vn')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot correlations between gdot and each term in the sum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_gdot_correlations
    outputFileNames = cell(2, 1) ;
    outputFileNames{1} = {fullfile(mKDir, 'correlation_alltime_div_gdot'), ...
                       fullfile(mKDir, 'correlation_earlytimes_div_gdot')} ;
    outputFileNames{2} = {fullfile(mKDir, 'correlation_alltime_2Hvn_gdot'), ...
                       fullfile(mKDir, 'correlation_earlytimes_2Hvn_gdot')} ;
    alphaVal = 0.6 ;
    sz = 10 ;
    cmap = parula ;
    close all
    set(gcf, 'visible', 'off')
    for pairIdx = 1:2
        for tspanIdx = 1:2
            fnout = outputFileNames{pairIdx}{tspanIdx} ;
            timeSpan_i = timeSpans{tspanIdx} ;
            ntspan = length(timeSpan_i) ;
            titles = {'left lateral', 'right lateral', 'dorsal', 'ventral'} ;
            markers = QS.plotting.markers ;
            colors = mapValueToColor(1:ntspan, [1, ntspan], cmap) ;
            close all
            cols = round(nV * [0.2, 0.85]) ;
            sphCollection = cell(4, 1) ;
            sposCollection = cell(4, 1) ;
            for qq = 1:4  % consider left, right, dorsal, ventral
                disp(['qq = ', num2str(qq), ': ', titles{qq}])
                if pairIdx == 1
                    divv = divvsK{qq + 1} ;
                    gdot = gdotsK{qq + 1} ;
                    sphCollection{qq} = subplot(2, 2, qq) ;
                    for row = 1:ntspan
                        disp(['row = ', num2str(row)])
                        scatter(divv(row, cols), gdot(row, cols), sz, ...
                            markers{qq}, 'MarkerFaceColor', 'none', ...
                            'MarkerEdgeColor', colors(row, :), ...
                            'MarkerEdgeAlpha', alphaVal) ;
                        hold on ;
                    end

                    % Label the x axis if on the bottom row
                    if qq > 2
                        xlabel(['$\nabla \cdot \bf{v}_\parallel$ ' unitstr], ...
                                'Interpreter', 'Latex') ;
                    end
                    axis equal
                    % Add dashed y=x line
                    xlims = get(gca, 'xlim') ;
                    ylims = get(gca, 'ylim') ;
                    leftdot = max(xlims(1), ylims(1)) ;
                    rightdot = min(xlims(2), ylims(2)) ;
                    plot([leftdot, rightdot], [leftdot, rightdot], 'k--')
                else
                    H2vn = H2vnsK{qq + 1} ;
                    gdot = gdotsK{qq + 1} ;
                    sphCollection{qq} = subplot(2, 2, qq) ;
                    for row = 1:ntspan
                        disp(['row = ', num2str(row)])
                        scatter(H2vn(row, cols), gdot(row, cols), sz, ...
                            markers{qq}, 'MarkerFaceColor', 'none', ...
                            'MarkerEdgeColor', colors(row, :), ...
                            'MarkerEdgeAlpha', alphaVal) ;
                        hold on ;
                    end

                    % Label the x axis if on the bottom row
                    if qq > 2
                        xlabel(['$2Hv_n$ ' unitstr], 'Interpreter', 'Latex') ;
                    end
                    axis equal
                    % Add dashed y=x line
                    xlims = get(gca, 'xlim') ;
                    ylims = get(gca, 'ylim') ;
                    leftdot = max(xlims(1), -ylims(2)) ;
                    rightdot = min(xlims(2), -ylims(1)) ;
                    plot([leftdot, rightdot], [-leftdot, -rightdot], 'k--')
                end

                % Label the y axis if on the left column
                if qq == 1 || qq == 3
                    ylabel(['$\textrm{Tr}[g^{-1} \dot{g}]$ ' unitstr], ...
                            'Interpreter', 'Latex')
                end
                title(titles{qq}, 'Interpreter', 'Latex')

                % Grab axis position
                sposCollection{qq} = get(sphCollection{qq}, 'Position');
            end

            % Move subplots left a bit for colorbar space
            for qq = 1:length(sphCollection)
                spos = sposCollection{qq} ;
                wh = min(spos(3)-0.05, spos(4)) ;
                if mod(qq, 2) == 1
                    set(sphCollection{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                else
                    set(sphCollection{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                end
            end

            % Add colorbar
            c = colorbar('Position',[.9 .333 .02 .333]) ;
            % Make colorbar share the alpha of the image
            % Manually flush the event queue and force MATLAB to render the colorbar
            % necessary on some versions
            drawnow
            % Get the color data of the object that correponds to the colorbar
            cdata = c.Face.Texture.CData;
            % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
            cdata(end,:) = uint8(alphaVal * cdata(end,:));
            % Ensure that the display respects the alpha channel
            c.Face.Texture.ColorType = 'truecoloralpha';
            % Update the color data with the new transparency information
            c.Face.Texture.CData = cdata;
            c.Label.Interpreter = 'Latex' ;
            c.Label.String = ['time [' QS.timeunits ']'] ;
            c.Ticks = [0, 1] ;
            c.TickLabels = [tps(1), max(timeSpan_i)] ;

            % Save figure
            saveas(gcf, [fnout '.png']) ;
            saveas(gcf, [fnout '.pdf']) ;
            close all
            set(gcf, 'visible', 'off')
        end
    end
    disp('done')
end


%% Metric Kinematics -- decompose into isotropic and other component
if plot_gdot_decomp
    % Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
    dvDir = fullfile(mKDir, 'avgDV') ;
    lDir = fullfile(mKDir, 'avgLeft') ;
    rDir = fullfile(mKDir, 'avgRight') ;
    dDir = fullfile(mKDir, 'avgDorsal') ;
    vDir = fullfile(mKDir, 'avgVentral') ;
    outdirs = {dvDir, lDir, rDir, dDir, vDir} ;
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    gdotK0 = gdotsK{1} ;
    isogrowth = sum(gdotK0, 2) ;
    
    %% Plot as 1d Curve
    plot(tps, isogrowth); 
    xlabel(['time [' QS.timeunits ']'], 'Interpreter', 'Latex')
    ylabel('Isotropic component of growth') 
    saveas(gcf, fullfile(odir, 'average_growth.png'))  
    
    %% Plot quadrant contributions
    for qq = 1:length(outdirs)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end
        % Unpack what to plot (averaged kymographs, vary averaging region)
        HHK = HHsK{qq} ;
        gdotK = gdotsK{qq} ;
        divvK = divvsK{qq} ;
        velnK = velnsK{qq} ;
        H2vnK = H2vnsK{qq} ;
        m2plot = {gdotK, HHK, divvK, velnK, H2vnK} ;
        titles = {'$\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$',...
            'mean curvature, $H$', ...
            'divergence of flow, $\nabla \cdot \mathbf{v}$', ...
            'normal velocity, $v_n$', ...
            'normal motion, $v_n 2 H$'} ;
        labels = {['$\textrm{Tr}[g^{-1}\dot{g}]$ ' unitstr], ...
            ['mean curvature, $H$ ' Hunitstr], ...
            ['$\nabla \cdot \mathbf{v}$ ' unitstr], ...
            ['normal velocity, $v_n$ ' vunitstr] , ...
            ['normal motion, $v_n 2 H $ ' unitstr]} ;
        names = {'gdot', 'HH', 'divv', 'veln', 'H2vn'} ;
        climits = [climit, climit_H, climit, climit_veln, climit_err] ;

        %% Plot gdot/HH/divv/veln/H2vn DV-averaged kymograph
        for pp = 1:length(m2plot)            
            fn = fullfile(odir, [ names{pp} '.png']) ;
            if ~exist(fn, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                colormap bwr
                imagesc((1:nU)/nU, tps, m2plot{pp})
                caxis([-climits(pp), climits(pp)])
                % Add folds to plot
                hold on;
                fons1 = max(1, fons(1)) ;
                fons2 = max(1, fons(2)) ;
                fons3 = max(1, fons(3)) ;
                plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
                plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
                plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

                % title and save
                title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeunits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
            end
        end
    end
end
