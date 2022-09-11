%% Load pullback pathline metric Kinematics
t0Pathline = 1; 
sresStr = '' ;
nU = tubi.nU ;
outdir = './TubULAR_Paper_Figures/Kymographs_Pathlines/' ;

ntp = length(tubi.xp.fileMeta.timePoints) ;
H2vn_apM = zeros(ntp-1, nU) ;
divv_apM = zeros(ntp-1, nU) ;
radius_apM = zeros(ntp-1, nU) ;
veln_apM = zeros(ntp-1, nU) ;
HH_apM = zeros(ntp-1, nU) ;
gdot_apM = zeros(ntp-1, nU) ;
for tidx = 1:ntp-1
    disp(['loading tidx = '  num2str(tidx)])
    tp = tubi.xp.fileMeta.timePoints(tidx) ;
    mKDir = fullfile(tubi.dir.metricKinematics.root, ...
    strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d'], ...
    tubi.smoothing.lambda, tubi.smoothing.lambda_mesh, ...
    tubi.smoothing.lambda_err, tubi.smoothing.nmodes, ...
    tubi.smoothing.zwidth), '.', 'p'));
    mKPDir = fullfile(mKDir, sprintf('pathline_%04dt0', t0Pathline)) ;
    outdir = fullfile(mKPDir, 'measurements') ;
    Hfn = fullfile(outdir, sprintf('HH_pathline%04d_%06d.mat', t0Pathline, tp))   ;
    efn = fullfile(outdir, sprintf('gdot_pathline%04d_%06d.mat', t0Pathline, tp)) ;
    dfn = fullfile(outdir, sprintf('divv_pathline%04d_%06d.mat', t0Pathline, tp)) ;
    nfn = fullfile(outdir, sprintf('veln_pathline%04d_%06d.mat', t0Pathline, tp)) ;
    rfn = fullfile(outdir, sprintf('radius_pathline%04d_%06d.mat', t0Pathline, tp)) ;
    H2vnfn = fullfile(outdir, sprintf('H2vn_pathline%04d_%06d.mat', t0Pathline, tp)) ;

    tmp1 = load(H2vnfn) ;
    H2vn_apM(tidx, :) = tmp1.H2vn_ap ;
    tmp2 = load(dfn) ;
    divv_apM(tidx, :) = tmp2.divv_ap ;
    tmp = load(rfn) ;
    radius_apM(tidx, :) = tmp.radius_ap ;
    tmp = load(nfn) ;
    veln_apM(tidx, :) = tmp.veln_ap ;
    tmp = load(Hfn) ;
    HH_apM(tidx, :) = tmp.HH_ap ;
    tmp = load(efn) ;
    gdot_apM(tidx, :) = tmp.gdot_ap ;
end


%% Some plotting params

climit = 0.1; % 0.2 ;
climit_err = 0.1; % 0.2 ;
climit_veln = climit * 10 ;
climit_H = climit * 2 ;
climit_radius = 45;

%% Generate Kymograph Plot ================================================
close all; clc;
% figure; 

% Unit definitions for axis labels
convert_to_period = true;
if convert_to_period
    
    % unitstr = '[1/T]' ;
    % Hunitstr = [ '[1/' tubi.spaceUnits ']' ];
    % vunitstr = [ '[' tubi.spaceUnits '/T]' ];
    
    unitstr = '[1/T]' ;
    Hunitstr = [ '[1/' char(956) 'm]' ];
    vunitstr = [ '[' char(956) 'm/T]' ];
    
else
    
    unitstr = [ '[1/' tubi.timeUnits ']' ];
    Hunitstr = [ '[1/' tubi.spaceUnits ']' ];
    vunitstr = [ '[' tubi.spaceUnits '/' tubi.timeUnits ']' ];
    
end

titleadd = ': circumferentially averaged';

plotTypes = {'gdot'} ; % , 'divv', 'H2vn', 'veln', 'HH', 'radius'} ;

for ii = 1:length(plotTypes)
    plotType= plotTypes{ii} ;
    if strcmpi(plotType, 'gdot')
        m2plot = gdot_apM;
        titles = '$\frac{1}{2}\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$';
        labels = ['$\frac{1}{2}\textrm{Tr}[g^{-1}\dot{g}]$ ' unitstr];
        climits = climit;
    elseif strcmpi(plotType, 'HH')
        m2plot = HH_apM;
        % titles = 'mean curvature, $H$';
        % labels = ['mean curvature, $H$ ' Hunitstr];
        titles = 'mean curvature, H';
        labels = ['mean curvature, H ' Hunitstr];
        climits = climit_H;
    elseif strcmpi(plotType, 'divv')
        m2plot = divv_apM;
        % titles = 'divergence of flow, $\nabla \cdot \mathbf{v}$';
        % labels = ['$\nabla \cdot \mathbf{v}$ ' unitstr];
        titles = 'flow divergence';
        labels = ['\nabla \cdot v_{||} ' unitstr];
        climits = climit;
    elseif strcmpi(plotType, 'veln')
        m2plot = veln_apM;
        % titles = 'normal velocity, $v_n$';
        % labels = ['normal velocity, $v_n$ ' vunitstr];
        titles = 'normal velocity, v_n';
        labels = ['normal velocity, v_n ' vunitstr];
        climits = climit_veln;
    elseif strcmpi(plotType, 'H2vn')
        m2plot = H2vn_apM;
        % titles = 'normal motion, $v_n 2 H$';
        % labels = ['normal motion, $v_n 2 H $ ' unitstr];
        titles = 'normal motion';
        labels = ['2 H v_n' unitstr];
        climits = climit_err;
    elseif strcmpi(plotType, 'radius')
        m2plot = radius_apM;
        % titles = 'radius';
        % labels = ['radius [' tubi.spaceUnits ']'];
        titles = 'radius';
        labels = ['radius [' char(956) ']'];
        climits = climit_radius;
    else
        error('Invalid plot type');
    end

    % We relate the normal velocities to the divergence / 2 * H.
    tps = tubi.xp.fileMeta.timePoints(1:end-1) - tubi.t0;

    fig = figure('Visible', 'on',  'units', 'centimeters') ;

    if convert_to_period

        T = 11; % Length of period
        tps = tps ./ T;

        if ismember(lower(plotType), ...
                {'gdot', 'divv', 'veln', 'h2vn'})
            m2plot = T * m2plot;
            climits = T * climits;
        end

        imagesc((1:nU)/nU, tps, m2plot)

        if strcmpi(plotType, 'radius')
            if climits > 0
                % caxis([min(radius_apM(:)), max(radius_apM(:))]);
                caxis([20 45]);
            end
            cmap = brewermap(512, 'RdBu');
            colormap(cmap((257:end).', :));

        else
            if climits > 0
                caxis([-climits, climits])
            end
            colormap(brewermap(256, '*RdBu'));
        end

        axis square

        xticks(0:0.2:1);
        yticks(0:0.5:2.5);

        % title and save
        % title([titles, titleadd], 'Interpreter', 'Latex')
        % ylabel(['time [t/T]'], 'Interpreter', 'Latex')
        % xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
        % cb = colorbar() ;
        % ylabel(cb, labels, 'Interpreter', 'Latex')
        % set(gcf, 'Color', [1 1 1]);

        % title([titles, titleadd], 'FontWeight', 'normal')
        if contains(titles, '$')
            title(titles, 'FontWeight', 'normal', 'interpreter', 'latex');
        else
            title(['\langle' titles char(9002)], 'FontWeight', 'normal');
        end
        ylabel('time [T]');
        xlabel('position [s/L]');

        cb = colorbar() ;
        if contains(titles, '$')
            cb.Label.Interpreter = 'latex';
        end
        cb.Label.String = labels;
        cb.Position(1) = 0.75;
        cb.Position(2) = 0.25;
        cb.Position(3) = 0.045;
        cb.Position(4) = 0.6;
        cb.Label.Position(1) = 3.5;

        axPos = get(gca, 'Position');
        axRatio = axPos(4) / axPos(3);
        axPos(1) = 0.10;
        axPos(2) = 0.20;
        axPos(3) = 0.675;
        axPos(4) = axRatio * axPos(3);
        set(gca, 'Position', axPos);



    else

        imagesc((1:nU)/nU, tps, m2plot)
        if climits > 0
            caxis([-climits, climits])
        end
        colormap(brewermap(256, '*RdBu'));
        axis equal

        % title and save
        title([titles, titleadd], 'Interpreter', 'Latex')
        ylabel(['time [' tubi.timeUnits ']'], 'Interpreter', 'Latex')
        xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
        cb = colorbar() ;
        ylabel(cb, labels, 'Interpreter', 'Latex')

    end

    set(gcf, 'Color', [1 1 1]);

    % Resize Figure for Paper -------------------------------------------------
    set(fig, 'Units', 'centimeters');

    set(fig, 'Position', [0,0,4,4]) ;
    % ratio = fig.Position(4) ./ fig.Position(3);
    % fig.Position(3) = 3.6;
    % fig.Position(4) = ratio * fig.Position(3);

    set(gca, 'FontSize', 10);
    set(gca, 'FontWeight', 'normal');

    set(fig, 'PaperPositionMode', 'auto');
    set(fig.Children, 'FontName', 'Helvetica');
    
    saveas(gcf, fullfile(outdir, [ plotType '.pdf']))
end



%% Generate Cross-Correlation Plot ========================================
close all; clc;

T = 11; % Period size

A = H2vn_apM;
B = divv_apM;

% A = ones(size(H2vn_apM)) ;
% B = ones(size(H2vn_apM)) ;
A = T * A; B = T* B;


close all; clc;

A = H2vn_apM;   
B = divv_apM;

A = T * A; B = T* B;
A = A(4:4+2*T, :); B = B(4:4+2*T, :);

% C = zeros(size(B));
% for i = 0:(size(B,1)-1)
%     for j = 0:(size(B,2)-1)
%         
%         C(i+1,j+1) = sum(A .* circshift(conj(B), [i j]), 'all');
%         
%     end
% end

C = zeros(size(B,1), 1);
shifts = (-size(B,1)+1)*0.5:(size(B,1)-1)*0.5 ;
for idx = 1:length(shifts)
    i = shifts(idx) ;
    C(idx) = sum(A .* circshift(conj(B), [i 0]), 'all');
end
[~, t0_idx] = min(abs(shifts)) ;
timeSteps = shifts / T ;

%% NEW VERSION for plotting

% Fit nonlinear model to correlation function
% fitFunc = @(p, t) p(1) * cos(2 * pi .* t - p(2)) ;
% param0 = [ max(C) 0.5 ];
% options = optimoptions( 'lsqcurvefit', ...
%     'Algorithm', 'levenberg-marquardt', ...
%     'MaxIterations', 1000, ...
%     'MaxFunctionEvaluations', 10000, ...
%     'Display', 'none');
% [param, resnorm, ~, exitflag, output] = ...
%     lsqcurvefit(fitFunc, param0, timeSteps, C', [], [], options);

fitOptions = fitoptions('Method', 'NonlinearLeastSquares', ...
    'Lower', [-Inf,2*pi-eps,-Inf], 'Upper', [Inf, 2*pi+eps, Inf], ...
    'StartPoint', [400, 2*pi, 0]) ;
modelFit = fit(timeSteps', C,'sin1', fitOptions) ;

% Generate Visualization --------------------------------------------------

plot_lw = 0.75;
scatter_size = 7.5;

fig = figure('Visible', 'on',  'units', 'centimeters') ;

% Assumes A and B are the same size
scatter(timeSteps, C, scatter_size, 'o', ...
    'MarkerEdgeColor', tubi.plotting.colors(1,:), ...
    'LineWidth', plot_lw );

hold on

tplot = linspace(timeSteps(1), timeSteps(end), 300);
%plot(tplot, fitFunc(param, tplot), '-', 'LineWidth', plot_lw, ...
%    'Color', tubi.plotting.colors(1,:) );
plot(tplot, modelFit.a1 * sin(modelFit.b1 * tplot + modelFit.c1), ...
    '-', 'LineWidth', plot_lw, 'Color', tubi.plotting.colors(1,:) )

% Zero point of the sin is = -modelFit.c1, so the MAX is 1/4 * T later.
% Multiply modelFit.c1 by 1/2pi to convert to period, then add 1/4 T, where
% T = 1 here because we put timeSteps in units of T. 
maxLoc = -modelFit.c1 / (2*pi) + 0.25 ;
ci = confint(modelFit) ;
ci = ci(:, 3) ;
maxLoc_ci = -ci / (2*pi) + 0.25 ;
assert(all(abs(maxLoc - maxLoc_ci) == abs(maxLoc - maxLoc_ci(1))))
maxLoc_unc = (maxLoc - maxLoc_ci) * 0.5 ;
maxLoc_unc = maxLoc_unc(maxLoc_unc > 0) ;
scatter(maxLoc, fitFunc(param, maxLoc), 1.5 * scatter_size, 'filled', 'k');
plot([maxLoc, maxLoc], [yLim(1) fitFunc(param, maxLoc)], '--k', ...
    'LineWidth', max(0.75 * plot_lw, 0.5))

hold off

% Restrict attention to near zero.
xlim([-1,1])
ylim(yLim);

% xlabel(['time [' tubi.timeUnits ']'], 'Interpreter', 'Latex');
xlabel(['time shift \Delta [T]']);
ylabel({'cross correlation', 'C [2Hv_n(t), \nabla\cdotv_{||}(t+\Delta)]'});


tanno = annotation('textarrow', ...
    (maxLoc-xLim(1))/diff(xLim) * [1.15 1.1], ...
    ( fitFunc(param, maxLoc)-yLim(1))/diff(yLim) * [1.025 0.975], ...
    'String', sprintf('\\Delta = %0.3f\\pm%0.3f T', maxLoc, maxLoc_unc), ...
    'LineWidth', plot_lw, ...
    'FontSize', 10, 'HeadLength', 2.5, 'HeadWidth', 2.5);

set(gcf, 'Color', [1 1 1]);

% Resize Figure for Paper -------------------------------------------------
set(fig, 'Units', 'centimeters');

% ratio = fig.Position(4) ./ fig.Position(3);
% fig.Position(3) = 3.6;
% fig.Position(4) = ratio * fig.Position(3);
set(fig, 'Position', [0, 0, 7, 4])

set(gca, 'FontSize', 10);
set(gca, 'FontWeight', 'normal');

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');
saveas(gcf, fullfile('./TubULAR_Paper_Figures', 'metricKinematics_phaseRelation.pdf'))


%% OLD VERSION

% % tileT = T*2;
% B = B(4:end, :) ;
% tileA1 = A(4:4+T-1,:) ;
% tileA2 = A(4+T:4+2*T-1, :) ;
% AA1 = cat(1, tileA1, tileA1, tileA1, tileA1, tileA1, tileA1, tileA1, tileA1) ;
% AA2 = cat(1, tileA2, tileA2, tileA2, tileA2, tileA2, tileA2, tileA2, tileA2) ;
% 
% corrAB2 = xcorr2(AA1, B);
% % Take the column that only does cross correlation in TIME, not space
% corrAB_1 = corrAB2(:, size(A,2));
% corrAB2 = xcorr2(AA1, B);
% % Take the column that only does cross correlation in TIME, not space
% corrAB_2 = corrAB2(:, size(A,2));
% 
% % Average the two periods together 
% corrAB = (corrAB_1 + corrAB_2)*0.5 ;
%
% timeSteps = -(size(A,1)-1):(size(A,1)-1);
% timeSteps = timeSteps / T;
%
% Find locations of maximum signal
% [maxCorr, maxLoc] = max(corrAB);
% 
% xLim = [timeSteps(1) timeSteps(end)];
% yLim = [-5, 5];
% yLim = 550 * [-1 1];

% 
% %%
% % Fit nonlinear model to correlation function
% % fitFunc = @(p, t) p(1) * cos(p(2) .* (t-p(4)) - p(3)) .* exp(-(t-p(4)).^2 / p(5).^2);
% % param0 = [ min(corrAB) 2*pi 0 0 1];
% fitFunc = @(p, t) p(1) * cos(2 * pi .* (t-p(2)) - p(3)) .* exp(-(t-p(2)).^2 / p(4).^2);
% fitFunc = @(p, t) p(1) * cos(2 * pi .* (t-p(2)) - p(3)) .* exp(-(t-p(2)).^2 / p(4).^2);
% param0 = [ min(corrAB) 0 0 1];
% options = optimoptions( 'lsqcurvefit', ...
%     'Algorithm', 'levenberg-marquardt', ...
%     'MaxIterations', 1000, ...
%     'MaxFunctionEvaluations', 10000, ...
%     'Display', 'none');
% [param, resnorm, ~, exitflag, output] = ...
%     lsqcurvefit(fitFunc, param0, timeSteps, corrAB.', [], [], options);
% 
% % Generate Visualization --------------------------------------------------
% 
% plot_lw = 0.75;
% scatter_size = 7.5;
% 
% fig = figure('Visible', 'on',  'units', 'centimeters') ;
% 
% % Assumes A and B are the same size
% scatter(timeSteps, corrAB, scatter_size, 'o', ...
%     'MarkerEdgeColor', tubi.plotting.colors(1,:), ...
%     'LineWidth', plot_lw );
% 
% hold on
% 
% tplot = linspace(timeSteps(1), timeSteps(end), 300);
% plot(tplot, fitFunc(param, tplot), '-', 'LineWidth', plot_lw, ...
%     'Color', tubi.plotting.colors(1,:) );
% 
% maxLoc = (2*pi*param(2) + param(3))  +1;
% scatter(maxLoc, fitFunc(param, maxLoc), 1.5 * scatter_size, 'filled', 'k');
% plot([maxLoc, maxLoc], [yLim(1) fitFunc(param, maxLoc)], '--k', ...
%     'LineWidth', max(0.75 * plot_lw, 0.5))
% 
% hold off
% 
% xlim(xLim);
% ylim(yLim);
% 
% % xlabel(['time [' tubi.timeUnits ']'], 'Interpreter', 'Latex');
% xlabel(['time shift ' char(916) ' [T]']);
% ylabel({'cross correlation', 'C [2Hv_n, \nabla\cdotv_{||}]'});
% 
% 
% tanno = annotation('textarrow', ...
%     (maxLoc-xLim(1))/diff(xLim) * [1.15 1.1], ...
%     ( fitFunc(param, maxLoc)-yLim(1))/diff(yLim) * [1.025 0.975], ...
%     'String', sprintf('\\Delta = %0.2f T', maxLoc), ...
%     'LineWidth', plot_lw, ...
%     'FontSize', 5, 'HeadLength', 2.5, 'HeadWidth', 2.5);
% 
% set(gcf, 'Color', [1 1 1]);
% 
% % Resize Figure for Paper -------------------------------------------------
% set(fig, 'Units', 'centimeters');
% 
% ratio = fig.Position(4) ./ fig.Position(3);
% fig.Position(3) = 3.6;
% fig.Position(4) = ratio * fig.Position(3);
% 
% set(gca, 'FontSize', 5);
% set(gca, 'FontWeight', 'normal');
% 
% set(fig, 'PaperPositionMode', 'auto');
% set(fig.Children, 'FontName', 'Helvetica');
% 
% 
% 
% %%
% 
% 
% plot(timeSteps, corrAB, 'LineWidth', 3);
% hold on
% scatter(timeSteps(maxLoc), maxCorr, 100, 'filled', 'k');
% plot([timeSteps(maxLoc), timeSteps(maxLoc)], [yLim(1) maxCorr], '--k', ...
%     'LineWidth', 3)
% hold off
% 
% 
% 
% 
% % tanno = annotation('textarrow', ...
% %     (timeSteps(maxLoc)-xLim(1))/diff(xLim) * [1.05 1], ...
% %     (maxCorr-yLim(1))/diff(yLim) * [1.025 0.975], ...
% %     'String', sprintf('$$\\Delta \\frac{t}{T} = %0.2f$$', timeSteps(maxLoc)), ...
% %     'LineWidth', 2, 'Interpreter', 'Latex' );
% tanno = annotation('textarrow', ...
%     (timeSteps(maxLoc)-xLim(1))/diff(xLim) * [1.075 1], ...
%     (maxCorr-yLim(1))/diff(yLim) * [1.00 0.975], ...
%     'String', sprintf('$$\\Delta = %0.2f$$ T', timeSteps(maxLoc)), ...
%     'LineWidth', 2, 'Interpreter', 'Latex' );
% 
% 
% 
% %% 
% % get rid of the first and last few timepoints and endcaps
% ecap = 3 ;
% H2vnM = H2vn_apM(3:end-3, ecap:end-ecap) ;
% divvM = divv_apM(3:end-3, ecap:end-ecap) ;
% binscatter(H2vn_apM(:), divv_apM(:), 100)

