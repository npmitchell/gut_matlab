% Compute psiN for some point sets on a ring
% choose initial angles
clear; close all
thetafix = [0 1 2 3 4 5] * 2*pi/6 ;
% global x axis 
theta0s = (0:0.02:1)*pi/3 ;

% Options
% -------
% disorder strength
deltas = 0:0.001:1 ;
nn = length(thetafix) ;
% weight ranges from 1 for <few lines to 0 for 1000 lines
weight = max(0, 1 - length(deltas) / 15000000) ;
cmap = colormap ;
    
for ii = 1:length(theta0s)
    disp(['ii = ' num2str(ii)])
    close all
    fig = figure('visible', 'off') ;
    set(get(fig,'Children'),'Visible','off');
    thetas = theta0s(ii) + thetafix ;
    for k = 1:length(deltas)
        delta = deltas(k) ;
        dt = [0 delta*(rand(1, nn-1)-0.5)*2*pi/nn ] ;
        theta = thetas + dt ;
        psi(k) = 1. / nn * sum(exp(1j * nn * theta)) ;

        sp1 = subplot(1,2,1, 'visible', 'off') ;    
        cind = max(uint8(k / length(deltas) * length(cmap)), 1) ;
        cind = min(cind, length(cmap)) ;
        cinds(k) = cind ;

        if mod(k, length(deltas)/100) < 1
            for qq = 1:length(theta)
                plot([0, cos(theta(qq))], [0, sin(theta(qq))], '-', ...
                    'color', [1 1 1]*(1-weight) + weight*cmap(cind, :), ...
                    'linewidth', 0.01)
                hold on;
            end
        end
    end
    
    % Plot the orientational order parameters
    plot([0, cos(theta0s(ii))], [0, sin(theta0s(ii))], 'k-')
    axis equal
    axis off
    suptitle(['\theta_0 = ' sprintf('%0.2f', theta0s(ii) / pi) '\pi'])
    sp2 = subplot(1,2,2, 'visible', 'off') ;
    scatter(real(psi), imag(psi), 5, cmap(cinds, :))
    xlim([-1.1, 1.1])
    ylim([-1.1, 1.1])
    axis equal    
    xlim([-1.1, 1.1])
    ylim([-1.1, 1.1])
    hold on;
    plot(cos(0:0.01:2*pi), sin(0:0.01:2*pi), 'k-')
    xlabel(['\Re \psi_' num2str(nn)])
    ylabel(['\Im \psi_' num2str(nn)])
    set(gcf, 'visible', 'off')
    saveas(gcf, sprintf(['psi' num2str(nn) '_%06d.png'], ii))
    
    % Alt figure
    subplot(1,2,2)
    cla
    if length(deltas) < 2000
        sc = scatter(deltas, abs(psi), 5, cmap(cinds, :),...
            'filled','markeredgecolor', 'none', 'MarkerFaceAlpha',.3) ;
    else
        sc = scatter(deltas, abs(psi), 5, cmap(cinds, :),...
            'filled','markeredgecolor', 'none', 'MarkerFaceAlpha',.1) ;
    end
    xlim([0, 1])
    ylim([0, 1])
    xlabel('bond angle disorder')
    ylabel(['\psi_' num2str(nn)])
    saveas(gcf, sprintf(['psialt' num2str(nn) '_%06d.png'], ii))
end

%% 
% Compute psiN for some point sets on a ring varying only one bond
% choose initial angles
clear; close all
thetafix = [0 1 2 3 4 5] * 2*pi/6 ;
% global x axis 
theta0s = (0:0.02:1)*pi/3 ;

% Options
% -------
% disorder strength
deltas = fliplr(0:0.001:1) ;
nn = length(thetafix) ; 
% weight ranges from 1 for <few lines to 0 for 1000 lines
weight = max(0, 1 - length(deltas) / 15000000) ;
cmap = flipud(colormap) ;

for ii = 1:length(theta0s)
    disp(['ii = ' num2str(ii)])
    close all
    fig = figure('visible', 'off') ;
    set(get(fig,'Children'),'Visible','off');
    thetas = theta0s(ii) + thetafix ;
    for k = 1:length(deltas)
        delta = deltas(k) ;
        dt = zeros(1, nn) ;
        dt(1) = delta*(rand(1, 1)-0.5)*2*pi/nn ;
        theta = thetas + dt ;
        psi(k) = 1. / nn * sum(exp(1j * nn * theta)) ;

        sp1 = subplot(1,2,1, 'visible', 'off') ;    
        cind = max(uint8(k / length(deltas) * length(cmap)), 1) ;
        cind = min(cind, length(cmap)) ;
        cinds(k) = cind ;

        if mod(k, length(deltas)/100) < 1
            for qq = 1:length(theta)
                plot([0, cos(theta(qq))], [0, sin(theta(qq))], '-', ...
                    'color', [1 1 1]*(1-weight) + weight*cmap(cind, :), ...
                    'linewidth', 0.01)
                hold on;
            end
        end
    end

    % Plot the orientational order parameters
    plot([0, cos(theta0s(ii))], [0, sin(theta0s(ii))], 'k-')
    axis equal
    axis off
    suptitle(['\theta_0 = ' sprintf('%0.2f', theta0s(ii) / pi) '\pi'])
    sp2 = subplot(1,2,2, 'visible', 'off') ;
    scatter(real(psi), imag(psi), 5, cmap(cinds, :), 'filled', 'markeredgecolor', 'none')
    xlim([-1.1, 1.1])
    ylim([-1.1, 1.1])
    axis equal    
    xlim([-1.1, 1.1])
    ylim([-1.1, 1.1])
    hold on;
    plot(cos(0:0.01:2*pi), sin(0:0.01:2*pi), 'k-')
    xlabel(['\Re \psi_' num2str(nn)])
    ylabel(['\Im \psi_' num2str(nn)])
    set(gcf, 'visible', 'off')
    saveas(gcf, sprintf(['psi' num2str(nn) '_%06d.png'], ii))
    
    subplot(1,2,2)
    cla
    sc = scatter(deltas, abs(psi), 5, cmap(cinds, :),...
        'filled','markeredgecolor', 'none', 'MarkerFaceAlpha',.1) ;
    xlim([0, 1])
    ylim([0, 1])
    xlabel('bond angle disorder')
    ylabel(['\psi_' num2str(nn)])
    saveas(gcf, sprintf(['psialt' num2str(nn) '_%06d.png'], ii))
end

%% 
% Compute a bastardized version of psiN for some point sets on a ring 
% where we use interior angles instead of bond angles
% choose initial angles
clear; close all
thetafix = [0 1 2 3 4 5] * 2*pi/6 ;
% global x axis 
theta0s = 0 % (0:0.02:1)*pi/3 ;

% Options
% -------
% disorder strength
deltas = 0:0.0001:1 ;
nn = length(thetafix) ;
% weight ranges from 1 for <few lines to 0 for 1000 lines
weight = max(0, 1 - length(deltas) / 15000000) ;
cmap = colormap ;
    
for ii = 1:length(theta0s)
    disp(['ii = ' num2str(ii)])
    close all
    fig = figure('visible', 'off') ;
    set(get(fig,'Children'),'Visible','off');
    thetas = theta0s(ii) + thetafix ;
    for k = 1:length(deltas)
        delta = deltas(k) ;
        dt = [0 delta*(rand(1, nn-1)-0.5)*2*pi/nn ] ;
        theta = thetas + dt ;
        dtheta = [diff(theta), (theta(1) + 2*pi) - theta(end)] ;
        psi(k) = 1. / nn * sum(exp(1j * nn * dtheta)) ;

        sp1 = subplot(1,2,1, 'visible', 'off') ;    
        cind = max(uint8(k / length(deltas) * length(cmap)), 1) ;
        cind = min(cind, length(cmap)) ;
        cinds(k) = cind ;

        if mod(k, length(deltas)/100) < 1
            for qq = 1:length(theta)
                plot([0, cos(theta(qq))], [0, sin(theta(qq))], '-', ...
                    'color', [1 1 1]*(1-weight) + weight*cmap(cind, :), ...
                    'linewidth', 0.01)
                hold on;
            end
        end
    end
    
    % Plot the orientational order parameters
    plot([0, cos(theta0s(ii))], [0, sin(theta0s(ii))], 'k-')
    axis equal
    axis off
    suptitle(['\theta_0 = ' sprintf('%0.2f', theta0s(ii) / pi) '\pi'])
    sp2 = subplot(1,2,2, 'visible', 'off') ;
    scatter(real(psi), imag(psi), 5, cmap(cinds, :))
    xlim([-1.1, 1.1])
    ylim([-1.1, 1.1])
    axis equal    
    xlim([-1.1, 1.1])
    ylim([-1.1, 1.1])
    hold on;
    plot(cos(0:0.01:2*pi), sin(0:0.01:2*pi), 'k-')
    xlabel(['\Re \psi_' num2str(nn)])
    ylabel(['\Im \psi_' num2str(nn)])
    set(gcf, 'visible', 'off')
    saveas(gcf, sprintf(['phi' num2str(nn) '_%06d.png'], ii))
    
    subplot(1,2,2)
    cla
    sc = scatter(deltas, real(psi), 5, cmap(cinds, :),...
        'filled','markeredgecolor', 'none', 'MarkerFaceAlpha',.1) ;xlim([0, 1])
    xlabel('bond angle disorder')
    ylabel(['\psi_' num2str(nn)])
    saveas(gcf, sprintf(['phialt' num2str(nn) '_%06d.png'], ii))
    
    
end