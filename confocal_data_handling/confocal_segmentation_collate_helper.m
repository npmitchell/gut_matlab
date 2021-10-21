function confocal_segmentation_collate_helper(datdirs, midX, outmatFn,...
    timeUnits, shadingStyle, color)
%confocal_segmentation_collate_helper(datdirs, midX, outmatFn)
%
% Compute and plot statistics for a/b, theta, Qxx and Qyy of cell
% segmentation in 3d 
%   
%
% NPMitchell 2021
nsE = 250 ;
for ii = 1:length(datdirs)
    datdir = datdirs{ii} ;
    load(fullfile(datdir, 'results.mat'), ...
        'stats', 'keep', 'timestamps', 'exciteIdx') ;
    ntp = length(timestamps) ;
    ar_ak = zeros(ntp, 1) ; 
    ar_sk = zeros(ntp, 1) ;
    ar_ek = zeros(ntp, 1) ;
    th_ak = zeros(ntp, 1) ; 
    th_sk = zeros(ntp, 1) ;
    th_ek = zeros(ntp, 1) ;
    Qa_ak = zeros(ntp, 1) ;
    Qa_sk = zeros(ntp, 1) ;
    Qa_ek = zeros(ntp, 1) ;
    Qt_ak = zeros(ntp, 1) ;
    Qt_sk = zeros(ntp, 1) ;
    Qt_ek = zeros(ntp, 1) ;
    Qct_ak = zeros(ntp, 1) ;
    Qct_sk = zeros(ntp, 1) ;
    Qct_ek = zeros(ntp, 1) ;
    Qst_ak = zeros(ntp, 1) ;
    Qst_sk = zeros(ntp, 1) ;
    Qst_ek = zeros(ntp, 1) ;
    idx2do = find(keep) ;
    fullResults = cell(ntp, 1) ;
    disp('collating stats for this dataset')
    for tidx = 1:ntp
        idx = idx2do(tidx) ;

        lobe12 = [0, midX(ii), Inf] ;
        m1 = stats{idx}.raw.moment1 ;
        m2 = stats{idx}.raw.moment2 ;
        aa = stats{idx}.raw.areas ;
        cx = stats{idx}.raw.centroids2d ;
        a1 = stats{idx}.raw.angle ;
        ar = sqrt(m2./m1) ;
        % Compute the means/std/ste weighted
        [~, meanR, stdR, nR, steR] = binDataMeanStdWeighted(cx, ...
            ar, lobe12, aa, nsE) ;
        qcos2 = (sqrt(m2./m1)-1) .* cos(2 * a1) ;
        qsin2 = (sqrt(m2./m1)-1) .* sin(2 * a1) ;
        
        strength = sqrt(m2./m1)-1 ;
        Q0 = QtensorFromAngle(a1) ;
        keep = cx < midX(ii) ;
        QQ = strength(keep) .* Q0(keep, :, :) ;
        [meanQ, stdQ, steQ] = QtensorStats(QQ, aa(keep)) ;
        [meanQLobeAspect, meanQLobeTheta, meanQLobeAspectStd, ...
            meanQLobeThetaStd] = QtensorAspectTheta(meanQ, stdQ) ;
        [~, ~, meanQLobeAspectSte, ...
            meanQLobeThetaSte] = QtensorAspectTheta(meanQ, steQ) ;
        
        [~, meanC, stdC, nC, steC] = binDataMeanStdWeighted(cx, ...
            qcos2, lobe12, aa, nsE) ;
        [~, meanS, stdS, nS, steS] = binDataMeanStdWeighted(cx, ...
            qsin2, lobe12, aa, nsE) ;
        [~, meanCt, stdCt, nCt, steCt] = binDataMeanStdWeighted(cx, ...
            qcos2, lobe12, aa, nsE) ;
        [~, meanSt, stdSt, nSt, steSt] = binDataMeanStdWeighted(cx, ...
            qsin2, lobe12, aa, nsE) ;
        meanT =  mod(0.5 * atan2(meanSt(1), meanCt(1)), pi) ;
        % stdTheta
        numerator1 = meanSt(1)^2 * (stdCt(1))^2 ;
        numerator2 = meanCt(1)^2 * (stdSt(1))^2 ;
        numerator = numerator1+numerator2 ;
        denominator = (meanCt(1)^2 + meanSt(1)^2)^2 ;
        stdT = 0.5 * sqrt(numerator / denominator) ;
        % steTheta
        numerator1 = meanSt(1)^2 * (steCt(1))^2 ;
        numerator2 = meanCt(1)^2 * (steSt(1))^2 ;
        numerator = numerator1+numerator2 ;
        denominator = (meanCt(1)^2 + meanSt(1)^2)^2 ;
        steT = 0.5 * sqrt(numerator / denominator) ;

        % Collate
        ar_ak(tidx) = meanR(1) ;
        ar_sk(tidx) = stdR(1) ;
        ar_ek(tidx) = steR(1) ;
        th_ak(tidx) = meanT ;
        th_sk(tidx) = stdT ;
        th_ek(tidx) = steT ;

        Qa_ak(tidx) = meanQLobeAspect ;
        Qa_sk(tidx) = meanQLobeAspectStd ;
        Qa_ek(tidx) = meanQLobeAspectSte ;
        Qt_ak(tidx) = mod(meanQLobeTheta, pi) ;
        Qt_sk(tidx) = meanQLobeThetaStd ;
        Qt_ek(tidx) = meanQLobeThetaSte ;
        
        Qct_ak(tidx) = meanC(1) ;
        Qct_sk(tidx) = stdC(1) ;
        Qct_ek(tidx) = steC(1) ;

        Qst_ak(tidx) = meanS(1) ;
        Qst_sk(tidx) = stdS(1) ;
        Qst_ek(tidx) = steS(1) ;

        fullResults{tidx} = stats{idx} ;
        fullResults{tidx}.keep = keep ;
        fullResults{tidx}.QQ = Q0 ;
        fullResults{tidx}.strength = strength ;
    end

    if contains(lower(timeUnits), 'hr') || contains(lower(timeUnits), 'hour')
        timestamps = timestamps /60 ;
    end
    % aR
    subplot(3, 2, 5)
    hold on;
    if shadingStyle == 1
        lineProps = {'.-', 'color', color} ;
        h1=shadedErrorBar(timestamps, ar_ak, ar_sk, 'lineProps', lineProps) ;
        plot(timestamps, ar_ak-ar_ek, '--', 'color', color) ;
        plot(timestamps, ar_ak+ar_ek, '--', 'color', color) ;
        pbaspect([1, 1.5, 1])
    else
        lineProps = {'.-', 'color', color} ;
        h1=shadedErrorBar(timestamps, ar_ak, ar_ek, 'lineProps', lineProps) ;
        pbaspect([1, 1, 1])
    end
    xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
    ylabel('$a/b$', 'interpreter', 'latex')

    % theta
    subplot(3, 2, 6)
    hold on;
    if shadingStyle == 1
        lineProps = {'.-', 'color', color} ;
        h1=shadedErrorBar(timestamps, th_ak, th_sk, 'lineProps', lineProps) ;
        plot(timestamps, th_ak-th_ek, '--', 'color', color) ;
        plot(timestamps, th_ak+th_ek, '--', 'color', color) ;
        pbaspect([1, 1.5, 1])
    else
        lineProps = {'.-', 'color', color} ;
        h1=shadedErrorBar(timestamps, th_ak, th_ek, 'lineProps', lineProps) ;
        pbaspect([1, 1, 1])
    end
    xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
    ylabel('$\theta$', 'interpreter', 'latex')
    ylim([0, pi])
    yticks([0, pi/2, pi])
    set(gca, 'yticklabels', {'0', '\pi/2', '\pi'})

    % Qxx
    subplot(3, 2, 1)
    hold on;
    if shadingStyle == 1
        lineProps = {'.-', 'color', color} ;
        h1=shadedErrorBar(timestamps, 2*Qct_ak, 2*Qct_sk, ...
            'lineProps', lineProps) ;
        plot(timestamps, 2*(Qct_ak-Qct_ek), '--', 'color', color) ;
        plot(timestamps, 2*(Qct_ak+Qct_ek), '--', 'color', color) ;
        pbaspect([1, 1.5, 1])
    else
        lineProps = {'.-', 'color', color} ;
        h1=shadedErrorBar(timestamps, 2*Qct_ak, 2*Qct_ek,...
            'lineProps', lineProps) ;
        pbaspect([1, 1, 1])
    end
    xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
    ylabel('$Q_{xx}$', 'interpreter', 'latex')

    % Qyy
    subplot(3, 2, 2)
    hold on;
    if shadingStyle == 1
        lineProps = {'.-', 'color', color} ;
        h1=shadedErrorBar(timestamps, 2*Qst_ak, 2*Qct_sk, ...
            'lineProps', lineProps) ;
        plot(timestamps, 2*(Qst_ak-Qst_ek), '--', 'color', color) ;
        plot(timestamps, 2*(Qst_ak+Qst_ek), '--', 'color', color) ;
        pbaspect([1, 1.5, 1])
    else
        % h1=plot(timestamps, Qst_ak, '.-', 'color', color) ;
        lineProps = {'.-', 'color', color} ;
        h1=shadedErrorBar(timestamps, 2*Qst_ak, 2*Qct_ek, ...
            'lineProps', lineProps) ;
        pbaspect([1, 1, 1])
    end
    xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
    ylabel('$Q_{yy}$', 'interpreter', 'latex')

    
    % Qaspect
    subplot(3, 2, 3)
    hold on;
    if shadingStyle == 1
        lineProps = {'.-', 'color', color} ;
        h1=shadedErrorBar(timestamps, 2*Qa_ak + 1, 2*Qa_sk + 1, ...
            'lineProps', lineProps) ;
        plot(timestamps, 2*(Qa_ak-Qa_ek) + 1, '--', 'color', color) ;
        plot(timestamps, 2*(Qa_ak+Qa_ek) + 1, '--', 'color', color) ;
        pbaspect([1, 1.5, 1])
    else
        lineProps = {'.-', 'color', color} ;
        h1=shadedErrorBar(timestamps, 2*Qa_ak, 2*Qa_ek,...
            'lineProps', lineProps) ;
        pbaspect([1, 1, 1])
    end
    xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
    ylabel('$(a/b)_{\langle Q \rangle}$', 'interpreter', 'latex')

    % Qtheta
    subplot(3, 2, 4)
    hold on;
    if shadingStyle == 1
        lineProps = {'.-', 'color', color} ;
        h1=shadedErrorBar(timestamps, Qt_ak, Qt_sk, ...
            'lineProps', lineProps) ;
        plot(timestamps, Qt_ak-Qt_ek, '--', 'color', color) ;
        plot(timestamps, Qt_ak+Qt_ek, '--', 'color', color) ;
        pbaspect([1, 1.5, 1])
    else
        % h1=plot(timestamps, Qst_ak, '.-', 'color', color) ;
        lineProps = {'.-', 'color', color} ;
        h1=shadedErrorBar(timestamps, Qt_ak, Qt_ek, ...
            'lineProps', lineProps) ;
        pbaspect([1, 1, 1])
    end
    ylim([0, pi])
    yticks([0, pi/2, pi])
    set(gca, 'yticklabels', {'0', '\pi/2', '\pi'})
    xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
    ylabel('$\theta_{\langle Q \rangle}$', 'interpreter', 'latex')

    
    save(sprintf(outmatFn, ii), ...
        'timestamps', 'ar_ak', 'ar_sk', 'ar_ek', ...
        'th_ak', 'th_sk', 'th_ek', ...
        'Qct_ak', 'Qct_sk', 'Qct_ek', ...
        'Qst_ak', 'Qst_sk', 'Qst_ek', ...
        'Qa_ak', 'Qa_sk', 'Qa_ek', ...
        'Qt_ak', 'Qt_sk', 'Qt_ek', ...
        'timestamps', 'exciteIdx', ...
        'fullResults')
    end
end
