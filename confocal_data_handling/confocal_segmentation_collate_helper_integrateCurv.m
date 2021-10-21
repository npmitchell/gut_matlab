function [intAr, intTh, intQc, intQs, intQa, intQt] = ...
    confocal_segmentation_collate_helper_integrateCurv(timebins, outmatFn)
%confocal_segmentation_collate_helper(datdirs, midX, outmatFn)
%
% Compute resampled area under curve for a/b, theta, Qxx and Qyy of cell
% segmentation in 3d 
%   
%
% NPMitchell 2021
nsE = 250 ;
for ii = 1:length(datdirs)
    
    load(sprintf(outmatFn, ii), ...
        'timestamps', 'ar_ak', 'ar_sk', 'ar_ek', ...
        'th_ak', 'th_sk', 'th_ek', ...
        'Qct_ak', 'Qct_sk', 'Qct_ek', ...
        'Qst_ak', 'Qst_sk', 'Qst_ek', ...
        'Qa_ak', 'Qa_sk', 'Qa_ek', ...
        'Qt_ak', 'Qt_sk', 'Qt_ek')
    
    % resample curve to uniform bins
    dt = diff(timebins) ;
    assert(all(dt == dt(1)))
    
    % preallocate first round
    if ii == 1
        %
        length(find(midx > -eps))
    end
    
    % aspect ratio (raw)
    weights = 1 ./ ar_ek.^2 ;
    [midx, meany, stdy] = ...
        binDataMeanStdWeighted(timestamps, ar_ak, timebins, ...
        weights, nsE, ar_ek) ;
    intAr(ii, :) = cumsum(meany(midx > -eps)) * dt ;
    
    % theta (raw)
    weights = 1 ./ ar_ek.^2 ;
    [midx, meany, stdy] = ...
        binDataMeanStdWeighted(timestamps, th_ak, timebins, ...
        weights, nsE, th_ek) ;
    intTh(ii, :) = cumsum(meany(midx > -eps)) * dt ;
    uncTh(ii, :) = cumsum(meany(midx > -eps)) * dt ;
    
    % Qcmean
    weights = 1 ./ Qc_ek.^2 ;
    [midx, meany, stdy] = ...
        binDataMeanStdWeighted(timestamps, Qc_ak, timebins, ...
        weights, nsE, Qc_ek) ;
    intQc(ii, :) = cumsum(meany(midx > -eps)) * dt ;
    
    % Qsmean
    weights = 1 ./ Qc_ek.^2 ;
    [midx, meany, stdy] = ...
        binDataMeanStdWeighted(timestamps, Qs_ak, timebins, ...
        weights, nsE, Qs_ek) ;
    intQs(ii, :) = cumsum(meany(midx > -eps)) * dt ;
    
    % Qa
    weights = 1 ./ Qc_ek.^2 ;
    [midx, meany, stdy] = ...
        binDataMeanStdWeighted(timestamps, Qa_ak, timebins, ...
        weights, nsE, Qa_ek) ;
    intQa(ii, :) = cumsum(meany(midx > -eps)) * dt ;
    
    % Qt
    weights = 1 ./ Qc_ek.^2 ;
    [midx, meany, stdy] = ...
        binDataMeanStdWeighted(timestamps, Qt_ak, timebins, ...
        weights, nsE, Qt_ek) ;
    intQt(ii, :) = cumsum(meany(midx > -eps)) * dt ;
end

% mean and stdev of curves
binDataMeanStd
