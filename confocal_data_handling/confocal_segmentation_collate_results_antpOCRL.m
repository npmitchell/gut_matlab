%% Collate confocal_segmentation_antpOCRL
% npmitchell 2021
clearvars
addpath('/mnt/data/code/gut_matlab/addpath_recurse/')
addpath_recurse('/mnt/data/code/gut_matlab/')
colors = define_colors ;

timeUnits = 'hr' ;
xlims = [-0.5, 1.125] ;
outdir = '/mnt/data/optogenetics_confocal/antpGAL4/' ;
OCRLdirs = {'/mnt/data/optogenetics_confocal/antpGAL4/huygens_deconvolution_noKlar/202103181730_antpG4OCRLGap43mCh_40x1p6x_5mpf_4pc3pc_to_12pc9pc_600ns_lav3_DC/cellSegmentation/', ...
    '/mnt/data/optogenetics_confocal/antpGAL4/huygens_deconvolution_withKlar/202105252111_AntpG4kOCRLgap43_0p75um_1p25x40x_lav3_3t6pc3t6pc_5mpf_480ns_LED4/cellSegmentation/ground_truth/',...
    '/mnt/data/optogenetics_confocal/antpGAL4/huygens_deconvolution_withKlar/202105311838_AntpG4kOCRL_0p75um_1p5x40_3t6pc_lav3_86s_5mpf/cellSegmentation/',...
  } ;
OCRLmidX = [665, 575, 570] ;
% OCRL: Left side - 202103181730_antpG4OCRLGap43mCh_40x1p6x_5mpf_4pc3pc_to_12pc9pc_600ns_lav3_DC
%   midConstrictionX = 
% Left side - 202105252111_AntpG4kOCRLgap43_0p75um_1p25x40x_lav3_3t6pc3t6pc_5mpf_480ns_LED4
%   midConstrictionX = 
% Right side - 202105311838_AntpG4kOCRL_0p75um_1p5x40_3t6pc_lav3_86s_5mpf
%   midConstrictionX = 


% poster frames: use 202105252111


WTdirs = {'/mnt/data/optogenetics_confocal/WTcontrol_48YG4kCAAXmCh/202106081815_48YG4kCAAXmCh_0p75um_1p5x40x_3t6pc3t6pc_lav3_615ns_5mpf_leftDorsoLateral_included/cellSegmentation/', ...
    } ; % '/mnt/data/optogenetics_confocal/WTcontrol_48YG4kCAAXmCh/'} ;
WTmidX = [345] ;
% Left -- 202106081815_48YG4kCAAXmCh_0p75um_1p5x40x_3t6pc3t6pc_lav3_615ns_5mpf_leftDorsoLateral_included
% Right -- 202106061730_48YCAAXmCh_0p75um_1p25x40x_2p5t5pc_lav3_615ns_rightLateral_include

ANTPdirs = {'/mnt/data/confocal_data/gut/AntpNSRC3_48YG4CAAXmCh/202108291112_antpNSRC3_48YG4CAAXmCh_2um_5pc_10mpf_lav3fac2_e1e2e3e5mutant/e1/cellSegmentation/', ...
    '/mnt/data/confocal_data/gut/AntpNSRC3_48YG4CAAXmCh/20210827438_Antp48YG4CAAXmCh_control_1p2x40x_2um_10mpf_2t5pc_lav4fac2/e2/cellSegmentation/'};
ANTPmidX = [422, 393] ;

% 202108291112_antpNSRC3_48YG4CAAXmCh_2um_5pc_10mpf_lav3fac2_e1e2e3e5mutant/e1
% 20210827438_Antp48YG4CAAXmCh_control_1p2x40x_2um_10mpf_2t5pc_lav4fac2/e2

%% ASPECT RATIO
for shadingStyle = [2, 1]
    close all
    % fig = figure('units', 'cm', 'outerposition', [0.5 0 0.5 1]);
    
    hf = figure('Position', [100 100 320 600], 'units', 'centimeters');
    hold on;
    outmatFn = fullfile(outdir, 'OCRL%03d.mat') ;
    confocal_segmentation_collate_helper(OCRLdirs, OCRLmidX, outmatFn, ...
        timeUnits, shadingStyle, colors(1, :))

    % Save figure
    saveas(gcf, fullfile(outdir, sprintf('aRQ_results_OCRL_style%d.png', shadingStyle)))
    saveas(gcf, fullfile(outdir, sprintf('aRQ_results_OCRL_style%d.pdf', shadingStyle)))
    
    

    %% Compare to WT
    outmatFn = fullfile(outdir, 'WT%03d.mat') ;
    confocal_segmentation_collate_helper(WTdirs, WTmidX, outmatFn, ...
        timeUnits, shadingStyle, colors(2, :))

    % Save figure
    saveas(gcf, fullfile(outdir, sprintf('aRQ_results_OCRL_WT_confocal_style%d.png', shadingStyle)))
    saveas(gcf, fullfile(outdir, sprintf('aRQ_results_OCRL_WT_confocal_style%d.pdf', shadingStyle)))
    
    
    
    %% Compare to ANTP
    outmatFn = fullfile(outdir, 'ANTP%03d.mat') ;
    confocal_segmentation_collate_helper(ANTPdirs, ANTPmidX, outmatFn, ...
        timeUnits, shadingStyle, colors(3, :))
    % Save figure
    saveas(gcf, fullfile(outdir, sprintf('aRQ_results_OCRL_WT_ANTP_confocal_style%d.png', shadingStyle)))
    saveas(gcf, fullfile(outdir, sprintf('aRQ_results_OCRL_WT_ANTP_confocal_style%d.pdf', shadingStyle)))
        
    %% Compare to WT Lightsheet
    % % WTfn =  '/mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/msls_output/cellSegmentation/seg3d_corrected/stats_summary.mat' ;
    WTfn =  '/mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/msls_output/cellSegmentation/seg3d_corrected/stats_summary_L12.mat' ;
    tmp = load(WTfn) ;
    
    if contains(lower(timeUnits), 'hr') || contains(lower(timeUnits), 'hour')
        tmp.timeStamps = tmp.timeStamps /60 ;
    end
    
    % Qxx
    subplot(3, 2, 1)
    lineProps = {'.-', 'color', colors(2, :)} ;
    h2 = shadedErrorBar(tmp.timeStamps, tmp.meanxsL12, ...
        0.5 * tmp.stdxs, 'lineProps', lineProps) ;
    plot(tmp.timeStamps, (tmp.meanxsL12 + tmp.stdmeanxs), '--', 'color', colors(2, :)) ;
    plot(tmp.timeStamps, (tmp.meanxsL12 - tmp.stdmeanxs), '--', 'color', colors(2, :)) ;
    xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
    ylabel('$Q_{xx}$', 'interpreter', 'latex')
    
    % Qyy
    subplot(3, 2, 2)
    lineProps = {'.-', 'color', colors(2, :)} ;
    h2 = shadedErrorBar(tmp.timeStamps, tmp.meanysL12, ...
        tmp.stdys, 'lineProps', lineProps) ;
    plot(tmp.timeStamps, (tmp.meanysL12 + tmp.stdmeanys), '--', 'color', colors(2, :)) ;
    plot(tmp.timeStamps, (tmp.meanysL12 - tmp.stdmeanys), '--', 'color', colors(2, :)) ;
    xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
    ylabel('$Q_{yy}$', 'interpreter', 'latex')
    
    % Save figure
    saveas(gcf, fullfile(outdir, 'aR_results_OCRL_WT_ANTP_wLightsheet.png'))
    saveas(gcf, fullfile(outdir, 'aR_results_OCRL_WT_ANTP_wLightsheet.pdf'))
    
end

%% Group similar timepoints together
% timebins = [-22.5:10:90] / 60 ;
timebins = [-30:20:90] / 60 ;
xlims = [-0.3, 1.125] ;
clf

%WT
Qa_ak = [] ;
Qt_ak = [] ;
Qa_ek = [] ;
Qt_ek = [] ;
Qct_ak = [] ;
Qst_ak = [] ;
Qct_ek = [] ;
Qst_ek = [] ;
ar_ak = [] ;
ar_ek = [] ;
th_ak = [] ;
th_ek = [] ;
tp = [] ;
for ii = 1:length(WTdirs)
    outmatFn = fullfile(outdir, sprintf('WT%03d.mat', ii)) ;
    tmp = load(outmatFn) ;
    Qa_ak = [Qa_ak; tmp.Qa_ak] ;
    Qt_ak = [Qt_ak; tmp.Qt_ak] ;
    Qa_ek = [Qa_ek; tmp.Qa_ek] ;
    Qt_ek = [Qt_ek; tmp.Qt_ek] ;
    Qct_ak = [Qct_ak; tmp.Qct_ak] ;
    Qst_ak = [Qst_ak; tmp.Qst_ak] ;
    Qct_ek = [Qct_ek; tmp.Qct_ek] ;
    Qst_ek = [Qst_ek; tmp.Qst_ek] ;
    ar_ak = [ar_ak; tmp.ar_ak] ;
    ar_ek = [ar_ek; tmp.ar_ek] ;
    th_ak = [th_ak; tmp.th_ak] ;
    th_ek = [th_ek; tmp.th_ek] ;
    tp = [tp, tmp.timestamps ] ;
    
end
weights = 1./ar_ek.^2 ;
[tps_wt, aRmean_wt, aRste_wt, n_wt] = ...
   binDataMeanStdWeighted(tp', ar_ak, timebins, weights, ar_ek) ;
weights = 1./th_ek.^2 ;
[tps_wt, thmean_wt, thste_wt, n_wt] = ...
   binDataMeanStdWeighted(tp', th_ak, timebins, weights, th_ek) ;
weights = 1./Qa_ek.^2 ;
[tps_wt, Qamean_wt, Qaste_wt, n_wt] = ...
   binDataMeanStdWeighted(tp', Qa_ak, timebins, weights, Qa_ek) ;
weights = 1./Qt_ek.^2 ;
[tps_wt, Qtmean_wt, Qtste_wt, n_wt] = ...
   binDataMeanStdWeighted(tp', Qt_ak, timebins, weights, Qt_ek) ;
weights = 1./Qct_ek.^2 ;
[tps_wt, Qcmean_wt, Qcste_wt, n_wt] = ...
   binDataMeanStdWeighted(tp', Qct_ak, timebins, weights, Qct_ek) ;
weights = 1./Qst_ek.^2 ;
[tps_wt, Qsmean_wt, Qsste_wt, n_wt] = ...
   binDataMeanStdWeighted(tp', Qst_ak, timebins, weights, Qst_ek) ;

%OCRL
Qa_ak = [] ;
Qt_ak = [] ;
Qa_ek = [] ;
Qt_ek = [] ;
Qct_ak = [] ;
Qst_ak = [] ;
Qct_ek = [] ;
Qst_ek = [] ;
ar_ak = [] ;
ar_ek = [] ;
th_ak = [] ;
th_ek = [] ;
tp = [] ;
for ii = 1:length(OCRLdirs)
    outmatFn = fullfile(outdir, sprintf('OCRL%03d.mat', ii)) ;
    tmp = load(outmatFn) ;
    Qa_ak = [Qa_ak; tmp.Qa_ak] ;
    Qt_ak = [Qt_ak; tmp.Qt_ak] ;
    Qa_ek = [Qa_ek; tmp.Qa_ek] ;
    Qt_ek = [Qt_ek; tmp.Qt_ek] ;
    Qct_ak = [Qct_ak; tmp.Qct_ak] ;
    Qst_ak = [Qst_ak; tmp.Qst_ak] ;
    Qct_ek = [Qct_ek; tmp.Qct_ek] ;
    Qst_ek = [Qst_ek; tmp.Qst_ek] ;
    ar_ak = [ar_ak; tmp.ar_ak] ;
    ar_ek = [ar_ek; tmp.ar_ek] ;
    th_ak = [th_ak; tmp.th_ak] ;
    th_ek = [th_ek; tmp.th_ek] ;
    tp = [tp, tmp.timestamps ] ;
    
end
weights = 1./ar_ek.^2 ;
[tps_ocrl, aRmean_ocrl, aRste_ocrl, n_ocrl] = ...
   binDataMeanStdWeighted(tp', ar_ak, timebins, weights, ar_ek) ;
weights = 1./th_ek.^2 ;
[tps_ocrl, thmean_ocrl, thste_ocrl, n_ocrl] = ...
   binDataMeanStdWeighted(tp', th_ak, timebins, weights, th_ek) ;
weights = 1./Qa_ek.^2 ;
[tps_ocrl, Qamean_ocrl, Qaste_ocrl, n_ocrl] = ...
   binDataMeanStdWeighted(tp', Qa_ak, timebins, weights, Qa_ek) ;
weights = 1./Qt_ek.^2 ;
[tps_ocrl, Qtmean_ocrl, Qtste_ocrl, n_ocrl] = ...
   binDataMeanStdWeighted(tp', Qt_ak, timebins, weights, Qt_ek) ;
weights = 1./Qct_ek.^2 ;
[tps_ocrl, Qcmean_ocrl, Qcste_ocrl, n_ocrl] = ...
   binDataMeanStdWeighted(tp', Qct_ak, timebins, weights, Qct_ek) ;
weights = 1./Qst_ek.^2 ;
[tps_ocrl, Qsmean_ocrl, Qsste_ocrl, n_ocrl] = ...
   binDataMeanStdWeighted(tp', Qst_ak, timebins, weights, Qst_ek) ;

%ANTP
% timebins = [-0.51, -0.17, 0.16, 0.49, 0.82, 1.16,  1.52]+0.05  ;
% timebins = timebins + 0.
Qa_ak = [] ;
Qt_ak = [] ;
Qa_ek = [] ;
Qt_ek = [] ;
Qct_ak = [] ;
Qst_ak = [] ;
Qct_ek = [] ;
Qst_ek = [] ;
ar_ak = [] ;
ar_ek = [] ;
th_ak = [] ;
th_ek = [] ;
tp = [] ;
for ii = 1:length(ANTPdirs)
    outmatFn = fullfile(outdir, sprintf('ANTP%03d.mat', ii)) ;
    tmp = load(outmatFn) ;
    Qa_ak = [Qa_ak; tmp.Qa_ak] ;
    Qt_ak = [Qt_ak; tmp.Qt_ak] ;
    Qa_ek = [Qa_ek; tmp.Qa_ek] ;
    Qt_ek = [Qt_ek; tmp.Qt_ek] ;
    Qct_ak = [Qct_ak; tmp.Qct_ak] ;
    Qst_ak = [Qst_ak; tmp.Qst_ak] ;
    Qct_ek = [Qct_ek; tmp.Qct_ek] ;
    Qst_ek = [Qst_ek; tmp.Qst_ek] ;
    ar_ak = [ar_ak; tmp.ar_ak] ;
    ar_ek = [ar_ek; tmp.ar_ek] ;
    th_ak = [th_ak; tmp.th_ak] ;
    th_ek = [th_ek; tmp.th_ek] ;
    tp = [tp, tmp.timestamps - 0.05 ] ;
end
weights = 1./ar_ek.^2 ;
[tps_antp, aRmean_antp, aRste_antp, n_antp] = ...
   binDataMeanStdWeighted(tp', ar_ak, timebins, weights, ar_ek) ; 
weights = 1./th_ek.^2 ;
[tps_antp, thmean_antp, thste_antp, n_antp] = ...
   binDataMeanStdWeighted(tp', th_ak, timebins, weights, th_ek) ;
weights = 1./Qa_ek.^2 ;
[tps_antp, Qamean_antp, Qaste_antp, n_antp] = ...
   binDataMeanStdWeighted(tp', Qa_ak, timebins, weights, Qa_ek) ;
weights = 1./Qt_ek.^2 ;
[tps_antp, Qtmean_antp, Qtste_antp, n_antp] = ...
   binDataMeanStdWeighted(tp', Qt_ak, timebins, weights, Qt_ek) ;
weights = 1./Qct_ek.^2 ;
[tps_antp, Qcmean_antp, Qcste_antp, n_antp] = ...
   binDataMeanStdWeighted(tp', Qct_ak, timebins, weights, Qct_ek) ;
weights = 1./Qst_ek.^2 ;
[tps_antp, Qsmean_antp, Qsste_antp, n_antp] = ...
   binDataMeanStdWeighted(tp', Qst_ak, timebins, weights, Qst_ek) ;


% PVALUES at 1 hr
pid = find(tps_wt == 1) ;

% T-test 
% tobs = (x-mu) / (std / sqrt(n)) 
% z score = (xmean1 - xmean2) / sqrt(std1^2 + std2^2)
% assert(sum(~isnan(areaCurvControl)) == nControlExpts)

% % % num = mean(areaCurv_ocrl) - nanmean(areaCurv_wt) ;
% % % denom = sqrt(nanvar(areaCurv_wt) / sum(~isnan(areaCurv_wt)) + ...
% % %     nanvar(areaCurv_ocrl) / sum(~isnan(areaCurv_wt))) ;
% % % zscore = num / denom ;
% % % pvalue_ocrlArea = normcdf(zscore);

% ALternative T-test at single time -- ocrl QA
num = -Qamean_ocrl(pid) + Qamean_wt(pid) ;
denom = sqrt(Qaste_ocrl(pid).^2 + Qaste_wt(pid).^2) ;
zscoreEndVals = num / denom ;
pvalue_ocrlEndVals_Qa = normcdf(zscoreEndVals);

% ALternative T-test at single time -- antp QA
num = -Qamean_antp(pid) + Qamean_wt(pid) ;
denom = sqrt(Qaste_antp(pid).^2 + Qaste_wt(pid).^2) ;
zscoreEndVals = num / denom ;
pvalue_antpEndVals_Qa = normcdf(zscoreEndVals);

% ALternative T-test at single time -- ocrl Qt
num = -Qtmean_ocrl(pid) + Qtmean_wt(pid) ;
denom = sqrt(Qtste_ocrl(pid).^2 + Qtste_wt(pid).^2) ;
zscoreEndVals = num / denom ;
pvalue_ocrlEndVals_Qt = normcdf(zscoreEndVals);

% ALternative T-test at single time -- antp Qt
num = -Qtmean_antp(pid) + Qtmean_wt(pid) ;
denom = sqrt(Qtste_antp(pid).^2 + Qtste_wt(pid).^2) ;
zscoreEndVals = num / denom ;
pvalue_antpEndVals_Qt = normcdf(zscoreEndVals);


% ALternative T-test at single time -- ocrl Qc
num = Qcmean_ocrl(pid) - Qcmean_wt(pid) ;
denom = sqrt(Qcste_ocrl(pid).^2 + Qcste_wt(pid).^2) ;
zscoreEndVals = num / denom ;
pvalue_ocrlEndVals_Qc = normcdf(zscoreEndVals);

% ALternative T-test at single time -- antp Qc
num = Qcmean_antp(pid) - Qcmean_wt(pid) ;
denom = sqrt(Qcste_antp(pid).^2 + Qcste_wt(pid).^2) ;
zscoreEndVals = num / denom ;
pvalue_antpEndVals_Qc = normcdf(zscoreEndVals);


% ALternative T-test at single time -- ocrl Qs
num = Qsmean_ocrl(pid) - Qsmean_wt(pid) ;
denom = sqrt(Qsste_ocrl(pid).^2 + Qsste_wt(pid).^2) ;
zscoreEndVals = num / denom ;
pvalue_ocrlEndVals_Qs = normcdf(zscoreEndVals);

% ALternative T-test at single time -- antp Qs
num = Qsmean_antp(pid) - Qsmean_wt(pid) ;
denom = sqrt(Qsste_antp(pid).^2 + Qsste_wt(pid).^2) ;
zscoreEndVals = num / denom ;
pvalue_antpEndVals_Qs = normcdf(zscoreEndVals);


% ALternative T-test at single time -- ocrl aR
num = -aRmean_ocrl(pid) + aRmean_wt(pid) ;
denom = sqrt(aRste_ocrl(pid).^2 + aRste_wt(pid).^2) ;
zscoreEndVals = num / denom ;
pvalue_ocrlEndVals_aR = normcdf(zscoreEndVals);

% ALternative T-test at single time -- antp aR
num = -aRmean_antp(pid) + aRmean_wt(pid) ;
denom = sqrt(aRste_antp(pid).^2 + aRste_wt(pid).^2) ;
zscoreEndVals = num / denom ;
pvalue_antpEndVals_aR = normcdf(zscoreEndVals);

% ALternative T-test at single time -- ocrl theta
num = -thmean_ocrl(pid) + thmean_wt(pid) ;
denom = sqrt(thste_ocrl(pid).^2 + thste_wt(pid).^2) ;
zscoreEndVals = num / denom ;
pvalue_ocrlEndVals_th = normcdf(zscoreEndVals);

% ALternative T-test at single time -- antp theta
num = -thmean_antp(pid) + thmean_wt(pid) ;
denom = sqrt(thste_antp(pid).^2 + thste_wt(pid).^2) ;
zscoreEndVals = num / denom ;
pvalue_antpEndVals_th = normcdf(zscoreEndVals);

% PLOT it
close all
hf = figure('Position', [100 100 320 600], 'units', 'centimeters');
hold on;

% Aspect ratio (raw)
subplot(3, 2, 5); hold on;
lineProps = {'.-', 'color', colors(1, :)} ;
h1 = shadedErrorBar(tps_wt, aRmean_wt, ...
    movmean(aRste_wt, 3, 'omitnan'), 'lineProps', lineProps) ;
lineProps = {'.-', 'color', colors(2, :)} ;
h2 = shadedErrorBar(tps_ocrl, aRmean_ocrl, ...
     movmean(aRste_ocrl, 3, 'omitnan'), 'lineProps', lineProps) ;
lineProps = {'.-', 'color', colors(3, :)} ;
h3 = shadedErrorBar(tps_antp, aRmean_antp, ...
     movmean(aRste_antp, 3, 'omitnan'), 'lineProps', lineProps) ;
xlim(xlims)
ylabel('$\langle a/b \rangle$ (raw)', 'interpreter', 'latex')
title(['$z = ' num2str(pvalue_ocrlEndVals_aR) ',' ...
    num2str(pvalue_antpEndVals_aR) '$'], ...
    'interpreter', 'latex')

% Theta (raw)
subplot(3, 2, 6); hold on;
lineProps = {'.-', 'color', colors(1, :)} ;
h1 = shadedErrorBar(tps_wt, thmean_wt, ...
    movmean(thste_wt, 3, 'omitnan'), 'lineProps', lineProps) ;
lineProps = {'.-', 'color', colors(2, :)} ;
h2 = shadedErrorBar(tps_ocrl, thmean_ocrl, ...
    movmean(thste_ocrl, 3, 'omitnan'), 'lineProps', lineProps) ;
lineProps = {'.-', 'color', colors(3, :)} ;
h3 = shadedErrorBar(tps_antp, thmean_antp, ...
    movmean(thste_antp, 3, 'omitnan'), 'lineProps', lineProps) ;
xlim(xlims)
yticks([0, pi/2, pi])
set(gca, 'yticklabels', {'0', '\pi/2', '\pi'})
ylim([0,pi])
ylabel('$\langle \theta \rangle$ (raw)', 'interpreter', 'latex')
title(['$z = ' num2str(pvalue_ocrlEndVals_th) ',' ...
    num2str(pvalue_antpEndVals_th) '$'], ...
    'interpreter', 'latex')

% Qaspect (averaged shapes)
subplot(3, 2, 3); hold on;
lineProps = {'.-', 'color', colors(1, :)} ;
h1 = shadedErrorBar(tps_wt, Qamean_wt, movmean(Qaste_wt, 3, 'omitnan'), ...
    'lineProps', lineProps) ;
lineProps = {'.-', 'color', colors(2, :)} ;
h2 = shadedErrorBar(tps_ocrl, Qamean_ocrl, ...
    movmean(Qaste_ocrl, 3, 'omitnan'), ...
    'lineProps', lineProps) ;
lineProps = {'.-', 'color', colors(3, :)} ;
h3 = shadedErrorBar(tps_antp, Qamean_antp, ...
    movmean(Qaste_antp, 3, 'omitnan'),...
    'lineProps', lineProps) ;
xlim(xlims)
ylabel('$(a/b)_{\langle Q \rangle}$', 'interpreter', 'latex')
title(['$z = ' num2str(pvalue_ocrlEndVals_Qa) ',' ...
    num2str(pvalue_antpEndVals_Qa) '$'], ...
    'interpreter', 'latex')

% Qtheta (averaged shapes)
subplot(3, 2, 4); hold on;
lineProps = {'.-', 'color', colors(1, :)} ;
h1 = shadedErrorBar(tps_wt, Qtmean_wt, movmean(Qtste_wt, 3, 'omitnan'), ...
    'lineProps', lineProps) ;
lineProps = {'.-', 'color', colors(2, :)} ;
h2 = shadedErrorBar(tps_ocrl, Qtmean_ocrl, ...
    movmean(Qtste_ocrl, 3, 'omitnan'), ...
    'lineProps', lineProps) ;
lineProps = {'.-', 'color', colors(3, :)} ;
h3 = shadedErrorBar(tps_antp, Qtmean_antp, ...
    movmean(Qtste_antp, 3, 'omitnan'),...
    'lineProps', lineProps) ;
xlim(xlims)
ylim([0, pi])
yticks([0, pi/2, pi])
set(gca, 'yticklabels', {'0', '\pi/2', '\pi'}) 
ylabel('$\theta_{\langle Q \rangle}$',...
    'interpreter', 'latex')
title(['$z = ' num2str(pvalue_ocrlEndVals_Qt) ',' ...
    num2str(pvalue_antpEndVals_Qt) '$'], ...
    'interpreter', 'latex')

% Qcos2theta (averaged shapes)
subplot(3, 2, 1); hold on;
lineProps = {'.-', 'color', colors(1, :)} ;
h1 = shadedErrorBar(tps_wt, -Qcmean_wt+1, movmean(Qcste_wt, 3, 'omitnan'), ...
    'lineProps', lineProps) ;
lineProps = {'.-', 'color', colors(2, :)} ;
h2 = shadedErrorBar(tps_ocrl, -Qcmean_ocrl+1, movmean(Qcste_ocrl, 3, 'omitnan'), ...
    'lineProps', lineProps) ;
lineProps = {'.-', 'color', colors(3, :)} ;
h3 = shadedErrorBar(tps_antp, -Qcmean_antp+1, movmean(Qcste_antp, 3, 'omitnan'),...
    'lineProps', lineProps) ;
xlim(xlims)
ylabel('$1-\left(a/b-1 \right)\cos 2\theta$', 'interpreter', 'latex')
title(['$z = ' num2str(pvalue_ocrlEndVals_Qc) ',' ...
    num2str(pvalue_antpEndVals_Qc) '$'], ...
    'interpreter', 'latex')

% Qsin2theta (averaged shapes)
subplot(3, 2, 2); hold on;
lineProps = {'.-', 'color', colors(1, :)} ;
h1 = shadedErrorBar(tps_wt, Qsmean_wt, movmean(Qsste_wt, 3, 'omitnan'), ...
    'lineProps', lineProps) ;
lineProps = {'.-', 'color', colors(2, :)} ;
h2 = shadedErrorBar(tps_ocrl, Qsmean_ocrl, movmean(Qsste_ocrl,3, 'omitnan'), ...
    'lineProps', lineProps) ;
lineProps = {'.-', 'color', colors(3, :)} ;
h3 = shadedErrorBar(tps_antp, Qsmean_antp, movmean(Qsste_antp, 3, 'omitnan'), ...
    'lineProps', lineProps) ;
xlim(xlims)
ylabel('$\left(a/b-1 \right)\sin 2\theta$', 'interpreter', 'latex')
title(['$z = ' num2str(pvalue_ocrlEndVals_Qs) ',' ...
    num2str(pvalue_antpEndVals_Qs) '$'], ...
    'interpreter', 'latex')

sgh = sgtitle('$p$ values for OCRL, Antp at 1 hr', 'interpreter', 'latex') ;
outfn = fullfile(outdir, 'stats_results_binned_OCRL_WT_ANTP.pdf') ;
disp(['Saving ' outfn])
saveas(gcf, outfn)

%% Get p value from all timepoints under curv from t=(0, 1) hr
% confocal_segmentation_collate_helper_integrateCurv()
assert(all(tps_wt == tps_ocrl))
assert(all(tps_wt == tps_antp))

% Qxx
tid = tps_wt > 0.2 & tps_wt < 1.01 ;
ntid = length(tid) ;
diffQc_ocrl = Qcmean_ocrl(tid) - Qcmean_wt(tid) ;
diffQs_ocrl = Qsmean_ocrl(tid) - Qsmean_wt(tid) ;
diffQa_ocrl = - Qamean_ocrl(tid) + Qamean_wt(tid) ;
diffQt_ocrl = Qtmean_ocrl(tid) - Qtmean_wt(tid) ;
diffaR_ocrl = - aRmean_ocrl(tid) + aRmean_wt(tid) ;
diffth_ocrl = thmean_ocrl(tid) - thmean_wt(tid) ;
denomQc_ocrl = sqrt( Qcste_ocrl(tid).^2 + Qcste_wt(tid).^2 ) ;
denomQs_ocrl = sqrt( Qsste_ocrl(tid).^2 + Qsste_wt(tid).^2 ) ;
denomQa_ocrl = sqrt( Qaste_ocrl(tid).^2 + Qaste_wt(tid).^2 ) ;
denomQt_ocrl = sqrt( Qtste_ocrl(tid).^2 + Qtste_wt(tid).^2 ) ;
denomaR_ocrl = sqrt( aRste_ocrl(tid).^2 + aRste_wt(tid).^2 ) ;
denomth_ocrl = sqrt( thste_ocrl(tid).^2 + thste_wt(tid).^2 ) ;
zQc_ocrl = mean( diffQc_ocrl ./ denomQc_ocrl) * sqrt(ntid) ;
zQs_ocrl = mean( diffQs_ocrl ./ denomQs_ocrl) * sqrt(ntid) ;
zQa_ocrl = mean( diffQa_ocrl ./ denomQa_ocrl) * sqrt(ntid) ;
zQt_ocrl = mean( diffQt_ocrl ./ denomQt_ocrl) * sqrt(ntid) ;
zaR_ocrl = mean( diffaR_ocrl ./ denomaR_ocrl) * sqrt(ntid) ;
zth_ocrl = mean( diffth_ocrl ./ denomth_ocrl) * sqrt(ntid) ;
pQc_ocrl = normcdf(zQc_ocrl) ;
pQs_ocrl = normcdf(zQs_ocrl) ;
pQa_ocrl = normcdf(zQa_ocrl) ;
pQt_ocrl = normcdf(zQt_ocrl) ;
paR_ocrl = normcdf(zaR_ocrl) ;
pth_ocrl = normcdf(zth_ocrl) ;

diffQc_antp = Qcmean_antp(tid) - Qcmean_wt(tid) ;
diffQs_antp = Qsmean_antp(tid) - Qsmean_wt(tid) ;
diffQa_antp = -Qamean_antp(tid) + Qamean_wt(tid) ;
diffQt_antp = Qtmean_antp(tid) - Qtmean_wt(tid) ;
diffaR_antp = -aRmean_antp(tid) + aRmean_wt(tid) ;
diffth_antp = thmean_antp(tid) - thmean_wt(tid) ;
denomQc_antp = sqrt( Qcste_antp(tid).^2 + Qcste_wt(tid).^2 ) ;
denomQs_antp = sqrt( Qsste_antp(tid).^2 + Qsste_wt(tid).^2 ) ;
denomQa_antp = sqrt( Qaste_antp(tid).^2 + Qaste_wt(tid).^2 ) ;
denomQt_antp = sqrt( Qtste_antp(tid).^2 + Qtste_wt(tid).^2 ) ;
denomaR_antp = sqrt( aRste_antp(tid).^2 + aRste_wt(tid).^2 ) ;
denomth_antp = sqrt( thste_antp(tid).^2 + thste_wt(tid).^2 ) ;
zQc_antp = mean(diffQc_antp ./ denomQc_antp) * sqrt(ntid) ;
zQs_antp = mean(diffQs_antp ./ denomQs_antp) * sqrt(ntid) ;
zQa_antp = mean(diffQa_antp ./ denomQa_antp) * sqrt(ntid) ;
zQt_antp = mean(diffQt_antp ./ denomQt_antp) * sqrt(ntid) ;
zaR_antp = mean(diffaR_antp ./ denomaR_antp) * sqrt(ntid) ;
zth_antp = mean(diffth_antp ./ denomth_antp) * sqrt(ntid) ;
pQc_antp = normcdf(zQc_antp) ;
pQs_antp = normcdf(zQs_antp) ;
pQa_antp = normcdf(zQa_antp) ;
pQt_antp = normcdf(zQt_antp) ;
paR_antp = normcdf(zaR_antp) ;
pth_antp = normcdf(zth_antp) ;


subplot(3, 2, 1) ;
title(['$z = ' num2str(pQc_ocrl) ',' ...
    num2str(pQc_antp) '$'], ...
    'interpreter', 'latex')
subplot(3, 2, 2) ;
title(['$z = ' num2str(pQs_ocrl) ',' ...
    num2str(pQs_antp) '$'], ...
    'interpreter', 'latex')
subplot(3, 2, 3) ;
title(['$z = ' num2str(pQa_ocrl) ',' ...
    num2str(pQa_antp) '$'], ...
    'interpreter', 'latex')
subplot(3, 2, 4) ;
title(['$z = ' num2str(pQt_ocrl) ',' ...
    num2str(pQt_antp) '$'], ...
    'interpreter', 'latex')
subplot(3, 2, 5) ;
title(['$z = ' num2str(paR_ocrl) ',' ...
    num2str(paR_antp) '$'], ...
    'interpreter', 'latex')
subplot(3, 2, 6) ;
title(['$z = ' num2str(pth_ocrl) ',' ...
    num2str(pth_antp) '$'], ...
    'interpreter', 'latex')
sgh.String = '$p$ values for OCRL, Antp from t=[20,60] min' ;
outfn = fullfile(outdir, 'stats_results_binned_OCRL_WT_ANTP_fullPVal.pdf') ;
disp(['Saving ' outfn])
saveas(gcf, outfn) 


%% 
error('here')
    % QQ
    in12 = cx < OCRLmidX(ii) ;
    nCells = length(find(in12)) ;
    mratio = m2(in12) ./ m1(in12) ;
    strength = zeros(nCells, 1) ;
    eccentricity = zeros(nCells, 1) ;
    QQ = zeros(nCells, 2, 2) ;
    for qq = 1:nCells
        if ~isempty(intersect(keep, qq))
            tt = mod(seg3d.qualities.ang1(qq), pi) ;
            nn = [cos(tt), sin(tt)] ;
            % Create traceless symmetric matrix using unit vec
            strength(qq) = abs(sqrt(mratio(qq))) - 1 ;
            QQ(qq, :, :) = nn' * nn - 0.5 * [1, 0; 0, 1] ;
        end
    end
    [meanQA, stdQA, steQA] = QtensorStats(aa .* QQ(keep, :, :), aa) ;
    [meanQAAspect, meanQATheta, ...
        meanQAAspectStd, meanQAThetaStd] = ...
        QtensorAspectTheta(meanQA, stdQA) ;
    [~,~, meanQAAspectSte, meanQAThetaSte] = ...
        QtensorAspectTheta(meanQA, steQA) ;


%% Qxx
for shadingStyle = 1:2
    close all
    % fig = figure('units', 'cm', 'outerposition', [0.5 0 0.5 1]);
    
    hf = figure('Position', [100 100 320 120], 'units', 'centimeters');
    for ii = 1:length(OCRLdirs)
        ocrldir = OCRLdirs{ii} ;
        load(fullfile(ocrldir, 'results.mat'), ...
            'Qct_ak', 'Qct_sk', 'Qct_ek', ...
            'Qst_ak', 'Qst_sk', 'Qst_ek', 'timestamps', 'exciteIdx') ;
        if contains(lower(timeUnits), 'hr') || contains(lower(timeUnits), 'hour')
            timestamps = timestamps /60 ;
        end
        
        % Qxx
        subplot(1, 2, 1)
        hold on;
        if shadingStyle == 1
            lineProps = {'.-', 'color', colors(1, :)} ;
            h1=shadedErrorBar(timestamps, Qct_ak, Qct_sk, 'lineProps', lineProps) ;
            plot(timestamps, Qct_ak-Qct_ek, '--', 'color', colors(1, :)) ;
            plot(timestamps, Qct_ak+Qct_ek, '--', 'color', colors(1, :)) ;
            pbaspect([1, 1.5, 1])
        else
            lineProps = {'.-', 'color', colors(1, :)} ;
            h1=shadedErrorBar(timestamps, Qct_ak, Qct_ek, 'lineProps', lineProps) ;
            pbaspect([1, 1, 1])
        end
        xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
        ylabel('$Q_{xx}$', 'interpreter', 'latex')

        % Qyy
        subplot(1, 2, 2)
        hold on;
        if shadingStyle == 1
            lineProps = {'.-', 'color', colors(1, :)} ;
            h1=shadedErrorBar(timestamps, Qst_ak, Qct_sk, 'lineProps', lineProps) ;
            plot(timestamps, Qst_ak-Qst_ek, '--', 'color', colors(1, :)) ;
            plot(timestamps, Qst_ak+Qst_ek, '--', 'color', colors(1, :)) ;
            pbaspect([1, 1.5, 1])
        else
            % h1=plot(timestamps, Qst_ak, '.-', 'color', colors(1, :)) ;
            lineProps = {'.-', 'color', colors(1, :)} ;
            h1=shadedErrorBar(timestamps, Qst_ak, Qct_ek, 'lineProps', lineProps) ;
            pbaspect([1, 1, 1])
        end
        xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
        ylabel('$Q_{yy}$', 'interpreter', 'latex')
    end

    % Save figure
    saveas(gcf, fullfile(outdir, sprintf('Q_results_OCRL_style%d.png', shadingStyle)))
    saveas(gcf, fullfile(outdir, sprintf('Q_results_OCRL_style%d.pdf', shadingStyle)))
    
    

    %% Compare to WT

    for ii = 1:length(WTdirs)
        wtdir = WTdirs{ii} ;
        load(fullfile(wtdir, 'results.mat'), ...
            'Qct_ak', 'Qct_sk', 'Qct_ek', ...
            'Qst_ak', 'Qst_sk', 'Qst_ek', 'timestamps', 'exciteIdx') ;

        if contains(lower(timeUnits), 'hr') || contains(lower(timeUnits), 'hour')
            timestamps = timestamps /60 ;
        end
        % Qxx
        subplot(1, 2, 1)
        hold on;
        if shadingStyle == 1
            lineProps = {'.-', 'color', colors(2, :)} ;
            h1=shadedErrorBar(timestamps, Qct_ak, Qct_sk, 'lineProps', lineProps) ;
            plot(timestamps, Qct_ak-Qct_ek, '--', 'color', colors(2, :)) ;
            plot(timestamps, Qct_ak+Qct_ek, '--', 'color', colors(2, :)) ;
            pbaspect([1, 1.5, 1])
        else
            % h1=plot(timestamps, Qct_ak, '.-', 'color', colors(2, :)) ;
            lineProps = {'.-', 'color', colors(2, :)} ;
            h1=shadedErrorBar(timestamps, Qct_ak, Qct_ek, 'lineProps', lineProps) ;
            pbaspect([1, 1, 1])
        end
        xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
        ylabel('$Q_{xx}$', 'interpreter', 'latex')

        % Qyy
        subplot(1, 2, 2)
        hold on;
        if shadingStyle == 1
            lineProps = {'.-', 'color', colors(2, :)} ;
            h1=shadedErrorBar(timestamps, Qst_ak, Qct_sk, 'lineProps', lineProps) ;
            plot(timestamps, Qst_ak-Qst_ek, '--', 'color', colors(2, :)) ;
            plot(timestamps, Qst_ak+Qst_ek, '--', 'color', colors(2, :)) ;
            pbaspect([1, 1.5, 1])
        else
            % h1=plot(timestamps, Qst_ak, '.-', 'color', colors(2, :)) ;
            lineProps = {'.-', 'color', colors(2, :)} ;
            h1=shadedErrorBar(timestamps, Qst_ak, Qst_ek, 'lineProps', lineProps) ;
            pbaspect([1, 1, 1])
        end
        xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
        ylabel('$Q_{yy}$', 'interpreter', 'latex')
    end

    % Save figure
    saveas(gcf, fullfile(outdir, sprintf('Q_results_OCRL_WT_confocal_style%d.png', shadingStyle)))
    saveas(gcf, fullfile(outdir, sprintf('Q_results_OCRL_WT_confocal_style%d.pdf', shadingStyle)))
    
    
    
    %% Compare to ANTP
    for ii = 1:length(ANTPdirs)
        Antpdir = ANTPdirs{ii} ;
        load(fullfile(Antpdir, 'results.mat'), ...
            'Qct_ak', 'Qct_sk', 'Qct_ek', ...
            'Qst_ak', 'Qst_sk', 'Qst_ek', 'timestamps', 'exciteIdx') ;

        if contains(lower(timeUnits), 'hr') || contains(lower(timeUnits), 'hour')
            timestamps = timestamps /60 ;
        end
        % Qxx
        subplot(1, 2, 1)
        hold on;
        if shadingStyle == 1
            lineProps = {'.-', 'color', colors(3, :)} ;
            h1=shadedErrorBar(timestamps, Qct_ak, Qct_sk, 'lineProps', lineProps) ;
            plot(timestamps, Qct_ak-Qct_ek, '--', 'color', colors(3, :)) ;
            plot(timestamps, Qct_ak+Qct_ek, '--', 'color', colors(3, :)) ;
            pbaspect([1, 1.5, 1])
        else
            % h1=plot(timestamps, Qct_ak, '.-', 'color', colors(2, :)) ;
            lineProps = {'.-', 'color', colors(3, :)} ;
            h1=shadedErrorBar(timestamps, Qct_ak, Qct_ek, 'lineProps', lineProps) ;
            pbaspect([1, 1, 1])
        end
        xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
        ylabel('$Q_{xx}$', 'interpreter', 'latex')

        % Qyy
        subplot(1, 2, 2)
        hold on;
        if shadingStyle == 1
            lineProps = {'.-', 'color', colors(3, :)} ;
            h1=shadedErrorBar(timestamps, Qst_ak, Qct_sk, 'lineProps', lineProps) ;
            plot(timestamps, Qst_ak-Qst_ek, '--', 'color', colors(3, :)) ;
            plot(timestamps, Qst_ak+Qst_ek, '--', 'color', colors(3, :)) ;
            pbaspect([1, 1.5, 1])
        else
            % h1=plot(timestamps, Qst_ak, '.-', 'color', colors(2, :)) ;
            lineProps = {'.-', 'color', colors(3, :)} ;
            h1=shadedErrorBar(timestamps, Qst_ak, Qst_ek, 'lineProps', lineProps) ;
            pbaspect([1, 1, 1])
        end
        xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
        ylabel('$Q_{yy}$', 'interpreter', 'latex')
    end

    % Save figure
    saveas(gcf, fullfile(outdir, sprintf('Q_results_OCRL_WT_ANTP_confocal_style%d.png', shadingStyle)))
    saveas(gcf, fullfile(outdir, sprintf('Q_results_OCRL_WT_ANTP_confocal_style%d.pdf', shadingStyle)))
    
    
    %% Compare to WT Lightsheet
    % % WTfn =  '/mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/msls_output/cellSegmentation/seg3d_corrected/stats_summary.mat' ;
    WTfn =  '/mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/msls_output/cellSegmentation/seg3d_corrected/stats_summary_L12.mat' ;
    tmp = load(WTfn) ;
    
    if contains(lower(timeUnits), 'hr') || contains(lower(timeUnits), 'hour')
        tmp.timeStamps = tmp.timeStamps /60 ;
    end
    
    % Qxx
    subplot(1, 2, 1)
    lineProps = {'.-', 'color', colors(2, :)} ;
    h2 = shadedErrorBar(tmp.timeStamps, 0.5 * tmp.meanxsL12, ...
        0.5 * tmp.stdxs, 'lineProps', lineProps) ;
    plot(tmp.timeStamps, 0.5 * (tmp.meanxsL12 + tmp.stdmeanxs), '--', 'color', colors(2, :)) ;
    plot(tmp.timeStamps, 0.5 * (tmp.meanxsL12 - tmp.stdmeanxs), '--', 'color', colors(2, :)) ;
    xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
    ylabel('$Q_{xx}$', 'interpreter', 'latex')
    
    % Qyy
    subplot(1, 2, 2)
    lineProps = {'.-', 'color', colors(2, :)} ;
    h2 = shadedErrorBar(tmp.timeStamps, 0.5 * tmp.meanysL12, ...
        0.5 * tmp.stdys, 'lineProps', lineProps) ;
    plot(tmp.timeStamps, 0.5 * (tmp.meanysL12 + tmp.stdmeanys), '--', 'color', colors(2, :)) ;
    plot(tmp.timeStamps, 0.5 * (tmp.meanysL12 - tmp.stdmeanys), '--', 'color', colors(2, :)) ;
    xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
    ylabel('$Q_{yy}$', 'interpreter', 'latex')
    
    % Save figure
    saveas(gcf, fullfile(outdir, 'Q_results_OCRL_WT_ANTP_wLightsheet.png'))
    saveas(gcf, fullfile(outdir, 'Q_results_OCRL_WT_ANTP_wLightsheet.pdf'))

end

