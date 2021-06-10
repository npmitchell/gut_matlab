%% Collate confocal_segmentation_antpOCRL
% npmitchell 2021
clearvars
addpath('/mnt/data/code/gut_matlab/addpath_recurse/')
addpath_recurse('/mnt/data/code/gut_matlab/')
colors = define_colors ;

timeUnits = 'min' ;
outdir = '/mnt/data/optogenetics_confocal/antpGAL4/' ;
OCRLdirs = {'/mnt/data/optogenetics_confocal/antpGAL4/huygens_deconvolution_noKlar/202103181730_antpG4OCRLGap43mCh_40x1p6x_5mpf_4pc3pc_to_12pc9pc_600ns_lav3_DC/', ...
    '/mnt/data/optogenetics_confocal/antpGAL4/huygens_deconvolution_withKlar/202105252111_AntpG4kOCRLgap43_0p75um_1p25x40x_lav3_3t6pc3t6pc_5mpf_480ns_LED4/',...
    '/mnt/data/optogenetics_confocal/antpGAL4/huygens_deconvolution_withKlar/202105311838_AntpG4kOCRL_0p75um_1p5x40_3t6pc_lav3_86s_5mpf/',...
    } ;

close all
for ii = 1:length(OCRLdirs)
    ocrldir = OCRLdirs{ii} ;
    load(fullfile(ocrldir, 'cellSegmentation/results.mat'), ...
        'Qct_ak', 'Qct_sk', 'Qct_ek', ...
        'Qst_ak', 'Qst_sk', 'Qst_ek', 'timestamps', 'exciteIdx') ;

    % Qxx
    subplot(1, 2, 1)
    hold on;
    lineProps = {'.-', 'color', colors(1, :)} ;
    h1=shadedErrorBar(timestamps, Qct_ak, Qct_sk, 'lineProps', lineProps) ;
    plot(timestamps, Qct_ak-Qct_ek, '--', 'color', colors(1, :)) ;
    plot(timestamps, Qct_ak+Qct_ek, '--', 'color', colors(1, :)) ;
    xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
    ylabel('$Q_{xx}$', 'interpreter', 'latex')

    % Qyy
    subplot(1, 2, 2)
    hold on;
    lineProps = {'.-', 'color', colors(1, :)} ;
    h1=shadedErrorBar(timestamps, Qst_ak, Qst_sk, 'lineProps', lineProps) ;
    plot(timestamps, Qst_ak-Qst_ek, '--', 'color', colors(1, :)) ;
    plot(timestamps, Qst_ak+Qst_ek, '--', 'color', colors(1, :)) ;
    xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
    ylabel('$Q_{yy}$', 'interpreter', 'latex')
end

% Save figure
saveas(gcf, fullfile(outdir, 'Q_results_OCRL.png'))
saveas(gcf, fullfile(outdir, 'Q_results_OCRL.pdf'))


% clf
% hold on;
% lineProps = {'.-', 'color', colors(1, :)} ;
% h1=shadedErrorBar(timestamps, Qct_ak, Qct_sk, 'lineProps', lineProps) ;
% lineProps = {'.-', 'color', colors(2, :)} ;
% h2=shadedErrorBar(timestamps, Qst_ak, Qst_sk, 'lineProps', lineProps) ;
% 
% plot(timestamps, Qct_ak-Qct_ek, '--', 'color', colors(1, :)) ;
% plot(timestamps, Qct_ak+Qct_ek, '--', 'color', colors(1, :)) ;
% plot(timestamps, Qst_ak-Qst_ek, '--', 'color', colors(2, :)) ;
% plot(timestamps, Qst_ak+Qst_ek, '--', 'color', colors(2, :)) ;
% 
% legend([h1.patch, h2.patch], {'$Q_{xx}$', '$Q_{yy}$'}, 'interpreter', 'latex')
% 
% xlabel('time [min]', 'interpreter', 'latex')
% ylabel('$Q_{xx}$, $Q_{yy}$', 'interpreter', 'latex')
% imfn = fullfile(segDir, 'stats_dynamics') ;
% disp(['Saving statistics (x,t): ' imfn])
% saveas(gcf, [imfn '.pdf']) ;
% saveas(gcf, [imfn '.png']) ;


%% Compare to WT
% WTfn =  '/mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/msls_output/cellSegmentation/seg3d_corrected/stats_summary.mat' ;
WTfn =  '/mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/msls_output/cellSegmentation/seg3d_corrected/stats_summary_L12.mat' ;
tmp = load(WTfn) ;

% Qxx
subplot(1, 2, 1)
lineProps = {'.-', 'color', colors(2, :)} ;
h2 = shadedErrorBar(tmp.timeStamps, tmp.meanxsL12, tmp.stdxs, 'lineProps', lineProps) ;
plot(tmp.timeStamps, tmp.meanxsL12 + tmp.stdmeanxs, '--', 'color', colors(2, :)) ;
plot(tmp.timeStamps, tmp.meanxsL12 - tmp.stdmeanxs, '--', 'color', colors(2, :)) ;
xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
ylabel('$Q_{xx}$', 'interpreter', 'latex')

% Qyy
subplot(1, 2, 2)
lineProps = {'.-', 'color', colors(2, :)} ;
h2 = shadedErrorBar(tmp.timeStamps, tmp.meanysL12, tmp.stdys, 'lineProps', lineProps) ;
plot(tmp.timeStamps, tmp.meanysL12 + tmp.stdmeanys, '--', 'color', colors(2, :)) ;
plot(tmp.timeStamps, tmp.meanysL12 - tmp.stdmeanys, '--', 'color', colors(2, :)) ;
xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
ylabel('$Q_{yy}$', 'interpreter', 'latex')

% Save figure
saveas(gcf, fullfile(outdir, 'Q_results_OCRL_WT.png'))
saveas(gcf, fullfile(outdir, 'Q_results_OCRL_WT.pdf'))

