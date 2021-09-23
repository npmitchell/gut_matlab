%% Collate confocal_segmentation_antpOCRL
% npmitchell 2021
clearvars
addpath('/mnt/data/code/gut_matlab/addpath_recurse/')
addpath_recurse('/mnt/data/code/gut_matlab/')
colors = define_colors ;

timeUnits = 'hr' ;
outdir = '/mnt/data/optogenetics_confocal/antpGAL4/' ;
OCRLdirs = {'/mnt/data/optogenetics_confocal/antpGAL4/huygens_deconvolution_noKlar/202103181730_antpG4OCRLGap43mCh_40x1p6x_5mpf_4pc3pc_to_12pc9pc_600ns_lav3_DC/cellSegmentation/', ...
    '/mnt/data/optogenetics_confocal/antpGAL4/huygens_deconvolution_withKlar/202105252111_AntpG4kOCRLgap43_0p75um_1p25x40x_lav3_3t6pc3t6pc_5mpf_480ns_LED4/cellSegmentation/ground_truth/',...
    '/mnt/data/optogenetics_confocal/antpGAL4/huygens_deconvolution_withKlar/202105311838_AntpG4kOCRL_0p75um_1p5x40_3t6pc_lav3_86s_5mpf/cellSegmentation/',...
  } ;

WTdirs = {'/mnt/data/optogenetics_confocal/WTcontrol_48YG4kCAAXmCh/202106081815_48YG4kCAAXmCh_0p75um_1p5x40x_3t6pc3t6pc_lav3_615ns_5mpf/cellSegmentation/'} ;

ANTPdirs = {'/mnt/data/confocal_data/gut/AntpNSRC3_48YG4CAAXmCh/202108291112_antpNSRC3_48YG4CAAXmCh_2um_5pc_10mpf_lav3fac2_e1e2e3e5mutant/e1/cellSegmentation/'};

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

