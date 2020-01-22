% Get autocorrelation in velocities
        disp('Obtaining autocorrelation in velocities...')
        acorr = zeros(size(vsmM, 2), size(vsmM, 3), 21) ;
        for j = 1:size(vsmM, 2)
            for k = 1:size(vsmM, 3)
                acorr(j, k, :) = autocorr(squeeze(vM(:, j, k))) ;
            end
        end
        mean_acorr = squeeze(mean(acorr, 1)) ;
        std_acorr = squeeze(std(acorr, 1)) ;
        % plot autocorrelation
        close all
        for nn=1:3
            errorbar(1:21, mean_acorr(nn, :), std_acorr(nn,:))
            hold on;
        end
        legend({'v_x', 'v_y', 'v_z'})
        title('raw correlations in velocities')
        xlabel('dt [min]')
        ylabel('correlation')
        saveas(gcf, fullfile(pivSimpleAvgDir, 'autocorr_velocities.png'))

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get autocorrelation in smoothed velocities
        disp('Obtaining autocorrelation in smoothed velocities...')
        acorr = zeros(size(vsmM, 2), size(vsmM, 3), 21) ;
        for j = 1:size(vsmM, 2)
            for k = 1:size(vsmM, 3)
                acorr(j, k, :) = autocorr(squeeze(vsmM(:, j, k))) ;
            end
        end
        mean_acorr = squeeze(mean(acorr, 1)) ;
        std_acorr = squeeze(std(acorr, 1)) ;
        % plot autocorrelation
        close all
        for nn=1:3
            errorbar(1:21, mean_acorr(nn, :), std_acorr(nn,:))
            hold on;
        end
        legend({'v_x', 'v_y', 'v_z'})
        title('raw correlations in velocities')
        xlabel('dt [min]')
        ylabel('correlation')
        saveas(gcf, fullfile(pivSimpleAvgDir, 'autocorr_smoothed_velocities.png'))

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get autocorrelation in normal vectors
        disp('Obtaining autocorrelation...')
        acorr = zeros(size(vsmM, 2), size(vsmM, 3), 21) ;
        for j = 1:size(vsmM, 2)
            for k = 1:size(vsmM, 3)
                acorr(j, k, :) = autocorr(squeeze(nM(:, j, k))) ;
            end
        end
        mean_acorr = squeeze(mean(acorr, 1)) ;
        std_acorr = squeeze(std(acorr, 1)) ;
        % plot autocorrelation
        close all
        for nn=1:3
            errorbar(1:21, mean_acorr(nn, :), std_acorr(nn,:))
            hold on;
        end
        legend({'v_x', 'v_y', 'v_z'})
        title('raw correlations in normal vectors')
        xlabel('dt [min]')
        ylabel('correlation')
        saveas(gcf, fullfile(pivSimpleAvgDir, 'autocorr_normals.png'))
        close all