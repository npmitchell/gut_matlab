cd /mnt/data/montell_surface/3D_surface_models/
clear all
dirs = {'septin_overexpression', 'wildtype_control', ...
    'septin_knockdown'} ;

% Do not use the normed signal!

colors = define_colors ;
addpath /mnt/data/code/gut_matlab/plotting/shadedErrorBar/
        
%% Default Options
% overwrite previous results if on disk
overwrite = true ;
overwriteImages = false ;
% The number of Laplacian eigenvectors to calculate
nModes = 200;
signalTypes = {'HH', 'radialu', 'dist'};

for jj = 1:length(dirs)
    fns = dir(fullfile(dirs{jj}, '*.ply')) ;

    outdir = fullfile(dirs{jj}, 'analysis') ;

    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end
    
    %%
    for ii = 1:length(fns)
        fn = strrep(fns(ii).name, '.ply', '') ;
        disp(['ii = ' num2str(ii) ': ' fn])

        outfn1 = fullfile(outdir, ...
            [fn '_powerSpectrum_' signalTypes{1} '.mat']) ;
        outfn2 = fullfile(outdir, ...
            [fn '_powerSpectrum_' signalTypes{2} '.mat']) ;
        if ~exist(outfn1, 'file') || ~exist(outfn2, 'file') || overwrite

            mesh = read_ply_mod(fullfile(fns(ii).folder, [fn '.ply'])) ;

            % On which fields do we compute LBS
            clf
            Ufn = fullfile(outdir, ...
                [fn '_conformalMappingToUnitSphere.mat']) ;
            if ~exist(Ufn, 'file') || overwrite
                U = conformalized_mean_curvature_flow(mesh.v,mesh.f, 'LaplacianType', 'cotangent') ;
                Urescaled = conformalized_mean_curvature_flow(mesh.v,mesh.f, 'LaplacianType', 'cotangent', 'RescaleOutput', true) ;

                [Ur2fit,ind2fit]= farthest_points(Urescaled, 200) ;
                offset= mean(Ur2fit) ;
                ptCloud = pointCloud(Ur2fit - offset) ;
                sphereModel = pcfitsphere(ptCloud,1) ;
                sphereCenter = sphereModel.Center + offset ;
                sphereRadius = sphereModel.Radius ;
                sphereParameters = sphereModel.Parameters ;
                radii = vecnorm(mesh.v - sphereCenter, 2, 2) ;
                radii0 = vecnorm(Urescaled - sphereCenter, 2, 2) ;

                % Check if we get a more uniform radius from using offset in addition
                tmp = mean([offset; sphereCenter]) ;
                if std(radii0) > std(vecnorm(Urescaled - tmp, 2, 2))
                    if std(vecnorm(Urescaled - tmp, 2, 2)) < std(vecnorm(Urescaled - tmp, 2, 2))
                        disp('Using FarthestPtSearch sampling mean only, since here it is most accurate')
                        sphereCenter = offset ;
                    else
                        disp('Using 1/2 * (sphereModel.Center + FarthestPtSearch sampling mean)')
                        sphereCenter = tmp ;
                    end
                    sphereRadius = mean(vecnorm(Ur2fit - sphereCenter, 2, 2)) ;
                end
                radii = vecnorm(mesh.v - sphereCenter, 2, 2) ;
                radii0 = vecnorm(Urescaled - sphereCenter, 2, 2) ;

                %% plot the result
                clf
                subplot(2, 2, 1)
                 trisurf(triangulation(mesh.f, mesh.v), 'facecolor', [0., 0., 0], ...
                     'edgecolor', 'none', 'facealpha', 0.2); 
                 hold on;
                 trisurf(triangulation(mesh.f, Urescaled), 'facecolor', [0.5, 0.5, 0], ...
                     'edgecolor', 'none', 'facealpha', 0.2); 
                 hold on; 
                 %  %h = plot(sphereModel) ;
                 %  set(h, 'FaceAlpha', 0.2) ;
                 %  set(h, 'EdgeColor', 'none')
                 %  set(h, 'faceColor',  [1, 1, 0.2]) ;
                 camlight
                 axis equal ;
                 grid off ;
                 trisurf(triangulation(mesh.f, mesh.v), 'facecolor', [0.2, 0.2, 1],...
                     'edgecolor', 'none', 'facealpha', 0.2)
                 axis equal; 
                 grid off ; % axis off ;
                 subplot(2, 2, 2)
                 trisurf(triangulation(mesh.f, Urescaled), ...
                     'edgecolor', 'none', 'facevertexCdata', radii)
                 axis equal ;
                 grid off ;
                 title('radial position of surface');
                 clims = caxis ;

                 subplot(2, 2, 3)
                 trisurf(triangulation(mesh.f, Urescaled), ...
                     'edgecolor', 'none', 'facevertexCdata', radii0)
                 axis equal ;
                 grid off
                 title('radius of mapped surface');
                 caxis(clims)

                 set(gcf, 'color', 'w')
                 figfn = fullfile(outdir,[ fn '_sphericalFit.png']) ;
                 saveas(gcf, figfn)

                %% save the result
                save(Ufn, 'U', 'Urescaled', 'radii', 'radii0', ...
                    'sphereCenter', 'sphereRadius', 'sphereParameters')
            else
                load(Ufn, 'U', 'Urescaled', 'radii', 'radii0', ...
                    'sphereCenter', 'sphereRadius', 'sphereParameters')
            end

            Hfn = fullfile(outdir, [fn '_meanCurvature.mat']) ;
            if ~exist(Hfn, 'file') || overwrite 
                DEC = DiscreteExteriorCalculus(mesh.f, mesh.v) ;
                H3d = sum(mesh.vn .* DEC.laplacian(mesh.v), 2) * 0.5 ;
                disp('saving mean curvature')
                save(Hfn, 'H3d')
            else
                load(Hfn, 'H3d')
            end


            %% Construct Mesh Laplace-Beltrami Operator ===============================

            % laplaceBeltrami = construct_laplace_beltrami( face, vertex );
            laplaceBeltrami = cotmatrix( Urescaled - sphereCenter, mesh.f );
            % DEC = DiscreteExteriorCalculus( face, vertex );

            %% Find Eigenfunctions of Laplace-Beltrami Operator =======================
            tic
            disp('computing spherical harmonics')
            [V,~] = eigs(laplaceBeltrami,nModes,0);

            toc

            %% View Results ===========================================================
            close all ;
            figure('units', 'centimeters', 'position', [0,0,10,5])
            for kk = 1:length(signalTypes)
                if strcmpi(signalTypes{kk}, 'HH')
                    ff = H3d ;
                    titleStr = 'Laplace-Beltrami power spectrum of H' ;
                elseif strcmpi(signalTypes{kk}, 'radialu')
                    ff =  radii - radii0 ;
                    titleStr = 'Laplace-Beltrami power spectrum of \deltar' ;
                elseif strcmpi(signalTypes{kk}, 'dist')
                    ff = vecnorm(mesh.v - Urescaled, 2, 2) ;
                    titleStr = 'Laplace-Beltrami power spectrum of \deltax' ;
                end
                % The 'Fourier Transform' of the signal
                rawPowers = V' * ff;
                rawPowersNormV = rawPowers ./ length(mesh.v) ;

                % Plot Results ------------------------------------------------------------
                % plot( abs(x), '.-', 'LineWidth', 0.5, 'MarkerFaceColor', 'b' );
                clf
                subplot(1, 2, 1)
                bar(abs(rawPowers), 'FaceColor',colors(1, :),'EdgeColor','none')
                xlabel('spherical harmonic index')
                ylabel('spectral power')
                
                subplot(1, 2, 2)
                bar(abs(rawPowersNormV), 'FaceColor',colors(1, :),'EdgeColor','none')
                sgtitle(titleStr);
                xlabel('spherical harmonic index')
                ylabel('spectral power - normed')
                figfn = fullfile(outdir, ...
                    [fn '_powerSpectrum_' signalTypes{kk} '.png']) ;
                saveas(gcf, figfn)


                %% Sort by \ell value
                llvals = [] ;
                lls = 0:13 ;
                dmyk = 1 ;
                powers = zeros(numel(lls),1) ;
                powersNormV = powers ;
                for Lind = 1:numel(lls)
                    ll = lls(Lind) ;
                    for qq = 1:(2*ll + 1)
                        llvals(dmyk)= ll ;
                        powers(Lind) = powers(Lind) + abs(rawPowers(dmyk)) ;
                        powersNormV(Lind) = powersNormV(Lind) + abs(rawPowersNormV(dmyk)) ;
                        % disp(['done accounting for dmyk=' num2str(dmyk) '-> L=' num2str(ll)])
                        dmyk = dmyk + 1 ;
                    end
                end
                
                clf
                subplot(1, 2, 1)
                bar(lls, powers, 'FaceColor',colors(1, :),'EdgeColor','none')
                xlabel('degree of roughness')
                ylabel('spectral power')
                subplot(1, 2, 2)
                bar(lls, powersNormV, 'FaceColor',colors(1, :),'EdgeColor','none')
                xlabel('degree of roughness')
                ylabel('spectral power - normed')
                sgtitle(titleStr)
                figfn = fullfile(outdir, ...
                     [fn '_powerSpectrumSummed_' signalTypes{kk} '.pdf']) ;
                saveas(gcf, figfn)
                % save results
                outfn = fullfile(outdir, ...
                    [fn '_powerSpectrum_' signalTypes{kk} '.mat']) ;
                save(outfn, 'powers', 'powersNormV', ...
                    'lls', 'llvals', 'rawPowers', ...
                    'rawPowersNormV')

            end
        end
    end

    %% Compare all surfaces
    for qq = 1:length(signalTypes)
        signalType = signalTypes{qq} ;
        for ii = 1:length(fns)
            fn = strrep(fns(ii).name, '.ply', '') ;
            outfn = fullfile(outdir, ...
                [fn '_powerSpectrum_' signalType '.mat']) ;
            load(outfn, 'powers', 'powersNormV', 'llvals')

            if ii == 1
                powersAll = powers ;
                powersAllNormV = powersNormV ;
            else
                powersAll = cat(2, powersAll, powers) ;
                powersAllNormV = cat(2, powersAllNormV, powersNormV) ;
                assert(~all(powersNormV == powersAllNormV(:, 1)))
            end
        end
        clf
        bar(lls, powersAll, 'edgecolor', 'none') 
        xlabel('degree of roughness')
        ylabel('spectral power')
        title('Comparison across surfaces')
        saveas(gcf, fullfile(outdir, ...
            ['comparison_of_powerSpectra_' signalType '.pdf']))
        
        clf
        bar(lls, powersAllNormV, 'edgecolor', 'none') 
        xlabel('degree of roughness')
        ylabel('spectral power')
        title('Comparison across surfaces -- normed')
        saveas(gcf, fullfile(outdir, ...
            ['comparison_of_powerSpectraNormV_' signalType '.pdf']))
    
        % save stats for these powers
        save(fullfile(outdir, [signalType '_spectralPowers.mat']), ...
            'powersAll', 'powersAllNormV', 'lls', 'dirs')
    end
end


%% Compare all experiment cases 
close all
colors = [colors(1, :); [0.2,0.2,0.2]; colors(3, :)] ;
clc
meanPowers = {} ;
stdPowers = {} ;
stePowers = {} ;
for normV = [true, false] 
    for qq = 1:length(signalTypes)
        signalType = signalTypes{qq} ;

        hs = {} ;
        for pp = 1:length(dirs)

            outdir = dirs{pp} ;
            load(fullfile(outdir, 'analysis', [signalType '_spectralPowers.mat']), ...
                'powersAll',  'powersAllNormV', 'lls')

            if normV 
                meanPower = mean(powersAllNormV, 2) ;
                stdPower = std(powersAllNormV, [], 2) ;
                stePower = stdPower ./ sqrt(size(powersAllNormV, 2)) ;
            else
                meanPower = mean(powersAll, 2) ;
                stdPower = std(powersAll, [], 2) ;
                stePower = stdPower ./ sqrt(size(powersAll, 2)) ;
            end
            if pp == 1
                lls0 = lls ;
            else
                assert(all(lls == lls0)) ;
            end
            meanPowers{pp} = meanPower ;
            stdPowers{pp} = stdPower ;
            stePowers{pp}= stePower ;

            lineProps = {'-','color', colors(pp, :)} ;
            means = meanPower' ;
            h =shadedErrorBar(lls, means, stdPower', ...
                'lineProps', lineProps, 'patchSaturation', 0.2) ;
            hold on;

            ylim([0, Inf])
        end

        legend(strrep(dirs, '_', ' '),'AutoUpdate','off')

        % Now add stes

        for pp = 1:length(dirs)
            lineProps = {'-','color', colors(pp, :)} ;
            errorbar(lls, meanPowers{pp}, stePowers{pp}, '.')
            % h =shadedErrorBar(lls, meanPowers{pp}, stePowers{pp}, ...
            %    'lineProps', lineProps, 'patchSaturation', 0.2) ;

            % h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            % h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end

        xlabel('degree of roughness')
        ylabel('spectral power')
        if normV
            title('Comparison across conditions -- normed')
            saveas(gcf, ['comparison_of_powerSpectraNormV_' signalType '.pdf'])
        else
            title('Comparison across conditions')
            saveas(gcf, ['comparison_of_powerSpectra_' signalType '.pdf'])
        end
       close all

       %% Save data
       if normV
           save(['statisticsNormV_' signalTypes{qq} '.mat'], ...
               'meanPowers', 'stdPowers', 'stePowers', 'dirs')
       else
           save(['statistics_' signalTypes{qq} '.mat'], ...
               'meanPowers', 'stdPowers', 'stePowers', 'dirs')
       end
    end
end
