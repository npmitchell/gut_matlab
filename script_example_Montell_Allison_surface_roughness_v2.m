cd /mnt/data/3D_surface_models/3D_surface_models
clear all
dirs = {'septin_knockdown', 'wildtype_control', ...
   'septin_overexpression' } ;
shorthand = {'KD', 'WT', 'OE'} ;
% Do not use the normed signal!

%% Testing code here
sumH = cell(length(dirs), 1) ;
sumHf = cell(length(dirs), 1) ;

colors = [   0.9000    0.5500    0.5500
    0.5000    0.5000    0.5000
    0.6200    0.7600    0.8400] ;
clf
for jj = 1:length(dirs)
    outdir = fullfile(dirs{jj}, 'analysis') ;
    fns = dir(fullfile(dirs{jj}, '*.ply')) ;
    for ii = 1:length(fns)
        
        % mesh = read_ply_mod(fullfile(fns(ii).folder, fns(ii).name)) ;
        fn = strrep(fns(ii).name, '.ply', '') ;
        Hfn = fullfile(outdir, [fn '_meanCurvature.mat']) ;
        % load(Hfn, 'H3d', 'areas')
        % [V2F, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
        % H3d_faces = (V2F * H3d) ;
        load(Hfn, 'H3d', 'areas', 'H3d_faces')
        sumH{jj}(ii) = sum(abs(H3d)) ;
        sumHf{jj}(ii) = sum(abs(H3d_faces) .* areas / sum(areas)) ;
    end
    
    % errorbar(jj, mean(sumH{jj}), std(sumH{jj}) / sqrt(length(sumH{jj})), 'o')
    plot(jj * ones(length(sumHf{jj})), sumHf{jj}, '.', 'color', colors(jj, :))
    hold on;
    errorbar(jj, mean(sumHf{jj}), std(sumHf{jj}) / sqrt(length(sumHf{jj})), 's', 'color', colors(jj, :))
    xticks([1,2,3])
    xticklabels(shorthand)
    ylabel('$\langle |H| \rangle$ [$\mu$m$^{-1}$]', 'interpreter', 'latex')
    xlim([0, 4])
    hold on
end
saveas(gcf, 'average_absMeanCurvature.png')

%% Testing code part 2, coloring a surface by radialu

for jj = 2
    fns = dir(fullfile(dirs{jj}, '*.ply')) ;
        
    outdir = fullfile(dirs{jj}, 'analysis') ;
    for ii = 1:length(fns)        
        fn = strrep(fns(ii).name, '.ply', '') ;
        disp(['ii = ' num2str(ii) ': ' fn])
        Ufn = fullfile(outdir, ...
            [fn '_conformalMappingToUnitSphere.mat']) ;
        load(Ufn, 'U', 'Urescaled', 'radii', 'radii0', ...
            'sphereCenter', 'sphereRadius', 'sphereParameters')

        mesh = read_ply_mod(fullfile(fns(ii).folder, fns(ii).name)) ;
        fn = strrep(fns(ii).name, '.ply', '') ;

        clf
        patch( 'Faces', mesh.f, 'Vertices', mesh.v, ...
            'FaceVertexCdata', radii-radii0, ...
            'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceLighting', 'gouraud' );
        % caxis([-1, 1])
        % colorbar
        colormap(bam)
        caxis([-max(abs(radii-radii0)), max(abs(radii-radii0))])
        axis equal
        % axis off
        set(gcf, 'color', 'w')

        camlight('left')
        xlabel('x [pix]')
        zlabel('z [pix]')
        saveas(gcf, fullfile(dirs{jj}, ['radialu_' fn '.png']))
        clf
    end
end

%%

colors = define_colors ;
addpath /mnt/data/code/gut_matlab/plotting/shadedErrorBar/
        
%% Default Options
% overwrite previous results if on disk
overwrite = false ;
overwrite2 = false ;
overwriteImages = false ;
% The number of Laplacian eigenvectors to calculate
nModes = 1000;
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
        if ~exist(outfn1, 'file') || ~exist(outfn2, 'file') || overwrite || overwrite2

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
                areas = 0.5 * doublearea(mesh.v, mesh.f) ;
                [V2F, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
                H3d_faces = (V2F * H3d) ;
                disp('saving mean curvature')
                save(Hfn, 'H3d', 'areas', 'H3d_faces')
            else
                load(Hfn, 'H3d')
            end

            % At this point we have performed mean curvature flow to the surface
            % according to "Can mean curvature flow be made non-singular?" 
            % [Kazhdan et al. 2012], giving a sphere of radius R centered at 
            % r0. We found the distance of the            % surface vertices from r0 and subtracted R, the mapped radial
            % distance from r0. This gives us a measure of how much the
            % surface is extended beyond R or retracted from R at each
            % point on the mesh. This defines the 'corrugation' field 
            % \deltar, which can be expressed either
            % as a function of position on the cell surface or as a function of
            % position on the mapped sphere found by mean curvature flow, with 
            % values \deltar(x) mapped to \deltar(.
            % Conveniently, scalar field patterns on the sphere can be compared
            % directly across samples. Therefore, we then decompose the
            % scalar field defined on the mapped (spherical) surface into 
            % spherical harmonics. Spherical harmonics are a set of 
            % functions which can reproduce the original field when multiplied by 
            % weights and summed together. More formally, they are 
            % eigenmodes of the Laplace-Beltrami operator defined on the sphere.
            % The inner product between a given eigenmode and the measured
            % corrugation field \delta r.

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
                
                outfn = fullfile(outdir, ...
                    [fn '_powerSpectrum_' signalTypes{kk} '.mat']) ;
                if ~exist(outfn, 'file') || overwrite || overwrite2
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
                    % rawPowersNormV = rawPowers ./ length(mesh.v) ;

                    % Plot Results ------------------------------------------------------------
                    % plot( abs(x), '.-', 'LineWidth', 0.5, 'MarkerFaceColor', 'b' );
                    clf
                    % subplot(1, 2, 1)
                    bar(abs(rawPowers), 'FaceColor',colors(1, :),'EdgeColor','none')
                    xlabel('spherical harmonic index')
                    ylabel('spectral power')
                    % subplot(1, 2, 2)
                    % bar(abs(rawPowersNormV), 'FaceColor',colors(1, :),'EdgeColor','none')
                    sgtitle(titleStr);
                    % xlabel('spherical harmonic index')
                    % ylabel('spectral power - normed')
                    figfn = fullfile(outdir, ...
                        [fn '_powerSpectrum_' signalTypes{kk} '.png']) ;
                    saveas(gcf, figfn)


                    %% Sort by \ell value
                    llvals = [] ;
                    lls = 0:30 ;
                    dmyk = 1 ;
                    powers = zeros(numel(lls),1) ;
                    % powersNormV = powers ;
                    for Lind = 1:numel(lls)
                        ll = lls(Lind) ;
                        for qq = 1:(2*ll + 1)
                            llvals(dmyk)= ll ;
                            powers(Lind) = powers(Lind) + abs(rawPowers(dmyk)) ;
                            % powersNormV(Lind) = powersNormV(Lind) + abs(rawPowersNormV(dmyk)) ;
                            % disp(['done accounting for dmyk=' num2str(dmyk) '-> L=' num2str(ll)])
                            dmyk = dmyk + 1 ;
                        end
                    end

                    clf
                    % subplot(1, 2, 1)
                    bar(lls, powers, 'FaceColor',colors(1, :),'EdgeColor','none')
                    xlabel('degree of roughness')
                    ylabel('spectral power')
                    % subplot(1, 2, 2)
                    % bar(lls, powersNormV, 'FaceColor',colors(1, :),'EdgeColor','none')
                    % xlabel('degree of roughness')
                    % ylabel('spectral power - normed')
                    sgtitle(titleStr)
                    figfn = fullfile(outdir, ...
                         [fn '_powerSpectrumSummed_' signalTypes{kk} '.pdf']) ;
                    saveas(gcf, figfn)
                    % save results
                    save(outfn, 'powers', ...
                        'lls', 'llvals', 'rawPowers')
                else
                    load(outfn, 'powers', 'lls', 'llvals', 'rawPowers')

                end
            end
        
        else
            load(outfn2, 'powers', 'lls', 'llvals', 'rawPowers')
            
        end
    end

    %% Compare all surfaces
    for qq = 1:length(signalTypes)
        signalType = signalTypes{qq} ;
        for ii = 1:length(fns)
            fn = strrep(fns(ii).name, '.ply', '') ;
            outfn = fullfile(outdir, ...
                [fn '_powerSpectrum_' signalType '.mat']) ;
            load(outfn, 'powers', 'llvals')

            if ii == 1
                powersAll = powers ;
                % powersAllNormV = powersNormV ;
            else
                powersAll = cat(2, powersAll, powers) ;
                % powersAllNormV = cat(2, powersAllNormV, powersNormV) ;
                % assert(~all(powersNormV == powersAllNormV(:, 1)))
            end
        end
        clf
        bar(lls, powersAll, 'edgecolor', 'none') 
        xlabel('degree of roughness')
        ylabel('spectral power')
        title('Comparison across surfaces')
        saveas(gcf, fullfile(outdir, ...
            ['comparison_of_powerSpectra_' signalType '.pdf']))
        
        % clf
        % bar(lls, powersAllNormV, 'edgecolor', 'none') 
        % xlabel('degree of roughness')
        % ylabel('spectral power')
        % title('Comparison across surfaces -- normed')
        % saveas(gcf, fullfile(outdir, ...
        %    ['comparison_of_powerSpectraNormV_' signalType '.pdf']))
    
        % save stats for these powers
        save(fullfile(outdir, [signalType '_spectralPowers.mat']), ...
            'powersAll',  'lls', 'dirs')
    end
end


%% Compare all experiment cases 
close all
colors = [   0.9000    0.5500    0.5500
    0.5000    0.5000    0.5000
    0.6200    0.7600    0.8400] ;
clc
meanPowers = {} ;
stdPowers = {} ;
stePowers = {} ;
% for normV = [false] 
    for qq = 1:length(signalTypes)
        signalType = signalTypes{qq} ;

        hs = {} ;
        for pp = 1:length(dirs)

            outdir = dirs{pp} ;
            load(fullfile(outdir, 'analysis', [signalType '_spectralPowers.mat']), ...
                'powersAll',  'lls')


            meanPower = mean(powersAll, 2) ;
            stdPower = std(powersAll, [], 2) ;
            stePower = stdPower ./ sqrt(size(powersAll, 2)) ;

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
            errorbar(lls, meanPowers{pp}, stePowers{pp}, '.', 'color', colors(pp, :))
            % h =shadedErrorBar(lls, meanPowers{pp}, stePowers{pp}, ...
            %    'lineProps', lineProps, 'patchSaturation', 0.2) ;

            % h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            % h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end

        xlabel('degree of roughness')
        ylabel('spectral power')
        % if normV
        %     title('Comparison across conditions -- normed')
        %     saveas(gcf, ['comparison_of_powerSpectraNormV_' signalType '.pdf'])
        % else
            title('Comparison across conditions')
            saveas(gcf, ['comparison_of_powerSpectra_' signalType '.pdf'])
        % end
        
        xlim([10, 20])
        ylim([200,1200])
        axis square
            saveas(gcf, ['comparison_of_powerSpectra_' signalType '_zoom2.pdf'])
       close all

       %% Save data
       
           save(['statistics_' signalTypes{qq} '.mat'], ...
               'meanPowers', 'stdPowers', 'stePowers', 'dirs', 'lls')
       
    end
% end


%% Bin powers
qq = 2 ;
load(['statistics_' signalTypes{qq} '.mat'], ...
               'meanPowers', 'stdPowers', 'stePowers', 'dirs', 'lls')
ind = 2 ;

%% low mode
kd = meanPowers{1}(ind) ;
wt = meanPowers{2}(ind) ;
oe = meanPowers{3}(ind) ;
kd_unc = stePowers{1}(ind) ;
wt_unc = stePowers{2}(ind) ;
oe_unc = stePowers{3}(ind) ;

    num = -kd+wt ;
    denom = sqrt(kd_unc.^2 + wt_unc.^2) ;
    zscore_kdwt = num / denom ;
    pval_kdwt = normcdf(zscore_kdwt);

    num = -kd+oe ;
    denom = sqrt(kd_unc.^2 + oe_unc.^2) ;
    zscore_kdoe = num / denom ;
    pval_kdoe = normcdf(zscore_kdoe);

    num = -wt+oe ;
    denom = sqrt(oe_unc.^2 + wt_unc.^2) ;
    zscore_wtoe = num / denom ;
    pval_wtoe = normcdf(zscore_wtoe);

    % superbar plot
    hf = figure('Position', [100 100 400 400], 'units', 'centimeters');
    clf;
    E = cat(3, kd_unc, wt_unc, oe_unc) ;

    Colors = [
        0.90    0.55    0.55
        0.5     0.5     0.5
        0.62    0.76    0.84
        0.89    0.10    0.11
        0.12    0.47    0.70
        ];
    % Colors = reshape(Colors, [2 2 3]);

    P = [0,  pval_kdwt, pval_kdoe; ...
        pval_kdwt, 0, pval_wtoe;...
        pval_kdoe, pval_wtoe, 0];
    
    E = [kd_unc, wt_unc, oe_unc] ;
    superbar([1,2,3],[kd,wt,oe], 'E', E, 'P', P,...
        'BarFaceColor', Colors) ;
    %, 'Orientation', 'v', ...
    %        'ErrorbarStyle', 'I', 'PLineOffset', 0.1, 'PStarShowGT', false)
    xticks([1,2,3])
    categories = {'KD', 'WT', 'OE'} ;
    xticklabels(categories)
    ylabel(sprintf('Spectral power in mode l=%d', lls(ind)))
    saveas(gcf, sprintf('low_order_mode_comparison_l%02d.pdf', lls(ind)))
    means = [kd,wt,oe] ;
    save(sprintf('low_oder_mode_comparison_l%02.mat', lls(ind)), ...
        'P', 'E', 'means', 'categories')
    
%% high order modes
upperThres = 2 ;
close all
for thres = 0
    % thres  = 4 ;

    kd = sum(meanPowers{1}(lls > thres & lls < upperThres, :)) ;
    wt = sum(meanPowers{2}(lls > thres & lls < upperThres, :)) ;
    oe = sum(meanPowers{3}(lls > thres & lls < upperThres, :)) ;
    kd_unc = sqrt(sum((stePowers{1}(lls > thres & lls < upperThres, :)).^2)) ;
    wt_unc = sqrt(sum((stePowers{2}(lls > thres & lls < upperThres, :)).^2)) ;
    oe_unc = sqrt(sum((stePowers{3}(lls > thres & lls < upperThres, :)).^2)) ;

    % pvalue
    num = kd-wt ;
    denom = sqrt(kd_unc.^2 + wt_unc.^2) ;
    zscore_kdwt = num / denom ;
    pval_kdwt = normcdf(zscore_kdwt);

    num = kd-oe ;
    denom = sqrt(kd_unc.^2 + oe_unc.^2) ;
    zscore_kdoe = num / denom ;
    pval_kdoe = normcdf(zscore_kdoe);

    num = wt-oe ;
    denom = sqrt(oe_unc.^2 + wt_unc.^2) ;
    zscore_wtoe = num / denom ;
    pval_wtoe = normcdf(zscore_wtoe);

    %% superbar plot
    hf = figure('Position', [100 100 400 400], 'units', 'centimeters');
    clf;
    E = cat(3, kd_unc, wt_unc, oe_unc) ;

    Colors = [
        0.90    0.55    0.55
        0.5     0.5     0.5
        0.62    0.76    0.84
        0.89    0.10    0.11
        0.12    0.47    0.70
        ];
    % Colors = reshape(Colors, [2 2 3]);

    P = [0,  pval_kdwt, pval_kdoe; ...
        pval_kdwt, 0, pval_wtoe;...
        pval_kdoe, pval_wtoe, 0];
    
    E = [kd_unc, wt_unc, oe_unc] ;
    superbar([1,2,3],[kd,wt,oe], 'E', E, 'P', P,...
        'BarFaceColor', Colors) ;
    %, 'Orientation', 'v', ...
    %        'ErrorbarStyle', 'I', 'PLineOffset', 0.1, 'PStarShowGT', false)
    xticks([1,2,3])
    categories = {'KD', 'WT', 'OE'} ;
    xticklabels(categories)
    ylabel(sprintf('Spectral power in modes l>%d', thres))
    outfn = sprintf('high_order_mode_comparison_lgt%02d_lt%02d.pdf', thres, upperThres) ;
    saveas(gcf, outfn)
    
    outfn = sprintf('high_order_mode_comparison_lgt%02d_lt%02d.mat', thres, upperThres) ;
    save(outfn, 'P', 'E', 'means', 'categories')
    
end