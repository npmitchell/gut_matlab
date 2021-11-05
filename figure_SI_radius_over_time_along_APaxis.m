%% Script for supplementary figure showing radius (s, t)

% Define QS from 48YGAL4-UASCAAXmCh
overwrite = true ;
t0 = QS.t0set() ;
rmax = 400 ;
lobeDir = QS.dir.lobe ;
timePoints = QS.xp.fileMeta.timePoints ;
timeInterval = QS.timeInterval ;
uvexten = QS.uvexten ;  % string like sprintf('_nU%04d_nV%04d', nU, nV) ;
foldImDir = '/mnt/data/analysis/gut_figSI_radius_over_time/' ;
spcutMeshBase = QS.fullFileBase.spcutMesh ;
sphiBase = spcutMeshBase ;
[~, ~, ~, xyzlim] = QS.getXYZLims() ;

foldfn = fullfile(lobeDir, ['fold_locations_sphi' uvexten '_avgpts.mat']) ;
if exist(foldfn, 'file') 
    disp('Loading lobes')
    % Save the fold locations as a mat file
    load(foldfn, 'ssfold', 'folds', 'ssfold_frac', 'ssmax', 'fold_onset', ...
        'rssfold', 'rssfold_frac', 'rssmax', 'rmax')
end

if ~exist(foldImDir, 'dir')
    mkdir(foldImDir)
end

methods = {'ringpath', 'avgpts'} ;
for qq = 1:2
    method = methods{qq} ;
    mexten = ['_' methods{qq} ] ;
    
    % Plot location of each fold superimposed on the mean radius for uniformly
    % sampled DVhoops
    blue = [0 0.4470 0.7410] ;
    red = [0.8500    0.3250    0.0980] ;
    maroon = [0.6350    0.0780    0.1840]; 
    yellow = [0.9290    0.6940    0.1250 ]; 
    maxrmax = max(rmax) ;
    maxss = max(ssmax) ;
    for kk = [114, 144, 174, 204, 234] % 1:length(timePoints)
        % Translate to which timestamp
        t = timePoints(kk) ;
        timestr = sprintf('_%04d', t) ;

        ofn = fullfile(foldImDir, ['radii_folds' dvexten mexten timestr '.pdf']) ;
        if ~exist(ofn, 'file') || overwrite      
            fig = figure('visible', 'off', 'units', 'centimeters',...
                'position', [0,0,7,7]) ;  
            load(sprintf(sphiBase, t), 'spcutMesh') ;

            % Find the average radius over the uniformly sampled hoops. 
            % (make radius 1d)
            rad = mean(spcutMesh.radii_from_mean_uniform_rs, 2) ;
            minrad = min(spcutMesh.radii_from_mean_uniform_rs, [], 2) ;
            maxrad = max(spcutMesh.radii_from_mean_uniform_rs, [], 2) ;

            stdradlo = rad - std(spcutMesh.radii_from_mean_uniform_rs, [], 2) ;
            stdradhi = rad + std(spcutMesh.radii_from_mean_uniform_rs, [], 2) ;

            if strcmp(method, 'avgpts')
                ss = spcutMesh.avgpts_ss ;
            elseif strcmp(method, 'ringpath')
                ss = spcutMesh.ringpath_ss ;
            end

            % Save plot of radius with fold locations marked
            fill([ss; flipud(ss)], [minrad; flipud(maxrad)], blue, ...
                'facealpha', 0.3, 'edgecolor', 'none')
            hold on;
            fill([ss; flipud(ss)], [stdradlo; flipud(stdradhi)], maroon, ...
                'facealpha', 0.3, 'edgecolor', 'none')
            plot(ss, rad, 'Color', 'k')
            plot(ss(folds(kk, :)), rad(folds(kk, :)), 'o', 'Color', yellow)
            xlim([0, maxss])
            ylim([0, maxrmax])
            % axis equal
            if strcmpi(method, 'ringpath')
                xlabel('geodesic position along AP axis [$\mu$m]', 'interpreter', 'latex')
            else
                xlabel('AP position of centerline position [$\mu$m]', 'interpreter', 'latex')
            end
            title(['$t=$' num2str(t-t0) ' min'], 'interpreter', 'latex')
            ylabel('radius [$\mu$m]', 'interpreter', 'latex')
            disp(['Saving radius vs ss plot for t=' num2str(t)])
            saveas(fig, ofn)
            close all
        end
        

        %% Show what we are measuring
        close all   
        fig = figure('visible', 'on', 'units', 'centimeters',...
            'position', [0,0,18,7]) ;  
        % tidx0 = find(timePoints-t0  == 60) ;
        % Translate to which timestamp
        % t = timePoints(tidx0) ;
        % timestr = sprintf('_%04d', t) ;

        load(sprintf(QS.fullFileBase.spcutMeshSmRSC, t), 'spcutMeshSmRSC') ;   
        load(sprintf(sphiBase, t), 'spcutMesh') ;
        mesh = spcutMeshSmRSC ;

        % geodesic and centerline length
        avgpts = spcutMesh.avgpts ;
        ssa = spcutMesh.avgpts_ss ;
        ssr = spcutMesh.ringpath_ss ;

        % check it
        test = cumsum(vecnorm(diff(avgpts), 2, 2)) ;
        assert(all(test == ssa(2:end)))

        cntrline_ss = ssa .* ones(100, 99) ;
        ringpath_ss = ssr .* ones(100, 99) ;
        xx = reshape(mesh.v(:, 1), [100, 99]) ;
        yy = reshape(mesh.v(:, 2), [100, 99]) ;
        zz = reshape(mesh.v(:, 3), [100, 99]) ;
        cmax = 400 ;

        % now plot each
        colorsV = viridis(256) ;
        subplot(1, 2, 1)
        for ii = 1:size(xx, 1)
            cidx = min(round(ssr(ii) / cmax * 255 + 1), 256) ;
            colorii = colorsV(cidx, :) ;
            plot3([xx(ii, :), xx(ii, 1)], ...
                [yy(ii, :), yy(ii, 1)], ...
                [zz(ii, :), zz(ii, 1)], ...
                'color', colorii)
            hold on;
        end
        caxis( [0, cmax])
        axis equal
        axis off
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        grid off
        %xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
        %zlabel('dv position [$\mu$m]', 'interpreter', 'latex')
        %ylabel('lateral position [$\mu$m]', 'interpreter', 'latex')
        cb = colorbar('location', 'southoutside') ;
        pos = get(cb, 'position') ;
        set(cb, 'position', [pos(1) + 0.2*pos(3), pos(2), pos(3)*0.5, pos(4)])
        pos = get(gca, 'position') ;
        set(gca, 'position', [pos(1), pos(2) + 0.0*pos(4), pos(3), pos(4)])
        caxis( [0, cmax])
        ylabel(cb, 'geodesic position along AP axis [$\mu$m]', ...
            'interpreter', 'latex')
        subplot(1, 2, 2)
        for ii = 1:size(xx, 1)
            cidx = min(round(ssa(ii) / cmax * 255 + 1), 256) ;
            colorii = colorsV(cidx, :) ;
            plot3([xx(ii, :), xx(ii, 1)], ...
                [yy(ii, :), yy(ii, 1)], ...
                [zz(ii, :), zz(ii, 1)], ...
                'color', colorii)
            hold on;
        end
        axis equal
        axis off
        caxis( [0, cmax])
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        grid off
        %xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
        %zlabel('dv position [$\mu$m]', 'interpreter', 'latex')
        %ylabel('lateral position [$\mu$m]', 'interpreter', 'latex')
        cb = colorbar('location', 'southoutside') ;
        pos = get(cb, 'position') ;
        set(cb, 'position', [pos(1) + 0.2*pos(3), pos(2), pos(3)*0.5, pos(4)])
        pos = get(gca, 'position') ;
        set(gca, 'position', [pos(1), pos(2) + 0.0*pos(4), pos(3), pos(4)])
        ylabel(cb, 'AP position of centerline position [$\mu$m]', ...
            'interpreter', 'latex')

        fig.Renderer='Painters';
        saveas(gcf, fullfile(foldImDir, sprintf('demo_ssr_ssa_%06d.pdf', t)))
    end
end






