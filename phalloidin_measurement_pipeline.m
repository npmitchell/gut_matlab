addpath('/mnt/data/code/gut_matlab/data_handling/')

% Load the image
optoStyle = 'OCRL' ;
bluecolor = [0    0.4470    0.7410] ;
redcolor = [0.8500    0.3250    0.0980] ;
datDir = ['/mnt/data/confocal_data/gut/phalloidin_OCRL/', ...
    '20210110_phalloidin647_48YG4kOCRL_PFA0108_5min470nm45mPFA_0p3um_4p5x63x/', ...
    'for_segmentation'] ;
fns = dir(fullfile(datDir, '*.png')) ;
membraneChannel = 1 ;
optoChannel = 2 ;
actinChannel = 3 ;
opto_thres = 15 ;
overwrite = false ;
pix2um = dlmread(fullfile(datDir, 'pix2micron.txt')) ;
fid = fopen(fullfile(datDir, 'embryo_stages.txt'));
stages = textscan(fid,'%s%s');
fclose(fid);

% Create folder for segmentation snaps
if ~exist(fullfile(fns(1).folder, 'seg_images'), 'dir')
    mkdir(fullfile(fns(1).folder, 'seg_images'))
end

for qq = 1:length(fns)
    name = fns(qq).name ;
    name = name(1:end-4) ;
    
    outfn = fullfile(fns(qq).folder, [fns(qq).name(1:end-4) '_cells.mat']) ;
    
    if ~exist(outfn, 'file') || overwrite
        fn = fullfile(fns(qq).folder, [name '.png']) ;
        im = imread(fn) ;
        im = permute(im, [2,1,3]) ;

        probfn = fullfile(fns(qq).folder, [name '_Probabilities.h5']) ;
        prob = h5read(probfn, '/exported_data') ;
        mem = squeeze(prob(1, :, :)) ;
        bg = squeeze(prob(2, :, :)) ;
        bw = mem > 0.5 ;

        % Skeletonize
        skel = bwskel(bw, 'MinBranchLength',50);
        imshow(labeloverlay(im, skel, 'Transparency',0))
        % Segment into regions
        CC = bwconncomp(imcomplement(skel), 4) ;

        % Choose good regions
        clickfn = fullfile(fns(qq).folder, [fns(qq).name(1:end-4) '_clicks.mat']) ;
        if ~exist(clickfn, 'file')
            clf
            opto_tmp = squeeze(im(:, :, optoChannel)) ;
            subplot(1, 2, 1)
            labels = bwlabel(imcomplement(skel), 4) ;
            imshow(labeloverlay(labels, skel, 'Transparency',0.5))
            subplot(1, 2, 2)
            if mean(opto_tmp(:)) > opto_thres
                imshow(labeloverlay(im * 2, skel, 'Transparency',0.5))
            else
                imshow(labeloverlay(im * 5, skel, 'Transparency',0.5))
            end
            title('Select regions to use, then press enter.')
            [xi, yi] = getpts() ;
            save(clickfn, 'xi', 'yi')
        else
            load(clickfn, 'xi', 'yi')
        end

        % Find which regions these points are inside
        labels = bwlabel(imcomplement(skel), 4) ;
        regs = zeros(length(xi), 1) ;
        for ptId = 1:length(xi)    
            regs(ptId) = labels(round(yi(ptId)), round(xi(ptId))) ;
        end
        % check that there are no duplicates
        try
            assert(numel(unique(regs)) == numel(sort(regs)))
        catch
            error('you selected duplicate regions. Should we remove those?')
        end

        % Measure intensities at membrane and cytoplasm of those regions
        results = {} ;
        for regId = 1:length(regs)
            regi = regs(regId) ;
            % what is the actin in membrane?
            regMask = 0*bw ;
            regMask(labels == regi) = 1 ;
            memMask = regMask .* bw ;
            cytMask = regMask .* imcomplement(bw) ;
            membrane = double(squeeze(im(:, :, membraneChannel))) ;
            mem_mem = membrane .* memMask ;
            mem_cyt = membrane .* cytMask ;
            mem_mem(~memMask) = NaN ;
            mem_cyt(~cytMask) = NaN ;

            % Get distance from boundary -- note we dilate by 1 pixel to
            % include boundary
            SE = strel('disk', 1) ;
            dists = bwdist(logical(~imdilate(regMask, SE))) ;
            dists(dists < 1) = NaN ;
            dist = dists(~isnan(dists)) ;
            tmp = membrane(~isnan(dists)) ;
            [dist, idx] = sort(dist) ;
            memsort = tmp(idx) ;
            mem_dist = [dist, memsort] ;
            [binx, meanmem, stdmem] = binDataMeanStd(dist, memsort) ;
            mem_meanstd = [binx, meanmem, stdmem] ;

            % actin channel
            actin = double(squeeze(im(:, :, actinChannel))) ;
            actin_mem = actin .* memMask ;
            actin_cyt = actin .* cytMask ;
            actin_mem(~memMask) = NaN ;
            actin_cyt(~cytMask) = NaN ;
            tmp = actin(~isnan(dists)) ;
            actinsort = tmp(idx) ;
            actin_dist = [dist, actinsort] ;
            [binx, meanact, stdact] = binDataMeanStd(actin_dist(:, 1), actin_dist(:, 2)) ;
            actin_meanstd = [binx, meanact, stdact] ;

            % Decide if opto
            optog = double(squeeze(im(:, :, optoChannel))) ;
            optog_mem = optog .* memMask ;
            optog_cyt = optog .* cytMask ;
            optog_mem(~memMask) = NaN ;
            optog_cyt(~cytMask) = NaN ;
            opto2cyt = nanmean(optog_mem(:)) / nanmean(optog_cyt(:)) ;
            optog_reg = optog .* regMask ;
            optog_reg(~regMask) = NaN ;
            optog_mean = nanmean(optog_reg(:))  ;


            % Plot to check
            clf
            cyan = [ zeros(256,1), linspace(0,1,256)', linspace(0,1,256)'] ;
            yellow = [ linspace(0,1,256)', linspace(0,1,256)', zeros(256,1)] ;
            red = [linspace(0,1,256)', zeros(256,2)] ;
            colormap(cyan); 
            subplot(3, 3, 1)
            imagesc(membrane) 
            axis equal; axis off
            subplot(3, 3, 2)
            imagesc(mem_mem)
            axis equal; axis off
            subplot(3, 3, 3)
            imagesc(mem_cyt)
            axis equal; axis off
            subplot(3, 3, 4)
            imagesc(optog) 
            axis equal; axis off
            colormap(gca, red); 
            subplot(3, 3, 5)
            imagesc(optog_mem)
            axis equal; axis off
            colormap(gca, red); 
            subplot(3, 3, 6)
            imagesc(optog_cyt)
            axis equal; axis off
            colormap(gca, red); 
            subplot(3, 3, 7)
            imagesc(actin) 
            axis equal; axis off
            colormap(gca, yellow); 
            subplot(3, 3, 8)
            imagesc(actin_mem)
            axis equal; axis off
            colormap(gca, yellow); 
            subplot(3, 3, 9)
            imagesc(actin_cyt)
            axis equal; axis off
            colormap(gca, yellow); 
            if optog_mean > opto_thres
                wtstr = 'Opto ' ;
            else
                wtstr = 'WT ' ;
            end
            sgtitle([wtstr 'cell ' num2str(regId)], ...
                ... % ' of ' strrep(fns(qq).name(1:end-4), '_', '\_')], ...
                'interpreter', 'latex')
            figfn = fullfile(fns(qq).folder, 'seg_images', ...
                [fns(qq).name(1:end-4) sprintf('_cells%05d.png', regId)]) ;
            saveas(gcf, figfn) 

            % Get ratio (holistic)
            mem2cyt = nanmean(actin_mem(:)) / nanmean(actin_cyt(:)) ;

            % Get stage
            embryostring = strsplit(name, 'e') ;
            embryostring = embryostring{2} ;
            embryostring = strsplit(embryostring, '_') ;
            embryostring = ['e' embryostring{1}] ;
            stage = stages{2}(find(strcmpi(stages{1}, embryostring))) ;
            stage = str2double(stage{1}) 
            
            results{regId} = struct('actin_mem', actin_mem, ...
                'actin_cyt', actin_cyt, ...
                'optog_mem', optog_mem, ...
                'optog_cyt', optog_cyt, ...
                'optog_mean', optog_mean, ...
                'mem_mem', mem_mem, ...
                'mem_cyt', mem_cyt, ...
                'mem2cyt', mem2cyt, ...
                'mem_dist', mem_dist, ...
                'actin_dist', actin_dist, ...
                'mem_meanstd', mem_meanstd, ...
                'actin_meanstd', actin_meanstd, ...
                'opto', opto2cyt, ...
                'stage', stage) ;
        end

        % Save results for this image
        save(outfn, 'results')
    else
        disp(['already done with image: ' outfn])
    end
end



%% All results on one plot
clf
acts = {} ;
for qq = 1:length(fns)
    name = fns(qq).name ;
    name = name(1:end-4) ;
    disp(name)
    outfn = fullfile(fns(qq).folder, [fns(qq).name(1:end-4) '_cells.mat']) ;
    load(outfn, 'results')

    for id = 1:length(results)
        disp(['cell ' num2str(id)])
        act = results{id}.actin_meanstd ;
        if results{id}.optog_mean > opto_thres
            p1 = scatter(act(:, 1) * pix2um, act(:, 2), 4, 'markeredgealpha', 0.3, ...
                'markeredgecolor', redcolor) ;
            plot(act(:, 1) * pix2um, act(:, 2), '--', 'color', redcolor)
        else
            p2 = scatter(act(:, 1) * pix2um, act(:, 2), 8, 'markeredgealpha', 0.3, ...
                'markeredgecolor', bluecolor) ;
            plot(act(:, 1) * pix2um, act(:, 2), '--', 'color', bluecolor)
        end
        hold on;
        
        % Check that this is a unique cell
        for pp = 1:length(acts)
            if numel(act) == numel(acts{pp})
                if all(all(act - acts{pp}) < 1)
                    error('here')
                end
            end
        end
        if isempty(pp)
            acts{1} = act ;
        else
            acts{pp+1} = act ;
        end
    end
end
legend([p1, p2], {optoStyle, 'WT'})
xlabel('distance from membrane [$\mu$m]', 'interpreter', 'latex')
ylabel('Intensity [a.u.]', 'interpreter', 'latex')
saveas(gcf, fullfile(datDir, 'seg_images', 'intensity_versus_distance.png'))


%% Normed plot
clf
acts = {} ;
for qq = 1:length(fns)
    name = fns(qq).name ;
    name = name(1:end-4) ;
    disp(name)
    outfn = fullfile(fns(qq).folder, [fns(qq).name(1:end-4) '_cells.mat']) ;
    load(outfn, 'results')

    for id = 1:length(results)
        disp(['cell ' num2str(id)])
        act = results{id}.actin_meanstd ;
        if results{id}.optog_mean > opto_thres
            p1 = scatter(act(:, 1) * pix2um, act(:, 2) / max(act(:, 2)), 4, 'markeredgealpha', 0.6, ...
                'markeredgecolor', redcolor) ;
            % plot(act(:, 1) * pix2um, act(:, 2) / max(act(:, 2)), 'r--') 
        else
            p2 = scatter(act(:, 1) * pix2um, act(:, 2) / max(act(:, 2)), 4, 'markeredgealpha', 0.6, ...
                'markeredgecolor', bluecolor) ;
            % plot(act(:, 1) * pix2um, act(:, 2) / max(act(:, 2)), 'b--')
        end
        hold on;
        
        % Check that this is a unique cell
        for pp = 1:length(acts)
            if numel(act) == numel(acts{pp})
                if all(all(act - acts{pp}) < 1)
                    error('here')
                end
            end
        end
        if isempty(pp)
            acts{1} = act ;
        else
            acts{pp+1} = act ;
        end
    end
end
legend([p1, p2], {optoStyle, 'WT'})
xlabel('distance from membrane [$\mu$m]', 'interpreter', 'latex')
ylabel('Intensity [a.u.]', 'interpreter', 'latex')
saveas(gcf, fullfile(datDir, 'seg_images', 'intensity_versus_distance_normed.png'))

%% Plot versus stage
clf
acts = [] ;
stages = [] ;
pp = 1;
for qq = 1:length(fns)
    name = fns(qq).name ;
    name = name(1:end-4) ;
    disp(name)
    outfn = fullfile(fns(qq).folder, [fns(qq).name(1:end-4) '_cells.mat']) ;
    load(outfn, 'results')
    for id = 1:length(results)
        disp(['cell ' num2str(id)])
        acts(pp) = results{id}.mem2cyt ;
        stages(pp) = results{id}.stage ;
        optog_means(pp) = results{id}.optog_mean ;
        if optog_means(pp) > opto_thres
            p1 = scatter(stages(pp), acts(pp), 10, 'filled', 'markerfacecolor', redcolor, ...
                'markerfacealpha', 0.5) ;
        else
            p2 = scatter(stages(pp), acts(pp), 15, 'markeredgecolor', bluecolor, ...
                'markerfacealpha', 0) ;
        end
        pp = pp + 1 ;
        hold on;
    end
end
legend([p1, p2], {optoStyle, 'WT'}, 'location', 'best')
xlim([min(stages) - 0.5, max(stages) + 0.5])
xlabel('morphological stage', 'interpreter', 'latex')
ylabel('$ {I_\textrm{membrane}} / {I_\textrm{cytoplasm}}$', ...
    'interpreter', 'latex')
saveas(gcf, fullfile(datDir, 'seg_images', 'intensity_versus_stage.png'))