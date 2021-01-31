addpath('/mnt/data/code/gut_matlab/data_handling/')
addpath('/mnt/data/code/gut_matlab/plotting/shadedErrorBar/')
addpath('/mnt/data/code/gut_matlab/plotting/')
addpath('/mnt/data/code/gut_matlab/plotting/violin/')
addpath('/mnt/data/code/gut_matlab/tiff_handling/')

% Load the image
optoStyle = 'OCRL' ;
bluecolor = [0    0.4470    0.7410] ;
redcolor = [0.8500    0.3250    0.0980] ;
optoColor = redcolor ;
WTcolor = bluecolor ;
memthres = 0.25 ;  % um, max dist from membrane to label as membrane-bound

% Dataset 1 -- jan 8th, 2021
% datDir = ['/mnt/data/confocal_data/gut/phalloidin_OCRL/', ...
%     '20210110_phalloidin647_48YG4kOCRL_PFA0108_5min470nm45mPFA_0p3um_4p5x63x/', ...
%     'for_segmentation'] ;
% membraneChannel = 1 ;
% optoChannel = 2 ;
% actinChannel = 3 ;
% opto_thres = 15 ;
minBranchLength = 250 ;

% Dataset 2 -- jan 10th, 2021
datDir = fullfile('/mnt/data/confocal_data/gut/phalloidin_OCRL/', ...
    '20210114b_phalloidin647_48YG4kOCRL_0109_10m470_10m470PFA_15mPFA_0p3um_1pc1pc0p07pc4x63x_aqP/', ...
    '4x', 'for_segmentation') ;
overDir = fullfile('/mnt/data/confocal_data/gut/phalloidin_OCRL/', ...
    '20210114b_phalloidin647_48YG4kOCRL_0109_10m470_10m470PFA_15mPFA_0p3um_1pc1pc0p07pc4x63x_aqP/', ...
    '0p75x') ;
membraneChannel = 1 ;
optoChannel = 3 ;
actinChannel = 2 ;
opto_thres = -1 ;
minBranchLength = 50 ;
segThres  = 0.6 ;

% Dataset 3
%


%% I/O
fns = dir(fullfile(datDir, '*.png')) ;
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
        bw = mem > segThres ;

        % Skeletonize
        skel = bwskel(bw, 'MinBranchLength',minBranchLength);
        
        imshow(labeloverlay(im, skel, 'Transparency',0))
        pause(1) 
        
        B = bwmorph(skel, 'branchpoints');
        E = bwmorph(skel, 'endpoints');
        [y,x] = find(E);
        B_loc = find(B);
        Dmask = false(size(skel));
        for k = 1:numel(x)
            D = bwdistgeodesic(skel,x(k),y(k));
            distanceToBranchPt = min(D(B_loc));
            Dmask(D < distanceToBranchPt) =true;
        end
        skel = skel - Dmask;
        
        imshow(labeloverlay(im, skel, 'Transparency',0))
        % Segment into regions
        CC = bwconncomp(imcomplement(skel), 4) ;

        % Choose good regions
        clickfn = fullfile(fns(qq).folder, [fns(qq).name(1:end-4) '_clicks.mat']) ;
        if ~exist(clickfn, 'file')
            clf
            opto_tmp = squeeze(im(:, :, optoChannel)) ;
            % subplot(1, 2, 1)
            % labels = bwlabel(imcomplement(skel), 4) ;
            % imshow(labeloverlay(labels, skel, 'Transparency',0.5))
            % subplot(1, 2, 2)
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
            hold on;
            scatter(xi, yi, 10, 'filled')
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
            figfn = fullfile(fns(qq).folder, 'seg_images', ...
                [fns(qq).name(1:end-4) sprintf('_cells%05d.png', regId)]) ;
            if ~exist(figfn, 'file') || true
                clf
                cyan = [ zeros(256,1), linspace(0,1,256)', linspace(0,1,256)'] ;
                yellow = [ linspace(0,1,256)', linspace(0,1,256)', zeros(256,1)] ;
                red = [linspace(0,1,256)', zeros(256,2)] ;
                colormap(cyan); 
                subplot(3, 3, 1)
                imagesc(membrane) 
                caxis([min(membrane(:)), max(membrane(:))])
                axis equal; axis off
                subplot(3, 3, 2)
                imagesc(mem_mem)
                caxis([min(membrane(:)), max(membrane(:))])
                axis equal; axis off
                subplot(3, 3, 3)
                imagesc(mem_cyt)
                caxis([min(membrane(:)), max(membrane(:))])
                axis equal; axis off
                subplot(3, 3, 4)
                imagesc(optog) 
                caxis([min(optog(:)), max(optog(:))])
                axis equal; axis off
                colormap(gca, red); 
                subplot(3, 3, 5)
                imagesc(optog_mem)
                caxis([min(optog(:)), max(optog(:))])
                axis equal; axis off
                colormap(gca, red); 
                subplot(3, 3, 6)
                imagesc(optog_cyt)
                caxis([min(optog(:)), max(optog(:))])
                axis equal; axis off
                colormap(gca, red); 
                subplot(3, 3, 7)
                imagesc(actin) 
                caxis([min(actin(:)), max(actin(:))])
                axis equal; axis off
                colormap(gca, yellow); 
                subplot(3, 3, 8)
                imagesc(actin_mem)
                caxis([min(actin(:)), max(actin(:))])
                axis equal; axis off
                colormap(gca, yellow); 
                subplot(3, 3, 9)
                imagesc(actin_cyt)
                caxis([min(actin(:)), max(actin(:))])
                axis equal; axis off
                colormap(gca, yellow); 

                % If there is a threshold, use it. Otherwise use name
                if opto_thres > 0
                    if optog_mean > opto_thres
                        wtstr = 'Opto ' ;
                        is_opto = true ;
                    else
                        wtstr = 'WT ' ;
                        is_opto = false ;
                    end
                else
                    if contains(name, 'WT')
                        wtstr = 'WT' ;
                        is_opto = false ;
                    elseif contains(name, 'OC')
                        wtstr = 'Opto' ;
                        is_opto = true ;
                    else
                        error('here')
                    end
                end
                sgtitle([wtstr 'cell ' num2str(regId)], ...
                    ... % ' of ' strrep(fns(qq).name(1:end-4), '_', '\_')], ...
                    'interpreter', 'latex')
                saveas(gcf, figfn) 
            end

            % Get ratio (holistic)
            mem2cyt = nanmean(actin_mem(:)) / nanmean(actin_cyt(:)) ;

            % Get stage
            embryostring = strsplit(name, 'e') ;
            embryostring = embryostring{2} ;
            embryostring = strsplit(embryostring, '_') ;
            embryostring = ['e' embryostring{1}] ;
            stage = stages{2}(find(strcmpi(stages{1}, embryostring))) ;
            stage = str2double(stage{1}) ;
            
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
                'stage', stage, ...
                'is_opto', is_opto, ...
                'embryoID', embryostring) ;
            is_opto
        end

        % Save results for this image
        save(outfn, 'results')
    else
        
        load(outfn, 'results')
        embryoIDs_done = true;
        for pp = 1:length(results)
            if ~isfield(results{pp}, 'embryoID')
                embryoIDs_done = false ;
            end
        end
        if embryoIDs_done
            disp(['already done with image: ' outfn])
        else        
            disp(['adding embryoID to: ' outfn])
            % Get stage
            embryostring = strsplit(name, 'e') ;
            embryostring = embryostring{2} ;
            embryostring = strsplit(embryostring, '_') ;
            embryostring = ['e' embryostring{1}] ;
            for pp = 1:length(results)
                results{pp}.embryoID = embryostring ;
            end
            
            % Save results for this image
            save(outfn, 'results')
        end
    end
end

%% Normalization factors
clf
bolis = [] ;
embryoIDs = stages{1} ;
for eID = 2:length(embryoIDs)
    embryoID = embryoIDs{eID} ;
    fns2find = dir(fullfile(overDir, ['MAX_' embryoID '*'])) ;
    outfn = fullfile(fns2find(1).folder, [embryoID '_meanActin.mat']);
    
    if ~exist(outfn, 'file')
        if length(fns2find) > 1
            error('Too many matching files!')
        else
            fn = fullfile(fns2find(1).folder, fns2find(1).name) ;
            im = readTiff4D(fn, 4) ;
            actin = im{actinChannel} ;
            mem = im{membraneChannel} ;
            bf = im{4} ;
            imshow(cat(3, bf, actin + mem, actin))
            mask = roipoly() ;
            % what is the full (bolometric) intensity of actin?
            meanActin = sum(uint8(mask(:)) .* actin(:)) / sum(uint8(mask(:))) ;
            save(outfn, 'mask', 'meanActin')
        end
    else
        load(outfn, 'meanActin')
    end
    bolis(eID) = meanActin ;
end

%% Cluster by embryo vs distance
for normalized = [0, 1]
    clf
    embryos = {} ;
    acts = [] ;
    embryoIDs = stages{1} ;
    for eID = 2:length(embryoIDs)
        % Embryo identification
        embryostr = embryoIDs{eID} ;

        first = true ;
        for qq = 1:length(fns)
            name = fns(qq).name ;
            name = name(1:end-4) ;
            disp(name)
            outfn = fullfile(fns(qq).folder, [fns(qq).name(1:end-4) '_cells.mat']) ;
            if contains(name, embryostr)
                load(outfn, 'results')

                % for each cell of this image
                for id = 1:length(results)
                    if strcmpi(results{id}.embryoID, embryostr)
                        disp(['cell ' num2str(id)])
                        act = results{id}.actin_meanstd ;
                        if normalized
                            act(:, 2:3) = act(:, 2:3) / bolis(eID) ;
                        end
                        if first
                            disp('instantiating acts')
                            is_opto = results{id}.is_opto ;
                            acts = act ;
                            first = false ;
                        else
                            disp('adding to acts')
                            assert(is_opto == results{id}.is_opto) ;
                            acts = [acts; act] ;
                        end
                    end
                end
            end
        end
        is_optos(eID) = is_opto ;

        [midx, amean, astd] = binDataMeanStd(acts(:, 1), acts(:, 2)) ;

        if is_opto
            disp('Opto')            
            if str2double(stages{2}{eID}) > 13
                lineProps = {'-', 'color', redcolor} ;
            else
                lineProps = {'--', 'color', redcolor} ;
            end
            p1 = shadedErrorBar(midx * pix2um,amean, astd, ...
                'patchSaturation', 0.1, 'lineProps', lineProps);
        else
            disp('WT')
            if str2double(stages{2}{eID}) > 13
                lineProps = {'-', 'color',bluecolor} ;
            else
                lineProps = {'--', 'color',bluecolor} ;
            end
            p2 = shadedErrorBar(midx * pix2um,amean, astd, ...
                'patchSaturation', 0.1, 'lineProps',  lineProps);
        end
        hold on;
        pause(0.001) ;
        drawnow
    end
    legend([p1.patch, p2.patch], {optoStyle, 'WT'})
    xlabel('distance from membrane [$\mu$m]', 'interpreter', 'latex')
    if normalized
        ylabel('Intensity relative to mean actin channel', 'interpreter', 'latex')
    else
        ylabel('Intensity [a.u.]', 'interpreter', 'latex')
    end
    if normalized
        saveas(gcf, fullfile(datDir, 'seg_images', 'intensity_versus_distance_embryos_normed.png'))
    else
        saveas(gcf, fullfile(datDir, 'seg_images', 'intensity_versus_distance_embryos.png'))
    end
end


%% Cluster by embryo vs stages -- mem/average or  mem/cyto
for normed = [1, 0 ]
    if normed 
        bw = 0.01 ;
    else
        bw = 0.05 ;
    end
    close all ; 
    clf;
    actsWT = {} ;
    stagesWT = [] ;
    collectWT = 1 ;
    actsOpto = {} ;
    stagesOpto = [] ;
    collectOpto = 1 ;
    % Consider each embryo
    for eID = 1:length(embryoIDs)
        embryoID = embryoIDs{eID} ;
        
        % Clear acts to build up collection
        acts = [] ; 
        if ~contains(embryoID, 'embryoID')
            first = true ;
            pp = 1 ;
            % Find all images of that embryo
            for qq = 1:length(fns)
                name = fns(qq).name ;
                name = name(1:end-4) ;
                disp(name)
                outfn = fullfile(fns(qq).folder, [fns(qq).name(1:end-4) '_cells.mat']) ;

                % Is this image a match to the embryoID under consideration?
                if contains(name, embryoID)
                    disp(['Loading embryoID file: ' name])
                    % grab the mem2cyt ratio for this embryo in this image
                    load(outfn, 'results')
                    for id = 1:length(results)
                        disp(['cell ' num2str(id)])
                        if normed
                            adi = results{id}.actin_dist(:, 1:2) ;
                            [midx, amean, astd] = binDataMeanStd(...
                                adi(:, 1)*pix2um, adi(:, 2), [0, memthres, max(adi(:,1))*pix2um]) ;
                            actin_mem_normed = amean(1) / bolis(eID) ;
                            acts(pp) = actin_mem_normed ;
                        else
                            acts(pp) = results{id}.mem2cyt ;
                        end
                        assert(acts(pp) > 0)
                        
                        % Check for consistency across embryos regarding
                        % opto label and staging
                        if first
                            stage = results{id}.stage ;
                            is_opto = results{id}.is_opto ;
                        else
                            assert(stage == results{id}.stage)
                            assert(is_opto == results{id}.is_opto)
                        end
                        
                        % add to temporary plot
                        hold on;
                        if is_opto
                            plot(stage, acts(pp), '.', 'color', optoColor)
                        else
                            plot(stage, acts(pp), '.', 'color', WTcolor)
                        end
                        pause(0.0001)
    
                        % prepare for next cell
                        pp = pp + 1;
                    end
                end
            end

            % Collect into cell
            if is_opto
                actsOpto{collectOpto} = acts ;
                stagesOpto(collectOpto) = stage ;
                medianOpto(collectOpto) = median(acts) ;
                meanOpto(collectOpto) = mean(acts) ;
                collectOpto = collectOpto + 1; 
            else
                actsWT{collectWT} = acts ;
                stagesWT(collectWT) = stage ;
                medianWT(collectWT) = median(acts) ;
                meanWT(collectWT) = mean(acts) ;
                collectWT = collectWT + 1; 
            end

        end
    end

    % Plot the distribution of mem2cyt for this embryo as violin plot
    clf
    p1 = violin(actsOpto, 'x', stagesOpto, ...
        'facecolor', optoColor,'edgecolor','none',...
        'bw',bw, 'mc',[],'medc', [], 'facealpha', 0.1) ;
    hold on;
    scatter(stagesOpto, medianOpto, 50, 's', 'markeredgeColor', optoColor)
    scatter(stagesOpto, meanOpto, 50, 'o', 'filled', ...
        'markerfaceColor', optoColor, 'markeredgeColor', optoColor)
    hold on;
    p2 = violin(actsWT, 'x', stagesWT, ...
        'facecolor', WTcolor,'edgecolor','none',...
        'bw', bw, 'mc',[],'medc', [], 'facealpha', 0.1) ;
    scatter(stagesWT, medianWT, 50, 's', 'markeredgeColor', WTcolor)
    scatter(stagesWT, meanWT, 50, 'o', 'filled', ...
        'markerfaceColor', WTcolor, 'markeredgeColor', WTcolor)

    % Save the figure
    legend([p1(1), p2(1)], {optoStyle, 'WT'}, 'location', 'northwest')
    xlim([min(min(stagesWT),min(stagesOpto)) - 0.5, max(max(stagesOpto), max(stagesWT)) + 0.5])
    ylims = ylim() ;
    ylim([0, ylims(2)])
    xticks([13,14,15,16])
    
    xlabel('morphological stage', 'interpreter', 'latex')
    title('Membrane-bound phalloidin intensities', 'interpreter', 'latex')
    if normed
        ylabel('$ {I_\textrm{membrane}} / \langle{I_\textrm{embryo}}\rangle$', ...
            'interpreter', 'latex')
        saveas(gcf, fullfile(datDir, 'seg_images', 'intensity_versus_stage_embryo_normed.png'))
    else
        ylabel('$ {I_\textrm{membrane}} / \langle{I_\textrm{cytoplasm}}\rangle$', ...
            'interpreter', 'latex')
        saveas(gcf, fullfile(datDir, 'seg_images', 'intensity_versus_stage_embryo.png'))
    end
    clf
end


%% Cluster by embryo vs stages -- mem/average or  mem/cyto -- LOBES
lobeClusters = {[1,2], [3,4]} ;
lobeStrings = {'12', '34'} ;

for normed = [1, 0 ]
    if normed 
        bw = 0.01 ;
    else
        bw = 0.2 ;
    end
    
    close all ; 
    clf;

    for lobeID = 1:length(lobeClusters)
        lobes = lobeClusters{lobeID} ;
        lobestr = lobeStrings{lobeID} ;
    
        % new arrays
        actsWT = {} ;
        stagesWT = [] ;
        collectWT = 1 ;
        actsOpto = {} ;
        stagesOpto = [] ;
        collectOpto = 1 ;
        % Consider each embryo
        for eID = 1:length(embryoIDs)
            embryoID = embryoIDs{eID} ;

            % Clear acts to build up collection
            acts = [] ; 
            if ~contains(embryoID, 'embryoID')
                first = true ;
                pp = 1 ;
                % Find all images of that embryo
                for qq = 1:length(fns)
                    name = fns(qq).name ;
                    name = name(1:end-4) ;
                    disp(name)
                    outfn = fullfile(fns(qq).folder, [fns(qq).name(1:end-4) '_cells.mat']) ;

                    % Get lobe
                    lobeSplit = strsplit(name, '_L') ;
                    lobeSplit = lobeSplit{2} ;
                    lobeSplit = strsplit(lobeSplit, '_') ;
                    lobestr = lobeSplit{1} ;
                    
                    % Is this image a match to the embryoID under consideration?
                    if contains(name, embryoID) && contains(lobeStrings{lobeID}, lobestr)
                        disp(['Loading embryoID file: ' name])
                        % grab the mem2cyt ratio for this embryo in this image
                        load(outfn, 'results')
                        for id = 1:length(results)
                            disp(['cell ' num2str(id)])
                            if normed
                                adi = results{id}.actin_dist(:, 1:2) ;
                                [midx, amean, astd] = binDataMeanStd(...
                                    adi(:, 1)*pix2um, adi(:, 2), [0, memthres, max(adi(:,1))*pix2um]) ;
                                actin_mem_normed = amean(1) / bolis(eID) ;
                                acts(pp) = actin_mem_normed ;
                            else
                                acts(pp) = results{id}.mem2cyt ;
                            end
                            assert(acts(pp) > 0)

                            % Check for consistency across embryos regarding
                            % opto label and staging
                            if first
                                stage = results{id}.stage ;
                                is_opto = results{id}.is_opto ;
                            else
                                assert(stage == results{id}.stage)
                                assert(is_opto == results{id}.is_opto)
                            end

                            % prepare for next cell
                            pp = pp + 1;
                        end
                    end
                end

                % Collect into cell
                if is_opto
                    actsOpto{collectOpto} = acts ;
                    stagesOpto(collectOpto) = stage ;
                    medianOpto(collectOpto) = median(acts) ;
                    meanOpto(collectOpto) = mean(acts) ;
                    collectOpto = collectOpto + 1; 
                else
                    actsWT{collectWT} = acts ;
                    stagesWT(collectWT) = stage ;
                    medianWT(collectWT) = median(acts) ;
                    meanWT(collectWT) = mean(acts) ;
                    collectWT = collectWT + 1; 
                end

            end
        end

        % Clean away empty entries in acts,stages,mean,median
        keepO = true(length(actsOpto), 1) ;
        dmyk = 1 ;
        actsOptoNew = {} ;
        for ijk = 1:length(actsOpto)
            if isempty(actsOpto{ijk})
                keepO(ijk) = false;
            else
                actsOptoNew{dmyk} = actsOpto{ijk} ; 
                dmyk = dmyk + 1;
            end
        end
        actsOpto = actsOptoNew ; 
        stagesOpto = stagesOpto(keepO) ;
        meanOpto = meanOpto(keepO) ;
        medianOpto = medianOpto(keepO) ;
        
        % Clean away in WT
        keepW = true(length(actsWT), 1) ;
        dmyk = 1 ;
        actsWTNew = {} ;
        for ijk = 1:length(actsWT)
            if isempty(actsWT{ijk})
                keepW(ijk) = false;
            else
                actsWTNew{dmyk} = actsWT{ijk} ; 
                dmyk = dmyk + 1;
            end
        end
        actsWT = actsWTNew ; 
        stagesWT = stagesWT(keepW) ;
        meansWT = meanWT(keepW) ;
        medianWT = medianWT(keepW) ;
        
        % Plot the distribution of mem2cyt for this embryo as violin plot
        subplot(length(lobeClusters), 1, lobeID)
        p1 = violin(actsOpto, 'x', stagesOpto, ...
            'facecolor', optoColor,'edgecolor','none',...
            'bw',bw, 'mc',[],'medc', [], 'facealpha', 0.1) ;
        hold on;
        scatter(stagesOpto, medianOpto, 50, 's', 'markeredgeColor', optoColor)
        scatter(stagesOpto, meanOpto, 50, 'o', 'filled', ...
            'markerfaceColor', optoColor, 'markeredgeColor', optoColor)
        hold on;
        p2 = violin(actsWT, 'x', stagesWT, ...
            'facecolor', WTcolor,'edgecolor','none',...
            'bw', bw, 'mc',[],'medc', [], 'facealpha', 0.1) ;
        scatter(stagesWT, medianWT, 50, 's', 'markeredgeColor', WTcolor)
        scatter(stagesWT, meanWT, 50, 'o', 'filled', ...
            'markerfaceColor', WTcolor, 'markeredgeColor', WTcolor)

        % Format the figure
        legend([p1(1), p2(1)], {optoStyle, 'WT'}, 'location', 'northwest')
        xlim([12.5, 16.5])
        ylims = ylim() ;
        ylim([0, ylims(2)])
        xticks([13,14,15,16])

        xlabel('morphological stage', 'interpreter', 'latex')
        if normed
            ylabel(['$ {I_\textrm{mem}} / \langle{I_\textrm{emb}}\rangle$, L' lobeStrings{lobeID}], ...
                'interpreter', 'latex')
        else
            ylabel(['$ {I_\textrm{mem}} / \langle{I_\textrm{cyto}}\rangle$, L' lobeStrings{lobeID}], ...
                'interpreter', 'latex')
        end
    end
    
    % Save the N-panel figure
    sgtitle('Membrane-bound phalloidin intensities', 'interpreter', 'latex')
    if normed
        saveas(gcf, fullfile(datDir, 'seg_images', ['intensity_versus_stage_embryoLobes_normed.png']))
    else
        saveas(gcf, fullfile(datDir, 'seg_images', ['intensity_versus_stage_embryoLobes.png']))
    end
    clf
end

%% All results on one plot -- raw vs distance
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
        is_opto = results{id}.is_opto ;
        
        
        if is_opto
            disp('Opto')
            p1 = scatter(act(:, 1) * pix2um, act(:, 2), 4, 'markeredgealpha', 0.3, ...
                'markeredgecolor', redcolor) ;
            plot(act(:, 1) * pix2um, act(:, 2), '--', 'color', redcolor)
        else
            disp('WT')
            p2 = scatter(act(:, 1) * pix2um, act(:, 2), 8, 'markeredgealpha', 0.3, ...
                'markeredgecolor', bluecolor) ;
            plot(act(:, 1) * pix2um, act(:, 2), '--', 'color', bluecolor)
        end
        hold on;
        
        % Check that this is a unique cell
        for prevCell = 1:length(acts)
            if numel(act) == numel(acts{prevCell})
                if all(all(abs(act - acts{prevCell}) < 0.001))
                    error(['actin signal from ' num2str(prevCell) ' is same as ' num2str(id)])
                end
            end
        end
        if isempty(prevCell)
            acts{1} = act ;
        else
            acts{prevCell+1} = act ;
        end
    end
end
legend([p1, p2], {optoStyle, 'WT'})
xlabel('distance from membrane [$\mu$m]', 'interpreter', 'latex')
ylabel('Intensity [a.u.]', 'interpreter', 'latex')
saveas(gcf, fullfile(datDir, 'seg_images', 'intensity_versus_distance.png'))

%% Plot versus stage
clf
acts = [] ;
stage = [] ;
is_opto = [] ;
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
        stage(pp) = results{id}.stage ;
        optog_means(pp) = results{id}.optog_mean ;
        is_opto(pp) = results{id}.is_opto ;
        if is_opto(pp)
            p1 = scatter(stage(pp), acts(pp), 10, 'filled', 'markerfacecolor', redcolor, ...
                'markerfacealpha', 0.5) ;
        else
            p2 = scatter(stage(pp), acts(pp), 15, 'markeredgecolor', bluecolor, ...
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

