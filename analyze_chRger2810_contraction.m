% Look at motion of chRger2810 muscle during optogenetic stimulation

% paths for reading data
datdir = ['/mnt/data/confocal_data/gut/chRger_mef2G4mcD8klar/', ...
    '202202071508_chRger2810mef2G4mcD8klar_5pc5pc_20pcArgon_0p5spf/', ...
    'e2_15a/chRger_posteriorFold_inset'] ;
cd(datdir)
fns = dir(fullfile(datdir, '*Probabilities.h5'));
imfnsBF = dir(fullfile(datdir, 'separated_channels', 'chRger_posteriorFold_inset_t*_c003.png')) ;
imfnsGC = dir(fullfile(datdir, 'separated_channels', 'chRger_posteriorFold_inset_t*_c001.png')) ;
pix2um = 1/3.5200 ; % um / pix
colors = define_colors ;
yellow = colors(3, :) ;
cyan = [0, 1, 1] ;

tt = 0 ;
dt = 0.5 ;
dts = [120, 120, 120, 0] ;

dmyk = 1 ;
for qq = 1:length(fns)
    prob = h5read(fullfile(fns(qq).folder, fns(qq).name), ...
        '/exported_data') ;
    coms = zeros(size(prob, 4), 2) ;
    tt0 = tt ;
    for tidx = 1:size(prob, 4) 
        page = prob(:, :, 1, tidx) ;
        page2 = page(200:end, :) ;
        coms(tidx, :) = com_region2d(page2) + [200, 0] ;
        tt = tt + dt ;
        timestamp(tidx) = tt ;
    end
    coms_seq{qq} = coms ; 
    timestamps{qq} = timestamp ;
    
    % Inspect result for all images in seq
    disp('overlaying images...')
    figdir = fullfile(datdir, 'overlayTracking') ;
    ensureDir(figdir) 
    h1 = figure(1) ;
    for tidx = 1:size(prob, 4)
        clf ;
        set(h1, 'visible', 'off')
        im = imread(fullfile(imfnsBF(dmyk).folder, imfnsBF(dmyk).name)) ;
        im2 = imread(fullfile(imfnsGC(dmyk).folder, imfnsGC(dmyk).name)) ;
        imO = imoverlay(im, im2 > 50, cyan) ;
        imshow(imO)
        set(h1, 'visible', 'off')
        hold on;
        scatter(coms(tidx, 1), coms(tidx, 2), 30, 'filled', 'markerfacecolor', yellow);
        title([' t = ' sec2hms(floor(timestamp(tidx)), '%02d:%02d:%02d')])
        outfn = fullfile(figdir, imfnsBF(dmyk).name) ;
        export_fig(outfn)
        dmyk = dmyk + 1 ;
    end
    
    % Prepare for next sequence
    tt = tt0 + dts(qq) ;
end

%% Plot all data now
clearvars recovery
clf
finalpt = zeros(length(fns), 2) ;
initialpt = zeros(length(fns), 2) ;
dist = zeros(length(fns), 1) ;
recovery = zeros(length(fns)-1, 1) ;

figure('Units', 'centimeters', 'position', [0, 0, 9,6]);

initPos = coms_seq{1}(1, :) ;
for qq = 1:length(fns)
    plot(timestamps{qq}, -(coms_seq{qq}(:, 2) - initPos(2)) * pix2um, '.-')
    xlabel('time [s]', 'Interpreter', 'latex')
    ylabel('DV contraction [$\mu$m]', 'Interpreter', 'latex')
    hold on;
    
    finalpt(qq, :) = coms_seq{qq}(end, :)  ;
    initialpt(qq, :) = coms_seq{qq}(1, :) ;
    distV = finalpt(qq, :) - initialpt(qq, :) ;
    dist(qq) = sqrt(sum(distV.^2)) ;
    
    %  o---------------------o
    %                  o<-----
    if qq > 1
        distFromPrevInit = initialpt(qq, :) - initialpt(qq-1, :) ;
        distFromPrevInit = sqrt(sum(distFromPrevInit.^2)) ;
        % distFromPrevFinal = initialpt(qq, :) - finalpt(qq-1, :) ;
        recovery(qq-1) = 1 - distFromPrevInit / dist(qq-1) ;
    end
end
saveas(gcf, fullfile(datdir, 'displacement_during_activation.png'))
saveas(gcf, fullfile(datdir, 'displacement_during_activation.pdf'))

%% Plot recovery
cla
bar(recovery)
ylim([0, 1])
xlabel('strobe iteration')
ylabel('fractional recovery')
saveas(gcf, fullfile(datdir, 'recovery_after_activation.png'))
saveas(gcf, fullfile(datdir, 'recovery_after_activation.pdf'))

