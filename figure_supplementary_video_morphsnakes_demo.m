% Load and label tiff made by:
% options = struct() ;
% options.plot_evolution = false ;
% options.growth_t0 = 85 ;
% QS.visualizeMeshEvolution(options)


fn0 = '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/msls_output/demo_morphsnakes_figs/'; 
fn = [fn0 'GrowEvolve_CombinedStacks_085_t0_3panel_labeled.tif'];
outdir = fullfile(fn0, 'CombinedStacks_for_label') ;

dat = dir(fullfile(outdir, 'preLabel/GrowEvolve*.png')) ;
nFrames = length(dat) ;
% dat = readTiff4D(fn, 1) ;
% size(dat) ;
% nFrames = size(dat, 4) ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

timePoints = zeros(nFrames, 1) ;
timePoints(81:103) = -38 ;
timePoints(103:end) = -38:140 ;
for frame = 81:nFrames
    disp(['frame = ', num2str(frame)])
    outfn = fullfile(outdir, sprintf('%06d.png', frame)) ;
    % im = squeeze(dat(:, :, :, frame)) ;
    im = imread(fullfile(dat(frame).folder, dat(frame).name)) ;
    
    % What to label
    position = [500, 500] ;
    stringTitle = ['$t = $', sprintf('%03d', timePoints(frame)) ' min'] ;
    figure('visible', 'off');
    plot([0, 1, 1, 0], [0, 0, 1, 1], 'k-')
    text(0.5, 0.5, stringTitle,'Interpreter','latex', 'fontsize', 48, 'Horizontalalignment', 'center')
    FF = getframe(gca);
    dF = FF.cdata ;
    dF = dF(120:end-120, 50:end-50, :) ;
    % dF = uint8(( dF0 / 255) .^ 1.2 * 255)
    % imshow(dF)
    dFs = imresize(dF, 0.37, 'bilinear');
    
    dFX = size(dFs, 1) ;
    dFY = size(dFs, 2) ;
    xmin = 125 ;
    ymin = 850 ;
    
    im2 = im ;
    im2(xmin+1:xmin+dFX, ymin+1:ymin+dFY, :) = dFs ;
    % imshow(im2)
    imwrite(im2, outfn)
    
    close all
end

disp('done')