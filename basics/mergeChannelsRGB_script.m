% Merge images to RGB for imagestack
im1s = dir('./endoderm_luminance/*.png') ;
im2s = dir('./muscle_luminance/*.png') ;
timePoints = 0:59 ;
timeInterval = 2 ;
outDir = './combined_blue_yellow/' ;

color1 = [100, 200, 255] / 255 ;  % orig: [123, 205, 240] / 255 ;
color2 = [255, 200, 100] / 255 ;  % orig: [255, 234, 157] / 255 

% area to merge
xymin = [204, 169] ;
xymax = [640, 1050] ;

% area to mask
maskmin = [140, 525] ;
maskmax = [195, 700] ;
maskTextXY = [mean(maskmin(1), maskmax(1)), mean(maskmin(2), maskmax(2))];

iptsetpref('ImshowBorder','tight');
for ii = 1:length(im1s)
    % Read in the images to merge/combine
    im1 = imread(fullfile(im1s(ii).folder, im1s(ii).name)) ;
    im2 = imread(fullfile(im2s(ii).folder, im2s(ii).name)) ;
    
    % Rescale by clim
    im1Sub = im1(xymin(1):xymax(1), xymin(2):xymax(2)) ;
    im2Sub = im2(xymin(1):xymax(1), xymin(2):xymax(2)) ;
    val1 = mean(double(im1Sub(:))) ;
    val2 = mean(double(im1Sub(:))) ;
    if ii == 1
        val01 = val1 ;
        val02 = val2 ;
    end
    
    im1r = mat2gray(im1, [0, 130 * val1 / val01]) ;
    im2r = mat2gray(im2, [0, 100 * val2 / val02]) ;
    im1 = mat2gray(im1, [0, 255]) ;
    im2 = mat2gray(im2, [0, 255]) ;
    
    % Put merged image in center of canvas (RGB section)
    outim = zeros(size(im1, 1), size(im1, 2), 3) ;
    outMerge = mergeChannelsRGB(im1r, im2r, color1, color2, 8) ;
    outim(xymin(1):xymax(1), xymin(2):xymax(2), :) = ...
        outMerge(xymin(1):xymax(1), xymin(2):xymax(2), :) ;
    
    % Border of canvas is just grayscale of original image
    stub = im2(:, 1:xymin(2)) ;
    outim(:, 1:xymin(2), :) = cat(3, stub,stub,stub) ;
    stub = im2(:, xymax(2):end) ;
    outim(:, xymax(2):end, :) = cat(3, stub,stub,stub) ;
    stub = im2(1:xymin(1), :) ;
    outim(1:xymin(1), :, :) = cat(3, stub,stub,stub) ;
    stub = im2(xymax(1):end, :) ;
    outim(xymax(1):end, :, :) = cat(3, stub,stub,stub) ;
    
    % Mask out mask area
    outim(maskmin(1):maskmax(1), maskmin(2):maskmax(2), :) = 0 ;
    
    % Add title
    % outim = insertText(outim, maskTextXY, sprintf('$t= %d$ min', timePoints(ii) * timeInterval)), 'interpreter', 'latex')
    close all
    imshow(outim)
    hold on;
    text(maskTextXY(2), maskTextXY(1), sprintf('$t= %d$ min', timePoints(ii) * timeInterval), ...
        'interpreter', 'latex', 'color', 'w', 'fontsize', 16, ...
        'horizontalalignment', 'left', 'verticalalignment', 'top')
    
    
    % Output the merged/combined image to disk
    % write canvas to disk
    outfn = fullfile(outDir, 'canvas', ...
        sprintf('combined_%06d_canvas.png', timePoints(ii))); 
    if ~exist( fullfile(outDir, 'canvas'), 'dir')
        mkdir( fullfile(outDir, 'canvas'))
    end
    if ~exist( fullfile(outDir, 'image'), 'dir')
        mkdir( fullfile(outDir, 'image'))
    end
    daspect([1,1,1])
    drawnow
    
    H = getframe(gcf) ;
    imwrite(H.cdata, outfn)
    
    % Write image only
    outfn = fullfile(outDir, 'image', ...
        sprintf('combined_%06d_image.png', timePoints(ii))); 
    imwrite(outim, outfn) 
end


%% Add title from canvas to image
if ~exist(fullfile(outDir, 'image_with_header'), 'dir')
    mkdir(fullfile(outDir, 'image_with_header'))
end

ims = dir(fullfile(outDir, 'image', 'combined*.png')) ;
tims = dir(fullfile(outDir, 'canvas', 'combined*.png')) ;
for ii = 1:length(ims)
    % Read in the images to merge/combine
    im = imread(fullfile(ims(ii).folder, ims(ii).name)) ;
    tim = imread(fullfile(tims(ii).folder, tims(ii).name)) ;
    
    im(maskmin(1):maskmax(1), maskmin(2):maskmax(2), :) = ...
        tim(maskmin(1):maskmax(1), maskmin(2):maskmax(2), :) ;
    
    outfn = fullfile(outDir, 'image_with_header', ...
        sprintf('combined_%06d_image.png', timePoints(ii))); 
    imwrite(im, outfn) 
end

