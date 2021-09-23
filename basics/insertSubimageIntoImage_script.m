% Merge images to RGB for imagestack
im1s = dir('./DVLR/*.png') ;
im1s = dir('./perspective_upsampled/*.png') ;
im2s = dir('./timestamper/*.png') ;
timePoints = 0:169 ;
% outDir = './DVLR_withHeader' ;
outDir = './perspective_withHeader' ;

% area to mask with subimage
Xc = 514 ;
Yc = 10 ;

%% Add subimage onto other image
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

for ii = 1:length(im1s)
    % Read in the images to merge/combine
    im1 = imread(fullfile(im1s(ii).folder, im1s(ii).name)) ;
    im2 = imread(fullfile(im2s(ii).folder, im2s(ii).name)) ;
    
    im1(Yc:Yc+size(im2, 1)-1, Xc:Xc+size(im2, 2)-1, :) = im2 ;
    outfn = fullfile(outDir, ...
        sprintf('combined_%06d_image.png', timePoints(ii))); 
    imwrite(im1, outfn) 
end


%% Add transition -- ENDING
Ntrans = 30 ;
infn = fullfile(outDir, ...
        sprintf('combined_%06d_image.png', timePoints(end))) ;
im1 = imread(infn) ;
for qq = 1:Ntrans
    
    im2 = uint8(im1 * (1 - qq / Ntrans)) ;
    outfn = fullfile(outDir, ...
        sprintf('combined_%06d_image.png', timePoints(end)+qq));
    imwrite(im2, outfn) 
end



%% Add transition -- Beginning
Ntrans = 30 ;
infn = fullfile(outDir, ...
        sprintf('combined_%06d_image.png', timePoints(1))) ;
im1 = imread(infn) ;
for qq = 1:Ntrans
    
    im2 = uint8(im1 * (1 - qq / Ntrans)) ;
    outfn = fullfile(outDir, ...
        sprintf('combined_%06d_image.png', timePoints(1) - Ntrans -1 + qq));
    imwrite(im2, outfn) 
end
