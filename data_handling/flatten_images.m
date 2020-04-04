

% datadir = '/mnt/crunch/embryo/48Ygal4pmCIBNGFP/201901061300_goodfold/Time12views_60sec_1.2um_25x_1/data/';
% datadir = '/mnt/crunch/embryo/48Ygal4-UAShistRFP/201901031740_folded_interruptoutofmem/Time12views_60sec_1.2um_25x_1/data/';
datadir = '/mnt/crunch/embryo/48Ygal4-UAShistRFP/201901021550_folded_2part/Time12views_60sec_1.2um_25x_3/data/' ;
files = dir(fullfile(datadir,'TP*.tif'));
outsubdir = 'flattened_imagestack_deconv/' ;

datadir = '/mnt/crunch/embryo/48Ygal4-UAShistRFP/201901021550_folded_2part/Time12views_60sec_1.2um_25x_4/data/' ;
files = dir(fullfile(datadir,'Time*_Angle_0_c1_ls_1.ome.tif'));
outsubdir = 'flattened_imagestack_angle0_90to160/';

% for deconv: (20, 390). For regular, (80, 120)
aind = 90
zind = 160 ; % 120

% load each tif image in a directory and flatten the array along the z axis
first = true ;
for kk = 1:numel(files)
    disp([num2str(kk) '/' num2str(numel(files))])
    imfn = fullfile(datadir, files(kk).name)
    %imfinfo(imfn)
    
    % build up 3d image
    % read one subimage
    im = imread(imfn, aind);
    dims = size(im) ;
    % imflat = zeros(dims(1), dims(2)) ;  % , zind-aind) ;
    for jj = 1:(zind - aind)
        disp([' > building image', num2str(kk), ':', num2str(aind+jj)])
        im = imread(imfn, aind + jj) ;
        if jj == 1
            imflat = im ;
        else
            imflat = imflat + im ;
        end
    end
    % imflat = sum(im, 3) ;
    % imagesc(imflat)
    
    if first
        imflat = im2double(imflat);
        scale = mean(imflat(:)) + 8 * std(imflat(:)) ;
        % scale = max(imflat(:))
        first = false ;
    end
    imflat = im2uint16(im2double(imflat) / scale) ;
    
    % max(imflat(:))
    % imagesc(imflat)
    % Save image 
    outname = fullfile([datadir outsubdir], ['T' sprintf('%06d', kk - 1) '.tif'])
    imwrite(imflat, outname)
end