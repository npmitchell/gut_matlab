%% Example for calibration of stabilizeImages.m
% Run this then stabilizeImagesCorrect.m to see the calibration on the
% example data created by this script.
% >> stabilizeImagesExampleScript
% >> stabilizeImagesCorrect

%
% NOTE:
% view1 = along third dimension, near half
% view2 = along third dimension, far half
% view11 = along first dimension, near half
% view12 = along first dimension, far half
% view21 = along second dimension, near half
% view22 = along second dimension, far half

%% Define/add paths
mipoutdir = 'mips' ;
addpath('/mnt/data/code/imsaneV1/external/bfmatlab/');
dataDir = './';

%% Make the subdirectories for the mips if not already existing
mipdirs = {mipoutdir, ...
    fullfile(mipoutdir, 'view1/'), ...
    fullfile(mipoutdir, 'view2/'), ...
    fullfile(mipoutdir, 'view11/'), ...
    fullfile(mipoutdir, 'view12/'),...
    fullfile(mipoutdir, 'view21/'),...
    fullfile(mipoutdir, 'view22/')} ;
for i = 1:length(mipdirs)
    if ~exist(mipdirs{i},'dir')
        mkdir(mipdirs{i})
    end
end

%% Make the data
N = 25 ;
W = 350 ;
stepX = 1 ;
stepY = 2 ;
t_off = 0 ;
x0 = 50 ;
y0 = 50 ;
tmp = load('spiralVol.mat');
im = double(tmp.spiralVol) ;
filenameFormat  = 'Time_%06d_c1.tif';  
for time = 1:(N-1)
    disp(['Considering time=' num2str(time)])
    
    % create a 3d volume of data
    % dstack = 0.1 * rand(N, N, N) ;
    % y = min(N, max(round(time*0.50 + N*0.5), 1)) ;
    % z = min(N, max(round(time*0.25 + N*0.25), 1)) ;
    xind = x0 + (time*stepX):(time*stepX+W) ;
    yind = y0 + (time*stepY):(time*stepY+W) ;
    dstack = im(xind, yind, :) ;
    % dstack(1:N, 1:N, 1:N) = 0 ;
    % dstack(end:end-N, end:end-N, end:end-N) = 0 ;
    % Convert to grayscale at 16 bit depth
    im2 = mat2gray(dstack,[0 max(dstack(:))]);
    im2 = uint16(2^16*im2);
    imSize = size(im2);
    
    % Define data filename
    fileName = sprintf(filenameFormat, time);
    fullFileName = fullfile(dataDir, fileName);
    
    % Write to disk. Note that if uint16 specified upon instantiation, then
    % imwrite will respect the class type.
    for z = 1:size(im2, 3)
        imwrite(im2(:,:,z), fullFileName,'tiff',...
            'Compression','none','WriteMode','append') ;
    end

    % Define mip filenames
    m1fn = fullfile(mipoutdir, sprintf('view1/mip_1_%03d_c1.tif',  time-t_off)) ;
    m2fn = fullfile(mipoutdir, sprintf('view2/mip_2_%03d_c1.tif',  time-t_off)) ;
    m11fn = fullfile(mipoutdir, sprintf('view11/mip_11_%03d_c1.tif',time-t_off)) ;
    m21fn = fullfile(mipoutdir, sprintf('view21/mip_21_%03d_c1.tif',time-t_off)) ;
    m12fn = fullfile(mipoutdir, sprintf('view12/mip_12_%03d_c1.tif',time-t_off)) ;
    m22fn = fullfile(mipoutdir, sprintf('view22/mip_22_%03d_c1.tif',time-t_off)) ;
    mexist = exist(m1fn, 'file') && exist(m2fn, 'file') && ...
        exist(m11fn, 'file') && exist(m21fn, 'file') && ...
        exist(m12fn, 'file') && exist(m22fn, 'file') ;

    % Creat MIPs (maximum intensity projections)
    disp(['creating mips for timepoint=' num2str(time)])
    disp(fullFileName)

    mip_1 = max(im2(:,:,1:round(imSize(end)/2)),[],3);
    mip_2 = max(im2(:,:,round(imSize(end)/2):end),[],3);
    mip_11 = squeeze(max(im2(1:round(imSize(1)/2),:,:),[],1));
    mip_21 = squeeze(max(im2(round(imSize(1)/2):end,:,:),[],1));
    mip_12 = squeeze(max(im2(:,1:round(imSize(2)/2),:),[],2));
    mip_22 = squeeze(max(im2(:,round(imSize(2)/2):end,:),[],2));

    % Write the mip files
    imwrite(mip_1, m1fn, 'tiff','Compression','none');
    imwrite(mip_2, m2fn, 'tiff','Compression','none');
    imwrite(mip_11,m11fn,'tiff','Compression','none');
    imwrite(mip_21,m21fn,'tiff','Compression','none');
    imwrite(mip_12,m12fn,'tiff','Compression','none');
    imwrite(mip_22,m22fn,'tiff','Compression','none');
    
end
