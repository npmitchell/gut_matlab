%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script to read 16 bit data images and output the mips 
% NPMitchell 2019, based on original script by SJS
%
% NOTE:
% view1 = along third dimension, near half
% view2 = along third dimension, far half
% view11 = along first dimension, near half
% view12 = along first dimension, far half
% view21 = along second dimension, near half
% view22 = along second dimension, far half
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

%% OPTIONS
overwrite_mips = true ;
mipoutdir = 'mips' ; %_stab
% tiff filenames of which to make mips 
% filenameFormat  = 'Time_%06d_c1_stab.tif';  %_stab
filenameFormat = 'TP%d_Ch0_Ill0_Ang0,60,120,180,240,300.tif' ;
scale = 300; % 0.002; % 0.02
% Offset for setting what timestep is t=0
t_off = 0;
% Where the data lives
dataDir = cd;
% which timepoints to make mips of
timePoints = 0:270;
use_scale = true ;

%% Create directories & add paths
if ~exist(mipoutdir, 'dir')
    mkdir(mipoutdir)
end

addpath('/mnt/data/code/imsaneV1/external/bfmatlab/');
cd(dataDir)
msgLevel = 1;
setpref('ImSAnE', 'msgLevel', msgLevel);

%% Make the subdirectories for the mips if not already existing
mipdirs = {fullfile(mipoutdir, 'view1/'), ...
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
% Define naming scheme for each half-volume MIP
% along dim 3
name1  = fullfile('view1', 'mip_1_%03d_c1.tif');
name2  = fullfile('view2', 'mip_2_%03d_c1.tif');
% along dim 1
name11 = fullfile('view11', 'mip_11_%03d_c1.tif');
name21 = fullfile('view21', 'mip_21_%03d_c1.tif');
% along dim 2 
name12 = fullfile('view12', 'mip_12_%03d_c1.tif');
name22 = fullfile('view22', 'mip_22_%03d_c1.tif');

%% Cycle through timepoints to make mips
for time = timePoints
    disp(['Considering time=' num2str(time)])
    fileName = sprintf(filenameFormat, time);
    fullFileName = fullfile(dataDir, fileName);
    m1fn = fullfile(mipoutdir, sprintf(name1,  time-t_off)) ;
    m2fn = fullfile(mipoutdir, sprintf(name2,  time-t_off)) ;
    m11fn = fullfile(mipoutdir, sprintf(name11, time-t_off)) ;
    m21fn = fullfile(mipoutdir, sprintf(name21, time-t_off)) ;
    m12fn = fullfile(mipoutdir, sprintf(name12, time-t_off)) ;
    m22fn = fullfile(mipoutdir, sprintf(name22, time-t_off)) ;
    mexist = exist(m1fn, 'file') && exist(m2fn, 'file') && ...
        exist(m11fn, 'file') && exist(m21fn, 'file') && ...
        exist(m12fn, 'file') && exist(m22fn, 'file') ;
    
    % Only consider this timepoint if mips don't exist, or overwrite
    if exist(fullFileName, 'file') && (~mexist || overwrite_mips)
        disp([ 'reading ' fullFileName]) 
        % data = readSingleTiff(fullFileName);
        data = bfopen(fullFileName) ;
        tmp = data{1} ;
        dstack = zeros([size(tmp{1}), length(data{1})]) ;
        for i = 1:length(data{1})
            dstack(:, :, i) = tmp{i};  
        end
        
        % Convert to grayscale at 16 bit depth
        if ~use_scale 
            scale = max(dstack(:)) ;
        end
        im2 = mat2gray(dstack,[0 scale]);
        im2 = uint16(2^16*im2);
        imSize = size(im2);

        % Creat MIPs (maximum intensity projections)
        disp(['creating mips for timepoint=' num2str(time)])
        disp(fullFileName)
        % take max intensity projections of half-space volumes
        mip_1 = max(im2(:,:,1:round(imSize(end)/2)),[],3);
        mip_2 = max(im2(:,:,round(imSize(end)/2):end),[],3);
        mip_11 = squeeze(max(im2(1:round(imSize(1)/2),:,:),[],1));
        mip_21 = squeeze(max(im2(round(imSize(1)/2):end,:,:),[],1));
        mip_12 = squeeze(max(im2(:,1:round(imSize(2)/2),:),[],2));
        mip_22 = squeeze(max(im2(:,round(imSize(2)/2):end,:),[],2));
        imwrite(mip_1, m1fn,'tiff','Compression','none');
        imwrite(mip_2, m2fn,'tiff','Compression','none');
        imwrite(mip_11,m11fn,'tiff','Compression','none');
        imwrite(mip_21,m21fn,'tiff','Compression','none');
        imwrite(mip_12,m12fn,'tiff','Compression','none');
        imwrite(mip_22,m22fn,'tiff','Compression','none');
    else
        if ~exist(fullFileName, 'file') 
            disp(['WARNING: file does not exist, skipping: ', fullFileName])
        elseif mexist
            disp(['MIPS already exist for ' fullFileName, ' -- skipping'])
        else 
            error('Somehow the file exists and mips do not, but failed to execute. Investigate here.')
        end
    end
end
disp('done')
