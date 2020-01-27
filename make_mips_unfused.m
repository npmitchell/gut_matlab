%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script to read 16 bit TIFF images (raw, not fused) and output the mips 
% NPMitchell 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% Main options
% whether to overwrite the mip images
overwrite = true ;  
% what are the files called, with time and angle as variable
filenameFormat  = 'Time_%06d_Angle_%d_c3_ls_1.ome.tif';  
% which timepoints to make
timePoints      = [0:87];
% which angles to make
angles = [0, 45, 90, 135, 180, 225, 270, 315] ;
% Maximum intensity cutoff
scale = 1000 ;
start_frac = 0.3 ; % what fraction of volume to start at
stop_early = 10 ; % how many frames 

%% Directory setup
mipoutdir = 'mips_raw' ;
if ~exist(mipoutdir, 'dir')
    mkdir(mipoutdir)
end

addpath('/mnt/data/code');
addpath('/mnt/data/code/imsaneV1/external/bfmatlab/');
dataDir    = cd;
cd(dataDir)

%% Other options
msgLevel = 1;
setpref('ImSAnE', 'msgLevel', msgLevel);

% Make the subdirectories for the mips if not already existing
for kk = 1:length(angles)
    angle = angles(kk) ;
    mipdirs{kk} = fullfile(mipoutdir, sprintf('angle%03d', angle)) ;
end
for i = 1:length(mipdirs)
    if ~exist(mipdirs{i},'dir')
        mkdir(mipdirs{i})
    end
end

% Consider each angle in angles, make mip of that view
for kk = 1:length(angles)
    angle = angles(kk) ;
    % Consider each timepoint in timePoints and make mip of that time in
    % this current view angle.
    for time = timePoints
        disp(['Considering time=' num2str(time)])
        fileName = sprintf(filenameFormat, time, angle);
        fullFileName = fullfile(dataDir, fileName);

        if exist(fullFileName, 'file') || overwrite
            disp([ 'reading ' fullFileName]) 
            % data = readSingleTiff(fullFileName);
            data = bfopen(fullFileName) ;
            tmp = data{1} ;
            dstack = zeros([size(tmp{1}), length(data{1})]) ;
            for i = 1:length(data{1})
                dstack(:, :, i) = tmp{i};  
            end

            % Convert to grayscale at 16 bit depth
            im2 = mat2gray(dstack,[0 scale]);
            im2 = uint16(2^16*im2);
            imSize = size(im2);

            % Creat MIPs (maximum intensity projections)
            disp(['creating mips for timepoint=' num2str(time)])
            disp(fullFileName)
            mip_1 = max(im2(:,:,round(imSize(end)* start_frac):end-stop_early),[],3);
%             mip_2 = max(im2(:,:,round(imSize(end)/2):end),[],3);
%             mip_11 = squeeze(max(im2(1:round(imSize(1)/2),:,:),[],1));
%             mip_21 = squeeze(max(im2(round(imSize(1)/2):end,:,:),[],1));
%             mip_12 = squeeze(max(im2(:,1:round(imSize(2)/2),:),[],2));
%             mip_22 = squeeze(max(im2(:,round(imSize(2)/2):end,:),[],2));

            imwrite(mip_1, fullfile(mipdirs{kk}, sprintf('mip_1_%03d_c1.tif',  time)),'tiff','Compression','none');
%             imwrite(mip_2, fullfile(mipoutdir, sprintf('view2/mip_2_%03d_c1.tif',  time-t_off)),'tiff','Compression','none');
%             imwrite(mip_11,fullfile(mipoutdir, sprintf('view11/mip_11_%03d_c1.tif',time-t_off)),'tiff','Compression','none');
%             imwrite(mip_21,fullfile(mipoutdir, sprintf('view21/mip_21_%03d_c1.tif',time-t_off)),'tiff','Compression','none');
%             imwrite(mip_12,fullfile(mipoutdir, sprintf('view12/mip_12_%03d_c1.tif',time-t_off)),'tiff','Compression','none');
%             imwrite(mip_22,fullfile(mipoutdir, sprintf('view22/mip_22_%03d_c1.tif',time-t_off)),'tiff','Compression','none');
        else
            disp(['WARNING: file does not exist, skipping: ', fullFileName])
        end
    end
    disp(['Done with angle ' num2str(angle)])
end
disp('done')