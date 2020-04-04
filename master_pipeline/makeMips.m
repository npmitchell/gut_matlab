function makeMips(timePoints, dir16bit, fileName, mipDir, Options)
% MAKEMIPS(timePoints, dir16bit, fileName, mipDir, Options)
% read 16 bit data images and output the mips 
%
% Parameters
% ----------
%  
%
% Outputs
% -------
%
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

overwrite_mips = true ;
use_scale = false ; 
if isfield(Options, 'overwrite_mips')
    overwrite_mips = Options.overwrite_mips ;
end
if isfield(Options, 'scale')
    if Options.scale > 0 
        use_scale = true ;
        scale = Options.scale ;
    end
end

msgLevel = 1;
setpref('ImSAnE', 'msgLevel', msgLevel);

%% Make the subdirectories for the mips if not already existing
mipdirs = {mipDir, fullfile(mipDir, 'view1/'), ...
    fullfile(mipDir, 'view2/'), ...
    fullfile(mipDir, 'view11/'), ...
    fullfile(mipDir, 'view12/'),...
    fullfile(mipDir, 'view21/'),...
    fullfile(mipDir, 'view22/')} ;
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
    fullFileName = fullfile(dir16bit, sprintf(fileName, time)) ;
    m1fn = fullfile(mipDir, sprintf(name1,  time)) ;
    m2fn = fullfile(mipDir, sprintf(name2,  time)) ;
    m11fn = fullfile(mipDir, sprintf(name11, time)) ;
    m21fn = fullfile(mipDir, sprintf(name21, time)) ;
    m12fn = fullfile(mipDir, sprintf(name12, time)) ;
    m22fn = fullfile(mipDir, sprintf(name22, time)) ;
    mexist = exist(m1fn, 'file') && exist(m2fn, 'file') && ...
        exist(m11fn, 'file') && exist(m21fn, 'file') && ...
        exist(m12fn, 'file') && exist(m22fn, 'file') ;
    
    if mexist
        disp(['MIPs already exist on disk for t=' num2str(time)])
        if exist(fullFileName, 'file') && overwrite_mips
            disp('Overwriting...')
        elseif ~exist(fullFileName, 'file') && overwrite_mips
            disp('cannot overwrite since tif does not exist')
        end
    end
    
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
        im2 = mat2gray(dstack, [0 scale]);
        im2 = uint16(2^16 * im2);
        imSize = size(im2);

        % Creat MIPs (maximum intensity projections)
        disp(['creating mips for timepoint=' num2str(time)])
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

        % declare done
        clearvars im2
        disp(['finished mips for ' fullFileName])
    else
        if ~exist(fullFileName, 'file') 
            disp(['WARNING: file does not exist, skipping: ', fullFileName])
        elseif mexist
            disp(['MIPs skipped for ' fullFileName, ' -- skipping'])
        else 
            error('Somehow the file exists and mips do not, but failed to execute. Investigate here.')
        end
    end
end
disp('done')
