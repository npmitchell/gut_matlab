function applyStabilizationToChannel(stabChannel, ...
    fileName, fileNameOut, rgbName, typename, ...
    timePoints, times_todo, t_ref, mipDir, mipStabDir, mipsRGBDir, Options)
% STABILIZEIMAGES()
% Stabilize images of a single channel by removing jitter taken from MIPs
% of that channel.
% 
% We compute the shifts to stabilize all volumes and write the 16
% bit stabilized volumes to disk. We generate the stabilized MIPs and show
% an RGB overlay of the reference timepoint data and stabilized data in
% cyan/magenta.
% 
% NOTE:
% view1 = along third dimension, near half
% view2 = along third dimension, far half
% view11 = along first dimension, near half
% view12 = along first dimension, far half
% view21 = along second dimension, near half
% view22 = along second dimension, far half
%
% Requirements
% ------------
% MIPs of the data specified by fileName must exist in:
%    mipDir/view1/mip_1_%03d_c1.tif
%    mipDir/view2/mip_2_%03d_c1.tif
%    mipDir/view11/mip_11_%03d_c1.tif
%    mipDir/view21/mip_21_%03d_c1.tif
%    mipDir/view12/mip_12_%03d_c1.tif
%    mipDir/view22/mip_22_%03d_c1.tif
% This is accomplished by running make_mips.m
%
% Parameters
% ----------
% mipDir : where MIPs are stored
% mips_stab_check : output dir for the stabilized RGB overlays 
% mipoutdir : output dir for the stabilized MIPs
% im_intensity : scale for MIP output intensity
% imref_intensity : scale for reference MIP intensity in RGB overlay
% t_ref : timestamp of reference, matches other timepoints to this one
% timePoints : all timestamps to consider in the sequence
% times_todo : timestamps to correct/make output
% typename : output volume datatype, as string ('uint8', 'uint16')
% fileName : input filename format string for data to stabilize
% fileNameOut : output filename format string for stabilized volumes
% rgbName : RGB overlay filename format string
% Options : struct with fields
%   stabChannel
%
% Returns
% -------
% - jitter/drift-corrected 3D volumes 
% - MIPs of the jitter/drift-corrected volumes within mipoutdir
% - colored overlays of one view with the reference time. 
%
% See also
% --------
% stabilizeImagesCorrect.m -- script version of this function
% 
% NPMitchell 2020

%% Default options
if count(fileName, '%') == 2
    timechannel = 'tc' ;
elseif count(fileName, '%') == 1
    timechannel = 't' ;
else
    error('handle case here')
end

%% Unpack Options
im_intensity = Options.im_intensity ;
imref_intensity = Options.imref_intensity ;
overwrite_mips = Options.overwrite_mips ;
overwrite_tiffs = Options.overwrite_tiffs ;

%% Make the subdirectories for the mips if not already existing
mipdirs = {mipStabDir, mipsRGBDir, ...
    fullfile(mipStabDir, 'view1/'), ...
    fullfile(mipStabDir, 'view2/'), ...
    fullfile(mipStabDir, 'view11/'), ...
    fullfile(mipStabDir, 'view12/'),...
    fullfile(mipStabDir, 'view21/'),...
    fullfile(mipStabDir, 'view22/')} ;
for i = 1:length(mipdirs)
    if ~exist(mipdirs{i},'dir')
        mkdir(mipdirs{i})
    end
end
name1  = fullfile('view1', 'mip_1_%03d_c%d.tif');
name2  = fullfile('view2', 'mip_2_%03d_c%d.tif');
name11 = fullfile('view11', 'mip_11_%03d_c%d.tif');
name21 = fullfile('view21', 'mip_21_%03d_c%d.tif');
name12 = fullfile('view12', 'mip_12_%03d_c%d.tif');
name22 = fullfile('view22', 'mip_22_%03d_c%d.tif');
t_ref_ind = find( timePoints == t_ref ) ;

%% Define shifts for each time point ======================================
disp('Defining shifts...')
shiftfn = fullfile(mipDir, 'shifts_stab.mat') ;
if exist(shiftfn, 'file')
    disp('Loading shifts from disk')
    load(shiftfn, 'shifts', 't_ref', 't_ref_ind')
    x_1 = cat(1,shifts.x_1); % rows        (x) 
    y_1 = cat(1,shifts.y_1); % columns     (y)
    x_2 = cat(1,shifts.x_2); % columns #2  (y)
    y_2 = cat(1,shifts.y_2); % leaves      (z) 
    x_3 = cat(1,shifts.x_3); % rows    #2  (x) 
    y_3 = cat(1,shifts.y_3); % leaves  #2  (z)   
    % Average contributions from different views on the same shift axis
    dx = round(0.5 * (x_1 + x_3)) ;
    dy = round(0.5 * (y_1 + x_2)) ;
    dz = round(0.5 * (y_2 + y_3)) ;
else
    error('no shifts for stabilization on disk!')
end

%% Clear the individual shift values
clearvars x_1 y_1 x_2 y_2 x_3 y_3
close('all')
% Define stackSize, which is the number of frames in the z dimension
done = false ;
stackSize = 0 ;
if strcmp(timechannel, 'tc')
    name_ref = sprintf(fileName, t_ref, stabChannel) ;
elseif strcmp(timechannel, 't')
    name_ref = sprintf(fileName, t_ref) ;
else
    error('handle case here')
end
while ~done
    stackSize = stackSize + 1 ;
    try 
        tmp = imread(name_ref, stackSize);
    catch
        done = true ;
    end
end    
stackSize = stackSize - 1 ;

%% Build reference MIP for RGB overlay
disp('building reference MIP...')
if strcmp(timechannel, 't')
    name_ref = sprintf(fileName, t_ref);
elseif strcmp(timechannel, 'tc')
    name_ref = sprintf(fileName, t_ref, stabChannel);
else
    error('handle case here')
end
% preallocate im_ref3D for speed
tmp = imread(name_ref, 1);
im_ref3D = zeros([size(tmp) stackSize]) ;
for z = 1 : stackSize
    im_ref3D(:,:,z)= imread(name_ref,z);
end
mip_ref = squeeze(max(im_ref3D,[],3));
clear tmp
disp('done creating reference MIP')

%% Build MIP image for each stabilized timepoint ==========================
disp('Running through timepoints to build ims...')
tidx_todoA = 1:10:length(timePoints) ;
tidx_todoB = setdiff(1:length(timePoints), tidx_todoA) ;
tidx_todo = [tidx_todoA, tidx_todoB] ;
for tid = tidx_todo
    disp(['considering tidx = ' num2str(tid)])
    time = timePoints(tid);
    if ismember(time, times_todo)
        % The original image in 3d
        if strcmp(timechannel, 't')
            im0fn = sprintf(fileName, time);
            name_out = sprintf(fileNameOut, time) ;
            rgb_outname = sprintf(rgbName, time) ;
        elseif strcmp(timechannel, 'tc')
            im0fn = sprintf(fileName, time, stabChannel);    
            name_out = sprintf(fileNameOut, time, stabChannel) ;
            rgb_outname = sprintf(rgbName, time, stabChannel) ;
        else
            error('handle case here')
        end
        
        % Check that we're not appending to an existing file
        tiff_exists = exist(name_out, 'file') ;
        if tiff_exists && ~overwrite_tiffs && ~overwrite_mips
            disp(['Output file already exists: ' name_out ])
        else
            disp('Creating mips for this timepoint')
            % preallocate im_2_3D for speed, and 
            % Specify typename for correct output bitdepth
            tmp = imread(im0fn, 1) ;
            im0 = zeros([size(tmp) stackSize], typename) ;
            clear tmp

            for z = 1:stackSize
                im0(:,:,z)= imread(im0fn,z);
            end

            % Initialize the image 
            im = im0 ;
            imx = 0 * im0 ;  % for shift in x
            imy = 0 * im0 ;  % for shift in y

            % Offset if there is a shift in X
            if dx(tid)~=0
                if dx(tid)>0
                    % imx(dx(tid):end,:,:) = im(1:end-dx(tid)+1,:,:);
                    imx(1+dx(tid):end,:,:) = im(1:end-dx(tid),:,:);
                else
                    % imx(1:(end+dx(tid)+1),:,:) = im(-dx(tid):end,:,:);
                    imx(1:(end+dx(tid)),:,:) = im((1-dx(tid)):end,:,:);
                end
            else
                imx = im ;
            end

            % Offset if there is a shift in Y
            if dy(tid)~=0
                if dy(tid)>0
                    % imy(:,dy(tid):end,:) = imx(:,1:(end-dy(tid)+1),:);
                    imy(:,1+dy(tid):end,:) = imx(:,1:(end-dy(tid)),:);
                else
                    % imy(:,1:(end+dy(tid)+1),:) = imx(:,-dy(tid):end,:);
                    imy(:,1:(end+dy(tid)),:) = imx(:,(1-dy(tid)):end,:);
                end
            else
                imy = imx ;
            end

            % Offset if there is a shift in Z
            if dz(tid)~=0
                if dz(tid)>0
                    %im(:,:,dz(tid):end) = imy(:,:,1:end-dz(tid)+1);
                    im(:,:,1+dz(tid):end) = imy(:,:,1:end-dz(tid));
                else
                    % im(:,:,1:(end+dz(tid)+1)) = imy(:,:,-dz(tid):end);
                    im(:,:,1:(end+dz(tid))) = imy(:,:,(1-dz(tid)):end);
                end
            else
                im = imy ;
            end

            % Write to disk. Note that if uint16 specified upon instantiation, then
            % imwrite will respect the class type.
            if ~tiff_exists || overwrite_tiffs
                disp('Saving stabilized TIFF')
                for z = 1:stackSize
                    if z == 1
                        imwrite(im(:,:,z), name_out,'tiff',...
                            'Compression','none') ;    
                    else
                        imwrite(im(:,:,z), name_out,'tiff',...
                            'Compression','none','WriteMode','append') ;
                    end
                end
                disp(['saved image ' name_out ])
            end
            
            % Make a color of the current mip wrt the reference
            mip_1 = squeeze(max(im,[],3));
            if strcmp(typename, 'uint8')
                rgb = zeros(size(mip_ref,1),size(mip_ref,2),3,'uint8');
                rgb(:,:,1)= uint8(mip_1 * im_intensity);
                rgb(:,:,2)= uint8(mip_1 * im_intensity);
                rgb(:,:,2)= uint8(mip_ref * im_intensity);
                rgb(:,:,3)= uint8(mip_ref * imref_intensity);
                rgb(rgb(:) > 255) = 255 ;
            elseif strcmp(typename, 'uint16')
                rgb = zeros(size(mip_ref,1),size(mip_ref,2),3,'uint16');
                rgb(:,:,1)= uint16(mip_1 * im_intensity);
                rgb(:,:,2)= uint16(mip_1 * im_intensity);
                rgb(:,:,2)= uint16(mip_ref * im_intensity);
                rgb(:,:,3)= uint16(mip_ref * imref_intensity);
                rgb(rgb(:) > 65535) = 65535 ;
            else
                error(['Have not coded for RGB overlay ', ...
                    'with given typename: ', typename, '. Do so here'])
            end
            
            % Make record of the mip drift
            imwrite(rgb, fullfile(mipsRGBDir, rgb_outname))

            % Creat MIPs (maximum intensity projections) of stabilized images
            disp(['creating mips for timepoint=' num2str(time)])
            imSize = size(im) ;
            mip_1 = max(im(:,:,1:round(imSize(end)/2)),[],3);
            mip_2 = max(im(:,:,round(imSize(end)/2):end),[],3);
            mip_11 = squeeze(max(im(1:round(imSize(1)/2),:,:),[],1));
            mip_21 = squeeze(max(im(round(imSize(1)/2):end,:,:),[],1));
            mip_12 = squeeze(max(im(:,1:round(imSize(2)/2),:),[],2));
            mip_22 = squeeze(max(im(:,round(imSize(2)/2):end,:),[],2));
            if strcmp(timechannel, 't')
                m1fn = fullfile(mipStabDir, sprintf(name1,  time)) ;
                m2fn = fullfile(mipStabDir, sprintf(name2,  time)) ;
                m11fn = fullfile(mipStabDir, sprintf(name11,time)) ;
                m21fn = fullfile(mipStabDir, sprintf(name21,time)) ;
                m12fn = fullfile(mipStabDir, sprintf(name12,time)) ;
                m22fn = fullfile(mipStabDir, sprintf(name22,time)) ;
            elseif strcmp(timechannel, 'tc')
                m1fn = fullfile(mipStabDir, sprintf(name1,  time, stabChannel)) ;
                m2fn = fullfile(mipStabDir, sprintf(name2,  time, stabChannel)) ;
                m11fn = fullfile(mipStabDir, sprintf(name11,time, stabChannel)) ;
                m21fn = fullfile(mipStabDir, sprintf(name21,time, stabChannel)) ;
                m12fn = fullfile(mipStabDir, sprintf(name12,time, stabChannel)) ;
                m22fn = fullfile(mipStabDir, sprintf(name22,time, stabChannel)) ;
            end

            % Convert to uint8
            % mip_1 = uint8(mip_1) ;
            % mip_2 = uint8(mip_2) ;
            % mip_11 = uint8(mip_11) ;
            % mip_21 = uint8(mip_21) ;
            % mip_12 = uint8(mip_12) ;
            % mip_2 = uint8(mip_22) ;
            mip_1 = mip_1 * im_intensity ;
            mip_2 = mip_2 * im_intensity ;
            mip_11 = mip_11 * im_intensity ;
            mip_21 = mip_21 * im_intensity ;
            mip_12 = mip_12 * im_intensity ;
            mip_22 = mip_22 * im_intensity ;

            imwrite(mip_1, m1fn,'tiff','Compression','none');
            imwrite(mip_2, m2fn,'tiff','Compression','none');
            imwrite(mip_11,m11fn,'tiff','Compression','none');
            imwrite(mip_21,m21fn,'tiff','Compression','none');
            imwrite(mip_12,m12fn,'tiff','Compression','none');
            imwrite(mip_22,m22fn,'tiff','Compression','none');
        end
    end
end
disp('done building ims and saving stab tifs...')
