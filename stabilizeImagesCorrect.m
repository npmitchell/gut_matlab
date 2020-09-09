%% Stabilize images by removing jitter taken from MIPs. 
% NPM 2019-2020
%
% Script version of stabilizeImages() function.
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
% MIPs of the data specified by dataFileName must exist in:
%    mipsDir/view1/mip_1_%03d_c1.tif
%    mipsDir/view2/mip_2_%03d_c1.tif
%    mipsDir/view11/mip_11_%03d_c1.tif
%    mipsDir/view21/mip_21_%03d_c1.tif
%    mipsDir/view12/mip_12_%03d_c1.tif
%    mipsDir/view22/mip_22_%03d_c1.tif
% This is accomplished by running make_mips.m
%
% Parameters
% ----------
% mipsDir : where MIPs are stored
% mips_stab_check : output dir for the stabilized RGB overlays 
% mipoutdir : output dir for the stabilized MIPs
% im_intensity : scale for MIP output intensity
% imref_intensity : scale for reference MIP intensity in RGB overlay
% t_ref : timestamp of reference, matches other timepoints to this one
% alltimes : all timestamps to consider in the sequence
% times_todo : timestamps to correct/make output
% typename : output volume datatype, as string ('uint8', 'uint16')
% dataFileName : input filename format string for data to stabilize
% dataFileNameOut : output filename format string for stabilized volumes
% previewName : RGB overlay filename format string
%
% Returns
% -------
% - jitter/drift-corrected 3D volumes 
% - MIPs of the jitter/drift-corrected volumes within mipoutdir
% - colored overlays of one view with the reference time. 

%%
clc ; clearvars
disp('defining options...')
overwrite = false ;
% addpath('/mnt/data/code/spimCode/');

%% Define directory structure
mipsDir = fullfile(cd, ['mips' filesep]) ;
% name of directory to check the stabilization of mips
mips_stab_check = [mipsDir 'stab_check' filesep] ;
mipoutdir = [mipsDir 'mips_stab' filesep] ;

%% Options for scaling the image intensity
im_intensity = 1 ; % 0.01 ;
imref_intensity = 1 ; % 0.005 ;
% Choose reference time for stabilization
t_ref = 11;  % timestamp (not index) of the reference image
alltimes = [0:169] ;
times_todo = [0:169];  % times to overwrite as tifs
% Choose bit depth as typename
typename = 'uint16' ;
% Give file names for I/O
dataFileName = 'Time_%06d_c1.tif';
dataFileNameOut = 'Time_%06d_c1_stab.tif';
previewName = 'Time_%06d_c1_stab.png' ;

disp('done defining options')

%% Make the subdirectories for the mips if not already existing
mipdirs = {mipoutdir, mips_stab_check, ...
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
name1  = fullfile('view1', 'mip_1_%03d_c1.tif');
name2  = fullfile('view2', 'mip_2_%03d_c1.tif');
name11 = fullfile('view11', 'mip_11_%03d_c1.tif');
name21 = fullfile('view21', 'mip_21_%03d_c1.tif');
name12 = fullfile('view12', 'mip_12_%03d_c1.tif');
name22 = fullfile('view22', 'mip_22_%03d_c1.tif');

%% Load MIP data into im_1 and im_2 for all times
disp('Loading MIP data for all times...')
NTimes = length(alltimes);
% preallocate im_1 for speed
tmp = imread([mipsDir, sprintf(name1,alltimes(1))]) ;
im_1 = zeros([size(tmp) length(alltimes)]) ;
% preallocate im_2 for speed
tmp = imread([mipsDir, sprintf(name11,alltimes(1))]) ;
im_2 = zeros([size(tmp) length(alltimes)]) ;
for tid = 1:length(alltimes)
    time = alltimes(tid) ;
    % weight the different views equally
    im_1(:,:,tid) = imread([mipsDir,sprintf(name1,time)])+imread([mipsDir,sprintf(name2,time)]);
    im_2(:,:,tid) = imread([mipsDir,sprintf(name11,time)])+imread([mipsDir,sprintf(name21,time)]);
    im_3(:,:,tid) = imread([mipsDir,sprintf(name12,time)])+imread([mipsDir,sprintf(name22,time)]);
end
t_ref_ind = find( alltimes == t_ref ) ;

disp('done loading data into im_1 and im_2')
%% Define shifts for each time point ======================================
disp('Defining shifts...')
shifts = struct('x_1',[],'y_1',[],'x_2',[],'y_2',[]);
for time = 1 :NTimes
    [shiftx,shifty,~] = xcorr2fft(im_1(:,:,time), im_1(:,:,t_ref_ind));
    shifts(time).x_1 = shiftx;
    shifts(time).y_1 = shifty; 
    
    [shiftx,shifty,~] = xcorr2fft(im_2(:,:,time), im_2(:,:,t_ref_ind));
    shifts(time).x_2 = shiftx;
    shifts(time).y_2 = shifty; 
    
    [shiftx,shifty,~] = xcorr2fft(im_3(:,:,time), im_3(:,:,t_ref_ind));
    shifts(time).x_3 = shiftx;
    shifts(time).y_3 = shifty; 
end
disp('done defining phase correlations')

%% Convert correlations to shifts
x_1 =  cat(1,shifts.x_1); % rows        (x) 
y_1 =  cat(1,shifts.y_1); % columns     (y)
x_2 =  cat(1,shifts.x_2); % columns #2  (y)
y_2 =  cat(1,shifts.y_2); % leaves      (z) 
x_3 =  cat(1,shifts.x_3); % rows    #2  (x) 
y_3 =  cat(1,shifts.y_3); % leaves  #2  (z) 
% Average contributions from different views on the same shift axis
dx = round(0.5 * (x_1 + x_3)) ;
dy = round(0.5 * (y_1 + x_2)) ;
dz = round(0.5 * (y_2 + y_3)) ;

%% Plot the shifts
disp('Plotting shifts...')
close('all')
hold all;
% Consider view 0
plot(alltimes, x_1, '.-', 'DisplayName', 'x (dim1 of view 0)')
plot(alltimes, y_1, '.-', 'DisplayName', 'y (dim2 of view 0)')
% Also consider view 1,2
plot(alltimes, y_2, '.-', 'DisplayName', 'z (dim2 of view 1)')
plot(alltimes, x_3, 's--', 'DisplayName', 'x (dim1 of view 2)')
plot(alltimes, x_2, 'o--', 'DisplayName', 'y (dim1 of view 1)')
plot(alltimes, y_3, '^--', 'DisplayName', 'z (dim2 of view 2)')
ylabel('shift [pixels]')
xlabel('timestamp')
legend('location', 'best')
title('Jitter stabilization')
saveas(gcf, fullfile(mipsDir, 'jitter_stabilization.png'))
disp('done plotting shifts, see Figure')
disp('Saving shifts to shifts_stab.mat in mipsDir')
save(fullfile(mipsDir, 'shifts_stab.mat'), 'shifts')

%% Clear the individual shift values
clearvars x_1 y_1 x_2 y_2 x_3 y_3

%% Make reference MIP
disp('Creating reference MIP...')
close('all')
% Define stackSize, which is the number of frames in the z dimension
disp('defining stackSize...')
done = false ;
stackSize = 0 ;
name_ref = sprintf(dataFileName, t_ref);
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
name_ref = sprintf(dataFileName, t_ref);
% preallocate im_ref3D for speed
im_ref3D = zeros([size(tmp) stackSize]) ;
for z = 1 : stackSize
    im_ref3D(:,:,z)= imread(name_ref,z);
end
mip_ref = squeeze(max(im_ref3D,[],3));
clear tmp
disp('done creating reference MIP')

%% Build image for each timepoint =========================================
disp('Running through timepoints to build ims...')
for tid = 1 : length(alltimes)
    time = alltimes(tid);
    if ismember(time, times_todo)
        % The original image in 3d
        im0fn = sprintf(dataFileName, time);

        % Check that we're not appending to an existing file
        name_out = sprintf(dataFileNameOut,time) ;
        if exist(name_out, 'file')
            if overwrite
                disp(['Output file already exists: ' name_out ])
            else
                error(['Output file already exists: ' name_out ])
            end
            overwrite_dat = true ;
        else
            overwrite_dat = false ;
        end

        % preallocate im_2_3D for speed, and 
        % Specify typename for correct output bitdepth
        tmp = imread(im0fn,z) ;
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
        for z = 1 : stackSize
            if overwrite_dat
                imwrite(im(:,:,z), name_out,'tiff',...
                    'Compression','none') ;
                overwrite_dat = false ;
            else
                imwrite(im(:,:,z), name_out,'tiff',...
                    'Compression','none','WriteMode','append') ;
            end
        end
        disp(['saved image ' sprintf(dataFileNameOut,time)])

        % Make a color of the current mip wrt the reference
        mip_1   = squeeze(max(im,[],3));
        rgb = zeros(size(mip_ref,1),size(mip_ref,2),3,'uint8');
        rgb(:,:,1)= uint8(mip_1 * im_intensity);
        rgb(:,:,2)= 0.5 * uint8(mip_1 * im_intensity);
        rgb(:,:,2)= 0.5 * uint8(mip_ref * im_intensity);
        rgb(:,:,3)= uint8(mip_ref * imref_intensity);
        %res(time).rgb = rgb;

        % Make record of the mip drift
        % imshow(rgb,[]);
        % drawnow
        imwrite(rgb, fullfile(mips_stab_check, sprintf(previewName,time)))
       
        % Creat MIPs (maximum intensity projections) of stabilized images
        disp(['creating mips for timepoint=' num2str(time)])
        imSize = size(im) ;
        mip_1 = max(im(:,:,1:round(imSize(end)/2)),[],3);
        mip_2 = max(im(:,:,round(imSize(end)/2):end),[],3);
        mip_11 = squeeze(max(im(1:round(imSize(1)/2),:,:),[],1));
        mip_21 = squeeze(max(im(round(imSize(1)/2):end,:,:),[],1));
        mip_12 = squeeze(max(im(:,1:round(imSize(2)/2),:),[],2));
        mip_22 = squeeze(max(im(:,round(imSize(2)/2):end,:),[],2));
        m1fn = fullfile(mipoutdir, sprintf(name1,  time)) ;
        m2fn = fullfile(mipoutdir, sprintf(name2,  time)) ;
        m11fn = fullfile(mipoutdir, sprintf(name11,time)) ;
        m21fn = fullfile(mipoutdir, sprintf(name21,time)) ;
        m12fn = fullfile(mipoutdir, sprintf(name12,time)) ;
        m22fn = fullfile(mipoutdir, sprintf(name22,time)) ;
        
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
disp('done building ims and saving stab tifs...')

%% check the result; 
disp('Checking the first and reference frame...')
mip_1   = squeeze(max(im,[],2));
mip_ref = squeeze(max(im_ref3D,[],2));
rgb = zeros(size(mip_ref,1),size(mip_ref,2),3,'uint8');
rgb(:,:,1)= uint8(mip_1 * im_intensity);
rgb(:,:,2)= uint8(mip_ref * imref_intensity);
imshow(rgb,[])
set(gcf, 'visible', 'on')

