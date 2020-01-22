%% Stabilize images by removing jitter taken from MIPs. 
% SJS and NPM 2019

% A directory 'mipsDir' should already exist and have the non-stab mips in
% there. We compute the shifts to stabilize all volumes and write the 16
% bit stabilized volumes to disk

%%
clc
clear all
disp('defining options...')
addpath('/mnt/data/code/spimCode/');

workDir = [pwd '/']; 

mipsDir = ['./mips/'];
% name of directory to check the stabilization of mips
mips_stab = [mipsDir 'stab_check/'] ;
if ~exist(mips_stab, 'dir')
    mkdir(mips_stab)
end

% Choose reference time for stabilization
t_ref = 28;  % INDEX, not timestamp, of the reference image
alltimes = 27:58 ;
times_todo = [27:58];  % times to overwrite as tifs
num2check = 50 ;
% Choose bit depth as typename
typename = 'uint16' ;
% Give file names for I/O
dummyName = 'Time_%06d_c1.tif';
dummyNameOut = 'Time_%06d_c1_stab.tif';
previewName = 'Time_%06d_c1_stab.png' ;

name_1  = 'view1/mip_1_%03d_c1.tif';
name_11 = 'view2/mip_2_%03d_c1.tif';
name_in  = 'view11/mip_11_%03d_c1.tif';
name_21 = 'view21/mip_21_%03d_c1.tif';
name_22 = 'view22/mip_22_%03d_c1.tif';
im_intensity = 0.01 ;
imref_intensity = 0.005 ;

disp('done defining options')

%% Load MIP data into im_1 and im_2 for all times
disp('Loading MIP data for all times...')
NTimes = length(alltimes);
% preallocate im_1 for speed
tmp = imread([mipsDir, sprintf(name_1,alltimes(1))]) ;
im_1 = zeros([size(tmp) length(alltimes)]) ;
% preallocate im_2 for speed
tmp = imread([mipsDir, sprintf(name_in,alltimes(1))]) ;
im_2 = zeros([size(tmp) length(alltimes)]) ;
for time_id = 1:length(alltimes)
    time = alltimes(time_id) ;
    im_1(:,:,time_id) = imread([mipsDir,sprintf(name_1,time)])+imread([mipsDir,sprintf(name_11,time)]);
    im_2(:,:,time_id) = imread([mipsDir,sprintf(name_in,time)])+imread([mipsDir,sprintf(name_21,time)]);
end
disp('done loading data into im_1 and im_2')
%% Define shifts for each time point ======================================
disp('Defining shifts...')
shifts = struct('x_1',[],'y_1',[],'x_2',[],'y_2',[]);
for time = 1 :NTimes
    [shiftx,shifty,~] = xcorr2fft(im_1(:,:,time), im_1(:,:,t_ref));
    shifts(time).x_1 = shiftx;
    shifts(time).y_1 = shifty; 
    
    [shiftx,shifty,~] = xcorr2fft(im_2(:,:,time), im_2(:,:,t_ref));
    shifts(time).x_2 = shiftx;
    shifts(time).y_2 = shifty; 
end
disp('done defining shifts')

%% Plot the shifts
disp('Plotting shifts...')
close('all')
x_1 =  cat(1,shifts.x_1); % rows in the 3D image; 
plot(alltimes, x_1, '.-', 'DisplayName', 'x1')
hold all;
y_1 =  cat(1,shifts.y_1); %columns in the 3D image
plot(alltimes, y_1, '.-', 'DisplayName', 'y1')
x_2 =  cat(1,shifts.x_2); %is the columns 
plot(alltimes, x_2, '.-', 'DisplayName', 'x2')
y_2 =  cat(1,shifts.y_2); %must then be the leaves ; 
plot(alltimes, y_2, '.-', 'DisplayName', 'y2')
ylabel('shift [pixels]')
xlabel('timestamp')
legend('location', 'best')
title('Jitter stabilization')
saveas(gcf, fullfile(mipsDir, 'jitter_stabilization.png'))
disp('done plotting shifts, see Figure')
disp('Saving shifts to shifts_stab.mat in mipsDir')
save(fullfile(mipsDir, 'shifts_stab.mat'), 'shifts')

%% Make reference MIP
disp('Creating reference MIP...')
close('all')
% Define stackSize, which is the number of frames in the z dimension
disp('defining stackSize...')
done = false ;
stackSize = 0 ;
name_ref = sprintf(dummyName,alltimes(t_ref));
while ~done
    stackSize = stackSize + 1 ;
    try
        tmp = imread(name_ref, stackSize);
    catch
        done = true ;
    end
end    
stackSize = stackSize - 1 ;

disp('building reference MIP...')
name_ref = sprintf(dummyName,alltimes(t_ref));
% preallocate im_ref3D for speed
im_ref3D = zeros([size(tmp) stackSize]) ;
for z = 1 : stackSize
    im_ref3D(:,:,z)= imread(name_ref,z);
end
mip_ref = squeeze(max(im_ref3D,[],2));
clear tmp
disp('done creating reference MIP')

%% Build image for each timepoint =========================================
disp('Running through timepoints to build ims...')
for time_id = 1 : length(alltimes)
    time = alltimes(time_id);
    if ismember(time, times_todo)
        name_in = sprintf(dummyName, time);

        % Check that we're not appending to an existing file
        name_out = sprintf(dummyNameOut,time) ;
        if exist(name_out, 'file')
            error(['Output file already exists: ' name_out ])
        end

        % preallocate im_2_3D for speed, and 
        % Specify typename for correct output bitdepth
        tmp = imread(name_in,z) ;
        im_2_3D = zeros([size(tmp) stackSize], typename) ;

        clear tmp
        for z = 1 : stackSize
            im_2_3D(:,:,z)= imread(name_in,z);
        end

        % Padding?
        im = im_2_3D;
        if x_1(time_id)~=0
            if x_1(time_id)>0
                im(x_1(time_id):end,:,:) = im(1:end-x_1(time_id)+1,:,:);
            else
                im(1:end+x_1(time_id)+1,:,:) = im(-x_1(time_id):end,:,:);
            end
        end

        % Padding in Y?
        %y_1(time)= 10;
        if y_1(time_id)~=0
            if y_1(time_id)>0
                im(:,y_1(time_id):end,:) = im(:,1:end-y_1(time_id)+1,:);
            else
                im(:,1:end+y_1(time_id)+1,:) = im(:,-y_1(time_id):end,:);
            end
        end

        if y_2(time_id)~=0
            if y_2(time_id)>0
                im(:,:,+y_2(time_id):end) = im(:,:,1:end-y_2(time_id)+1);
            else
                im(:,:,1:end+y_2(time_id)+1) = im(:,:,-y_2(time_id):end);
            end
        end

        % Write to disk. Note that if uint16 specified upon instantiation, then
        % imwrite will respect the class type.
        for z = 1 : stackSize
            imwrite(im(:,:,z), name_out,'tiff',...
                'Compression','none','WriteMode','append') ;
        end
        disp(['saved image ' sprintf(dummyNameOut,time)])

        % Make a color of the current mip wrt the reference
        mip_1   = squeeze(max(im,[],2));
        rgb = zeros(size(mip_ref,1),size(mip_ref,2),3,'uint8');
        rgb(:,:,1)= uint8(mip_1 * im_intensity);
        rgb(:,:,2)= uint8(mip_ref * imref_intensity);
        %res(time).rgb = rgb;

        % Make record of the mip drift
        % imshow(rgb,[]);
        % drawnow
        imwrite(rgb, fullfile(mips_stab, sprintf(previewName,time)))
    end
end
disp('done building ims and saving stab tifs...')

%% check the result; 
disp('Checking the first and last frame...')
mip_1   = squeeze(max(im,[],2));
mip_ref = squeeze(max(im_ref3D,[],2));
rgb = zeros(size(mip_ref,1),size(mip_ref,2),3,'uint8');
rgb(:,:,1)= uint8(mip_1/100);
rgb(:,:,2)= uint8(mip_ref/100);
imshow(rgb,[])
waitfor(fig)
disp('done checking, see Figure')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check on the first images; 
disp('Check on each time point for first mips...')
for time = 1:num2check
    im = im_1(:,:,time);
    rgb = zeros(size(im,1),size(im,2),3,'uint8');

    im_ref = im_1(:,:,t_ref);
    if x_1(time)~=0
        if x_1(time)>0
            im(x_1(time):end,:) = im(1:end-x_1(time)+1,:);
        else
            im(1:end+x_1(time)+1,:) = im(-x_1(time):end,:);
        end
    end

    %y_1(time)= 10;
    if y_1(time)~=0
        if y_1(time)>0
            im(:,y_1(time):end) = im(:,1:end-y_1(time)+1);
        else
            im(:,1:end+y_1(time)+1) = im(:,-y_1(time):end);
        end
    end
    %imshow(double(im)+1*double(im_ref),[])
    rgb(:,:,1)= uint8(im * im_intensity);
    rgb(:,:,2)= uint8(im_ref * imref_intensity);
    imshow(rgb,[])
    mm(time)= getframe;
end
disp('Done checking each timepoint MIP #1')

%% check on the second images; 
disp('Check on each time point for second mips...')
for time = 1:num2check
    im = im_2(:,:,time);
    rgb = zeros(size(im,2),size(im,1),3,'uint8');

    im_ref = im_2(:,:,t_ref);
    if x_2(time)~=0
        if x_2(time)>0
            im(x_2(time):end,:) = im(1:end-x_2(time)+1,:);
        else
            im(1:end+x_2(time)+1,:) = im(-x_2(time):end,:);
        end
    end

    %y_1(time)= 10;
    if y_2(time)~=0
        if y_2(time)>0
            im(:,y_2(time):end) = im(:,1:end-y_2(time)+1);
        else
            im(:,1:end+y_2(time)+1) = im(:,-y_2(time):end);
        end
    end
    %imshow(double(im)+1*double(im_ref),[])
    rgb(:,:,1)= uint8(im' * im_intensity);
    rgb(:,:,2)= uint8(im_ref' * imref_intensity);
    imshow(rgb,[])
    mm2(time)= getframe;
end
disp('Done checking each timepoint MIP #2')