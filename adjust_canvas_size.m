%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Go through many images that are different sizes in pwd, plot them on a 
% canvas that is as large as the biggest one (in each dim) so that one can 
% make a single tiff of all the images.
% SJS and NPM 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Options
option = 2; % 1 == select timepoints, 2 == adjust all name matches in dir
save_stack = false ; % true for stack output, false for individual images
% outdir = '../PullbackImages_010step_extended_shifted_adjusted/' ;
% outdir = '../PullbackImages_010step_relaxed_extended_adjusted/' ;
outdir = fullfile(cd, 'adjusted') ;
bitdepth = 'uint16' ;

%% Naming structure
% For option 1
chan  = 1;
nameFormat = ['cmp_',num2str(chan),'_1_T%04d.tif']; 

% For option 2
nameFormat = 'Time_*_c1_stab_corr.png' ;
nameFormat = 'Time_*_c1_stab.png' ;

%% Create outdir
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

%% Option 1: select timepoints
if option == 1
    
    % Select timepoints
    tInit = 110 ; 
    tEnd  = 242 ;
    step  = 1;
    % determine the maximum with and height in the image series; 
    w = zeros(1,tEnd-tInit+1);
    h = zeros(1,tEnd-tInit+1);
    s = struct('info',[]);
    tInd = tInit : step: tEnd;
    for time = 1:length(tInd)
        name = sprintf(nameFormat,tInd(time));
        %name = ['Cyl2_p_',num2str(time),'.tif'];
        info = imfinfo(name);
        s(time).info = info; 

        w(time) = s(time).info.Width;
        h(time) = s(time).info.Height;
    end

    % Now take the max of both width and height
    mw = max(w);
    mh = max(h);

    % adjust the canvas size by initializing large empty image
    im = zeros(mh,mw,tEnd-tInit+1,'uint16');
    for time = 1:length(tInd)
        name = sprintf(nameFormat,tInd(time));
        %name = ['Cyl2_p_',num2str(time),'.tif'];
        temp = imread(name);

        st = size(temp);
        stD = floor(([mh,mw]-st)/2)+1;

        im( stD(1):(stD(1)+st(1)-1) , stD(2):(stD(2)+st(2)-1) ,time) = temp;
    end
end

%%  Option 2: all matching files
if option == 2
    fns = dir(nameFormat) ;
    % determine the maximum with and height in the image series; 
    w = zeros(1, length(fns));
    h = zeros(1, length(fns));
    s = struct('info',[]);
    for i = 1:length(fns)
        name = fns(i).name ;
        %name = ['Cyl2_p_',num2str(time),'.tif'];
        info = imfinfo(name);
        s(i).info = info; 
        w(i) = s(i).info.Width;
        h(i) = s(i).info.Height;
    end
    
    % Now take the max of both width and height
    mw = max(w);
    mh = max(h);

    % adjust the canvas size by initializing large empty image
    im = zeros(mh, mw, length(fns), bitdepth);
    for i = 1:length(fns)
        name = fns(i).name ;
        %name = ['Cyl2_p_',num2str(time),'.tif'];
        temp = imread(name);
        % convert to 2d
        if length(size(temp)) > 2
            temp = rgb2gray(temp);
        end

        st = size(temp);
        stD = floor(([mh, mw] - [st(1), st(2)]) * 0.5) + 1;

        im( stD(1):(stD(1)+st(1)-1), stD(2):(stD(2)+st(2)-1), i) = temp;
    end
end

%% Save the adjusted images
if save_stack
    % Write the result as a tiff stack
    for k = 1 : size(im,3)    
        imwrite(im(:,:,k),fullfile(outdir, ['adjust',num2str(chan),'.tif']), ...
            'tiff','Compression','none','WriteMode','append');
    end
else
    % Write the result as a tiff (single frame)
    for k = 1 : size(im,3)    
        imwrite(im(:,:,k),fullfile(outdir, [fns(k).name(1:end - 4),'_adjust.tif']), ...
            'tiff', 'Compression', 'none');
    end    
end
disp('done')
