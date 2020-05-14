%% Walking average of frames for increasing time resolution

% 1/4  1/2  1/4
%      1/2  1/2
%           1/4  1/2  1/4
%                1/2  1/2
%                     1/4  1/2  1/4
%                          1/2  1/2
%                               1/4  1/2  1/4

fns = './Time_*.png' ;
outdir = './walking_average' ;
outname = 'walking_interp_%06d.png' ;

% build output prep
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
outFileName = fullfile(outdir, outname) ;
fns = dir(fns) ;
dmy = 0 ;
for ii = 1:length(fns)
    disp(num2str(ii)) 
    im = imread(fullfile(fns(ii).folder, fns(ii).name)) ;
    imwrite(im, sprintf(outFileName, dmy)) ;
    dmy = dmy + 1 ;
    
    % 50/50 split with next frame
    im = 0.5 * im + 0.5 * imread(fullfile(fns(ii+1).folder, fns(ii+1).name)) ;
    imwrite(im, sprintf(outFileName, dmy)) ;
    dmy = dmy + 1 ;
end
disp('done')