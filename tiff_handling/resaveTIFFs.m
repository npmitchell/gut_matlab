%% Go through all tiffs in a directory and resave them

addpath('/mnt/data/code/tubular/utility/bfmatlab/')

datdir = '/mnt/crunch/48YGAL4klarGFPnlsCAAXmCh/' ;
datdir = [datdir '202308101028_180s_1p4um_2mW2mW_48YG4knlsGFPCAAXmCh_0p25_3p0msexposure/'] ;
datdir = [datdir 'data_testReconfiguredToFix/'] ;
outdir = [datdir 'data_resavedViaMATLAB/'] ;
mkdir(outdir)

fns = dir(fullfile(datdir, 'Time_*.ome.tif')) ;

for ii = 1:length(fns)
    
    fn = fullfile(fns(ii).folder, fns(ii).name) ;
    bla = imfinfo(fn);
    filesize = length(bla);
    % image = readtimetable
    % for k = 1:filesize
    %     temp =  imread(name,k ); 
    %     image(:,:,k)  = temp(1:end,1:end);
    % end
    disp(['loading tif: ' fn])
    im = loadtiff(fn) ;

    % Reshape and save file
    outname = fullfile(outdir, fns(ii).name) ;
    disp(['Saving tif: ' outname])
    im = reshape(im,[size(im,1),size(im,2),1,size(im,3),1]); 
    bfsave(im, outname, 'dimensionOrder', 'XYTZC')
end


disp('done')