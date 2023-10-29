%% Go through all tiffs in a directory and resave them

% addpath('/mnt/data/code/tubular/utility/bfmatlab/')
addpath('~/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/tubular/utility/bfmatlab/')

scale = 0.25 ;
datdir = '/Volumes/LuxData/MBL/' ;
outdir = '/Volumes/LuxData/MBL_downsampled/' ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

fns = dir(fullfile(datdir, 'Time_*.ome.tif')) ;

for ii = 1:length(fns)
    
    fn = fullfile(fns(ii).folder, fns(ii).name) ;
    
    outname = fullfile(outdir, fns(ii).name) ;
    
    if ~exist(outname, 'file')
        bla = imfinfo(fn);
        filesize = length(bla);

        % read the tiff
        disp(['loading tif: ' fn])
        im = loadtiff(fn) ;

        % resize the tiff
        im = imresize(im, scale) ;

        % Reshape and save file
        disp(['Saving tif: ' outname])
        im = reshape(im,[size(im,1),size(im,2),1,size(im,3),1]); 
        bfsave(im, outname, 'dimensionOrder', 'XYTZC')
    end

end

disp('done')