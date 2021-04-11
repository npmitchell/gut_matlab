% make mips of data (for example, to see saggital plane of data over time)
% using tiffMIP

addpath('/mnt/data/code/gut_matlab/tiff_handling')
slices = [] ;
outDir = './mipPreviews/' ;
for ch = [1,2]
    for tt = 0:10:110
        fileName = sprintf('./TP%d_Ch%d_Ill0_Ang0,60,120,180,240,300.tif', ...
            tt, ch) ;
        
        disp([ 'reading ' fileName]) 
        
        mipfns = {fullfile(outDir, sprintf('mip_t%04d_c%d_slice1.png', tt, ch)),...
            fullfile(outDir, sprintf('mip_t%04d_c%d_slice2.png', tt, ch)),...
            fullfile(outDir, sprintf('mip_t%04d_c%d_slice3.png', tt, ch))};
               
        mip = tiffMIP(fileName, mipfns, 'middle', 0.5, [1,2,3]) ;
    end
end