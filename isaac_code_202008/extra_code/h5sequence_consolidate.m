addpath_recurse('C:/Streichan/code/gut_matlab/') ;

cd('C:\Streichan\data\confocal_data\mef2GAL4klarUASCAAXmChHiFP\202003111830_mef2gal4klarUASCAAXmChHiFP_wo_63x_e4e5\midfold_e4\overlays\masked_ch2')


fn = 'midfold_e4_T%02d_Ch2_overlayMIP_masked_Probabilities.h5' ;

array = zeros(2, 512,512,33) ;

for tt = 1:33
    slice = h5read(fullfile(cd, sprintf(fn,tt)), '/exported_data') ;
    array(:,:,:,tt) = slice ;
end

h5create( fullfile(cd, 'midfold_e4_Ch2_overlayMIP_masked_Probabilities.h5'), '/data', [2 512 512 33]) ;
h5write( fullfile(cd, 'midfold_e4_Ch2_overlayMIP_masked_Probabilities.h5'), '/data', array) ;
    