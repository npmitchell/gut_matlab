
clear all; clc; close all

addpath_recurse('L:/Streichan/code/gut_matlab/') ;

cd L:\Streichan\data\relMotion_test\Ch1_MIPs
tiff = loadtiff('Ch1_testMIP.tiff') ;

for tt = 1:25 
    rstiff = tiff((0+10*tt):(255+10*tt), 128:383) ;
    saveastiff(rstiff, ['test_Ch1_' sprintf('T%02d', tt) '_MIP.tiff']) ;
end

cd L:\Streichan\data\relMotion_test\Ch2_MIPs
tiff = loadtiff('Ch2_testMIP.tiff') ;

for tt = 1:25 
    rstiff = tiff((10):(265), 128:383) ;
    saveastiff(rstiff, ['test_Ch2_' sprintf('T%02d', tt) '_MIP.tiff']) ;
end
