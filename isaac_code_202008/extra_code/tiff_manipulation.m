%% TIFF MANIPULATION
% Isaac B. Breinyn 2020
%
% This is a script that can seperate a multi-channel 4D tiff into multiple
% TP tifs. This script would be unecessary if the bioformat plugin actually
% worked
%
% Have since changed script to include various modules that are convenient
% for miscellaneous tiff manipulation schemes

clear all; clc; close all

addpath_recurse('L:/Streichan/code/gut_matlab/') ;

%% Trim individual TIFFS in a sequence

% define time and z range for tif
ltp = 59 ;
dir = 'L:\Streichan\data\202007142100_mef2GAL4klarH2AGFPCAAXmCh_e1e2\1006_e2\' ;

for cc = 0:1
parfor tt = 0:ltp
    fn = [dir, '1006_e2_C', num2str(cc), '_', sprintf('T%02d', tt)] ;
    tiff = loadtiff([fn '.tiff']) ;
    tiff = tiff(:,:,5:16) ;
    if cc == 0
        ch = 2 ;
    elseif cc == 1
        ch = cc ;
    end
    saveastiff(tiff, [dir, '1006_e2_Ch', num2str(ch), '_', sprintf('T%02d', tt+1), '.tiff']) ;
end
end

%% Turn a sequence of Tiffs into MIPs

% define time and z range for tif
ltp = 33 ;
dir = 'L:\Streichan\data\202003111630_mef2gal4klarUASCAAXmChHiFP_wo_great\1644_folding_30s_2um_la8_4p5zoom_63x\' ;

parfor cc = 1:2
    for tt = 1:ltp
        fn = [dir, '1644_folding_Ch', num2str(cc), '_', sprintf('T%02d', tt)] ;
        tiff = loadtiff([fn '.tiff']) ;
        
        tiff = max(tiff, [], 3) ;
        tiff = imadjust(tiff) ;
        saveastiff(tiff, [dir, '1644_folding_Ch', num2str(cc), '_', sprintf('T%02d', tt), '_MIP.tiff']) ;
    end
end

%% Turn a sequence of Individual tiffs to a stack

% define time and z range for tif
ltp = 35 ;
dir = 'L:\Streichan\data\202003111830_mef2gal4klarUASCAAXmChHiFP_wo_63x_e4e5\2038_e5\Ch1_MIPs\' ;
fnstack = [dir, '2038_e5_Ch1_MIPs'] ;

for tt = 1:ltp
    fn = [dir, '2038_e5_Ch', num2str(1), '_', sprintf('T%02d', tt), '_MIP'] ;
    tiff = loadtiff([fn '.tiff']) ;
    values = unique(tiff) ;
    thresh = values(round(length(values)/5)) ;
    tiff(tiff<thresh) = 0 ;
    tiff = imadjust(tiff) ;
    tiffstack(:, :, tt) = tiff ;
end

saveastiff(tiffstack, [fnstack '.tiff']) ;

%% Turn a Sequuence of tiffs into a sequence of H5s

dir = 'L:\Streichan\data\202007142100_mef2GAL4klarH2AGFPCAAXmCh_e1e2\1006_e2\Ch2_MIPs\' ;
ltp = 60 ;

for tt = 1:ltp
    fn = [dir, '1006_e2_Ch', num2str(2), '_', sprintf('T%02d', tt), '_maskedMIP_smoothed'] ;
    tiff = loadtiff([fn '.tiff']) ;
%     values = unique(tiff) ;
%     thresh = values(round(length(values)/5)) ;
%     tiff(tiff<thresh) = 0 ;
%     tiff = imadjust(tiff) ;
    h5create([fn '.h5'], '/myData', size(tiff)) ;
    h5write([fn '.h5'], '/myData', tiff) ;
end

%% Turn a sequence of h5s into a single h5

dir = 'L:\Streichan\data\202007142100_mef2GAL4klarH2AGFPCAAXmCh_e1e2' ;
fnstack = [dir, '1006_e2_Ch2_stack'] ;
ltp = 60 ;
fsize = [256 256] ;
h5create([fnstack '.h5'], '/Data', [ltp fsize], 'Datatype', 'uint8') ;

for tt = 1:ltp
    fn = [dir, '\1006_e2\1006_e2_Ch', num2str(2), '_', sprintf('T%02d', tt)] ;
    dset = '/inputData' ;
    h5 = uint8(h5read([fn '.h5'], dset)) ;
    h5stack(tt, :,:) = h5 ;
end

h5write([fnstack '.h5'], '/myData', h5stack) ;

%% Make Colored MIPs
% Make MIPs that are colored based on the z-location of that pixel before
% projection along the z-axis

dir = 'L:\Streichan\data\202007142100_mef2GAL4klarH2AGFPCAAXmCh_e1e2\1006_e2\' ;
ltp = 60 ;

parfor cc = 1:2
    for tt = 1:ltp
        fn = [dir, '1006_e2_Ch', num2str(cc), '_', sprintf('T%02d', tt)] ;
        nfn  = [dir, 'Ch' num2str(cc) '_MIPs\' '1006_e2_Ch', num2str(cc), '_', sprintf('T%02d', tt) ] ;
        tiff = loadtiff([fn '.tiff']) ;
        tiff = max(tiff, [], 3) ;
        saveastiff(tiff, [nfn '_MIP.tiff']) ;
    end
end

