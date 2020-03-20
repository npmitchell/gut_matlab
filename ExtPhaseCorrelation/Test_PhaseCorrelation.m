%Sample to perform image registration using 
%   ExtPhaseCorrelation.m 
%   PhaseCorrelationPotentialBound.m
%
%                                               eL-ardiansyah
%                                               January, 2010
%                                                       CMIIW
%                               modifications NPMitchell 2020
%                                                        KITP
%============================================================

xshift = 0 ;
yshift = 1 ;

clc
img1 = imread('2.bmp');
img1 = double(rgb2gray(img1)) / double(max(img1(:))) ;
% translate the image by a known amount
img2 = ones(size(img1)) ;
img2((1+yshift):end, (1+xshift):end) = img1(1:(end-yshift), 1:(end-xshift)) ;
montage({img1, img2})
[deltaX, deltaY] = ExtPhaseCorrelation(img2,img1)
[deltaX, deltaY, x, y] = ExtPhaseCorrelationPotentialBound(img2,img1, 'none', 2, 'verbose', true, 'sigmaFilter', 0)

% Add noise to img2
img1 = imread('2.bmp');
img1 = double(rgb2gray(img1)) / double(max(img1(:))) ;
% translate the image by a known amount
img2 = ones(size(img1)) ;
img2((1+yshift):end, (1+xshift):end) = img1(1:(end-yshift), 1:(end-xshift)) ;
q = 0.2 ;
img2 = img2 * (1-q) + q * (rand(size(img2)) - 0.5);
img2(img2 > 1) = 1 ;
img2(img2 < 0) = 0 ;
img2 = imgaussfilt(img2, 1) ;
montage({img1, img2})
[deltaX, deltaY] = ExtPhaseCorrelation(img2,img1)
[deltaX, deltaY, x, y] = ExtPhaseCorrelationPotentialBound(img2,img1, 'none', 2, 'verbose', true, 'sigmaFilter', 3)


% img1 = imread('1.bmp');
% img1 = rgb2gray(img1);
% img2 = imread('3.bmp');
% img2 = rgb2gray(img2);
% [deltaX, deltaY] = ExtPhaseCorrelation(img1,img2)
% img1 = imread('1.bmp');
% img1 = rgb2gray(img1);
% img2 = imread('4.bmp');
% img2 = rgb2gray(img2);
% [deltaX, deltaY] = ExtPhaseCorrelation(img1,img2)