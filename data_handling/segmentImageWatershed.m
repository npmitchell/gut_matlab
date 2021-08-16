function [bgm, DL] = segmentImageWatershed(img, adaphisteqClip, strelRadius)
% [bgm, DL] = segmentImageWatershed(img, adaphisteqClip, strelRadius)
% Uses Lin's algorithm to segment cells/polygons from image img. 
%
% Parameters
% ----------
% img : 2d image
% adaphisteqClip : float (default=0.2 if empty)
%   threshold for clipping
% strelRadius : int (default=3 if empty)
%   radius for structured element in local minimum filter
%
% Returns
% -------
% bgm : 2d NxM logical array
%   segmentation of img from watershed transform of morphological
%   reconstruction of img
% DL : 2d NxM float array
%   watershed transform of morphological reconstruction of img
%
% Yuzeng Lin & NPMitchell 2021

if nargin < 2
    adaphisteqClip = 0.2 ;
elseif isempty(adaphisteqClip)
    adaphisteqClip = 0.2 ;
end
if nargin < 3
    strelRadius = 3 ;
elseif isempty(adaphisteqClip)
    strelRadius = 3 ;
end

% First adaptively histogram the image
% roi = drawrectangle();
% rectangle_position = round(roi.Position);
% cropped1 = imcrop(img, rectangle_position);
%histeq(img)
clipped = adapthisteq(img,'clipLimit',adaphisteqClip,'Distribution','rayleigh');

% Now perform morphological reconstruction
% see https://www.mathworks.com/help/images/understanding-morphological-reconstruction.html
Iflip = max(clipped(:))-clipped;
se = strel('disk', strelRadius);
Ie = imerode(Iflip,se); % eroded image by strel disk
Iobr = imreconstruct(Ie,Iflip);

% Perform second morphological reconstruction 
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(Iobrd, Iobr);
bw = imbinarize(Iobrcbr);

% Now watershed the binary image
DL = watershed(bwdist(bw));
bgm = DL == 0;
