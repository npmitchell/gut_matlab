function [rgb, individuals] = stack2rgb(stack, colors, maxI, pcTiles, outname)
% Assign each channel to a color and adds them together
%
% Parameters
% ----------
% stack : either a filename for a tiff stack or a NxMxQ grayscale stack 
% colors (optional)
%   colors to assign to each channel of the stack
% maxI (optional) 
%   rescaling of the intensities (for each channel)
% outname (optional)
%   filename to output the rgb overlay
%
% Returns
% ------
% rgb : NxMx3 RGB image
%
% Example usage
% -------------
% stack = 'Stack.tif' ;
% colors = random_color(8) ;
% maxI = [] ;
% outname = 'composite.tif'
% stack2rgb(stack, colors, maxI, outname)

% if the stack is specified as a string ('char'), then load it assuming it
% is a tiff!
if ischar(stack)
    stack = readTiff4D(stack) ;
end

nchannels = size(stack, 3) ;
 
% load the colors
if nargin < 2 || isempty(colors)
    colors = define_colors(nchannels) ;
end

% Define the maximum intensity in the LUT
if nargin < 3 || isempty(maxI)
   maxI = ones(nchannels, 1) ./ nchannels ;
end

if nargin < 4
    pcTiles = [1, 99] ;
end

% We want a 2D RGB image. Preallocate Red, Green, and Blue 2d image
rr = zeros(size(stack(:, :, 1))) ;
gg = rr ;
bb = gg ;

% preallocate individuals as a cell array
if nargout > 1
    individuals = cell(nchannels, 1) ;
end

% Consider each channel in the stack
for ii = 1:nchannels
    % normalize the current channel. This is the image we add to RGB,
    % weighted by the current color.
    if strcmpi(class(stack), 'uint16')
        dstack = double(stack(:,:,ii)) ;
        im = (dstack - prctile(dstack(:), pcTiles(1))) ./ ...
            (maxI(ii) * prctile(dstack(:), pcTiles(2))) ;
    else
        error('code for other bit depth here')
    end
    
    % what is the current color, as an RGB array
    ccolor = colors(ii, :) ;
    
    % add in the Red, green, and blue contributions to the
    % composite/overlay based on the current image intensities
    rr = rr + im * ccolor(1) ;
    gg = gg + im * ccolor(2) ;
    bb = bb + im * ccolor(3) ;
    
    if nargout > 1
        individuals{ii} = cat(3, ...
            im * ccolor(1), im * ccolor(2), im * ccolor(3)) ;
    end
    
end

rgb = cat(3, rr, gg, bb) ;

% output the file if desired
if nargin > 4
    imwrite(rgb, outname) ;
end

