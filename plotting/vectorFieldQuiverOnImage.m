function [h1, h2] = vectorFieldQuiverOnImage(im, xx, yy, vx, vy, vscale, ...
    options)
%VECTORFIELDQUIVERONIMAGE(im, xx, yy, vx, vy, vscale, options)
%   Plot a vector field (vx,vy) evaluated at grid[xx, yy] on an image im as
%   quiverplot (subsampled quiver)
%
% xx : N x 1 float array
%   x values of PIV grid evaluation points
% yy : M x 1 float array
%   y values of PIV grid evaluation points
% vx : N*M x 1 float array
%   velocity in x direction
% vy : N*M x 1 float array
%   velocity in y direction
% vscale : float
%   magnitude associated with maximum color/intensity in velocity image
% qopts : struct with fields
%   outfn : str
%       path to save image if given
%   label : str
%       colorbar label. Default is '$|v|$ [$\mu$m / min]' 
%   qsubsample : int
%       subsampling factor of the quiver field
%   qscale : float
%       overall scale of the quivers
%   outfn : str
%       output filename for figure as png 
%   colorVectors : bool
%       color the vectors by qopts.colors (default parula colors)
%
%
% Returns
% -------
% h1 : handle for imshow
% h3 : handle for quiverplot
%
% Example Usage 
% -------------
% 
%
%
% NPMitchell 2020

% Default options
labelstr = '$|v|$ [$\mu$m / min]' ;
overlay_quiver = true ;
qsubsample = 10 ;
qscale = 10 ;
transposeVectorShape = false ;
lw = 1.2 ;
colorVectors = false ;

% Unpack options
if isfield(options, 'label')
    labelstr = options.label ;
end
if isfield(options, 'overlay_quiver')
    overlay_quiver = options.overlay_quiver ;
end
if isfield(options, 'qsubsample')
    qsubsample = options.qsubsample ;
end
if isfield(options, 'qscale') 
    qscale = options.qscale ;
end
if isfield(options, 'transposeVectorShape') 
    transposeVectorShape = options.transposeVectorShape ;
end
if isfield(options, 'fig') 
    fig = options.fig ;
end
if isfield(options, 'colorVectors') 
    colorVectors = options.colorVectors ;
end
if isfield(options, 'lw') 
    lw = options.lw ;
end

if colorVectors
    mags = sqrt(vx.^2 + vy.^2) ; 
    
    if isfield(options, 'climv')
        climv = options.climv ;
    else
        maxv = max(mags(:)) ;
        climv = [0,maxv] ;
    end
    if isfield(options, 'colors')
        colors = options.colors ;
    else
        colors = parula(100) ;
    end
    
end


% 
% vangle = reshape(mod(atan2(vy, -vx), 2* pi), gridsz) ;
% speed = reshape(vecnorm([v2dsm_ii(:, 1), v2dsm_ii(:, 2)], 2, 2), gridsz);
ww = length(xx) ;
hh = length(yy) ;
vangle = mod(atan2(vy, -vx), 2* pi) ;
speed = reshape(vecnorm([vx(:), vy(:)], 2, 2), [ww, hh]);

% Compute angle of the velocity vector
if ~all(size(vangle) == [ww, hh])
    vangle = reshape(vangle, [ww, hh]) ;
end

% Set up the figure
if ~ishandle(fig)
    close all
    fig = figure('units', 'normalized', ...
        'outerposition', [0 0 1 1], 'visible', 'off') ;
end
if max(im(:)) <= 1
    h1 = imshow(im) ;
else
    h1 = imagesc(im) ;
    axis equal 
    axis tight
end
hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUIVER 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if transposeVectorShape
    vxR = reshape(vx, [hh, ww]) ;
    vyR = reshape(vy, [hh, ww]) ;
    QX = imresize(vxR, [hh / qsubsample, ww / qsubsample], 'bicubic') ;
    QY = imresize(vyR, [hh / qsubsample, ww / qsubsample], 'bicubic') ;
else
    vxR = reshape(vx, [ww, hh]) ;
    vyR = reshape(vy, [ww, hh]) ;
    QX = imresize(vxR, [ww / qsubsample, hh / qsubsample], 'bicubic') ;
    QY = imresize(vyR, [ww / qsubsample, hh / qsubsample], 'bicubic') ;
end
xq = 1:qsubsample:ww ;
yq = 1:qsubsample:hh ;
[xg, yg] = meshgrid(xx(xq), yy(yq)) ;

if colorVectors
    mags = sqrt(QX.^2 + QY.^2) ; 
    h2 = quiverColorVectors2D(xg(:), yg(:), QX(:)*  qscale, QY(:)*qscale, mags(:), colors, climv, lw) ;
else
    h2 = quiver(xg(:), yg(:), qscale * QX(:), qscale * QY(:), 0, 'k', 'LineWidth', lw) ;   
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phasemap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error('finish this function: make arrow to show scale')

if isfield(options, 'outfn')
    saveas(fig, options.outfn) ;   
    close all
end
