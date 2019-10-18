function h2 = heatmap_on_alphaimage(im, xfield, yfield, cfield, options)
% HEATMAP_ON_IMAGE Plot a scalar field as heatmap on image
% This function has some unresolved bugs 2019-10-06: please fix
%
% Parameters
% ----------
% im : 2d image 
% xfield : x values for heatmap, can be ndgrid, linspace, or scattered (not meshgrid)
% yfield : y values for heatmap, can be ndgrid, linspace, or scattered (not meshgrid)
% cfield : scalar field values evaluated at xfield, yfield
% options : struct optionally containing the following
%   flipy      -- whether to flip the heatmap field in Y
%   alpha      -- opacity value (0, 1) 
%   clims      -- color limits [min max]
%   cmap       -- colormap
%
% Returns
% -------
% fig : figure instance
% ax1 : the image axis
% ax2 : the heatmap axis
% 
% NPMitchell 2019

if isfield(options, 'flipy')
    flipy = options.flipy ;
else
    flipy = true ;
end
if isfield(options, 'clims')
    clims = options.clims ;
else
    clims = [-1 1];
end
if isfield(options, 'cmap')
    cmap = options.cmap ;
else
    cmap = [] ;
end
if isfield(options, 'gridded')
    if options.gridded
        scattered = false ;
    else
        scattered = true ;
    end
    check_for_grid = false ;
else
    check_for_grid = true ;
end


% Get x and y linspace for the image
xx = 1:size(im, 1) ;
yy = 1:size(im, 2) ;
[xg, yg] = ndgrid(xx, yy) ;

% Check if the input data is gridded
if check_for_grid
    dx = gradient(xfield) ;
    dy = gradient(yfield) ;
    if all(dx == dx(1)) && all(dy == dy(1))
        scattered = false ;
    else 
        scattered = true ;
    end
end

% Background image is in imshow then used as alpha
if scattered
    gF = scatteredInterpolant(xfield, yfield, cfield);
    vF = gF(xx, yy) ;
else
    if ~all(size(xfield) == size(cfield))
        % Here we assume since the size is not right, we make a grid
        % This ensures that the image is not transposed. 
        % Note that meshgrid would be transposed.
        [xfield, yfield] = ndgrid(xfield, yfield) ;
    end
    gF = griddedInterpolant(xfield, yfield, cfield);
    vF = gF(xg, yg) ;
end

% PLOT IT
% h1 = imshow(im) ;
h2 = imagesc(xx, yy, vF) ;
set(h2, 'AlphaData', im)

caxis(clims)

if ~isempty(cmap)
    colormap(cmap);
end

if flipy
    set(gca,'ydir','normal');
end


end

