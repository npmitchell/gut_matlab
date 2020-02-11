function [h1, h2, h3] = scalarVectorFieldsOnImage(im, xxs, yys, sf, ...
    xxv, yyv, vx, vy, options)
%SCALARVECTORFIELDSONIMAGE(im, xx, yy, vx, vy, options)
%   Plot both a scalar field (sf) and vector field (vx,vy) on an image im
%
% xx : N x 1 float array
%   x values of PIV grid evaluation points
% yy : M x 1 float array
%   y values of PIV grid evaluation points
% vx : N*M x 1 float array
%   velocity in x direction
% vy : N*M x 1 float array
%   velocity in y direction
% qopts : struct with fields
%   style : str ('diverging' 'phase')
%       style of overlaid scalar field, default is diverging
%   sscale : float (optional, used only if style == 'phase')
%       scalar field maximum for opacity/alpha
%   alpha : float (optional, used if style ~= 'phase')
%       opacity of overlaid scalar field
%   outfn : str
%       path to save image if given
%   label : str
%       colorbar label. Default is '$|v|$ [$\mu$m / min]' 
%   qsubsample : int
%       subsampling factor of the quiver field
%   overlay_quiver : bool (default=true)
%       whether to show the quiverplot overlay
%   qscale : float
%       overall multiplication of the quiver magnitude for overlay
%       (then plotted at true value, "scale=0" in quiver())
%   outfn : str
%       output filename for figure as png 
%
% Returns
% -------
% h1 : handle for imshow
% h2 : handle for imagesc
% h3 : handle for quiverplot
%
% NPMitchell 2020

% Default options
labelstr = '' ;
overlay_quiver = true ;
qsubsample = 10 ;
qscale = 10 ;
sscale = 0 ;
alphaVal = 0.8 ;
style = 'diverging' ;

% Unpack options
if isfield(options, 'style') 
    style = options.style ;
end
if isfield(options, 'sscale') 
    sscale = options.sscale ;
end
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
if isfield(options, 'alpha') 
    alphaVal = options.alpha ;
end

% 
% vangle = reshape(mod(atan2(vy, -vx), 2* pi), gridsz) ;
% speed = reshape(vecnorm([v2dsm_ii(:, 1), v2dsm_ii(:, 2)], 2, 2), gridsz);
wws = length(xxs) ;
hhs = length(yys) ;

% Set up the figure
close all
fig = figure('units', 'normalized', ...
    'outerposition', [0 0 1 1], 'visible', 'off') ;
h1 = imshow(im) ;
hold on;
h2 = imagesc(xxs, yys, sf) ;
if strcmp(style, 'phase')
    set(h2, 'AlphaData', sf / sscale)
elseif strcmp(style, 'diverging')
    alpha(alphaVal) ;
    if sscale > 0
        caxis([-sscale, sscale])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUIVER 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if overlay_quiver    
    wwv = length(xxv) ;
    hhv = length(yyv) ;
    vx = reshape(vx, [wwv, hhv]) ;
    vy = reshape(vy, [wwv, hhv]) ;
    QX = imresize(vx, [wwv / qsubsample, hhv / qsubsample], 'bicubic') ;
    QY = imresize(vy, [wwv / qsubsample, hhv / qsubsample], 'bicubic') ;
    xq = 1:qsubsample:wwv ;
    yq = 1:qsubsample:hhv ;
    [xg, yg] = meshgrid(xxv(xq), yyv(yq)) ;

    h3 = quiver(xg(:), yg(:), qscale * QX(:), qscale * QY(:), 0, 'k', 'LineWidth', 1.2) ;
else
    h3 = [] ;
end

% Add title (optional)
if isfield(options, 'title')
    title(options.title, 'Interpreter', 'latex')
end

% Add the colorbar in the style set in options struct
if strcmp(style, 'phase')
    %%%%%%%%%%%%%%%%%%%
    % Phasemap
    %%%%%%%%%%%%%%%%%%%
    colormap phasemap
    caxis([0, 2*pi])
    ylim([size(im, 2) * 0.25, size(im, 2) * 0.75])
    set(gca, 'Position', [0 0.11 0.85 0.8]) ;
    % Add phasebar
    phasebar('location', [0.87, 0.7, 0.1, 0.1]) ;
    % Add colorbar
    cax = axes('Position',[.9 .3 .02 .3]) ;
    [~, yyq] = meshgrid(0:4, 0:100) ;
    imshow(fliplr(yyq/max(yyq(:))))
    axis on
    yyaxis right
    ylabel(labelstr, 'color', 'k', ...
        'Interpreter', 'Latex')
    yticks([0 1])
    yticklabels({'0', num2str(sscale)})
    xticks([])
    yyaxis left
    yticks([])
    cax.YAxis(1).Color = 'k';
    cax.YAxis(2).Color = 'k';
elseif strcmp(style, 'diverging')
    colormap bwr
    if isfield(options, 'ylim')
        ylim(options.ylim)
    end
    set(gca, 'Position', [0 0.11 0.85 0.8]) ;
    % Add colorbar
    c = colorbar('Position',[.9 .333 .02 .333]) ;
    % ylabel(cax, labelstr, 'color', 'k', ...
    %     'Interpreter', 'Latex')
    
    % Make colorbar share the alpha of the image
    % Manually flush the event queue and force MATLAB to render the colorbar
    % necessary on some versions
    drawnow
    % Get the color data of the object that correponds to the colorbar
    cdata = c.Face.Texture.CData;
    % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
    cdata(end,:) = uint8(alphaVal * cdata(end,:));
    % Ensure that the display respects the alpha channel
    c.Face.Texture.ColorType = 'truecoloralpha';
    % Update the color data with the new transparency information
    c.Face.Texture.CData = cdata;
    c.Label.String = labelstr ;
    c.Label.Interpreter = 'latex' ;

else
    error('have not coded for this style yet')
end

% Save the image if outfn is supplied
if isfield(options, 'outfn')
    disp(['scalarVectorFieldsOnImage: saving ' options.outfn])
    saveas(fig, options.outfn) ;   
    close all
end
