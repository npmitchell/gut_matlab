function [h1, h2, h3] = scalarVectorFieldsOnSurface(faces, vertices, sf, ...
    xxv, yyv, zzv, vx, vy, vz, options)
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
interpreter = 'latex' ;
qscale = 10 ;
sscale = 0 ;
alphaVal = 0.8 ;
style = 'diverging' ;
qsubsample = 1 ;
xlabelstr = '' ;
ylabelstr = '' ;
zlabelstr = '' ;
titlestr = '' ;

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
if isfield(options, 'qscale') 
    qscale = options.qscale ;
end
if isfield(options, 'qsubsample')
    qsubsample = options.qsubsample ;
end
if isfield(options, 'alpha') 
    alphaVal = options.alpha ;
end


% Set up the figure
close all
fig = figure('units', 'normalized', ...
    'outerposition', [0 0 1 1], 'visible', 'off') ;

% Add the scalar field defined on faces
h1 = patch( 'Faces', faces, 'Vertices', vertices, 'FaceVertexCData', sf, ...
    'FaceColor', 'flat', 'EdgeColor', 'none', ...
    'SpecularStrength', 0.1, 'DiffuseStrength', 0.1, ...
    'AmbientStrength', 0.8 );
hold on;
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUIVER 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if qsubsample > 1
    error('have not implemented here')
end
h3 = quiver3(xxv(:), yyv(:), zzv(:), ...
    qscale * vx(:), qscale * vy(:), qscale * vz(:), 0, 'k', 'LineWidth', 1.2) ;

% View and axis limits
if isfield(options, 'view') 
    view(options.view(1), options.view(2)) ;
end
if isfield(options, 'xlim') 
    xlim(options.xlim) ;  
end
if isfield(options, 'ylim') 
    ylim(options.ylim) ;
end
if isfield(options, 'zlim')
    zlim(options.zlim) ;
end

% Labels and plot title
if ~isempty(xlabelstr)
    xlabel(xlabel, 'Interpreter', interpreter)
end
if ~isempty(ylabelstr)
    ylabel(ylabel, 'Interpreter', interpreter)
end
if ~isempty(zlabelstr)
    zlabel(zlabel, 'Interpreter', interpreter)
end
if ~isempty(titlestr)
    title(titlestr, 'Interpreter', interpreter) 
end
hold on;


% Add title (optional)
if isfield(options, 'title')
    title(options.title, 'Interpreter', interpreter)
end

% Add the colorbar in the style set in options struct
if strcmp(style, 'phase')
    disp('setting scalar field to phase style')
    %%%%%%%%%%%%%%%%%%%
    % Phasemap
    %%%%%%%%%%%%%%%%%%%
    colormap phasemap
    caxis([0, 2*pi])
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
    disp('setting scalar field to diverging style')
    colormap bwr
    if isfield(options, 'ylim')
        ylim(options.ylim)
    end
    set(gca, 'Position', [0 0.11 0.85 0.8]) ;
    
    % Set color axis limits
    if sscale > 0
        caxis([-sscale, sscale])
    end
    
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
