function [h1, h2, h3] = vectorFieldHeatPhaseOnImage(im, xx, yy, vx, vy, vscale, ...
    options)
%PLOTVECTORFIELDHEATPHASEONIMAGE(im, xy, vangle, speed, options)
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
%   overlay_quiver : bool
%       whether to show the quiverplot overlay or not
%   qscale : float
%       overall scale of the quivers
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

labelstr = '$|v|$ [$\mu$m / min]' ;
if isfield('label', options)
    labelstr = options.label ;
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
close all
fig = figure('units', 'normalized', ...
    'outerposition', [0 0 1 1], 'visible', 'off') ;
h1 = imshow(im) ;
hold on;
h2 = imagesc(xx, yy, vangle) ;
set(h2, 'AlphaData', speed / vscale)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUIVER 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.overlay_quiver
    qsubsample = options.qsubsample ;
    qscale = options.qscale ;
    
    vx = reshape(vx, [ww, hh]) ;
    vy = reshape(vy, [ww, hh]) ;
    QX = imresize(vx, [ww / qsubsample, hh / qsubsample], 'bicubic') ;
    QY = imresize(vy, [ww / qsubsample, hh / qsubsample], 'bicubic') ;
    xq = 1:qsubsample:ww ;
    yq = 1:qsubsample:hh ;
    [xg, yg] = meshgrid(xx(xq), yy(yq)) ;

    h3 = quiver(xg(:), yg(:), qscale * QX(:), qscale * QY(:), 0) ;
else
    h3 = [] ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phasemap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
yticklabels({'0', num2str(vscale)})
xticks([])
yyaxis left
yticks([])
cax.YAxis(1).Color = 'k';
cax.YAxis(2).Color = 'k';

% folds
% plot([foldx; foldx], [0, 0, 0; yesz, yesz, yesz], 'k--')
if isfield(options, 'outfn')
    saveas(fig, options.outfn) ;   
    close all
end
