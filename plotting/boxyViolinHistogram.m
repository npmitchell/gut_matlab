function boxyViolinHistogram(xvals, yvals, nBins_or_edges, maxBoxWidth, ...
    style, faceColors, edgeColors, lineWidths)
% boxyViolinHistogram(xvals, yvals, boxWidths, ...
%    boxHeights, style, faceColors, edgeColors, lineWidths)
%
% Parameters
% ----------
% xvals : length(N) numeric
%   x positions of each violin box
% yvals : length(N) numeric (default = [1])
%   y positions of each violin box
% nBins_or_edges: int or Px1 numeric (default = 100)
%   number of bins into which to histogram the data or else edges around
%   which to histogram
% maxBoxWidth : optional 1x1 numeric (default = 1)
%   maximum width of any box in the plot
% style : optional string specifier ('left', 'center', default ='center')
%   if left, uses bottom left for (xvals, yvals)
%   if center, interprets (xvals, yvals) as center of box
% faceColors : optional #unique(xvals) x 3 numeric or #unique(xvals) x 1 cell array
% 	RGB triplets to color each violin histogram   
% edgeColors : optional #unique(xvals) x 3 numeric or #unique(xvals) x 1 cell array
% 	RGB triplets to color each violin histogram    
% LineWidths : optional #unique(xvals) x 1 numeric (default=0)
% 	linewidth for each histogram
%
% Example usage
% -------------
% N = 1000 ;
% x = ones(N, 1) ;
% y = rand(N, 1) ;
% boxyViolinPlot(x, y)
% 
% NPMitchell 2021


% Argument parsing
if nargin < 8 || isempty(lineWidths)
    lineWidths = 1 ;
end

if nargin < 7 || isempty(edgeColors)
    edgeColors = 'k' ;
end

if nargin < 6 || isempty(faceColors)
    faceColors = [  0        0.447        0.741] ;
elseif any(faceColors > 1+1e-14)
    faceColors = double(faceColors) ./ 255.0 ;
end

if nargin < 5 || isempty(style)
    style = 'center' ;
end

if nargin < 4  || isempty(maxBoxWidth)
    % try something reasonable for widths
    maxBoxWidth = 1 ;
end


% Input handling -- allow for single values in any of the supplied inputs
% if there are multiple xvals or multiple yvals
assert(numel(maxBoxWidth) == 1)
if numel(xvals) == 1 && numel(yvals) == 1
    disp('Only one x,y value given to boxyViolinHistogram')
elseif numel(xvals) == 1
    xvals = xvals * ones(size(yvals)) ;
elseif numel(yvals) == 1
    yvals = yvals * ones(size(xvals)) ;
end


% done input handling/parsing

% Discretize into histograms for each xvalue
ux = unique(xvals) ;

% Match size of unique(xvals) for faceColors, edgeColors, lineWidths
if ischar(faceColors) && ~iscell(faceColors) 
    faceColors = repmat({faceColors}, numel(ux),1) ;
elseif numel(faceColors) == 3
    faceColors = faceColors .* ones(numel(ux), 3) ;
end
if ischar(edgeColors) && ~iscell(edgeColors) 
    edgeColors = repmat({edgeColors}, numel(ux),1) ;
elseif numel(edgeColors) == 3
    edgeColors = edgeColors .* ones(numel(ux), 3) ;
end
if numel(lineWidths) == 1
    lineWidths = lineWidths .* ones(size(ux)) ;
end

% Discretize into histograms for each xvalue
if all(numel(xvals) == numel(yvals))
    Ncount = cell(length(ux), 1) ;
    edges = cell(length(ux), 1) ;
    Nmax = 0 ;
    for kk = 1:length(unique(xvals))
        [Ncount{kk}, edges{kk}] = histcounts(yvals(xvals==ux(kk)), nBins_or_edges) ;
        Nmax = max(Nmax, max(Ncount{kk})) ;
    end
elseif numel(xvals) < 2
    disp('There is only one x value -- there should be only one yvalue too')
    Ncount = cell(1, 1) ;
    edges = cell(1, 1) ;
    [Ncount{1}, edges{1}] = histcounts(yvals, nBins_or_edges) ;
    Nmax = max(Ncount{1}) ;
else
    error('Could not parse xvals')
end

for kk = 1:length(ux)
    % interpret current face color
    if iscell(faceColors)
        fColor = faceColors{kk} ;
    else
        fColor = faceColors(kk, :) ;
    end
    if iscell(edgeColors)
        eColor = edgeColors{kk} ;
    else
        eColor = edgeColors(kk, :) ;
    end
    
    boxWidths = Ncount{kk} * maxBoxWidth / Nmax ;
    boxHeights = diff(edges{kk}) ;
    yv = 0.5 * (edges{kk}(1:end-1) + edges{kk}(2:end)) ;
    
    
    % Plot the box
    for qq = 1:length(yv)
        switch style
            case 'left' 
                rectangle('Position',[ux(kk), yv(qq), ...
                    boxWidths(qq), boxHeights(qq)],...
                    'FaceColor', fColor,...
                    'EdgeColor', eColor,...
                    'LineWidth',lineWidths(kk))
            case 'center'
                xc = ux(kk) -0.5*boxWidths(qq) ;
                yc =  yv(qq) - 0.5*boxHeights(qq) ;
                rectangle('Position',...
                    [xc, yc, ...
                    boxWidths(qq), boxHeights(qq)],...
                    'FaceColor', fColor,...
                    'EdgeColor', eColor,...
                    'LineWidth',lineWidths(kk))
        end
    end
end

