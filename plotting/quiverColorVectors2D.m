function [midx, hh] = quiverColorVectors2D(x, y, u, v, mags, colors, clims, lw)
%QUIVERCOLORVECTORS2D(x, y, u, v, mags, numBins)
% Color vectors of quiver plot in 2D.
% Splits mags (which are used to determine color) into #colors bins. 
% Then we make #colors different quiver objects, each with one color. 
% If mags is empty, then #colors must = length(x) = length(y) = length(u) =
% length(v), and we plot each arrow as a separate quiver object in the
% active figure.
% 
% Parameters
% ----------
% x : #vectors x 1 float array
%   vector magnitudes in y direction
% y : #vectors x 1 float array
%   vector magnitudes in y direction
% u : #vectors x 1 float array
%   vector magnitudes in x direction
% v : #vectors x 1 float array
%   vector magnitudes in y direction
%
% Returns
% -------
% midx : #colors x 1 float array
%   values of mags at which we bin and plot each set of colored arrows
% hh : quiver handle
%   handle for ONE set of vectors (the last one plotted)
%
%
% NPMitchell 2022


% Default method is to bin mags into #colors
ncolors = size(colors, 1) ;

if nargin < 7 || isempty(clims)
    clims = [min(mags), max(mags)] ;
end

if isempty(mags) 
    mags = sqrt(u.^2 + v.^2) ;
    mags = mags(:) ;
end

xedges = linspace(clims(1), clims(2), ncolors) ;
xedges(1) = -Inf ;
xedges(end) = Inf ;
[~,~,loc] = histcounts(mags, xedges);
midx = 0.5*(xedges(1:end-1)+xedges(2:end));

if nargin < 8
    lw = 3 ;
end

% Prep output if requested
if nargout > 1
   hh = cell(ncolors, 1) ;
end

% Plot each batch of vectors
hold on;
for i=1:ncolors
    idx = find(loc == i) ;
    if nargout > 1
        hh{i} = quiver(x(idx),y(idx),u(idx),v(idx),0,'color', colors(i, :), 'linewidth', lw) ;
        set(hh{i},'MaxHeadSize',1e4)
    else
        h = quiver(x(idx),y(idx),u(idx),v(idx), 0, 'color', colors(i, :), 'linewidth', lw) ;
        set(h,'MaxHeadSize',1e4)
    end
end
hold off;

   