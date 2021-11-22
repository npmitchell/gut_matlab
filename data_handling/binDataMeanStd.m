function [midx, meany, stdy, ny, stey] = binDataMeanStd(x, y, xedges)
% binDataMeanStd(x, y, xedges)
% bin data in x, take means and stdevs of y data and output binned mean and
% stdev curves
%
%
% Parameters
% ----------
% x : values for binning into xedges
% y : data values to take stats of
% xedges : bin edges for binning x values and grouping y values
%
% Returns
% -------
% midx : middle value of each bin (average between adjacent xedges)
% meany : mean y value in each bin
% stdy : standard deviations of y values in each bin
%
% NPMitchell 2021

% Default edges are positive integers
if nargin < 3
    xedges = linspace(floor(min(x)), ceil(max(x)), ...
        ceil(max(x)) - floor(min(x)) + 1)' ;
end

[~,~,loc]=histcounts(x,xedges);
meany = accumarray(loc(:),y(:))./accumarray(loc(:),1);
midx = 0.5*(xedges(1:end-1)+xedges(2:end));
stdy = accumarray(loc(:), y(:), [], @std);


% Stanadard error on the mean
if nargout > 3
    nBins = length(xedges) -1 ;
    ny = zeros(nBins, 1) ;
    for bin = 1:nBins
        idx = find(loc == bin) ;
        if ~isempty(idx)
            % see also: weightedMeanStdSte.m

            ny(bin) = length(idx) ;
        end
    end
    stey = stdy ./ sqrt(ny) ;
end


% Check it
% scatter(x, y, 5, 'filled')
% hold on; 
% plot(midx, meany, '.-')
% plot(midx, meany-stdy); 
% plot(midx, meany+stdy); 
