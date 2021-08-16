function [midx, mediany, meany, stdy, ny] = binDataMedianStd(x, y, xedges, weights)
% binDataMeanStd(x, y, xedges)
% bin data in x, take medians and stdevs of y data and output binned mean 
% and stdev curves
%
%
% Parameters
% ----------
% x : values for binning into xedges
% y : data values to take stats of
% xedges : bin edges for binning x values and grouping y values
% weights : optional, for scaling means and stdevs
%
% Returns
% -------
% midx : middle value of each bin (average between adjacent xedges)
% mediany : median y value in each bin
% meany : mean y value in each bin
% stdy : standard deviations of y values in each bin
%
% NPMitchell 2021

% Default edges are positive integers
if nargin < 3
    xedges = linspace(floor(min(x)), ceil(max(x)), ...
        ceil(max(x)) - floor(min(x)) + 1)' ;
end

if nargin < 4
    weights = ones(size(x)) ;
end

[~, ~, loc] = histcounts(x,xedges);
% Consider each bin
nBins = length(xedges) -1 ;
mediany = nan(nBins, 1) ;
meany = nan(nBins, 1) ; 
stdy = nan(nBins, 1) ;
ny = zeros(nBins, 1) ;

for bin = 1:nBins
    idx = find(loc == bin & ~isnan(weights)) ;
    if ~isempty(idx)
        ny(bin) = length(idx) ;
        ww = weights(idx) / nansum(weights(idx)) ;
        mediany(bin) = nanmedian(y(idx)) ;
        meany(bin) = nansum(ww .* y(idx)) ;

        % The var function divides by N instead of (N-1)
        % stdy(bin) = sqrt(var(x(idx), ww)) ;

        % We choose to divide by (N-1)
        num = nansum(ww .* (y(idx)-meany(bin)).^2) ;
        % denom is very close to 1, slightly smaller because of N-1
        denom = (length(idx) -1) / length(idx) * nansum(ww) ;
        stdy(bin) = sqrt(num / denom) ;
        
        if isnan(stdy(bin))
            assert(length(idx) == 1)
        end
    end
end

% Unweighted for reference:
% meany = accumarray(loc(:),y(:))./accumarray(loc(:),1);
% stdy = accumarray(loc(:), y(:), [], @std);

midx = 0.5*(xedges(1:end-1)+xedges(2:end));

% Check it
% scatter(x, y, 5, 'filled')
% hold on; 
% plot(midx, meany, '.-')
% plot(midx, meany-stdy); 
% plot(midx, meany+stdy); 
