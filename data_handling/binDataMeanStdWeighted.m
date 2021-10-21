function [midx, meany, stdy, ny, stey] = binDataMeanStdWeighted(...
    x, y, xedges, weights, nsamplesForStE, minstds)
% binDataMeanStdWeighted(x, y, xedges, weights)
% bin data in x, take weighted means and stdevs of y data and output binned mean and
% weighted stdev curves
% Weighted stdev can be defined in various ways. See, for ex,
% https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf
% or 
% https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Mathematical_definition
% Here we use MATLAB's var(x, w) with weights. Note that differs from std()
% because MATLAB seems not to divide by N-1 but instead by N for N
% measurements
% 
% 
% Parameters
% ----------
% x : values for binning into xedges
% y : data values to take stats of
% xedges : bin edges for binning x values and grouping y values
% weights : weights for each observation x that bias the importance of each
% observation. Note that these need not be normalized. Each bin's weights
% will be normalized
% nsamplesForStE : optional (default=500)
% minstds : if only one observation k lies in a bin, ascribe a stdev to
% that bin of minstds(k).
% 
% Returns
% -------
% midx : middle value of each bin (average between adjacent xedges)
% meany : mean y value in each bin
% stdy : standard deviations of y values in each bin
% ny : number in each bin
% stey : standard error on the mean determined via 'simple' bootstrapping
%
% See also
% --------
% weightedMeanStdSte.m
%
% NPMitchell 2021

% Default edges are positive integers
if nargin < 5
    nsamplesForStE = 500 ;
end
if nargin < 3
    xedges = linspace(floor(min(x)), ceil(max(x)), ...
        ceil(max(x)) - floor(min(x)) + 1)' ;
end
if nargin < 4 
    weights = 1 / length(x) ;
end

[~, ~, loc] = histcounts(x,xedges);
% Consider each bin
nBins = length(xedges) -1 ;
meany = nan(nBins, 1) ; 
stdy = nan(nBins, 1) ;
stey = nan(nBins, 1) ;
ny = zeros(nBins, 1) ;

for bin = 1:nBins
    idx = find(loc == bin & ~isnan(weights)) ;
    if ~isempty(idx)
        % see also: weightedMeanStdSte.m
        
        ny(bin) = length(idx) ;
        ww = weights(idx) / nansum(weights(idx)) ;
        meany(bin) = nansum(ww .* y(idx)) ;

        % The var function divides by N instead of (N-1)
        % stdy(bin) = sqrt(var(x(idx), ww)) ;

        % We choose to divide by (N-1)
        num = nansum(ww .* (y(idx)-meany(bin)).^2) ;
        % denom is very close to 1, slightly smaller because of N-1
        denom = (length(idx) -1) / length(idx) * nansum(ww) ;
        stdy(bin) = sqrt(num / denom) ;
        
        % Lower limit to std if supplied
        if ny(bin) == 1 && exist('minstds', 'var')
            stdy(bin) = minstds(idx) ;
        end
        
        if isnan(stdy(bin))
            assert(length(idx) == 1)
            stey(bin) = NaN ;
        else
            % bootstrap for sterror if desired as output var
            if nargout > 4
                opts = struct('simpleOrFull', 'simple', 'nsamples', nsamplesForStE) ;
                
                stey(bin) = bootstrapErrorWithWeights(y(idx), ww, opts) ;
            end
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
