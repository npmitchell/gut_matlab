function [midx, meany, stdy] = binDataMeanStd(x, y, xedges)
% binDataMeanStd(x, y, xedges)
% bin data in x, take means and stdevs of y data and output binned mean and
% stdev curves
%
%
% Parameters
% ----------
%
% Returns
% -------
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

% Check it
% scatter(x, y, 5, 'filled')
% hold on; 
% plot(midx, meany, '.-')
% plot(midx, meany-stdy); 
% plot(midx, meany+stdy); 
