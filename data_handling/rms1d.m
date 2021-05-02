function out = rms1d(arr, include_nan)
% out = rms1d(arr) 
% Find root mean square value in array arr.
% 
% todo: handle complex variable case
%
% NPM 2021

if nargin < 2
    include_nan = false;
end

if include_nan
    out = sqrt((1/length(arr)) * sum(arr.^2)) ;
else
    out = sqrt((1/sum(~isnan(arr))) * nansum(arr.^2)) ;    
end