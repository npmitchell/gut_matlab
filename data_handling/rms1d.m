function out = rms1d(arr)
% out = rms1d(arr) 
% Find root mean square value in array arr.
% 
% todo: handle complex variable case
%
% NPM 2021

out = sqrt((1/length(arr)) * sum(arr.^2)) ;
end