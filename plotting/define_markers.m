function markers = define_markers(nmarkers)
%DEFINE_MARKERS return cell array of markerstyles
%   
% Parameters
% ----------
% nmarkers : int
%   The number of markers to hold in the array
%
% Returns
% -------
% markers : 1 x nmarkers cell array 
%   The symbols for linespec for each kind of marker
%
% NPMitchell 2019

markers = {'+','o','*','x','v','d','^','s','>','<'} ;
markers = markers(1:nmarkers) ;
% todo: if nmarkers > length(markers), set up a cycle here

end

