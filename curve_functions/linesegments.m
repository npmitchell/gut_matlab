function lsegs = linesegments(xyz, closed)
%LINESEGMENTS Collate linesegments from an ordered list of points
%   The ith row is the segment from i to i+1, as xi yi zi xi+1 yi+1 zi+1.
%   The argument can be any dimensionality (N x D).
%
% Parameters
% ----------
% xyz : N x D float array
%   The D-dimensional curve whose linesegments we find
% closed : bool
%   Whether the curve is closed (true) or open (false)
% 
% Returns
% -------
% lsegs : N x 2D or (N-1) x 2D float array
%   Each row is a linesegment. For ex, in 3D, a row is x1 y1 z1 x2 y2 z2
%
% NPMitchell 2019

if closed
    lsegs = [xyz, circshift(xyz, -1, 1)] ;
else
    lsegs = [xyz(1:end-1, :), xyz(2:end, :)] ;
end
end

