function QQ = QtensorFromAngle(angle1, keep)
% QtensorFromAngle(angle1, keep)
%
% Parameters
% ----------
% angle1 : nCells x 1 float array
%   angles of nematic director (elongation axis of ellipses, for ex)
% keep : optional nCells2Keep x 1 int array
%   indices of which cells to keep in the computation
%   Absent indices qq will have QQ(qq, :, :) = 0.
%
% Returns
% -------
% QQ : nCells x 2 x 2 float array
%   traceless symmetric matrices with nematic orientation along angle1
%
% See also
% --------
% QtensorStats.m 
% QtensorAspectTheta.m
%
%
% NPMitchell 2021


nCells = length(angle1) ;
if nargin < 2
    keep = 1:nCells ;
end
QQ = zeros(nCells, 2, 2) ;
for qq = 1:nCells
    if ~isempty(intersect(keep, qq))
        tt = mod(angle1(qq), pi) ;
        nn = [cos(tt), sin(tt)] ;
        % Create traceless symmetric matrix using unit vec
        QQ(qq, :, :) = nn' * nn - 0.5 * [1, 0; 0, 1] ;
    end
end