function [mss, mcline, radii_from_mean, avgpts_ss, avgpts] = srFromDVCurves(curvesDV)
%SRFROMDVCURVES(curvesDV) Find new centerline, pathlength, and radii for 3D DV hoops 
%   Given the array of 1d DV curves, compute the new centerline from
%   averaging each hoop and resampling the string of averages. Also returns
%   the pathlength of the resampling and distance from the average points.
%
% Parameters
% ----------
% curvesDV : N*M x 3 float
%   The 3D embedding square grid (UV) of the 2D DV curves
% 
% Returns
% -------
% mss : 1000 x 1 float
%   pathlength values of centerline constructed via mean of hoops
% mcline : 1000 x 3 float
%   Resampled centerline created from avgpts (which is average 3d positions 
%   of hoops)
% radii_from_mean : N*M x 1 float
%   distance of each gridpt (on hoops) to averaged hoop points (avgpts)
% avgpts_ss : N x 1 float
%   pathlength values of the centerline made by connecting average of each
%   hoop
% avgpts : N x 3 float
%   positions of the average of each hoop (3d positions of each hoop)
%
% NPMitchell 2019
                   
% Iterate over u coordinate
% Compute radius of curves3d from new centerline
avgpts = zeros(size(curvesDV, 1), 3) ;
radii_from_mean = zeros(size(curvesDV, 1), size(curvesDV, 2)) ;
for jj = 1:size(curvesDV, 1)
    % Consider this hoop
    hoop = squeeze(curvesDV(jj, :, :)) ;
    avgpts(jj, :) = mean(hoop) ; 
    radii_from_mean(jj, :) = vecnorm(hoop - avgpts(jj, :), 2, 2) ;
end
avgpts_ss = ss_from_xyz(avgpts) ;

% Now define new centerline built from avgpts
fprintf('New centerline from avgpts...\n')
mcline = curvspace(avgpts, 1000) ;
mss = ss_from_xyz(mcline) ;
% avgpts_ss = mss(pointMatch(avgpts, mcline)) ;

end

