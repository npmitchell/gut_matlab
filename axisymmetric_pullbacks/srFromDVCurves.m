function [mss, mcline, radii_from_mean, avgpts_ss, avgpts] = srFromDVCurves(curvesDV)
%srFromDVCurves Find pathlength match and radius for hoops in 3d 
%   Given the centerline pathlength cseg_ss, the centerline cseg, and the
%   array of 1d DV curves, compute the matching centerline pathlength for
%   each 
%
% Parameters
% ----------
% curvesDV : N*M x 3 float
%   The 3D embedding square grid (UV) of the 2D DV curves
% 
% Returns
% -------
% ssv : N*M x 1 float
%   pathlength values of supplied centerline
% radii : N*M x 1 float
%   distance of each gridpt (on hoops) to pointmatch on supplied
%   centerline
% avgpts : N x 3 float
%   positions of the average of each hoop
% cids : N x 1 int
%   indices into supplied centerline that are closest on average to
%   points in each hoop
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

