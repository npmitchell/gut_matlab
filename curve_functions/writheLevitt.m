function [Wr, wr] = writheLevitt(xyz, closed)
%WRITHELEVITT(xyz, s) Compute the writhe by considering linesegments.
%   Wr = sum_i sum_j Omega_ij / 4pi = 2 sum_i=2..N sum_j<i Omega_ij / 4pi
%   Note that if the curve is closed, the first and last point in xyz need
%   not be identical.
%   
% Parameters
% ----------
% xyz : N x 3 float array
%   The 3D curve whose writhe we compute
% ss : N x 1 float array, optional
%   The cumulative pathlength from the endpoint of the curve xyz(1, :) to a
%   given curve point
% 
% Returns
% -------
% Wr : float 
%   The writhe of the curve as computed using the Gauss Integral
% wr : N x 1 float array 
%   The writhe "density" obtained by integrating only over ds', not over ds
% 
% NPMitchell 2019
lsegs = linesegments(xyz, closed) ;

% Consider each linesegment and compare to all others
for ii = 2:size(lsegs, 1)
    
    s1 = lseg(ii, 1:3) ;
    s2 = lseg(ii, 4:6) ;
    
    for jj = 1:(ii - 1)
        if ii == jj
            error('Should not enter this situation -> double-counting.')
        else
            s3 = lseg(jj, 1:3) ;
            s4 = lseg(jj, 4:6) ;
            r12 = s2 - s1 ;
            r13 = s3 - s1 ;
            r14 = s4 - s1 ;
            r23 = s3 - s2 ;
            r24 = s4 - s2 ;
            r34 = s4 - s3 ;
            
            % Define normalized cross products
            n1 = cross(r13, r14) ;
            n1 = n1 / vecnorm(n1) ;
            n2 = cross(r14, r24) ;
            n2 = n2 / vecnorm(n2) ;
            n3 = cross(r24, r23) ;
            n3 = n3 / vecnorm(n3) ;
            n4 = cross(r23, r13) ;
            n4 = n4 / vecnorm(n4) ;
            
            % Sum 
            ostar1 = arcsin(dot(n1, n2)) ;
            ostar2 = arcsin(dot(n2, n3)) ;
            ostar3 = arcsin(dot(n3, n4)) ;
            ostar4 = arcsin(dot(n1, n4)) ;
            ostar = ostar1 + ostar2 + ostar3 + ostar4 ;
            r34xr12dr13 = dot(cross(r34, r12), r13) ;
            Omega = ostar * sign(r34xr12dr13) ;
        end
    end
end
Wr = nansum(Omega) / (2 * pi) ; % Note that there is a factor of 2 / 4pi
end

