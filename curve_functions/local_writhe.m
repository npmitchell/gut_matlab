function wr = local_writhe(ss, xyz)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% Based on Berger & Prior, The writhe of open and closed curves, 2006
% 
% Note that in the paper, prime denotes d/dz.
lamb = dz /ds ;
% take cross product of tangent with change of tangent
txtprime = cross(tangent, tanprime) ;
wr = (1 / (2 * pi)) * (1 / (1 + abs(lamb))) * txtprime(:, 3) ; 

end

