% npm 2018
%
% TRI is a (#tris x 3) array with one row for each triangle in the
%     tessalation. The elements of the rows are the inds of the particles in
%     that tri, in ascending order. If particle N is in 6 tris, it will most
%     likely appear in each column twice, once as the lowest index, once as the
%     middle, and once as the highest.
%
% xyz is a (# nodes x 3) float array, 
%     each row is a position vector of a node
%
% sign is a (#tris x 1) float array or a single float (+1 or -1). 
%     If sign is 1, then normals of each triangle in TRI are individually 
%     flipped so that their z component is positive. 
%     If sign is -1, then z components are forced to be negative.
%     If sign is an array, then each element determines whether or not each 
%     triangle is flipped.
%
% normed is a boolean
%     If true, the output normals have unit length. Otherwise, the length
%     is determined by the area of each triangle --- and that is determined
%     by the cross product magnitude, which can be useful for
%     applying a pressure (force per unit area of the triangle).
%
% normals is a (#particles x 3) array giving the normal vectors at each
%     node of a network in xyz. If a particle in xyz has no bonds, its normal
%     will be zero.
%

function [normals] = TRI2normals(TRI, xyz, sign, normed)

normals = 0. * xyz;
col1 = TRI(:, 1) ;    
col2 = TRI(:, 2) ;
col3 = TRI(:, 3) ;
crossprod = cross(xyz(col2, :) - xyz(col1, :), xyz(col3, :) - xyz(col1, :)) ;

if sign == -1
    crossprod(crossprod(:, 3) > 0, :) = crossprod(crossprod(:, 3) > 0, :) * -1;
elseif sign == 1
    crossprod(crossprod(:, 3) < 0, :) = crossprod(crossprod(:, 3) < 0, :) * -1;
else
    crossprod = transpose(sign) .* crossprod;
end

% nns are normals for each triangle
% Note that here vecnorm takes arguments of p=2 for p-norm, and dim=2
if normed
    nns = crossprod ./ vecnorm(crossprod, 2, 2) ;
  
    % Translate nns into normals (change size of array to fit xyz)
    normals(col1, :) = normals(col1, :) + nns ;
    normals(col2, :) = normals(col2, :) + nns ;
    normals(col3, :) = normals(col3, :) + nns ;
else
    % Translate nns into normals (change size of array to fit xyz)
    normals(col1, :) = normals(col1, :) + 0.5 * crossprod ;
    normals(col2, :) = normals(col2, :) + 0.5 * crossprod ;
    normals(col3, :) = normals(col3, :) + 0.5 * crossprod ;

end

return

