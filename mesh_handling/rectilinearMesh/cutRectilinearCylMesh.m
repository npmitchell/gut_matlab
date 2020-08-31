function cutMesh = cutRectilinearCylMesh(mesh, options)
%cutRectilinearMesh(mesh, options)
%
% Parameters
% ----------
% mesh : struct, closed cylinder mesh with fields
%   nU : int
%   v : (nU*(nV-1)) x 3 float array
%       3d vertices of the mesh embedding
%   u : (nU*(nV-1)) x 2 float array
%       2d vertices of the rectilinear mesh in pullback space
%   f : #faces x 3 int array
%       indices into v (or equivalently into u) of mesh connectivity
%       (faces)
%   and optional fields
%       vn : (nU*nV) x 3 float array
%           vertex normals 
% options : optional struct with fields
%   
% 
% Returns 
% -------
% cutMesh : struct with fields
%   nU : int
%   v : (nU*nV) x 3 float array
%       3d vertices of the mesh embedding
%   u : (nU*nV) x 2 float array
%       2d vertices of the rectilinear mesh in pullback space
%   f : #faces x 3 int array
%       indices into v (or equivalently into u) of mesh connectivity
%       (faces)
% If input mesh contains field vn, then output will also:
%   vn : (nU*nV) x 3 float array
%       vertex normals 
%
% NPMitchell 2020

if nargin > 1
    if isfield(options, 'vmax')
        vmax = options.vmax ;
    else
        vmax = 1 ;
    end
else
    vmax = 1 ;
end

nU = mesh.nU ;
nV = length(mesh.v(:, 1)) / nU + 1;
cutMesh = mesh ;

% Duplicate the first row as the last row
cutMesh.v(nU*(nV-1) + 1:nU*nV, :) = mesh.v(1:nU, :) ;

% Duplicate the first row as last in pullback space
cutMesh.u(nU*(nV-1) + 1:nU*nV, :) = [0, vmax] + mesh.u(1:nU, :) ;
cutMesh.f = defineFacesRectilinearGrid(mesh.u, nU, nV) ;

if isfield(mesh, 'vn')
    cutMesh.vn((nV-1)*nU + 1:nU*nV, :) = cutMesh.vn(1:nU, :) ;
end
