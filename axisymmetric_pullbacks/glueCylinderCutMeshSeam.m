function [cutMeshClosed, glue2cut] = glueCylinderCutMeshSeam(cutMesh)
%GLUECUTMESHSEAM(cutMesh)
% close a rectilinear cutMesh by gluing the seam back together 
% 
% Parameters
% ----------
% cutMesh : struct
%   the cylinderCutMesh object with fields
%       nU : int
%       nV : int
%       v : (nU*nV) x 3 float array
%       f : (nU-1)*(nV-1) x 3 int array
%       vn : (nU*nV) x 3 float array
%       u : nU*nV
% 
% Returns
% -------
% cutMeshClosed : struct
% glue2cut : #vertices in cutMesh x 1 int
%   glue2cut(i) gives the index in output glued mesh of the ith vertex in 
%   the cut mesh, so cutVtx = glueVtx(glue2cut)
% 
%
% NPMitchell 2020
%

% todo: assert that the vertices have a rectilinear structure with nU, nV
nU = cutMesh.nU ;
nV = cutMesh.nV ;

cutMeshClosed.v = cutMesh.v(1:end-nU, :) ;
cutMeshClosed.f = mod(cutMesh.f, (nV-1)*nU + 1) ;
cutMeshClosed.f(cutMesh.f > (nV-1)*nU) = cutMeshC.f(cutMesh.f > (nV-1)*nU) + 1 ;
cutMeshClosed.vn = cutMesh.vn(1:end-nU, :) ;
cutMeshClosed.u = cutMesh.u(1:end-nU, :) ;

if length(varargout) > 1
    % produce a map from the glued mesh back to the cut mesh indices, 
    % so that cutVtx = glueVtx(glue2cut)
    glue2cut = [ 1:(length(cutMesh.u(:, 1))-nU), 1:nU ] ;
end
