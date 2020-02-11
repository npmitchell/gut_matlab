function [v0n, v0t, v0t2d, jac, facenormals] = ...
    resolveTangentNormalVelocities(faces, vertices, v0, vertices2d, ...
    fieldfaces, varargin)
%RESOLVETANGENTNORMALVELOCITIES
% Resolve a 3d vector field into the tangential and normal components of a
% 2d mesh in its embedding and in its 2d image (flattened mesh)
%
% Parameters
% ----------
% faces : #faces x 3 int array
%   indices into vertices2d of the mesh faces
% vertices : #vertices x 3 float array
%   mesh embedding in 3d
% v0 : N x 3 float array
%   velocities in 3d evaluated at points living in fieldfaces
% vertices2d : #vertices x 2 float array
%   mesh embedding in 2d
% fieldfaces : N x 1 int array 
%   indices into faces in which velocities v0 are evaluated
% varargin : keyword arguments
%   'facenormals' : #faces x 3 float array, normal vectors for each face
%
% Returns
% -------
% v0n : 
% v0t : 
%   original scale velocities in tangent plane
% vfield2d : #fieldfaces x 2 float array
%   The vector field mapped to the 2d mesh
% jac : length(ff) x 1 cell array
%   A cell array containing the jacobian for each face as each element -- 
%   transformation from 3d to 2d
%
% NPMitchell 2020

% Obtain face normals, either from varargin or compute them
if ~isempty(varargin)
    for i = 1:length(varargin)
        if isa(varargin{i},'double') 
            continue;
        end
        if isa(varargin{i},'logical')
            continue;
        end

        if ~isempty(regexp(varargin{i},'^[Ff]ace[Nn]ormals','match'))
            facenormals = varargin{i+1} ;
        end
    end
else
    % Option 1 : compute using matlab built-in
    facenormals = faceNormal( triangulation(faces, vertices) );
    % Option 2 : project faces onto the already-computed normals
    % aux_alternate_velocity_projection
end


%% Take dot product of flow fields with normals in 3D
v0n = dot(facenormals(fieldfaces, :), v0, 2) ;
% Subtract the normal velocity to obtain tangential velocity
v0t = v0 - v0n .* facenormals(fieldfaces, :) ;

%% Compute the tangential velocities in plane
% u is 3d, w is 2d. jac takes u->w, jjac takes w->u
[v0t2d, jac] = pullVectorField3Dto2DMesh(v0t, vertices2d, vertices, faces, fieldfaces) ;
