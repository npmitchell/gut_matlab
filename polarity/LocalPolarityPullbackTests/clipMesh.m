function [newFace, newFaceIDx, newVertex, newVertexIDx] = ...
    clipMesh(face, vertex, rmAllTypes)
% CLIPMESH Remove dangling triangles, ears, and unreferenced vertices from 
% a mesh triangulation
% 
% INPUT PARAMETERS:
%
%   face:         #Fx3 input face connectivity list
%
%   vertex:       #VxD input vertex coordinate list
%
%   rmAllTypes:   If true, this algorithm will remove dangling triangles,
%                 ears, and unreferenced vertices.  If false, it will only
%                 remove unreferenced vertices.
%
% OUTPUT PARAMETERS:
%
%   newFace:      #F'x3 output face connectivity list
%
%   newFaceIDx:   #F'x1 indices of new faces in old triangulation
%
%   newVertex:    #V'xD output vertex coordinate list
%
%   newVertexIDx: #V'x1 indices of new vertices in old triangulation
% 
% By Dillon Cislo 02/04/2020

%--------------------------------------------------------------------------
% Validate Inputs
%--------------------------------------------------------------------------

validateattributes( vertex, {'numeric'}, ...
    {'2d', 'finite', 'real', 'nonnan'} );
validateattributes( face, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'real', 'integer', 'positive'} );

assert(all(ismember(face, (1:size(vertex,1)).'), 'all'), ...
    'Face connectivity list contains unreferenced vertices!');

if nargin < 3
    rmAllTypes = true;
end

%--------------------------------------------------------------------------
% Remove Undesired Vertices
%--------------------------------------------------------------------------

% Lone points are defined as vertices belonging to ears (boundary triangles
% sharing only one edge with the mesh) or vertices belonging to triangles
% sharing only one vertex with the edge of the mesh or vertices that are
% not referenced by the input triangulation
num_lone_points = Inf;

% Make temporary copies of the input parameters
f = face;
v = vertex;

while( num_lone_points ~= 0 )
    
    % Determine the number of faces attached to each vertex
    vertexFaceCount = full(sparse(f, 1, 1, size(v,1), 1));
    
    % The number and ID of vertices to remove
    if rmAllTypes
        lone_points = vertexFaceCount < 2;
    else
        lone_points = vertexFaceCount < 1;
    end
    num_lone_points = sum(lone_points);
    lone_points = find(lone_points);
    
    % Remove the lone vertices at the specified index values
    newVertex = v;
    newVertex(lone_points, :) = [];
    
    % Create a map from old vertex indices to update vertex indices
    % mapVIDx(i) gives the new vertex index of the ith vertex in the
    % original mesh ("0" if the vertex was not kept)
	% NOTE: This assumes there are no duplicate vertices/faces
    [~, mapVIDx] = ismember(v, newVertex, 'rows');
    
    % Find any faces that contained the lone vertices and remove them
    newFace = f;
    lone_faces = any( ismember( newFace, lone_points ), 2 );
    newFace( lone_faces, : ) = [];
    
    % Now update the vertex indices in the connectivity list
    newFace = [ mapVIDx(newFace(:,1)), mapVIDx(newFace(:,2)), ...
        mapVIDx(newFace(:,3)) ];
    
    % Update the temporary copies
    f = newFace;
    v = newVertex;
    
end

% if isempty(newVertex)
%     warning('All vertices have been removed!');
% end
% 
% if isempty(newFace)
%     warning('All faces have been removed!');
% end

% Find the indices of the new vertices in the old triangulation
% NOTE: This assumes there are no duplicate vertices
[~, newVertexIDx] = ismember( newVertex, vertex, 'rows' );

% Find the indices of the new faces in the old triangulation
% NOTE: This assumes there are no duplicate faces
newFaceIDx = zeros( size(newFace, 1), 1 );
for fID = 1:size(newFace,1)
    
    f = newVertexIDx( newFace(fID, : ) ).';
    oldFID = find(ismember(f, face, 'rows'));
    
    newFaceIDx(fID) = oldFID;
    
end

end


    
