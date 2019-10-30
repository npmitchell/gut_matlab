function [ path ] = ...
    findShortestPath(faceIn, vertexIn, cp1, cp2)
%CYLINDERCUTMESH creates a cut mesh structure from a an input mesh.  
% Input mesh should be a topological cylinder.  
% Output mesh will be a topological disk.  
% The mesh is cut along an edge-based mesh geodesic between two input
% points.  Currently this method can only support a single cut path
%   INPUT PARAMETERS:
%       - faceIn:         #Fx3 face connectivity list
%       - vertexIn:       #VxD vertex coordinate list
%       - normalIn:       #VxD vertex normal list
%       - cp1:            Vertex ID of the cut path origin
%       - cp2:            Vertex ID of the cut path termination
%   
%   OUTPUT PARAMETERS:
%       - cutMesh:      An ImSAnE-style mesh struct with additional fields
%                       to describe the cut mesh properties
%       - cp1Out:       Vertex ID of the final cut path origin
%       - cp2Out:       Vertex ID of the final cut path termination
%       - P:            shortest path, indices of original vertices

%==========================================================================
% Calculate the Shortest Path Between Input Points Along Mesh Edges
%==========================================================================

% MATLAB-style triangulation
meshTri = triangulation( faceIn, vertexIn );

% The vertex IDs of vertices on the mesh boundary
bdyIDx = meshTri.freeBoundary;
bdyIDx = unique(bdyIDx(:));

% The #Ex3 edge connectivity list of the mesh
edgeIn = edges( meshTri );

% Check that the input mesh is a topological cylinder
if ( length(vertexIn) - length(edgeIn) + length(faceIn) ) ~= 0
    error( 'Input mesh is NOT a topological cylinder!' );
end

% Check that the origin/terminal points lie on the mes
if ~ismember(cp1, bdyIDx) || ~ismember(cp2, bdyIDx)
    error('One or more input point does not lie on the mesh boundary!');
end

% Find edge lengths
L = vertexIn( edgeIn(:,2), : ) - vertexIn( edgeIn(:,1), : );
L = sqrt( sum( L.^2, 2 ) );

% Construct an #V x #V vertex adjacency matrix
A = sparse( [ edgeIn(:,1); edgeIn(:,2) ], ...
    [ edgeIn(:,2), edgeIn(:,1) ], ...
    [ L; L ], size(vertexIn, 1), size(vertexIn, 1) );

% A MATLAB-style weighted, undirected graph representation of the mesh
% triangulation.
G = graph(A);

% The shortest path
path = shortestpath( G, cp1, cp2 )' ;