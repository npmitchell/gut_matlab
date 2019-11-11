function [ Graph ] = makeGraph(faceIn, vertexIn)
%CYLINDERCUTMESH creates a MATLAB-style weighted, undirected graph 
% representation of the mesh triangulation.
% 
%   INPUT PARAMETERS:
%       - faceIn:         #Fx3 face connectivity list
%       - vertexIn:       #VxD vertex coordinate list
%       - cp1:            Vertex ID of the cut path origin
%       - cp2:            Vertex ID of the cut path termination
%   
%   OUTPUT PARAMETERS:
%       - Graph:          undirected graph rep of mesh triangulation

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

% % Check that the input mesh is a topological cylinder
% if ( length(vertexIn) - length(edgeIn) + length(faceIn) ) ~= 0
%     error( 'Input mesh is NOT a topological cylinder!' );
% end

% % Check that the origin/terminal points lie on the mesh
% if ~ismember(cp1, bdyIDx) || ~ismember(cp2, bdyIDx)
%     error('One or more input point does not lie on the mesh boundary!');
% end

% Find edge lengths
L = vertexIn( edgeIn(:,2), : ) - vertexIn( edgeIn(:,1), : );
L = sqrt( sum( L.^2, 2 ) );

% Construct an #V x #V vertex adjacency matrix
A = sparse( [ edgeIn(:,1); edgeIn(:,2) ], ...
    [ edgeIn(:,2), edgeIn(:,1) ], ...
    [ L; L ], size(vertexIn, 1), size(vertexIn, 1) );

% A MATLAB-style weighted, undirected graph representation of the mesh
% triangulation.
Graph = graph(A);

% Can now compute, for ex, The shortest path
% path = shortestpath( G, cp1, cp2 )' ;
