function [ Habs, Hn ] = calculate_mean_curvature( face, vertex )
%CALCULATE_MEAN_CURVATURE Calculates the discrete mean curvature of a
%triangulated surface.  Here the mean curvature is defined on a mesh vertex
%with ID v from:
%   L * x(v) = -2 H(v) n(v)
% Where x(v) is the 3D vertex position, n(v) is the vertex mean curvature
% unit normal, L is the Laplace-Beltrami operator, and H(v) is the vertex
% mean curvature
%
%   Input Parameters:
%       - face:             The triangulation connectivity list
%       - vertex:           The vertex coordinate list
%
%   Ouput Parameter:
%       - Habs:             The absolute mean curvature at each vertex
%       - Hn:               The mean curvature vector ( -2 H n )
%
%   by Dillon Cislo 02/26/2019

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------

if ( numel(size(face)) ~= 2 )
    error( 'Connectivity list is improperly sized!' );
end

if ( numel(size(vertex)) ~= 2 )
    error( 'Vertex coordinate list is improperly sized!' );
end

if ( size(face,2) ~= 3 )
    face = face';
    if ( size(face,2) ~= 3 )
        error('Mesh facets must be triangles!');
    end
end

if ( size(vertex,2) ~= 3 )
    vertex = vertex';
    if ( size(vertex,2) ~= 3 )
        error('Vertex coordinates must be 3D!');
    end
end

%--------------------------------------------------------------------------
% Calculate mean curvature
%--------------------------------------------------------------------------

% The Laplace-Beltrami operator
L = construct_laplace_beltrami( face, vertex );

% The mean curvature vector
Hn = L * vertex;

% The absolute mean curvature
Habs = sqrt( sum( Hn.^2, 2 ) ) ./ 2;

end

