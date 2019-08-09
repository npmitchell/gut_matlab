function [ K, Kr ] = calculate_gaussian_curvature( face, vertex )
%CALCULATE_GAUSSIAN_CURVATURE Calculates the discrete Gaussian curvature of
%a triangulated surface. Here the Gaussian curvature is defined on a mesh
%vertex with ID v as:
%
%   K(v) = ( 2*pi - sum( theta ) ) / A(v) ( Bulk Vertex )
%        = (   pi - sum( theta ) ) / A(v) ( Boundary Vertex )
%
% Where the sum( theta ) denotes a sum of the angles of the incident
% triangles at vertex v and A(v) denotes the mixed Voronoi averaging area.
%
%   Input Parameters:
%       - face:             The triangulation connectivity list
%       - vertex:           The vertex coordinate list
%
%   Ouput Parameter:
%       - K:                The Gaussian curvature at each vertex
%       - Kr:               The 'raw' Gaussian curvature ( i.e. no
%                           normalization by A(v)
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

if ( size(vertex,2) == 3 )
    numVertex = size(vertex,1);
else
    vertex = vertex';
    if ( size(vertex,2) == 3 )
        numVertex = size(vertex,1);
    else
        error('Vertex coordinates must be 3D!');
    end
end

%--------------------------------------------------------------------------
% Construct triangulation tools
%--------------------------------------------------------------------------

% Triangulation structure
surfTri = triangulation( face, vertex );

% Vertex IDs of boundary vertices
bdyID = unique( surfTri.freeBoundary );

% Face IDs of facets attached to each vertex
vertexFace = surfTri.vertexAttachments;

% The circumcenters of each facet
CoM = surfTri.circumcenter;

%--------------------------------------------------------------------------
% Calculate interior angles of triangle faces
%--------------------------------------------------------------------------

% Face edges
e1 = vertex( face(:,3), : ) - vertex( face(:,2), : );
e2 = vertex( face(:,1), : ) - vertex( face(:,3), : );
e3 = vertex( face(:,2), : ) - vertex( face(:,1), : );

ee = cat( 3, e1, e2, e3 );

% Face edge lengths
l1 = sqrt( sum( e1.^2, 2 ) );
l2 = sqrt( sum( e2.^2, 2 ) );
l3 = sqrt( sum( e3.^2, 2 ) );

z1 = cross( e3, -e2, 2 ); z1 = z1 ./ sqrt( sum( z1.^2, 2 ) );
z2 = cross( e1, -e3, 2 ); z2 = z2 ./ sqrt( sum( z2.^2, 2 ) );
z3 = cross( e2, -e1, 2 ); z3 = z3 ./ sqrt( sum( z3.^2, 2 ) );

theta1 = 2 .* atan2( dot( cross( e3, -e2, 2 ), z1, 2 ), ...
    l3 .* l2 + dot( e3, -e2, 2 ) );
theta2 = 2 .* atan2( dot( cross( e1, -e3, 2 ), z2, 2 ), ...
    l1 .* l3 + dot( e1, -e3, 2 ) );
theta3 = 2 .* atan2( dot( cross( e2, -e1, 2 ), z3, 2 ), ...
    l2 .* l1 + dot( e2, -e1, 2 ) );

theta = [ theta1, theta2, theta3 ];

% Face IDs of obtuse faces
obtuseFace = any( theta > (pi/2), 2 );

%--------------------------------------------------------------------------
% Calculate the Gaussian curvature of each vertex
%--------------------------------------------------------------------------
K = zeros( numVertex, 1 );
Kr = zeros( numVertex, 1 );

% Used to construct mixed Voronoi averaging area of each vertex
EE1_ID = [3,1,2];
EE2_ID = [2,3,1];

for vID = 1:numVertex
    
    angSum = 0;
    vArea = 0;
    vFIDx = vertexFace{vID};
    
    for f = 1:length(vFIDx)
        
        % Find the vID of the current vertex in the attached face ---------
        fID = vFIDx( f );
        curFace = face( fID, : );
        
        [~, loc] = ismember( vID, curFace );
        
        if ( loc == 0 )
            error('Vertex ID not found!');
        end
        
        % Sum the interior face angles around the current vertex ----------
        angSum = angSum + theta( fID, loc );
        
        % Find the mixed Voronoi averaging area for the current vertex ----       
        if obtuseFace(fID)
            v2COM = ( ee(fID,:,loc) ./ 2 ) - vertex(vID,:);
        else
            v2COM = CoM(fID,:) - vertex(vID,:);
        end
        
        EE1 = ee( fID, :, EE1_ID(loc) ) ./ 2;
        EE2 = -ee( fID, :, EE2_ID(loc) ) ./ 2;
        
        A1 = sqrt( sum( cross( EE1, v2COM, 2 ).^2, 2 ) ) ./ 2;
        A2 = sqrt( sum( cross( EE2, v2COM, 2 ).^2, 2 ) ) ./ 2;
        
        vArea = vArea + A1 + A2;
        
    end

    % Combine contributions to the local Gaussian curvature ---------------
    if ( ismember( vID, bdyID ) )
        Kr(vID) = pi - angSum;
    else
        Kr(vID) = 2 .* pi - angSum;
    end
    
    K(vID) = Kr(vID) ./ vArea;
    
end


end

