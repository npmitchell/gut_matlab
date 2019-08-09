function laplaceBeltrami = construct_laplace_beltrami( face, vertex )
%CONSTRUCT_LAPLACE_BELTRAMI Constructs the discrete Laplace-Beltrami
%operator on a mesh triangulation as a sparse matrix.  Here we use the
%cotangent discretization:
%
%   lB * f(vi) = sum( ( cot(aij) + cot(bij) ) * ( fj - fi ) ) / ( 2 A(vi) )
%
%Where vi denotes a vertex, f(vi) = fi, fj are scalar functions defined on
%vertices, A(vi) is the mixed Voronoi averaging region of the vertex and
%aij and bij are the angles opposite the edge defined by vertices vi and
%vj.
%
%   Input Parameters:
%       - face:             The triangulation connectivity list
%       - vertex:           The vertex coordinate list
%
%   Ouput Parameter:
%       - laplaceBeltrami:  A sparse matrix operator which acts on vectors
%                           with the same length as the vertex list
%
%   by Dillon Cislo 02/26/2019

%--------------------------------------------------------------------------
% INPUT PROCESSING
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

if ( ( size(vertex,2) == 3 ) || ( size(vertex,2) == 2 ) )
    numVertex = size(vertex,1);
else
    vertex = vertex';
    if ( ( size(vertex,2) == 3 ) || ( size(vertex,2) == 2 ) )
        numVertex = size(vertex,1);
    else
        error('Vertex coordinates must be 2D or 3D!');
    end
end

%--------------------------------------------------------------------------
% Construct triangulation tools
%--------------------------------------------------------------------------

% Triangulation structure
surfTri = triangulation( face, vertex );

% Vertex IDs defining triangulation edges
eIDx = surfTri.edges;

% Edge IDs of boundary edges
bdyID = surfTri.freeBoundary;

if ~isempty(bdyID)
    bdyID = ismember(sort(eIDx,2),sort(bdyID,2),'rows');
else
    bdyID = false( size(eIDx,1), 1 );
end

% Face IDs of facets attached to each vertex
vertexFace = surfTri.vertexAttachments;

% Face IDx of the facets attached to each edge
edgeFace = surfTri.edgeAttachments( eIDx );

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
% Re-structure interior angles to be stored on edges
%--------------------------------------------------------------------------
thetaEdge = zeros( size( eIDx, 1 ), 2 );

for eID = 1:size( eIDx, 1 )
    
    eFIDx = edgeFace{eID};
    
    for i = 1:length(eFIDx) % Either 2 or 1 (boundary edge)
        
        % Find out which  edge in the face the current edge corresponds to
        fID = eFIDx( i );
        loc = find( ~ismember( face(fID,:), eIDx(eID,:) ) );
        
        thetaEdge( eID, i ) = theta( fID, loc );
          
    end
     
end

%--------------------------------------------------------------------------
% Construct mixed Voronoi vertex averaging regions
%--------------------------------------------------------------------------
vA = zeros( numVertex, 1 );

EE1_ID = [3,1,2];
EE2_ID = [2,3,1];

for vID = 1:numVertex
    
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
    
    vA(vID) = vArea;
    
end

%--------------------------------------------------------------------------
% Construct Laplace-Beltrami Operator
%--------------------------------------------------------------------------

% Pre-allocate storage
II = zeros( 2 * size( eIDx, 1 ), 1 );
JJ = zeros( 2 * size( eIDx, 1 ), 1 );
VV = zeros( 2 * size( eIDx, 1 ), 1 );

index = 0;
for vID = 1:numVertex
    
    % The edge IDs of all edges attached to the current vertex
    vEIDx = find( any( eIDx == vID, 2 ) );
    
    for i = 1:length( vEIDx )
        
        % Find the vID of the other vertex in the edge --------------------
        eID = vEIDx( i );
        curEdge = eIDx( eID, : );
        
        if curEdge(1) == vID
            v2ID = curEdge(2);
        else
            v2ID = curEdge(1);
        end
        
        % Construct the Laplace-Beltrami weight ---------------------------
        if bdyID(eID)
            w = 0; % cot( thetaEdge(eID,1) );
        else
            w = cot( thetaEdge(eID,1) ) + cot( thetaEdge(eID,2) );
        end
        
        w = w ./ ( 2 .* vA(vID) );
        
        % Add entries to sparse matrix construction vectors ---------------
        index = index + 1;
        II(index) = vID;
        JJ(index) = vID;
        VV(index) = -w;
        
        index = index + 1;
        II(index) = vID;
        JJ(index) = v2ID;
        VV(index) = w;
        
    end
    
end

% The Laplace-Beltrami Operator
laplaceBeltrami = sparse( II, JJ, VV, numVertex, numVertex );

end

