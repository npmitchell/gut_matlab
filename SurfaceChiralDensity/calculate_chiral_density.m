function chiral = calculate_chiral_density( face, vertex )
%CALCULATE_CHIRAL_DENSITY Calculates the surface chiral density
%pseudotensor of a mesh triangulation
%   Input Parameters:
%       - face:         The connectivity list of the triangulation
%
%       - vertex:       The vertex coordinate list
%
%   Output Parameters:
%       - chiral:       The chiral density pseudotensor defined on each
%                       face
%
% By Dillon Cislo 02/2019
% Edits by Noah Mitchell 02/2019

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

if ( numel(size(face)) ~= 2 )
    error( 'Connectivity list is improperly sized!' );
end

if ( numel(size(vertex)) ~= 2 )
    error( 'Vertex coordinate list is improperly sized!' );
end

if ( size(face,2) == 3 )
    numFaces = size(face,1);
else
    face = face';
    if ( size(face,2) == 3 )
        numFaces = size(face,1);
    else
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
% GEOMETRY PROCESSING
%--------------------------------------------------------------------------

surfTri = triangulation( face, vertex );

% Calculate Basic Quantities ---------------------------------------------- 

% The face edges
e1 = vertex(face(:,3),:) - vertex(face(:,2),:);
e2 = vertex(face(:,1),:) - vertex(face(:,3),:);

% The edge lengths
l1 = sqrt( sum( e1.^2, 2 ) );
l2 = sqrt( sum( e2.^2, 2 ) );

% The face unit normal vector
n = cross( e1, e2, 2 );
n = n ./ sqrt( sum( n.^2, 2 ) );

% Caluclate Mid-Edge Unit Normal Vectors ----------------------------------

% Vertex IDs defining each edge
eIDx = surfTri.edges;

% Edge attachment list
edgeFace = surfTri.edgeAttachments( eIDx );
% If there are boundaries, these edges have only one element in edgeFace.
% Fill in their 'second' face in edgeFace to be another copy of the same
% face with which they already are paired.
resizeCell = @(x) repmat( x, 1, 1+mod(numel(x),2) );
try
    edgeFace = cell2mat( cellfun( resizeCell, edgeFace, ...
        'UniformOutput', false ) );
catch
    edgeFace = cellfun( resizeCell, edgeFace, ...
        'UniformOutput', false );
    
    % Build up warning message
    badindices = find(cellfun('length', edgeFace) ~= 2) ;
    nfaces = [];
    for jj=1:length(badindices)
        bindex = badindices(jj);
        nfaces = [nfaces; length(edgeFace{bindex})] ;
        disp(['bad faces = ', num2str(edgeFace{bindex})])
    end
    msg = ['Some elements were not included: # ignored =  ', ...
            num2str(nnz(cellfun('length', edgeFace) ~= 2))];
    msg2 = ['The ignored edges share ', mat2str(nfaces), ' faces'];
    disp(msg)
    disp(msg2)

    % Inspect the offending faces
    trisurf(face, vertex(:,1), vertex(:,2), vertex(:,3));
    badinds = cellfun('length', edgeFace) ~= 2 ;
    badfaces = edgeFace(badinds) ;
    % highlight the vertices of the bad faces
    badf = [] ;
    for jj=1:length(badfaces)
        badf = [badf; badfaces{jj}];
    end
    badf = unique(badf(:)) ;
    badpts = face(badf, :) ;

    % Check which ones are problematic
    edgeFace(cellfun('length', edgeFace) ~= 2)
    % Filter these rows from eIDx
    eIDx = eIDx(cellfun('length', edgeFace) == 2, :);
    % Filter these rows from edgeFace
    edgeFace = edgeFace(cellfun('length', edgeFace) == 2) ;
    edgeFace = cell2mat(edgeFace);

end

% The mid-edge unit normal
nEdge = ( n( edgeFace(:,1), : ) + n( edgeFace(:,2), : ) ) ./ 2;
nEdge = nEdge ./ sqrt( sum( nEdge.^2, 2 ) );

% Re-structure Mid-Edge Normals to Lie on Facet Edges ---------------------
e1IDx = sort( [ face(:,3), face(:,2) ], 2 );
e2IDx = sort( [ face(:,1), face(:,3) ], 2 );
e3IDx = sort( [ face(:,2), face(:,1) ], 2 );

eIDx = sort( eIDx, 2 );

[~, e1IDx] = ismember( e1IDx, eIDx, 'rows' );
[~, e2IDx] = ismember( e2IDx, eIDx, 'rows' );
[~, e3IDx] = ismember( e3IDx, eIDx, 'rows' );

% The mid-edge unit normal corresponding to each edge in each face
n1 = nEdge( e1IDx, : );
n2 = nEdge( e2IDx, : );
n3 = nEdge( e3IDx, : );

%--------------------------------------------------------------------------
% CONSTRUCT CHIRAL DENSITY PSEUDOTENSOR
%--------------------------------------------------------------------------
chiral = cell( numFaces, 1 );

for fID = 1:numFaces
    
    % The unit-normal of the current face ---------------------------------
    fN = n( fID, : )';
    
    % The tangent space projection operator -------------------------------
    Pp = eye(3) - fN * fN';
    
    % Construct the coordinate change operator ----------------------------
    
    li = l1(fID); % Length of first edge
    lj = l2(fID); % Length of second edge
    
    ei = e1(fID,:)' ./ li; % Unit vector along first edge
    ej = e2(fID,:)' ./ lj; % Unit vector along second edge
    
    CC = [ cross( ej, fN )'; -cross( ei, fN )' ];
    CC = CC ./ dot( fN, cross( ei, ej ) );
    
    % Construct normal variation operator ---------------------------------
    ni = n1(fID,:)';
    nj = n2(fID,:)';
    nk = n3(fID,:)';
    
    NV = 2 .* [ ( nj - nk ) ./ li, ( nk - ni ) ./ lj ];   
    
    % Construct matrix cross-product operator -----------------------------
    cxN = [ 0 -fN(3) fN(2); fN(3) 0 -fN(1); -fN(2) fN(1) 0 ];
    
    % The chiral density tensor -------------------------------------------
    chiral{fID} = cxN * NV * CC * Pp;
    
   
end


end

