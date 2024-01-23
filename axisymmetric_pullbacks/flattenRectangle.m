function newCutMesh = flattenRectangle( cutMesh, LType, leftCornerID, rightCornerID )
% [ THIS IS NOT FINISHED  -- npmitchell 2020]
%FLATTENANNULUS This function homotopically flattens a 3D topological
%annulus by embedding it within a 2D annular Euclidean orbifold
%   INPUT PARAMETERS:
%       - cutMesh:      A struct defining the cut 3D annulus.
%                       See 'cylinderCutMesh.m'.
%
%       - LType:        The type of mesh Laplacian operator used to
%                       construct the linear system solved to find the
%                       embedding coordinates.  This method currently
%                       supports the 'Dirichlet' Laplacian and the 'Mean
%                       Value Coordinates' operator
%   
%   OUTPUT PARAMETERS:
%       - newCutMesh:   An updated struct with the additional field
%                       'newCutMesh.u' which contains the (u,v)-coordinates
%                       of the mesh vertices in the 2D orbifold
%
% by Dillon Cislo & Noah Mitchell 2019-2020 [ THIS IS NOT FINISHED]

%==========================================================================
% THE GEOMETRY OF THE OUTPUT:
%
% The field 'cutMesh.pathPairs' is a (#CV)x2 array that holds the
% correspondences between the original vertex IDs of the cut vertices and
% the vertex IDs of their duplicates.  
%
% Its columns are given by: [ (S1->E1), (S2->E2) ]
%
% Where (S1,S2) are the IDs of the left edge path
% and (E1,E2) are the IDs of the right edge path.
% The final output of the embedding procedure will have the following 
% geometry
%
%           (4)
%    (S1)--------(E1)
%     |            |
%     |            |
% (1) |            | (3)
%     |            |
%     |            |
%    (S2)--------(E2)
%           (2)
%
% Where segments (2)+(4) are identified to create the topological annulus.
% The template for the output region in the (u,v) plane is the domain:
% [0 1] X [0 1]
%
% FINDING THE EMBEDDING:
%
% The embedding coordinates are found as the solution to the linear system
%
%   L u = 0 subject to A u = b
%
% Where L is some mesh Laplacian operator and u is a vector corresponding
% to the (u,v) coordinates of each vertex.  (A,b) is a linear set of
% equality constraints that specify the boundary conditions of the problem.
% For mor information see 'Orbifold Tutte Embeddings' by Noam Aigerman and
% Yaron Lipman.
%
%==========================================================================

if nargin < 2
    LType = 'Dirichlet';
end

S1 = leftCornerID(1);
S2 = leftCornerID(2);
E1 = rightCornerID(1);
E2 = rightCornerID(2);

%==========================================================================
% Treating the cut mesh as a topologcal disk, extract the boundary of the
% cut mesh and partition into segments corresponding to the sides of the
% annular orbifold
%==========================================================================

% Extract an ordered list of all boundary vertex IDs ----------------------
bdyIDx = freeBoundary( triangulation( cutMesh.f, cutMesh.v ) );
assert( isequal( bdyIDx(:,1), [ bdyIDx(end,2); bdyIDx(1:end-1,2) ] ) );
bdyIDx = bdyIDx(:,1);

% Shift the list so that the original start point is in the first position
ind = find( bdyIDx == S1 );
bdyIDx = bdyIDx( [ ind:end, 1:ind-1 ] );

% Ensure that the list is consistently ordered (see diagram above) --------
if find( bdyIDx == E1 ) < find( bdyIDx == S2 )
    bdyIDx = [ bdyIDx(1); flipud( bdyIDx(2:end) ) ];
end

% Find the location of the corner vertices in the updated list ------------
ind = [ find( bdyIDx == S1 ); find( bdyIDx == S2 );
    find( bdyIDx == E2 ); find( bdyIDx == E1 ) ];

% Partition the boundary into segments ------------------------------------
bdyPaths = cell(4,1);
for i = 1:length(ind)-1
    bdyPaths{i} = bdyIDx( ind(i):ind(i+1) );
end
bdyPaths{end} = [ bdyIDx(ind(end):end); bdyIDx(1) ];

%==========================================================================
% Construct the linear constraint sub-system
%==========================================================================

A = sparse( 0, 2 * size( cutMesh.v, 1 ) );
b = [];

%--------------------------------------------------------------------------
% NOTE: I chose here to set hard positional constraints to force the true
% boundary vertices to lie on the lines supporting the vertical edges of
% the square at fixed positions.  Originally, I tried to allow the
% vertices to slide along these lines (see the 'freesquare' functionality
% of Noam Aigerman's 'euclidean_orbifold' package).  For particular meshes
% this led to issues with self-inersections in the embedding.  Perhaps a
% future version can figure out the pathology of these problems and improve
% this method
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Constrain boundary segment (1) to lie on the line (0,1)->(0,0)
%--------------------------------------------------------------------------

% The vertex IDs along segment (1)
segIDx = bdyPaths{1};

v1 = [0,1]; % The start coordinates of the target line
v2 = [0,0]; % The final coordinates of the target line

% Calculate distances between consecutive vertices
d = sqrt( sum( (cutMesh.v(segIDx(1:end-1),:) ...
    - cutMesh.v(segIDx(2:end),:)).^2, 2 ) );

% Convert these to a fractional distance along the lenght of the line
d = [ 0; cumsum(d)/sum(d) ];

% Interpolate to find the (u,v)-coordinates of each vertex
interpXY = [ v1(1) .* (1-d) + v2(1) .* d,  v1(2) .* (1-d) + v2(2) .* d ];

% Set constraints
for i = 1:length(segIDx)
    A(end+1, 2*segIDx(i)-1) = 1; b = [ b; interpXY(i,1) ];
    A(end+1, 2*segIDx(i)) = 1;   b = [ b; interpXY(i,2) ];
end

%--------------------------------------------------------------------------
% Constrain boundary segment (3) to lie on the line (1,0)->(1,1)
%--------------------------------------------------------------------------

% The vertex IDs along segment (3)
segIDx = bdyPaths{3};

v1 = [1,0]; % The start coordinates of the target line
v2 = [1,1]; % The final coordinates of the target line

% Calculate distances between consecutive vertices
d = sqrt( sum( (cutMesh.v(segIDx(1:end-1),:) ...
    - cutMesh.v(segIDx(2:end),:)).^2, 2 ) );

% Convert these to a fractional distance along the lenght of the line
d = [ 0; cumsum(d)/sum(d) ];

% Interpolate to find the (u,v)-coordinates of each vertex
interpXY = [ v1(1) .* (1-d) + v2(1) .* d,  v1(2) .* (1-d) + v2(2) .* d ];

% Set constraints
for i = 1:length(segIDx)
    A(end+1, 2*segIDx(i)-1) = 1; b = [ b; interpXY(i,1) ];
    A(end+1, 2*segIDx(i)) = 1;   b = [ b; interpXY(i,2) ];
end

%--------------------------------------------------------------------------
% Set constraints to identify corresponding vertices in segments (2)+(4)
%--------------------------------------------------------------------------

bdy2 = bdyPaths{2}; bdy2 = bdy2(2:end-1);
bdy4 = bdyPaths{4}; bdy4 = flipud(bdy4(2:end-1));

for i = 1:length(bdy2)
    
    % The u-components are identical
    A(end+1, 2*[ bdy2(i) bdy4(i) ]-1) = [1,-1];
    b = [ b; 0 ];
    
    % The v-components are identical up to a fixed translation
    A(end+1, 2*[ bdy2(i) bdy4(i) ]) = [1,-1];
    b = [ b; -1 ];
    
end

%==========================================================================
% Construct the mesh Laplacian operator
%==========================================================================

if strcmp( LType, 'Dirichlet' ) == 1
    
    L = cotmatrix( cutMesh.v, cutMesh.f );
    
    m = min( min( tril( L, -1 ) ) );
    if m < 0
        
        warning('Mesh is NOT Delaunay!');
        clamp = 1e-2;
        L(L<0) = clamp;
        
        inds = sub2ind( size(L), 1:length(L), 1:length(L) );
        L(inds) = 0;
        L(inds) = -sum(L);
        
    end
    
elseif strcmp( LType, 'MVC' )
    
    
    L = mean_value_laplacian( cutMesh.v, cutMesh.f );
    
else
    
    error('Invalid Laplacian Type Supplied!');
    
end

%--------------------------------------------------------------------------
% Duplicate the Laplacian to work on both planar coordinates
%--------------------------------------------------------------------------

RealL = sparse( size(L,1)*2, size(L,2)*2 );
RealL(1:2:end, 1:2:end) = L;
RealL(2:2:end, 2:2:end) = L;

L = RealL;

%==========================================================================
% Solve the Full Constrained System
%==========================================================================

n_vars = size(A,2);
n_eq = size(A,1);

M = [ L A'; A sparse( n_eq, n_eq ) ];
rhs = [ zeros( n_vars, 1 ); b ];

u_lambda = M \ rhs;

e = max( abs( M*u_lambda -rhs ) );
if e > 1e-6
    error('Linear system solution failed!');
end

u = u_lambda(1:n_vars);
u = [ u(1:2:end), u(2:2:end) ];

%==========================================================================
% Format Output Results
%==========================================================================

newCutMesh = cutMesh;
newCutMesh.u = u;

end

