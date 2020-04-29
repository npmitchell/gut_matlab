function tesMesh = generateFPSSTessellation(F, V, varargin)
%GENERATEFPSSTESSELLATION Generates a Voronoi tessellation of a mesh
%triangulation using a 'Farthest Point Sampling Set' as the seed points
%
%   INPUT ARGUMENTS:
%
%       - F:        #Fx3 connectivity list
%       - V:        #VxD vertex coordinate list
%
%   (Name, Value)-pair Options:
%
%       - 'NumPoints':      The number of seed points (100)
%       - 'Overlap':        The overlap distance between adjacent cells
%                           If a number [0 1) is provided it is considered
%                           a percentage of the largest dimension of each
%                           cell.  A number >= 1 is a fixed distance for
%                           all cells. (0.1)
%       - 'SeedIDx':        A set of vertex indices used as a seed point
%                           for the farthest point sampling set (Default
%                           is a single randomly chosen vertex ID)
%       - 'ClipEars':       If true, ears will be clipped from the
%                           submeshes
%       - 'VertexNormals':  #VxD array of vertex unit normals
%       - 'Display':        Set to true for verbose output
%       - 'MinFaces':       The minimum number of allowed faces in a
%                           submesh
%       - 'MinVertices':    The minimum number of allowed vertices in a
%                           submesh
%
%   OUTPUT ARGUMENTS:
%
%       - tesMesh:          A collection of Vornoi cells stored as a
%                           (NumPoints)x1 struct array with the following
%                           fields
%
%                           - 'f': The face connectivity list of each
%                           submesh
%                           - 'v': The 3D vertex coordinate list of each
%                           submesh
%                           - 'vn': The 3D vertex normal list of each
%                           submesh
%                           - 'u': The 2D vertex coordinate list of each
%                           submesh
%                           - 'Center': The face ID in the full mesh AND
%                           the submesh of the face that contains the point
%                           (x,y) = (0,0) in the domain of parameterization
%                           - 'VinFullMesh': The vertex ID in the full mesh
%                           of the corresponding point in the submesh
%                           - 'FinFullMesh': The face ID in the full mesh
%                           of the corresponding face in the submesh
%
% by Dillon Cislo 02/04/2020

%==========================================================================
% INPUT PROCESSING
%==========================================================================

% Validate mandatory inputs -----------------------------------------------
validateattributes(F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'positive', 'finite'});
validateattributes(V, {'numeric'},  {'2d', 'finite', 'nonnan'});

numVertex = size(V,1);

% The IDs of vertices on the mesh boundary
% bdyIDx = unique(T.freeBoundary);

% Process optional inputs -------------------------------------------------
numPoints = 100;
overlap = 0.1;
seedIDx = randi( numVertex, 1, 1 );
clipEars = true;
VN = [];
dispType = true;
minFaces = 6;
minVertices = 7;

for i = 1:length(varargin)
    if isa(varargin{i},'double') 
        continue;
    end
    if isa(varargin{i},'logical')
        continue;
    end
    if ~isempty(regexp(varargin{i},'^[Nn]um[Pp]oints','match'))
        numPoints = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Oo]verlap','match'))
        overlap = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Ss]eed[Ii][Dd][Xx]','match'))
        seedIDx = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Cc]lip[Ee]ars','match'))
        clipEars = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Vv]ertex[Nn]ormals','match'))
        VN = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Dd]isplay','match'))
        dispType = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Mm]in[Ff]aces','match'))
        minFaces = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Mm]in[Vv]ertices','match'))
        minVertices = varargin{i+1};
    end
end

validateattributes(numPoints, {'numeric'}, ...
    {'scalar', 'positive', 'integer', '<=', numVertex});
validateattributes(overlap, {'numeric'}, ...
    {'scalar', 'positive', 'finite', 'nonnan'});
validateattributes(seedIDx, {'numeric'}, ...
    {'vector', 'positive', 'integer', '<=', numVertex});
validateattributes(minFaces, {'numeric'}, ...
    {'scalar', 'positive', 'finite', 'nonnan'});
validateattributes(minVertices, {'numeric'}, ...
    {'scalar', 'positive', 'finite', 'nonnan'});

assert(numel(seedIDx) <= numPoints, ...
    'Number of supplied seeds exceeds the number of desired submeshes!');

if isempty(VN)    
    VN = per_vertex_normals(V, F, 'Weighting', 'angle');   
else   
    validateattributes(VN, {'numeric'}, {'2d', 'finite', 'nonnan'});
    assert(isequal(size(VN), size(V)), ...
    'Input vertex normal list is improperly sized!');
end

%==========================================================================
% DETERMINE REMAINING SAMPLE POINTS
%==========================================================================
if dispType, fprintf('Finding seed points: '); end

% The vertex ID of each sample point
sampleIDx = zeros(numPoints, 1);
sampleIDx(1:numel(seedIDx)) = seedIDx;

% Will hold the distances of all mesh vertices to the final sample set
% D(i,j) is the distance of the ith mesh vertex to the jth seed point
D = inf(numVertex, numPoints);

for i = 1:numel(seedIDx)
    D(:,i) = perform_fast_marching_mesh(V, F, seedIDx(i));
end

% The minimum distance of each vertex in the mesh to any of the extant
% sample points
Dmin = min(D, [], 2);

if dispType, progressbar(i, numPoints); end
while (i < numPoints)
    
    i = i+1;
    
    if dispType, progressbar(i, numPoints); end
    
    [~, newID] = max(Dmin);
    assert(~ismember(newID, sampleIDx), ...
        'Update vertex is already a sample point!');
    sampleIDx(i) = newID;
    
    % Update the minimum distances
    D(:,i) = perform_fast_marching_mesh(V, F, newID);
    Dmin = min(Dmin, D(:,i));
    
end

%==========================================================================
% GENERATE A VORONOI TESSELLATION OF THE 3D MESH
%==========================================================================
if dispType, fprintf('Subdividing input mesh: '); end

% Will hold the vertex IDs in the original mesh associated 
% to each Voronoi cell
submVIDx = cell(numPoints, 1);

% Will hold the face IDs in the original mesh associated
% to each Voronoi cell
submFIDx = cell(numPoints, 1);

% The ID of the seed point nearest to each mesh vertex
[~, nearestSeed] = min(D, [], 2);

% The distance from each vertex to its tentative seed point
seedVtxDist = V(nearestSeed, :) - V;
seedVtxDist = sqrt( sum( seedVtxDist.^2, 2 ) );

% The maximum distance of each seed point to any of its cell's member
% vertices
maxDist = zeros(numPoints, 1);
for i = 1:numPoints
    maxDist(i) = max( seedVtxDist( nearestSeed == i ) );
end

% Calculate the overlap distance based on the maximum distance
if overlap < 1
    overlapDist = overlap .* maxDist;
else
    overlapDist = overlap;
end

% Determine the member vertices of each Voronoi cell given the cell's
% overlap distance
for i = 1:numPoints
    
    if dispType, progressbar(i, numPoints); end
    
    % The distances from mesh vertices to seed points including overlap
    DD = D;
    DD(:, setdiff(1:numPoints, i)) = DD(:, setdiff(1:numPoints, i)) + ...
        overlapDist(i);
    
    % The member vertices of the current cell
    [~, VinSubm] = min(DD, [], 2);
    VinSubm = find(VinSubm == i);
    
    submVIDx{i} = VinSubm;
    
    % Find the faces in the whole mesh comprised only of vertices in the
    % current submesh
    FinSubm = find( all( ismember( F, VinSubm ), 2 ) );
    submFIDx{i} = FinSubm;
    
    assert( numel(VinSubm) >= minVertices, ...
        sprintf(['Number of vertices in submesh %d is less than the ' ...
        'the minimum. Please choose a different seed point, ' newline ...
        'increase the overlap distance, decrease the number of ' ...
        'desired points, ' newline 'or refine the input mesh ' ...
        'and try again.'], i) );
    
    assert( numel(FinSubm) >= minFaces, ...
        sprintf(['Number of faces in submesh %d is less than the ' ...
        'the minimum. Please choose a different seed point, ' newline ...
        'increase the overlap distance, decrease the number of ' ...
        'desired points, ' newline 'or refine the input mesh ' ...
        'and try again.'], i) );

end

%==========================================================================
% MAP EACH VORONOI CELL TO THE PLANE
%==========================================================================
if dispType, fprintf('Mapping submeshes to the plane: '); end

% For each Voronoi submesh:
submF = cell(numPoints, 1); % The face connectivity list
submV = cell(numPoints, 1); % The 3D vertex coordinate list
submU = cell(numPoints, 1); % The 2D vertex coordinate list
submVN = cell(numPoints, 1); % The 3D vertex normal list

% The ID of the face containing the point (x,y) = (0,0) in the pullback of
% the submesh to the unit disk
submCenter = cell(numPoints, 1);

for i = 1:numPoints
    
    if dispType, progressbar(i, numPoints); end
    
    %----------------------------------------------------------------------
    % Clean/Process Submesh Triangle Soup Arrays
    %----------------------------------------------------------------------
    
    vIDx_OM = submVIDx{i};
    fIDx_OM = submFIDx{i};
    
    % The 3D coordinates of the submesh vertices
    VCell = V(vIDx_OM, :);
    
    % The faces in the original mesh contained in the current submesh
    FCell_OM = F(fIDx_OM, :);
    
    % The updated face list with the local vertex IDs
    FCell = FCell_OM;
    for k = 1:numel(vIDx_OM)
        FCell(FCell == vIDx_OM(k)) = k;
    end
    
    if clipEars
        
        [ newFCell, ~, newVCell, ~ ] = clipMesh( FCell, VCell, true );
        
        if ( (size(newFCell, 1) < minFaces) || ...
                (size(newVCell, 1) < minVertices) )
            
            % Only remove unreferenced vertices
            [ newFCell, ~, newVCell, ~ ] = clipMesh( FCell, VCell, false );
            
        end
        
        
    else
        
        % Only remove unreferenced vertices
        [ newFCell, ~, newVCell, ~ ] = clipMesh( FCell, VCell, false );
        
    end
    
    assert( size(newVCell, 1) >= minVertices, ...
        sprintf(['Number of vertices in submesh %d is less than the ' ...
        'the minimum. Please choose a different seed point, ' newline ...
        'increase the overlap distance, decrease the number of ' ...
        'desired points, ' newline 'or refine the input mesh ' ...
        'and try again.'], i) );
    
    assert( size(newFCell, 1) >= minFaces, ...
        sprintf(['Number of faces in submesh %d is less than the ' ...
        'the minimum. Please choose a different seed point, ' newline ...
        'increase the overlap distance, decrease the number of ' ...
        'desired points, ' newline 'or refine the input mesh ' ...
        'and try again.'], i) );
    
    % Determine which submesh vertices were removed
    [ rmVIDx, ~ ] = ismember( VCell, newVCell, 'rows' );
    rmVIDx = find(~rmVIDx);
    
    % Determine which faces were removed
    rmFIDx = any( ismember( FCell, rmVIDx ), 2 );
    
    % Update current submesh triangle soup
    VCell = newVCell;
    FCell = newFCell;
    
    % Update whole mesh/submesh vertex correspondence
    vIDx_OM( rmVIDx ) = [];
    submVIDx{i} = vIDx_OM;
    
    % Update whole mesh/submesh face correspondence
    fIDx_OM( rmFIDx ) = [];
    submFIDx{i} = fIDx_OM;
    
    %----------------------------------------------------------------------
    % Map the Submesh to the Unit Disk
    %----------------------------------------------------------------------
    [UCell, FParam] = surface_parameterization( FCell, VCell, struct() );
    UCell = 2 .* ( UCell - 0.5 );
    
    assert( isequal(FParam, FCell) && size(UCell, 1) == size(VCell, 1), ...
        [ 'Mesh topology for submesh %d was changed during ' ...
        'paremeterization.' ], i );
    
    % Correct for possible reflections ------------------------------------
    e12 = UCell( FCell(:,2), : ) - UCell( FCell(:,1), : );
    e13 = UCell( FCell(:,3), : ) - UCell( FCell(:,1), : );
    
    e12 = [ e12, zeros( size(e12,1), 1 ) ];
    e13 = [ e13, zeros( size(e13,1), 1 ) ];
    
    pm = cross( e12, e13, 2 );
    pm = sign( pm(3) );
    
    if numel(unique(pm)) > 1     
        warning( 'Faces of submesh %d are inconsistently ordered!\n', i );     
    elseif unique(pm) < 0      
        UCell = [ UCell(:,1), -UCell(:,2) ];      
    end
    
    % Clip the boundary points of the pullback to the unit disk -----------
    z = complex(UCell(:,1), UCell(:,2));
    outIDx = abs(z) > 1;
    z( outIDx ) = exp( 1i .* angle( z( outIDx ) ) );
    
    UCell = [ real(z), imag(z) ];

    % Minimize area distortion in the mapping -----------------------------
    a = minimizeIsoarealMobiusEnergy( FCell, VCell, UCell );
    
    % Make sure that the point being mapped to the origin is not a vertex
    if any( ismember( UCell, a, 'rows' ) )
        a = a + 10 * eps;
    end
    
    % Convert to complex coordinates
    a = complex(a(1), a(2)); 

    % Apply Mobius mapping
    w = (z - a) ./ ( 1 - conj(a) .* z );
    
    UCellIso = [ real(w), imag(w) ];
    
    % Check for self-intersections ----------------------------------------
    % If the Mobius optimized mesh contains self-intersections, we simply
    % return the non-optimized mesh instead
    
    [ intersects, ~ ] = mesh_self_intersection_3d( FCell, ...
        [ UCellIso, zeros(size(UCellIso,1), 1) ] );
    
    if intersects
        
        % Check for intersections on the original pullback mesh
        [ intersects, ~ ] = mesh_self_intersection_3d( FCell, ...
            [ UCell, zeros(size(UCell,1), 1) ] );
        
        if intersects
            warning( [ 'The pullback for mesh %d contains ' ...
                'self-intersections!\n' ], i ); 
        end
        
    else
        
        UCell = UCellIso;
        
    end
    
    % Determine which triangle contains the point (0,0) -------------------
    centerFace = pointLocation( triangulation( FCell, UCell ), 0, 0 );

    %----------------------------------------------------------------------
    % Update Output Arrays
    %----------------------------------------------------------------------
    submF{i} = FCell;
    submV{i} = VCell;
    submU{i} = UCell;
    submVN{i} = VN(vIDx_OM, :);
    submCenter{i} = [ fIDx_OM(centerFace) centerFace ];

end

%==========================================================================
% ASSEMBLE TESSELATED MESH
%==========================================================================

tesMesh = struct( 'f', submF, 'v', submV, 'vn', submVN, ...
    'u', submU, 'Center', submCenter, ...
    'VinFullMesh', submVIDx, 'FinFullMesh', submFIDx );

end



