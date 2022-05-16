function [U, Uall, nnIDx] = laplacian_boundary_smooth( F, V, nnOrder, ...
    L_method, fixIDx, lambda, method, max_iter, tol )
% LAPLACIAN_BOUNDARY_SMOOTH Attempts to smooth the boundaries of a
% triangular mesh via a modified Laplacian smoothing method designed to
% mitigate boundary shrinkage
%
%	INPUT PARAMETERS:
%
%       - F:            #Fx3 face connectivity list
%
%       - V:            #Vx3 input vertex coordinate list
%
%       - nnOrder:      The size of the natural neighborhood around the
%                       boundary to consider for smoothing
%
%       - L_method:     Method for basic Laplacian construction
%                       ('uniform', 'cotan')
%
%       - fixIDx:       #FVx1 list of fixed vertex indices
%
%       - lambda:       diffusion speed parameter
%
%       - method:       Method for smoothing
%                       ('implicitx', 'explicit')
%
%       - max_iter:     The maximum number of smoothing iterations
%
%       - tol:          The threshold for insufficient change used for
%                       smoothing termination
%
%   OUTPUT PARAMETERS:
%
%       - U:            #Vx3 output vertex coordinate list
%
%       - Uall:         #Vx3xiters list of vertex coordinate lists for each
%                       smoothing iteration
%
%   by Dillon Cislo 12/09/2021


%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

if (~exist('F', 'var'))
    error('Please supply face connectivity list');
end

if (~exist('V', 'var'))
    error('Please supply vertex coordinate list');
end

validateattributes(V, {'numeric'}, {'2d', 'real', 'finite', 'nonnan'});
validateattributes(F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'real', 'positive', 'integer', '<=', size(V,1)});

TR = triangulation(F,V);
E = edges(TR);

% #Ex2 array of fIDs of the faces attached to a particular edge.
% If an edge is a border edge (i.e., only attached to a single face), then
% that fID is listed twice for dimensional consistency
resizeCell = @(x) repmat( x, 1, 1+mod(numel(x),2) );
edgeFace = edgeAttachments( TR, E );
edgeFace = cell2mat( cellfun( resizeCell, edgeFace, ...
    'UniformOutput', false ) );

% The average edge length in the mesh
h = avgedge(V,F);

if (~exist('nnOrder', 'var')), nnOrder = 3; end
if (~exist('L_method','var')), L_method = 'uniform'; end
if (~exist('fixIDx','var')), fixIDx = []; end
if (~exist('lambda','var')), lambda = 0.1; end
if (~exist('method','var')), method = 'implicit'; end
if (~exist('tol','var')), tol = 0.001; end
if (~exist('max_iter','var')), max_iter = 10; end

% Compute/validate the natural neighborhood of the boundary components ----

% The boundary components of the mesh
allBdyIDx = compute_boundaries(F);

% Generate an adjacency matrix from mesh edges
A = sparse( [E(:,1); E(:,2)], [E(:,2); E(:,1)], 1, size(V,1), size(V,1) );

% Create a graph object from this adjacency matrix
G = graph(A);

% Find distances between pairs of vertices
D = distances(G, 'Method', 'unweighted');

allNN = cell(numel(allBdyIDx), 1);
for i = 1:numel(allBdyIDx)
    
    % Assemble the current boundary neighborhood
    curBdyIDx = allBdyIDx{i};
    curNN = [];
    for j = 1:numel(curBdyIDx)
        
        curID = curBdyIDx(j);
        curNN = [curNN; find( (D(:,curID) <= nnOrder) )];
        % curNN = [curNN; ...
        %     find( (D(:,curID) <= nnOrder) & (D(:,curID) > 0) )];
            
    end
    
    % Ensure that boundary neighborhoods are non-intersecting
    curNN = unique(curNN(:));
    for j = 1:(i-1)
        assert(~any(ismember(curNN, allNN{j})), ...
            'Intersecting boundary neighborhoods detected');
    end
    
    allNN{i} = curNN;

end

bdyIDx = [allBdyIDx{:}].';
nnIDx = vertcat(allNN{:});

assert(all(ismember(bdyIDx, nnIDx)), ...
    'Some boundary vertices are not properly included');

% Fix all non-included vertices
allIDx = (1:size(V,1)).';
fixIDx = [fixIDx; allIDx(~ismember(allIDx, nnIDx))];

%--------------------------------------------------------------------------
% Perform Boundary Relaxation
%--------------------------------------------------------------------------
        
% Only construct the uniform Laplacian once
if strcmpi(L_method, 'uniform')
    L = uniform_boundary_laplacian(E, bdyIDx, nnIDx);
end

% Create a copy of the input vertex coordinate list
S = V;

% Build sparse identity matrix
I = speye(size(V,1), size(V,1));

% Place for factorization and symmtery flag used by min_quad_with_fixed
P = [];
% sym = [];

iter = 0;
U = S;
U_prev = S;
if nargout >= 2
    Uall = [];
end

while( iter < max_iter && (iter == 0 || max(abs(U(:)-U_prev(:)))>tol*h))
    
    U_prev = U; % Store a copy of the previous iteration
    
    % Update geometric Laplacian operator
    if strcmpi(L_method, 'cotan')
        L = cotan_boundary_laplacian(F, U, bdyIDx, nnIDx);
    end
    
    switch method
        
        case 'implicit'
            Q = (I-lambda*L);
            % could prefactor Q for 'uniform' case
            for d = 1:size(S,2)
                [U(:,d),P] = ...
                    min_quad_with_fixed( Q*0.5, -U(:,d), fixIDx, ...
                    S(fixIDx,d), [], [], P);
            end
            
        case 'explicit'
            Q = (I+lambda*L);
            U = Q * U;
            % enforce boundary
            U(fixIDx,:) = S(fixIDx,:);
            
        otherwise
            error(['' method ' is not a supported smoothing method']);
            
    end
    
    if (nargout >= 2)
        Uall = cat(3, Uall, U);
    end
    
    iter = iter + 1;
    
end
    
end

function L = uniform_boundary_laplacian(E, bdyIDx, nnIDx)

I = [E(:,1); E(:,2)];
J = [E(:,2); E(:,1)];

rmIDx = ismember(I, bdyIDx) & (~ismember(J, bdyIDx));
% rmIDx = rmIDx | ~ismember(I, nnIDx);

I(rmIDx) = [];
J(rmIDx) = [];

L = sparse(I, J, 1, max(E(:)), max(E(:)));
% Lnorm = sum(L,2);
L = L - diag(sum(L,2));
% L = repmat(1./Lnorm, 1, size(L,2)) .* L;

end

function L = cotan_boundary_laplacian(F, V, bdyIDx, nnIDx)

C = cotangent(V, F); C = [ C(:); C(:) ];

I = F(:,[2 3 1]);
J = F(:,[3 1 2]);

tmp = I(:);
I = [ I(:); J(:) ];
J = [ J(:); tmp(:) ];

rmIDx = ismember(I, bdyIDx) & (~ismember(J, bdyIDx));
% rmIDx = rmIDx | ~ismember(I, nnIDx);

I(rmIDx) = [];
J(rmIDx) = [];
C(rmIDx) = [];

L = sparse(I, J, C, size(V,1), size(V,1));
% Lnorm = sum(L,2);
L = L - diag(sum(L,2));
% L = repmat(1./Lnorm, 1, size(L,2)) .* L;

end


function boundaries = compute_boundaries(F)
%COMPUTE_BOUNDARIES Finds the boundary vertices of a potientially multiply
%connected mesh. Re-creates the functionality of 'compute_boundaries.m' in
%the ImSAnE package.
%
%   INPUT PARAMETERS:
%
%       - F:            #FxP polygonal face connectivity list (P=3 for a
%                       triangulation).  Faces do not need to have
%                       consistent orientation (even though they should and
%                       everyone's life would be better if they did)
%
%   OUTPUT PARAMETERS:
%
%       - boundaries:   1x#B cell array of boundary components.
%                       boundaries{i} is a 1x#BV row vector of the vertex
%                       IDs in the ith boundary componenent
%
% by Dillon Cislo 11/18/2019

% Build the unordered vertex adjacency list
A = sparse( F, F(:, [2:end 1]), 1 );
A = A+A';

% Find vertices in the edges that occur only once
[ bdyIDx, ~ ] = find(A == 1);
bdyIDx = unique(bdyIDx);

numBdy = 0; % The number of boundary components
boundaries = {}; % The output parameter
while ~isempty(bdyIDx)
    
    numBdy = numBdy + 1;
    
    % Choose the starting point of the current boundary component
    i = bdyIDx(1);
    
    % If a boundary vertex is connected to more than 2 other boundary
    % vertices, something is wrong (i.e. a triangle sharing only one
    % vertex with the rest of the mesh)
    u = find(A(i,:) == 1);
    if (numel(u) ~= 2), warning('Problem in boundary'); end
    
    boundary = [i u(1)];
    
    % Iterate over the rest of the boundary
    s = boundary(2);
    i = 2;
    while i <= size(A,2)
        
        u = find(A(s,:) == 1);
        if (numel(u) ~= 2), warning('Problem in boundary'); end
        
        if u(1) == boundary(i-1)
            s = u(2);
        else
            s = u(1);
        end
        
        if s ~= boundary(1)
            boundary = [ boundary s ];
        else
            break;
        end
        
        i = i+1;
       
    end
    
    bdyIDx = setdiff( bdyIDx, boundary );
    boundaries{numBdy} = boundary;
    
end

end
