function [curves, curves2d] = meshConstantCoordSurfaceCurves( F, V2D, V3D, ...
    constCoord, nCurves, coordLim )
%MESHCONSTANTCOORDSURFACECURVES Given a 3D mesh triangulation with a
%corresponding Cartesian 2D parameterization this function calculates
%face-spanning surface curves along the 3D mesh that correspond to lines
%of constant x- or y- in the 2D parameterization.
%
%   Input Parameters:
%       - F:            #Fx3 face connectivity list
%       - V2D:          #Vx2 2D vertex coordinate list
%       - V3D:          #Vx3 3D vertex coordinate list
%       - constCoord:   The coordinate to hold constant along the curves
%                       ('X' or 'Y')
%       - nCurves:      The number of curves to calculate
%       - coordLim:     The range over which to space the curves
%
%   Output Parameters:
%       - curves:       nCurves x 1 cell array.  curves{i} is a list of the
%                       the 3D coordinates where the surface curve
%                       intersects mesh edges/vertices
%
%   by Dillon Cislo 11/13/2019

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

% Validate mandatory triangulation input parameters
validateattributes( F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'positive'} );
validateattributes( V2D, {'numeric'}, ...
    {'2d', 'ncols', 2, 'finite', 'real', 'nonnan'} );
validateattributes( V3D, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'real', 'nonnan'} );

% Generate a MATLAB-triangulation representation of the input
TR = triangulation( F, V2D );

% By default, the vertical coordinate will be held constant
if (nargin < 4), constCoord = 'Y'; end

% Re-format constant coordinate input
if any(strcmp(constCoord, {'x', 'X'}))
    ccID = 1; % Column index of constant coordinate
    ncID = 2; % Column index of non-constant coordinate
elseif any(strcmp(constCoord, {'y', 'Y'}))
    ccID = 2;
    ncID = 1;
else
    error('Invalid user-supplied constant coordinate');
end

% Default settings for remaining optional parameters
if (nargin < 5), nCurves = 100; end
if (nargin < 6), coordLim = [ min(V2D(:,ccID)) max(V2D(:,ccID)) ]; end

% Calculate the locations of the output curves along the constant
% coordinate
curveCoords = coordLim(1):diff(coordLim)/(nCurves-1):coordLim(2);

%--------------------------------------------------------------------------
% CREATE SORTED EDGE LIST STRUCTURES
%--------------------------------------------------------------------------

% Vertex IDs defining edges in the triangulation
eIDx = TR.edges;

% The constant coordinates of the edge vertices
CC = [ V2D(eIDx(:,1), ccID) V2D(eIDx(:,2), ccID) ];

% The non-constant coordinates of the edge vertices
NC = [ V2D(eIDx(:,1), ncID) V2D(eIDx(:,2), ncID) ];

% Sort the constant coordinate array so that the minimum constant
% coordinate value of each edge is in the first column
[ CC, I ] = sort( CC, 2 );

% Determine the indexing that will appropriately map the remaining arrays
% to match the constant coordinate array
I = I.';
I( I == 2 ) = (size(I,2)+1):(2*size(I,2));
I( I == 1 ) = 1:size(I,2);
I = I.';

eIDx = eIDx(I); % Update the edge list
NC = NC(I); % Udpate the non-constant coordinate list

%==========================================================================
% NOTE: For whatever ridiculous reason, MATLAB does not have a built-in
% binary search function.  Some Googling reveals that writing one yourself
% does not even result in a noticable speed-up over using 'find' on an
% unsorted vector of less than ~10^6 elements due to the fact 
% that 'find' is most implemented in pre-compiled C-code...
%
% So until I feel like writing a MEX-function binary
% search there is no reason to actually sort any of these arrays.  I will
% try to make obvious any line in subsequent code blocks that should be
% replaced when a superior search function is implemented
%==========================================================================

% A sorted list of the minimum constant coordinates in each edge
% [ CCmin, Imin ] = sort( CC(:,1) );

% A sorted list of the maximum constant coordinates in each edge
% [ CCmax, Imax ] = sort( CC(:,2) );

%--------------------------------------------------------------------------
% DETERMINE EDGE COLLISIONS WITH CONSTANT COORDINATE CURVES
%--------------------------------------------------------------------------

curves = cell( nCurves, 1 );
if nargout > 1
    curves2d = cell( nCurves, 1) ;
end

for i = 1:nCurves
    
    % The constant coordinate value of the current curve
    cs = curveCoords(i);
    
    % Find all edges intersecting this line (REPLACE FOR BINARY SEARCH)
    collisions = find( ( CC(:,1) <= cs ) & ( CC(:,2) >= cs ) );
    
    % The minimum constant coordinate value of each intersecting edge
    c1 = CC(collisions, 1);
    
    % The maximum constant coordinate value of each interesecting edge
    c2 = CC(collisions, 2);
    
    % The interpolation constant corresponding to the intersection
    t = ( cs - c1 ) ./ ( c2 - c1 );
    
    % Calculate the non-constant coordinate value at each collision
    ncs = t .* NC(collisions, 2) + (1-t) .* NC(collisions, 1);
    
    % Sort the intersections by their non-constant coordinate value
    [ ~, Incs ] = sort( ncs );
    
    % The sorted list of interpolation constants
    t = t(Incs);
    
    % The sorted list of intersecting edge IDs
    collisions = collisions(Incs);
    
    % Interpolate to find the coordinates of each collision in the 3D mesh
    curves{i} = t .* V3D( eIDx(collisions, 2), : ) + ...
        (1-t) .* V3D( eIDx(collisions, 1), : );
    
    % If the 2D positions are requested, supply those here
    if nargout > 1
        curves2d{i} = t .* V2D( eIDx(collisions, 2), : ) + ...
                      (1-t) .* V2D( eIDx(collisions, 1), : );
    end
    
end

end

