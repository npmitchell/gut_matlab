function [ TF, TV2D, TV3D ] = tileAnnularCutMesh( cutMesh, tileCount )
%TILEANNULARCUTMESH This function vertically tiles the orbifold pullback of
%an annular cutMesh and returns the parameters of a single triangulation
%   INPUT PARAMETERS:
%       - cutMesh:          A struct defining the cut 3D annulus.
%                           See 'cylinderCutMesh.m'.
%
%       - tileCount:        The vertical tiling parameters.
%                           tileCount(1) tiles above the basic tile.
%                           tileCount(2) tiles below the basic tile.
%
%   OUTPUT PARAMETERS:
%       - TF:               #Fx3 face connectivity list of the combined
%                           triangulation.
%       - TV2D:             #Vx2 2D pullback coordinate list of the
%                           combined triangulation.
%       - TV3D:             #Vx3 3D embedding coordinate list of the
%                           combined triangulation.
%
% by Dillon Cislo

%==========================================================================
% THE GEOMETRY OF THE CUT MESH:
%
% The field 'cutMesh.pathPairs' is a (#CV)x2 array that holds the
% correspondences between the original vertex IDs of the cut vertices and
% the vertex IDs of their duplicates.  
%
% Its columns are given by: [ (S1->E1), (S2->E2) ]
%
% Where (S1,E1) are the IDs of the original start and end points of the cut
% path and (S2,E2) are the IDs of the duplicated start and end points. The
% final output of the embedding procedure will have the following geometry
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
%==========================================================================

% Default tiling creates three stacked tiles
if nargin < 2
    tileCount = [1 1];
end

% Verify input cut mesh
if ~isfield( cutMesh, 'u' )
    error('Cut mesh must contain 2D pullback coordinates');
end

pathPairs = cutMesh.pathPairs;

% The vertex IDs defining the current top seam of the combined tile
topSeamIDx = pathPairs(:,1);

% Find the location of the basic bottom seam vertices in the basic face
% connectivity list as a #(BF)x3 logical array
bottomSeamLoc = ismember( cutMesh.f, pathPairs(:,2) );

% Combined triangulation is just the basic tile to start
TF = cutMesh.f;
TV2D = cutMesh.u;
TV3D = cutMesh.v;

% Find the vertical shift between tiles (should just be 1)
shift = cutMesh.u( pathPairs(1,1), 2 ) - cutMesh.u( pathPairs(1,2), 2 );

% Due to the structure of the cutMesh generation process it is easiest to
% add all new tiles to the top of the basic tile and then shift to reflect
% the desired tiles below the basic tile
for i = 1:sum(abs(tileCount))
    
    % The parameters of the basic tile
    face = cutMesh.f(:);
    V2D = cutMesh.u;
    V3D = cutMesh.v;
    
    % Sew a shifted basic tile to the top of the current combined
    % triangulation -------------------------------------------------------
    
    % The shifted pullback vertices
    V2D(:,2) = V2D(:,2) + i * shift;
    
    % Remove the bottom seam from the vertex lists
    V2D( pathPairs(:,2), : ) = [];
    V3D( pathPairs(:,2), : ) = [];
    
    % Update the basic tile face list to reflect the fact that vertices
    % will be added at the end of the combined vertex coordinate list
    face( ~bottomSeamLoc(:) ) = face( ~bottomSeamLoc(:) ) + size(TV2D,1);
    
    % Reshape the face connectivity list
    face = reshape( face, size( cutMesh.f ) );
    
    % Replace basic bottom seam IDs with the current top seam IDs
    for j = 1:length(topSeamIDx)
        face( ismember(face, pathPairs(j,2) ) ) = topSeamIDx(j);
    end
    
    % Update the current top seam vertex IDs
    topSeamIDx = pathPairs(:,1) + size(TV2D,1);
    
    % Update combined lists
    TF = [ TF; face ];
    TV2D = [ TV2D; V2D ];
    TV3D = [ TV3D; V3D ];
    
end

% Shift the coordinates of the combined triangulation to reflect the
% desired number of tiles below the basic tile
TV2D(:,2) = TV2D(:,2) - abs(tileCount(2)) * shift;

end

