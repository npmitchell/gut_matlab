function [ patchIm, imref, zeroID, MIP, SIP ] = ...
    texture_patch_to_image( FF, VV, TF, TV, I, Options )
%TEXTURE_PATCH_TO_IMAGE This function generates an image from the texture
%mapping of a mesh triangulation with respect to an texture volume I. One 
%input is the "physical" or "real-space" mesh triangulation. This
%is the mesh triangulation you actually wish to display.  The other input
%is the "virtual" or "texture-space" mesh triangulation.  This input maps
%the physical mesh into texture space.  The virtual mesh must have the same
%number of facets as the physical mesh (each face needs a texture
%mapping!), but may have different connectivities (i.e. the physical mesh
%may be simply connected while the virtual mesh is multiply connected) or a
%different number of vertices (i.e. re-using the texturing from a single
%template face for multiple physical faces).
%
%   texture_patch_to_image( FF, VV, TF, TV, I, Options )
%
%   Input Parameters:
%       - FF:       #Fx3 face connectivity list of the physical mesh
%       - VV:       #Vx2 vertex coordinate list of the physical mesh.
%                   NOTE: VV is provided in (X,Y) format
%       - TF:       #Fx3 face connectivity list of the virtual mesh
%       - TV:       #KxD texture coordinate list. Coordinates can be given
%                   as real pixel positions in the texture image or as a
%                   range [0..1]. D = (2,3). NOTE: TV is provided in
%                   (row, column)/(row, column, page) format
%       - I:        The texture image volume. This function supports:
%                   Grayscale:  [I1 x I2] (2D)
%                   Grayscale:  [I1 x I2 x I3] (3D)
%
%       - Options:  Structure containing the standard options for a
%                   textured surface patch, such as EdgeColor, EdgeAlpha,
%                   etc.  See MATLAB documentation for more information.
%
%       - Options.imSize:       The size of the output image
%       - Options.baseSize:     The side length in pixels of the smallest
%                               side of the output image when using the
%                               tight mesh bounding box
%       - Options.xLim:         The x-bounds of the output image
%       - Options.yLim:         The y-bounds of the output image
%       - Options.pixelSearch:  The method of searching for the faces
%                               containing pixel centers
%                                   - 'AABB' (requires GPToolBox)
%                                   - 'Default' (MATLAB built-ins)
%       - Options.numLayers:    The number of onion layers to create
%                               Format is [ (num +), (num -) ]
%       - Options.layerSpacing: The spacing between adjacent onion layers
%       - Options.smoothIter:   Number of iterations of Laplacian mesh
%                               smoothing to run on the mesh prior to
%                               vertex normal displacement (requires
%                               GPToolBox)
%       - Options.vertexNormal: User supplied vertex unit normals to the
%                               texture triangulation
%       - Options.Interpolant:  A pre-made texture image volume interpolant
%
%   Output Parameters:
%       - patchIm:  The output image stack
%       - imref:    The image spatial referencing object
%       - zeroID:   patchIm(:,:,zeroID) is data layer zero
%       - MIP:      The maximum intensity projection of the output stack
%       - SIP:      The summed intensity projection of the output stack
%
%   by Dillon Cislo 10/25/2019

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

% Verify input triangulations ---------------------------------------------

% Check that the number of faces is consistent
if ( size(FF,1) ~= size(TF,1) )
    error('texture_patch_to_image:inputs', ...
        [ 'Number of real-space faces must match', ...
        'the number of texture-space faces' ] );
end

% Check that the faces are triangles
if ( size(FF,2) ~= 3 ) || ( size(TF,2) ~= 3 )
    error('texture_patch_to_image:inputs', ...
        'All mesh faces must be elements of a triangulation' );
end

% Check the dimensions of the physical vertices
if( size(VV,2) ~= 2 )
    error('texture_patch_to_image:inputs', ...
        'Real-space vertex coordinates must be 2D');
end

% Check the dimensions of the texture vertices/image
if ~isempty(I) % Allow users to supply pre-made interpolant
    if ~( (size(TV,2) == 3) && (ndims(I) == 3) )
        if ~( (size(TV,2) == 2) && ismatrix(I) )
            error('texture_patch_to_image:inputs',...
                'Invalid texture vertex/image input');
        end
    end
end

% Re-scale texture vertices to true pixel positions if necessary
if max(TV(:)) < 2 
    for i = size(TV,2)
        TV(:,i) = (size(I,i)-1) .* TV(:,i) + 1;
    end
end

% Create triangulation objects
realTri = triangulation( FF, VV );
textureTri = triangulation( TF, TV );

% Extract input options ---------------------------------------------------

% Use default options if none are supplied
if nargin < 6, Options = struct(); end

% Extract the output image x-bounds
if isfield( Options, 'xLim' )
    xLim = Options.xLim;
else
    xLim = [ min(VV(:,1)) max(VV(:,1)) ];
end

% Extract the output image y-bounds
if isfield( Options, 'yLim' )
    yLim = Options.yLim;
else
    yLim = [ min(VV(:,2)) max(VV(:,2)) ];
end

% Extract base image size
if isfield( Options, 'baseSize')
    baseSize = Options.baseSize;
else
    baseSize = 512;
end

% Extract the output image size
if isfield( Options, 'imSize' )
    imSize = Options.imSize;
else
    bndRatio = [ yLim(2)-yLim(1), xLim(2)-xLim(1) ];
    bndRatio = bndRatio ./ min(bndRatio);
    imSize = ceil( baseSize .* bndRatio );
end

% Construct the image spatial referencing object
if nargout > 1
    imref = imref2d( imSize, xLim, yLim );
end

% Determine the pixel center search method
if isfield( Options, 'pixelSearch' )
    pixelSearch = Options.pixelSearch;
    % Check that GPToolBox is on the path
    if strcmp(pixelSearch, 'AABB')
        if ~exist('in_element_aabb', 'file')
            warning('GPToolbox was not found. Using default settings.');
            pixelSearch = 'Default';
        end
    end
else
    pixelSearch = 'Default';
end

% Determine if any onion layers are to be produced
makePosLayers = false;
makeNegLayers = false;
makeOnion = false;
if isfield( Options, 'numLayers' )
    numLayers = Options.numLayers;
    if (abs(numLayers(1)) > 0), makePosLayers = true; end
    if (abs(numLayers(2)) > 0), makeNegLayers = true; end
    if ( makePosLayers || makeNegLayers ), makeOnion = true; end
else
    numLayers = [0 0];
end

% Determine the location of the data layer zero
if nargout > 2
    zeroID = numLayers(2) + 1;
end

% Determine the onion layer spacing
if isfield( Options, 'layerSpacing' )
    layerSpacing = Options.layerSpacing;
else
    layerSpacing = 5;
end

% Determine how many iterations of Laplacian mesh smoothing to run on the
% input mesh prior to vertex displacement along normal vectors
if isfield( Options, 'smoothIter' )
    smoothIter = Options.smoothIter;
    % Check that GPToolBox is on the path
    if smoothIter > 0
        if ~exist('laplacian_smooth', 'file')
            warning('GPToolBox was not found. Using default settings.');
            smoothIter = 0;
            smoothMesh = false;
        else
            smoothMesh = true;
        end
    else
        smoothMesh = false;
    end
else
    smoothIter = 0;
    smoothMesh = false;
end

% Process texture vertex normals
if isfield( Options, 'vertexNormal' )
    vN = Options.vertexNormal;
else
    vN = [];
end

% Create the interpolant object from the input image object ---------------
if isfield( Options, 'Interpolant' )
    II = Options.Interpolant;
else
    II = griddedInterpolant(single(I), 'cubic');
end

%--------------------------------------------------------------------------
% CREATE TEXTURE PATCH IMAGE FOR DATA LAYER ZERO
%--------------------------------------------------------------------------

% Calculate the centers of the pixels in real-space ----------------------

% The real-space spacing between adjacent columns
dx = ( xLim(2) - xLim(1) ) / imSize(2);

% The real-space spacing between adjacent rows
dy = ( yLim(2) - yLim(1) ) / imSize(1);

% The coordinates of the centers of pixels in real-space
[X, Y] = meshgrid( (xLim(1)+dx/2):dx:(xLim(2)-dx/2), ...
    (yLim(1)+dy/2):dy:(yLim(2)-dy/2) );

XY = [ X(:), Y(:) ];

% Find the faces containing the specified pixel centers in real-space and
% the barycentric interpolation coordinates of the pixel centers with
% respect to those containing faces ---------------------------------------

if strcmp( pixelSearch, 'AABB' )
    
    % NOTE: Points not located within any face return 0
    [ pixFaces, ~, ~, ~ ] = in_element_aabb( VV, FF, XY, [], [], [] );
    pixFaces( pixFaces == 0 ) = NaN;
    
    pixBaryReal = nan( size(XY,1), 3 );
    pixBaryReal(~isnan(pixFaces), :) = cartesianToBarycentric( realTri, ...
        pixFaces(~isnan(pixFaces)), XY(~isnan(pixFaces), :) );

else
    
    % NOTE: Points not located within any face return NaN
    [ pixFaces, pixBaryReal ] = pointLocation( realTri, XY );
    
end

% Find the texture coordinates of the pixel centers -----------------------
pixTexture = nan( size(XY,1), ndims(I) );
pixTexture(~isnan(pixFaces), :) = barycentricToCartesian( textureTri, ...
    pixFaces(~isnan(pixFaces)), pixBaryReal(~isnan(pixFaces), :) );

% Interpolate to find the texture intensity at the pixel center -----------
patchIm = zeros( size(XY,1), 1 );

if ismatrix(I)
    patchIm(~isnan(pixFaces)) = II(pixTexture(~isnan(pixFaces), 1), ...
        pixTexture(~isnan(pixFaces), 2));
else
    patchIm(~isnan(pixFaces)) = II(pixTexture(~isnan(pixFaces), 1), ...
        pixTexture(~isnan(pixFaces), 2), pixTexture(~isnan(pixFaces), 3 ));
end

% Resize to desired output dimensions
patchIm = reshape( patchIm, imSize );

%--------------------------------------------------------------------------
% GENERATE ONION LAYERS
%--------------------------------------------------------------------------

if makeOnion
    
    % Smooth the texture triangulation if desired
    if smoothMesh
        % Hold the positions of the boundary vertices fixed under smoothing
        % to avoid mesh collapse
        bdyIDx = unique(realTri.freeBoundary);
        smoothV = laplacian_smooth( TV, TF, 'cotan', ...
            bdyIDx, 0.1, 'implicit', TV, smoothIter );
    else
        smoothV = textureTri.Points;
    end
    
    % Calculate the vertex unit normals of the texture triangulation
    if isempty(vN)
        vN = vertexNormal( triangulation( TF, smoothV ) );
    end
    
    %----------------------------------------------------------------------
    % Generate negative layers
    %----------------------------------------------------------------------
    
    if makeNegLayers
        
        % The stack holding the negative layers
        mStack = zeros( [ imSize abs(numLayers(2)) ] );
        
        for i = 1:abs(numLayers(2))
            
            % Determine the texture coordinates for the current layer
            OV = smoothV - i .* layerSpacing .* vN;
            
            % Construct a triangulation representation
            smoothTri = triangulation( TF, OV );
            
            % Find the texture coordinates of the pixel centers
            pixTexture = nan( size(XY,1), ndims(I) );
            pixTexture(~isnan(pixFaces), :) = barycentricToCartesian( ...
                smoothTri, pixFaces(~isnan(pixFaces)), ...
                pixBaryReal(~isnan(pixFaces), :) );
            
            % Interpolate to find the texture intensity at the pixel center
            layerIm = zeros( size(XY,1), 1 );
            
            if ismatrix(I)
                
                layerIm(~isnan(pixFaces)) = ...
                    II( pixTexture(~isnan(pixFaces), 1), ...
                    pixTexture(~isnan(pixFaces), 2) );
                
            else
                
                layerIm(~isnan(pixFaces)) = ...
                    II( pixTexture(~isnan(pixFaces), 1), ...
                    pixTexture(~isnan(pixFaces), 2), ...
                    pixTexture(~isnan(pixFaces), 3) );
                
            end
            
            % Resize to desired output dimensions
            layerIm = reshape( layerIm, imSize );
            
            % Add the current layer image to the stack
            mStack(:,:,i) = layerIm;
            
        end
        
        % Flip stack to reflect proper spatial ordering
        mStack = flip( mStack, 3 );
        
        % Concatenate the negative layers and data layer zero
        patchIm = cat(3, mStack, patchIm);
        
    end
    
    if makePosLayers
        
        % The stack holding the negative layers
        pStack = zeros( [ imSize abs(numLayers(1)) ] );
        
        for i = 1:abs(numLayers(1))
            
            % Determine the texture coordinates for the current layer
            OV = smoothV + i .* layerSpacing .* vN;
            
            % Construct a triangulation representation
            smoothTri = triangulation( TF, OV );
            
            % Find the texture coordinates of the pixel centers
            pixTexture = nan( size(XY,1), ndims(I) );
            pixTexture(~isnan(pixFaces), :) = barycentricToCartesian( ...
                smoothTri, pixFaces(~isnan(pixFaces)), ...
                pixBaryReal(~isnan(pixFaces), :) );
            
            % Interpolate to find the texture intensity at the pixel center
            layerIm = zeros( size(XY,1), 1 );
            
            if ismatrix(I)
                
                layerIm(~isnan(pixFaces)) = ...
                    II( pixTexture(~isnan(pixFaces), 1), ...
                    pixTexture(~isnan(pixFaces), 2) );
                
            else
                
                layerIm(~isnan(pixFaces)) = ...
                    II( pixTexture(~isnan(pixFaces), 1), ...
                    pixTexture(~isnan(pixFaces), 2), ...
                    pixTexture(~isnan(pixFaces), 3) );
                
            end
            
            % Resize to desired output dimensions
            layerIm = reshape( layerIm, imSize );
            
            % Add the current layer image to the stack
            pStack(:,:,i) = layerIm;
            
        end
        
        % Concatenate the negative layers and data layer zero
        patchIm = cat(3, patchIm, pStack);
        
    end
    
end

%--------------------------------------------------------------------------
% FORMAT OUTPUT
%--------------------------------------------------------------------------

% Calculate output stack MIP ----------------------------------------------
if nargin > 3
    MIP = mat2gray(max( patchIm, [], 3 ));
end

% Calculate output stack SIP ----------------------------------------------
if nargin > 4
    SIP = mat2gray(sum( patchIm, 3 ));
end

% Re-scale output image stack ---------------------------------------------
limits = [ min(patchIm(:)) max(patchIm(:)) ];

if limits(2) ~= limits(1)
    delta = 1 / (limits(2) - limits(1));
    patchIm = delta .* ( patchIm - limits(1) );
end

% Make sure all values are on the range [0, 1]
patchIm( patchIm < 0 ) = 0;
patchIm( patchIm > 1 ) = 1;
    

end

