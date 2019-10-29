function patchIm = texture_patch_to_image( FF, VV, TF, TV, I, Options )
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
%       - I:        The texture image volume.This function supports:
%                   Grayscale:  [I1 x I2] (2D)
%                   Grayscale:  [I1 x I2 x I3] (3D)
%
%       - Options:  Structure containing the standard options for a
%                   textured surface patch, such as EdgeColor, EdgeAlpha,
%                   etc.  See MATLAB documentation for more information.
%
%       - Options.imSize:       The size of the output image
%       - Options.xLim:         The x-bounds of the output image
%       - Options.yLim:         The y-bounds of the output image
%
%   Output Parameters:
%       - patchIm:  The output image
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
if ~( (size(TV,2) == 3) && (ndims(I) == 3) )
    if ~( (size(TV,2) == 2) && ismatrix(I) )
        error('texture_patch_to_image:inputs',...
            'Invalid texture vertex/image input');
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

% Create the interpolant object from the input image object ---------------
II = griddedInterpolant(single(I), 'cubic');

% Extract input options ---------------------------------------------------

% Use default options if none are supplied
if nargin < 6, Options = struct(); end

% Extract the output image size
if isfield( Options, 'imSize' )
    imSize = Options.imSize;
else
    imSize = [ 512 512 ];
end

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

%--------------------------------------------------------------------------
% CREATE TEXTURE PATCH IMAGE
%--------------------------------------------------------------------------

% Calculate the centers of the pixels in real-space ----------------------

% The real-space spacing between adjacent columns
dx1 = ( xLim(2) - xLim(1) ) / ( imSize(2) - 1 );
dx2 = ( xLim(2) - xLim(1) - dx1 ) / ( imSize(2) -1 );

% The real-space spacing between adjacent rows
dy1 = ( yLim(2) - yLim(1) ) / ( imSize(1) - 1 );
dy2 = ( yLim(2) - yLim(1) - dy1 ) / ( imSize(1) - 1 );

% The coordinates of the centers of pixels in real-space
[X, Y] = meshgrid( (xLim(1)+dx1/2):dx2:(xLim(2)-dx1/2), ...
    (yLim(1)+dy1/2):dy2:(yLim(2)-dy1/2) );

XY = [ X(:), Y(:) ];

% Find the faces containing the specified pixel centers in real-space -----
% NOTE: Points not located within any face return NaN
pixFaces = pointLocation( realTri, XY );

% Find the barycentric interpolation coordinates of the pixel centers with
% respect to their containing face ----------------------------------------
pixBaryReal = nan( size(XY,1), 3 );
pixBaryReal(~isnan(pixFaces), :) = cartesianToBarycentric( realTri, ...
    pixFaces(~isnan(pixFaces)), XY(~isnan(pixFaces), :) );

% Find the texture coordinates of the pixel centers -----------------------
pixTexture = nan(  size(XY,1), ndims(I) );
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
patchIm = mat2gray(reshape( patchIm, imSize ));

end

