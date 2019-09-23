function texture_patch_3d( FF, VV, TF, TV, IV, Options)
%TEXTURE_PATCH_3D This function will display a triangulated mesh, like the
%MATLAB function "patch", with a texture mapping given by the 3D image
%volume I. One set of inputs to
%this function is the "physical" or "real-space" mesh triangulation. This
%is the mesh triangulation you actually wish to display.  The other input
%is the "virtual" or "texture-space" mesh triangulation.  This input maps
%the physical mesh into texture space.  The virtual mesh must have the same
%number of facets as the physical mesh (each face needs a texture
%mapping!), but may have different connectivities (i.e. the physical mesh
%may be simply connected while the virtual mesh is multiply connected) or a
%different number of vertices (i.e. re-using the texturing from a single
%template face for multiple physical faces).
%
%   texture_patch_3d( FF, VV, TF, TV, IV, Options );
%
%   Input Parameters:
%       - FF:       #Fx3 face connectivity list of the physical mesh
%       - VV:       #Vx3 vertex coordinate list of the physical mesh. NOTE:
%                   VV is provided in (X,Y,Z) format!
%       - TF:       #Fx3 face connectivity list of the virtual mesh
%       - TV:       #Kx3  texture coordinate list. Coordinates can be given
%                   as real pixel positions in the texture image or as a
%                   range [0..1]. NOTE: TV is provided in
%                   (row, column, page) format!
%       - IV:        The 3D texture image volume.This function supports:
%                   Grayscale:  [I1 x I2 x I3]
%       - Options:  Structure containing the standard options for a
%                   textured surface patch, such as EdgeColor, EdgeAlpha,
%                   etc.  See MATLAB documentation for more information.
%
%       - Options.PSize:    Special option, defines the image texture size
%                           for each individual polygon.  A lower number
%                           gives a more pixelated texture. Defaults to 64.
%       - Options.Rotation: Special option, defines the rotation matrix of 
%                           all surface patches.
%       - Options.Translation: Special option, defines the translation 
%                           vector applied after rotation on all surfaces.
%       - Options.Rotation: Special option, defines the dilation factor 
%                           applied to all surfaces after rotation and 
%                           translation. 
%
%
%   by Dillon Cislo 08/14/2019
%   NPMitchell added Rotation, Translation, & Dilation options 09/2019
%   NPMitchell grouped surfaces into a Parent container for speedup 09/2019

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

% Patch FaceColor MUST be a texture.
Options.FaceColor = 'texturemap';

% Size of the 2D texture image used for each triangle
if isfield( Options, 'PSize' )
    sizep = round(Options.PSize(1));
    Options = rmfield(Options, 'PSize');
else
    sizep = 64;
end

% Also get rotation for all patches
if isfield( Options, 'Rotation' )
    rot = Options.Rotation ;
    Options = rmfield(Options, 'Rotation');
    rotate = true ;
else
    rotate = false ;
end

% Also get translation for all patches
if isfield( Options, 'Translation' )
    trans = Options.Translation ;
    Options = rmfield(Options, 'Translation');
    translate = true ;
else
    translate = false ;
end

% Also get dilation for all patches
if isfield( Options, 'Dilation' )
    dilation = Options.Dilation ;
    Options = rmfield(Options, 'Dilation');
    dilate = true ;
else
    dilate = false ;
end

% Check that the number of faces is consistent
if ( size(FF,2) ~= size(TF,2) )
    error('texture_patch_3d:inputs', ...
        [ 'Number of real-space faces must match', ...
        'the number of texture-space faces' ] );
end

% Check the dimensions of the physical vertices
if ( size(VV,2) ~= 3 )
    if( size(VV,2) == 2 )
        VV = [ VV zeros( size(VV,1), 1 ) ];
    else
        error('texture_patch_3d:inputs', ...
            'Invalid real-space vertex input');
    end
end

% Check the dimensions of the texture vertices
if ( size(TV,2) ~= 3 )
    error('texture_patch_3d:inputs', 'Invalid texture vertex input');
end

% Check the dimensions of the texure image
if ( ndims(IV) ~= 3 )
    error('texture_patch_3d:inputs', 'Invalid texture image input');
end

% Re-scale texture vertices to true pixel positions if necessary
if ( max(TV(:)) < 2 )
    
    TV = [ (size(IV,1)-1) .* TV(:,1) + 1, ...
        (size(IV,2)-1) .* TV(:,2) + 1, ...
        (size(IV,3)-1) .* TV(:,2) + 1 ];
    
end

% Create the interpolant object from the input image object
IVI = griddedInterpolant(single(IV), 'cubic');

%--------------------------------------------------------------------------
% FIND PATCH INTERPOLATION VALUES
%--------------------------------------------------------------------------

% The patch used for every triangle
% For simplicity these will be the same for each triangle regardless of
% size in either the physical or texture space. A smarter algorithm would
% probably take these features into account
J = zeros( (sizep+1), (sizep+1), class(IV) );

% Linear indices of the 2D image associated with the triangle patch
jind = (sizep+1)^2:-1:1;

% Barycentric interpolation values for query points in the triangle patch
[ lambda1, lambda2, lambda3 ] = ...
    calculateBarycentricInterpolationValues( sizep );

%--------------------------------------------------------------------------
% CREATE THE SURFACE PLOT
%--------------------------------------------------------------------------

hold on

% Create container for all surface objects 
container = hggroup(gca) ;

% Loop through all of the faces of the mesh
for i = 1:size(FF,1)
    
    % The current face's physical vertex coordinates
    V = VV( FF(i,:), : );
    
    % The current face's texture vertex coordinates
    tV = TV( TF(i,:), : );
    
     % Define the physical face in 'surface' format
    x=[V(1,1) V(2,1); V(3,1) V(3,1)];
    y=[V(1,2) V(2,2); V(3,2) V(3,2)];
    z=[V(1,3) V(2,3); V(3,3) V(3,3)];
    
    % Define the texture coordinates of the 'surface'
    xyz = [ tV(1,:); tV(2,:); tV(3,:); tV(3,:) ];
    
    % Calculate the texture interpolation coordinatex ---------------------
    
    pos(:,1) = xyz(1,1)*lambda1 + xyz(2,1)*lambda2 + xyz(3,1)*lambda3;
    pos(:,2) = xyz(1,2)*lambda1 + xyz(2,2)*lambda2 + xyz(3,2)*lambda3;
    pos(:,3) = xyz(1,3)*lambda1 + xyz(2,3)*lambda2 + xyz(3,3)*lambda3;
    
    % Map texture to surface image ----------------------------------------
    J(jind) = IVI( pos(:,1), pos(:,2), pos(:,3) ) ;
    
    % Show the surface ----------------------------------------------------
    if rotate
        xyz = [x(1) y(1) z(1); x(3) y(3) z(3); x(2) y(2) z(2)] ;
        xyzp = (rot * xyz')' ;
        x = [xyzp(1, 1) xyzp(2, 1); xyzp(3, 1) xyzp(3, 1)];
        y = [xyzp(1, 2) xyzp(2, 2); xyzp(3, 2) xyzp(3, 2)];
        z = [xyzp(1, 3) xyzp(2, 3); xyzp(3, 3) xyzp(3, 3)];
    end
    if translate
        x = x + trans(1) ;
        y = y + trans(2) ;
        z = z + trans(3) ;
    end
    if dilate
        x = x * dilation ;
        y = y * dilation ;
        z = z * dilation ;
    end
    surface( container, x, y, z, J, Options );
    
end

% t1 = hgtransform('Parent',gca);
% Txy = makehgtform('translate',Options.translate);
% set(t1,'Matrix',Txy)

hold off

end


%==========================================================================
% CALCULATEBARYCENTRICINTERPOLATIONVALUES
%==========================================================================

function [ lambda1, lambda2, lambda3 ] = ...
    calculateBarycentricInterpolationValues( sizep )
% CALCULATEBARYCENTRICINTERPOLATIONVALUES Calculates the barycentric
% interpolation values of a set of query points in a square of size
% (sizep+1)X(sizep+1) relative to an abstract triangle in the upper region
% of that square

% Define the triangle in the upperpart of the square
x1 = sizep; y1 = sizep;
x2 = sizep; y2 = 0;
x3 = 0; y3 = sizep;

% Calculate the barycentric coordinates for each query point on the square
% (These will correspond to the linear indices corresponding to the 2D
% image associated with the triangle)
[x,y] = ndgrid(0:sizep,0:sizep); % NOTE THE USE OF NDGRID
x = x(:); y = y(:);

detT = (x1-x3)*(y2-y3) - (x2-x3)*(y1-y3);
lambda1 = ( (y2-y3).*(x-x3) + (x3-x2).*(y-y3) ) ./ detT;
lambda2 = ( (y3-y1).*(x-x3) + (x1-x3).*(y-y3) ) ./ detT;
lambda3 = 1 - lambda1 - lambda2;

end

