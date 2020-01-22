function [pts, fieldfaces, tr0] = interpolate2Dpts_3Dmesh(faces, v2d, v3d, uv)
%INTERPOLATE2DPTS_3DMESH Map points in 2D to 3D using mesh faces  
%   Map points in 2D space living in 2D representation of a mesh to 3D
%   space living on the 3D representation of the same mesh, using
%   barycentric coordinates.
%
% Parameters
% ----------
% faces : M x 3 int array
%   The mesh connectivity list, indexing into the vertex array(s)
% v2d : P x 2 float array
%   The mesh vertex locations in 2d
% v3d : P x 3 float array
%   The mesh vertex locations in 3d
% uv : N x 2 float array
%   The points to map to 3D using barycentric coordinates
%
% Returns
% -------
% pts : N x 3 float array
%   the 3d point coordinates, living on 3D-representation of the mesh
% fieldfaces : N x 1 int array
%   indices into faces of the faces on which uv live. These are where 
%   the supplied field uv is defined on the mesh, and where the output
%   pts live on the 3d mesh.
%
% NPMitchell 2019

tr0 = triangulation(faces, v2d) ;
[fieldfaces, baryc0] = pointLocation(tr0, uv) ; 

% Handle case where there are NaNs
bad = find(isnan(fieldfaces)) ;
baryc0(isnan(fieldfaces), :) = 0 ; 
fieldfaces(isnan(fieldfaces)) = 1 ;  

% Interpolate the position in 3D given relative position within 2D
% triangle.
% x123(i) is the x coords of the elements of triangle t_contain(i)
vxa = v3d(:, 1) ;
vya = v3d(:, 2) ;
vza = v3d(:, 3) ;

assert(size(vxa, 1) == size(v2d, 1))

% Map to the faces
tria = tr0.ConnectivityList(fieldfaces, :) ;

% Modulo the vertex IDs: trisa are the triangle vertex IDs
% trisa = mod(tria, size(vxa, 1)) ;
% trisa(trisa == 0) = size(vxa, 1) ;
% x123a = vxa(tria) ;
% y123a = vya(tria) ;
% z123a = vza(tria) ;

% Multiply the vertex positions by relative weights.
% Note that baryc gives weights of the three vertices of triangle
% t_contain(i) for pointLocation x0(i), y0(i)
pts = [sum(baryc0 .* vxa(tria), 2), ...
    sum(baryc0 .* vya(tria), 2), ...
    sum(baryc0 .* vza(tria), 2) ] ;

% Handle case where there are NaNs
pts(bad, :) = NaN  ;
fieldfaces(bad) = NaN ;


%%%
% Note: an alternative method is to interpolate
% Xai = scatteredInterpolant(tm0X, tm0Y, tm0v3d(:, 1)) ;
% Yai = scatteredInterpolant(tm0X, tm0Y, tm0v3d(:, 2)) ;
% Zai = scatteredInterpolant(tm0X, tm0Y, tm0v3d(:, 3)) ;
% pt0 = [Xai(x0(:), y0(:)), Yai(x0(:), y0(:)), Zai(x0(:), y0(:))] ;
% disp('Finding barycentric coordinates')
% % Compute barycentric coordinates for later
% tr0 = triangulation(tm0f, [tm0X, tm0Y]) ;
% [fieldfaces, baryc0] = pointLocation(tr0, [x0(:), y0(:)]) ;
%%%

end

