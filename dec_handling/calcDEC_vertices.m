function [div_grid, rot_grid, dec, divv, rotv, bc2d] = ...
    calcDEC_vertices(mesh2d, mesh3d, x0, y0, vx, vy, dt)
% calcDEC(mesh2d, mesh3d, xx, yy, vx, vy)
% 
% Parameters
% ----------
% mesh2d : 
% mesh3d : 
% xx : 
% yy : 
% vx : 
% vy : 
% dt : float (default=1)
%   time interval between t0 and t1 
%
% Returns 
% -------
% div_grid :
% rot_grid :
% dec : 
% divv : 
% rotv : 
% bc2d :

% grab the advected points in 2d
x1 = x0 + vx ;
y1 = y0 + vy ;

% now we have 2d mesh, and 2d vectors

%%
% Perform DEC analysis of these fields
% FF = mesh.f ;  % face connectivity list of the mesh
% VV = mesh.v ;  % vertex positions in 3D of the CLOSED mesh

%for rectilinear grid
% cutmesh = struct() ;
% [cutmesh.f, ~] = defineFacesRectilinearGrid(pixelXY, size(image,1), size(image,2)) ; 
% [xx,yy] = meshgrid(1:size(image,1), 1:size(image,2)) ; % maybe the dimensions are swapped here, not sure
% cutmesh.v = [xx(:), yy(:)];

% % optional(?): Tile the cut mesh (topological disk, so endcaps are cut AND
% % a seam is cut) and then map heads and tails of velocity vectors in PB
% % space into 3d to get vel3D.
% tileCount = [1, 1] ;  % #Tiles in periodic dimension to generate from cut mesh
% [tmF, tmV2d, tmV3d, ~] = tileAnnularCutMesh(cutmesh, tileCount);
% [pt0, fieldfaces0, tr0] = interpolate2Dpts_3Dmesh(tmF, tmV2d, tmV3d, [x0(:), y0(:)]) ;
% % Advect to the next poisitions along 2D velocity field
% x1 = x0(:) + velocity2d(:, 1) ;
% y1 = y0(:) + velocity2d(:, 2) ;

%% Push to 3D
FF = mesh_ss.f ; 
VV = mesh_ss.v ;
assert(all(all(FF == mesh_rect_ss.f)))


% x0 and y0 have the same units as mesh_rect --> the pullback space untis!
% y0 = y0(:) * sdf ;
% x0 = x0(:) * sdf ;
% y1 = y1(:) * sdf ;
% x1 = x1(:) * sdf ;
[pt0, fieldfaces0, ~] = interpolate2Dpts_3Dmesh(FF, mesh_rect_ss.v(:, 1:2), VV, [y0(:), x0(:)]) ;
[pt1, fieldfaces1, ~] = interpolate2Dpts_3Dmesh(FF, mesh_rect_ss.v(:, 1:2), VV, [y1(:), x1(:)]) ;
% Convert displacement into velocity
dt = 1;
vel3D = (pt1 - pt0) / dt ;
% now we have 3d mesh, 3d velocities on some faces fieldfaces0

%% Interpolate vel3D onto facebarycenters so we have one vector for each
% face -- we need one vector per face! To do so, get bc's in 2D pb space.
bc = barycenter(mesh_ss.v, mesh_ss.f) ;
xx = X0(1, :)' ;
yy = Y0(:,1) ;
vel3Dx = reshape(vel3D(:, 1), [size(yy, 1), size(xx, 1)]) ;
vel3Dy = reshape(vel3D(:, 2), [size(yy, 1), size(xx, 1)]) ;
vel3Dz = reshape(vel3D(:, 3), [size(yy, 1), size(xx, 1)]) ;
v3Dxi = griddedInterpolant(Y0 * sdf, X0 * sdf, vel3Dx) ;
v3Dyi = griddedInterpolant(Y0 * sdf, X0 * sdf, vel3Dy) ;
v3Dzi = griddedInterpolant(Y0 * sdf, X0 * sdf, vel3Dz) ;

% % check orientation
% trisurf(triangulation(FF, mesh_rect.v), 'EdgeColor', 'none', 'FaceVertexCData', [1:length(FF)*0.5, 1:length(FF)*0.5]', 'facealpha', 0.2)
% hold on;
% scatter3(Y0(:) * sdf, X0(:)*sdf, 0*X0(:), 20, piv_temp.VY(:), 'markeredgealpha', 0.4)
% quiver(y0(:), x0(:), piv_temp.VY(:), piv_temp.VX(:), 1) 
% view(2)
% %waitfor(gcf)
% scatter(y0,x0, 200, vel3Dz(:), 'filled')
% maxPBxy = max(mesh_rect.v) ;
% maxPBx = maxPBxy(1) ;  % this is max of z in embedding 
% maxPBy = maxPBxy(2) ; % this is the circumference 
% % scale factor is the ratio of circumferential direction in mesh to PB
% % image
% % sdf2 = maxPBy / size(im1, 1);

%%
% Initialize a face-based velocity field in 3D, with one vector per face

bc2d = barycenter(mesh_rect_ss.v, mesh_rect_ss.f) ;
vel3Df = zeros(size(FF, 1), 3) ;
vel3Df(:, 1) = v3Dxi(bc2d(:,1), bc2d(:, 2)) ;
vel3Df(:, 2) = v3Dyi(bc2d(:,1), bc2d(:, 2)) ;
vel3Df(:, 3) = v3Dzi(bc2d(:,1), bc2d(:, 2)) ;

% check the result
% clf
% hold on;
% scatter(bc2d(:, 1), bc2d(:, 2), 20, vel3Df(:, 3), 'markeredgealpha', 0.05)
% scatter(y0, x0, 30, vel3Dz(:), 'filled')

%% Construct DEC class instance
% You have a fibration , where every face in 2D has a 3D vector associated
% with it. DEC treats a 3D mesh with the SAME topology so that the same mesh
% faces are now treated in embedding with the 3d vectors in hand.

% make assertion that the faces of the 2d mesh which you just used to get
% the face-based 3d vectors are the same as the faces of the 3d mesh which
% you feed DEC. 

assert(all(all(mesh_rect_ss.f == mesh_ss.f)))

DEC = DiscreteExteriorCalculus(mesh_ss.f, mesh_ss.v) ;

% Now resolve the vector field for decomposition
[~, v0t, ~, ~, ~, ~, ~] = ...
    resolveTangentNormalVelocities(mesh_rect_ss.f, mesh_ss.v, vel3Df, 1:size(mesh_ss.f, 1), mesh_rect_ss.v) ;

% Compute divergence and rotational flow
%divv = DEC.divergence(v0t);
rotv = DEC.curl(v0t);
divv = DEC.div(v0t);

%% visualize div and curl
% opts = struct() ;
% opts.visible = 'on' ;
% opts.clims = {2*std(divv(:)), std(rotv(:))} ;
% opts.labels = {'$\nabla \cdot v$', '$\star \textrm{d} (v_t^\flat)$'} ;
% nFieldsOnSurface(mesh, {divv, rotv}, opts)

%%
% figure(1)
% opts_2 = struct();
% opts_2.visible = 'on';
% opts_2.clims = {std(rotv(:))} ;
% opts_2.labels = {'$\star \textrm{d} (v_t^\flat)$'} ;
% nFieldsOnSurface(mesh_rect_ss, {rotv}, opts_2)
% caxis([-0.1, 0.1])
% 
% pause(.1)
% 
% DEC_2d_Curl_mov(t) = getframe(gcf);

%% interpolate the div/curl results into the 2d mesh's vertices 

%2d face barycenters
bc2d = barycenter(mesh_rect_ss.v, mesh_rect_ss.f);

curl_interp = scatteredInterpolant(bc2d(:,1), bc2d(:,2), rotv);
div_interp = scatteredInterpolant(bc2d(:,1), bc2d(:,2), divv);

rot_grid = curl_interp(Y0*sdf, X0*sdf);
div_grid = div_interp(Y0*sdf, X0*sdf);

