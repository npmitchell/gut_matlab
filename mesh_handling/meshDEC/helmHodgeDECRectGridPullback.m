function [divs, rots, harms, glueMesh] = helmHodgeDECRectGridPullback(cutM, facevels)
% helmHodgeDECRectGridPullback(cutM)
% 
%
% Parameters
% ----------
%
% Returns
% -------
% 
%
% NPMitchell 2020

% Unpack the cutMesh
FF = cutM.f ;
V2D = cutM.u ;
v3drs = cutM.v ;

% % compute COM for each triangle in 2d and 3d --> note that we do this
% % before gluing so that the 2d case is not messed up, and the 3d case is
% % the same either way
% bc = cat( 3, v3drs(FF(:,1), :), v3drs(FF(:,2), :), v3drs(FF(:,3), :) );
% bc = mean( bc, 3 ) ;
% bc2d = cat( 3, V2D(FF(:,1), :), V2D(FF(:,2), :), V2D(FF(:,3), :) );
% bc2d = mean( bc2d, 3 ) ;

% num faces in each row, col is nfU, nfV
% faces are arranged as (nU-1)*(nV-1) * 2
% [~, faceIDgrid] = defineFacesRectilinearGrid(V2D, nU, nV) ;

% % Check that this makes sense
% trisurf(FF, v3drs(:, 1), v3drs(:, 2), v3drs(:, 3), divv, ...
%     'FaceAlpha', 0.2, 'EdgeColor', 'none')

% hold on;
% inds = 1:21:length(bc) ;
% quiver3(bc(inds, 1), bc(inds, 2), bc(inds, 3), ...
%     vfsm(inds, 1), vfsm(inds, 2), vfsm(inds, 3), 1) 
% % quiver3(bc(inds, 1), bc(inds, 2), bc(inds, 3), ...
% %     v0t(inds, 1), v0t(inds, 2), v0t(inds, 3), 1) 
% axis equal

% Take divergence and curl
        
% TODO: glue the mesh back together, FF will change
% cutMC is cutM that is Closed at the seam
[glueMesh, glue2cut] = glueCylinderCutMeshSeam(cutM) ;  
DEC = DiscreteExteriorCalculus( glueMesh.f, glueMesh.v ) ;

% Now resolve the vector field for decomposition
[v0n, v0t, v0t2d, jac3d_to_2d, ~, ~, dilation] = ...
    resolveTangentNormalVelocities(FF, v3drs, facevels, V2D, 1:length(FF)) ;

divv = DEC.divergence(v0t) ;
rotv = DEC.curl(v0t) ;

% Perform Helmholtz-Hodge decomposition
[divU, rotU, harmU, scalarP, vectorP] = ...
    DEC.helmholtzHodgeDecomposition(vfsm) ;

% DO LAPLACIAN SMOOTHING HERE
niter = 3 ;
fixed_verts = [] ;  % note: could use boundaries here, seems unnecessary
divvsm = laplacian_smooth(glueMesh.v, glueMesh.f, 'uniform', fixed_verts, ...
    0.1, 'explicit', divv, niter) ;

% weight curl by face area
% Note: hodge dual (star) is applied at end of Curl(), so no area exists in
% the curl/rotational field, so just average with any weighting we like (we
% think).
faceAreas = 0.5 * doublearea(glueMesh.v, glueMesh.f) ;
tri = triangulation(glueMesh.f, glueMesh.v) ;
% compute face neighbors 
faceNeighbors = [ (1:size(glueMesh.f,1)).', tri.neighbors() ];
% Handle the NaNs, which don't have four adjacent faces 
fNfix = faceNeighbors ;
badn = find(isnan(faceNeighbors)) ;
fNfix(badn) = 1 ;
rotvsm = rotv(fNfix) .* faceAreas(fNfix) ./ ...
    sum(faceAreas(fNfix), 2);
rotvsm(badn) = 0 ;


% add nU points back to divv from phi=0 by duplication --> convert back to
% cut mesh indices
divvCut = divvsm(glue2cut) ;
rotvCut = rotvsm(glue2cut) ;

% Coarse grain
% bcx = bc(:, 1) ;
% bcy = bc(:, 2) ;
% bcz = bc(:, 3) ;
% bcgridx = imresize(bcx(faceIDgrid), 1/qsub) ;
% bcgridy = imresize(bcy(faceIDgrid), 1/qsub) ;
% bcgridz = imresize(bcz(faceIDgrid), 1/qsub) ;
% divUgrid = divU(faceIDgrid) ;
% rotUgrid = rotU(faceIDgrid) ;
% harmUgrid = harmU(faceIDgrid) ;

% Pullback Vector Fields to Domain of Parameterization ====================
divU2d = zeros(size(FF, 1), 2);
for f = 1:size(FF,1)
    divU2d(f,:) = jac3d_to_2d{FF(f)} * divU(f,:)' / dilation(f) ;
end
rotU2d = zeros(size(FF, 1), 2);
for f = 1:size(FF,1)
    rotU2d(f,:) = jac3d_to_2d{FF(f)} * rotU(f,:)'  / dilation(f) ;
end
harmU2d = zeros(size(FF, 1), 2);
for f = 1:size(FF,1)
    harmU2d(f,:) = jac3d_to_2d{FF(f)} * harmU(f,:)' / dilation(f) ;
end

divs = {divvCut, divU, divU2d} ;
rots = {rotvCut, rotU, rotU2d} ;
harms = {harmU, harmU2D} ;

