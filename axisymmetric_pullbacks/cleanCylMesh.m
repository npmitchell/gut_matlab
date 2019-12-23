function [mesh] = cleanCylMesh(mesh)
%CLEANCYLMESH(mesh) Orient faces, smooth normals, remove ears 
%   Detailed explanation goes here
%
% Parameters
% ----------
% mesh : cylinderMesh with fields v, f
% 
% Returns
% -------
% mesh : cleaned cylinderMesh with fields v, f, vn
% 
% NPMitchell 2019

% Consistently orient mesh faces
disp('orienting faces')
mesh.f = bfs_orient( mesh.f );

% Compute vertex normals by weighting the incident triangles by
% their incident angles to the vertex
% This uses gptoolbox.
mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
% Average normals with neighboring normals
disp('Averaging normals with neighboring normals')
mesh.vn = average_normals_with_neighboring_vertices(mesh, 0.5) ;

% Clip the ears of the triangulation and update the AD/PD points if
% necessary -----------------------------------------------------------
disp('Clipping ears...') ;
[ ff, ~, vv, newind] = clip_mesh_mod( mesh.f, mesh.v );
mesh.f = ff; mesh.v = vv;
% Remove zeros from newind in order to assign the new vtx normals
newind_keep = newind(newind > 0) ;
if all(diff(newind_keep) > 0)
    mesh.vn = mesh.vn(newind > 0, :) ;
else
    error('Index array is non-monotonic, handle slightly fancier case here')
end


% % Check it and the normals
% tmp = mesh ;
% gone = find(newind == 0) ;
% inds = 1:min(gone)
% ind = min(gone) ;
% tmp.v = tmp.v + tmp.vn * (130) ;
% trisurf(tmp.f, tmp.v(:, 1), tmp.v(:, 2), tmp.v(:, 3), 'EdgeColor', 'none')
% hold on
% axis equal
%
% Show which points were removed
% plot3(tmp0.v(gone, 1), tmp0.v(gone, 2), tmp0.v(gone, 3), 'ro')
% axis equal
% 
% tmp2 = mesh ;
% % plot3(tmp2.v(:, 1), tmp2.v(:, 2), tmp2.v(:, 3), '.')
% trisurf(tmp2.f, tmp2.v(:, 1), tmp2.v(:, 2), tmp2.v(:, 3), 'EdgeColor', 'none')
% hold on;
% quiver3(tmp2.v(:, 1), tmp2.v(:, 2), tmp2.v(:, 3), tmp2.vn(:, 1), tmp2.vn(:, 2), tmp2.vn(:, 3), 10)
% %quiver3(tmp2.v(inds, 1), tmp2.v(inds, 2), tmp2.v(inds, 3), tmp2.vn(inds, 1), tmp2.vn(inds, 2), tmp2.vn(inds, 3), 10)
% quiver3(tmp2.v(ind, 1), tmp2.v(ind, 2), tmp2.v(ind, 3), tmp2.vn(ind, 1), tmp2.vn(ind, 2), tmp2.vn(ind, 3), 10, 'color', 'r')
% axis equal
% mags = vecnorm(mesh.vn - mesh2.vn(1:length(mesh.vn), :), 2, 1)


end

