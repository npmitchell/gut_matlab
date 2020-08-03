function cmesh = closeRectilinearCylMesh(mesh)
%CLOSERECTILINEARCYLMESH(mesh)
% Given a cut mesh (topologically, a square with implied periodicity in Y),
% glue the mesh back together into a topological cylinder.
%
% Parameters
% ----------
% mesh : struct with fields v,u,vn
%   mesh whose vertices are defined in u plane by 1...nU cols, 1...nV rows
%
% Returns
% -------
% cmesh : struct with fields v, u, vn
%   closed mesh whose vertices are defined as before plus two endcap
%   vertices
% 
% NPMitchell 2020

nU = mesh.nU ;
nV = mesh.nV ;

% Make the u and v fields are list of vectors not gridded structures 
if ~size(mesh.u, 2) == 2 || any(size(mesh.u) == nU)
    mesh.u = reshape(mesh.u, [nU * nV, 2]) ;
end
if ~size(mesh.v, 2) == 3 || any(size(mesh.v) == nU)
    mesh.v = reshape(mesh.v, [nU * nV, 3]) ;
end

% check that uv has increasing u, then increasing v
assert(mesh.u(1, 1) ~= mesh.u(2, 1))
assert(mesh.u(1, 2) == mesh.u(2, 2))

% Close the seam connecting u(:, 2) == 0 to u(:, 2) == 1
cmesh.v = mesh.v(1:end-nU, :) ;
% Replace the final row with the first row. First send the last row from
% (nV-1)*nU+1:nV*nU to 0:nU-1. Then send 0:nU-1 to 1:nU.
cmesh.f = mod(mesh.f, (nV-1)*nU + 1) ;
cmesh.f(mesh.f > (nV-1)*nU) = cmesh.f(mesh.f > (nV-1)*nU) + 1 ;
cmesh.u = mesh.u(1:end-nU, :) ;

% Adjust normals field if it exists
if isfield(cmesh, 'vn')
    cmesh.vn = mesh.vn(1:end-nU, :) ;
end

apt = mean(mesh.v(1:nV:end, :), 1) ;
ppt = mean(mesh.v(nU:nV:end, :), 1) ;

% add front and rear faces
for pp = 1:2
    if pp == 1
        % front faces
        endpts = 1:nV:length(mesh.v) ;
        cmesh.v = cat(1, cmesh.v, apt) ;
    else
        % back faces
        endpts = nU:nV:length(mesh.v) ;
        cmesh.v = cat(1, cmesh.v, ppt) ;
    end
    aidx = length(cmesh.v);
    frontf = zeros(length(endpts), 3) ;
    for qq = 1:length(endpts)
        nid = mod(qq + 1, length(endpts)) ;
        if nid == 0
            nid = length(endpts) ;
        end
        frontf(qq, :) = [endpts(qq), aidx, endpts(nid)];
    end
    cmesh.f = cat(1, cmesh.f, frontf) ;
end

% Provide normals to the added points at endcaps. If normals were included
% in the original, then only update the endcap normals. Otherwise define
% per_vertex_normals.
if isfield(cmesh, 'vn')
    assert(all(size(cmesh.v) == [nU*(nV-1)+2, 3]))
    vn2 = per_vertex_normals(cmesh.v, cmesh.f, 'Weighting', 'angle') ;
    cmesh.vn(nU*(nV-1)+1:end, :) = vn2(nU*(nV-1)+1:end, :) ;
else
    cmesh.vn = per_vertex_normals(cmesh.v, cmesh.f, 'Weighting', 'angle') ;
end
