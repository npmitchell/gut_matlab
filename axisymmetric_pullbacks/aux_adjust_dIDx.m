function [adIDx, pdIDx] = aux_adjust_dIDx(mesh, t, dpFile, ADBase, PDBase, cylinderMeshCleanBase, outadIDxfn, outpdIDxfn, xp) 

% Auxilliary function for adjusting adIDx and pdIDx in
% Generate_Axisymmetric_Pullbacks_Orbifold.m script

% Load the AD/PD vertex IDs
disp('Loading ADPD vertex IDs...')
if t == xp.fileMeta.timePoints(1)
    adIDx = h5read( dpFile, sprintf( ADBase, t ) );
    pdIDx = h5read( dpFile, sprintf( PDBase, t ) );

    ad3D = mesh.v( adIDx, : );
    pd3D = mesh.v( pdIDx, : );
else
    % Load previous mesh and previous adIDx, pdIDx
    prevcylmeshfn = sprintf( cylinderMeshCleanBase, t -1 ) ;
    prevmesh = read_ply_mod( prevcylmeshfn ); 
    prevadIDx = h5read(outadIDxfn, ['/' sprintf('%06d', t-1) ]) ;
    % read previous pdIDx with new indices
    prevpdIDx = h5read(outpdIDxfn, ['/' sprintf('%06d', t-1) ]) ;
    ad3D = prevmesh.v(prevadIDx, :) ;
    pd3D = prevmesh.v(prevpdIDx, :) ;
end

trngln = triangulation(mesh.f, mesh.v) ;
boundary = trngln.freeBoundary ;
adIDx = boundary(pointMatch( ad3D, mesh.v(boundary(:, 1), :) ), 1);
pdIDx = boundary(pointMatch( pd3D, mesh.v(boundary(:, 1), :) ), 1);