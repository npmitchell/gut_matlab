%% NOTES -- Dillon

% Add path
addpath('/mnt/data/code/RicciFlow_MATLAB/')

% Ricci flow
% gut -> annulus
tp = 150 ;
cutMesh = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
cutMesh = cutMesh.spcutMeshSmRS ;
glueMesh = glueCylinderCutMeshSeam(cutMesh) ;

% Check it
trisurf(triangulation(glueMesh.f, glueMesh.v))
axis equal

% Generate conformal parameterization in the unit disk
[~, U, ~] = DiscreteRicciFlow.EuclideanRicciFlow(glueMesh.f, glueMesh.v, ...
    'BoundaryType', 'Fixed', 'BoundaryShape', 'Circles'); 
