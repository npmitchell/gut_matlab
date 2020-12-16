%% NOTES -- Dillon

% Add path
addpath('../../RicciFlow_MATLAB/')

% Ricci flow
% gut -> annulus
tp = 160 ;
% cutMesh = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
cutMesh = load('./spcutMeshSmRS_000160.mat')
cutMesh = cutMesh.spcutMeshSmRS ;
glueMesh = glueCylinderCutMeshSeam(cutMesh) ;

% Check it
trisurf(triangulation(glueMesh.f, glueMesh.v))
axis equal

% Generate conformal parameterization in the unit disk
[~, U, ~] = DiscreteRicciFlow.EuclideanRicciFlow(glueMesh.f, glueMesh.v, ...
    'BoundaryType', 'Fixed', 'BoundaryShape', 'Circles', ...
    'MaxCircIter', 50);
% [labels, dbonds, topStructTools] = labelRectilinearMeshBonds(cutMesh) ;
% Let MaxCircIter >= 50
triplot(triangulation(glueMesh.f, U))


%% Beltrami
mu = bc_metric(glueMesh.f, U, glueMesh.v, 3) ;

% plot beltrami
mesh2d = triangulation(glueMesh.f, [U(:, 1), U(:, 2)]) ; %  ; , 0*U(:, 1)]) ;
m2plot = triangulation(glueMesh.f, [U(:, 1), U(:, 2), 0*U(:, 1)]) ;
options.view = [0, 90] ;
nFieldsOnSurface(m2plot, {real(mu), imag(mu), abs(mu)}, ...
    options)
set(gcf, 'visible', 'on')

%% Plot just the bare triangulation
close all
triplot(triangulation(glueMesh.f, U))
axis equal
scatter(0,0,50, 'filled', 'r')

%% Push barycenter of inner hole to origin via Mobius transformation
boundaries = DiscreteRicciFlow.compute_boundaries(glueMesh.f) ;
% which boundary is the inner one?
LL = [0, 0] ;
for qq = 1:length(boundaries)
    % make periodic 
    LL(qq) = sum(sqrt(sum((U(boundaries{qq},:)- ...
        circshift(U(boundaries{qq},:), 1, 1)).^2, 2))) ;
end
% which is smaller?
[~, innerID] = min(LL) ;
inner = boundaries{innerID} ;
outer = boundaries{mod(innerID+1, 2)} ;


% barycenter of inner -- THIS IS BIASED
% baryc = mean(U(inner, :), 1) ;
% baryc = complex(baryc(1), baryc(2)) ;

% Take center (not barycenter) of circle
% note: https://www.mathworks.com/matlabcentral/fileexchange/5557-circle-fit
[xc,yc, innerR] = circfit(U(inner, 1), U(inner, 2)); 
baryc = xc + 1j * yc ;

% Enforce circularity in inner Boundary
error('todo')

% Enforce circularity in outer boundary ==> radius=1
error('todo')

% Covert U to complex
zz = complex(U(:, 1), U(:, 2)) ;

% Mobius
zz = (zz - baryc) ./ (1 - conj(baryc) .* zz) ;
UU = [real(zz), imag(zz) ] ;

% inspect centered mesh
triplot(triangulation(glueMesh.f, UU))
hold on;
plot(UU(boundaries{1}, 1), UU(boundaries{1}, 2), 'ko-')
plot(UU(boundaries{2}, 1), UU(boundaries{2}, 2), 'co-')
scatter(0,0, 'r', 'filled')

%% Branch cut

% THIS DOES NOT WORK
% branch = [0,0, 2, 0] ;
% E = triangulation(glueMesh.f, UU).edges ; 
% xy2 = [UU(E(:,1), :), UU(E(:,2), :)] ;
% intersections = lineSegmentIntersect(branch, xy2) ;
% % which lines intersect
% linesThatIntersect = find(intersections.intAdjacencyMatrix) ;
% % order them by increasing x coordinate along branch 
% xs = intersections.intMatrixX(linesThatIntersect) ;
% [xsort, sortIDs] = sort(xs) ;
% cutPath = E(linesThatIntersect(sortIDs), 1) ;
% 
% close all
% triplot(triangulation(glueMesh.f, UU), 'Color', [0.8,0.8,0.8]) ;
% hold on;
% plot(UU(cutPath, 1), UU(cutPath, 2), 'ro-') ;
% plot([0, 1], [0, 0], 'k--')

% USE MATLAB graph 

%% Take log
zz = complex(UU(:, 1), UU(:, 2)) ;
rr = real(log(zz)) ;
phi = imag(log(zz)) ;
triplot(triangulation(glueMesh.f, [rr, phi])) ;

%% Think about variation
[labels, dbonds, topStructTools] = labelRectilinearMeshBonds(cutMesh) ;
tmp = vecnorm(dbonds.realSpace.v, 2, 2) ;
trisurf(triangulation(cutMesh.f, cutMesh.v), 'FaceColor', tmp)
axis equal




