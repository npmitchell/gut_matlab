%% DISCRETE CURVATURE TEST ================================================
% by Dillon Cislo 02/26/2019

clear; close all; clc

% Construct hemispherical mesh
sphereTri = sphericalTriangulation( 'NumIterations', 5 );

face = sphereTri.ConnectivityList;
vertex = sphereTri.Points;


%% CALCULATE GAUSSIAN CURVATURE ===========================================

[ K, Kr ] = calculate_gaussian_curvature( face, vertex );


%% CALCULATE MEAN CURVATURE ===============================================
[ Habs, Hn ] = calculate_mean_curvature( face, vertex );

%% EXAMPLE ================================================================
% Here we show to to compute the signed mean curvature given an input set
% of vertex normals
%--------------------------------------------------------------------------

% Outward pointing unit normals on the sphere
vn = vertex ./ sqrt( sum( vertex.^2, 2 ) );

% The mean curvature unit normals
nn = -Hn ./ sqrt( sum( Hn.^2, 2 ) );

% The sign of the mean curvature
sgnH = sign( dot( vn, nn, 2 ) );

% The signed mean curvature
H = sgnH .* Habs;

% View mean curvature normals ---------------------------------------------
trisurf(sphereTri,'FaceColor',[0.8 0.8 1.0]);
axis equal
hold on
quiver3( vertex(:,1), vertex(:,2), vertex(:,3), ...
    nn(:,1), nn(:,2), nn(:,3), ...
    0.5, 'Color', 'b' );

