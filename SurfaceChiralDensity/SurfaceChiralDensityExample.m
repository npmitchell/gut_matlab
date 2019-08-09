%% SURFACE CHIRAL DENSITY EXAMPLE =========================================
% An example detailing the functionality of the discrete surface chiral
% density pseudotensor.
% By Dillon Cislo 02/26/2019
%==========================================================================

clear; close all; clc;

%% CONSTRUCT HELICAL STRIP MESH ===========================================

M = 50; % Number of rows (should always be even)
N = 100; % Number of columns

xlim = [0 20]; xdist = xlim(2)-xlim(1); 
ylim = [-2.5 2.5]; ydist = ylim(2)-ylim(1);

% Construct 2D vertices----------------------------------------------------
[X,Y] = meshgrid( xlim(1):xdist/(N-1):xlim(2), ...
    ylim(1):ydist/(M-1):ylim(2) );


v2D = [ X(:), Y(:) ];

clear X Y xlim ylim xdist ydist

% Construct face connectivity list ----------------------------------------
rhv = (1:(M-1))';
rhv = repmat( rhv, 1, (N-1) );
rhv = rhv + M .* meshgrid(0:(N-2),1:(M-1));
rhv = rhv(:);
rhf = [ rhv, rhv+1, rhv+M ];

lhv = (2:M)';
lhv = repmat( lhv, 1, (N-1) );
lhv = lhv + M .* meshgrid(0:(N-2),1:(M-1));
lhv = lhv(:);
lhf = [ lhv, lhv+M, lhv+(M-1) ];

face = [rhf;lhf];

clear rhv rhf lhv lhf

% Construct 3D vertices ---------------------------------------------------
vertex = [ v2D, zeros(size(v2D,1),1) ];

% We rotate the vertices around the X-axis by an x-dependent rotation angle
c = -0.5;
thetaX = c .* vertex(:,1);

xHat = [ ones(size(v2D,1),1), ...
    zeros(size(v2D,1),1), zeros(size(v2D,1),1) ];

vertex = vertex .* cos( thetaX ) + ...
    cross( xHat, vertex, 2 ) .* sin( thetaX ) + ...
    ( 1 - cos( thetaX ) ) .* vertex(:,1) .* xHat;

helixTri = triangulation( face, vertex );

clear thetaX xHat

% View Result -------------------------------------------------------------
trisurf( helixTri ); axis equal
xlabel('x')
ylabel('y')
zlabel('z')

%% CONSTRUCT ANALYTIC CHIRAL DENSITY TENSOR ===============================

% Construct the tensor on each vertex -------------------------------------
trueChiralVertex = cell(size(v2D,1),1);
for vID = 1:size(v2D,1)
    
    u = v2D(vID,1); v = v2D(vID,2);
    
   chi = [ 1, -c*v*sin(c*u), c*v*cos(c*u); ...
       -c*v*sin(c*u), -(1+(1+2*c^2*v^2)*cos(2*c*u))/2, -(1+2*c^2*v^2)*sin(2*c*u)/2; ...
       c*v*cos(c*u), -(1+2*c^2*v^2)*sin(2*c*u)/2, -(1-(1+2*c^2*v^2)*cos(2*c*u))/2 ];
   
   chi = c .* chi ./ ( (1 + c^2 * v^2 ).^2 );
   
   trueChiralVertex{vID} = chi;
    
end

% Average the tensor onto each face ---------------------------------------
trueChiral = cell( size(face,1), 1 );
for fID = 1:size(face,1)
    
    chi = trueChiralVertex{ face(fID,1) } + ...
        trueChiralVertex{ face(fID,2) } + ...
        trueChiralVertex{ face(fID,3) };
    
    trueChiral{fID} = chi ./ 3;
    
end

clear chi u v

%% CONSTRUCT CHIRAL DENSITY PSEUDOTENSOR ==================================

chiral = calculate_chiral_density( face, vertex );

%% VIEW RESULTS ===========================================================

% Create surface vector field based on faces ------------------------------
rhv = vertex( face( 1:((M-1)*(N-1)), 2 ), : ) - ...
    vertex( face( 1:((M-1)*(N-1)), 1 ), : );

lhv = vertex( face( ((M-1)*(N-1)+1):size(face,1), 2 ), : ) - ...
    vertex( face( ((M-1)*(N-1)+1):size(face,1), 3 ), : );

fvf = [ rhv; lhv ]; 
% clear rhv lhv

e1 = vertex(face(:,3),:) - vertex(face(:,2),:);
e2 = vertex(face(:,1),:) - vertex(face(:,3),:);

% The face unit normal vector
n = cross( e1, e2, 2 );
n = n ./ sqrt( sum( n.^2, 2 ) );
clear e1 e2

fvf = fvf ./ sqrt( sum( fvf.^2, 2 ) );
fvf = fvf + n;
fvf = fvf ./ sqrt( sum( fvf.^2, 2 ) );

% View surface vector field -----------------------------------------------
CoM = cat( 3, vertex(face(:,1),:), ...
    vertex(face(:,2),:), vertex(face(:,3),:) );
CoM = mean( CoM, 3 ); % Face centroids

trisurf(helixTri,'FaceColor',[1.0 1.0 1.0]);
axis equal
hold on
% quiver3( CoM(:,1), CoM(:,2), CoM(:,3), ...
%     fvf(:,1), fvf(:,2), fvf(:,3), ...
%     0.5, 'Color', 'b' );
quiver3( CoM(1:((M-1)*(N-1)), 1), ...
    CoM(1:((M-1)*(N-1)), 2), ...
    CoM(1:((M-1)*(N-1)), 3), ...
    rhv(:,1), rhv(:,2), rhv(:,3), ...
    0.2, 'Color', 'b' );
quiver3( CoM(((M-1)*(N-1)+1):size(face,1), 1), ...
    CoM(((M-1)*(N-1)+1):size(face,1), 2), ...
    CoM(((M-1)*(N-1)+1):size(face,1), 3), ...
    lhv(:,1), lhv(:,2), lhv(:,3), ...
    0.2, 'Color', 'r' );

quiver3( CoM(:, 1), ...
    CoM(:, 2), ...
    CoM(:, 3), ...
    fvf(:,1), fvf(:,2), fvf(:,3), ...
    0.2, 'Color', 'g' );

clear CoM

%%
% Act the chirality operators on each face vector -------------------------
trueChi = zeros( size(face,1), 1 );
chi = zeros( size(face,1), 1 );
for i = 1:size(face,1)
    
    trueChi(i) = fvf(i,:) * trueChiral{i} * fvf(i,:)';
    chi(i) = fvf(i,:) * chiral{i} * fvf(i,:)';
    
    % test for using embedding space vector rather than surface vector
    % tmp = chiral{i} ;
    %chi(i) = tmp(2, 2);
    
end

% The relative error ------------------------------------------------------
chiError = sign(chi) .* sign(trueChi) .* ...
    abs( ( chi - trueChi ) ./ trueChi );

if any( chiError(:) < 0 )
    warning("A sign discrepancy occurred!");
end

% Get the relative error without the boundary faces -----------------------
faceNeighbors = helixTri.neighbors;
bdyFaceIDx = any(isnan(faceNeighbors),2);

maxErrNoBdy = max(chiError(~bdyFaceIDx));
chiErrNoBdy = chiError;
chiErrNoBdy( bdyFaceIDx ) = maxErrNoBdy;

clear faceNeighbors maxErrNoBdy

% View the error ----------------------------------------------------------

% Error with boundary faces
% patch( 'Faces', face, 'Vertices', vertex, ...
%     'FaceVertexCData', chiError, 'FaceColor', 'flat' );
% axis equal
% colorbar

% Error not including boundary faces
patch( 'Faces', face, 'Vertices', vertex, ...
    'FaceVertexCData', chiErrNoBdy, 'FaceColor', 'flat' );
axis equal
colorbar
    