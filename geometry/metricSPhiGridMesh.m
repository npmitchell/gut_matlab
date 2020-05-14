function gg = metricSPhiGridMesh(mesh, varargin)
% metricSPhiGridMesh(mesh, varargin)
%   Compute the metric tensor for each face of a rectilinear grid mesh
%   whose bonds are oriented along the s and phi 
%
% todo: allow passing topologicalStructureTools instead of recompute
%
% Parameters
% ----------
% mesh : struct with fields
%   v : N x 3 float
%       mesh vertices
%   f : M x 3 int
%       mesh faces, as defined by indexing into vertices. Assumes
%       consistent face orientation throughout the mesh
%   nU : int
%   nV : int
% 
% Returns
% -------
% gg : #faces x 3 float array
%   the metric tensor for each face, as 3-column array, with cols
%   g_ss, g_sphi, g_phiphi. Assumes g_sphi = g_phis
%
% NPMitchell 2020

% Store mesh vertices and faces
VV = mesh.v ;
FF = mesh.f ;
nU = mesh.nU ;
nV = mesh.nV ;

% Construct Triangulation
tri = triangulation(FF, VV) ;

[eIDx, feIDx, ~] = topologicalStructureTools(tri) ;

%% Construct Physical/Target Configurations and Geometries ================
% Construct Initial Configuration Metric -------------------------------
% Directed edge vectors
eij = VV(eIDx(:,2), :) - VV(eIDx(:,1), :);

% Initial bonds define (s,phi) coordinate system
% Are they shat (1), phihat (2), or diagonal (0)?
% Build sorp==1 is shat, sorp==2 is phihat, sorp==0 is diagonal
sorp = zeros(size(eIDx(:, 1))) ;
% shat will differ by one
bondDX = abs(eIDx(:, 2) - eIDx(:, 1)) ;
sorp(bondDX == 1) = 1 ;
% phihat will differ by multiple of nU
sorp(mod(bondDX, nU) == 0) = 2 ;

% Current edge lengths
eL0 = sqrt( sum( eij.^2, 2 ) );

%% Current metric 
% build metric for each triangle: [ss sphi phiphi]
gg = zeros(size(feIDx)) ;

% fe_is_s is boolean, true where feIDx element is along s
fe_is_s = sorp(feIDx) == 1 ;    
% fe_is_phi is boolean, true where feIDx element is along phi
fe_is_phi = sorp(feIDx) == 2 ;

% check that there is one shat vector in each triangle
assert(all(any(fe_is_s, 2)))  
% check that there is one phihat vector in each triangle
assert(all(any(fe_is_phi, 2)))  

% The bond index along s = sum(fe_is_s .* feIDx, 2) for each triangle
gg(:, 1) = eL0(sum(fe_is_s .* feIDx, 2)) ;
gg(:, 3) = eL0(sum(fe_is_phi .* feIDx, 2)) ;

% Now for the cross-term, grab directed bonds in s and phi
dbond_s = eij(sum(fe_is_s .* feIDx, 2), :) ;
dbond_phi = eij(sum(fe_is_phi .* feIDx, 2), :) ;
% note there is a minus sign so that the tails of both vectors
% share a vertex.
gg(:, 2) = dot(dbond_s', -dbond_phi')' ;

