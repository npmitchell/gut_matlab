function [labels, dbonds, topStructTools] = labelRectilinearMeshBonds(mesh, varargin)
% metricSPhiGridMesh(mesh, varargin)
%   Identify the bond on each face of a rectilinear grid mesh that is along
%   x and y (s and phi). Assumes input mesh is rectilinear with indices
%   increasing first along x, then along y. Mesh may be periodic in either
%   or both directions.
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
%       number of mesh vertices along each row of the rectilinear mesh
%   nV : int
%       number of mesh vertices along each column of the mesh
% vargin : positional arguments (optional)
%   eIDx : #bonds x 2 int
%       bond vertex IDs, bond vertices from topologicalStructureTools()
%   feIDx : #faces x 3 int
%       face bond IDs, face bonds from topologicalStructureTools()
%   bulkEdgeIDx : #bonds x 1 int
%       label of whether in bulk (0) or on edge (1)
%
% Returns
% -------
% labels : struct with fields
%   fe_is_u : #faces x 3 int array
%       whether a bond in the triangulation is along u (row)
%   fe_is_v : #faces x 3 int array
%       whether a bond in the triangulation is along v (column)
% dbonds : struct with fields
%   If mesh has no baseSpace (no mesh.u field), then output has fields:
%       u : #faces x dim float array
%           directed bond vector oriented along the 'u' direction of the
%           rectilinear mesh (rows, connecting (i, i+1) vertices)
%       v : #faces x dim float array
%           directed bond vector oriented along the 'v' direction of the
%           rectilinear mesh (columns, connecting (i, i+nU) vertices)
%   Otherwise, if mesh does have base space, there is a separate field for
%   each set of directed bond vectors
%       baseSpace : struct with fields u,v as above
%           directed bonds in domain of parameterization
%       realSpace : struct with fields u,v as above
%           directed bonds in coordinate space of mesh vertices
% topStructTools : optional cell array output, with elements
%   eIDx : #bonds x 2 int
%       bond vertex IDs, bond vertices from topologicalStructureTools()
%   feIDx : #faces x 3 int
%       face bond IDs, face bonds from topologicalStructureTools()
%   bulkEdgeIDx : #bonds x 1 int
%       label of whether in bulk (0) or on edge (1)
%
% NPMitchell 2020

% Store mesh vertices and faces
VV = mesh.v ;
FF = mesh.f ;
nU = mesh.nU ;
% nV = mesh.nV ;

% Construct Triangulation
tri = triangulation(FF, VV) ;

% Unpack or compute topological structure tools
load_tstools = true ;
if ~isempty(varargin)
    try
        eIDx = varargin{1} ;
        feIDx = varargin{2} ;
        load_tstools = false ;
    catch
        disp(['labelRectilinearMeshBonds: ', ...
            'topological structure tools could not be parsed'])
    end
end
if load_tstools
    [eIDx, feIDx, ~] = topologicalStructureTools(tri) ;
    if nargout > 2
        topStructTools = {eIDx, feIDx, } ;
    end
end

%% Count which bonds connect distant/adjacent/modular indices for 0/1/2
% Initial bonds define (s,phi) coordinate system
% Are they shat (1), phihat (2), or diagonal (0)?
% Build sorp==1 is shat, sorp==2 is phihat, sorp==0 is diagonal
sorp = zeros(size(eIDx(:, 1))) ;
% shat will differ by one
bondDX = abs(eIDx(:, 2) - eIDx(:, 1)) ;
sorp(bondDX == 1) = 1 ;
% phihat will differ by multiple of nU
sorp(mod(bondDX, nU) == 0) = 2 ;

% fe_is_s is boolean, true where feIDx element is along s
fe_is_u = sorp(feIDx) == 1 ;    
% fe_is_phi is boolean, true where feIDx element is along phi
fe_is_v = sorp(feIDx) == 2 ;

% check that there is one shat vector in each triangle
assert(all(any(fe_is_u, 2)))  
% check that there is one phihat vector in each triangle
assert(all(any(fe_is_v, 2)))  

labels.fe_is_u = fe_is_u ;
labels.fe_is_v = fe_is_v ;

if nargout > 1
    % Directed edge vectors in mesh configuration space
    eij = VV(eIDx(:,2), :) - VV(eIDx(:,1), :);
    eij = VV(eIDx(:,2), :) - VV(eIDx(:,1), :);

    % Return directed bond vectors for bonds along u and along v
    dbond_u = eij(sum(fe_is_u .* feIDx, 2), :) ;
    dbond_v = eij(sum(fe_is_v .* feIDx, 2), :) ;
    dbonds.u = dbond_u ;
    dbonds.v = dbond_v ;
    
    % Check for a base space of the mesh
    if isfield(mesh, 'u')
        UU = mesh.u ;
        
        % Directed edge vectors in mesh base space
        eij = UU(eIDx(:,2), :) - UU(eIDx(:,1), :);
        eij = UU(eIDx(:,2), :) - UU(eIDx(:,1), :);

        % Return directed bond vectors for bonds along u and along v
        dbond_u = eij(sum(fe_is_u .* feIDx, 2), :) ;
        dbond_v = eij(sum(fe_is_v .* feIDx, 2), :) ;
        dbonds2.u = dbond_u ;
        dbonds2.v = dbond_v ;
        
        % Return as struct with fields realSpace and baseSpace
        dbonds.realSpace = dbonds ;
        dbonds.baseSpace = dbonds2 ;
    end
end

