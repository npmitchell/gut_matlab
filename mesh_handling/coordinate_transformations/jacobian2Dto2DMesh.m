function J_2D_To_2D = jacobian2Dto2DMesh(u2d, x2d, ff)
%JACOBIAN2DTO2DMESH(x2d, u3d, ff) 
%   Construct jacobian for mapping 2d->2d, from X->U.
%
% Parameters
% ----------
% u3d : #vertices x 2 float array
%   The target vertex locations in 2d (image space)
% x2d : #vertices x 2 float array
%   The source vertex locations in 2d (domain space)
% ff : #faces x 3 int array
%   connectivity list for triangulated meshes
% 
% Returns
% -------
% J_2D_To_2D : length(ff) x 1 cell array
%   A cell array containing the jacobian for each face as each element
%   The components of this matrix are the partial derivatives: du^i/dx^j
%   i.e., the Jacobian used to transform the CONTRAVARIANT components of a
%   tensor
% 
% Dillon Cislo, NPMitchell 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Construct Topological Structure Tools
%--------------------------------------------------------------------------
% Construct edge ID list --------------------------------------------------
tri2d = triangulation(ff, x2d) ;
EE = tri2d.edges;

% Construct face-edge correspondence tool ---------------------------------
% Given a list of scalar edge quantities, 'EQ', the output of
% 'EQ(feIDx(f,i))' is that quantity corresponding to the edge opposite the
% ith vertex in face f

% Look for the edge (32, 13, or 21) and be agnostic to the ordering (21=12)
e1IDx = sort( [ ff(:,3), ff(:,2) ], 2 );
e2IDx = sort( [ ff(:,1), ff(:,3) ], 2 );
e3IDx = sort( [ ff(:,2), ff(:,1) ], 2 );

[~, e1IDx] = ismember( e1IDx, EE, 'rows' );
[~, e2IDx] = ismember( e2IDx, EE, 'rows' );
[~, e3IDx] = ismember( e3IDx, EE, 'rows' );

feIDx = [ e1IDx e2IDx e3IDx ];

clear e1IDx e2IDx e3IDx

%% Calculate edge lengths in the domain space -----------------------------
L_E = x2d(EE(:,2),:) - x2d(EE(:,1),:);
L_E = sqrt(sum(L_E.^2, 2));

L_F = L_E(feIDx);

%% Calculate internal angles in the domain space --------------------------

% Some convenience variables to vectorize the cosine law calculation
Gi = L_F; Gj = circshift(L_F, [0 -1]); Gk = circshift(L_F, [0 -2]);

% The cosine of the internal angles
cosAng = ( Gj.^2 + Gk.^2 - Gi.^2 ) ./ ( 2 .* Gj .* Gk );

% The sine of the internal angles
sinAng = sin(acos(cosAng));

%% Calculate the (e1Hat, t1Hat)-basis in the domain space =================
% Calculate face unit normals  --------------------------------------------
normals = faceNormal( triangulation(ff, [x2d, zeros(size(x2d,1), 1)]) );

if (any(normals(:,1) ~= 0) || any(normals(:,2) ~= 0))
    error('Invalid domain triangulation supplied');
end

if (numel(unique(normals(:,3))) ~= 1)
    warning('Mesh faces are inconsistently ordered');
end

% Calculate the (e1Hat, t1Hat)-basis --------------------------------------
e1Hat = (x2d(ff(:,3),:) - x2d(ff(:,2),:)) ./ L_F(:,1);

t1Hat = cross(normals, [e1Hat, zeros(size(e1Hat,1), 1)], 2);
t1Hat = t1Hat(:, [1 2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct Jacobian Operator on Faces ===================================
% Calculate the components of the 2D Jacobian matrix ----------------------

% U-coordinates of triangles in 2D image space
u1 = u2d(ff(:,1), 1); u2 = u2d(ff(:,2), 1); u3 = u2d(ff(:,3), 1);
u23 = u3-u2; u21 = u1-u2;

% V-coordinates of triangles in 2D image space
v1 = u2d(ff(:,1), 2); v2 = u2d(ff(:,2), 2); v3 = u2d(ff(:,3), 2);
v23 = v3-v2; v21 = v1-v2;

J11 = u23 ./ L_F(:,1);
J21 = v23 ./ L_F(:,1);

J12 = u21 ./ L_F(:,3) - u23 .* cosAng(:,2) ./ L_F(:,1);
J12 = J12 ./ sinAng(:,2);

J22 = v21 ./ L_F(:,3) - v23 .* cosAng(:,2) ./ L_F(:,1);
J22 = J22 ./ sinAng(:,2);

% clear u1 u2 u3 u23 u21 u1 u2 u3 u23 u21 cosTheta2 sinTheta2

% Construct full Jacobian from 3D to 2D on each face ----------------------
J_2D_To_2D = cell(size(ff,1),1);

for f = 1:size(ff,1)
    
    % The coordinate change operator into the (e1Hat, t1Hat)-basis
    CC = [ e1Hat(f, :); t1Hat(f, :) ];
    
    % The 2D Jacobian -matrix
    JF = [ J11(f) J12(f); J21(f) J22(f) ];
    
    % Combine to construct complete Jacobian operator
    J_2D_To_2D{f} = JF * CC;
    
end

% clear J11 J21 J12 J22 
% clear CC JF f