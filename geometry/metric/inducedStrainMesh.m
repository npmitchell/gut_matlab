function [strain, tre, dev, theta, outputStruct] = ...
    inducedStrainMesh(deformedMesh, refMesh, options)
%strainRateMesh(cutMesh, dxdy, options)
%   Compute strain tensor on faces of a mesh given the rest lengths of the 
%   mesh bonds in the fixed/reference Lagrangian frame.
%       ==> strain = 0.5 * (inv(g0) * (g1 - g0))
%
%   DO NOT tile the mesh before making computations. 
%   Returns strain tensors on each face, the trace
%   (note this is the full trace, not 1/2 trace) of the strain rate in the
%   embedding Tr[inv(g) epsilon], the deviatoric component magnitude
%   (Frobenius norm of the deviator, which is
%   A = epsilon - (1/2)Tr[inv(g) epsilon] g, the angle of elongation of the
%   deviator in the embedding coordinates relative to the projection of the
%   zeta direction in the embedding coordinates, and additional output. 
%   Additional output returns angles theta_pb in the coordinates
%   of the pullback space (instead of in the embedding space relative
%   to the embedding-space projection of zeta_hat), the scaling factors
%   bondDxDy used to scale the eigenvectors in the theta determination, the
%   fundamental forms fundForms, and the symmetrized gradient of the
%   velocities.
%
% 
% Parameters
% ----------
% deformedMesh : struct with fields 
%   Mesh with pullback space u and embedding space v. 
%   f : #faces x 3 int array
%       mesh face connectivity list 
%   v : #vertices x 3 float array
%       embedding space mesh vertices of the cut mesh
%   u : #vertices x 2 float array
%       pullback mesh vertices of the cut mesh
% refMesh : struct with fields 
%   Mesh with pullback space u and embedding space v. 
%   f : #faces x 3 int array
%       mesh face connectivity list 
%   v : #vertices x 3 float array
%       embedding space mesh vertices of the cut mesh
%   u : #vertices x 2 float array
%       pullback mesh vertices of the cut mesh
% options : optional struct with fields
%   passed to labelRectilinearMeshBonds()
%
% Returns
% -------
% strainrate : #faces x 1 cell array of 2x2 float matrices
% tre : #faces x 1 float array
% dev : #faces x 1 float array
% theta : #faces x 1 float array
% outputStruct
%   fundForms : struct with fields
%       gg : first fundamental form on each face for tiled mesh
%       bb : second fundamental form on each face for tiled mesh
%   bondDxDy : struct with fields
%       dx : length of dx in embedding space / length of dx in pullback space
%       dy : length of dy in embedding space / length of dy in pullback space
%   theta_pb : #tiledFaces x 1 float array
%       elongation angle, theta, in pullback space (differs from
%       theta by scaling of the eigenvectors by (dx, dy) returned in
%       bondDxDy
%
% See Also
% --------
% inducedStrainPeriodicMesh -- for 1d-periodic boundary conditions / cylindrical-like meshes
%
% NPMitchell 2022

%% Default options
debug = false ;
mesh = [] ;

%% Unpack options
if nargin > 2
    if isfield(options, 'debug')
        debug = options.debug ;
    end
else
    options = struct() ;
end

%% Compute the fundamental forms
[g0cell, b0cell] = constructFundamentalForms(refMesh.f, refMesh.v, refMesh.u) ;
% Use tiled, open mesh (glued seam) to compute the second fundamental form
% [~, b0cell] = constructFundamentalForms(TF0, TV3D0, TV2D0) ;

[g1cell, b1cell] = constructFundamentalForms(deformedMesh.f, deformedMesh.v, deformedMesh.u) ;
% Use tiled, open mesh (glued seam) to compute the second fundamental form
% [~, b1cell] = constructFundamentalForms(TF1, TV3D1, TV2D1) ;

%% Strain tensor
% Metric for reference mesh
% Could instead use different function: 
% g0cell = inducedMetric(refCutMesh.f, refCutMesh.v, refCutMesh.u) ;
% Convert metric tensor to #Faces x 2 x 2
g0 = zeros(size(g0cell, 1), 2, 2) ;
for qq = 1:size(g0cell, 1) 
    g0(qq, :, :) = g0cell{qq} ;
end
b0 = zeros(size(b0cell, 1), 2, 2) ;
for qq = 1:size(b0cell, 1) 
    b0(qq, :, :) = b0cell{qq} ;
end

% Metric for deformed mesh
% g1cell = inducedMetric(deformedCutMesh.f, deformedCutMesh.v, refCutMesh.u) ;
% Convert metric tensor to #Faces x 3
g1 = zeros(size(g1cell, 1), 2, 2) ;
for qq = 1:size(g1cell, 1)
    g1(qq, :, :) = g1cell{qq} ;
end
b1 = zeros(size(b1cell, 1), 2, 2) ;
for qq = 1:size(b1cell, 1)
    b1(qq, :, :) = b1cell{qq} ;
end

% Metric strain
strain = zeros(size(g0)) ;
gdiffs = zeros(size(g0)) ;
for qq = 1:size(g0cell, 1)
    g0q = g0cell{qq} ;
    g1q = g1cell{qq} ;
    strain(qq, :, :) = 0.5 * (inv(g0q) * (g1q - g0q)) ;
    gdiffs(qq, :, :) = 0.5 * (g1q - g0q) ;
    % tmp = inv(g0q) * g0q ;
    % assert(abs(tmp(1) - 1)  < 1e-7)
    % assert(abs(tmp(2))  < 1e-7)
    
end

% disp('compute strain here')

%% Find dx and dy -- lengths of projected / lengths of pullback
[~, dbonds_ref, ~] = labelRectilinearMeshBonds(refMesh, options) ;
% [~, dbonds_def, ~] = labelRectilinearMeshBonds(defTiled, options) ;
dx = (vecnorm(dbonds_ref.realSpace.u, 2, 2)) ./ vecnorm(dbonds_ref.baseSpace.u, 2, 2) ;
dy = (vecnorm(dbonds_ref.realSpace.v, 2, 2)) ./ vecnorm(dbonds_ref.baseSpace.v, 2, 2) ;

%% Metric strain -- separate trace and deviatoric strain comp, angle
tre = zeros(size(strain, 1), 1) ;  % traceful dilation
dev = zeros(size(strain, 1), 1) ;  % deviatoric magnitude
theta = zeros(size(strain, 1), 1) ;  % angle of elongation
theta_pb = zeros(size(strain, 1), 1) ;  % angle of elongation in PB space
for qq = 1:size(strain, 1)
    eq = squeeze(gdiffs(qq, :, :)) ;
    gq = g0cell{qq} ;
    
    %% Trace / deviator / theta
    % take angle of deviator to be in embedding space, angle relative
    % to the embedding-space projection of zeta_hat. 
    try
        [tre(qq), dev(qq), theta(qq), theta_pb(qq)] = ...
            traceDeviatorPullback(eq, gq, dx(qq), dy(qq)) ;
    catch
        disp('what is going wrong?')
    end
end
% Modulo is not necessary
% theta = mod(theta, pi) ;
% theta_pb = mod(theta_pb, pi) ;

%% Output the fundamental forms, embedding bond lengths, and theta_pb
if nargout > 4
    %% Find dx and dy -- lengths of projected / lengths of pullback
    [~, dbonds_ref, ~] = labelRectilinearMeshBonds(refMesh, options) ;
    [~, dbonds_def, ~] = labelRectilinearMeshBonds(deformedMesh, options) ;
    % dx_ref = vecnorm(dbonds_ref.realSpace.u, 2, 2) ./ vecnorm(dbonds_ref.baseSpace.u, 2, 2) ;
    % dy_ref = vecnorm(dbonds_ref.realSpace.v, 2, 2) ./ vecnorm(dbonds_ref.baseSpace.v, 2, 2) ;
    % dx_def = vecnorm(dbonds_def.realSpace.u, 2, 2) ./ vecnorm(dbonds_def.baseSpace.u, 2, 2) ;
    % dy_def = vecnorm(dbonds_def.realSpace.v, 2, 2) ./ vecnorm(dbonds_def.baseSpace.v, 2, 2) ;
    dx_strain = (vecnorm(dbonds_def.realSpace.u, 2, 2) - vecnorm(dbonds_ref.realSpace.u, 2, 2)) ./ vecnorm(dbonds_ref.realSpace.u, 2, 2) ;
    dy_strain = (vecnorm(dbonds_def.realSpace.v, 2, 2) - vecnorm(dbonds_ref.realSpace.v, 2, 2)) ./ vecnorm(dbonds_ref.realSpace.v, 2, 2) ;

    % fundamental forms
    fundForms = struct() ;
    fundForms.g0 = g0 ;
    fundForms.g1 = g1 ;
    fundForms.b0 = b0 ;
    fundForms.b1 = b1 ;
    
    % ratio of embedding bond lengths to reference bond lengths
    bondDxDy = struct() ;
    bondDxDy.dbonds_ref = dbonds_ref ;
    bondDxDy.dbonds_def = dbonds_def ;
    bondDxDy.dx_strain = dx_strain ;
    bondDxDy.dy_strain = dy_strain ;
    bondDxDy.dx = dx ;
    bondDxDy.dy = dy ;

    % theta in pullback space -- theta_pb    
    % Pack them all into output struct
    outputStruct.fundForms = fundForms ;
    outputStruct.bondDxDy = bondDxDy ;
    outputStruct.theta_pb = theta_pb ;
    
end


