function generateCurrentCutMesh(QS)
%

% Unpack parameters
tt = QS.currentTime ;
cutMeshfn = sprintf(QS.fullFileBase.cutMesh, tt) ;
QS.getCleanCntrlines() ;
mesh = QS.currentMesh.cylinderMesh ;
if isempty(mesh)
    QS.loadCurrentCylinderMesh() ;
    mesh = QS.currentMesh.cylinderMesh ;
end

% Grab ad and pd indices for cylinder mesh
adIDx = h5read(QS.fileName.aBoundaryDorsalPtsClean,...
    ['/' sprintf('%06d', tt)]) ;
pdIDx = h5read(QS.fileName.pBoundaryDorsalPtsClean,...
    ['/' sprintf('%06d', tt)]) ;

% Output names
outcutfn = sprintf(QS.fullFileBase.cutPath, tt) ;

% try geodesic if first timepoint
if tt == QS.xp.fileMeta.timePoints(1)
    cutOptions.method = 'fastest' ;
    disp(['Cutting mesh using method ' cutOptions.method])
    cutMesh = cylinderCutMesh( mesh.f, mesh.v, mesh.vn, adIDx, pdIDx, cutOptions );
    cutP = cutMesh.pathPairs(:, 1) ;
    adIDx = cutP(1) ;
    pdIDx = cutP(end) ;
else 
    % If a previous Twist is not held in RAM, compute it
    % if ~exist('prevTw', 'var')
    % Load previous mesh and previous cutP
    prevcylmeshfn = sprintf( cylinderMeshCleanBase, tt-1) ;
    prevmesh = read_ply_mod( prevcylmeshfn ); 
    prevcutP = dlmread(sprintf(outcutfn, tt-1), ',', 1, 0) ;
    previousP = prevmesh.v(prevcutP, :) ;
    % Load previous centerline in raw units
    prevcline = cntrlines{tt-1} ; % use previous CORRECTED centerline (non-anomalous)
    % Compute Twist for this previous timepoint
    prevTw = twist(previousP, prevcline) ;

    % [edgelen, annulusv2d] = DiscreteRicciFlow.EuclideanRicciFlow(mesh.f, mesh.v, ...
    %     'BoundaryType', 'fixed', 'BoundaryShape', 'Circles', ...
    %     'MaxIter', 25, 'MaxCircIter', 21, ...
    %     'Tolerance', 1e-6, 'CircTolerance', 1e-4, 'PCGTolerance', 1e-4) ;
    % findAnnularPathZeroWindingNumber(mesh.f, annulusv2d, adIDx, pdIDx)

    % Which path to match this one to: choose previous timepoint
    % Load previous mesh and previous cutP
    prevcylmeshfn = sprintf( cylinderMeshCleanBase, tt-1) ;
    prevmesh = read_ply_mod( prevcylmeshfn ); 
    prevcutP = dlmread(sprintf(outcutfn, tt-1), ',', 1, 0) ;
    previousP = prevmesh.v(prevcutP, :) ;

    [cutMesh, adIDx, pdIDx, cutP, ~] = ...
        generateCutMeshFixedTwist(mesh, adIDx, pdIDx, ...
        cntrlines{tt},...  % supply the current corrected centerline
        nsegs4path, prevTw, previousP, ...
        'MaxTwChange', maxTwChange, 'MaxJitter', maxJitter, ...
        'PrevCntrline', prevcline) ;
end

% Store this path for the next one to be nearby
% Limit the number of segments to nsegs4path
% previousP = cutMesh.v(cutP, :) ;
% pstep = round(length(cutP) / nsegs4path ) ;
% previousP = previousP(1:pstep:end, :) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Done with generating initial 3D CutMesh with cutPath\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the cutPath to txt file
header = 'cutP (path of cut), indexing into vertices' ;
write_txt_with_header(sprintf(outcutfn, tt), cutP, header)  
clearvars header

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate pullback to rectangular domain ---------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The surface parameterization algorithm (optionally) takes four vertex IDs
% as input to specify the corners of the square parameterization domain.
% Maddeningly, the order in which these points are specified to not seem to
% effect the output. For consistency, we perform a post-hoc correction so
% that the final output has the following geometric ordering
%
%   (AD1)-------(PD1)
%     |           |
%     |           |
%     |           |
%     |           |
%   (AD2)-------(PD2)
%
% Note that the pathPairs variable has the following columns:
%   [ ( AD1 -> PD1 ), ( AD2 -> PD2 ) ]
%--------------------------------------------------------------------------

% View results --------------------------------------------------------
% P = cutMesh.pathPairs(:,1);
% 
% trisurf( triangulation( mesh.f, mesh.v ) );
% 
% hold on
%
% line( mesh.v(P,1), mesh.v(P,2), mesh.v(P,3), ...
%     'Color', 'c', 'LineWidth',2);
% 
% scatter3( mesh.v(adIDx,1), mesh.v(adIDx,2), mesh.v(adIDx,3), ...
%     'filled', 'r' );
% scatter3( mesh.v(pdIDx,1), mesh.v(pdIDx,2), mesh.v(pdIDx,3), ...
%     'filled', 'm' );
% 
% hold off
% 
% axis equal
% 
% clear P

%----------------------------------------------------------------------
% Generate Pullback to Annular Orbifold Domain
%----------------------------------------------------------------------
fprintf('Relaxing network via Affine transformation... ');
cutMesh = flattenAnnulus( cutMesh );

% Find lateral scaling that minimizes spring network energy
ar = minimizeIsoarealAffineEnergy( cutMesh.f, cutMesh.v, cutMesh.u );
% Assign scaling based on options: either a0 or a_fixed
% if tidx == 1 && ~a_fixed
%     a_fixed = ar ;
% end      
% a = a_fixed ;

% Scale the x axis by a or ar
uvtx = cutMesh.u ;
cutMesh.u = [ QS.a_fixed .* uvtx(:,1), uvtx(:,2) ];
cutMesh.ar = ar ;
cutMesh.umax = QS.a_fixed ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Done flattening cutMesh. Now saving.\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save cutMesh
disp(['Saving cutMesh to ' cutMeshfn]) 
save(cutMeshfn, 'cutMesh', 'adIDx', 'pdIDx', 'cutP')

% Displacing mesh along normal direction
disp(['Evolving mesh along normal shift for pullback ',...
    'images: shift=' num2str(QS.normalShift)])
cutMesh.v = cutMesh.v + cutMesh.vn * QS.normalShift ;

disp('Assigning current cutMesh to self')
QS.currentMesh.cutMesh = cutMesh ;
disp('done with plotting & saving cut')