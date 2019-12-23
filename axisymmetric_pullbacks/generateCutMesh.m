function [cutMesh, adIDx, pdIDx, cutP, success] = generateCutMeshFixedTwist(mesh, adIDx, pdIDx, centerline, nsegs4path, previousP, prevTw)
% GENERATECUTMESH(mesh, adIDx, pdIDx, centerline, nsegs4path)
%    Generating cutMesh from a mesh and the start/endpoints a/pdIDx
% 
% Parameters
% -----------
% mesh
% adIDx 
% pdIDx
% centerline
% nsegs4path
%
% Returns
% -------
% prevTw : 
% previousP
% compute_pullback = success
%
% NPMitchell 2019 


success = true ;  % mark as successful unless otherwise caught

% First try geodesic approach
cutOptions.method = 'fastest' ;
disp(['Cutting mesh using method ' cutOptions.method])

try
    [ cutMesh, adIDx, pdIDx, cutP ] = ...
        cylinderCutMesh( mesh.f, mesh.v, mesh.vn, adIDx, pdIDx, cutOptions );

    % If twist about centerline has jumped, find nearest piecewise geodesic
    Tw = twist(mesh.v(cutP, :), centerline) ;
    disp(['Found Twist = ' num2str(Tw)])

    % Initialize some book-keeping variables
    trykk = 0 ;
    nsegs4path_kk = nsegs4path ;
    while abs(Tw - prevTw) > 0.5
        disp('Twist out of range. Forcing nearest curve.')
        % If this is the first timepoint we are considering
        % or if we have to iteratively refine the target 
        % curve, load previous cutP.
        if new_run || trykk > 0
            % Load previous mesh and previous cutP
            prevcylmeshfn = sprintf( cylinderMeshCleanBase, t-trykk-1) ;
            prevmesh = read_ply_mod( prevcylmeshfn ); 
            prevcutP = dlmread(sprintf(outcutfn, t-1-trykk), ',', 1, 0) ;
            previousP = prevmesh.v(prevcutP, :) ;
            % subsample the path to get previousP
            % Note: avoid oversampling by taking min, with
            % fudge factor to avoid perfect sampling (which
            % is too dense also).
            pstep = min(round(length(prevcutP) / nsegs4path ), round(length(prevcutP) * 0.7)) ;
            previousP = previousP(1:pstep:end, :) ;
        end

        cutOptions.method = 'nearest' ;
        cutOptions.path = previousP ;
        disp(['Cutting mesh using method ' cutOptions.method])
        [ cutMesh, adIDx, pdIDx, cutP ] = ...
            cylinderCutMesh( mesh.f, mesh.v, mesh.vn, adIDx, pdIDx, cutOptions );                
        Tw = twist(mesh.v(cutP, :), centerline) ;
        disp(['Found new Twist = ' num2str(Tw)])
        if abs(Tw - prevTw) > 0.5
            nsegs4path_kk = nsegs4path_kk * 2 ;
            trykk = trykk + 1 ;
        end
    end
    disp('Twist within range. Continuing...')
    prevTw = Tw ;
catch
    disp('Could not cut this timepoint: Input mesh probably NOT a cylinder')
    success = false ;
end

