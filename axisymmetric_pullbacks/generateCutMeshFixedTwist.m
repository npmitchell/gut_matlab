function [cutMesh, adIDx, pdIDx, cutP, Tw] = generateCutMeshFixedTwist(mesh, ...
    adIDx, pdIDx, centerline, nsegs4path, prevTw, previousP, varargin)
% GENERATECUTMESH(mesh, adIDx, pdIDx, centerline, nsegs4path)
%    Generating cutMesh from a mesh and the start/endpoints a/pdIDx
% 
% Parameters
% -----------
% mesh
% adIDx 
% pdIDx
% centerline : N x 3 float
%   approx centerline of the mesh in original pixel coordinates
% nsegs4path : int
%   how many segments to divide the path into initially to compute a nearby
%   cut
% prevTw : float
%   The twist value or winding number of previous mesh
% outcutfn : str
%   The output Base filename for the cutPath (evaluatable string with
%   sprintf)
% cylinderMeshCleanBase : str
%   The output Base filename for the cylinderMesh 
% t : int
%   time index of current mesh, whos path is to be compared to previous
%
% Optional Input Arguments
% ------------------------
% MaxTwChange : float
%   Maximum twist change allowed in the path
% MaxJitter : float 
%
% Returns
% -------
% prevTw : 
% previousP
% compute_pullback = success
% prevTw = Tw : the current Twist value, to be stored for next time
%
% NPMitchell 2019 


% success = true ;  % mark as successful unless otherwise caught

max_Tw_change = 0.25 ;
max_jitter = 10 ;
max_nsegs4path = 10 ;
jitter_growth_rate = 10 ;  %linear growth rate
nsegs_growth_rate = 2 ; % linear growth rate
for i = 1:length(varargin)
    if isa(varargin{i},'double') 
        continue;
    end
    if isa(varargin{i},'logical')
        continue;
    end
    
    if ~isempty(regexp(varargin{i},'^[Mm]ax[Tt]w[Cc]hange','match'))
        max_Tw_change = varargin{i+1} ;
    end
    if ~isempty(regexp(varargin{i},'^[Mm]ax[Jj]itter','match'))
        max_jitter = varargin{i+1} ;
    end
    if ~isempty(regexp(varargin{i},'^[Pp]rev[Cc]enter[Ll]ine','match'))
        prevcntrline = varargin{i+1} ;
    elseif ~isempty(regexp(varargin{i},'^[Pp]rev[Cc]ntr[Ll]ine','match'))
        prevcntrline = varargin{i+1} ;
    else
        prevcntrline = [] ;
    end
end


% First try geodesic approach
cutOptions.method = 'fastest' ;
disp(['Cutting mesh using method ' cutOptions.method])

% try
    cutMesh = cylinderCutMesh( mesh.f, mesh.v, mesh.vn, adIDx, pdIDx, cutOptions );
    cutP = cutMesh.pathPairs(:, 1) ;
    adIDx = cutP(1) ;
    pdIDx = cutP(end) ;

    % If twist about centerline has jumped, find nearest piecewise geodesic
    Tw = twist(mesh.v(cutP, :), centerline) ;
    disp(['Found Twist = ' num2str(Tw) '; previous twist = ' num2str(prevTw)])

    % Initialize some book-keeping variables
    trykk = 0 ;
    jitter_amp = 0 ;
    nsegs4path_kk = nsegs4path ;
    while abs(Tw - prevTw) > max_Tw_change
        disp('Twist out of range. Forcing nearest curve.')
        % If this is the first timepoint we are considering
        % or if we have to iteratively refine the target 
        % curve, load previous cutP.
        
        %%%%%%%%%%%%%%%%%%%
        % subsample the previous path to get previousP_kk
        % Note: avoid oversampling by taking min, with
        % fudge factor to avoid perfect sampling (which
        % is too dense also).
        disp(['nsegs4path = ' num2str(nsegs4path_kk)])
        pstep = round(length(previousP) / min(max_nsegs4path, nsegs4path_kk)) ;
        
        disp(['pstep = ' num2str(pstep)])
        disp(['jitter_amp = ', num2str(min(jitter_amp, max_jitter))])
        previousP_kk = previousP(1:pstep:end, :);
        jitter = min(jitter_amp, max_jitter) * (rand(size(previousP_kk)) - 0.5);
        previousP_kk = previousP_kk + jitter ;
        %%%%%%%%%%%%%%%%%%%%%

        cutOptions.method = 'nearest' ;
        cutOptions.path = previousP_kk;
        disp(['Cutting mesh using method ' cutOptions.method])
        try
            cutMesh = cylinderCutMesh( mesh.f, mesh.v, mesh.vn, adIDx, pdIDx, cutOptions );  
            cutP = cutMesh.pathPairs(:,1) ;
            Tw = twist(mesh.v(cutP, :), centerline) ;
            disp(['Found new Twist = ' num2str(Tw)])
            if abs(Tw - prevTw) > max_Tw_change
                nsegs4path_kk = round(nsegs4path_kk + nsegs_growth_rate) ;
                jitter_amp = jitter_amp + jitter_growth_rate ;
                trykk = trykk + 1 ;

                % show progress
                figure(1);
                set(gcf, 'visible', 'on')
                trisurf(cutMesh.f, cutMesh.v(:, 1), cutMesh.v(:, 2), ...
                    cutMesh.v(:, 3), cutMesh.v(:, 3), ...
                    'EdgeColor', 'none', 'FaceAlpha', 0.5)
                hold on;
                plot3(cutMesh.v(cutP, 1), cutMesh.v(cutP, 2), cutMesh.v(cutP, 3), 'k-')
                plot3(previousP_kk(:, 1), previousP_kk(:, 2), previousP_kk(:, 3), 'o-')
                plot3(centerline(:, 1), centerline(:, 2), centerline(:, 3), '-')
                if ~isempty(prevcntrline)
                    plot3(prevcntrline(:, 1), prevcntrline(:, 2), prevcntrline(:, 3), '--')
                end
                title(['sampling = ' num2str(length(previousP_kk))])
                hold off 
                axis equal
                pause(0.0000000000001)
            end
        catch
            disp('Could not generate cutMesh, likely not a topological disk')
            pause(0.0001)
            close all
            nsegs4path_kk = round(nsegs4path_kk + nsegs_growth_rate) ;
            jitter_amp = jitter_amp + jitter_growth_rate ;
            trykk = trykk + 1 ;
        end
        
        % Give up on perturbing the curve after 10 tries, instead construct
        % the solution by projecting onto an annulus
        % if trykk > 10
        %     % Glue the previous mesh back together again
        %     % prevcmesh = closeRectilinearCylMesh(prevmesh) ;
        %     % Glue the current mesh back together again
        %     % cmesh = closeRectilinearCylMesh(mesh) ;
        %     wN = annularPathWindingNumber(cmesh.f, cmesh.v, cutMesh.pathPairs) ;
        % end
    end
    disp('Twist within range. Continuing...')
    error('break')
% catch
%     disp('Could not cut this timepoint: Input mesh probably NOT a cylinder')
%     success = false ;
% end

% Redefine the endpoints
adIDx = cutP(1) ;
pdIDx = cutP(end) ;
