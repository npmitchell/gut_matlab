function plotRelativeMotionTracks(QS, Options)
% plotRelativeMotionTracks(QS, Options)
%   Load tracked nuclei/objects in two layers of the evolving surface and
%   plot their pathlines in 3D on the surface.
%
% Parameters
% ----------
% QS : QuapSlap object
% Options : struct with fields
%   layerLabel : str
%   relMotionFn : str
%       default is fullfile(QS.dir.tracking, 'relative_motion_tracks.mat')
%   lwTrail : float (default=2)
%   uniqueCorrespondence : bool (default = false)
%       plot only one endoderm cell for each muscle cell AND only one
%       muscle cell for each endoderm cell
%   jitterCorrespondence : float >=0 (default = 0) 
%       If more than one match in object A list for a given object B,
%       jitter matching object As by this amount in 3d
%
% Saves to disk
% -------------
% pngs: fullfile(QS.dir.tracking, 'relativeMotionImages',
%                'tracks_%06d.png')
%
%
% NPMitchell 2021



%% Prepare plot settings
if isfield(Options, 'layerLabel')
    subdir = Options.layerLabel ;
    exten = [ '_' subdir ] ;
else
    exten = '' ;
end
opts = load(fullfile(QS.dir.texturePatchIm, ['metadat' exten '.mat'])) ;
xyzlims = opts.metadat.xyzlim ;
smoothing_lambda = opts.metadat.smoothing_lambda ;
normal_shift = opts.metadat.normal_shift ;
backdrop_normal_shift = opts.metadat.normal_shift ;
% Rotation = opts.Options.Rotation ;
% Translation = opts.Options.Translation ;
% Dilation = opts.Options.Dilation ;
flipy = QS.flipy ;
meshFileBase = QS.fullFileBase.mesh ;
initialDistanceThres = Inf ;
overwrite = false ;
lwTrail= 1 ;
uniqueCorrespondence = false ;
jitterCorrespondence = 0 ;

%% Unpack options
relMotionFn = fullfile(QS.dir.tracking, 'relative_motion_tracks.mat') ;
timePoints = QS.xp.fileMeta.timePoints ;
if isfield(Options, 'relMotionFn')  
    relMotionFn = Options.relMotionFn ;
end
if isfield(Options, 'specifyTracks')
    tracks2plot = Options.specifyTracks ;
end
if isfield(Options, 'timePoints')
    timePoints = Options.timePoints ;
end
if isfield(Options, 'normal_shift')
    normal_shift = Options.normal_shift ;
elseif isfield(Options, 'normalShift')
    normal_shift = Options.normalShift ;
end
if isfield(Options, 'initialDistanceThres')
    initialDistanceThres = Options.initialDistanceThres ;
end
if isfield(Options, 'overwrite')  
    overwrite = Options.overwrite ;
end
if isfield(Options, 'lwTrail')  
    lwTrail = Options.lwTrail ;
end
if isfield(Options, 'jitterCorrespondence')
    jitterCorrespondence = Options.jitterCorrespondence ;
end
if isfield(Options, 'uniqueCorrespondence')
    uniqueCorrespondence = Options.uniqueCorrespondence ;
end

% Load relative tracks
load(relMotionFn, 'dusEuclidean', 'dusGeodesic', ...
        'tracks1', 'tracks2', 'pairIDs', 'U0s', 'V0s', ...
        'geodesicPaths', 'ptBarycenters', 'ptFaceLocations', ...
        'euclideanDistanceTraveledU', 'euclideanDistanceTraveledV', ...
        'v3d_u', 'v3d_v', 'nSaved') 
        
nTracks = size(v3d_u, 1) ;
assert(nSaved == nTracks) ;
if ~exist('tracks2plot', 'var')
    tracks2plot = find(dusEuclidean(:, 1) < initialDistanceThres) ;
    if uniqueCorrespondence
        disp('Removing duplicate correspondences in candidate tracks... (uniqueCorrespondence)')
        [~, w] = unique( pairIDs, 'stable' );
        duplicate_indices = setdiff( 1:numel(pairIDs), w ) ;
        duplicate_indices = intersect(tracks2plot, duplicate_indices) ;
        tracks2plot = setdiff(tracks2plot, duplicate_indices) ;
    end
end

% Prepare I/O
if isfield(Options, 'specifyTracks')
    subdir = 'specifiedTracks' ;
elseif isfinite(initialDistanceThres)
    subdir = sprintf('initialDistThres%0.2f', initialDistanceThres) ;
else
    subdir = 'allTracks' ;
end
if jitterCorrespondence > 0
    subdir = [subdir '_jitterCorrespondences'] ;
elseif uniqueCorrespondence > 0
    try
        assert(~isfield(Options, 'specifyTracks'))
    catch
        error('Cannot both specify tracks and demand unique Correspondences, for intelligibility')
    end
    subdir = [subdir '_uniqueCorrespondences'] ;
end
lateralDir = fullfile(QS.dir.tracking, 'relativeMotionImages', subdir, 'lateral1') ;
ventralDir = fullfile(QS.dir.tracking, 'relativeMotionImages', subdir, 'ventral') ;
lateralVentralDir = fullfile(QS.dir.tracking, 'relativeMotionImages', subdir, 'lateralAndVentral') ;
if ~exist(lateralDir, 'dir')
    mkdir(lateralDir)
end
if ~exist(ventralDir, 'dir')
    mkdir(ventralDir)
end
if ~exist(lateralVentralDir, 'dir')
    mkdir(lateralVentralDir)
end

% Title and timing
t0 = QS.t0set() ;

close all
% Pick colors maximally distinguishable from both black and white
trackColors = distinguishable_colors(nTracks, [0,0,0; 1,1,1]) ;
colormap(trackColors)
for tidx = 1:length(timePoints)
    % Set the current time ------------------------------------------------
    tp = timePoints(tidx) ;
    QS.setTime(tp) 
    
    % Obtain mesh for this timepoint as blackdrop (black backdrop)---------
    % Read in the mesh file 
    disp('Reading mesh...')
    % Specfiy the mesh file to load
    meshfn = sprintf( meshFileBase, tp );
    mesh = read_ply_mod( meshfn );

    % If we smooth before pushing along the normal
    if smoothing_lambda > 0 
        disp('smoothing mesh via laplacian filter')
        mesh.v = laplacian_smooth(...
            mesh.v, mesh.f, 'cotan', [], smoothing_lambda, 'implicit') ;

    end

    % Make sure vertex normals are normalized -----------------------------
    mesh.vn = mesh.vn ./ sqrt( sum( mesh.vn.^2, 2 ) );
    mesh0 = mesh ;
    % Normally evolve vertices
    mesh.v = mesh.v + backdrop_normal_shift .* mesh.vn;
    
    % transform into APDV coordsys
    % Allow for overall flip
    % --> apply rotation and translation and dilation BEFORE flipping
    VV = mesh.v ;
    V0 = mesh0.v ;
    if isfield(opts.Options, 'Rotation')
        disp('rotating...')
        VV = (opts.Options.Rotation * VV')' ;
        V0 = (opts.Options.Rotation * V0')' ;
    end
    if isfield(opts.Options, 'Translation')
        disp('translating...')
        VV = VV + opts.Options.Translation ;
        V0 = V0 + opts.Options.Translation ;
    end
    if isfield(opts.Options, 'Dilation')
        disp('dilating...')
        VV = VV * opts.Options.Dilation ;
        V0 = V0 * opts.Options.Dilation ;
        dilation = opts.Options.Dilation ;
    else
        dilation = 1 ;
    end
    if flipy
        VV(:, 2) = -VV(:, 2) ;
        V0(:, 2) = -V0(:, 2) ;
    end
    
    
    % Push muscle layer tracks by normal_shift
    vnormals = normals(V0, mesh.f) ;
    vnormals = vnormals ./ vecnorm(vnormals, 2, 2) ;
    faceIndices = ptFaceLocations(tracks2plot, tidx, 1) ;
    nanID = find(isnan(faceIndices)) ;
    faceIndices(nanID) = 1 ;
    norms2add = vnormals(faceIndices, :) ;
    norms2add(nanID, :) = 0 ;
    
    uus = squeeze(v3d_u(tracks2plot, tidx, :)) ;
    vvs = squeeze(v3d_v(tracks2plot, tidx, :)) + ...
        (normal_shift * dilation) .* norms2add ;
    
    % If jitterCorrespondence is nonzero then jitter any points that have
    % multiple pairIDs
    if jitterCorrespondence > 0 
        [~, w] = unique( pairIDs, 'stable' );
        duplicate_indices = setdiff( 1:numel(pairIDs), w ) ;
        duplicate_indices = intersect(tracks2plot, duplicate_indices) ;
        uind2jitter = ismember(tracks2plot, duplicate_indices) ;
        uus(uind2jitter, :) = uus(uind2jitter, :) + ...
            jitterCorrespondence * ...
            rand(size( uus(uind2jitter, :), 1), size( uus(uind2jitter, :), 2)) ;
    elseif uniqueCorrespondence
        % We should have no duplicates by construction. Assert!
        [~, w] = unique( pairIDs, 'stable' );
        duplicate_indices = setdiff( 1:numel(pairIDs), w ) ;
        duplicate_indices = intersect(tracks2plot, duplicate_indices) ;
        assert(isempty(duplicate_indices))
    end
    
    bcs = barycenter(V0, mesh.f) ;
    
    % check it
    % quiver3(bcs(:, 1), bcs(:, 2), bcs(:, 3), ...
    %     vnormals(:, 1), vnormals(:, 2), vnormals(:, 3), 1)

    
    % Alternative mesh
    % mesh = QS.getCurrentSPCutMeshSmRSC() ;
    
    % Which tracks are initially on the left side (near from view(0,0))?
    nearSide = v3d_u(tracks2plot, 1, 2) < 0;
    
    % Plot the tracks at this time ----------------------------------------
    styles = {'backdrop', 'surf', 'nearOnly'} ;
    for styleID = 1:length(styles)
        exten = styles{styleID} ;
        latDir = fullfile(lateralDir, styles{styleID}) ;
        venDir = fullfile(ventralDir, styles{styleID}) ;
        latVenDir = fullfile(lateralVentralDir, styles{styleID}) ;
        if ~exist(latDir, 'dir')
            mkdir(latDir)
        end
        if ~exist(venDir, 'dir')
            mkdir(venDir)
        end
        if ~exist(latVenDir, 'dir')
            mkdir(latVenDir)
        end
        outputfnBase = fullfile(latVenDir, 'tracks_LateralAndVentral_%s_%06d.png') ;
        
        if ~exist(sprintf(outputfnBase, exten, tp), 'file')   || overwrite

            yzlims = [min(xyzlims(2:3, 1)), max(xyzlims(2:3, 2))] ;
            
            clf; 
            sploth = cell(2, 1) ;
            for subplotID = 1:2
                sploth{subplotID} = subplot(2, 1, subplotID) ;
                hold on;

                % make a backdrop mesh for effect
                if styleID ~= 2
                    trisurf(triangulation(mesh.f, VV + [0, 200, 0]), 'facecolor', 'k')
                    trisurf(triangulation(mesh.f, VV + [0, 0, 200]), 'facecolor', 'k')
                else
                    % make a foreground mesh for opacity/depth perception
                    trisurf(triangulation(mesh.f, V0), 'facecolor', 'k', ...
                        'edgecolor', 'none', 'facealpha', 0.2)
                end

                for ii = 1:length(tracks2plot)
                    if styleID < 3 || nearSide(ii)
                        trackii = tracks2plot(ii) ;
                        scatter3(uus(ii, 1), uus(ii, 2), uus(ii, 3), ...
                            50, 'filled', 'markerfacecolor', ...
                            trackColors(trackii, :), ...
                            'markeredgecolor', 'none') ;
                        scatter3(vvs(ii, 1), vvs(ii, 2), vvs(ii, 3), ...
                            50, '^', 'filled', 'markerfacecolor',...
                            trackColors(trackii, :), ...
                            'markeredgecolor', 'none') ;
                    end
                end

                % Plot trailing track
                if tidx > 1
                    for ii = 1:length(tracks2plot)
                        if styleID < 3 || nearSide(ii)
                            trackid = tracks2plot(ii) ;
                            colorii = trackColors(trackid, :) ;
                            plot3(v3d_u(trackid, 1:tidx, 1), ...
                                v3d_u(trackid, 1:tidx, 2), ...
                                v3d_u(trackid, 1:tidx, 3), '-', ...
                                'color', colorii, 'linewidth', lwTrail)
                            plot3(v3d_v(trackid, 1:tidx, 1), ...
                                v3d_v(trackid, 1:tidx, 2), ...
                                v3d_v(trackid, 1:tidx, 3), '-', ...
                                'color', colorii, 'linewidth', lwTrail)
                        end
                    end
                end
            
                if subplotID == 1
                    % Lateral view
                    view(0, 0)

                    % Format the figure
                    axis equal
                    xlim(xyzlims(1, :))
                    ylim(yzlims + [0,  200])
                    zlim(yzlims)
                    ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
                    zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
                    timeStamp = num2str((timePoints(tidx) - t0) * QS.timeInterval) ;
                    title(['nuclear tracks, $t=$' timeStamp ' ' QS.timeUnits], 'interpreter', 'latex')
                    grid off
                    
                    sh1 = scatter3(NaN, NaN, NaN, ...
                        50, 'filled', 'markerfacecolor', ...
                        'k', ...
                        'markeredgecolor', 'none') ;
                    sh2 = scatter3(NaN, NaN, NaN, ...
                        50, '^', 'filled', 'markerfacecolor','k', ...
                        'markeredgecolor', 'none') ;
                    lh = legend([sh1, sh2], {'endoderm', 'muscle'},...
                        'position', [0.80    0.7182    0.1742    0.0724], ...
                        'interpreter', 'latex') ;
                    % drawnow
                    % axpos = get(sploth, 'position') ;
                else

                    % Ventral view
                    view(0, 270)

                    % Format the figure
                    axis equal
                    xlim(xyzlims(1, :))
                    ylim(yzlims)
                    zlim(yzlims + [0, 200])
                    xlabel(['ap position [' QS.spaceUnits ']'], 'interpreter', 'latex')
                    ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
                    zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
                    grid off
                    % drawnow
                    % ax2pos = get(sploth, 'position') ;
                    % try
                    %     assert(axpos(1) == ax2pos(1))
                    %     assert(axpos(3) == ax2pos(3))
                    %     assert(axpos(4) == ax2pos(4))
                    % catch
                    %     disp('Axes positions unequal, adjusting...')
                    %     set(gca, 'position', [axpos(1),  ax2pos(2), axpos(3), ax2pos(4)])
                    % end
                    % drawnow
                end
                
            end
            
            % Align axes
            axFixPos = get(sploth{1}, 'position') ;
            for subplotIDj = 1:2
                axpos = get(sploth{subplotIDj}, 'position') ;
                set(sploth{subplotIDj}, 'position', ...
                    [axFixPos(1), axpos(2), axFixPos(3), axFixPos(4)])
            end
            
            % Save the figure
            set(gcf,'color','w');
            outfn = sprintf(outputfnBase, exten, tp) ;
            disp(['Saving figure: ' outfn])
            export_fig(outfn, '-nocrop', '-r150')
            
        end
        
        % ventral and lateral separately, for overlay with fluor data.
        outputfnLateral = fullfile(latDir, 'tracks_lateral_%s_%06d.png') ;
        outputfnVentral = fullfile(venDir, 'tracks_ventral_%s_%06d.png') ;
        if ~exist(sprintf(outputfnLateral, exten, tp), 'file') || ...
            ~exist(sprintf(outputfnVentral, exten, tp), 'file') || overwrite

            clf; hold on;

            % make a backdrop mesh for effect
            if styleID ~= 2
                trisurf(triangulation(mesh.f, VV + [0, 200, 0]), 'facecolor', 'k')
                trisurf(triangulation(mesh.f, VV + [0, 0, 200]), 'facecolor', 'k')
            else
                % make a foreground mesh for opacity/depth perception
                trisurf(triangulation(mesh.f, V0), 'facecolor', 'k', ...
                    'edgecolor', 'none', 'facealpha', 0.2)
            end

            for ii = 1:length(tracks2plot)
                if styleID < 3 || nearSide(ii)
                    trackii = tracks2plot(ii) ;
                    scatter3(uus(ii, 1), uus(ii, 2), uus(ii, 3), ...
                        50, 'filled', 'markerfacecolor', ...
                        trackColors(trackii, :), ...
                        'markeredgecolor', 'none')
                    scatter3(vvs(ii, 1), vvs(ii, 2), vvs(ii, 3), ...
                        50, '^', 'filled', 'markerfacecolor', ...
                        trackColors(trackii, :), ...
                        'markeredgecolor', 'none')
                end
            end

            % Plot trailing track
            if tidx > 1
                for ii = 1:length(tracks2plot)
                    if styleID < 3 || nearSide(ii)
                        trackid = tracks2plot(ii) ;
                        colorii = trackColors(trackid, :) ;
                        plot3(v3d_u(trackid, 1:tidx, 1), ...
                            v3d_u(trackid, 1:tidx, 2), ...
                            v3d_u(trackid, 1:tidx, 3), '-', ...
                            'color', colorii, 'linewidth', lwTrail)
                        plot3(v3d_v(trackid, 1:tidx, 1), ...
                            v3d_v(trackid, 1:tidx, 2), ...
                            v3d_v(trackid, 1:tidx, 3), '-', ...
                            'color', colorii, 'linewidth', lwTrail)
                    end
                end
            end
            % Lateral view
            view(0, 0)

            % Format the figure
            set(gcf,'color','w');
            axis equal
            xlim(xyzlims(1, :))
            ylim([xyzlims(2, 1), xyzlims(2, 2) + 200])
            zlim(xyzlims(3, :))
            xlabel(['ap position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            timeStamp = num2str((timePoints(tidx) - t0) * QS.timeInterval) ;
            title(['$t=$' timeStamp ' ' QS.timeUnits], 'interpreter', 'latex')
            grid off
                    
            
            % Save the figure
            outfn = sprintf(outputfnLateral, exten, tp) ;
            disp(['Saving figure: ' outfn])
            export_fig(outfn, '-nocrop', '-r150')

            % Ventral view
            view(0, 270)

            % Format the figure
            axis equal
            xlim(xyzlims(1, :))
            ylim(xyzlims(2, :))
            zlim(xyzlims(3, :) + [0, 200])
            xlabel(['ap position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            timeStamp = num2str((timePoints(tidx) - t0) * QS.timeInterval) ;
            title(['$t=$' timeStamp ' ' QS.timeUnits], 'interpreter', 'latex')
            grid off
            
            % Save the figure
            outfn = sprintf(outputfnVentral, exten, tp) ;
            disp(['Saving figure: ' outfn])
            export_fig(outfn, '-nocrop', '-r150')
            
        end
    end
end
