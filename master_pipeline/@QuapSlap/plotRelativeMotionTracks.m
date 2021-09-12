function plotRelativeMotionTracks(QS, Options)
% plotRelativeMotionTracks(QS, Options)
%
% Parameters
% ----------
% QS : QuapSlap object
% Options : struct with fields
%   layerLabel : str
%   relMotionFn : str
%       default is fullfile(QS.dir.tracking, 'relative_motion_tracks.mat')
% NPMitchell 2021

%%
relMotionFn = fullfile(QS.dir.tracking, 'relative_motion_tracks.mat') ;
if isfield(Options, 'layerLabel')
    subdir = Options.layerLabel ;
    exten = [ '_' subdir ] ;
else
    exten = '' ;
end
if isfield(Options, 'relMotionFn')  
    relMotionFn = Options.relMotionFn ;
end

% Load relative tracks
load(relMotionFn, 'v3d_u', 'v3d_v')

% Prepare plot settings
colormap jetshuffle
opts = load(fullfile(QS.dir.texturePatchIm, ['metadat' exten '.mat'])) ;
xyzlims = opts.metadat.xyzlim ;
smoothing_lambda = opts.metadat.smoothing_lambda ;
normal_shift = opts.metadat.normal_shift ;
% Rotation = opts.Options.Rotation ;
% Translation = opts.Options.Translation ;
% Dilation = opts.Options.Dilation ;
flipy = QS.flipy ;
meshFileBase = QS.fullFileBase.mesh ;
outputfn = fullfile(QS.dir.tracking, 'relativeMotionImages', 'tracks_%06d.png') ; 

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
    % Normally evolve vertices
    mesh.v = mesh.v + normal_shift .* mesh.vn;

    % transform into APDV coordsys
    % Allow for overall flip
    % --> apply rotation and translation and dilation BEFORE flipping
    VV = mesh.v ;
    if isfield(opts.Options, 'Rotation')
        disp('rotating...')
        VV = (Options.Rotation * VV')' ;
    end
    if isfield(opts.Options, 'Translation')
        disp('translating...')
        VV = VV + opts.Options.Translation ;
    end
    if isfield(opts.Options, 'Dilation')
        disp('dilating...')
        VV = VV * opts.Options.Dilation ;
    end
    if flipy
        VV(:, 2) = -VV(:, 2) ;
    end
    
    % Alternative mesh
    % mesh = QS.getCurrentSPCutMeshSmRSC() ;
    
    % Plot the tracks at this time ----------------------------------------
    hold off ;
    trisurf(triangulation(mesh.f, VV + [0, 200, 0]), 'facecolor', 'k')
    hold on;
    trisurf(triangulation(mesh.f, VV), 'facecolor', 'k', ...
        'edgecolor', 'none', 'facealpha', 0.2)
    scatter3(v3d_u(:, tidx, 1), v3d_u(:, tidx, 2), v3d_u(:, tidx, 3), ...
        50, (1:nTracks), 'filled', 'markeredgecolor', 'none')
    scatter3(v3d_v(:, tidx, 1), v3d_v(:, tidx, 2), v3d_v(:, tidx, 3), ...
        50, 1:nTracks, 's')
    
    % Plot trailing track
    if tidx > 1
        plot3(v3d_u(:, 1:tidx, 1), v3d_u(:, 1:tidx, 2), v3d_u(:, 1:tidx, 3), '-')
        plot3(v3d_v(:, 1:tidx, 1), v3d_v(:, 1:tidx, 2), v3d_v(:, 1:tidx, 3), '-')
    end
    
    % Lateral view
    view(0, 0)
    
    % Format the figure
    axis equal
    xlim(xyzlims(1, :))
    ylim(xyzlims(2, :))
    zlim(xyzlims(3, :))
    
    export_fig(sprintf(outputfn, tp), '-nocrop', '-r200')
end
