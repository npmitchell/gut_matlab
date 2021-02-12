function plotSeriesOnSurfaceTexturePatch(QS,...
   options, TexturePatchOptions)
% PLOTDATAONSURFACETEXTUREPATCH(xp, xyzlim, rot, trans)
%   Plot intensity data timeseries on 3d mesh timeseries
%
% Parameters
% ----------
% QS : QuapSlap class instance
% overwrite : bool
% metadat : struct with fields
% TexturePatchOptions: struct with fields
%
%
% Returns
% -------
% <none>
%
% Saves to disk
% -------------
% metadat.mat : mat file with variables 'metadat' and 'Options'
%   options for how to plot and compute texturepatches
%
% NPMitchell 2020

%% Default options
overwrite = false ;
plot_dorsal = true ;
plot_ventral = true ;
plot_left = true ;
plot_right = true ;
plot_perspective = true ;
channel = [] ;  % by default, plot all channels

%% Unpack Options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'channel')
    channel = options.channel ;
end
if isfield(options, 'plot_dorsal')
    plot_dorsal = options.plot_dorsal ;
end
if isfield(options, 'plot_ventral')
    plot_ventral = options.plot_ventral ;
end
if isfield(options, 'plot_dorsal')
    plot_left = options.plot_left ;
end
if isfield(options, 'plot_dorsal')
    plot_right = options.plot_right ;
end
if isfield(options, 'plot_perspective')
    plot_perspective = options.plot_perspective ;
end
% Collate boolean plot indicators to decide which views to plot
plot_view = [plot_dorsal, plot_ventral, plot_left, ...
    plot_right, plot_perspective] ;

%% Unpack QS
meshFileBase = QS.fullFileBase.mesh ;
figoutdir = QS.dir.texturePatchIm ;
normal_shift = QS.normalShift ;
flipy = QS.flipy ;
% texture_axis_order = QS.data.axisOrder ;

try
    t0 = QS.t0set() ;
catch
    t0 = QS.xp.fileMeta.timePoints(1) ;
end

%% Load metadat and TexturePatchOptions if not supplied
metafn = fullfile(figoutdir, 'metadat.mat') ;
resave_metadat = false ;
if nargin < 4
    try
        load(metafn, 'Options')
    catch
        [rot, trans] = QS.getRotTrans() ;
        resolution = QS.APDV.resolution ;
        
        % Psize is the linear dimension of the grid drawn on each triangular face
        Options.PSize = 5;
        Options.EdgeColor = 'none';
        Options.Rotation = rot ;
        Options.Translation = trans ;
        Options.Dilation = resolution ;
        Options.numLayers = [1, 1];  % [2, 2] marches ~0.5 um in either dir
        Options.layerSpacing = 2 ;
        resave_metadat = true ;
    end
else
    Options = TexturePatchOptions ;
end

if nargin < 3
    try
        load(metafn, 'metadat')
    catch
        % Define & Save metadata
        [~, ~, ~, ~, xyzbuff] = QS.getXYZLims() ;
        xyzbuff(:, 1) = xyzbuff(:, 1) - 10 ;
        xyzbuff(:, 2) = xyzbuff(:, 2) + 10 ;
        metadat.normal_shift = QS.normalShift ;             % normal push, in pixels, along normals defined in data XYZ space
        metadat.xyzlim = xyzbuff ;                          % xyzlimits
        metadat.texture_axis_order = QS.data.axisOrder ;    % texture space sampling
        metadat.reorient_faces = false ;                    % set to true if some normals may be inverted
        resave_metadat = true ;
    end
else
    metadat = options ;
end

% Save it
if resave_metadat
    save(metafn, 'metadat', 'Options')
end

% pass metadat to options
options.normal_shift = metadat.normal_shift ;
options.xyzlim = metadat.xyzlim ;
options.texture_axis_order = metadat.texture_axis_order ;
texture_axis_order = options.texture_axis_order ;
options.reorient_faces = metadat.reorient_faces ;

%% Name output directories
figdDir = fullfile(figoutdir, 'dorsal') ;
figvDir = fullfile(figoutdir, 'ventral') ;
figlat1Dir = fullfile(figoutdir, 'lateral1') ;
figlat2Dir = fullfile(figoutdir, 'lateral2') ;
figPerspDir = fullfile(figoutdir, 'perspective') ;
dirs = {figoutdir, figdDir, figvDir, figlat1Dir, figlat2Dir, figPerspDir} ;
for i = 1:length(dirs)
    if ~exist(dirs{i}, 'dir')
        mkdir(dirs{i}) ;
    end
end

%% Define output filenames
fns = {fullfile(figdDir, 'patch_dorsal_%06d.png'), ...
    fullfile(figvDir, 'patch_ventral_%06d.png'), ...
    fullfile(figlat1Dir, 'patch_lateral1_%06d.png'), ...
    fullfile(figlat2Dir, 'patch_lateral2_%06d.png'), ...
    fullfile(figPerspDir, 'patch_persp_%06d.png') };

% Unpack metadat and save 
xyzlim = options.xyzlim ;
reorient_faces = options.reorient_faces ;
timeinterval = QS.timeInterval ;
timeunits = QS.timeUnits ;
timePoints = QS.xp.fileMeta.timePoints ;

% Unpack xyzlim
xmin = xyzlim(1, 1); xmax = xyzlim(1, 2) ;
ymin = xyzlim(2, 1); ymax = xyzlim(2, 2) ;
zmin = xyzlim(3, 1); zmax = xyzlim(3, 2) ;

%% Now draw for all TPs
% figure parameters
xwidth = 16 ; % cm
ywidth = 10 ; % cm

tidx_todoA = 1:20:length(timePoints) ;
tidx_todoB = setdiff(1:10:length(timePoints), tidx_todoA) ;
tidx_todo = [tidx_todoA, tidx_todoB] ;
tidx_todoC = setdiff(1:length(timePoints), tidx_todo) ;
tidx_todo = [tidx_todo, tidx_todoC] ;

for tidx = tidx_todo
    tp = timePoints(tidx) ;
    ondisk = true ;
    for ii = 1:length(fns)
        ondisk = ondisk && exist(sprintf(fns{ii}, tp), 'file') ;
    end
    
    if (overwrite || ~ondisk) && any(plot_view)
        tic 
        close all
        % Copy passed Options argument for unpacking
        Options = TexturePatchOptions ; 

        QS.setTime(tp)
        QS.getCurrentData()
        IV = QS.currentData.IV ;
        
        % % Get the data in 3d for this timepoint
        % tidx = QS.xp.tIdx(tp) ;
        % % if t ~= xp.currentTime
        % xp.loadTime(tp) ;
        % xp.rescaleStackToUnitAspect() ;
        % % end
        % 
        % disp('Applying image...')
        % IV = xp.stack.image.apply();
        % IV = QS.adjustIV(IV, adjustlow, adjusthigh) ;

        % Read in the mesh file -----------------------------------------------
        disp('Reading mesh...')
        % Specfiy the mesh file to load
        meshfn = sprintf( meshFileBase, tp );
        mesh = read_ply_mod( meshfn );

        % Make sure vertex normals are normalized
        mesh.vn = mesh.vn ./ sqrt( sum( mesh.vn.^2, 2 ) );
        % Normally evolve vertices
        mesh.v = mesh.v + normal_shift .* mesh.vn;
        % Re-orient faces
        if reorient_faces
            disp('reorienting faces...')
            mesh.f = reorient_facets( mesh.v, mesh.f );
        end

        % First args are physical vertices, then texture faces (same number as 
        % physical faces, but connectivity can be different), then texture
        % vertices, which can also be different. The vertices are in row, column, 
        % page format, so x and y get flipped. IV is the texture volume.
        % Options.PSize 
        % mesh.f = f{tidx} ;
        % mesh.v = v{tidx} ;
        % mesh.vn = vn{tidx} ;

        fig = figure('Visible', 'Off') ;
        % set(gcf, 'Visible', 'Off') 
        disp(['creating texture patch ' num2str(tp, '%06d')])

        % Allow for axis swapping
        disp(['texture axis order: [', ...
            num2str(texture_axis_order(1)), ...
            ' ', num2str(texture_axis_order(2)), ...
            ' ', num2str(texture_axis_order(3)) ']'])
        TV = mesh.v(:, texture_axis_order) ;
        
        
        % Allow for overall flip
        % --> apply rotation and translation and dilation BEFORE flipping
        VV = mesh.v ;
        if isfield(Options, 'Rotation')
            disp('rotating...')
            VV = (Options.Rotation * VV')' ;
            Options = rmfield(Options, 'Rotation') ;
        end
        if isfield(Options, 'Translation')
            disp('translating...')
            VV = VV + Options.Translation ;
            Options = rmfield(Options, 'Translation') ;
        end
        if isfield(Options, 'Dilation')
            disp('dilating...')
            VV = VV * Options.Dilation ;
            Options = rmfield(Options, 'Dilation') ;
        end
        if flipy
            VV(:, 2) = -VV(:, 2) ;
        end

        % Create the texture patch
        if ~isempty(channel)
            IV2plot = cell(1) ;
            IV2plot{1} = IV{channel} ;
            texture_patch_3d( mesh.f, VV, mesh.f, TV, IV2plot, Options );
        else
            texture_patch_3d( mesh.f, VV, mesh.f, TV, IV, Options );
        end
        
        % format the figure
        disp('formatting figure...')
        axis equal
        xlim([xmin, xmax])
        ylim([ymin, ymax])
        zlim([zmin, zmax])
        colormap bone
        titlestr = ['$t = $' num2str(tp*timeinterval-t0) ' ' timeunits] ;
        title(titlestr, 'Interpreter', 'Latex', 'Color', 'white') 
        xlabel('AP position [$\mu$m]', 'Interpreter', 'Latex')
        ylabel('lateral position [$\mu$m]', 'Interpreter', 'Latex')
        zlabel('DV position [$\mu$m]', 'Interpreter', 'Latex')

        % Rotate the camera angle using rotation and translation 
        % camorbit(theta, phi)
        % camup() --> ax.CameraUpVector = [sin(45) cos(45) 1]
        % camroll(dtheta)
        % This is an active rotation
        % rotate(hSurface,direction,25)
        set(fig, 'PaperUnits', 'centimeters');
        set(fig, 'PaperPosition', [0 0 xwidth ywidth]);

        % Make background black & Make tick labels white
        set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w', 'zcol', 'w')
        set(gcf, 'InvertHardCopy', 'off');
        set(gcf, 'Color', 'k')
        set(gcf, 'color', 'k')

        % Capture all four views
        disp(['saving figure...' num2str(tp, '%06d')])
        % Save each figure
        for ii = 1:length(fns)
            % Only plot this view if plot_view(ii) is true
            if plot_view(ii)
                if ii == 1
                    % dorsal
                    disp('saving dorsal image...')
                    view(0, 90)
                elseif ii == 2
                    % ventral
                    disp('saving ventral image...')
                    view(0, 270)
                elseif ii == 3
                    % Lateral views
                    disp('saving lateral image...')
                    view(0, 0)
                elseif ii == 4
                    % lateral view 2
                    disp('saving second lateral image...')
                    view(0, 180)
                elseif ii == 5
                    % perspective view
                    disp('saving perspective image...')
                    view(-20, 20)
                else
                    error(['Exhausted DorsalVentralLeftRight indices. ',...
                        'What is going on here?'])
                end

                % Use export_fig instead, from plotting/export_fig/
                % saveas(fig, fullfile(figvdir, fnv))
                export_fig(sprintf(fns{ii}, tp), '-nocrop', '-r200')
            end
        end
        close all
        toc
    end
end

clear Options IV