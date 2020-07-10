function plotSeriesOnSurfaceTexturePatch(QS,...
    overwrite, metadat, TexturePatchOptions)
% PLOTDATAONSURFACETEXTUREPATCH(xp, xyzlim, rot, trans)
%   Plot intensity data timeseries on 3d mesh timeseries
%
% Parameters
% ----------
% QS :
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

%% Unpack QS
meshFileBase = QS.fullFileBase.mesh ;
figoutdir = QS.dir.texturePatchIm ;
normal_shift = QS.normalShift ;
flipy = QS.flipy ;
texture_axis_order = QS.data.axisOrder ;

try
    t0 = QS.t0set() ;
catch
    t0 = QS.xp.fileMeta.timePoints(1) ;
end

% Load metadat and TexturePatchOptions if not supplied
metafn = fullfile(figoutdir, 'metadat.mat') ;
resave_metadat = false ;
if nargin < 4
    try
        load(metafn, 'Options')
    catch
        % Psize is the linear dimension of the grid drawn on each triangular face
        Options.PSize = 5;
        Options.EdgeColor = 'none';
        Options.Rotation = rot ;
        Options.Translation = trans ;
        Options.Dilation = resolution ;
        Options.numLayers = [1, -1];  % at layerSpacing 2, 2 marches ~0.5 um 
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
        metadat.normal_shift = 10 ;                 % pixels
        metadat.xyzlim = xyzbuff ;                  % xyzlimits
        metadat.texture_axis_order = [1, 2, 3] ;    % texture space sampling
        metadat.reorient_faces = false ;            % if some normals are inverted
        resave_metadat = true ;
    end
end

% Save it
if resave_metadat
    save(metafn, 'metadat', 'Options')
end

%% Name output directories
figddir = fullfile(figoutdir, 'dorsal') ;
figvdir = fullfile(figoutdir, 'ventral') ;
figlat1dir = fullfile(figoutdir, 'lateral1') ;
figlat2dir = fullfile(figoutdir, 'lateral2') ;
dirs = {figoutdir, figddir, figvdir, figlat1dir, figlat2dir} ;
for i = 1:length(dirs)
    if ~exist(dirs{i}, 'dir')
        mkdir(dirs{i}) ;
    end
end

%% Define output filenames
fns = {fullfile(figddir, 'patch_dorsal_%06d.png'), ...
    fullfile(figvdir, 'patch_ventral_%06d.png'), ...
    fullfile(figlat1dir, 'patch_lateral1_%06d.png'), ...
    fullfile(figlat2dir, 'patch_lateral2_%06d.png') };

% Unpack metadat and save 
xyzlim = metadat.xyzlim ;
reorient_faces = metadat.reorient_faces ;
timeinterval = QS.timeinterval ;
timeunits = QS.timeunits ;

% Unpack xyzlim
xmin = xyzlim(1, 1); xmax = xyzlim(1, 2) ;
ymin = xyzlim(2, 1); ymax = xyzlim(2, 2) ;
zmin = xyzlim(3, 1); zmax = xyzlim(3, 2) ;

%% Now draw for all TPs
% figure parameters
xwidth = 16 ; % cm
ywidth = 10 ; % cm

for tp = QS.xp.fileMeta.timePoints
    ondisk = true ;
    for ii = 1:length(fns)
        ondisk = ondisk && exist(sprintf(fns{ii}, tp), 'file') ;
    end
    
    if overwrite || ~ondisk
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
        texture_patch_3d( mesh.f, VV, mesh.f, TV, IV, Options );

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
            if ii == 1
                % dorsal
                view(0, 90)
            elseif ii == 2
                % ventral
                view(0, 270)
            elseif ii == 3
                % Lateral views
                view(0, 0)
            elseif ii == 4
                % lateral view 2
                view(0, 180)
            else
                error('Exhausted DVLR indices. What is going on here?')
            end

            % Use export_fig instead, from plotting/export_fig/
            % saveas(fig, fullfile(figvdir, fnv))
            export_fig(sprintf(fns{ii}, tp), '-nocrop', '-r200')
        end
        close all
        toc
    end
end

clear Options IV