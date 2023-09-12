% For Figure 4 panel a in Tubular paper

%% ************************************************************************
% *************************************************************************
%               GENERATE TEXTURES SURFACE IMAGES
% *************************************************************************
% *************************************************************************

%% Generate Textured Surface Images =======================================
close all; clc;

% Set metadata options
overwrite = false ;
texture_shift = 0 ;
reorient_faces = false ;
% plot_dorsal = true ;
% plot_ventral = true ;
% plot_left = true ;
% plot_right = true ;
% plot_perspective = true ;
blackFigure = false ;
makeColorbar = false ;
smoothing_lambda = 0.0 ;
channel = [] ;  % by default, plot all channels
% figoutdir = tubi.dir.texturePatchIm ;
normal_shift = 0 ;
% plot_time_points = [];
texture_axis_order = [1 2 3] ;
use_sphere_mesh = false;

% Choose view
viewType = 'dorsal';

% Set texture patch options
Options.PSize = 5;          % Psize is the linear dimension of the grid drawn on each triangular face. Set PSize > 1 for refinement of texture on each triangle of the surface triangulation. Higher numbers are slower but give more detailed images.
Options.numLayers = [1, 1];  % how many layers to MIP over/bundle into stack, as [outward, inward]
Options.layerSpacing = 2;   % Distance between layers over which we take MIP, in pixels, 
Options.edgeColor = 'none';
Options.faceLighting = 'gouraud';
Options.SpecularStrength = 0.5;
Options.DiffuseStrength = 0.9; % 0.2;
Options.AmbientStrength = 0.9; % 0.2;
Options.SpecularExponent = 0.5;
Options.BackFaceLighting = 'lit';

% figure parameters
xwidth = 16 ; % cm
ywidth = 10 ; % cm

timePoints = tubi.xp.fileMeta.timePoints ;

% Choose the time point to be visualized
% tp = 1; original
% tp = 4; original 
tp = 5 ;
% tp = 8; original

for tp = 1:30
    tubi.setTime(tp)
    
    try
        t0 = tubi.t0set() ;
    catch
        t0 = tubi.xp.fileMeta.timePoints(1) ;
    end
    
    % Unpack metadat and save 
    timeinterval = tubi.timeInterval ;
    timeunits = tubi.timeUnits ;
    
    % Unpack xyzlim
    [~, ~, ~, xyzlim] = tubi.getXYZLims() ;
        metadat.xyzlim = xyzlim ;
    buff = 4 ;
    xmin = xyzlim(1, 1) - buff; 
    xmax = xyzlim(1, 2) + buff;
    ymin = xyzlim(2, 1) - buff; 
    ymax = xyzlim(2, 2) + buff;
    zmin = xyzlim(3, 1) - buff; 
    zmax = xyzlim(3, 2) + buff;
    
    % Extract image data
    tubi.getCurrentData()
    IV = tubi.currentData.IV ;
    
    % Extract mesh data
    flipy = tubi.flipy ;
    if use_sphere_mesh
        meshFileBase = tubi.fullFileBase.mesh ;
        meshfn = sprintf( meshFileBase, tp );
        mesh = read_ply_mod( meshfn );
    else
        mesh = tubi.getCurrentSPCutMeshSmRSC;
        
        % Switch mesh face normals to point outward for visualization purposes
        mesh.f = mesh.f(:, [2 1 3]);
        % mesh.vn = -mesh.vn;
    end
    
    % If we smooth before pushing along the normal
    if smoothing_lambda > 0
        disp('smoothing mesh via laplacian filter')
        try
            disp('attemping smoothing using gptoolbox, if configured...')
            mesh.v = laplacian_smooth(...
                mesh.v, mesh.f, 'cotan', [], smoothing_lambda, 'implicit') ;
        catch
            disp('Could not smooth mesh because no gptoolbox')
        end
    end
    
    % Make sure vertex normals are normalized
    mesh.vn = mesh.vn ./ sqrt( sum( mesh.vn.^2, 2 ) );
    % Normally evolve vertices
    mesh.v = mesh.v + normal_shift .* mesh.vn;
    % Re-orient faces
    if reorient_faces && use_sphere_mesh
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
    
    fig = figure('Visible', 'on',  'units', 'centimeters', ...
        'position', [0,0,xwidth,ywidth]) ;
    % set(gcf, 'Visible', 'Off')
    disp(['creating texture patch ' num2str(tp, '%06d')])
    
    % Allow for axis swapping
    disp(['texture axis order: [', ...
        num2str(texture_axis_order(1)), ...
        ' ', num2str(texture_axis_order(2)), ...
        ' ', num2str(texture_axis_order(3)) ']'])
    
    % Shift the texture vertices by a custom amount without affecting
    % the mesh embedding
    if any(abs(texture_shift) > 0)
        if numel(texture_shift) == numel(mesh.v(:, 1)) || numel(texture_shift(:)) == 1
            TV = mesh.v(:, texture_axis_order) + texture_shift .* mesh.vn(:, texture_axis_order) ;
        elseif numel(texture_shift) == numel(mesh.v)
            TV = mesh.v(:, texture_axis_order) + texture_shift(:, texture_axis_order) ;
        else
            error('options.texture_shift should be scalar, #vertices x 1, or #vertices x dimension')
        end
    else
        TV = mesh.v(:, texture_axis_order)  ;
    end
    
    
    % Allow for overall flip
    % --> apply rotation and translation and dilation BEFORE flipping
    VV = mesh.v ;
    if use_sphere_mesh
        if isfield(Options, 'Rotation')
            disp('rotating...')
            VV = (Options.Rotation * VV')' ;
            Options = rmfield(Options, 'Rotation') ;
        else
            disp('WARNING: no rotation supplied, using APDV frame')
            tubi.getRotTrans()
            VV = (tubi.APDV.rot * VV')' ;
        end
        if isfield(Options, 'Translation')
            disp('translating...')
            VV = VV + Options.Translation ;
            Options = rmfield(Options, 'Translation') ;
        else
            disp('WARNING: no translation supplied, using APDV frame')
            tubi.getRotTrans()
            VV = VV + tubi.APDV.trans ;
        end
        if isfield(Options, 'Dilation')
            disp('dilating...')
            VV = VV * Options.Dilation ;
            Options = rmfield(Options, 'Dilation') ;
        else
            disp('WARNING: no dilation supplied, using APDV frame')
            VV = VV * tubi.APDV.resolution ;
        end
    else
        if isfield(Options, 'Dilation')
            disp('dilating...')
            TV = TV / Options.Dilation ;
            Options = rmfield(Options, 'Dilation') ;
        else
            disp('WARNING: no dilation supplied, using APDV frame')
            TV = TV / tubi.APDV.resolution ;
        end
        if isfield(Options, 'Translation')
            disp('translating...')
            TV = TV - Options.Translation ;
            Options = rmfield(Options, 'Translation') ;
        else
            disp('WARNING: no translation supplied, using APDV frame')
            tubi.getRotTrans()
            TV = TV - tubi.APDV.trans ;
        end
        if isfield(Options, 'Rotation')
            disp('rotating...')
            TV = (inv(Options.Rotation) * TV')' ;
            Options = rmfield(Options, 'Rotation') ;
        else
            disp('WARNING: no rotation supplied, using APDV frame')
            tubi.getRotTrans()
            TV = (inv(tubi.APDV.rot) * TV')' ;
        end
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
    
    % Use for testing lighting
    % trisurf(triangulation(mesh.f, VV), ...
    %     'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
    %     'SpecularStrength', Options.SpecularStrength, ...
    %     'DiffuseStrength', Options.DiffuseStrength, ...
    %     'AmbientStrength', Options.AmbientStrength, ...
    %     'SpecularExponent', Options.SpecularExponent );
    
    
    if makeColorbar
        cb = colorbar() ;
    end
    
    % format the figure
    disp('formatting figure...')
    axis equal
    grid off
    xlim([xmin, xmax])
    ylim([ymin, ymax])
    zlim([zmin, zmax])
    
    
    % 2023-09-05 edits here:
    % colormap bone
    colormap gray_r
    
    %
    titlestr = ['$t = $ ' num2str(tp*timeinterval-t0) ' ' timeunits] ;
    if blackFigure
        title(titlestr, 'Interpreter', 'Latex', 'Color', 'white')
    else
        title(titlestr, 'Interpreter', 'Latex', 'Color', 'k')
    end
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
    if blackFigure
        set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w', 'zcol', 'w')
        set(gcf, 'InvertHardCopy', 'off');
        set(gcf, 'Color', 'k')
        set(gcf, 'color', 'k')
    else
        set(gcf, 'Color', 'w')
    end
    
    % Check that mesh is oriented correctly
    % trisurf(triangulation(mesh.f, VV), 'edgecolor', 'none')
    % xlabel('x'); ylabel('y'); zlabel('z')
    % axis equal
    
    get(gcf, 'position')
    set(gcf, 'position', [0, 0, xwidth, ywidth])
    
    if strcmpi(viewType, 'dorsal')
        view(0, 90);
    elseif strcmpi(viewType, 'ventral')
        view(0, 270);
    elseif strcmpi(viewType, 'lateral1')
        view(0, 0);
    elseif strcmpi(viewType, 'lateral2')
        view(0, 180);
    elseif strcmpi(viewType, 'perspective')
        view(-20, 20);
    else
        error('Invalid view orientation');
    end
    
    
    % Optional: camlight
    % cl = camlight;
    % cl.Position = [134.887298762337,22.8013031584538,1562.27254887447];
    
    % export_fig(sprintf('TubULAR_Paper_Figures/Textured_Images/Textured_Surface_Dorsal_T%06d.pdf', timePoints(tp)), '-eps', '-r200', '-depsc')
    if blackFigure
        export_fig(sprintf('TubULAR_Paper_Figures/Textured_Images_black/Textured_Surface_Dorsal_T%06d_r600.png', timePoints(tp)), '-png', '-r600')
    else
        export_fig(sprintf('TubULAR_Paper_Figures/Textured_Images/Textured_Surface_Dorsal_T%06d_r600.png', timePoints(tp)), '-png', '-r600')
    end
    close all
end