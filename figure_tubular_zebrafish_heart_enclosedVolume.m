%% Start by running example_zebrafish_heart.m until tubi definition.

cd /mnt/data/tubular_test/zebrafish_heart/analysis

%% ************************************************************************
% *************************************************************************
%               GENERATE ENCLOSED VOLUME PLOT
% *************************************************************************
% *************************************************************************

%% Generate Enclosed Volume Plot ==========================================
close all; clc;

% Set options
smoothing_lambda = 0.0;
normal_shift = 0.0;
flipy = false;
reorient_faces = false;
Options = struct();

% Time points for which to calculate the enclosed volume
timePoints = tubi.xp.fileMeta.timePoints ;

% Unpack metadat and save
timeinterval = tubi.timeInterval ;
timeunits = tubi.timeUnits ;

% Unpack xyzlim
[~, ~, ~, xyzlim] = tubi.getXYZLims() ;
metadat.xyzlim = xyzlim ;
xmin = xyzlim(1, 1); xmax = xyzlim(1, 2) ;
ymin = xyzlim(2, 1); ymax = xyzlim(2, 2) ;
zmin = xyzlim(3, 1); zmax = xyzlim(3, 2) ;

% Extract mesh volumes and areas
meshVolumes = zeros(size(timePoints));
meshAreas = zeros(size(timePoints));
for tp = timePoints
    
    tidx = tubi.xp.tIdx(tp);
    
    % Extract (spherical) mesh data
    meshFileBase = tubi.fullFileBase.mesh ;
    meshfn = sprintf( meshFileBase, tp );
    mesh = read_ply_mod( meshfn );
    
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
    if reorient_faces
        disp('reorienting faces...')
        mesh.f = reorient_facets( mesh.v, mesh.f );
    end
    
    VV = mesh.v ;
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
    
    if flipy
        VV(:, 2) = -VV(:, 2) ;
    end
    
    [meshVolumes(tidx), meshAreas(tidx)] = ...
        meshVolumeArea(mesh.v, mesh.f);
    
end

meshVolumes = meshVolumes ./ meshVolumes(1);
meshAreas = meshAreas ./ meshAreas(1);

%%
close all; clc;
T = 11; % Period size
timeSteps = tubi.xp.fileMeta.timePoints(1:end) - tubi.t0;
timeSteps = timeSteps / T;

xLim = [timeSteps(1) timeSteps(end)];
% yLim = [-5, 5];

fig = figure('Visible', 'on',  'units', 'centimeters') ;

% Assumes A and B are the same size
plot_lw = 1;
plot(timeSteps, meshVolumes, ...
    'LineWidth', plot_lw, 'Color', tubi.plotting.colors(1,:));
hold on
plot(timeSteps, meshAreas, ...
    'LineWidth', plot_lw, 'Color', tubi.plotting.colors(2,:));
hold off

xlim(xLim);
% ylim(yLim);

% xlabel(['time [' tubi.timeUnits ']'], 'Interpreter', 'Latex');
xlabel('time [T]');

% ylabel( ...
%     {['\color[rgb]{' num2str(tubi.plotting.colors(1,:)) '} Volume $$V/V_0$$'], ...
%     ['\color[rgb]{' num2str(tubi.plotting.colors(2,:)) '} Surface area $$A/A_0$$']}, ...
%     'Interpreter', 'Latex');
ylabel( ...
    {['\color[rgb]{0, 0.4470, 0.7410} volume_{}V/V_0'], ...
    ['\color[rgb]{0.85, 0.3250, 0.0980} surface_{}area_{}A/A_0']}, ...
    'Interpreter', 'tex');

set(gcf, 'Color', [1 1 1]);

% Resize Figure for Paper -------------------------------------------------
set(fig, 'Units', 'centimeters');
set(fig, 'position', [0,0,4,4]) ;
%ratio = fig.Position(4) ./ fig.Position(3);
%fig.Position(3) = 3.5;
%fig.Position(4) = ratio * fig.Position(3);

set(gca, 'FontSize', 10);
set(gca, 'FontWeight', 'normal');

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');

% save data
save('./TubULAR_Paper_Figures/Enclosed_Volume.mat',...
    'timeSteps', 'meshAreas', 'meshVolumes')

% Save figure
saveas(fig, './TubULAR_Paper_Figures/Enclosed_Volume_01.pdf')

