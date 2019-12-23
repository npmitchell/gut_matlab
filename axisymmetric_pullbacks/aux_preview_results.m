%% Show 3D Texture Patch (Optional) ---------------------------------------

% Texture patch options
Options.PSize = 5 ;
Options.EdgeColor = 'none';

% Raw stack data
IV = xp.stack.image.apply();
IV = imadjustn(IV{1});

% Make a full screen image
figure('units', 'normalized', 'outerposition', [0 0 1 1] )

% Plot texture patch cut mesh in 2D
texture_patch_3d( mesh.f, mesh.v, ...
    mesh.f, mesh.v(:, [2 1 3]), IV, Options );

axis equal

colormap bone

clear Options IV

%% Show 2D Texture Patch --------------------------------------------------
% Texture patch options
Options.PSize = 5;
Options.EdgeColor = 'none';

% Raw stack data
IV = xp.stack.image.apply();
IV = imadjustn(IV{1});

% Make a full screen image
figure('units', 'normalized', 'outerposition', [0 0 1 1])

% Plot texture patch cut mesh in 2D
texture_patch_3d( cutMesh.f, cutMesh.u, ...
    cutMesh.f, cutMesh.v(:, [2 1 3]), IV, Options );

axis equal

% Format axes
% xlim([0 1]); ylim([0 1]);
% set(gca, 'xtick', []);
% set(gca, 'ytick', []);

colormap bone

clear Options IV

%% Show Tiled 2D Texture Patch --------------------------------------------

% Texture patch options
Options.PSize = 5;
Options.EdgeColor = 'none';

% Raw stack data
IV = xp.stack.image.apply();
IV = imadjustn(IV{1});

% Make a full screen image
figure('units', 'normalized', 'outerposition', [0 0 1 1])

% Plot texture patch cut mesh in 2D
texture_patch_3d( cutMesh.f, cutMesh.u, ...
    cutMesh.f, cutMesh.v(:, [2 1 3]), IV, Options );

hold on

texture_patch_3d( cutMesh.f, [ cutMesh.u(:,1), cutMesh.u(:,2)+1 ], ...
    cutMesh.f, cutMesh.v(:, [2 1 3]), IV, Options );

texture_patch_3d( cutMesh.f, [ cutMesh.u(:,1), cutMesh.u(:,2)-1 ], ...
    cutMesh.f, cutMesh.v(:, [2 1 3]), IV, Options );

axis equal

% Format axes
xlim([0 a]); ylim([0 1]);
set(gca, 'xtick', []);
set(gca, 'ytick', []);

colormap bone
clear Options IV