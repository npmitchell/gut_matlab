%% Color lagrangian surface
datdir = '/mnt/data/analysis/tubular/lagrangianDemo/'
mkdir(datdir)
trueLagrangian = true ;
useAngles = false ;

[~,~,~, xyzlims] = QS.getXYZLims() ;
% get face connectivity
QS.setTime(QS.t0set()) ;
refMesh = QS.getCurrentRicciMesh() ;
mesh0 = QS.getCurrentSPCutMeshSmRS();

rr = mesh0.u(:, 1) ;
th = atan2(refMesh.annulus.u(:, 2), refMesh.annulus.u(:, 1)) ;
th(QS.nU*(QS.nV-1)+1:QS.nU*QS.nV) = th(1:QS.nU) ;

viewAngles = [-0.75, -1.0, 0.7] ;

% get pathlines
plines = QS.getPullbackPathlines() ;

tps = [123:10:243] ;
for tp = tps

    tidx = QS.xp.tIdx(tp) ;
    xx = squeeze( plines.vertices3d.vXrs(tidx, :,:) ) ;
    yy = squeeze( plines.vertices3d.vYrs(tidx, :,:) ) ;
    zz = squeeze( plines.vertices3d.vZrs(tidx, :,:) ) ;
    

    t0Pathlines = QS.t0set() ;
    QS.setTime(QS.xp.fileMeta.timePoints(tidx)) ;
    
    clf
    set(gcf, 'color', 'w')
    
    mesh = QS.getCurrentSPCutMeshSmRS();
    
    if trueLagrangian
        mesh.v = [xx(:), yy(:), zz(:)] ;
    end    
    
    %% PLOT BASED ON MU
    % strain = QS.getCurrentPathlineStrain( t0Pathlines) ;
    % mu = strain.beltrami.mu_material_filtered ;
    % opts = struct('mesh', mesh, 'clim_mag', 0.8) ;
    % mag = abs(mu) ;
    % theta = atan2(imag(mu), real(mu));
    % plotNematicField(mag, theta, opts)
    
    
    
    %% PLOT BASED ON POSITION
    % opts = struct('mesh', mesh, 'clim_mag', 0.8) ;
    % plotNematicField(rr, th, opts)
    
    
    %% PLOT BASED ON POSITION -- QUALITATIVE
    
    load('romaO.mat')
    load('vikO.mat')
    load('bamO.mat')
    bamO = circshift(bamO,115) ;
    mesh.u(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) ;
    pulse1 = 1-2*abs((mesh.u(:, 1)-0.5))  ;
    pulse2 = max(0,1-2*mesh.u(:, 1)) ;
    pulse3 = max(0,2*mesh.u(:, 1)-1) ;
    colors = pulse1 .* bamO(round(mesh.u(:, 2) * 255+1), :) + ... % middle
        pulse2.* vikO(256-round(mesh.u(:, 2) * 254), :) + ... % start
        pulse3 .* romaO(256-round(mesh.u(:, 2) * 254), :) ;    % end
        
    h = trisurf(triangulation(mesh.f, mesh.v), ...
        'FaceVertexCData', colors, 'edgecolor', 'none') ;
    axis equal
    
    %% handling
    axis equal
    axis tight
    grid off
    if useAngles
        view(viewAngles)
    end
    
    % lgt = camlight('headlight') ;
    lighting gouraud

    shading interp
    lightangle(-45,30)
    h.FaceLighting = 'gouraud';
    h.AmbientStrength = 0.9;
    h.DiffuseStrength = 0.8;
    h.SpecularStrength = 0.1;
    h.SpecularExponent = 25;
    % h.BackFaceLighting = 'unlit';
    xlim(xyzlims(1, :))
    ylim(xyzlims(2, :))
    zlim(xyzlims(3, :))
    
    
    export_fig(fullfile(datdir, sprintf('snap_%06d.png', tp)), '-nocrop', '-r600')
end



    h = trisurf(triangulation(mesh.f, [rr(:).*cos(th(:)), rr(:).*sin(th(:)), 0*rr(:)]), ...
        'FaceVertexCData', colors, 'edgecolor', 'none') ;
    axis equal
    
    

    h = trisurf(triangulation(mesh.f, cat(2, mesh.u, mesh.u(:, 1)*0)), ... %[rr(:)/max(rr(:)), th(:)/(2*pi), 0*rr(:)]), ...
        'FaceVertexCData', colors, 'edgecolor', 'none') ;
    axis equal
    
view(2)
xlabel('x')
axis off
export_fig(fullfile(datdir, sprintf('parameterization_crop.png', tp)), '-r600')


