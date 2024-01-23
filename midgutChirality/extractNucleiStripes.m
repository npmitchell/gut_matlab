% Draw curves of muscle along AP axis on pullback images
clc; close all; 
clearvars -except tubi

addpath('/mnt/data/code/gut_matlab/travelling_salesman/tsp_ga')

outdir = fullfile(tubi.dir.tracking, 'muscleStripes');
trackfn = fullfile(tubi.dir.tracking, 'muscle_tracks.mat') ;
tmp = load(trackfn) ;
tracks = tmp.tracks ;
imFileBase = tubi.fullFileBase.im_sp_sm ;
preview = false ;

ncurves = 4 ;
timePoints = tubi.xp.fileMeta.timePoints ;

for curv = 0:ncurves-1 
        
    % First grab all tracks available for this curv by looking for any
    % point in a chosen polygon enclosing the curve at t0
    tidx0 = tubi.xp.tIdx(tubi.t0) ;
    % Gather all tracks at tidx0
    ids = [] ;
    for track = 1:length(tracks)
        ids = cat(1, ids, tracks{track}(tidx0, :)) ;
    end

    polyfn = fullfile(outdir, sprintf('polygon%d_t0.mat', curv)) ;
    if ~exist(polyfn, 'file')
        % Gather all tracks for this stripe
        disp('Gathering all tracks for this stripe at t=t0')
    
        % Load image
        close all
        im = imread(sprintf(imFileBase, tubi.t0)) ;
        imshow(im) ;
        hold on;
        plot(ids(:, 1), ids(:, 2), 'o')
    
        % Select tracks that are part of this curve
        h = drawpolygon;
        
        % Wait for the user to finish drawing the polygon
        wait(h);
        
        % Get the polygon's vertices
        polygonVertices = h.Position;
    
        % Save the polygon
        save(polyfn, 'polygonVertices') 
    else
        load(polyfn, 'polygonVertices')
    end
    
    % Check which points are inside the drawn polygon
    inPolygon = inpolygon(ids(:, 1), ids(:, 2), ...
        polygonVertices(:, 1), polygonVertices(:, 2));
    
    % Select the xy points that are inside the polygon
    sids = ids(inPolygon, :);
    
    % Display the selected points on the image
    if preview
        clf
        im = imread(sprintf(imFileBase, tubi.t0)) ;
        imshow(im) ;
        hold on;
        plot(ids(:, 1), ids(:, 2), 'o')
        hold on;
        plot(sids(:, 1), sids(:, 2), 'ro', 'MarkerSize', 10);
        title('Tracks that are part of this stripe')
        waitfor(gcf)
    end

    % Load the curve positions that are already defined on disk
    curvfn = fullfile(outdir, sprintf('curvePositions_stripe%d.mat', curv)) ;
    if exist(curvfn, 'file')
        load(curvfn, 'curvePositions')
    else
        curvePositions = cell(length(timePoints), 1) ;
    end

    % Now add to these tracks to fill in the stripe
    allPts = {} ;
    for tidx = 1:length(timePoints)
        tp = timePoints(tidx) ;

        % Get known tracks for this timepoint
        disp(['Gathering tracks for this stripe at t=' num2str(tp)])
        ids = [] ;
        for track = 1:length(tracks)
            ids = cat(1, ids, tracks{track}(tidx, :)) ;
        end
        % Keep only selected tracks for this curve
        ids = ids(inPolygon, 1:2) ; 


        % Add to the current ids
        if isempty(curvePositions{tidx})
    
            % Add more points as needed to draw curves along rostral/caudal axis
            clf
            im = imread(sprintf(imFileBase, tp)) ;
            imshow(im) ;
            hold on;
            plot(ids(:, 1), ids(:, 2), 'o')

            h = drawfreehand;
            % Wait for the user to modify the curve
            wait(h);
            % Store the modified coordinates for the current image
            curvePositions{tidx} = h.Position;
            % Save updated curvePositions
            disp('saving updates...')
            save(curvfn, 'curvePositions')
        else
            % h = drawpolygon(gca, 'Position', curvePositions{tidx});
            % title(['t=' num2str(tp)])
            % pause(0.1)
        end
        
        % Adding to all points
        % ids = [ids; curvePositions{tidx}] ;
        % allPts{tidx} = curvePositions{tidx} ;

        % Connect these points in order using traveling salesman
        % userConfig = struct('xy', ids);
        % res = tsp_ga(userConfig);

        % error('here')
    end
end

%% Process the curves into regularly-sampled 3D curves
close all
clc

% Load single image to get dimensions
im = imread(sprintf(imFileBase, tubi.t0)) ;
[~,~,~,xyzlims] = tubi.getXYZLims() ;
Xmax = size(im, 2) ;
colors = define_colors ;


% Push stripes to 3d
dat = cell(length(timePoints), 1) ;
for tidx = 1:length(timePoints)
    tp = timePoints(tidx) ;
    tubi.setTime(tp) ;

    % Obtain mesh and centerline
    mesh = tubi.getCurrentSPCutMeshSmRS() ;
    meshc = tubi.getCurrentSPCutMeshSmRSC() ;
        
    % preallocate
    xyz = cell(4, 1) ;
    uvs = cell(4,1) ;
    
    % Push each curve to 3D
    for curv = 0:ncurves-1
    
        % Load the curve positions that are already defined on disk
        curvfn = fullfile(outdir, sprintf('curvePositions_stripe%d.mat', curv)) ;
        load(curvfn, 'curvePositions')
    
        XY = curvePositions{tidx} ;
        [XX, uids] = unique(XY(:, 1)) ; 
        YY = XY(uids, 2) ;

        % interpolate and resample
        Xi = linspace(1, Xmax, nUi) ;
        Yi = interp1(XX, YY, Xi, 'linear', 'extrap') ;

        % check it
        % clf
        % plot(XY(:, 1), XY(:, 2), '.-') ; hold on;
        % plot(Xi, Yi, '.-')

        % Push to 3D
        uv = tubi.XY2uv(im, [Xi;Yi]', false, 1, 1) ;

        % Change cover if we roll into negative values
        push = uv(:, 2) < 0 ;
        uv(push, 2) = uv(push, 2) + 1 ;

        % Push uv away from boundaries a little
        push = uv(:, 1) == 0 ;
        uv(push, 1) = eps ;
        uvs{curv+1} = uv ;
        xyz{curv+1} = tubi.uv2APDV(uv, 'spsmrs', 1, 1) ;
    
        assert(~any(any(isnan(xyz{curv+1}))))
    end

    % Plot the result
    if preview 
        close all;
        figure('position', [0 0 875 656])
        h = trisurf(triangulation(meshc.f, meshc.v), 'edgecolor', 'none', ...
             'facecolor', 0.8 * [1,1,1], 'facealpha', 0.7) ;
            % 'faceVertexCData', meshc.u(:, 1), 'facealpha', 0.5) ;
        lighting gouraud
        camlight
        hold on;
        plot3(xyz{1}(:, 1), xyz{1}(:, 2), xyz{1}(:, 3), '-', 'linewidth', 2) ;
        plot3(xyz{2}(:, 1), xyz{2}(:, 2), xyz{2}(:, 3), '-', 'linewidth', 2) ;
        plot3(xyz{3}(:, 1), xyz{3}(:, 2), xyz{3}(:, 3), '-', 'linewidth', 2) ;
        plot3(xyz{4}(:, 1), xyz{4}(:, 2), xyz{4}(:, 3), '-', 'linewidth', 2) ;
        axis equal
        axis off
        xlim(xyzlims(1,:))
        ylim(xyzlims(2,:))
        zlim(xyzlims(3,:))
        set(gcf, 'color', 'w')
        drawnow
        title(['t=' num2str(tp) ' ' tubi.timeUnits])
        saveas(gcf, fullfile(outdir, 'images', 'perspective',...
            sprintf('muscleStripes_%06d.png', tp)))
        view(2)
        saveas(gcf, fullfile(outdir, 'images', 'dorsal',...
            sprintf('muscleStripes_dorsal_%06d.png', tp)))
    end

    % Save the result
    dat{tidx} = struct('mesh', mesh, 'meshc', meshc, ...
        'curv1', xyz{1}, 'curv2', xyz{2}, ...
        'curv3', xyz{3}, 'curv4', xyz{4}, ...
        'uv1', uvs{1}, 'uv2', uvs{2}, 'uv3', uvs{3}, 'uv4', uvs{4}) ;

    disp(['done with tp=' num2str(tp)])

end

% Save dat
datfn = fullfile(outdir, 'muscleNucleiStripes.mat') ;
save(datfn, 'dat')

%% Measure twist from this coarse grained curves
clc

nUi = 150 ;
% Load dat
datfn = fullfile(outdir, 'muscleNucleiStripes.mat') ;
load(datfn, 'dat')

omegas = zeros(length(timePoints)-1, 4, nUi) ;
twists = omegas(:, :, 1:end-1) ;
for tidx = 1:length(timePoints)-1
    tp = timePoints(tidx) ;
    disp(['tp= ' num2str(tp)])
    tubi.setTime(tp) ;
    d0 = dat{tidx} ;
    d1 = dat{tidx+1} ;
    for curv = 1:4 
        if curv == 1
            v0 = d1.curv1 - d0.curv1;
            curv0 = d0.curv1 ;
            uv = d0.uv1 ;
        elseif curv == 2
            v0 = d1.curv2 - d0.curv2;
            curv0 = d0.curv2 ;
            uv = d0.uv2 ;
        elseif curv == 3
            v0 = d1.curv3 - d0.curv3;
            curv0 = d0.curv3 ;
            uv = d0.uv3 ;
        elseif curv == 4
            v0 = d1.curv4 - d0.curv4;
            curv0 = d0.curv4 ;
            uv = d0.uv4 ;
        end
        % Project onto the mesh
        nearestVtx = pointMatch(curv0, d0.meshc.v) ;
        % % Subtract the normal component
        % vt = v0 - dot(v0, d0.meshc.vn(nearestVtx)) ;
        % To isolate the phihat direction, use 2d representation
        % First get fieldfaces in 2D image space
        uv(uv(:, 1) < 1e-6, 1) = uv(uv(:, 1) < 1e-6, 1) + 1e-6 ;
        fieldfaces = pointLocation( triangulation(d0.mesh.f, d0.mesh.u), uv) ;
        [v0n, v0t, v0t2d] = ...
            resolveTangentNormalVelocities(d0.mesh.f, d0.mesh.v, v0, fieldfaces, ...
            d0.mesh.u) ;
        v0t2d(:, 1) = 0 ;
        % push back into 3d
        jac2d_to_3d = jacobian2Dto3DMesh(d0.mesh.u, d0.mesh.v, d0.mesh.f) ;
        vphi3d = zeros(size(v0t2d, 1), 3) ;
        for qq = 1:length(v0t2d)
            vphi3d(qq, :) = jac2d_to_3d{qq} * v0t2d(qq, :)' ;
        end
        
        if preview % || mod(tidx, 10) == 0 
            if curv == 1
                close all;
                h = trisurf(triangulation(d0.meshc.f, d0.meshc.v), 'edgecolor', 'none', ...
                     'facecolor', 0.8 * [1,1,1], 'facealpha', 0.7) ;
                    % 'faceVertexCData', meshc.u(:, 1), 'facealpha', 0.5) ;
                lighting gouraud
                camlight
                hold on;
            end
            quiver3(curv0(:, 1), curv0(:, 2), curv0(:, 3), ...
                vphi3d(:, 1), vphi3d(:, 2), vphi3d(:, 3), 0)
            if curv == 4
                axis equal
                axis off
                xlim(xyzlims(1,:))
                ylim(xyzlims(2,:))
                zlim(xyzlims(3,:))
                set(gcf, 'color', 'w')
                drawnow
                title(['t=' num2str(tp) ' ' tubi.timeUnits])
                saveas(gcf, fullfile(outdir, 'images', 'vphi_perspective',...
                    sprintf('muscleStripes_%06d.png', tp)))
                view(2)
                saveas(gcf, fullfile(outdir, 'images', 'vphi_dorsal',...
                    sprintf('muscleStripes_dorsal_%06d.png', tp)))
            end
        end

        % Convert to angular velocity
        omega = sign(v0t2d(:, 2)) .* ...
            vecnorm(vphi3d, 2, 2) ./ d0.meshc.radius_um(nearestVtx) ;

        % store angular velocity
        omegas(tidx, curv, :) = omega ;
        % Twist is gradient of angular velocity wrt fractional length
        ds = 1/nUi ;
        twists(tidx, curv, :) = diff(omega) / ds ;
    end
end

%% Plot omegas and twists
close all
tps = timePoints * tubi.timeInterval ; % timepoints in proper units
labels = {'right ventral', 'right dorsal', 'left dorsal', 'left ventral'} ;
endCut = 10 ;
endFrac = endCut / nUi ;
soL = linspace(endFrac, 1-endFrac, nUi-2*endCut+1) ; % s/L
clims = [0.5, 5 ] ;
Figlabels = {'omega', 'twist'} ;

for ot = 1:2
    for qq = 1:4
        subplot(2, 2, qq)
        if ot == 1
            kymo = squeeze(omegas(:, qq, :)) ;
            kymo = imgaussfilt(kymo, 1) ;
        else
            kymo = squeeze(twists(:, qq, :)) ;
            kymo = imgaussfilt(kymo, 2) ;
        end
        kymo = kymo(:, endCut:end-endCut) ;
        imagesc(soL, tps, kymo)
        hold on;
        caxis(clims(ot) * [-1,1])
        if ismember(qq, [3,4])
            xlabel('RC position, s/L'); 
        end
        if ismember(qq, [1,3])
            ylabel(['time [' tubi.timeUnits ']']); 
        end
        title(labels{qq}, 'FontWeight', 'normal')

        % Show fold locations
        tubi.getFeatures() ;
        folds = tubi.features.folds ;
        for foldID = 1:3
            fx = double(folds(tubi.features.fold_onset(foldID):end, foldID)) / tubi.nU ;
            ftps = (tubi.features.fold_onset(foldID):length(timePoints)) * tubi.timeInterval ;
            plot(fx, ftps, '-')
        end
    end
    colormap bwr
    colorbar(gca, 'Position', [0.93, 0.2, 0.01, 0.6])
    set(gcf, 'color', 'w')
        
    if ot == 1
        sgtitle(['\omega [rad/' tubi.timeUnits ']'])
    else
        sgtitle(['Tw [rad/' tubi.timeUnits ' L]'])
    end
    saveas(gcf, fullfile(outdir, 'images', ['kymo_stripes_' Figlabels{ot} '.pdf']))
end

%% average the omegas for all bands
clf
omegaAvg = imgaussfilt(squeeze(mean(omegas, 2)), 1) ;
twistAvg = imgaussfilt(squeeze(mean(twists, 2)), 2) ;
omegaAvg = omegaAvg(:, endCut:end-endCut) ;
twistAvg = twistAvg(:, endCut:end-endCut) ;
imagesc(soL, tps, omegaAvg)
caxis(clims(1) * [-1,1])
xlabel('RC position, s/L'); 
ylabel(['time [' tubi.timeUnits ']']); 
colormap bwr
colorbar(gca, 'Location', 'eastOutside') 
set(gcf, 'color', 'w')
sgtitle(['\langle\omega\rangle_{\phi} [rad/' tubi.timeUnits ']'])
% Show fold locations
tubi.getFeatures() ;
folds = tubi.features.folds ;
for foldID = 1:3
    hold on;
    fx = double(folds(tubi.features.fold_onset(foldID):end, foldID)) / tubi.nU ;
    ftps = (tubi.features.fold_onset(foldID):length(timePoints)) * tubi.timeInterval ;
    plot(fx, ftps, '-')
end
saveas(gcf, fullfile(outdir, 'images', ['kymo_mean_' Figlabels{1} '.pdf']))

imagesc(soL, tps, twistAvg)
caxis(clims(2) * [-1,1])
xlabel('RC position, s/L'); 
ylabel(['time [' tubi.timeUnits ']']); 
colormap bwr
colorbar(gca, 'Location', 'eastOutside') 
set(gcf, 'color', 'w')
sgtitle(['\langleTw\rangle_{\phi} [rad/' tubi.timeUnits ' L]'])
tubi.getFeatures() ;
folds = tubi.features.folds ;
for foldID = 1:3
    hold on;
    fx = double(folds(tubi.features.fold_onset(foldID):end, foldID)) / tubi.nU ;
    ftps = (tubi.features.fold_onset(foldID):length(timePoints)) * tubi.timeInterval ;
    plot(fx, ftps, '-')
end
saveas(gcf, fullfile(outdir, 'images', ['kymo_mean_' Figlabels{2} '.pdf']))


%% Geodesic distance between DR-VR, VR-VL, VL-DL, DL-DR

% Here we first find the distance along isolines of s rather than true
% geodesics.
for tidx = 1:length(timePoints)-1
    tp = timePoints(tidx) ;
    disp(['tp= ' num2str(tp)])
    tubi.setTime(tp) ;
    d0 = dat{tidx} ;
    d1 = dat{tidx+1} ;
    for curv = 1:4 
        % [D,u,X,div_X,phi,pre,B,t] = heat_geodesic(V) ;
        % [path,vlist,plist] = compute_geodesic_mesh(D, vertex, face, x, options) 
        
        curv0 = curv ;
        curv1 = mod(curv + 1, 4);
        if curv1 == 0 
            curv1 = 4 ;
        end


        % sample the curve
        % Load the curve positions that are already defined on disk
        curvfn = fullfile(outdir, sprintf('curvePositions_stripe%d.mat', curv0)) ;
        load(curvfn, 'curvePositions')
    
        XY0 = curvePositions{tidx} ;
        [XX0, uids] = unique(XY0(:, 1)) ; 
        YY0 = XY0(uids, 2) ;

        % Load the adjacent curve
        curvfn = fullfile(outdir, sprintf('curvePositions_stripe%d.mat', curv1)) ;
        load(curvfn, 'curvePositions')
        XY1 = curvePositions{tidx} ;
        [XX1, uids] = unique(XY0(:, 1)) ; 
        YY1 = XY1(uids, 2) ;

        % interpolate and resample
        Xi = linspace(1, Xmax, nUi) ;
        Yi = interp1(XX, YY, Xi, 'linear', 'extrap') ;

        % Push to 3D
        uv = tubi.XY2uv(im, [Xi;Yi]', false, 1, 1) ;

        % Change cover if we roll into negative values
        push = uv(:, 2) < 0 ;
        uv(push, 2) = uv(push, 2) + 1 ;

        % Push uv away from boundaries a little
        push = uv(:, 1) == 0 ;
        uv(push, 1) = eps ;
        uvs{curv+1} = uv ;
        xyz{curv+1} = tubi.uv2APDV(uv, 'spsmrs', 1, 1) ;

    end
end


%% Measure tracked nuclei dynamics 