function getSkelCurvature(SE, overwriteData, overwriteImages)
% 
if nargin < 2
    overwriteData = false ;
    overwriteImages = false ;
elseif nargin < 3
    overwriteImages = false ;
end
if ~exist(SE.zDir.curvature, 'dir')
    mkdir(SE.zDir.curvature)
end
curvatureFn = fullfile(SE.zDir.curvature, 'curvature.mat') ;
if ~exist(curvatureFn, 'file') || overwriteData
    kymo = SE.getKymographs() ;
    curves = kymo.resampled ;
    landmarks = kymo.landmarks ;
    xyv = cat(3, curves.x_ventral, curves.y_ventral) ;
    xyd = cat(3, curves.x_dorsal, curves.y_dorsal) ;
    xyv = xyv * SE.resolution ;
    xyd = xyd * SE.resolution ;
    kappa_ventral = zeros(size(curves.x_ventral)) ;
    kappa_dorsal = zeros(size(curves.x_dorsal)) ;
    for tidx = 1:length(SE.timePoints)
        % Smooth these curves
        degree = 2 ;
        vx = smooth(squeeze(xyv(tidx, :, 1)), 50, 'sgolay',degree) ;
        vy = smooth(squeeze(xyv(tidx, :, 2)), 50, 'sgolay',degree) ;
        dx = smooth(squeeze(xyd(tidx, :, 1)), 50, 'sgolay',degree) ;
        dy = smooth(squeeze(xyd(tidx, :, 2)), 50, 'sgolay',degree) ;
        vv = [vx, vy] ;
        dd = [dx, dy] ;

        kappa_ventral(tidx, :) = -LineCurvature2D(vv) ;
        kappa_dorsal(tidx, :) = -LineCurvature2D(dd) ;
    end
    save(curvatureFn, 'kappa_ventral', 'kappa_dorsal', ...
        'curves', 'landmarks')
else
    load(curvatureFn, 'kappa_ventral', 'kappa_dorsal', ...
        'curves', 'landmarks')
end            

% Plot it
VkappaFigFn = fullfile(SE.zDir.curvature, ...
    'curvature_ventral.png') ;
VkappaFigFn_zoom = fullfile(SE.zDir.curvature, ...
    'curvature_ventral_zoom.png') ;
DkappaFigFn = fullfile(SE.zDir.curvature, ...
    'curvature_dorsal.png') ;
DkappaFigFn_zoom = fullfile(SE.zDir.curvature, ...
    'curvature_dorsal_zoom.png') ;

if ~exist(VkappaFigFn, 'file') || overwriteImages
    v1 = landmarks.v1 ;
    v2 = landmarks.v2 ;
    v3 = landmarks.v3 ;
    d1 = landmarks.d1 ;
    d2 = landmarks.d2 ;
    d3 = landmarks.d3 ;
    fns = {VkappaFigFn, DkappaFigFn} ;
    fns_zoom = {VkappaFigFn_zoom, DkappaFigFn_zoom} ;
    kymos = {kappa_ventral, kappa_dorsal} ;
    ps = max(curves.sdorsal, [], 2) ;
    xx = {curves.sventral, ps-curves.sdorsal} ;
    yy = {SE.dt * (SE.timePoints .* ones(1000, length(SE.timePoints)))', ...
        SE.dt * (SE.timePoints .* ones(1000, length(SE.timePoints)))'} ;
    lms = {{v1, v2, v3}, {ps - d1, ps - d2, ps - d3}} ;
    ktitleAdd = {'ventral side', 'dorsal side'} ;
    for kId = 1:2
        close all
        p1 = surf(xx{kId}, yy{kId}, kymos{kId},... 
            'FaceColor','interp',...
            'edgecolor','none');
        cb = colorbar() ;
        colormap(parula)
        ylabel(cb, ['curvature [' SE.spaceUnits '$^{-1}$]'], 'Interpreter', 'Latex') ;
        hold on ;
        view(2)
        grid off
        set(gca,'Ydir','reverse')
        ylabel(['time, [' SE.timeUnits ']'], 'interpreter', 'latex')
        xlabel(['pathlength, $s$ [' SE.spaceUnits ']'], 'interpreter', 'latex')
        title(['curvature: ' ktitleAdd{kId}])
        % Plot landmarks
        lm = lms{kId} ;
        for lid = 1:length(lm) 
            % Plot in 3d to be above the surface data
            plot3(lm{lid}, SE.timePoints * SE.dt, ...
                (max(kymos{kId}(:))+1) * ones(size(SE.timePoints)), 'o', ...
                'color', 'k'); 
            hold on;
        end
        ylim([min(SE.timePoints-0.5)*SE.dt, max(SE.timePoints+0.5)*SE.dt])
        caxis([-2, 2])
        colormap(bluewhitered(256)) ;
        saveas(gcf, fns{kId})  
        caxis([-0.5, 0.5])
        colormap(bluewhitered(256)) ;
        saveas(gcf, fns_zoom{kId}) 
        close all 
    end
end