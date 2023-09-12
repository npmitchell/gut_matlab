%% Create minimal data for figure 3

addpath('/mnt/data/code/gut_matlab/plotting/')

%% Clear workspace ========================================================
% We start by clearing the memory and closing all figures
clear; close all; clc;
cd /mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/
zwidth =1 ;
% cd /mnt/data/48Ygal4-UAShistRFP/201904031830_great/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2/data/deconvolved_16bit/
% zwidth = 1 ;
% cd /mnt/data/handGAL4klarHandGFPhistGFP/202105072030_1mWGFP/deconvolved_16bit/
% % cd /mnt/data/mef2GAL4klarUASCAAXmChHiFP/202003151700_1p4um_0p5ms3msexp/
% zwidth = 2 ;


dataDir = cd ;

%% ADD PATHS TO THIS ENVIRONMENT ==========================================
origpath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(origpath))
addpath(fileparts(origpath))
addpath(genpath('../'))
addpath(genpath('../utility'))
addpath(genpath('/mnt/data/code/gptoolbox'))
addpath('../TexturePatch')
addpath('../DECLab')
addpath(fullfile('../utility','plotting'))
addpath(fullfile('../utility','plotting'))
% go back to the data
cd(dataDir)

%% load xp
disp('loading xp struct from disk')
load(fullfile(dataDir, 'xp.mat'), 'xp', 'opts')

%% TubULAR class instance
disp('defining TubULAR class instance (tubi= tubular instance)')
tubi = TubULAR(xp, opts) ;
disp('done defining TubULAR instance')

%% Metric Kinematics Kymographs & Correlations -- Bandwidth Filtered
% options = struct() ;
% tubi.plotMetricKinematics(options)

top = [0,90] ;
cmap1 = brewermap(256, '*RdBu') ;
cmap2 = bwr ;
cmaps = {cmap1, cmap2} ;
cmapname = {'brewermapRdBu', 'bwr'} ;

for cmapID = 1:2
    cmap = cmaps{cmapID} ;
    
    lambda = tubi.smoothing.lambda ;
    lambda_mesh = tubi.smoothing.lambda_mesh ;
    lambda_err = tubi.smoothing.lambda_err ;
    nmodes = tubi.smoothing.nmodes ;
    zwidth = tubi.smoothing.zwidth ;
    nU = tubi.nU ;
    nV = tubi.nV ;
    climit = 0.2 ;
    climit_veln = climit * 10 ;
    
    mKDir = fullfile(tubi.dir.metricKinematics.root, ...
        strrep(sprintf(['lambda%0.3f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d'], ...
        lambda, lambda_mesh, lambda_err, nmodes, zwidth), '.', 'p'));
    outdir = fullfile(mKDir, 'measurements') ;
    
    outputdir = fullfile(tubi.dir.data, 'TubULAR_Figures', ['Figure03_' cmapname{cmapID}]) ;
    if ~exist(outputdir, 'dir')
        mkdir(outputdir)
    end
    
    % Compute or load all timepoints
    for tp = tubi.t0 + [30,60,90] 
        close all
        disp(['t = ' num2str(tp)])
        tidx = tubi.xp.tIdx(tp) ;
    
        tubi.setTime(tp) ;
        mesh = tubi.getCurrentSPCutMeshSmRSC ;
        cutMesh = tubi.getCurrentSPCutMeshSmRS ;
    
        outfn = fullfile(outputdir, sprintf('mesh_%06d.ply', tp)) ;
        plywrite(outfn, mesh.f, mesh.v)
    
        % Check for timepoint measurement on disk
        Hfn = fullfile(outdir, sprintf('HH_vertices_%06d.mat', tp))   ;
        efn = fullfile(outdir, sprintf('gdot_vertices_%06d.mat', tp)) ;
        dfn = fullfile(outdir, sprintf('divv_vertices_%06d.mat', tp)) ;
        nfn = fullfile(outdir, sprintf('veln_vertices_%06d.mat', tp)) ;
        H2vnfn = fullfile(outdir, sprintf('H2vn_vertices_%06d.mat', tp)) ;
        rfn = fullfile(outdir, sprintf('radius_vertices_%06d.mat', tp)) ;
    
        % Load timeseries measurements
        load(Hfn, 'HH_filt', 'HH_ap', 'HH_l', 'HH_r', 'HH_d', 'HH_v')
        load(efn, 'gdot_filt', 'gdot_ap', 'gdot_l', 'gdot_r', 'gdot_d', 'gdot_v')
        load(dfn, 'divv_filt', 'divv_ap', 'divv_l', 'divv_r', 'divv_d', 'divv_v')
        load(nfn, 'veln_filt', 'veln_ap', 'veln_l', 'veln_r', 'veln_d', 'veln_v') 
        load(H2vnfn, 'H2vn_filt', 'H2vn_ap', 'H2vn_l', 'H2vn_r', 'H2vn_d', 'H2vn_v') 
        load(rfn, 'radius_filt', 'radius_ap', 'radius_l', 'radius_r', 'radius_d', 'radius_v') 
    
        % Blank out the first and last ring
        gdot_filt(1,:) = gdot_filt(1,:) - divv_filt(1, :) ;
        gdot_filt(end,:) = gdot_filt(end,:) - divv_filt(end, :) ;
        divv_filt(1,:) = 0 ;
        divv_filt(end,:) = 0 ;
    
        % separate 2d/3d data
        gdot2d = gdot_filt ;
        gdot3d = gdot_filt(:, 1:nV-1);
        divv2d = divv_filt ;
        divv3d = divv_filt(:, 1:nV-1) ;
        veln2d = veln_filt ;
        veln3d = veln_filt(:, 1:nV-1) ;
        H2vn2d = H2vn_filt ;
        H2vn3d = H2vn_filt(:, 1:nV-1) ;
        
        % set(gcf, 'visible', 'off') ;
        colormap(cmap);
    
        % plot divv in 3d
        tubi.getXYZLims ;
        xyzlim = tubi.plotting.xyzlim_um ;
        trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
            divv3d, 'edgecolor', 'none')  
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        set(gcf, 'Units', 'centimeters')
        set(gcf, 'Position', [0,0,40,20])
        set(gcf, 'color', 'w')
        caxis([-climit, climit]) 
        view([0,0])
        axis off ;
        fillfig
        F = getframe() ;
        outfn = fullfile(outputdir, sprintf('divv3d_%06d.png', tp)) ;
        imwrite(F.cdata, outfn)
    
        % plot divv in 2d
        clf
        set(gcf, 'Units', 'centimeters')
        set(gcf, 'Position', [0,0,40,20])
        set(gcf, 'color', 'w')
        ss = cutMesh.u(:, 1) ;
        ssvals = tubi.a_fixed * ss / max(ss) ;
        trisurf(cutMesh.f, ssvals, ...
            cutMesh.u(:, 2), zeros(size(cutMesh.u(:,1))),...
            divv2d, 'edgecolor', 'none')
        xlim([0, tubi.a_fixed])
        ylim([0, 1]) 
        caxis([-climit, climit]) 
        axis off
         view(top)
        fillfig
        F = getframe() ;
        outfn = fullfile(outputdir, sprintf('divv2d_%06d.png', tp)) ;
        imwrite(F.cdata, outfn)
        outfn = fullfile(outputdir, sprintf('divv2d_%06d.txt', tp)) ;
        header = 'in-plane divergence of the tangential velocity field [1/min]' ;
        write_txt_with_header(outfn, divv3d, header)
    
        % plot 2Hvn in 3d
        tubi.getXYZLims ;
        xyzlim = tubi.plotting.xyzlim_um ;
        trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
            H2vn3d, 'edgecolor', 'none')  
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        set(gcf, 'Units', 'centimeters')
        set(gcf, 'Position', [0,0,40,20])
        set(gcf, 'color', 'w')
        caxis([-climit, climit]) 
        view([0,0])
        axis off ;
        fillfig
        F = getframe() ;
        outfn = fullfile(outputdir, sprintf('H2vn3d_%06d.png', tp)) ;
        imwrite(F.cdata, outfn)
    
        % plot H2vn in 2d
        clf
        set(gcf, 'Units', 'centimeters')
        set(gcf, 'Position', [0,0,40,20])
        set(gcf, 'color', 'w')
        ss = cutMesh.u(:, 1) ;
        ssvals = tubi.a_fixed * ss / max(ss) ;
        trisurf(cutMesh.f, ssvals, ...
            cutMesh.u(:, 2), zeros(size(cutMesh.u(:,1))),...
            H2vn2d, 'edgecolor', 'none')
        xlim([0, tubi.a_fixed])
        ylim([0, 1]) 
        caxis([-climit, climit]) 
        axis off
         view(top)
        fillfig
        F = getframe() ;
        outfn = fullfile(outputdir, sprintf('H2vn2d_%06d.png', tp)) ;
        imwrite(F.cdata, outfn)
        outfn = fullfile(outputdir, sprintf('H2vn2d_%06d.txt', tp)) ;
        header = 'in-plane divergence of the tangential velocity field [1/min]' ;
        write_txt_with_header(outfn, H2vn3d, header)
    
        % plot gdot in 3d
        tubi.getXYZLims ;
        xyzlim = tubi.plotting.xyzlim_um ;
        trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
            gdot3d, 'edgecolor', 'none')  
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        set(gcf, 'Units', 'centimeters')
        set(gcf, 'Position', [0,0,40,20])
        set(gcf, 'color', 'w')
        caxis([-climit, climit]) 
        view([0,0])
        axis off ;
        fillfig
        F = getframe() ;
        outfn = fullfile(outputdir, sprintf('gdot3d_%06d.png', tp)) ;
        imwrite(F.cdata, outfn)
    
        % plot gdot in 2d
        clf
        set(gcf, 'Units', 'centimeters')
        set(gcf, 'Position', [0,0,40,20])
        set(gcf, 'color', 'w')
        ss = cutMesh.u(:, 1) ;
        ssvals = tubi.a_fixed * ss / max(ss) ;
        trisurf(cutMesh.f, ssvals, ...
            cutMesh.u(:, 2), zeros(size(cutMesh.u(:,1))),...
            gdot2d, 'edgecolor', 'none')
        xlim([0, tubi.a_fixed])
        ylim([0, 1]) 
        caxis([-climit, climit]) 
        axis off
         view(top)
        fillfig
        F = getframe() ;
        outfn = fullfile(outputdir, sprintf('gdot2d_%06d.png', tp)) ;
        imwrite(F.cdata, outfn)
        outfn = fullfile(outputdir, sprintf('gdot2d_%06d.txt', tp)) ;
        header = 'in-plane divergence of the tangential velocity field [1/min]' ;
        write_txt_with_header(outfn, gdot3d, header)
    
        %%%%%%%%%%%%%%%%
    
        % plot vn in 3d
        tubi.getXYZLims ;
        xyzlim = tubi.plotting.xyzlim_um ;
        trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
            veln3d, 'edgecolor', 'none')  
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        set(gcf, 'Units', 'centimeters')
        set(gcf, 'Position', [0,0,40,20])
        set(gcf, 'color', 'w')
        caxis([-climit_veln, climit_veln]) 
        view([0,0])
        axis off ;
        fillfig
        F = getframe() ;
        outfn = fullfile(outputdir, sprintf('veln3d_%06d.png', tp)) ;
        imwrite(F.cdata, outfn)
    
        % plot veln in 2d
        clf
        set(gcf, 'Units', 'centimeters')
        set(gcf, 'Position', [0,0,40,20])
        set(gcf, 'color', 'w')
        ss = cutMesh.u(:, 1) ;
        ssvals = tubi.a_fixed * ss / max(ss) ;
        trisurf(cutMesh.f, ssvals, ...
            cutMesh.u(:, 2), zeros(size(cutMesh.u(:,1))),...
            veln2d, 'edgecolor', 'none')
        xlim([0, tubi.a_fixed])
        ylim([0, 1]) 
        caxis([-climit_veln, climit_veln]) 
        axis off
        view(top)
        fillfig
        F = getframe() ;
        outfn = fullfile(outputdir, sprintf('veln2d_%06d.png', tp)) ;
        imwrite(F.cdata, outfn)
        outfn = fullfile(outputdir, sprintf('veln2d_%06d.txt', tp)) ;
        header = 'normal velocity [um/min]' ;
        write_txt_with_header(outfn, veln3d, header)
    
    
    
        % %% Plot results
        % % % operational plotting options
        % pOptions.overwrite = overwrite_timePoints ;
        % pOptions.plot_flows = plot_flows ;
        % pOptions.plot_Hgdot = plot_Hgdot ;
        % pOptions.plot_factors = plot_factors ;
        % % parameter plotting options
        % pOptions.doubleResolution = doubleResolution; 
        % pOptions.lambda = lambda ;
        % pOptions.lambda_err = lambda_err ;
        % pOptions.lambda_mesh = lambda_mesh ;
        % pOptions.nmodes = nmodes ;
        % pOptions.zwidth = zwidth ;
        % pOptions.H2vn2d = H2vn2d ;
        % pOptions.divv2d = divv2d ;
        % pOptions.gdot2d = gdot2d ;
        % pOptions.veln2d = veln2d ;
        % pOptions.radi2d = radi2d ;
        % pOptions.H2d = H2d ;
        % pOptions.H2vn3d = H2vn3d ;
        % pOptions.divv3d = divv3d ;
        % pOptions.gdot3d = gdot3d ;
        % pOptions.veln3d = veln3d ;
        % pOptions.radi3d = radi3d ;
        % pOptions.H3d = H3d ;
        % pOptions.cutMesh = [] ;
        % pOptions.mesh = [] ;
        % pOptions.climit = climit ;
        % pOptions.climit_err = climit ;
        % pOptions.climit_veln = climit_veln ;
        % pOptions.climit_H = climit_H ;
        % tubi.plotMetricKinematicsTimePoint(tp, pOptions)
    end
end

close all
clc
disp('done')

%% Panel e-f: Show accumulated strain as magnitude
% tubi.plotPathlineStrain

clc
close all

% Default options 
overwrite = false ;
plot_kymographs = true ;
plot_strain3d = true ;
t0 = tubi.t0set() ;
t0Pathline = t0 ;
viewAngles = [-20,25] ;
clim_tre = 2 ;
clim_dev = 4 ;

strain_trace_label = '$\frac{1}{2}\mathrm{Tr} [\bf{g}^{-1}\varepsilon]$';
strain_trace_label_norm = '$||\mathrm{Tr} [\varepsilon]||$';
strain_deviator_label = ...
    '$||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\varepsilon\right]\bf{g}||$' ;
strain_deviator_label_short = ...
    '$||$Dev$\left[\varepsilon\right]||$' ;
% strain_theta_label = '$\theta_{\mathrm{Dev}[{\varepsilon}]}$'; 

% Parameter options
lambda = tubi.smoothing.lambda ;
lambda_mesh = tubi.smoothing.lambda_mesh ;
nmodes = tubi.smoothing.nmodes ;
zwidth = tubi.smoothing.zwidth ;
climit = 0.2 ;

% Unpack tubi
tubi.getXYZLims ;
xyzlim = tubi.plotting.xyzlim_um ;

% Colormap
close all
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;
bbr256 = blueblackred(256) ;
close all
pm256 = phasemap(256) ;

% load from tubi
nU = tubi.nU ;
nV = tubi.nV ;    

% Unit definitions for axis labels
unitstr = [ '[unitless]' ];
tidx0 = tubi.xp.tIdx(tubi.t0set()) ;



cmap = brewermap(256, '*RdBu') ;

% Output dir for figure panels
outputdir = fullfile(tubi.dir.data, 'TubULAR_Figures', ['Figure03ef_' cmapname{cmapID}]) ;
if ~exist(outputdir, 'dir')
    mkdir(outputdir)
end

% Plot strains in 3d
pb = tubi.getPullbackPathlines([], 'vertexPathlines3d') ;
refMesh = pb.refMesh ;
% Set the rotated scaled vertices as field "v" for
% inducedStrainPeriodicMesh() call
refMesh.v = refMesh.vrs ;
fa = doublearea(refMesh.vrs, refMesh.f) * 0.5 ;

tp = tubi.t0 + 90 ;
tubi.setTime(tp) ;

% EDIT HERE FROM ORIGINAL
outfn = fullfile(outputdir, ...
   [sprintfm('t0_%04d', t0Pathline)  sprintf('smoothed_strain_%06d.png', tp)]);

strain = tubi.getCurrentPathlineStrain(tubi.t0, 'strain') ;
tre = strain.strain.tre_sm ;
fra = strain.strain.fractionalAreaChange_sm ;
dev = strain.strain.dev_sm ;

% positions of the pathline vertices
xx = pb.vertices3d.vXrs(tidx, :, :) ;
yy = pb.vertices3d.vYrs(tidx, :, :) ;
zz = pb.vertices3d.vZrs(tidx, :, :) ;
v3d = [ xx(:), yy(:), zz(:) ] ;
mesh = struct() ;
mesh.f = refMesh.f ;
mesh.v = v3d ;
mesh.u = refMesh.u ;
mesh.pathPairs = refMesh.pathPairs ;
mesh.nU = tubi.nU ;
mesh.nV = tubi.nV ;

ax2 = subtightplot(2, 1, 1) ;
trisurf(triangulation(mesh.f, mesh.v), 'facevertexcdata', fra, 'edgecolor', 'none')
colorbar()
axis equal
view(viewAngles)
caxis(clim_tre*[-1,1])
colormap(ax2, brewermap(256, '*RdBu'))
title('\deltaA/A_0')
grid off
xlim(xyzlim(1, :))
ylim(xyzlim(2, :))
zlim(xyzlim(3, :))

ax4 = subtightplot(2, 1, 2) ;
trisurf(triangulation(mesh.f, mesh.v), 'facevertexcdata', dev, 'edgecolor', 'none')
colorbar()
axis equal
view(viewAngles)
caxis([0,clim_dev])
colormap(ax4, batlowk)
grid off
xlim(xyzlim(1, :))
ylim(xyzlim(2, :))
zlim(xyzlim(3, :))

set(gcf, 'color', 'w')
export_fig(outfn, '-nocrop', '-r600')

% Save data
dataPlotted = struct('mesh', mesh, 'fra', fra, 'dev', dev) ;
save(fullfile(outputdir, 'dataPlotted'), 'dataPlotted')

close all
clc
