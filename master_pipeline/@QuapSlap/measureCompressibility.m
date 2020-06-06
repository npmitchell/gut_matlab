function [cumerr, HHseries, divvseries, velnseries] = ...
    measureCompressibility(QS, lambda, lambda_err, climit, climit_err)
% Measure degree of incompressibility of the flow on the evolving surface
% 
% NPMitchell 2020

if nargin < 2
    lambda = 0.01 ;
end
if nargin < 3
    lambda_err = 0.01 ;
end
if nargin < 4
    climit = 0.2 ;
end
if nargin < 5
    climit_err = 0.03 ;
end
climit_veln = climit * 10 ;
climit_H = climit * 2 ;

% Unpack QS
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load vertex-based velocity measurements
vvsmMfn = fullfile(QS.dir.pivSimAvg, 'vvM_simpletimeavg.mat')  ;
tmp = load(vvsmMfn) ;
vertex_vels = tmp.vvsmM ;
vfsmMfn = fullfile(QS.dir.pivSimAvg, 'vfM_simpletimeavg.mat') ;
tmp = load(vfsmMfn) ;
face_vels = tmp.vfsmM ;

%% Load time offset for first fold, t0
QS.t0set() ;
tfold = QS.t0 ;

%% load from QS
nU = QS.nU ;
nV = QS.nV ;

%% Test incompressibility of the flow on the evolving surface
% We relate the normal velocities to the divergence / 2 * H.
tps = QS.xp.fileMeta.timePoints(1:end-1) - tfold;

% preallocate for cumulative error
cumerr = zeros(length(QS.xp.fileMeta.timePoints(1:end-1)), nU) ;
HHseries = zeros(length(QS.xp.fileMeta.timePoints(1:end-1)), nU) ;
divvseries = zeros(length(QS.xp.fileMeta.timePoints(1:end-1)), nU) ;
velnseries = zeros(length(QS.xp.fileMeta.timePoints(1:end-1)), nU) ;
for tp = QS.xp.fileMeta.timePoints(1:end-1)
    disp(['t = ' num2str(tp)])
    tidx = QS.xp.tIdx(tp) ;

    % Load current mesh
    tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp)) ;
    mesh = tmp.spcutMeshSmRSC ;
    clearvars tmp
    
    % Load cutMesh
    tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
    cutMesh = tmp.spcutMeshSmRS ;
    clearvars tmp
    
    % Compute mean curvature
    DEC = DiscreteExteriorCalculus(mesh.f, mesh.v) ;
    HH = sum(mesh.vn .* DEC.laplacian(mesh.v), 2) * 0.5 ;
    fvel = squeeze(face_vels(tidx, :, :)) ;
    divv = DEC.divergence(fvel) ;
    divv = laplacian_smooth(mesh.v, mesh.f, 'cotan', [],...
                            lambda, 'implicit', divv) ;
    
    % veln = divv / (2H) ; 
    % veln_pred = DEC.divergence(fvel) ./ (2.0 * H) ;
    % veln is normal velocity field (velocities along vertex normals)
    % mesh.vn are vertex normals
    
    % Smooth the velocities in space using gptoolbox
    vx = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 1)) ;
    vy = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 2)) ;
    vz = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 3)) ;
    vxs = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
        lambda, 'implicit', vx') ;
    vys = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
        lambda, 'implicit', vy') ;
    vzs = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
        lambda, 'implicit', vz') ;
    
    % Actual normal velocity
    veln = sum(cat(2, vxs, vys, vzs) .* mesh.vn, 2) ;
    
    % Predict the divergence
    div_pred = 2 * HH .* veln ;
    
    % Extend to have another row for 2d map
    divv2d = divv ;
    divv2d(nU*(nV-1) + 1:nU*nV) = divv(1:nU) ;
    HH2d = HH ;
    HH2d(nU*(nV-1)+1:(nU*nV)) = HH(1:nU) ;
    veln2d = veln ;
    veln2d(nU*(nV-1)+1:(nU*nV)) = veln(1:nU) ;
    % div_pred2d = div_pred ;
    % div_pred2d(nU*(nV-1)+1:(nU*nV)) = div_pred(1:nU) ;
    div_pred2d = 2 * HH2d .* veln2d ;
    
    % The difference
    if lambda_err > 0 
        err3d = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
            lambda_err, 'implicit', divv-div_pred) ;
        % expand error to 2d map
        err2d = err3d ;
        err2d(nU*(nV-1)+1:nU*nV) = err3d(1:nU) ;
    else
        err2d = divv2d - div_pred2d ;
        err3d = divv - div_pred ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the prediction, the measurement, and the difference
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    colors2d = {div_pred2d, divv2d, err2d} ; 
    colors3d = {div_pred, divv, err3d} ;
    dimDirs = {QS.dir.compressibility2d, QS.dir.compressibility3d} ;
    set(gcf, 'visible', 'off') ;
    colormap bwr
    plot_flows = false ;
    if plot_flows
        for dim = 2:3    
            dimDir = dimDirs{dim-1} ;
            for row = 1:2  % row
                for col = 1:3  % column
                    % create panel
                    subplot(3, length(colors2d), col + (row - 1) * 3)

                    % If 2d, plot in pullback space
                    % if 3d, plot in embedding space
                    if dim == 2 && row == 2
                        ss = cutMesh.u(:, 1) ;
                        ssvals = QS.a_fixed * ss / max(ss) ;
                        trisurf(cutMesh.f, ssvals, ...
                            cutMesh.u(:, 2), zeros(size(cutMesh.u(:,1))),...
                            colors2d{col}, 'edgecolor', 'none')
                        xlim([0, QS.a_fixed])
                        ylim([0, 1]) 
                        caxis([-climit, climit]) 
                        axis off
                    else
                        trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
                            colors3d{col}, 'edgecolor', 'none')
                        axis equal
                        axis off
                        caxis([-climit, climit]) 

                        % xlabel('AP position, [$\mu$m]', 'Interpreter', 'Latex')
                        % ylabel('lateral position, [$\mu$m]', 'Interpreter', 'Latex')
                        % zlabel('DV position, [$\mu$m]', 'Interpreter', 'Latex')
                        xlim(xyzlim(1, :))
                        ylim(xyzlim(2, :))
                        zlim(xyzlim(3, :))
                    end

                    % Set title and colorbars
                    if col == 1 && row == 1
                        % title(['$v_n 2 H$', ...
                        %     ', $t=$' sprintf('%03d', tp - tfold)], ...
                        %     'Interpreter', 'Latex')
                        view(0, 0)
                        caxis([-climit, climit]) 
                    elseif col == 2 && row == 1
                        % title(['$\nabla \cdot \bf{v}_\parallel$', ...
                        %     ', $t=$' sprintf('%03d', tp - tfold)], ...
                        %     'Interpreter', 'Latex')
                        title(['$t=$' sprintf('%03d', tp - tfold)], ...
                             'Interpreter', 'Latex')
                        view(0, 0)
                        caxis([-climit, climit]) 
                    elseif col == 3 && row == 1
                        % title(['$\textrm{Tr}[g^{-1} \dot{g}]$', ...
                        %     ', $t=$' sprintf('%03d', tp - tfold)], ...
                        %     'Interpreter', 'Latex')
                        view(0, 0)
                        caxis([-climit_err, climit_err]) 
                    elseif col == 1 && row == 2
                        view(0, 270)
                        cb = colorbar('south') ;
                        % set(cb, 'position',[.17 .1 .2 .03])
                        set(cb, 'position',[.165 .2 .15 .03])     
                        ylabel(cb, '$v_n 2H$', 'Interpreter', 'Latex')                    
                        caxis([-climit, climit]) 
                    elseif col == 2 && row == 2
                        view(0, 270)
                        cb = colorbar('south') ;
                        % set(cb, 'position',[.63 .1 .2 .03])
                        set(cb, 'position',[.445 .2 .15 .03])
                        ylabel(cb, '$\nabla \cdot \bf{v}_\parallel$', ...
                            'Interpreter', 'Latex')
                        caxis([-climit, climit]) 
                    elseif col == 3 && row == 2
                        view(0, 270)
                        cb = colorbar('south') ;
                        set(cb, 'position',[.725 .2 .15 .03], ...
                            'XTick', [-climit_err, 0, climit_err])
                        ylabel(cb, '$\textrm{Tr}[g^{-1} \dot{g}]$', ...
                            'Interpreter', 'Latex')
                        % xticks(cb, [-climit_err, climit_err]) 
                        caxis([-climit_err, climit_err]) 
                   end
                end
            end
            % check it
            % set(gcf, 'visible', 'on')
            % error('here')

            fn = fullfile(dimDir, ...
                sprintf('incompr_%1dd_%06d.png', dim, tp)) ;
            disp(['saving ', fn])
            export_fig(fn, '-png', '-nocrop', '-r200') 

        end
    end
    
    %% Plot the residual with separate factors of vn, H, div, and gdot
    plot_Hgdot = true ;
    if plot_Hgdot
        % Plot mean curvature, normal velocity, divv and error
        colors2d = {veln2d, HH2d, divv2d, err2d} ; 
        colors3d = {veln, HH, divv, err3d} ; 
        climits = [climit_veln, climit_H, climit, climit_err] ;
        for row = 1:2  % row
            for col = 1:4  % column
                % create panel
                subplot(3, length(colors2d), col + (row - 1) * 4)

                % plot in pullback space
                if row == 2
                    ss = cutMesh.u(:, 1) ;
                    ssvals = QS.a_fixed * ss / max(ss) ;
                    trisurf(cutMesh.f, ssvals, ...
                        cutMesh.u(:, 2), zeros(size(cutMesh.u(:,1))),...
                        colors2d{col}, 'edgecolor', 'none')
                    xlim([0, QS.a_fixed])
                    ylim([0, 1]) 
                    caxis([-climits(col), climits(col) ]) 
                    axis off
                else
                    trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
                        colors3d{col}, 'edgecolor', 'none')
                    axis equal
                    axis off
                    caxis([-climits(col), climits(col)])

                    % xlabel('AP position, [$\mu$m]', 'Interpreter', 'Latex')
                    % ylabel('lateral position, [$\mu$m]', 'Interpreter', 'Latex')
                    % zlabel('DV position, [$\mu$m]', 'Interpreter', 'Latex')
                    xlim(xyzlim(1, :))
                    ylim(xyzlim(2, :))
                    zlim(xyzlim(3, :))
                end

                % Set title and colorbars
                if col == 2 && row == 1
                    % 3d view with title
                    title(['$t=$' sprintf('%03d', tp - tfold)], ...
                         'Interpreter', 'Latex')
                    view(0, 0)
                elseif row == 1
                    % 3d view
                    view(0, 0)
                elseif col == 1 && row == 2
                    view(0, 270)
                    cb = colorbar('south') ;
                    set(cb, 'position',[.15 .2 .12 .03], 'Xtick', ...
                        [-climits(col), 0, climits(col)]) 
                    ylabel(cb, '$v_n$', 'Interpreter', 'Latex')
                elseif col == 2 && row == 2
                    view(0, 270)
                    cb = colorbar('south') ;
                    set(cb, 'position',[.355 .2 .12 .03], 'Xtick', ...
                        [-climits(col), 0, climits(col)])
                    ylabel(cb, '$H$', ...
                        'Interpreter', 'Latex')
                elseif col == 3 && row == 2
                    view(0, 270)
                    cb = colorbar('south') ;
                    set(cb, 'position',[.565 .2 .12 .03], 'Xtick', ...
                        [-climits(col), 0, climits(col)])
                    ylabel(cb, '$\nabla \cdot \mathbf{v}$', ...
                        'Interpreter', 'Latex')
                elseif col == 4 && row == 2
                    view(0, 270)
                    cb = colorbar('south') ;
                    set(cb, 'position',[.765 .2 .12 .03], 'Xtick', ...
                        [-climits(col), 0, climits(col)])
                    ylabel(cb, '$\textrm{Tr}[g^{-1} \dot{g}]$', ...
                        'Interpreter', 'Latex')
               end
               caxis([-climits(col), climits(col)])
            end
        end

        fn = fullfile(QS.dir.compressibility, 'gdot_vs_H', ...
            sprintf('gdot_vn_H_2d_%06d.png', tp)) ;
        disp(['saving ', fn])
        export_fig(fn, '-png', '-nocrop', '-r200')    

    end
    
    % Plot the factors separately of vn, 2H, and vn*2H
    plot_factors = false ;
    if plot_factors
        % Plot mean curvature, normal velocity, and product
        colors2d = {veln2d, HH2d, div_pred2d} ; 
        colors3d = {veln, HH, div_pred} ; 
        climits = [climit_veln, climit_H, climit] ;
        for row = 1:2  % row
            for col = 1:3  % column
                % create panel
                subplot(3, length(colors2d), col + (row - 1) * 3)

                % plot in pullback space
                if row == 2
                    ss = cutMesh.u(:, 1) ;
                    ssvals = QS.a_fixed * ss / max(ss) ;
                    trisurf(cutMesh.f, ssvals, ...
                        cutMesh.u(:, 2), zeros(size(cutMesh.u(:,1))),...
                        colors2d{col}, 'edgecolor', 'none')
                    xlim([0, QS.a_fixed])
                    ylim([0, 1]) 
                    caxis([-climits(col), climits(col) ]) 
                    axis off
                else
                    trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
                        colors3d{col}, 'edgecolor', 'none')
                    axis equal
                    axis off
                    caxis([-climits(col), climits(col)])

                    % xlabel('AP position, [$\mu$m]', 'Interpreter', 'Latex')
                    % ylabel('lateral position, [$\mu$m]', 'Interpreter', 'Latex')
                    % zlabel('DV position, [$\mu$m]', 'Interpreter', 'Latex')
                    xlim(xyzlim(1, :))
                    ylim(xyzlim(2, :))
                    zlim(xyzlim(3, :))
                end

                % Set title and colorbars
                if col == 2 && row == 1
                    % 3d view with title
                    title(['$t=$' sprintf('%03d', tp - tfold)], ...
                         'Interpreter', 'Latex')
                    view(0, 0)
                elseif row == 1
                    % 3d view
                    view(0, 0)
                elseif col == 1 && row == 2
                    view(0, 270)
                    cb = colorbar('south') ;
                    set(cb, 'position',[.165 .2 .15 .03])     
                    ylabel(cb, '$v_n$', 'Interpreter', 'Latex')
                elseif col == 2 && row == 2
                    view(0, 270)
                    cb = colorbar('south') ;
                    set(cb, 'position',[.445 .2 .15 .03])
                    ylabel(cb, '$H$', ...
                        'Interpreter', 'Latex')
                elseif col == 3 && row == 2
                    view(0, 270)
                    cb = colorbar('south') ;
                    set(cb, 'position',[.725 .2 .15 .03])
                    ylabel(cb, '$v_n 2H$', ...
                        'Interpreter', 'Latex') 
               end
               caxis([-climits(col), climits(col)])
            end
        end

        fn = fullfile(QS.dir.compressibility, 'factors', ...
            sprintf('factors_2d_%06d.png', tp)) ;
        disp(['saving ', fn])
        export_fig(fn, '-png', '-nocrop', '-r200')    

    end
    
    % cumulative error
    HHseries(tidx, :) = mean(reshape(HH2d, [nU,nV]), 2) ;
    cumerr(tidx, :) = mean(reshape(err2d, [nU,nV]), 2) ;
    divvseries(tidx, :) = mean(reshape(divv2d, [nU,nV]), 2) ;
    velnseries(tidx, :) = mean(reshape(veln2d, [nU,nV]), 2) ;

end

% write lambdas to disk
lambs = [lambda, lambda_err] ;
header = ['laplacian smoothing parameters lambda for vertex ', ...
    'velocities and div(v)-vn*2*H'] ;
filename = fullfile(QS.dir.compressibility, 'lambdas.txt') ;
write_txt_with_header(filename, lambs, header)

%% Now plot cumulative error as function of x over time
close all
colormap bwr
imagesc((1:nU)/nU, tps, cumerr)
caxis([-climit_err, climit_err])
% Add folds to plot
folds = load(QS.fileName.fold) ;
hold on;
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;
fons1 = fons(1) ;
fons2 = fons(2) ;
fons3 = fons(3) ;
plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

% title and save
title('$\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$',...
    'Interpreter', 'Latex')
ylabel('time [min]', 'Interpreter', 'Latex')
xlabel('AP position [$x/L$]', 'Interpreter', 'Latex')
cb = colorbar() ;
ylabel(cb, '$\textrm{Tr}[g^{-1}\dot{g}]$', 'Interpreter', 'Latex')  
fn = fullfile(QS.dir.compressibility, ...
    sprintf('cumerr.png')) ;
disp(['saving ', fn])
export_fig(fn, '-png', '-nocrop', '-r200')   

% Zoom in on small values
caxis([-climit_err/3, climit_err/3])
fn = fullfile(QS.dir.compressibility, ...
    sprintf('cumerr_zoom.png')) ;
disp(['saving ', fn])
export_fig(fn, '-png', '-nocrop', '-r200')   
% Zoom in on early times
ylim([min(tps), max(fons) + 10])
caxis([-climit_err/3, climit_err/3])
fn = fullfile(QS.dir.compressibility, ...
    sprintf('cumerr_zoom_early.png')) ;
disp(['saving ', fn])
export_fig(fn, '-png', '-nocrop', '-r200')   

%% Now plot timeseries of mean curvature as function of x over time
close all
colormap bwr
imagesc((1:nU)/nU, tps, HHseries)
caxis([-climit_H, climit_H])
% Add folds to plot
hold on;
plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

% title and save
title('mean curvature, $H$',...
    'Interpreter', 'Latex')
ylabel('time [min]', 'Interpreter', 'Latex')
xlabel('AP position [$x/L$]', 'Interpreter', 'Latex')
cb = colorbar() ;
ylabel(cb, 'mean curvature, $H$', 'Interpreter', 'Latex')  
fn = fullfile(QS.dir.compressibility, ...
    sprintf('HHseries.png')) ;
disp(['saving ', fn])
export_fig(fn, '-png', '-nocrop', '-r200')   

caxis([-climit_H/4, climit_H/4])
fn = fullfile(QS.dir.compressibility, ...
    sprintf('HHseries_zoom.png')) ;
disp(['saving ', fn])
export_fig(fn, '-png', '-nocrop', '-r200')   
 
ylim([min(tps), max(fons) + 10])
fn = fullfile(QS.dir.compressibility, ...
    sprintf('HHseries_zoom_early.png')) ;
disp(['saving ', fn])
export_fig(fn, '-png', '-nocrop', '-r200')   


%% Now plot timeseries of div(v) as function of x over time
close all
colormap bwr
imagesc((1:nU)/nU, tps, divvseries)
caxis([-climit, climit])
% Add folds to plot
hold on;
plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

% title and save
title('divergence of flow, $\nabla \cdot \mathbf{v}$',...
    'Interpreter', 'Latex')
ylabel('time [min]', 'Interpreter', 'Latex')
xlabel('AP position [$x/L$]', 'Interpreter', 'Latex')
cb = colorbar() ;
ylabel(cb, '$\nabla \cdot \mathbf{v}$', 'Interpreter', 'Latex')  
fn = fullfile(QS.dir.compressibility, ...
    sprintf('divvseries.png')) ;
disp(['saving ', fn])
export_fig(fn, '-png', '-nocrop', '-r200')   

caxis([-climit/4, climit/4])
fn = fullfile(QS.dir.compressibility, ...
    sprintf('divvseries_zoom.png')) ;
disp(['saving ', fn])
export_fig(fn, '-png', '-nocrop', '-r200')   
 
ylim([min(tps), max(fons) + 10])
fn = fullfile(QS.dir.compressibility, ...
    sprintf('divvseries_zoom_early.png')) ;
disp(['saving ', fn])
export_fig(fn, '-png', '-nocrop', '-r200')   


%% Now plot timeseries of veln as function of x over time
close all
colormap bwr
imagesc((1:nU)/nU, tps, velnseries)
caxis([-climit_veln, climit_veln])
% Add folds to plot
hold on;
plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

% title and save
title('normal velocity, $v_n$',...
    'Interpreter', 'Latex')
ylabel('time [min]', 'Interpreter', 'Latex')
xlabel('AP position [$x/L$]', 'Interpreter', 'Latex')
cb = colorbar() ;
ylabel(cb, 'normal velocity, $v_n$', 'Interpreter', 'Latex')  
fn = fullfile(QS.dir.compressibility, ...
    sprintf('velnseries.png')) ;
disp(['saving ', fn])
export_fig(fn, '-png', '-nocrop', '-r200')   

caxis([-climit_veln/4, climit_veln/4])
fn = fullfile(QS.dir.compressibility, ...
    sprintf('velnseries_zoom.png')) ;
disp(['saving ', fn])
export_fig(fn, '-png', '-nocrop', '-r200')   
 
ylim([min(tps), max(fons) + 10])
fn = fullfile(QS.dir.compressibility, ...
    sprintf('velnseries_zoom_early.png')) ;
disp(['saving ', fn])
export_fig(fn, '-png', '-nocrop', '-r200')   

disp('done')

