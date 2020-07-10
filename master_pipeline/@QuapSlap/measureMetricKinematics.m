function [gdot_apM, HH_apM, divv_apM, veln_apM] = ...
    measureMetricKinematics(QS, options)
% Measure degree of incompressibility of the flow on the evolving surface
% 
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%   timePoints : numeric 1D array
%       the timepoints to consider for the measurement. For ex, could
%       choose subset of the QS experiment timePoints
%   alphaVal : float
%       the opacity of the heatmap to overlay
%   invertImage : bool
%       invert the data pullback under the velocity map
%
% Returns 
% -------
%
% Saves to disk
% -------------
% 
% NPMitchell 2020

%% Default options 
overwrite = false ;
plot_Hgdot = true ;
plot_flows = true ;
plot_factors = true ;
plot_kymographs = true ;
plot_correlations = true ;
plot_gdot_decomp = true ;
lambda = 0.01 ; 
% by default, lambda_mesh = lambda, whether defined here or in options
lambda_err = 0.01 ;
climit = 0.1 ;
climit_err = 0.1 ;
climit_veln = climit * 10 ;
climit_H = climit * 2 ;

%% Unpack options & assign defaults
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'plot_Hgdot')
    plot_Hgdot = options.plot_Hgdot ;
end
if isfield(options, 'plot_flows')
    plot_flows = options.plot_flows ;
end
if isfield(options, 'plot_factors')
    plot_factors = options.plot_factors ;
end
if isfield(options, 'plot_kymographs')
    plot_kymographs = options.plot_kymographs ;
end
if isfield(options, 'plot_correlations')
    plot_correlations = options.plot_correlations ;
end
if isfield(options, 'plot_gdot_decomp')
    plot_gdot_decomp = options.plot_gdot_decomp ;
end
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_err')
    lambda_err = options.lambda_err ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
else
    % default lambda_mesh is equal to lambda 
    lambda_mesh = lambda ;
end
if isfield(options, 'climit')
    climit = options.climit ;
end
if isfield(options, 'climit_err')
    climit_err = options.climit_err ;
end
if isfield(options, 'climit_veln')
    climit_veln = options.climit_veln ;
end
if isfield(options, 'climit_H')
    climit_H = options.climit_H ;
end

%% Unpack QS
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;
mKDir = QS.dir.metricKinematics ;
folds = load(QS.fileName.fold) ;
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;

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
ntps = length(QS.xp.fileMeta.timePoints(1:end-1)) ;
HH_apM   = zeros(ntps, nU) ;   % dv averaged
divv_apM = zeros(ntps, nU) ;
veln_apM = zeros(ntps, nU) ;
gdot_apM = zeros(ntps, nU) ;
HH_lM   = zeros(ntps, nU) ;    % left averaged
divv_lM = zeros(ntps, nU) ;
veln_lM = zeros(ntps, nU) ;
gdot_lM = zeros(ntps, nU) ;
HH_rM   = zeros(ntps, nU) ;    % right averaged
divv_rM = zeros(ntps, nU) ;
veln_rM = zeros(ntps, nU) ;
gdot_rM = zeros(ntps, nU) ;
HH_dM   = zeros(ntps, nU) ;    % dorsal averaged
divv_dM = zeros(ntps, nU) ;
veln_dM = zeros(ntps, nU) ;
gdot_dM = zeros(ntps, nU) ;
HH_vM   = zeros(ntps, nU) ;    % ventral averaged
divv_vM = zeros(ntps, nU) ;
veln_vM = zeros(ntps, nU) ;
gdot_vM = zeros(ntps, nU) ;

% Build timepoint list so that we first do every 10, then fill in details
lastIdx = length(QS.xp.fileMeta.timePoints) - 1 ;
coarseIdx = 1:10:lastIdx ;
fineIdx = setdiff(1:lastIdx, coarseIdx) ;
allIdx = [coarseIdx, fineIdx ] ;
tp2do = QS.xp.fileMeta.timePoints(allIdx) ;

% Output directory is inside metricKinematics dir
outdir = fullfile(mKDir, 'measurements') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Unit definitions for axis labels
unitstr = [ '[1/' QS.timeunits ']' ];
Hunitstr = [ '[1/' QS.spaceunits ']' ];
vunitstr = [ '[' QS.spaceunits '/' QS.timeunits ']' ];
    
% Compute or load all timepoints
for tp = tp2do
    close all
    disp(['t = ' num2str(tp)])
    titlestr = ['$t=$' sprintf('%03d', tp - tfold) ' ' QS.timeunits ] ;
    tidx = QS.xp.tIdx(tp) ;

    % Check for timepoint measurement on disk
    Hfn = fullfile(outdir, sprintf('HH_series_%06d.mat', tp))   ;
    efn = fullfile(outdir, sprintf('gdot_series_%06d.mat', tp)) ;
    dfn = fullfile(outdir, sprintf('divv_series_%06d.mat', tp)) ;
    nfn = fullfile(outdir, sprintf('veln_series_%06d.mat', tp)) ;
    H2vnfn = fullfile(outdir, sprintf('H2vn_series_%06d.mat', tp)) ;

    redo_comp = overwrite || ~exist(Hfn, 'file') || ~exist(efn, 'file') ;
    redo_comp = redo_comp || ~exist(dfn, 'file') || ~exist(dfn, 'file') ;

    if redo_comp
        tic 
        % Load current mesh
        tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp)) ;
        mesh = tmp.spcutMeshSmRSC ;

        % Load cutMesh
        tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
        cutMesh = tmp.spcutMeshSmRS ;
        clearvars tmp

        % Compute mean curvature
        % Smooth the mesh with lambda_mesh
        if lambda_mesh > 0 
            disp('smoothing mesh vertices before computations')
            tri = triangulation(mesh.f, mesh.v) ;
            fbndy = tri.freeBoundary ;
            fbndy = fbndy(:, 1) ;
            mesh.v = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
                lambda_mesh, 'implicit', mesh.v) ;
        end
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


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This should not be necessary, since vertex_vels should already
        % account for APDV framing! Debug this issue earlier
        if QS.flipy
            vys = -vys ;
        end
        % DEBUG THE FLIP
        % vv = squeeze(vertex_vels(tidx, 1:nU*(nV-1), :)) ;
        % trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
        %     vv(:, 2), 'edgecolor', 'none', 'facealpha', 0.6)
        % pause(1);
        % hold on;
        % quiver3(mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), vv(:, 1), vv(:, 2), vv(:, 3))
        % tmp2 = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp+1)) ;
        % m2 = tmp2.spcutMeshSmRSC ;
        % trisurf(m2.f, m2.v(:, 1), m2.v(:, 2), m2.v(:, 3), ...
        %     vv(:, 2), 'edgecolor', 'k', 'facealpha', 0.6)
        % axis equal
        % figure ; 
        % plot(mesh.v(:, 2) - m2.v(:, 2), vv(:, 2), '.')
        % clearvars tmp

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        dimDirs = {QS.dir.metricKinematics2d, QS.dir.metricKinematics3d} ;
        set(gcf, 'visible', 'off') ;
        colormap bwr
        % Check if files already exist
        fn2 = fullfile(dimDirs{1}, sprintf('incompr_%1dd_%06d.png', 2, tp)) ;
        fn3 = fullfile(dimDirs{2}, sprintf('incompr_%1dd_%06d.png', 3, tp)) ;
        fns = {fn2, fn3} ;
        redo_prediction = ~exist(fn2, 'file') || ...
                          ~exist(fn3, 'file') || overwrite ;
        if plot_flows && redo_prediction
            for dim = 2:3    
                % Create all panels in 2d or 3d
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
                            title(titlestr, 'Interpreter', 'Latex')
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
                            ylabel(cb, ['$v_n 2H$ ' unitstr], ...
                                'Interpreter', 'Latex')                    
                            caxis([-climit, climit]) 
                        elseif col == 2 && row == 2
                            view(0, 270)
                            cb = colorbar('south') ;
                            % set(cb, 'position',[.63 .1 .2 .03])
                            set(cb, 'position',[.445 .2 .15 .03])
                            ylabel(cb, ...
                                ['$\nabla \cdot \bf{v}_\parallel$ ' unitstr ], ...
                                'Interpreter', 'Latex')
                            caxis([-climit, climit]) 
                        elseif col == 3 && row == 2
                            view(0, 270)
                            cb = colorbar('south') ;
                            set(cb, 'position',[.725 .2 .15 .03], ...
                                'XTick', [-climit_err, 0, climit_err])
                            ylabel(cb, ...
                                ['$\textrm{Tr}[g^{-1} \dot{g}]$ ' unitstr], ...
                                'Interpreter', 'Latex')
                            % xticks(cb, [-climit_err, climit_err]) 
                            caxis([-climit_err, climit_err]) 
                       end
                    end
                end

                % Save the plot
                fn = fns{dim-1} ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200') 
            end
        end

        %% Plot the residual with separate factors of vn, H, div, and gdot
        gdot_H_dir = fullfile(mKDir, 'gdot_vs_H') ;
        fn = fullfile(gdot_H_dir, sprintf('gdot_vn_H_2d_%06d.png', tp)) ;
        if plot_Hgdot && (~exist(fn, 'file') || overwrite)
            % Plot mean curvature, normal velocity, divv and error
            colors2d = {veln2d, HH2d, divv2d, err2d} ; 
            colors3d = {veln, HH, divv, err3d} ; 
            climits = [climit_veln, climit_H, climit, climit_err] ;
            if ~exist(gdot_H_dir, 'dir')
                mkdir(gdot_H_dir) ;
            end
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
                        title(titlestr, 'Interpreter', 'Latex')
                        view(0, 0)
                    elseif row == 1
                        % 3d view
                        view(0, 0)
                    elseif col == 1 && row == 2
                        view(0, 270)
                        cb = colorbar('south') ;
                        set(cb, 'position',[.15 .2 .12 .03], 'Xtick', ...
                            [-climits(col), 0, climits(col)]) 
                        ylabel(cb, ['$v_n$ ' vunitstr], 'Interpreter', 'Latex')
                    elseif col == 2 && row == 2
                        view(0, 270)
                        cb = colorbar('south') ;
                        set(cb, 'position',[.355 .2 .12 .03], 'Xtick', ...
                            [-climits(col), 0, climits(col)])
                        ylabel(cb, ['$H$ ' Hunitstr], ...
                            'Interpreter', 'Latex')
                    elseif col == 3 && row == 2
                        view(0, 270)
                        cb = colorbar('south') ;
                        set(cb, 'position',[.565 .2 .12 .03], 'Xtick', ...
                            [-climits(col), 0, climits(col)])
                        ylabel(cb, ['$\nabla \cdot \mathbf{v}$ ' unitstr], ...
                            'Interpreter', 'Latex')
                    elseif col == 4 && row == 2
                        view(0, 270)
                        cb = colorbar('south') ;
                        set(cb, 'position',[.765 .2 .12 .03], 'Xtick', ...
                            [-climits(col), 0, climits(col)])
                        ylabel(cb, ['$\textrm{Tr}[g^{-1} \dot{g}]$ ' unitstr], ...
                            'Interpreter', 'Latex')
                   end
                   caxis([-climits(col), climits(col)])
                end
            end
            disp(['saving ', fn])
            export_fig(fn, '-png', '-nocrop', '-r200')    

        end

        %% Plot the factors separately of vn, 2H, and vn*2H
        factorsDir = fullfile(mKDir, 'factors') ;
        fn = fullfile(factorsDir, sprintf('factors_2d_%06d.png', tp)) ;
        if plot_factors && (~exist(fn, 'file') || overwrite)
            % ensure output directory
            if ~exist(factorsDir, 'dir')
                mkdir(factorsDir)
            end

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
                        title(titlestr, 'Interpreter', 'Latex')
                        view(0, 0)
                    elseif row == 1
                        % 3d view
                        view(0, 0)
                    elseif col == 1 && row == 2
                        view(0, 270)
                        cb = colorbar('south') ;
                        set(cb, 'position',[.165 .2 .15 .03])     
                        ylabel(cb, ['$v_n$ ' vunitstr], 'Interpreter', 'Latex')
                    elseif col == 2 && row == 2
                        view(0, 270)
                        cb = colorbar('south') ;
                        set(cb, 'position',[.445 .2 .15 .03])
                        ylabel(cb, ['$H$ ' Hunitstr ], ...
                            'Interpreter', 'Latex')
                    elseif col == 3 && row == 2
                        view(0, 270)
                        cb = colorbar('south') ;
                        set(cb, 'position',[.725 .2 .15 .03])
                        ylabel(cb, ['$v_n 2H$ ' unitstr], ...
                            'Interpreter', 'Latex') 
                   end
                   caxis([-climits(col), climits(col)])
                end
            end

            % Save figure
            disp(['saving ', fn])
            export_fig(fn, '-png', '-nocrop', '-r200')    
        end

        % Store data on disk
        HH = reshape(HH2d, [nU,nV]) ;
        gdot = reshape(err2d, [nU,nV]) ;
        divv = reshape(divv2d, [nU,nV]) ;
        veln = reshape(veln2d, [nU,nV]) ;
        H2vn = reshape(HH2d .* veln2d, [nU, nV]) ;
        
        % Average along DV
        HH_ap = mean(HH, 2) ;
        gdot_ap = mean(gdot, 2) ;
        divv_ap = mean(divv, 2) ;
        veln_ap = mean(veln, 2) ;
        H2vn_ap = mean(H2vn, 2) ;
        
        % quarter bounds
        q0 = round(nV * 0.125) ;
        q1 = round(nV * 0.375) ;
        q2 = round(nV * 0.625) ;
        q3 = round(nV * 0.875) ;
        left = q0:q1 ;
        ventral = q1:q2 ;
        right = q2:q3 ;
        dorsal = [q3:nV, 1:q1] ;
        
        % left quarter
        HH_l = mean(HH(:, left), 2) ;
        gdot_l = mean(gdot(:, left), 2) ;
        divv_l = mean(divv(:, left), 2) ;
        veln_l = mean(veln(:, left), 2) ;
        H2vn_l = mean(H2vn(:, left), 2) ;
        
        % right quarter
        HH_r = mean(HH(:, right), 2) ;
        gdot_r = mean(gdot(:, right), 2) ;
        divv_r = mean(divv(:, right), 2) ;
        veln_r = mean(veln(:, right), 2) ;
        H2vn_r = mean(H2vn(:, right), 2) ;
        
        % dorsal quarter
        HH_d = mean(HH(:, dorsal), 2) ;
        gdot_d = mean(gdot(:, dorsal), 2) ;
        divv_d = mean(divv(:, dorsal), 2) ;
        veln_d = mean(veln(:, dorsal), 2) ;
        H2vn_d = mean(H2vn(:, dorsal), 2) ;
        
        % ventral quarter
        HH_v = mean(HH(:, ventral), 2) ;
        gdot_v = mean(gdot(:, ventral), 2) ;
        divv_v = mean(divv(:, ventral), 2) ;
        veln_v = mean(veln(:, ventral), 2) ;
        H2vn_v = mean(H2vn(:, ventral), 2) ;
        
        %% Save timeseries measurements
        save(Hfn, 'HH', 'HH_ap', 'HH_l', 'HH_r', 'HH_d', 'HH_v')
        save(efn, 'gdot', 'gdot_ap', 'gdot_l', 'gdot_r', 'gdot_d', 'gdot_v')
        save(dfn, 'divv', 'divv_ap', 'divv_l', 'divv_r', 'divv_d', 'divv_v')
        save(nfn, 'veln', 'veln_ap', 'veln_l', 'veln_r', 'veln_d', 'veln_v') 
        save(H2vnfn, 'H2vn', 'H2vn_ap', 'H2vn_l', 'H2vn_r', 'H2vn_d', 'H2vn_v')
        
        %% Save lambdas used to disk
        lambs = [lambda, lambda_mesh, lambda_err] ;
        header = ['laplacian smoothing parameters: lambda for ', ...
            'velocity & div(v), lambda for mesh, lambda for residual (gdot)'] ;
        filename = fullfile(mKDir, 'lambdas.txt') ;
        write_txt_with_header(filename, lambs, header)
        toc
    else
        % Load timeseries measurements
        load(Hfn, 'HH_ap', 'HH_l', 'HH_r', 'HH_d', 'HH_v')
        load(efn, 'gdot_ap', 'gdot_l', 'gdot_r', 'gdot_d', 'gdot_v')
        load(dfn, 'divv_ap', 'divv_l', 'divv_r', 'divv_d', 'divv_v')
        load(nfn, 'veln_ap', 'veln_l', 'veln_r', 'veln_d', 'veln_v') 
        load(H2vnfn, 'H2vn_ap', 'H2vn_l', 'H2vn_r', 'H2vn_d', 'H2vn_v') 
    end
    
    %% Store in matrices
    % dv averaged
    HH_apM(tidx, :) = HH_ap ;
    gdot_apM(tidx, :) = gdot_ap ;
    divv_apM(tidx, :) = divv_ap ;
    veln_apM(tidx, :) = veln_ap ;
    H2vn_apM(tidx, :) = H2vn_ap ;

    % left quarter
    HH_lM(tidx, :) = HH_l ;
    gdot_lM(tidx, :) = gdot_l ;
    divv_lM(tidx, :) = divv_l ;
    veln_lM(tidx, :) = veln_l ;
    H2vn_lM(tidx, :) = H2vn_l ;

    % right quarter
    HH_rM(tidx, :) = HH_r ;
    gdot_rM(tidx, :) = gdot_r ;
    divv_rM(tidx, :) = divv_r ;
    veln_rM(tidx, :) = veln_r ;
    H2vn_rM(tidx, :) = H2vn_r ;

    % dorsal quarter
    HH_dM(tidx, :) = HH_d ;
    gdot_dM(tidx, :) = gdot_d ;
    divv_dM(tidx, :) = divv_d ;
    veln_dM(tidx, :) = veln_d ;
    H2vn_dM(tidx, :) = H2vn_d ;

    % ventral quarter
    HH_vM(tidx, :) = HH_v ;
    gdot_vM(tidx, :) = gdot_v ;
    divv_vM(tidx, :) = divv_v ;
    veln_vM(tidx, :) = veln_v ;
    H2vn_vM(tidx, :) = H2vn_v ;
end

%% Store kymograph data in cell array
HHsK = {HH_apM, HH_lM, HH_rM, HH_dM, HH_vM} ;
gdotsK = {gdot_apM, gdot_lM, gdot_rM, gdot_dM, gdot_vM} ;
divvsK = {divv_apM, divv_lM, divv_rM, divv_dM, divv_vM} ;
velnsK = {veln_apM, veln_lM, veln_rM, veln_dM, veln_vM} ;
H2vnsK = {H2vn_apM, H2vn_lM, H2vn_rM, H2vn_dM, H2vn_vM} ;

%% Now plot different measured quantities as kymographs
if plot_kymographs
    % Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
    dvDir = fullfile(mKDir, 'avgDV') ;
    lDir = fullfile(mKDir, 'avgLeft') ;
    rDir = fullfile(mKDir, 'avgRight') ;
    dDir = fullfile(mKDir, 'avgDorsal') ;
    vDir = fullfile(mKDir, 'avgVentral') ;
    outdirs = {dvDir, lDir, rDir, dDir, vDir} ;
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    for qq = 1:length(outdirs)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end

        % Unpack what to plot (averaged kymographs, vary averaging region)
        HHK = HHsK{qq} ;
        gdotK = gdotsK{qq} ;
        divvK = divvsK{qq} ;
        velnK = velnsK{qq} ;
        H2vnK = H2vnsK{qq} ;
        m2plot = {gdotK, HHK, divvK, velnK, H2vnK} ;
        titles = {'$\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$',...
            'mean curvature, $H$', ...
            'divergence of flow, $\nabla \cdot \mathbf{v}$', ...
            'normal velocity, $v_n$', ...
            'normal motion, $v_n 2 H$'} ;
        labels = {['$\textrm{Tr}[g^{-1}\dot{g}]$ ' unitstr], ...
            ['mean curvature, $H$ ' Hunitstr], ...
            ['$\nabla \cdot \mathbf{v}$ ' unitstr], ...
            ['normal velocity, $v_n$ ' vunitstr] , ...
            ['normal motion, $v_n 2 H $ ' unitstr]} ;
        names = {'gdot', 'HH', 'divv', 'veln', 'H2vn'} ;
        climits = [climit, climit_H, climit, climit_veln, climit_err] ;

        %% Plot gdot/HH/divv/veln/H2vn DV-averaged kymograph
        for pp = 1:length(m2plot)
            close all
            set(gcf, 'visible', 'off')
            colormap bwr
            imagesc((1:nU)/nU, tps, m2plot{pp})
            caxis([-climits(pp), climits(pp)])
            % Add folds to plot
            hold on;
            fons1 = fons(1) ;
            fons2 = fons(2) ;
            fons3 = fons(3) ;
            plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
            plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
            plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

            % title and save
            title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
            ylabel(['time [' QS.timeunits ']'], 'Interpreter', 'Latex')
            xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
            cb = colorbar() ;
            ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
            fn = fullfile(odir, [ names{pp} '.png']) ;
            disp(['saving ', fn])
            export_fig(fn, '-png', '-nocrop', '-r200')   

            % Zoom in on small values
            caxis([-climits(pp)/3, climits(pp)/3])
            fn = fullfile(odir, [names{pp} '_zoom.png']) ;
            disp(['saving ', fn])
            export_fig(fn, '-png', '-nocrop', '-r200')   
            % Zoom in on early times
            ylim([min(tps), max(fons) + 10])
            caxis([-climits(pp)/3, climits(pp)/3])
            fn = fullfile(odir, [names{pp} '_zoom_early.png']) ;
            disp(['saving ', fn])
            export_fig(fn, '-png', '-nocrop', '-r200')   
        end
    end
end

%% Metric Kinematic Correlations
% Plot both all time and select times
timeSpans = {tps, tps(tps < max(fons) + 11)} ;
if plot_correlations
    outputFileNames = cell(2, 1) ;
    outputFileNames{1} = {fullfile(mKDir, 'correlation_alltime_div_gdot'), ...
                       fullfile(mKDir, 'correlation_earlytimes_div_gdot')} ;
    outputFileNames{2} = {fullfile(mKDir, 'correlation_alltime_2Hvn_gdot'), ...
                       fullfile(mKDir, 'correlation_earlytimes_2Hvn_gdot')} ;
    alphaVal = 0.6 ;
    sz = 10 ;
    cmap = parula ;
    close all
    set(gcf, 'visible', 'off')
    for pairIdx = 1:2
        for tspanIdx = 1:2
            fnout = outputFileNames{pairIdx}{tspanIdx} ;
            timeSpan_i = timeSpans{tspanIdx} ;
            ntspan = length(timeSpan_i) ;
            titles = {'left lateral', 'right lateral', 'dorsal', 'ventral'} ;
            markers = QS.plotting.markers ;
            colors = mapValueToColor(1:ntspan, [1, ntspan], cmap) ;
            close all
            cols = round(nV * [0.2, 0.85]) ;
            sphCollection = cell(4, 1) ;
            sposCollection = cell(4, 1) ;
            for qq = 1:4  % consider left, right, dorsal, ventral
                disp(['qq = ', num2str(qq), ': ', titles{qq}])
                if pairIdx == 1
                    divv = divvsK{qq + 1} ;
                    gdot = gdotsK{qq + 1} ;
                    sphCollection{qq} = subplot(2, 2, qq) ;
                    for row = 1:ntspan
                        disp(['row = ', num2str(row)])
                        scatter(divv(row, cols), gdot(row, cols), sz, ...
                            markers{qq}, 'MarkerFaceColor', 'none', ...
                            'MarkerEdgeColor', colors(row, :), ...
                            'MarkerEdgeAlpha', alphaVal) ;
                        hold on ;
                    end

                    % Label the x axis if on the bottom row
                    if qq > 2
                        xlabel(['$\nabla \cdot \bf{v}_\parallel$ ' unitstr], ...
                                'Interpreter', 'Latex') ;
                    end
                    axis equal
                    % Add dashed y=x line
                    xlims = get(gca, 'xlim') ;
                    ylims = get(gca, 'ylim') ;
                    leftdot = max(xlims(1), ylims(1)) ;
                    rightdot = min(xlims(2), ylims(2)) ;
                    plot([leftdot, rightdot], [leftdot, rightdot], 'k--')
                else
                    H2vn = H2vnsK{qq + 1} ;
                    gdot = gdotsK{qq + 1} ;
                    sphCollection{qq} = subplot(2, 2, qq) ;
                    for row = 1:ntspan
                        disp(['row = ', num2str(row)])
                        scatter(H2vn(row, cols), gdot(row, cols), sz, ...
                            markers{qq}, 'MarkerFaceColor', 'none', ...
                            'MarkerEdgeColor', colors(row, :), ...
                            'MarkerEdgeAlpha', alphaVal) ;
                        hold on ;
                    end

                    % Label the x axis if on the bottom row
                    if qq > 2
                        xlabel(['$2Hv_n$ ' unitstr], 'Interpreter', 'Latex') ;
                    end
                    axis equal
                    % Add dashed y=x line
                    xlims = get(gca, 'xlim') ;
                    ylims = get(gca, 'ylim') ;
                    leftdot = max(xlims(1), -ylims(2)) ;
                    rightdot = min(xlims(2), -ylims(1)) ;
                    plot([leftdot, rightdot], [-leftdot, -rightdot], 'k--')
                end

                % Label the y axis if on the left column
                if qq == 1 || qq == 3
                    ylabel(['$\textrm{Tr}[g^{-1} \dot{g}]$ ' unitstr], ...
                            'Interpreter', 'Latex')
                end
                title(titles{qq}, 'Interpreter', 'Latex')

                % Grab axis position
                sposCollection{qq} = get(sphCollection{qq}, 'Position');
            end

            % Move subplots left a bit for colorbar space
            for qq = 1:length(sphCollection)
                spos = sposCollection{qq} ;
                wh = min(spos(3)-0.05, spos(4)) ;
                if mod(qq, 2) == 1
                    set(sphCollection{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                else
                    set(sphCollection{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                end
            end

            % Add colorbar
            c = colorbar('Position',[.9 .333 .02 .333]) ;
            % Make colorbar share the alpha of the image
            % Manually flush the event queue and force MATLAB to render the colorbar
            % necessary on some versions
            drawnow
            % Get the color data of the object that correponds to the colorbar
            cdata = c.Face.Texture.CData;
            % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
            cdata(end,:) = uint8(alphaVal * cdata(end,:));
            % Ensure that the display respects the alpha channel
            c.Face.Texture.ColorType = 'truecoloralpha';
            % Update the color data with the new transparency information
            c.Face.Texture.CData = cdata;
            c.Label.Interpreter = 'Latex' ;
            c.Label.String = ['time [' QS.timeunits ']'] ;
            c.Ticks = [0, 1] ;
            c.TickLabels = [tps(1), max(timeSpan_i)] ;

            % Save figure
            saveas(gcf, [fnout '.png']) ;
            saveas(gcf, [fnout '.pdf']) ;
            close all
            set(gcf, 'visible', 'off')
        end
    end
    disp('done')
end


%% Metric Kinematics -- decompose into isotropic and other component
if plot_gdot_decomp
    % Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
    dvDir = fullfile(mKDir, 'avgDV') ;
    lDir = fullfile(mKDir, 'avgLeft') ;
    rDir = fullfile(mKDir, 'avgRight') ;
    dDir = fullfile(mKDir, 'avgDorsal') ;
    vDir = fullfile(mKDir, 'avgVentral') ;
    outdirs = {dvDir, lDir, rDir, dDir, vDir} ;
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    gdotK0 = gdotsK{1} ;
    isogrowth = sum(gdot0, 2) ;
    
    
    for qq = 1:length(outdirs)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end
        % Unpack what to plot (averaged kymographs, vary averaging region)
        HHK = HHsK{qq} ;
        gdotK = gdotsK{qq} ;
        divvK = divvsK{qq} ;
        velnK = velnsK{qq} ;
        H2vnK = H2vnsK{qq} ;
        m2plot = {gdotK, HHK, divvK, velnK, H2vnK} ;
        titles = {'$\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$',...
            'mean curvature, $H$', ...
            'divergence of flow, $\nabla \cdot \mathbf{v}$', ...
            'normal velocity, $v_n$', ...
            'normal motion, $v_n 2 H$'} ;
        labels = {['$\textrm{Tr}[g^{-1}\dot{g}]$ ' unitstr], ...
            ['mean curvature, $H$ ' Hunitstr], ...
            ['$\nabla \cdot \mathbf{v}$ ' unitstr], ...
            ['normal velocity, $v_n$ ' vunitstr] , ...
            ['normal motion, $v_n 2 H $ ' unitstr]} ;
        names = {'gdot', 'HH', 'divv', 'veln', 'H2vn'} ;
        climits = [climit, climit_H, climit, climit_veln, climit_err] ;

        %% Plot gdot/HH/divv/veln/H2vn DV-averaged kymograph
        for pp = 1:length(m2plot)
            close all
            set(gcf, 'visible', 'off')
            colormap bwr
            imagesc((1:nU)/nU, tps, m2plot{pp})
            caxis([-climits(pp), climits(pp)])
            % Add folds to plot
            hold on;
            fons1 = fons(1) ;
            fons2 = fons(2) ;
            fons3 = fons(3) ;
            plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
            plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
            plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

            % title and save
            title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
            ylabel(['time [' QS.timeunits ']'], 'Interpreter', 'Latex')
            xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
            cb = colorbar() ;
            ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
            fn = fullfile(odir, [ names{pp} '.png']) ;
            disp(['saving ', fn])
            export_fig(fn, '-png', '-nocrop', '-r200')   
        end
    end
end


