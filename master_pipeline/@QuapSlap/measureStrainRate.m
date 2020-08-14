function measureStrainRate(QS, options)
%measureStrainRate(QS, options)
%   Compute epsilon = 1/2 (\nabla_i v_j + \nabla_j v_i) - vN b_{ij} 
%   
% Parameters
% ----------
%
% Returns 
% -------
%
% NPMitchell 2020

%% Default options
lambda = 0.01 ;
lambda_mesh = 0.002 ;
overwrite = false ;
preview = true ;
clim_trace = 1 ;
clim_deviatoric = 0.5 ;
averagingStyle = 'Lagrangian' ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'averagingStyle')
    averagingStyle = options.averagingStyle ;
end


%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
    sresStr = '' ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
    sresStr = 'doubleRes_' ;
else 
    error("Could not parse samplingResolution: set to '1x' or '2x'")
end

%% Unpack QS
egDir = QS.dir.strainRate ;
egImDir = fullfile(QS.dir.strainRate, 'images') ;
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;

%% Prepare for plots
colors = define_colors ;
blue = colors(1, :) ;
red = colors(2, :) ;
green = colors(3, :) ;

%% Colormap
close all
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;
clf
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([0, 1])
pos256 = bluewhitered(256) ;
close all
pm256 = phasemap(256) ;

%% Prepare both metric styles
dirs2make = { egImDir, fullfile(egImDir, 'strainRate3d'), ...
    fullfile(egImDir, 'strainRate2d') } ;
for ii = 1:length(dirs2make)
    dir2make = dirs2make{ii} ;
    if ~exist(dir2make, 'dir')
        mkdir(dir2make)
    end
end

%% ------------------------------------------------------------------------
% Construct 2D mesh corresponding to the planar domain of
% parameterization
%--------------------------------------------------------------------------
% % Load from file
% path = fullfile(NESpath, 'NES_Examples') ;
% mesh = read_ply_mod(fullfile(path, 'tube_simple_h1p00_R0p00_w1p00.ply')) ;
% rmID = [length(mesh.v)-1, length(mesh.v)] ;
% [F, V] = remove_vertex_from_mesh(mesh.f, mesh.v, rmID) ;
% % Center the mesh around the x axis
% midx = 0.5 * (max(mesh.v(:, 1)) + min(mesh.v(:, 1))) ;
% V(:, 1) = V(:, 1) - midx ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Consider each kind of metric strain measurement
nU = QS.nU ;
nV = QS.nV ;
QS.t0set() ;
tfold = QS.t0 ;

% Load vertex-based velocity measurements
QS.getVelocityAverage('vv', 'vf')
vertex_vels = QS.velocityAverage.vv ;
face_vels = QS.velocityAverage.vf ;

% Pre-assign timepoints with velocities
tpts = QS.xp.fileMeta.timePoints(1:end-1) ;
tp2do = [60, tpts(1:10:end), setdiff(tpts, tpts(1:10:end))] ;

% Build metric from mesh
for tp = tp2do
    disp(['t = ' num2str(tp)])
    tidx = QS.xp.tIdx(tp) ;

    % Load current mesh
    tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp)) ;
    mesh = tmp.spcutMeshSmRSC ;
    % Normalize the zeta axis to unity or ar (aspectratio relaxed)
    mesh.u(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) * mesh.ar ;
    clearvars tmp

    % Define metric strain filename        
    estrainFn = sprintf(QS.fullFileBase.strainRate, tp) ;
    
    % Compute the metric strain if not on disk
    if ~exist(estrainFn, 'file') || overwrite
        if exist(estrainFn, 'file')
            disp('Overwriting strain rate on disk')
        else
            disp('Computing strain rate anew')
        end
        
        % Smooth the mesh vertices
        if lambda_mesh > 0 
            tri = triangulation(mesh.f, mesh.v) ;
            fbndy = tri.freeBoundary ;
            fbndy = fbndy(:, 1) ;
            mesh.v = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
                lambda_mesh, 'implicit', mesh.v) ;
            % check smoothed mesh
            % trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none')
        end
        
        % Obtain mean curvature H for checking against trace(b_ij)
        DEC = DiscreteExteriorCalculus(mesh.f, mesh.v) ;
        H3d = sum(mesh.vn .* DEC.laplacian(mesh.v), 2) * 0.5 ;
        H2d = H3d ;
        H2d(nU*(nV-1)+1:(nU*nV)) = H3d(1:nU) ;
        
        % Smooth the velocities in space using gptoolbox
        vraw = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, :)) ;
        vs = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
            lambda, 'implicit', vraw) ;
        
%         
%         % Smooth the velocities in space using gptoolbox
%         vx = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 1)) ;
%         vy = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 2)) ;
%         vz = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 3)) ;
%         vxs = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
%             lambda, 'implicit', vx') ;
%         vys = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
%             lambda, 'implicit', vy') ;
%         vzs = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
%             lambda, 'implicit', vz') ;
%         vs2 = [vxs, vys, vzs] ;
%         
%         assert(~any(any(vs2-vs))) 
        
        % Push vectors onto faces
        [V2F, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
        vf = V2F * vs ;
        
        %% Checking -- debug
        Hf = V2F * H3d ;
        divv = DEC.divergence(vf) ;
        divf = V2F * divv ;
        
        %% Convert to 2D mesh
        mesh.nU = QS.nU ;
        cutMesh = cutRectilinearCylMesh(mesh) ;
        
        % Compute the strain rate tensor
        disp('Computing covariantDerivative')
        [~, dvi] = vectorCovariantDerivative(vf, cutMesh.f, cutMesh.v, cutMesh.u) ;
        dvij = cell(size(dvi)) ;
        for qq = 1:length(dvi)
            dvij{qq} = 0.5 * ( dvi{qq} + dvi{qq}' ) ;
        end
        
        % Compute the second fundamental form
        [gg, ~] = constructFundamentalForms(cutMesh.f, cutMesh.v, cutMesh.u) ;
        [~, bb] = constructFundamentalForms(mesh.f, mesh.v, mesh.u) ;
        
        % Decompose velocity into components
        [v0n, ~] = resolveTangentNormalVector(cutMesh.f, cutMesh.v, vf) ;
        
        % Strain rate tensor
        strainrate = cell(size(dvi)) ;
        tre = zeros(size(dvi)) ;
        checkH = zeros(size(dvi)) ;
        checkdiv = zeros(size(dvi)) ;
        for qq = 1:length(dvi)
            strainrate{qq} = dvij{qq} - v0n(qq) .* bb{qq} ;
            tre(qq) = trace(inv(gg{qq}) * dvij{qq}) - v0n(qq) .* Hf(qq) ;
            checkH(qq) = trace(inv(gg{qq}) * bb{qq}) ;        % == 2 * Hf(qq) ; 
            checkdiv(qq) = trace(inv(gg{qq}) * dvij{qq}) ;    % == divf(qq) ; 
        end
        
        clf;
        figure;
        plot(checkdiv, divf, '.')
        hold on;
        plot(checkH, 2* Hf, 'o') 
        xlabel('Tr$b_{ij}$, Tr$\nabla_i v_j$', 'Interpreter', 'Latex')
        ylabel('$2H$, Tr$\nabla \cdot \mathbf{v}$', 'Interpreter', 'Latex')
        
        if preview
            %% Check Mean curvature
            subplot(2, 2, 1)
            trisurf(triangulation(mesh.f, mesh.v), checkH, 'edgecolor', 'none')
            axis equal; caxis([-.1, .1]); colorbar() ;
            title('Tr$\left[g^{-1} b\right]$', 'interpreter', 'latex')
            subplot(2, 2, 2)
            trisurf(triangulation(mesh.f, mesh.v), 2 * Hf, 'edgecolor', 'none')
            axis equal; caxis([-.1, .1]); colorbar() ;
            title('$2H$', 'interpreter', 'latex')
            subplot(2, 1, 2)
            trisurf(triangulation(mesh.f, mesh.v), checkH - 2 * Hf, 'edgecolor', 'none')
            axis equal; caxis([-.1, .1]); colorbar() ;
            title('Tr$[g^{-1}b] - 2H$', 'interpreter', 'latex')
            colormap(bwr)
        end

        %% Check div(v)
        clf
        clim = max(mean(abs(checkdiv)), mean(abs(divf))) ;
        subplot(2, 2, 1)
        trisurf(triangulation(mesh.f, mesh.v), checkdiv, ...
            'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        checkDivStr = '$\frac{1}{2}$Tr$\left[g^{-1} \left(\nabla_i v_j + \nabla_j v_i \right)\right]$' ;
        title(checkDivStr, 'interpreter', 'latex')
        subplot(2, 2, 2)
        trisurf(triangulation(mesh.f, mesh.v), divf, 'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        title('$\nabla \cdot \mathbf{v}_{\parallel}$', 'interpreter', 'latex')
        subplot(2, 1, 2)
        trisurf(triangulation(mesh.f, mesh.v), checkdiv - 2*pi*divf, ...
            'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        title([checkDivStr ' $-2\pi \nabla \cdot \mathbf{v}_{\parallel}$'], ...
            'interpreter', 'latex')
        colormap(bwr)
        
        %% Construct Topolgical Structure Tools ===============================
        % eg = metricStrainSPhiGridMesh(mesh, mesh2) ;
        % 
        % % Cut the mesh into a cutMesh to open up the pullback mesh & also
        % % grab du and dv for each face --  bond vecs along u and along v
        % mesh.nU = nU ;
        % tmp_options.preview = false ;
        % [~, dbonds] = labelRectilinearMeshBonds(mesh, tmp_options) ;
        % cutMesh = dbonds.cutMesh ; 
        
        % Metric strain -- separate trace and deviatoric strain comp, angle
        treps = zeros(size(strainrate, 1), 1) ;  % traceful dilation
        dvtre = zeros(size(strainrate, 1), 1) ;  % deviatoric magnitude
        theta = zeros(size(strainrate, 1), 1) ;  % angle of elongation
        for qq = 1:size(strainrate, 1)
            eq = strainrate{qq} ;
            gq = gg{qq} ;
            % traceful component -- 1/2 Tr[g^{-1} gdot] = Tr[g^{-1} eps] 
            treps(qq) = trace(inv(gq) * (eq)) ;
            % deviatoric component -- 
            % || epsilon - 1/2 Tr[g^{-1} epsilon] g|| = sqrt(Tr[A A^T]),
            % where A = epsilon - 1/2 Tr[g^{-1} epsilon] g.
            AA = eq - 0.5 * treps(qq) * gq ;
            dvtre(qq) = sqrt(trace(inv(gq) * (AA * (inv(gq) * AA)))) ;
            % angle of elongation -- first take eigvectors
            [evec, eval] = eig(AA) ;
            [evals, idx] = sort(diag(eval)) ;
            evec = evec(:, idx) ;
            pevec = evec(:, end) ;
            theta(qq) = atan2(pevec(2), pevec(1)) ;
        end
        theta = mod(theta, 2 * pi) ;
        
        % save the metric strain
        readme.eg = 'strain rate epsilon=1/2(nabla_i v_j + nabla_j v_i) - vn b_ij' ;
        readme.treps = 'Tr[g^{-1} epsilon]';
        readme.dvtre = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] )';
        readme.theta = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector';
        readme.lambda = 'Laplacian smoothing on velocities' ;
        readme.lambda_mesh = 'Laplacian smoothing on mesh vertices' ;
        readme.note = 'The pullback space is taken to range from zeta=[0, 1] and phi=[0, 1]' ; 
        disp(['saving ', estrainFn])
        save(estrainFn, 'strainrate', 'treps', 'dvtre', 'theta', ...
            'lambda', 'lambda_mesh', 'readme')
    else
        % load the metric strain
        load(estrainFn, 'strainrate', 'treps', 'dvtre', 'theta')
        % Convert to 2D mesh
        mesh.nU = QS.nU ;
        cutMesh = cutRectilinearCylMesh(mesh) ;
    end        

    %% Plot the metric components on trisurf
    % denom = sqrt(tg(:, 1, 1) .* tg(:, 2, 2)) ;
    % NOTE: \varepsilon --> ${\boldmath${\varepsilon}$}$
    labels = {'$\mathrm{Tr} [\bf{g}^{-1}\varepsilon] $', ...
        '$||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\varepsilon\right]\bf{g}||$'} ;
    time_in_units = (tp - tfold) * QS.timeInterval ;
    tstr = [': $t=$', sprintf('%03d', time_in_units), QS.timeUnits ];
    
    %% consider each metric element & plot in 3d
    fn = fullfile(egImDir, 'strainRate3d', sprintf([QS.fileBase.spcutMeshSmRSC '.png'], tp));
    if ~exist(fn, 'file') || overwrite
        clf
        set(gcf, 'visible', 'off') ;
        for qq = 1:2
            % For each view (dorsal, ventral, left, right)
            % for pp = 1:4
            subplot(1, 2, qq)
            if qq == 1
                trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
                    'FaceVertexCData', treps, 'edgecolor', 'none')
                caxis([-clim_trace, clim_trace])
                colormap(gca, bwr256)
                cb = colorbar('location', 'southOutside') ;      
                % ylabel(cb, labels{qq}, 'Interpreter', 'Latex')

            else
                % Intensity from dvtre and color from the theta
                indx = max(1, round(mod(2*theta, 2*pi)*size(pm256, 1)/(2 * pi))) ;
                colors = pm256(indx, :) ;
                colors = min(dvtre / clim_deviatoric, 1) .* colors ;
                trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
                    'FaceVertexCData', colors, 'edgecolor', 'none')

                % Colorbar and phasewheel
                colormap(gca, phasemap)
                phasebar('colormap', phasemap, ...
                    'location', [0.82, 0.1, 0.1, 0.135], 'style', 'nematic')
                ax = gca ;
                get(gca, 'position')
                cb = colorbar('location', 'southOutside') ;
                drawnow
                axpos = get(ax, 'position') ;
                cbpos = get(cb, 'position') ;
                set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3)*0.6, cbpos(4)])
                set(ax, 'position', axpos) 
                hold on;
                caxis([0, clim_deviatoric])
                colormap(gca, gray)
            end
            
            axis equal
            xlim(xyzlim(1, :))
            ylim(xyzlim(2, :))
            zlim(xyzlim(3, :))
            axis off
            title(labels{qq}, 'Interpreter', 'Latex')      
            
            % Save images;
            % dorsal
            % if pp == 1
            %     view(0, 90)
            % ventral
            % elseif pp == 2
            %     view(0, -90)
            % left
            % elseif pp == 3
            %     view(0, 0)
            % right
            % elseif pp == 4
            %     view(0, 180)
            % end
            % end
            
            % left view
            view(0, 0) 
        end
        sgtitle(['strain rate, ', tstr], 'Interpreter', 'latex') 

        % Save the image
        saveas(gcf, fn) ;
        clf

    end
    
    %% Now plot in 2d
    close all
    set(gcf, 'visible', 'off') ;
    fn = fullfile(egImDir, 'strainRate2d', ...
            sprintf([QS.fileBase.spcutMeshSm '.png'], tp));
    if ~exist(fn, 'file') || overwrite
        % Panel 1
        subplot(1, 2, 1)
        trisurf(cutMesh.f, ...
            cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
            cutMesh.u(:, 2), 0 * cutMesh.u(:, 2), ...
            'FaceVertexCData', treps, 'edgecolor', 'none')
        daspect([1,1,1])
        cb = colorbar('location', 'southOutside') ;

        caxis([-clim_trace, clim_trace])
        title(labels{1}, 'Interpreter', 'Latex')   
        colormap(bwr256)
        axis off
        view(2)
            
        % Panel 2 
        subplot(1, 2, 2)
        % Intensity from dvtre and color from the theta
        indx = max(1, round(mod(2*theta, 2*pi)*size(pm256, 1)/(2 * pi))) ;
        colors = pm256(indx, :) ;
        colors = min(dvtre / clim_deviatoric, 1) .* colors ;
        trisurf(cutMesh.f, cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
            cutMesh.u(:, 2), 0*cutMesh.u(:, 1), ...
            'FaceVertexCData', colors, 'edgecolor', 'none')
        daspect([1,1,1])
        title(labels{2}, 'Interpreter', 'Latex')   
        
        % Colorbar and phasewheel
        colormap(gca, phasemap)
        phasebar('colormap', phasemap, ...
            'location', [0.82, 0.12, 0.1, 0.135], 'style', 'nematic')
        axis off
        view(2)
        ax = gca ;
        get(gca, 'position')
        cb = colorbar('location', 'southOutside') ;
        drawnow
        axpos = get(ax, 'position') ;
        cbpos = get(cb, 'position') ;
        set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3)*0.6, cbpos(4)])
        set(ax, 'position', axpos) 
        hold on;
        caxis([0, clim_deviatoric])
        colormap(gca, gray)
        
        % Save the image
        sgtitle(['strain rate, ', tstr], 'Interpreter', 'latex') 
        saveas(gcf, fn) ;
        clf
    end    
    close all
    
    
    %% Compare trace to trace of gdot determined via kinematics
    close all
    set(gcf, 'visible', 'off') ;
    fn = fullfile(egImDir, 'strainRate2d', ...
            sprintf(['compare_' QS.fileBase.spcutMeshSm '.png'], tp));
    if ~exist(fn, 'file') || overwrite
        % Load gdot trace from kinematics
        fn_gdot = sprintf(QS.fullFileBase.metricKinematics.gdot, lambda, ...
            lambda, lambda_mesh, tidx) ;
        load([strrep(fn_gdot, '.', 'p'), '.mat'], 'gdot')
        
        % Panel 1
        subplot(1, 2, 1)
        trisurf(cutMesh.f, ...
            cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
            cutMesh.u(:, 2), 0 * cutMesh.u(:, 2), ...
            'FaceVertexCData', treps, 'edgecolor', 'none')
        daspect([1,1,1])
        cb = colorbar('location', 'southOutside') ;
        caxis([-clim_trace, clim_trace])
        title(labels{1}, 'Interpreter', 'Latex')   
        colormap(bwr256)
        axis off
        view(2)
            
        % Panel 2 
        subplot(1, 2, 2)
        % Comparison 1/2 * Tr[g^{-1}gdot]
        trisurf(cutMesh.f, cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
            cutMesh.u(:, 2), 0*cutMesh.u(:, 1), ...
            gdot, 'edgecolor', 'none')
        daspect([1,1,1])
        cb = colorbar('location', 'southOutside') ;
        caxis([-clim_trace * 0.5, clim_trace * 0.5])
        title('$\frac{1}{2}$Tr$[g^{-1}\dot{g}]$', 'Interpreter', 'Latex')   
        colormap(bwr256)
        axis off
        view(2)
                    
        % Save the image
        sgtitle(['comparison $\frac{1}{2}$Tr$[g^{-1}\dot{g}]$, ', tstr], 'Interpreter', 'latex') 
        saveas(gcf, fn) ;
        clf
    end    
    close all
    
    
    %% Compare to freshly computed
    close all
    set(gcf, 'visible', 'off') ;
    fn = fullfile(egImDir, 'strainRate2d', ...
            sprintf(['compare_fresh_' QS.fileBase.spcutMeshSm '.png'], tp));
    if ~exist(fn, 'file') || overwrite 
        DEC = DiscreteExteriorCalculus(mesh.f, mesh.v) ;
        H3d = sum(mesh.vn .* DEC.laplacian(mesh.v), 2) * 0.5 ;
        
        %% OPTION1: Could compute divergence anew on smoothed mesh
        % fvel = squeeze(face_vels(tidx, :, :)) ;
        % divv3d_test = DEC.divergence(fvel) ;
        
        %% OPTION2: Load divv from disk and convert using smoothed mesh
        % [~, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
        if strcmp(averagingStyle, 'Lagrangian')
            if doubleResolution
                dec_tp = load(sprintf(QS.fullFileBase.decAvg2x, tp)) ;
            else
                dec_tp = load(sprintf(QS.fullFileBase.decAvg, tp)) ;
            end
        else
            if doubleResolution
                dec_tp = load(sprintf(QS.fullFileBase.decSimAvg2x, tp)) ;
            else
                dec_tp = load(sprintf(QS.fullFileBase.decSimAvg, tp)) ;
            end
        end
        
        % Smooth divergence(v) [divv]
        divv3d = laplacian_smooth(mesh.v, mesh.f, 'cotan', [],...
                                lambda, 'implicit', dec_tp.divs.raw) ;

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

        % Actual normal velocity -- currently using ORIGINAL mesh normals
        veln3d = sum(cat(2, vxs, vys, vzs) .* mesh.vn, 2) ;

        % Predict the divergence
        H2vn3d = 2 * H3d .* veln3d ;
        
        % Extend to have another row for 2d map
        divv2d = divv3d ;
        divv2d(nU*(nV-1) + 1:nU*nV) = divv3d(1:nU) ;
        H2d = H3d ;
        H2d(nU*(nV-1)+1:(nU*nV)) = H3d(1:nU) ;
        veln2d = veln3d ;
        veln2d(nU*(nV-1)+1:(nU*nV)) = veln3d(1:nU) ;
        H2vn2d = 2 * H2d .* veln2d ;

        % The difference
        gdot2d = divv2d - H2vn2d ;
        gdot3d = divv3d - H2vn3d ;
        
        
        % Panel 1
        subplot(1, 2, 1)
        trisurf(cutMesh.f, ...
            cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
            cutMesh.u(:, 2), 0 * cutMesh.u(:, 2), ...
            'FaceVertexCData', treps, 'edgecolor', 'none')
        daspect([1,1,1])
        cb = colorbar('location', 'southOutside') ;
        caxis([-clim_trace, clim_trace])
        title(labels{1}, 'Interpreter', 'Latex')   
        colormap(bwr256)
        axis off
        view(2)
            
        % Panel 2 
        subplot(1, 2, 2)
        % Comparison 1/2 * Tr[g^{-1}gdot]
        trisurf(cutMesh.f, cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
            cutMesh.u(:, 2), 0*cutMesh.u(:, 1), ...
            gdot2d(:), 'edgecolor', 'none')
        daspect([1,1,1])
        cb = colorbar('location', 'southOutside') ;
        caxis([-clim_trace * 0.5, clim_trace * 0.5])
        title('$\frac{1}{2}$Tr$[g^{-1}\dot{g}]$', 'Interpreter', 'Latex')   
        colormap(bwr256)
        axis off
        view(2)
                    
        % Save the image
        sgtitle(['comparison fresh $\frac{1}{2}$Tr$[g^{-1}\dot{g}]$, ', tstr], 'Interpreter', 'latex') 
        saveas(gcf, fn) ;
        clf
    end
    error('stopping')
end


