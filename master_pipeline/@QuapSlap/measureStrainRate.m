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
clim_trgdot = 0.2 ;
clim_tg = 1 ;

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

%% Prepare both metric styles
dirs2make = { fullfile(egImDir) } ;
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
    clearvars tmp

    % Define metric strain filename        
    estrainFn = sprintf(QS.fullFileBase.strainRate, tp) ;
    
    % Compute the metric strain if not on disk
    if ~exist(estrainFn, 'file') || overwrite
        if exist(estrainFn, 'file')
            disp('Overwriting metric strain on disk')
        else
            disp('Computing metric strain anew')
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

        % Smooth the velocities in space using gptoolbox
        vraw = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, :)) ;
        vs = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
            lambda, 'implicit', vraw) ;
        
        % Push vectors onto faces
        [V2F, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
        vf = V2F * vs ;
        
        % Convert to 2D mesh
        mesh.nU = QS.nU ;
        cutMesh = cutRectilinearCylMesh(mesh) ;
        
        % Compute the strain rate tensor
        [~, dvi] = vectorCovariantDerivative(vf, cutMesh.f, cutMesh.v, cutMesh.u) ;
        dvij = cell(size(dvi)) ;
        for qq = 1:length(dvi)
            dvij{qq} = 0.5 * ( dvi{qq} + dvi{qq}' ) ;
        end
        
        % Compute the second fundamental form
        [gg, bb] = constructFundamentalForms(cutMesh.f, cutMesh.v, cutMesh.u) ;
        
        % Decompose velocity into components
        [v0n, ~] = resolveTangentNormalVector(cutMesh.f, cutMesh.v, vf) ;
        
        % Strain rate tensor
        epsilon = cell(size(dvi)) ;
        for qq = 1:length(dvi)
            epsilon{qq} = dvij{qq} - v0n(qq) .* bb{qq} ;
        end
        
        %% Construct Topolgical Structure Tools ===============================
        % eg = metricStrainSPhiGridMesh(mesh, mesh2) ;

        % Cut the mesh into a cutMesh to open up the pullback mesh & also
        % grab du and dv for each face --  bond vecs along u and along v
        mesh.nU = nU ;
        tmp_options.preview = false ;
        [~, dbonds] = labelRectilinearMeshBonds(mesh, tmp_options) ;
        cutMesh = dbonds.cutMesh ; 
        
        % Metric strain -- separate trace and deviatoric strain comp, angle
        tre = zeros(size(epsilon, 1), 1) ;  % traceful dilation
        dve = zeros(size(epsilon, 1), 1) ;  % deviatoric magnitude
        the = zeros(size(epsilon, 1), 1) ;  % angle of elongation
        for qq = 1:size(g0cell, 1)
            eq = epsilon{qq} ;
            gq = gg{qq} ;
            % traceful component
            tre(qq) = trace(inv(gq) * (eq)) ;
            % deviatoric component
            dve(qq) = sqrt(trace(inv(gq) * (eq * (inv(gq) * eq)))) ;
            % angle of elongation -- first take eigvectors
            [evec, eval] = eig() ;
            [evals, idx] = sort(diag(eval)) ;
            evec = evec(:, idx) ;
            pevec = evec(:, end) ;
            assert(prod(evals) == dve(qq)) ;
            the(qq) = atan2(pevec(2), pevec(1)) ;
        end

        % save the metric strain
        readme.eg = 'strain rate epsilon=1/2(nabla_i v_j + nabla_j v_i) - vn b_ij' ;
        readme.tre = 'Tr[g^{-1} epsilon]';
        readme.dve = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] )';
        readme.the = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector';
        readme.lambda = 'Laplacian smoothing on velocities' ;
        readme.lambda_mesh = 'Laplacian smoothing on mesh vertices' ;
        disp(['saving ', estrainFn])
        save(estrainFn, 'eg', 'tre', 'dve', 'the', ...
            'lambda', 'lambda_mesh', 'readme')
    else
        % load the metric strain
        load(estrainFn, 'eg', 'tre', 'dve', 'the')
    end        

    %% Plot the metric components on trisurf
    % denom = sqrt(tg(:, 1, 1) .* tg(:, 2, 2)) ;
    labels = {'$\mathrm{Tr} [g^{-1} \varepsilon] $', ...
        '$\sqrt\left( \mathrm{Tr}[g^{-1} \varepsilon g^{-1} \varepsilon] \right)$'} ;
    time_in_units = (tp - tfold) * QS.timeInterval ;
    tstr = [': $t=$', sprintf('%03d', time_in_units), QS.timeUnits ];
    es = [tre, dve, the] ;
    
    %% consider each metric element & plot in 3d
    fn = fullfile(egImDir, metric_style, ...
            sprintf([QS.fileBase.spcutMeshSmRSC '.png'], tp));
    if ~exist(fn, 'file') || overwrite
        clf
        set(gcf, 'visible', 'off') ;
        for qq = 1:2
            % For each view (dorsal, ventral, left, right)
            % for pp = 1:4
            subplot(1, 2, qq)
            colors = es(:, gelem(qq)) ;
            
            trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
                colors, 'edgecolor', 'none')
            axis equal
            cb = colorbar() ;
            caxis([-clim, clim])
            title(labels{qq}, 'Interpreter', 'Latex')            
            ylabel(cb, labels{qq}, 'Interpreter', 'Latex')
            
            colormap(bwr256)
            xlim(xyzlim(1, :))
            ylim(xyzlim(2, :))
            zlim(xyzlim(3, :))
            axis off
            
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
    fn = fullfile(egImDir, ['metricstrain_' metric_style '_2d'], ...
            sprintf([QS.fileBase.spcutMeshSm '.png'], tp));
    if ~exist(fn, 'file') || overwrite
        for qq = 1:4
            subplot(2, 2, qq) 
            if qq < 4
                disp(['coloring by tg ' num2str(gelem(qq))])
                colors = tg(:, gelem(qq)) ;
            else
                disp('coloring by trace')
                colors = dilation ;
            end
            trisurf(cutMesh.f, ...
                cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
                cutMesh.u(:, 2), 0 * cutMesh.u(:, 2), ...
                colors, 'edgecolor', 'none')
            daspect([1,1,1])
            cb = colorbar() ;

            if strcmp(metric_style, 'mesh')                 
                caxis([-0.5, 0.5])
                title(['surface deformation rate, ', labels{qq}, tstr], ...
                    'Interpreter', 'Latex')            
                ylabel(cb, labels{qq}, 'Interpreter', 'Latex')
            elseif strcmp(metric_style, 'strain')  
                if qq < 4
                    caxis([-clim_tg, clim_tg])
                else
                    caxis([-clim_trgdot, clim_trgdot])
                end
                title(['strain rate, ', strainlabels{qq}, tstr], ...
                    'Interpreter', 'Latex')            
                ylabel(cb, strainlabels{qq}, 'Interpreter', 'Latex')
            end
            % xlabel('AP position, [$\mu$m]', 'Interpreter', 'Latex')
            % ylabel('lateral position, [$\mu$m]', 'Interpreter', 'Latex')
            % zlabel('DV position, [$\mu$m]', 'Interpreter', 'Latex')
            colormap(bwr256)
            axis off
            view(2)
        end
        % Save the image
        saveas(gcf, fn) ;
        clf
    end    
    close all
end


