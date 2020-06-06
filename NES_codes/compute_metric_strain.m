%% Compute metric strain for evolving mesh ===============================
% Compute the metric strain for topological equivalent meshes over time
% (same connectivity, different geometry in embedding coordinates)

clear; close all; clc;

% Figure Options
figWidth = 20 ; 
figHeight = 12 ;
check = false ;
overwrite = false ;
lambda = 0.1 ; % smoothing parameter for velocity tensor strain

% Add paths
gutDir = '/mnt/data/code/gut_matlab/' ;
addpath(fullfile(gutDir, 'master_pipeline')) ;
addpath(fullfile(gutDir, 'addpath_recurse')) ;
addpath_recurse(fullfile(gutDir, 'basics')) ;
addpath_recurse(fullfile(gutDir, 'geometry')) ;
addpath_recurse(fullfile(gutDir, 'mesh_handling')) ;
addpath_recurse(fullfile(gutDir, 'plotting')) ;
addpath_recurse('/mnt/data/code/gptoolbox/') ;
addpath('/mnt/data/code/DEC/') ;

% Path options
NESpath = '/mnt/data/code/NonEuclideanShells/NES/' ;
addpath_recurse(NESpath)

dataDir = '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit_smallbox/' ;
meshDir = fullfile(dataDir, 'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/') ;

% Output dir defined by QS
stackResolution = [.2619 .2619 .2619] ;
nChannels = 1 ;
channelsUsed = 1 ;
timePoints = 0:169 ;
ssfactor = 4 ;

% Create xp
fn = 'Time_%06d_c1_stab';
xp = project.Experiment(dataDir, dataDir);
% A filename base template - to be used throughout this script
fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fn, '.tif'];
fileMeta.nChannels          = nChannels;
fileMeta.timePoints         = 110:263 ;
fileMeta.stackResolution    = stackResolution;
fileMeta.swapZT             = 1;
% first_tp is also required, which sets the tp to do individually.
first_tp = 1 ;
expMeta                     = struct();
expMeta.channelsUsed        = channelsUsed;
expMeta.channelColor        = 1;
expMeta.description         = 'Apical membrane in Drosophila gut';
expMeta.dynamicSurface      = 1;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.fitTime             = fileMeta.timePoints(first_tp);
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';
% Now set the meta data in the experiment.
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();

%% Create QuapSlap
opts.meshDir = meshDir ;
opts.flipy = false ;
opts.timeinterval = 1 ;
opts.timeunits = 'min';
opts.nV = 100 ;
opts.nU = 100 ;
opts.normalShift = 10 ;
opts.a_fixed = 2.0 ;
opts.adjustlow = 1.0 ;
disp('defining QS')
QS = QuapSlap(xp, opts) ;
disp('done')

%% Unpack QS
egDir = QS.dir.gstrain ;
engEgDir = QS.dir.gstrainRate ;
engEgImDir = QS.dir.gstrainRateIm ;
velEgDir = QS.dir.gstrainVel ;
velEgImDir = QS.dir.gstrainVelIm ;
geoEgDir = QS.dir.gstrainMesh ;
geoEgImDir = QS.dir.gstrainMeshIm ;
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;

%% Prepare for plots
colors = define_colors ;
blue = colors(1, :) ;
red = colors(2, :) ;
green = colors(3, :) ;

%% Prepare both metric styles
metric_styles = {'mesh', 'velocity', 'strain'} ;
imdirs = {geoEgImDir, velEgImDir, engEgImDir};
for qq = 1:3
    imdir = imdirs{qq} ;
    dirs2make = {fullfile(imdir, ['gss_ventral_' metric_styles{qq}]), ...
        fullfile(imdir, ['gss_dorsal_' metric_styles{qq}]), ...
        fullfile(imdir, ['gss_lateral_left_' metric_styles{qq}]), ...
        fullfile(imdir, ['gss_lateral_right_' metric_styles{qq}]), ...
        fullfile(imdir, ['gsphi_ventral_' metric_styles{qq}]), ...
        fullfile(imdir, ['gsphi_dorsal_' metric_styles{qq}]), ...
        fullfile(imdir, ['gsphi_lateral_left_' metric_styles{qq}]), ...
        fullfile(imdir, ['gsphi_lateral_right_' metric_styles{qq}]), ...
        fullfile(imdir, ['gphiphi_ventral_' metric_styles{qq}]), ...
        fullfile(imdir, ['gphiphi_dorsal_' metric_styles{qq}]), ...
        fullfile(imdir, ['gphiphi_lateral_left_' metric_styles{qq}]), ...
        fullfile(imdir, ['gphiphi_lateral_right_' metric_styles{qq}]) } ;
    for ii = 1:length(dirs2make)
        dir2make = dirs2make{ii} ;
        if ~exist(dir2make, 'dir')
            mkdir(dir2make)
        end
    end
end

%% -------------------------------------------------------------------------
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

%% Incompressibility
[cumerr, HHseries, divvseries, velnseries] = QS.measureCompressibility(0.1, 0.3) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Consider each kind of metric strain measurement
nU = QS.nU ;
nV = QS.nV ;
QS.t0set() ;
tfold = QS.t0 ;

% Load vertex-based velocity measurements
vvsmMfn = fullfile(QS.dir.pivSimAvg, 'vvM_simpletimeavg.mat')  ;
tmp = load(vvsmMfn) ;
vertex_vels = tmp.vvsmM ;
vfsmMfn = fullfile(QS.dir.pivSimAvg, 'vfM_simpletimeavg.mat') ;
tmp = load(vfsmMfn) ;
face_vels = tmp.vfsmM ;

for metric_style_kk = [3]
    metric_style = metric_styles{metric_style_kk} ;
    
    % Build metric from mesh
    for tp = xp.fileMeta.timePoints(1:end-1)
        disp(['t = ' num2str(tp)])
        tidx = xp.tIdx(tp) ;
        
        % Load current mesh
        tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp)) ;
        mesh = tmp.spcutMeshSmRSC ;
        clearvars tmp

        % Store mesh vertices and faces
        V0 = mesh.v ;
        F0 = mesh.f ;
        U0 = mesh.u ;
        
        % Save as ply if not already there
        fnply = sprintf(QS.fullFileBase.spcutMeshSmRSCPLY, tp) ;
        if ~exist(fnply, 'file')
            meshu1 = mesh.u(:, 1) / max(mesh.u(:, 1)) ;
            urgb = cat(2, meshu1, mesh.u(:, 2), 0*mesh.u(:, 1)) ;
            plywrite_with_normals(fnply, mesh.f, mesh.v, mesh.vn, urgb)
        end    
        
        % define metric strain filename        
        if strcmp(metric_style, 'mesh')
            gstrainFn = sprintf(QS.fullFileBase.gstrainMesh, tp) ;
        elseif strcmp(metric_style, 'velocity')
            gstrainFn = sprintf(QS.fullFileBase.gstrainVel, tp) ;
        elseif strcmp(metric_style, 'strain')
            gstrainFn = sprintf(QS.fullFileBase.gstrainRate, tp) ;
        end
       
        if ~exist(gstrainFn, 'file')
            if strcmp(metric_style, 'mesh')
                % Load next mesh also
                tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp + 1)) ;
                mesh2 = tmp.spcutMeshSmRSC ;
                clearvars tmp
            elseif strcmp(metric_style, 'velocity') || strcmp(metric_style, 'strain')
                mesh2 = mesh ;

                % Smooth the velocities in space using gptoolbox
                vx = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 1)) ;
                vy = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 2)) ;
                vz = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 3)) ;
                vxs = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], lambda, 'implicit', vx') ;
                vys = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], lambda, 'implicit', vy') ;
                vzs = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], lambda, 'implicit', vz') ;

                mesh2.v = mesh.v + vxs ;

                % Check the vels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if check
                    vxM = reshape(vx, [nU, nV-1]) ;
                    vxsM = reshape(vxs, [nU, nV-1]) ;
                    vysM = reshape(vys, [nU, nV-1]) ;
                    vzsM = reshape(vzs, [nU, nV-1]) ;
                    fig1=figure;
                    colormap bwr; 
                    trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3),...
                        vxM, 'edgecolor', 'none')
                    colorbar(); caxis([-4, 4]); 
                    axis equal; xlabel('x'); ylabel('y'); zlabel('z')
                    title('AP velocity snapshot -- pre-smoothing')
                    fig2=figure;
                    trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3),...
                        vxsM, 'edgecolor', 'none')
                    colorbar(); caxis([-4, 4]); 
                    axis equal; xlabel('x'); ylabel('y'); zlabel('z')
                    title('AP velocity snapshot')
                    waitfor(fig1)
                    waitfor(fig2)

                    trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3),...
                        vysM, 'edgecolor', 'none')
                    colorbar(); caxis([-4, 4]); 
                    axis equal; xlabel('x'); ylabel('y'); zlabel('z')
                    title('lateral velocity snapshot')
                    waitfor(gcf)
                    trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3),...
                        vzsM, 'edgecolor', 'none')
                    colorbar(); caxis([-4, 4]); 
                    axis equal; xlabel('x'); ylabel('y'); zlabel('z')
                    title('DV velocity snapshot')
                    waitfor(gcf)
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end

            V1 = mesh2.v ;
            F1 = mesh2.f ;
            U1 = mesh2.u ;

            %% Construct Topolgical Structure Tools ===============================
            % eg = metricStrainSPhiGridMesh(mesh, mesh2) ;

            % Metric for mesh
            g0cell = inducedMetric(F0, V0, U0) ;
            % Convert metric tensor to #Faces x 2 x 2
            g0 = zeros(size(g0cell, 1), 2, 2) ;
            for qq = 1:size(g0cell, 1) 
                g0(qq, :, :) = g0cell{qq} ;
            end

            % Metric for mesh2
            g1cell = inducedMetric(F1, V1, U1) ;
            % Convert metric tensor to #Faces x 3
            g1 = zeros(size(g1cell, 1), 2, 2) ;
            for qq = 1:size(g1cell, 1)
                g1(qq, :, :) = g1cell{qq} ;
            end

            % Metric strain
            if strcmp(metric_style, 'strain')
                eg = zeros(size(g0)) ;
                for qq = 1:size(g0cell, 1)
                    g0q = g0cell{qq} ;
                    g1q = g1cell{qq} ;
                    eg(qq, :, :) = 0.5 * (inv(g0q) * (g1q - g0q)) ;
                end
            else
                eg = g1 - g0 ;
            end
            
            % Build deformation tensor over mesh, whose elements are the change in
            % squared length in each direction (u.u,u.v,v.v) 
            tg = 0 * eg ;
            % Grab which bonds are along u and which along v
            mesh.nU = nU ;
            [labels, dbonds, tST] = labelRectilinearMeshBonds(mesh) ;
            du = dbonds.baseSpace.u ;
            dv = dbonds.baseSpace.v ;
            for qq = 1:size(g1cell, 1)
                gq = squeeze(eg(qq, :, :)) ;
                
                % NORMALIZATION 
                denom = [1 1 1 1] ;
                % % Frobenius norm of the metric
                % ggqq = squeeze(g0(qq, :, :)) ;
                % % denom = trace(gq .* ggqq') ;
                % 
                % denom(1) = du(qq, :) * ggqq * du(qq, :)' ;
                % denom(4) = dv(qq, :) * ggqq * dv(qq, :)' ;
                % % use fact that du and dv are orthogonal
                % denom(2) = sqrt(denom(1) * denom(4));
                % denom(3) = denom(2) ;
                % % Could consider norming by Euclidean norm of the three
                % % components so they are weighted equally
                
                if strcmp(metric_style, 'strain')
                    tg(qq, 1, 1) = du(qq, :) * gq * du(qq, :)' ;
                    tg(qq, 1, 2) = du(qq, :) * gq * dv(qq, :)' ;
                    tg(qq, 2, 1) = dv(qq, :) * gq * du(qq, :)' ;
                    tg(qq, 2, 2) = dv(qq, :) * gq * dv(qq, :)' ;
                else
                    tg(qq, 1, 1) = (du(qq, :) * gq * du(qq, :)') ;
                    tg(qq, 1, 2) = (du(qq, :) * gq * dv(qq, :)') ;
                    tg(qq, 2, 1) = (dv(qq, :) * gq * du(qq, :)') ;
                    tg(qq, 2, 2) = (dv(qq, :) * gq * dv(qq, :)') ;
                end
                
            end

            % save the metric strain
            readme.eg = ['strain of map from 2D to 3D measured ', ...
                'by ' metric_style ];
            readme.tg = ['strain of 3D bond lengths along s,phi ', ...
                'directions measured by ' metric_style];
            disp(['saving ', gstrainFn])
            save(gstrainFn, 'eg', 'tg', 'readme')
        else
            % load the metric strain
            load(gstrainFn, 'tg')
        end        

        %% Plot the metric components on trisurf
        % denom = sqrt(tg(:, 1, 1) .* tg(:, 2, 2)) ;
        labels = {'$\partial_t g_{ss}$', ...
            '$\partial_t g_{s\phi}$', ...
            '$\partial_t g_{\phi\phi}$'} ;
        strainlabels = {'$\dot{\varepsilon}_{ss}$', ...
            '$\dot{\varepsilon}_{s\phi}$', ...
            '$\dot{\varepsilon}_{\phi\phi}$'} ;
        glab = {'gss', 'gsphi', 'gphiphi'} ;
        gelem = [1 2 4] ;
    
        % consider each metric element
        for qq = 1:3
            colors = tg(:, gelem(qq)) ;
            set(gcf, 'visible', 'off') ;
            colormap bwr
            trisurf(F0, V0(:, 1), V0(:, 2), V0(:, 3), ...
                colors, 'edgecolor', 'none')
            axis equal
            cb = colorbar() ;
                
            if strcmp(metric_style, 'mesh')                 
                caxis([-0.5, 0.5])
                title(['surface deformation rate, ', labels{qq}, ...
                    ': $t=$' sprintf('%03d', tp - tfold)], ...
                    'Interpreter', 'Latex')            
                ylabel(cb, labels{qq}, 'Interpreter', 'Latex')
            elseif strcmp(metric_style, 'velocity')                
                caxis([-0.5, 0.5])
                title(['deformation rate, ', labels{qq}, ...
                    ': $t=$' sprintf('%03d', tp - tfold)], ...
                    'Interpreter', 'Latex')
                ylabel(cb, labels{qq}, 'Interpreter', 'Latex')
            elseif strcmp(metric_style, 'strain')  
                caxis([-0.01, 0.01])      
                title(['strain rate, ', strainlabels{qq}, ...
                    ': $t=$' sprintf('%03d', tp - tfold)], ...
                    'Interpreter', 'Latex')            
                ylabel(cb, strainlabels{qq}, 'Interpreter', 'Latex')
            end
            xlabel('AP position, [$\mu$m]', 'Interpreter', 'Latex')
            ylabel('lateral position, [$\mu$m]', 'Interpreter', 'Latex')
            zlabel('DV position, [$\mu$m]', 'Interpreter', 'Latex')
            xlim(xyzlim(1, :))
            ylim(xyzlim(2, :))
            zlim(xyzlim(3, :))

            % Save images
            imdir = imdirs{metric_style_kk} ;
            % dorsal
            view(0, 90)
            fn = fullfile(imdir, [glab{qq} '_dorsal_' metric_style], ...
                sprintf([glab{qq} '_' QS.fileBase.spcutMeshSmRSC '.png'], tp));
            if ~exist(fn, 'file') || overwrite
                saveas(gcf, fn) ;
            end
            % ventral
            view(0, -90)
            fn = fullfile(imdir, [glab{qq} '_ventral_' metric_style], ...
                sprintf([glab{qq} '_' QS.fileBase.spcutMeshSmRSC '.png'], tp));
            if ~exist(fn, 'file') || overwrite
                saveas(gcf, fn) ;
            end
            % lateral left
            view(0, 0)
            fn = fullfile(imdir, [glab{qq} '_lateral_left_' metric_style], ...
                sprintf([glab{qq} '_' QS.fileBase.spcutMeshSmRSC '.png'], tp));
            if ~exist(fn, 'file') || overwrite
                saveas(gcf, fn) ;
            end
            % dorsal
            view(0, 180)
            fn = fullfile(imdir, [glab{qq} '_lateral_right_' metric_style], ...
                sprintf([glab{qq} '_' QS.fileBase.spcutMeshSmRSC '.png'], tp));
            if ~exist(fn, 'file') || overwrite
                saveas(gcf, fn) ;
            end
            clf
        end
        close all
    end
end


%% Compare size of each contribution visually
% Consider each kind of metric strain measurement
for metric_style_kk = 3 
    metric_style = metric_styles{metric_style_kk} ;
    if strcmp(metric_style, 'mesh')
        gstrainBase = QS.fullFileBase.gstrainMesh ;
    elseif strcmp(metric_style, 'velocity')
        gstrainBase = QS.fullFileBase.gstrainVel;
    elseif strcmp(metric_style, 'strain')
        gstrainBase = QS.fullFileBase.gstrainRate;
    end
    tps = QS.xp.fileMeta.timePoints(1:end-1) - tfold;
    
    % build arrays of average magnitude of g11, g12, g22
    g11avg = zeros(1, length(tps)) ;
    g12avg = zeros(1, length(tps)) ;
    g22avg = zeros(1, length(tps)) ;
    
    g11std = zeros(1, length(tps)) ;
    g12std = zeros(1, length(tps)) ;
    g22std = zeros(1, length(tps)) ;
    
    g11min = zeros(1, length(tps)) ;
    g12min = zeros(1, length(tps)) ;
    g22min = zeros(1, length(tps)) ;
    
    g11max = zeros(1, length(tps)) ;
    g12max = zeros(1, length(tps)) ;
    g22max = zeros(1, length(tps)) ;
    
    % Build metric from mesh
    for tp = xp.fileMeta.timePoints(1:end-1)
        disp(['t = ' num2str(tp)])
        tidx = xp.tIdx(tp) ;
        
        % define metric strain filename        
        gstrainFn = sprintf(gstrainBase, tp) ;
        
        % load it
        load(gstrainFn, 'tg')
        
        % Store in linear array for this timepoint
        g11avg(tidx) = mean(sqrt(tg(:, 1, 1).^2)) ;
        g12avg(tidx) = mean(sqrt(tg(:, 1, 2).^2)) ;
        g22avg(tidx) = mean(sqrt(tg(:, 2, 2).^2)) ;
        g11std(tidx) = std(sqrt(tg(:, 1, 1).^2)) ;
        g12std(tidx) = std(sqrt(tg(:, 1, 2).^2)) ;
        g22std(tidx) = std(sqrt(tg(:, 2, 2).^2)) ;
        g11min(tidx) = min(sqrt(tg(:, 1, 1).^2)) ;
        g12min(tidx) = min(sqrt(tg(:, 1, 2).^2)) ;
        g22min(tidx) = min(sqrt(tg(:, 2, 2).^2)) ;
        g11max(tidx) = max(sqrt(tg(:, 1, 1).^2)) ;
        g12max(tidx) = max(sqrt(tg(:, 1, 2).^2)) ;
        g22max(tidx) = max(sqrt(tg(:, 2, 2).^2)) ;
    end
    
    % Draw each magnitude with stdevs
    tps2 = [tps, fliplr(tps)] ; 
    alph = 0.1 ;
    close all 
    fill(tps2, [max(0, g11avg - g11std), fliplr(g11avg + g11std)], ...
        red, 'facealpha', alph, 'edgecolor', 'none') ;
    hold on;
    fill(tps2, [max(0, g12avg - g12std), fliplr(g12avg + g12std)], ...
        blue, 'facealpha', alph, 'edgecolor', 'none') ;
    fill(tps2, [max(0, g22avg - g22std), fliplr(g22avg + g22std)], ...
        green, 'facealpha', alph, 'edgecolor', 'none') ;
    plot(tps, g11avg, '-', 'color', red) ;
    plot(tps, g12avg, '-', 'color', blue) ;
    plot(tps, g22avg, '-', 'color', green) ;    
    legend({'$g_{ss}$', '$g_{s\phi}$', '$g_{\phi\phi}$'},...
        'Interpreter', 'latex')
    xlabel('time [min]')
    ylabel('$|\partial_t g_{ij}|$', 'Interpreter', 'Latex')
    
    if strcmp(metric_style, 'mesh')
        title('Contribution of surface metric strain components')
        fn = fullfile(QS.dir.gstrain, ...
            ['gstrain_' metric_style '_magnitudes']) ;
    else
        title('Contribution of strain components')
        fn = fullfile(QS.dir.gstrain, ...
            ['gstrain_' metric_style '_magnitudes']) ;
    end
    saveas(gcf, [fn '.png']) ;
    ylim([0, 0.2])
    saveas(gcf, [fn '_zoom.png']) ;
    close all
    
end

%% Fit each hoop to Fourier series and add fit to saved .mat file
% Consider each kind of metric strain measurement
for metric_style_kk = 3 
    metric_style = metric_styles{metric_style_kk} ;
    if strcmp(metric_style, 'mesh')
        gstrainBase = QS.fullFileBase.gstrainMesh ;
    elseif strcmp(metric_style, 'velocity')
        gstrainBase = QS.fullFileBase.gstrainVel;
    elseif strcmp(metric_style, 'strain')
        gstrainBase = QS.fullFileBase.gstrainRate;
    end
    tps = QS.xp.fileMeta.timePoints(1:end-1) - tfold;
    
    % Build metric from mesh
    for tp = xp.fileMeta.timePoints(1:end-1)
        disp(['t = ' num2str(tp)])
        tidx = xp.tIdx(tp) ;
        
        % define metric strain filename        
        gstrainFn = sprintf(gstrainBase, tp) ;
        
        % load it
        load(gstrainFn, 'tg')
        
        % Fit each hoop to series
        phi = 2*pi * (0:(nV-2)) / (nV - 1) ;
        
        fitres = struct() ;
        for tcomp = [1,2,4]
            % preallocate fit coefficients for saving
            fitcoeffs = zeros(nU, 3) ;
            for ihoop = 1:nU
                ids = ihoop:nV:nU*(nV-1) ;
                thoop = squeeze(tg(ids, tcomp))' ;

                % Option 1 is curvfitting toolbox
                % ft = fittype('a + b*sin((x - shift))', ...
                %     'coefficients', {'a', 'b', 'shift'}) ;
                % mdl = fit(X,Y,ft,'startpoint',[shiftguess,xscaleguess,yscaleguess]);

                % Instead use raw MATLAB 
                bguess = [0, 0.01, 0] ;
                % Function to fit
                fit = @(b,phi)  b(1) + b(2)*sin(phi + b(3)) ;  
                % Least-Squares cost function
                fcn = @(b) sum((fit(b,phi) - thoop).^2);
                % Minimise Least-Squares
                fitcoeffs(ihoop, :) = fminsearch(fcn, bguess) ;
            end
            if tcomp == 1
                fitres.ss = fitcoeffs ;
            elseif tcomp == 2
                fitres.sphi = fitcoeffs ;
            elseif tcomp == 4
                fitres.phiphi = fitcoeffs ;
            end
        end
        
        % Save fitres
        save(gstrainFn, 'tg', 'fitres') ;
                
    end
end

%% Plot the fits

tps = QS.xp.fileMeta.timePoints(1:end-1) - tfold;

% Build metric from mesh
fss = zeros(length(tps), 3) ;
fsv = zeros(length(tps), 3) ;
fvv = zeros(length(tps), 3) ;
for tp = xp.fileMeta.timePoints(1:end-1)
    disp(['t = ' num2str(tp)])
    tidx = xp.tIdx(tp) ;

    % define metric strain filename        
    gstrainFitFn = sprintf(QS.fullFileBase.gstrainRate, tp) ;
    
    % load it
    load(gstrainFitFn, 'fitres')
    fss = fitres.ss(:, 1:2) ;
    fsv = fitres.sphi(:, 1:2) ;
    fvv = fitres.phiphi(:, 1:2) ;
    
    % domain of parameterization 
    xx = (0:(nU-1))/(nU-1) ;
    
    fn = fullfile(QS.dir.gstrainRateIm, sprintf('ess_fit_%06d.png', tp)) ;
    plot(xx, fss(:, 1), '-'); hold on;
    plot(xx, fss(:, 2), '-')
    legend({'$\langle \varepsilon_{\zeta\zeta} \rangle$', ...
        '$M_1(\varepsilon_{\zeta\zeta})$'}, ...
        'Location', 'northeastoutside', 'Interpreter', 'Latex')
    title('$\varepsilon_{\zeta\zeta}$', 'Interpreter', 'latex')
    xlabel('AP position, $\zeta/L$', 'Interpreter', 'latex') 
    ylabel('moments', 'Interpreter', 'latex')
    ylim([-0.025, 0.025])
    saveas(gcf, fn)
    clf
    
    fn = fullfile(QS.dir.gstrainRateIm, sprintf('esphi_fit_%06d.png', tp)) ;
    plot(xx, fsv, '-'); hold on;
    plot(xx, fsv(:, 2), '-');
    legend({'$\langle \varepsilon_{\zeta\phi} \rangle$', ...
        '$M_1(\varepsilon_{\zeta\phi})$'}, ...
        'Location', 'northeastoutside', 'Interpreter', 'Latex')
    title('$\varepsilon_{\zeta\phi}$', 'Interpreter', 'latex')
    xlabel('AP position, $\zeta/L$', 'Interpreter', 'latex') 
    ylabel('moments', 'Interpreter', 'latex')
    ylim([-0.025, 0.025])
    saveas(gcf, fn)
    clf
        
    fn = fullfile(QS.dir.gstrainRateIm, sprintf('ephiphi_fit_%06d.png', tp)) ;
    plot(xx, fvv(:, 1), '-'); hold on;
    plot(xx, fvv(:, 2), '-');
    legend({'$\langle \varepsilon_{\phi\phi} \rangle$', ...
        '$M_1(\varepsilon_{\phi\phi})$'}, ...
        'Location', 'northeastoutside', 'Interpreter', 'Latex')
    title('$\varepsilon_{\phi\phi}$', 'Interpreter', 'latex')
    xlabel('AP position, $\zeta/L$', 'Interpreter', 'latex') 
    ylabel('moments', 'Interpreter', 'latex')
    ylim([-0.025, 0.025])
    saveas(gcf, fn)
    clf
    
end


