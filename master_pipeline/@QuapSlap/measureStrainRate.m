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
lambda_mesh = 0.0 ;
overwrite = false ;
overwriteImages = false ;
preview = true ;
averagingStyle = 'Lagrangian' ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
debug = false ;

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
elseif isfield(options, 'overwrite_ims')
    overwriteImages = options.overwrite_ims ;
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
if isfield(options, 'debug')
    debug = options.debug ;
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
nU = QS.nU ;
nV = QS.nV ;

% Load vertex-based velocity measurements
QS.getVelocityAverage('vv', 'vf')
vertex_vels = QS.velocityAverage.vv ;
face_vels = QS.velocityAverage.vf ;

% Pre-assign timepoints with velocities
tpts = QS.xp.fileMeta.timePoints(1:end-1) ;
tp2do = [tpts(1:10:end), setdiff(tpts, tpts(1:10:end))] ;

% Build metric from mesh
for tp = tp2do
    disp(['t = ' num2str(tp)])
    tidx = QS.xp.tIdx(tp) ;

    % Load current mesh
    tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp)) ;
    mesh = tmp.spcutMeshSmRSC ;
    
    % DEBUG
    % Normalize the zeta to fixed aspect ratio (ar=aspectratio relaxed)
    % mesh.u(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) * mesh.ar ;
    clearvars tmp

    % Define metric strain filename        
    estrainFn = fullfile(strrep(sprintf( ...
        QS.dir.strainRate.measurements, lambda, lambda_mesh), '.', 'p'), ...
        sprintf('strainRate_%06d.mat', tp));
        
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
                
        % Smooth the velocities in space using gptoolbox
        vraw = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, :)) ;
        vs = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
            lambda, 'implicit', vraw) ;   
        % Push vectors onto faces
        [V2F, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
        vf = V2F * vs ;
        
        %% Checking -- debug
        if debug
            % Obtain mean curvature H for checking against trace(b_ij)
            DEC = DiscreteExteriorCalculus(mesh.f, mesh.v) ;
            H3d = sum(mesh.vn .* DEC.laplacian(mesh.v), 2) * 0.5 ;
            H2d = H3d ;
            H2d(nU*(nV-1)+1:(nU*nV)) = H3d(1:nU) ;

            Hf = V2F * H3d ;
            divv = DEC.divergence(vf) ;
            divf = V2F * divv ;
            vss = F2V * vf ;
        end
        
        %% Convert to 2D mesh
        mesh.nU = QS.nU ;
        cutMesh = cutRectilinearCylMesh(mesh) ;
        
        % Compute the strain rate tensor
        disp('Computing covariantDerivative')        
        % Decompose velocity into components and compute cov. derivative
        [v0n, v0t] = resolveTangentNormalVector(cutMesh.f, cutMesh.v, vf) ;
        [~, dvi] = vectorCovariantDerivative(v0t, cutMesh.f, cutMesh.v, cutMesh.u) ;
        dvij = cell(size(dvi)) ;
        for qq = 1:length(dvi)
            dvij{qq} = 0.5 * ( dvi{qq} + dvi{qq}' ) ;
        end
        
        % Compute the second fundamental form
        [gg, ~] = constructFundamentalForms(cutMesh.f, cutMesh.v, cutMesh.u) ;
        [~, bb] = constructFundamentalForms(mesh.f, mesh.v, mesh.u) ;
        
        % Strain rate tensor
        strainrate = cell(size(dvi)) ;
        tre = zeros(size(dvi)) ;
        checkH = zeros(size(dvi)) ;
        checkdiv = zeros(size(dvi)) ;
        for qq = 1:size(dvi,1)
            strainrate{qq} = dvij{qq} - v0n(qq) .* bb{qq} ;
        end
        
        %% Debug -- check results against DEC
        if debug
            for qq = 1:size(dvi,1)
                tre(qq) = trace(inv(gg{qq}) * dvij{qq}) - v0n(qq) .* Hf(qq) ;
                checkH(qq) = trace(inv(gg{qq}) * bb{qq}) ;        % == 2 * Hf(qq) ; 
                checkdiv(qq) = trace(inv(gg{qq}) * dvij{qq}) ;    % == divf(qq) ; 
            end
            
            % Compare directly
            clf;
            subplot(2, 1, 1)
            plot(checkdiv, divf, '.')
            hold on;
            plot(checkdiv, checkdiv, 'k--')
            axis equal
            xlabel('Tr$\nabla_i v_j$', 'Interpreter', 'Latex')
            ylabel('$\nabla \cdot \mathbf{v}$', 'Interpreter', 'Latex')
            subplot(2, 1, 2)
            plot(checkH, 2* Hf, 'o') 
            hold on;
            plot(checkH, checkH, 'k--') 
            axis equal
            xlabel('Tr$b_{ij}$', 'interpreter', 'latex')
            ylabel('$2H$', 'interpreter', 'latex')

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

                %% Check div(v)
                clf
                % set color limit, clim
                clim = max(2*std(abs(checkdiv)), 2*std(abs(divf))) ;
                subplot(2, 2, 1)
                trisurf(triangulation(mesh.f, mesh.v), checkdiv, ...
                    'edgecolor', 'none')
                axis equal; caxis([-clim, clim]); colorbar() ;
                checkDivStr = '$\frac{1}{2}$Tr$\left[g^{-1} \left(\nabla_i v_j + \nabla_j v_i \right)\right]$' ;
                title(checkDivStr, 'interpreter', 'latex')
                subplot(2, 2, 2)
                trisurf(triangulation(mesh.f, mesh.v), divv, 'edgecolor', 'none')
                axis equal; caxis([-clim, clim]); colorbar() ;
                title('$\nabla \cdot \mathbf{v}_{\parallel}$', 'interpreter', 'latex')
                subplot(2, 1, 2)
                trisurf(triangulation(mesh.f, mesh.v), checkdiv - divf, ...
                    'edgecolor', 'none')
                axis equal; caxis([-clim, clim]); colorbar() ;
                title([checkDivStr ' $-\nabla \cdot \mathbf{v}_{\parallel}$'], ...
                    'interpreter', 'latex')
                colormap(bwr)

                %% Check velocities (raw vs averaged)
                clf
                clim = 0.1 ;
                subplot(2, 2, 1)
                trisurf(triangulation(mesh.f, mesh.v), vss(:,1)-vs(:, 1), ...
                    'edgecolor', 'none')
                axis equal; caxis([-clim, clim]); colorbar() ;
                checkDivStr = '$\frac{1}{2}$Tr$\left[g^{-1} \left(\nabla_i v_j + \nabla_j v_i \right)\right]$' ;
                title(checkDivStr, 'interpreter', 'latex')
                xlabel('x'); ylabel('y'); zlabel('z')
                title('$v_x$')
                subplot(2, 2, 2)
                trisurf(triangulation(mesh.f, mesh.v), vss(:,2)-vs(:, 2), 'edgecolor', 'none')
                axis equal; caxis([-clim, clim]); colorbar() ;
                title('$\nabla \cdot \mathbf{v}_{\parallel}$', 'interpreter', 'latex')
                xlabel('x'); ylabel('y'); zlabel('z')
                title('$v_y$')
                subplot(2, 1, 2)
                trisurf(triangulation(mesh.f, mesh.v), vss(:,3)-vs(:, 3), ...
                    'edgecolor', 'none')
                axis equal; caxis([-clim, clim]); colorbar() ;
                title([checkDivStr ' $-\nabla \cdot \mathbf{v}_{\parallel}$'], ...
                    'interpreter', 'latex')
                xlabel('x'); ylabel('y'); zlabel('z')
                title('$v_z$')
                colormap(bwr)
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
        % Metric strain -- separate trace and deviatoric strain comp, angle
        epsilon_zz = zeros(size(strainrate, 1), 1) ;  % strain rate zeta zeta
        epsilon_zp = zeros(size(strainrate, 1), 1) ;  % strain rate zeta phi
        epsilon_pz = zeros(size(strainrate, 1), 1) ;  % strain rate phi zeta
        epsilon_pp = zeros(size(strainrate, 1), 1) ;  % strain rate phi phi
        g_zz = zeros(size(strainrate, 1), 1) ;  % metric tensor zeta zeta
        g_zp = zeros(size(strainrate, 1), 1) ;  % metric tensor zeta phi
        g_pz = zeros(size(strainrate, 1), 1) ;  % metric tensor phi zeta
        g_pp = zeros(size(strainrate, 1), 1) ;  % metric tensor phi phi
        b_zz = zeros(size(strainrate, 1), 1) ;  % 2nd fund form zeta zeta
        b_zp = zeros(size(strainrate, 1), 1) ;  % 2nd fund form zeta phi
        b_pz = zeros(size(strainrate, 1), 1) ;  % 2nd fund form phi zeta
        b_pp = zeros(size(strainrate, 1), 1) ;  % 2nd fund form phi phi
        treps = zeros(size(strainrate, 1), 1) ;  % traceful dilation
        dvtre = zeros(size(strainrate, 1), 1) ;  % deviatoric magnitude
        theta = zeros(size(strainrate, 1), 1) ;  % angle of elongation
        % eigv1 = zeros(size(strainrate, 1), 1) ;  % eigenvalue of smaller direction
        % eigv2 = zeros(size(strainrate, 1), 1) ;  % eigenvalue of larger deformation direction (along theta)
        for qq = 1:size(strainrate, 1)
            eq = strainrate{qq} ;
            gq = gg{qq} ;
            bq = bb{qq} ;
            
            %% Full strain rate for mesh averaging onto vertices
            epsilon_zz(qq) = eq(1, 1) ;
            epsilon_zp(qq) = eq(1, 2) ;
            epsilon_pz(qq) = eq(2, 1) ;
            epsilon_pp(qq) = eq(2, 2) ;
            %% Full metric tensor for mesh averaging onto vertices
            g_zz(qq) = gq(1, 1) ;
            g_zp(qq) = gq(1, 2) ;
            g_pz(qq) = gq(2, 1) ;
            g_pp(qq) = gq(2, 2) ;
            %% 2nd fundamental form for mesh averaging onto vertices
            b_zz(qq) = bq(1, 1) ;
            b_zp(qq) = bq(1, 2) ;
            b_pz(qq) = bq(2, 1) ;
            b_pp(qq) = bq(2, 2) ;
            
            %% Trace / deviator / theta
             [treps(qq), dvtre(qq), theta(qq)] = trace_deviator(eq, gq) ;
            
            % eigensystem for strain rate
            % [evec_e, evals_e] = eig(eq) ;
            % [evals_e, idx] = sort(diag(evals_e)) ;
            % evec_e = evec_e(:, idx) ;
            % pevec = evec_e(:, end) ;
            % theta(qq) = atan2(pevec(2), pevec(1)) ;
            % eigv1(qq) = evals_e(1) ;
            % eigv2(qq) = evals_e(2) ;
            
            % NOTE: I have checked that theta determined via full strain
            % rate tensor is identical to theta determined from deviatoric
            % component
        end
        theta = mod(theta, pi) ;
                
        %% Collate results as DVavg, L, R, D, V
        % Find trace and deviator on vertices instead of faces
        epsilon_zz_vtx = F2V * epsilon_zz ;
        epsilon_zp_vtx = F2V * epsilon_zp ;
        epsilon_pz_vtx = F2V * epsilon_pz ;
        epsilon_pp_vtx = F2V * epsilon_pp ;
        g_zz_vtx = F2V * g_zz ;
        g_zp_vtx = F2V * g_zp ;
        g_pz_vtx = F2V * g_pz ;
        g_pp_vtx = F2V * g_pp ;
        b_zz_vtx = F2V * b_zz ;
        b_zp_vtx = F2V * b_zp ;
        b_pz_vtx = F2V * b_pz ;
        b_pp_vtx = F2V * b_pp ;
        
        for qq = 1:size(epsilon_zz_vtx, 1)
            %% Traceful dilation
            eq = [epsilon_zz_vtx(qq), epsilon_zp_vtx(qq); ...
                  epsilon_pz_vtx(qq), epsilon_pp_vtx(qq)] ;
            gq = [g_zz_vtx(qq), g_zp_vtx(qq); ...
                  g_pz_vtx(qq), g_pp_vtx(qq)] ;
            
            % traceful component -- 1/2 Tr[g^{-1} gdot] = Tr[g^{-1} eps] 
            [treps_vtx(qq), dvtre_vtx(qq), theta_vtx(qq)] = ...
                trace_deviator(eq, gq) ;
        end
        theta_vtx = mod(theta_vtx, pi) ;
        
        %% Store measurements on vertices in grouped arrays
        strainrate_vtx = [epsilon_zz_vtx, epsilon_zp_vtx, ...
            epsilon_pz_vtx, epsilon_pp_vtx] ;
        gg_vtx = [g_zz_vtx, g_zp_vtx, ...
            g_pz_vtx, g_pp_vtx] ;
        bb_vtx = [b_zz_vtx, b_zp_vtx, ...
            b_pz_vtx, b_pp_vtx] ;
        treps_vtx((nU * (nV-1) + 1):nU*nV) = treps_vtx(1:nU) ;
        dvtre_vtx((nU * (nV-1) + 1):nU*nV) = dvtre_vtx(1:nU) ;
        theta_vtx((nU * (nV-1) + 1):nU*nV) = theta_vtx(1:nU) ;
        treps_vtx = reshape(treps_vtx, [nU,nV]) ;
        dvtre_vtx = reshape(dvtre_vtx, [nU,nV]) ;
        theta_vtx = reshape(theta_vtx, [nU,nV]) ;
        
        % Average along DV -- ignore last redudant row at nV
        [dvtre_ap, theta_ap] = ...
            QS.dvAverageNematic(dvtre_vtx(:, 1:nV-1), theta_vtx(:, 1:nV-1)) ;
        treps_ap = mean(treps_vtx(:, 1:nV-1), 2) ;
        
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
        [dvtre_l, theta_l] = ...
            QS.dvAverageNematic(dvtre_vtx(:, left), theta_vtx(:, left)) ;
        treps_l = mean(treps_vtx(:, left), 2) ;
        
        % right quarter
        [dvtre_r, theta_r] = ...
            QS.dvAverageNematic(dvtre_vtx(:, right), theta_vtx(:, right)) ;
        treps_r = mean(treps_vtx(:, right), 2) ;
        
        % dorsal quarter
        [dvtre_d, theta_d] = ...
            QS.dvAverageNematic(dvtre_vtx(:, dorsal), theta_vtx(:, dorsal)) ;
        treps_d = mean(treps_vtx(:, dorsal), 2) ;
        
        % ventral quarter
        [dvtre_v, theta_v] = ...
            QS.dvAverageNematic(dvtre_vtx(:, ventral), theta_vtx(:, ventral)) ;
        treps_v = mean(treps_vtx(:, ventral), 2) ;
        
        % save the metric strain
        readme.strainrate = 'strain rate on faces, epsilon=1/2(nabla_i v_j + nabla_j v_i) - vn b_ij' ;
        readme.treps = 'Tr[g^{-1} epsilon]';
        readme.dvtre = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] )';
        readme.theta = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector';
        readme.dvij = 'symmetrized covariant derivative 0.5 * (d_i v_j + d_j v_i)';
        readme.gg = 'metric tensor on faces';
        readme.bb = 'second fundamental form on faces';
        readme.lambda = 'Laplacian smoothing on velocities' ;
        readme.lambda_mesh = 'Laplacian smoothing on mesh vertices' ;
        readme.strainrate_vtx = 'strain rate on vertices, epsilon=1/2(nabla_i v_j + nabla_j v_i) - vn b_ij' ;
        readme.gg_vtx = 'metric tensor on vertices';
        readme.bb_vtx = 'second fundamental form on vertices';
        readme.treps_vtx = 'Tr[g^{-1} epsilon], on mesh vertices' ;
        readme.dvtre_vtx = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) on mesh vertices';
        readme.theta_vtx = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector on mesh vertices';
        readme.dvtre_ap = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged circumferentially';
        readme.dvtre_l = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on left quarter, on vertices';
        readme.dvtre_r = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on right quarter, on vertices';
        readme.dvtre_d = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on dorsal quarter, on vertices';
        readme.dvtre_v = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on ventral quarter, on vertices';
        readme.theta_ap = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged circumferentially, on vertices';
        readme.theta_l = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on left quarter, on vertices';
        readme.theta_r = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on right quarter, on vertices';
        readme.theta_d = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on dorsal quarter, on vertices';
        readme.theta_v = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on ventral quarter, on vertices';
        readme.treps_ap = 'Tr[g^{-1} epsilon], averaged circumferentially, on vertices';
        readme.treps_l = 'Tr[g^{-1} epsilon], averaged on left quarter, on vertices';
        readme.treps_r = 'Tr[g^{-1} epsilon], averaged on right quarter, on vertices';
        readme.treps_d = 'Tr[g^{-1} epsilon], averaged on dorsal quarter, on vertices';
        readme.treps_v = 'Tr[g^{-1} epsilon], averaged on ventral quarter, on vertices';
        readme.note = 'The pullback space is taken to range from zeta=[0, 1] and phi=[0, 1]' ; 
        disp(['saving ', estrainFn])
        save(estrainFn, 'strainrate', 'treps', 'dvtre', 'theta', ...
            'dvij', 'gg', 'bb', 'lambda', 'lambda_mesh', 'readme', ...
            'dvtre_ap', 'dvtre_l', 'dvtre_r', 'dvtre_d', 'dvtre_v', ...
            'theta_ap', 'theta_l', 'theta_r', 'theta_d', 'theta_v', ...
            'treps_ap', 'treps_l', 'treps_r', 'treps_d', 'treps_v', ...
            'strainrate_vtx', 'treps_vtx', 'dvtre_vtx', 'theta_vtx', ...
            'gg_vtx', 'bb_vtx')
    else
        % Convert to 2D mesh
        mesh.nU = QS.nU ;
        cutMesh = cutRectilinearCylMesh(mesh) ;
    end        
    
    % Plot the result
    options.overwrite = overwriteImages ;
    options.mesh = mesh ;
    options.cutMesh = cutMesh ;
    options.lambda = lambda ;
    options.lambda_mesh = lambda_mesh ;
    options.debug = debug ;
    QS.plotStrainRateTimePoint(tp, options)
end


