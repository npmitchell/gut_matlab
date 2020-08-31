function measurePathlineStrainRate(QS, options)
% measurePathlineStrainRate(QS, options)
%   Query the metric strain rate along lagrangian pathlines.
%   Plot results as kymographs.
%   
%  Coordinate transformation: if mesh shrinks in size, then accumulated xx
%  strain of previous timepoint is diminished via
%       inv(J01f) * strain0 * (inv(J01f).')
%  In the limit that the previous mesh was infinite in U and fixed size in
%  V, then the resulting contribution of the previously accumulated strain
%  to the current timepoint is dominated by e_phiphi, seen as 
%       imagesc(1:nU, fliplr(1:nV), reshape(epp ./ gpp, [nU, nV])')
%  If the mesh grows along U (Zeta) and stays fixed along V (phi), 
%  as is the case for the fly midgut, then the accumulated strain in xx 
%  will be magnified by transforming into the current coordinates. 
%  Given no rotation, then 
%       J01f = Jacobian_{t-1 to t} = [<1, 0; 0, 1] ;
%  If the mesh grows in U by a factor R, then 
%       strain0 = (1, 0; 0, 0) -->
%         inv(J01f) * strain0 * (inv(J01f).') = [1/R^2, 0; 0, 0] ;
%       strain0 = (0, 1; 1, 0) -->
%         inv(J01f) * strain0 * (inv(J01f).') = [0, 1/R; 1/R, 0] ;
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   plot_kymographs : bool
%   plot_kymographs_cumsum : bool
%   plot_gdot_correlations : bool
%   plot_gdot_decomp : bool
% 
% NPMitchell 2020

%% Default options 
overwrite = false ;
overwriteImages = false ;

%% Parameter options
lambda_mesh = 0.002 ;
lambda = 0.01 ; 
debug = false ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
averagingStyle = "Lagrangian" ;
% Load time offset for first fold, t0 -- default pathline t0
QS.t0set() ;
t0 = QS.t0 ;

%% Unpack options & assign defaults
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
end
%% parameter options
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
else
    % default lambda_mesh is equal to lambda 
    lambda_mesh = lambda ;
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
if isfield(options, 't0Pathline')
    t0Pathline = options.t0Pathline ;
else
    t0Pathline = t0 ;
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
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;
if strcmp(averagingStyle, 'Lagrangian')
    sKDir = fullfile(QS.dir.strainRate.root, ...
        strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f'], ...
        lambda, lambda_mesh), '.', 'p'));
else
    error('Have not implemented strain rate measurements based on simple averaging')
end
folds = load(QS.fileName.fold) ;
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;

%% Colormap
bwr256 = bluewhitered(256) ;

%% load from QS
if doubleResolution
    nU = QS.nU * 2 - 1 ;
    nV = QS.nV * 2 - 1 ;
else
    nU = QS.nU ;
    nV = QS.nV ;    
end

% We relate the normal velocities to the divergence / 2 * H.
tps = QS.xp.fileMeta.timePoints(1:end-1) - t0;

% Unit definitions for axis labels
unitstr = [ '[1/' QS.timeUnits ']' ];
vunitstr = [ '[' QS.spaceUnits '/' QS.timeUnits ']' ];
    
% DONE WITH PREPARATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load pathlines to build Kymographs along pathlines
QS.loadPullbackPathlines(t0Pathline, 'vertexPathlines')
vP = QS.pathlines.vertices ;

% Output directory is inside StrainRate dir
sKPDir = fullfile(sKDir, sprintf('pathline_%04dt0', t0Pathline)) ;
outdir = fullfile(sKPDir, 'measurements') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
% Data for kinematics on meshes (defined on vertices)
mdatdir = fullfile(sKDir, 'measurements') ;

% Load Lx, Ly by loadingPIV. 
QS.loadPIV()
Xpiv = QS.piv.raw.x ;
Ypiv = QS.piv.raw.y ;

% Also need velocities to advect mesh
% QS.loadVelocityAverage('vv')
% vPIV = QS.velocityAverage.vv ;

% Discern if piv measurements are done on a double covering or the meshes
if strcmp(QS.piv.imCoords(end), 'e')
    doubleCovered = true ;
end

% Compute or load all timepoints
strain = zeros(size(vP.vX, 2) * size(vP.vX, 3), 4) ; 
strain(:, 2) = 1e4 ;
strain(:, 3) = 1e4 ;
for tp = QS.xp.fileMeta.timePoints(1:end-1)
    % tp = 75 ;
    close all
    disp(['t = ' num2str(tp)])
    tidx = QS.xp.tIdx(tp) ;
    QS.setTime(tp) ;
    
    % Check for timepoint measurement on disk, on mesh vertices 
    estrainFn = fullfile(outdir, sprintf('strainRate_%06d.mat', tp)) ;
    
    if overwrite || ~exist(estrainFn, 'file') 
        % Load timeseries measurements defined on mesh vertices 
        srfnMesh = fullfile(mdatdir, sprintf('strainRate_%06d.mat', tp)) ;
        try
            load(srfnMesh, 'strainrate_vtx', 'gg_vtx', 'dx_vtx', 'dy_vtx') 
        catch
            msg = 'Run QS.measurePathlineStrainRate() ' ;
            msg = [msg 'with lambdas=(mesh,lambda,err)=('] ;
            msg = [msg num2str(lambda_mesh) ','] ;
            msg = [msg num2str(lambda) ','] ;
            msg = [msg ' before running ', ...
                    'QS.measurePathlineStrainRate()'] ;
            error(msg)
        end
        %% Interpolate StrainRate from vertices onto pathlines
        xx = vP.vX(tidx, :, :) ;
        yy = vP.vY(tidx, :, :) ;
        XY = [xx(:), yy(:)] ;
        Lx = vP.Lx(tidx) ;
        Ly = vP.Ly(tidx) ;
        options.Lx = Lx ;
        options.Ly = Ly ;
        XY = QS.doubleToSingleCover(XY, Ly) ;
        
        %% Recall strain rate at gridded vertices
        % strainrate from vertices to pathlines
        exx = strainrate_vtx(:, 1) ;
        exy = strainrate_vtx(:, 2) ;
        eyx = strainrate_vtx(:, 3) ;
        eyy = strainrate_vtx(:, 4) ;
        exx(nU*(nV-1)+1:nU*nV) = exx(1:nU) ;
        exy(nU*(nV-1)+1:nU*nV) = exy(1:nU) ;
        eyx(nU*(nV-1)+1:nU*nV) = eyx(1:nU) ;
        eyy(nU*(nV-1)+1:nU*nV) = eyy(1:nU) ;
        ezz = QS.interpolateOntoPullbackXY(XY, exx, options) ;
        ezp = QS.interpolateOntoPullbackXY(XY, exy, options) ;
        epz = QS.interpolateOntoPullbackXY(XY, eyx, options) ;
        assert(all(ezp == epz))
        epp = QS.interpolateOntoPullbackXY(XY, eyy, options) ;
        strainrate = [ezz, ezp, epz, epp] ;
        
        % Metric from vertices to pathlines
        gxx = gg_vtx(:, 1) ;
        gxy = gg_vtx(:, 2) ;
        gyx = gg_vtx(:, 3) ;
        gyy = gg_vtx(:, 4) ;
        gxx(nU*(nV-1)+1:nU*nV) = gxx(1:nU) ;
        gxy(nU*(nV-1)+1:nU*nV) = gxy(1:nU) ;
        gyx(nU*(nV-1)+1:nU*nV) = gyx(1:nU) ;
        gyy(nU*(nV-1)+1:nU*nV) = gyy(1:nU) ;
        gzz = QS.interpolateOntoPullbackXY(XY, gxx, options) ;
        gzp = QS.interpolateOntoPullbackXY(XY, gxy, options) ;
        gpz = QS.interpolateOntoPullbackXY(XY, gyx, options) ;
        gpp = QS.interpolateOntoPullbackXY(XY, gyy, options) ;
        assert(all(gzp == gpz))
        gg = [gzz, gzp, gpz, gpp] ;
        
        % dx and dy scalings
        dx = dx_vtx ;
        dy = dy_vtx ;
        dx(nU*(nV-1)+1:nU*nV) = dx(1:nU) ;
        dy(nU*(nV-1)+1:nU*nV) = dy(1:nU) ;
        dz = QS.interpolateOntoPullbackXY(XY, dx, options) ;
        dp = QS.interpolateOntoPullbackXY(XY, dy, options) ;
        
        %% Compute the strain rate trace and deviator at the pathlines
        for qq = 1:size(exx, 1)
            eq = [ezz(qq), ezp(qq); ...
                  epz(qq), epp(qq)] ;
            gq = [gzz(qq), gzp(qq); ...
                  gpz(qq), gpp(qq)] ;
            
            % traceful component is 1/2 Tr[g^{-1} gdot] = Tr[g^{-1} eps] 
            % deviatoric component is eps - 1/2 (traceful component) g
            [tre(qq), dev(qq), theta(qq)] = ...
                traceDeviatorPullback(eq, gq, dz(qq), dp(qq)) ;
        end
        
        %% Accumulate strain rate into STRAIN        
        % Transform the previous strain rate into current basis
        if tp > QS.xp.fileMeta.timePoints(1) 
            % Load previous mesh
            % DEBUG
            tp_prev = QS.xp.fileMeta.timePoints(tidx-1) ;
            
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
            mesh1 = tmp.spcutMeshSmRS ;
            clearvars tmp
            
            % Time increment
            % DEBUG
            dt = QS.timeInterval * (tp_prev - tp) ;
            
            % DEBUG
            % Normalize the zeta to fixed aspect ratio (ar=aspectratio relaxed)
            mesh1.u(:, 1) = mesh1.u(:, 1) / max(mesh1.u(:, 1)) * mesh1.ar ;
            
            %--------------------------------------------------------------------------
            % Find where mesh1.u would be at previous timepoint by tracing
            % back along the flow. 
            % If we want to define a transformation between the frame of time point 1
            % and into the frame of time point 1 that ends up recapitulating the mesh
            % at time point 1, we must ask "Where were the vertices of the mesh at time
            % point 1 back during time point 0?".  We COULD extract this information
            % using the motion of the material points, but it would be more accurate to
            % use the PIV fields from which the material points were built.
            % The PIV consists of velocities vPIV and evaluation positions (Xpiv,Ypiv).
            %--------------------------------------------------------------------------
            % qq is the previous timepoint
            qq = tidx-1 ;

            % The 'inverse map' 
            % 1. Interpolate velocities of time t=0 at their advected 
            %    locations in t=1. 
            %
            %   .---->*     .<----*
            %     qq=.       qq+1=*
            % 
            %    Since the advected x0+u_qq,y0+v_qq are spatially unstructured, we
            %    interpolate t=0 velocity at the advected positions,
            %    ie what velocities took XY0 to XY1, evaluate the velocities 
            %    that pushed to the next positions AT the next positions XY1,
            %    and pull back XY1 along those velocities to get displaced coordinates 
            %    DXY0 which will land on XY1 when moved along t0's flow.
            %
            uu = QS.piv.raw.u_filtered{qq} ;
            vv = QS.piv.raw.v_filtered{qq} ;

            % Load Lx, Ly for t=1, which are the extents of the domain (needed for
            % periodic boundary consitions and to keep all the points in the domain of
            % interpolation). 
            x01 = Xpiv{tidx}(:) + uu(:) ;
            y01 = Ypiv{tidx}(:) + vv(:) ;
            [xa, ya] = QS.clipXY(x01, y01, Lx, Ly) ;
            
            ui = scatteredInterpolant(xa, ya, uu(:), 'natural', 'nearest') ;
            vi = scatteredInterpolant(xa, ya, vv(:), 'natural', 'nearest') ;
            
            % 2. Evaluate at mesh1 vertices to pull them back along velocity to t=0.
            dx = ui(mesh1.u(:, 1), mesh1.u(:, 2)) ;
            dy = vi(mesh1.u(:, 1), mesh1.u(:, 2)) ;
            
            % 3. Pull mesh vertices back to t=0
            Xqq = mesh1.u(:, 1) - dx(:) ;
            Yqq = mesh1.u(:, 2) - dy(:) ;

            % The locations of the rectilinear vertices at time 0
            % 'Deformed Vertices from t1 at t0'
            DXY10 = [Xqq, Yqq ] ;

            % Check the deformed mesh that will advect into the gridded
            % mesh at t+1
            % clf
            % scatter(x01(:), y01(:), 1, uu(:), 'filled');
            % hold on;
            % quiver(Xqq, Yqq, dx(:), dy(:), 0, 'c')
            % triplot(triangulation(F1, DXY10));
            % axis equal            

            % Construct the Jacobian matrix on each mesh face
            % J01 = jacobian2Dto2DMesh(mesh1.u, mesh0.u, mesh1.f) ;
            J10 = jacobian2Dto2DMesh(mesh0.u, mesh1.u, F);

            %% ------------------------------------------------------------
            % Find which jacobians to use for each pathline point
            tri = triangulation(mesh1.f, mesh1.u) ;
            umax = max(mesh1.u(:, 1)) ;
            vmax = max(mesh1.u(:, 2)) ;
            if strcmp(QS.piv.imCoords, 'sp_sme')
                im = imfinfo(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
                im = [im.Height, im.Width] ;
                if im(1) ~= im(2)
                    error('check that dimensions are in correct order here')
                end
            else
                error('Handle this case for imCoords here')
            end
            uv = QS.XY2uv(im, XY, doubleCovered, umax, vmax) ;
            uv(:, 2) = mod(uv(:, 2), vmax) ; 
            uv(:, 1) = max(uv(:, 1), 1e-13) ; 
            uv(:, 1) = min(uv(:, 1), umax-1e-13) ; 
            fieldfaces = pointLocation(tri, uv) ;
            assert(all(min(fieldfaces) > 0))
            assert(length(fieldfaces) == length(ezz))
            
            %% ------------------------------------------------------------------------
            % Consider each pathline, add the strainRate * dt to the
            % accumulated strain from previous timepoint, which is called 
            % strain0. 
            % Accumulate the strain via implicit Euler (backward Euler
            % scheme)
            % Transform as a (0,2)-tensor (NOTICE THE MATRIX INVERSES)
            for qq = 1:size(fieldfaces,1)
                strain0 = [strain(qq, 1), strain(qq, 2); ...
                           strain(qq, 3), strain(qq, 4)] ;
                
                strainrateP = [ezz(qq), ezp(qq); epz(qq), epp(qq)] ;
                try
                    J01f = J01{fieldfaces(qq)} ;
                    strainM{qq} = inv(J01f) * strain0 * (inv(J01f).') + ...
                        dt * strainrateP ;
                catch
                    error('Ensure that all uv lie in the mesh.u')
                end
            end
            
            % Convert strain from cell (needed for face Jacobians) to array
            for qq = 1:size(ezz, 1)
                strain(qq, 1) = strainM{qq}(1, 1) ;
                strain(qq, 2) = strainM{qq}(1, 2) ;
                strain(qq, 3) = strainM{qq}(2, 1) ;
                strain(qq, 4) = strainM{qq}(2, 2) ; 
            end
        else
            % Load mesh for this timepoint
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
            mesh1 = tmp.spcutMeshSmRS ;
            clearvars tmp
            
            % Normalize the zeta to fixed aspect ratio (ar=aspectratio relaxed)
            mesh1.u(:, 1) = mesh1.u(:, 1) / max(mesh1.u(:, 1)) * mesh1.ar ;
            
            dt = QS.timeInterval ;
            % Note: strainrate is already interpolated onto pathlines
            strain = dt * strainrate ;
        end
        
        %% Trace/Determinant of strain 
        strain_tr = zeros(size(strain, 1), 1) ;
        strain_dv = zeros(size(strain, 1), 1) ;
        strain_th = zeros(size(strain, 1), 1) ;
        for qq = 1:size(strain, 1)
            eq = [strain(qq, 1), strain(qq, 2); ...
                  strain(qq, 3), strain(qq, 4)] ;
            gq = [gzz(qq), gzp(qq); ...
                  gpz(qq), gpp(qq)] ;
            [strain_tr(qq), strain_dv(qq), strain_th(qq)] = ...
                traceDeviatorPullback(eq, gq, dz(qq), dp(qq));
        end
        
        %% CHECK integration 
        if debug
            % Map intensity from dev and color from the theta
            close all
            pm256 = phasemap(256) ;
            indx = max(1, round(mod(2*strain_th(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
            colors = pm256(indx, :) ;
            devKclipped = min(strain_dv / 0.1, 1) ;
            colorsM = devKclipped(:) .* colors ;
            colorsM = reshape(colorsM, [nU, nV, 3]) ;
            imagesc(1:nU, 1:nV, permute(colorsM, [2, 1, 3]))
            caxis([0, 0.1])
            pause(1)
        end
                
        %% OPTION 1: simply reshape, tracing each XY pathline pt to its t0
        % % grid coordinate
        tre = reshape(tre, [nU, nV]) ;
        dev = reshape(dev, [nU, nV]) ;
        theta = reshape(theta, [nU, nV]) ;
        strain_tr = reshape(strain_tr, [nU, nV]) ;
        strain_dv = reshape(strain_dv, [nU, nV]) ;
        strain_th = reshape(strain_th, [nU, nV]) ;
        
        %% OPTION 2: the following regrids onto original XY coordinates,
        % rendering the process of following pathlines moot. 
        % Average into AP bins and take mean along 1/4 DV hoop arcs
        % if doubleCovered
        %     vminmax = [0.25 * Ly, 0.75 * Ly] ;
        % else
        %     vminmax = [1, Ly] ;
        % end
        %
        % Note the transposition: to plot as APDV, imshow(m')
        % HH = binData2dGrid([XY, HH], [1,Lx], vminmax, nU, nV) ;
        % gdot = binData2dGrid([XY, gdot], [1,Lx], vminmax, nU, nV) ;
        % divv = binData2dGrid([XY, divv], [1,Lx], vminmax, nU, nV) ;
        % veln = binData2dGrid([XY, veln], [1,Lx], vminmax, nU, nV) ;
        % H2vn = binData2dGrid([XY, H2vn], [1,Lx], vminmax, nU, nV) ;
        
        %% Average strainRATE along pathline DV hoops
        % Average along DV -- ignore last redudant row at nV
        [dev_ap, theta_ap] = ...
            QS.dvAverageNematic(dev(:, 1:nV-1), theta(:, 1:nV-1)) ;
        tre_ap = mean(tre(:, 1:nV-1), 2) ;
        
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
        [dev_l, theta_l] = ...
            QS.dvAverageNematic(dev(:, left), theta(:, left)) ;
        tre_l = mean(tre(:, left), 2) ;
        
        % right quarter
        [dev_r, theta_r] = ...
            QS.dvAverageNematic(dev(:, right), theta(:, right)) ;
        tre_r = mean(tre(:, right), 2) ;
        
        % dorsal quarter
        [dev_d, theta_d] = ...
            QS.dvAverageNematic(dev(:, dorsal), theta(:, dorsal)) ;
        tre_d = mean(tre(:, dorsal), 2) ;
        
        % ventral quarter
        [dev_v, theta_v] = ...
            QS.dvAverageNematic(dev(:, ventral), theta(:, ventral)) ;
        tre_v = mean(tre(:, ventral), 2) ;
        
        %% Average STRAIN (accumulated strain) along DV 
        % Average along DV -- ignore last redudant row at nV
        [strain_dv_ap, strain_th_ap] = ...
            QS.dvAverageNematic(strain_dv(:, 1:nV-1), strain_th(:, 1:nV-1)) ;
        strain_tr_ap = mean(strain_tr(:, 1:nV-1), 2) ;
        
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
        [strain_dv_l, strain_th_l] = ...
            QS.dvAverageNematic(strain_dv(:, left), strain_th(:, left)) ;
        strain_tr_l = mean(strain_tr(:, left), 2) ;
        % right quarter
        [strain_dv_r, strain_th_r] = ...
            QS.dvAverageNematic(strain_dv(:, right), strain_th(:, right)) ;
        strain_tr_r = mean(strain_tr(:, right), 2) ;
        % dorsal quarter
        [strain_dv_d, strain_th_d] = ...
            QS.dvAverageNematic(strain_dv(:, dorsal), strain_th(:, dorsal)) ;
        strain_tr_d = mean(strain_tr(:, dorsal), 2) ;
        % ventral quarter
        [strain_dv_v, strain_th_v] = ...
            QS.dvAverageNematic(strain_dv(:, ventral), strain_th(:, ventral)) ;
        strain_tr_v = mean(strain_tr(:, ventral), 2) ;
        
        %% Check result
        if debug
            close all
            scatter(1:nU, strain_th_ap, 20 * strain_dv_ap / max(strain_dv_ap),...
                strain_dv_ap, 'filled')
            pause(2)
        end
        
        % save the metric strain
        readme.strainrate = 'strain rate tensor interpolated onto pathlines' ;
        readme.gg = 'metric tensor interpolated onto pathlines' ;
        readme.tre = 'Tr[g^{-1} epsilon]';
        readme.dev = 'sqrt( Tr[g^{-1} deviator[epsilon] g^{-1} deviator[epsilon]] ) on mesh vertices';
        readme.theta = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector on mesh vertices';
        readme.dev_ap = 'sqrt( Tr[g^{-1} deviator[epsilon] g^{-1} deviator[epsilon]] ) averaged circumferentially';
        readme.dev_l = 'sqrt( Tr[g^{-1} deviator[epsilon] g^{-1} deviator[epsilon]] ) averaged on left quarter, on vertices';
        readme.dev_r = 'sqrt( Tr[g^{-1} deviator[epsilon] g^{-1} deviator[epsilon]] ) averaged on right quarter, on vertices';
        readme.dev_d = 'sqrt( Tr[g^{-1} deviator[epsilon] g^{-1} deviator[epsilon]] ) averaged on dorsal quarter, on vertices';
        readme.dev_v = 'sqrt( Tr[g^{-1} deviator[epsilon] g^{-1} deviator[epsilon]] ) averaged on ventral quarter, on vertices';
        readme.theta_ap = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged circumferentially, on vertices';
        readme.theta_l = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on left quarter, on vertices';
        readme.theta_r = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on right quarter, on vertices';
        readme.theta_d = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on dorsal quarter, on vertices';
        readme.theta_v = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on ventral quarter, on vertices';
        readme.tre_ap = 'Tr[g^{-1} epsilon], averaged circumferentially, on vertices';
        readme.tre_l = 'Tr[g^{-1} epsilon], averaged on left quarter, on vertices';
        readme.tre_r = 'Tr[g^{-1} epsilon], averaged on right quarter, on vertices';
        readme.tre_d = 'Tr[g^{-1} epsilon], averaged on dorsal quarter, on vertices';
        readme.tre_v = 'Tr[g^{-1} epsilon], averaged on ventral quarter, on vertices';
        
        readme.strain = 'integrated strain tensor' ;
        readme.strain_tr = 'integrated strain trace Tr[g^{-1} epsilon]';        
        readme.strain_dv = 'integrated strain deviator magnitude sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) on mesh vertices';
        readme.strain_th = 'integrated strain deviator angle -- arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector on mesh vertices';
        readme.strain_dv_ap = 'integrated sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged circumferentially';
        readme.strain_dv_l = 'integrated sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on left quarter, on vertices';
        readme.strain_dv_r = 'integrated sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on right quarter, on vertices';
        readme.strain_dv_d = 'integrated sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on dorsal quarter, on vertices';
        readme.strain_dv_v = 'integrated sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on ventral quarter, on vertices';
        readme.strain_th_ap = 'integrated strain deviator angle -- arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged circumferentially, on vertices';
        readme.strain_th_l = 'integrated strain deviator angle -- arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on left quarter, on vertices';
        readme.strain_th_r = 'integrated strain deviator angle -- arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on right quarter, on vertices';
        readme.strain_th_d = 'integrated strain deviator angle -- arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on dorsal quarter, on vertices';
        readme.strain_th_v = 'integrated strain deviator angle -- arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on ventral quarter, on vertices';
        readme.strain_tr_ap = 'integrated strain trace Tr[g^{-1} epsilon], averaged circumferentially, on vertices';
        readme.strain_tr_l = 'integrated strain trace Tr[g^{-1} epsilon], averaged on left quarter, on vertices';
        readme.strain_tr_r = 'integrated strain trace Tr[g^{-1} epsilon], averaged on right quarter, on vertices';
        readme.strain_tr_d = 'integrated strain trace Tr[g^{-1} epsilon], averaged on dorsal quarter, on vertices';
        readme.strain_tr_v = 'integrated strain trace Tr[g^{-1} epsilon], averaged on ventral quarter, on vertices';
        
        readme.note = 'Evaluated for Lagrangian paths. The pullback space is taken to range from zeta=[0, 1] and phi=[0, 1]' ; 
        disp(['saving ', estrainFn])
        save(estrainFn, 'strainrate', 'gg', 'strain', 'readme', ...
            'tre', 'dev', 'theta', ...
            'tre_ap', 'tre_l', 'tre_r', 'tre_d', 'tre_v', ...
            'dev_ap', 'dev_l', 'dev_r', 'dev_d', 'dev_v', ...
            'theta_ap', 'theta_l', 'theta_r', 'theta_d', 'theta_v', ...
            'strain_tr', 'strain_dv', 'strain_th', ...
            'strain_dv_ap', 'strain_dv_l', 'strain_dv_r', ...
            'strain_dv_d', 'strain_dv_v', ...
            'strain_tr_ap', 'strain_tr_l', 'strain_tr_r', ...
            'strain_tr_d', 'strain_tr_v', ...
            'strain_th_ap', 'strain_th_l', 'strain_th_r', ...
            'strain_th_d', 'strain_th_v')
    else
        load(estrainFn, 'strain')
        % Load mesh 
        tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
        mesh1 = tmp.spcutMeshSmRS ;
        clearvars tmp
        % Normalize the zeta to fixed aspect ratio (ar=aspectratio relaxed)
        mesh1.u(:, 1) = mesh1.u(:, 1) / max(mesh1.u(:, 1)) * mesh1.ar ;
    end
    
    % Plot the result
    options.overwrite = overwriteImages ;
    options.cutMesh = mesh1 ;
    options.lambda = lambda ;
    options.lambda_mesh = lambda_mesh ;
    options.debug = debug ;
    options.t0Pathline = t0Pathline ;
    QS.plotPathlineStrainRateTimePoint(tp, options)
    
end
disp('done with measuring pathline strain rate')


%% Combine DV-averaged profiles into kymographs
apKymoFn = fullfile(outdir, 'apKymographPathlineStrainRate.mat') ;
lKymoFn = fullfile(outdir, 'leftKymographPathlineStrainRate.mat') ;
rKymoFn = fullfile(outdir, 'rightKymographPathlineStrainRate.mat') ;
dKymoFn = fullfile(outdir, 'dorsalKymographPathlineStrainRate.mat') ;
vKymoFn = fullfile(outdir, 'ventralKymographPathlineStrainRate.mat') ;
files_exist = exist(apKymoFn, 'file') && ...
    exist(lKymoFn, 'file') && exist(rKymoFn, 'file') && ...
    exist(dKymoFn, 'file') && exist(vKymoFn, 'file') ;
if ~files_exist || true
    disp('Compiling kymograph data to save to disk...')
    for tp = QS.xp.fileMeta.timePoints(1:end-1)
        close all
        tidx = QS.xp.tIdx(tp) ;

        % Check for timepoint measurement on disk
        srfn = fullfile(outdir, sprintf('strainRate_%06d.mat', tp))   ;

        % Load timeseries measurements
        load(srfn, 'tre_ap', 'tre_l', 'tre_r', 'tre_d', 'tre_v', ...
            'dev_ap', 'dev_l', 'dev_r', 'dev_d', 'dev_v', ...
            'theta_ap', 'theta_l', 'theta_r', 'theta_d', 'theta_v', ...
            'strain_tr_ap', 'strain_tr_l', 'strain_tr_r', ...
            'strain_tr_d', 'strain_tr_v', ...
            'strain_th_ap', 'strain_th_l', 'strain_th_r', ....
            'strain_th_d', 'strain_th_v', ...
            'strain_dv_ap', 'strain_dv_l', 'strain_dv_r', ...
            'strain_dv_d', 'strain_dv_v') ;

        %% Store in matrices
        % dv averaged
        tr_apM(tidx, :) = tre_ap ;
        dv_apM(tidx, :) = dev_ap ;
        th_apM(tidx, :) = theta_ap ;

        % left quarter
        tr_lM(tidx, :) = tre_l ;
        dv_lM(tidx, :) = dev_l ;
        th_lM(tidx, :) = theta_l ;

        % right quarter
        tr_rM(tidx, :) = tre_r ;
        dv_rM(tidx, :) = dev_r ;
        th_rM(tidx, :) = theta_r ;

        % dorsal quarter
        tr_dM(tidx, :) = tre_d ;
        dv_dM(tidx, :) = dev_d ;
        th_dM(tidx, :) = theta_d ;

        % ventral quarter
        tr_vM(tidx, :) = tre_v ;
        dv_vM(tidx, :) = dev_v ;
        th_vM(tidx, :) = theta_v ;

        %% Store accumulated strain in matrices
        % dv averaged
        str_apM(tidx, :) = strain_tr_ap ;
        sdv_apM(tidx, :) = strain_dv_ap ;
        sth_apM(tidx, :) = strain_th_ap ;

        % left quarter
        str_lM(tidx, :) = strain_tr_l ;
        sdv_lM(tidx, :) = strain_dv_l ;
        sth_lM(tidx, :) = strain_th_l ;

        % right quarter
        str_rM(tidx, :) = strain_tr_r ;
        sdv_rM(tidx, :) = strain_dv_r ;
        sth_rM(tidx, :) = strain_th_r ;

        % dorsal quarter
        str_dM(tidx, :) = strain_tr_d ;
        sdv_dM(tidx, :) = strain_dv_d ;
        sth_dM(tidx, :) = strain_th_d ;

        % ventral quarter
        str_vM(tidx, :) = strain_tr_v ;
        sdv_vM(tidx, :) = strain_dv_v ;
        sth_vM(tidx, :) = strain_th_v ;
    end
    
    %% Save kymographs
    disp('Saving kymograph data files for Lagrangian pathlines')
    save(apKymoFn, 'tr_apM', 'dv_apM', 'th_apM', ...
        'str_apM', 'sdv_apM', 'sth_apM')
    disp(['Saved kymograph data to: ' apKymoFn])
    save(lKymoFn, 'tr_lM', 'dv_lM', 'th_lM', ...
        'str_lM', 'sdv_lM', 'sth_lM')
    disp(['Saved kymograph data to: ' lKymoFn])
    save(rKymoFn, 'tr_rM', 'dv_rM', 'th_rM', ...
        'str_rM', 'sdv_rM', 'sth_rM')
    disp(['Saved kymograph data to: ' rKymoFn])
    save(dKymoFn, 'tr_dM', 'dv_dM', 'th_dM', ...
        'str_dM', 'sdv_dM', 'sth_dM')
    disp(['Saved kymograph data to: ' dKymoFn])
    save(vKymoFn, 'tr_vM', 'dv_vM', 'th_vM', ...
        'str_vM', 'sdv_vM', 'sth_vM')
    disp(['Saved kymograph data to: ' vKymoFn])
end
disp('done')
