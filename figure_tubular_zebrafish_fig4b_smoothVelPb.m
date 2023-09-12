% script for making Fig 4 b panels


%% ************************************************************************
% *************************************************************************
%               GENERATE SMOOTHED PULLBACK VELOCITY PLOTS
% *************************************************************************
% *************************************************************************

%% Generate Smoothed Velocity Field Plots =================================
close all; clc;

plot_vxyz = false ;      
pivimCoords = tubi.piv.imCoords ;  % coordinate system of the pullback images used in PIV
averagingStyle = 'Lagrangian' ;  % Lagrangian or simple, how velocities are averaged over time
samplingResolution = '1x' ;      % 1x (2x was never implemented)
qscale = 20 ;                    %
vtscale = 5 ;                    % if zero, default is used
vnscale = 2 ;                    % if zero, default is used
gdotscale = 1/5.5 ;               % 
vscale = 0 ;                     % if zero, default is used
invertImage = true ;            % invert the data underneath velocity heatmaps
washout2d = 0. ;                % washout the data image under velocity heatmaps
alphaVal = 1.0 ; % 0.7 ;                % alpha for normal velocity heatmap
qsubsample = 14 ;               % quiver subsampling in pullback space 

% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
else 
    error("Could not parse samplingResolution: set to '1x' or '2x'")
end

% Unpack TubULAR
if strcmp(averagingStyle, 'Lagrangian')
    if doubleResolution
        pivDir = tubi.dir.piv.avg2x ;
    else
        pivDir = tubi.dir.piv.avg ;
    end
elseif strcmp(averagingStyle, 'simple')
    if doubleResolution
        pivDir = tubi.dir.pivSimAvg2x ;
    else
        pivDir = tubi.dir.pivSimAvg ;
    end    
end

% Load simple/Lagrangian average piv results
if strcmp(averagingStyle, 'Lagrangian')
    if doubleResolution
        disp('Loading double Resolution velocity sampling')
        tubi.getVelocityAverage2x()
        velstruct = tubi.velocityAverage2x ;
    else
        disp('Loading single Resolution velocity sampling')
        tubi.getVelocityAverage()
        velstruct = tubi.velocityAverage ;
    end
elseif strcmp(averagingStyle, 'Simple')
    if doubleResolution
        disp('Loading double Resolution simple velocity sampling')
        tubi.getVelocitySimpleAverage2x()
        velstruct = tubi.velocitySimpleAverage2x ;
    else
        disp('Loading single Resolution simple velocity sampling')
        tubi.getVelocitySimpleAverage()
        velstruct = tubi.velocitySimpleAverage ;
    end
end
vsmM = velstruct.v3d ;
v2dsmMum = velstruct.v2dum ;
vnsmM = velstruct.vn ;

% Extract PIV and Smoothed Velocities for a Particular Time Point ---------
timePoints = tubi.xp.fileMeta.timePoints ;

% Choose time point to view
if convert_to_period
    vtscale = T * vtscale;
    qscale = qscale / T;
    vnscale = T * vnscale;
    gdotscale = T * gdotscale;
end

tps = 1:length(timePoints)-1 ;
for i = [1,5,8,12] 
    % i = 1;
    % i = 4;
    % i = 8;
    tp = timePoints(i) ;
    tidx = tubi.xp.tIdx(tp) ;
    disp(['t = ', num2str(tp)])

    % Grab the tangential velocity for this timestep
    vsm_ii = squeeze(vsmM(i, :, :)) ;
    v2dsmum_ii = squeeze(v2dsmMum(i, :, :)) ;
    vnsm_ii = squeeze(vnsmM(i, :, :)) ;

    % Load piv results
    disp('Obtaining raw piv to get field size')
    tubi.getPIV(struct()) ; % NOTE: OPTIONS IS NOT ACTUALLY USED IN THIS FUNCTION
    piv = tubi.piv.raw ;
    % get size of images to make
    gridsz = size(tubi.piv.raw.x{1}) ;
    % Define Nx1 and Mx1 float arrays for xspace and yspace
    xx = tubi.piv.raw.x{tidx}(1, :) ;
    yy = tubi.piv.raw.y{tidx}(:, 1) ;

    % Load metric kinematics (gdot, 2Hvn, divv)
    nV = tubi.nV ;
    mKDir = fullfile(tubi.dir.metricKinematics.root, ...
        strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d'], ...
        tubi.smoothing.lambda, tubi.smoothing.lambda_mesh, tubi.smoothing.lambda_err, tubi.smoothing.nmodes, tubi.smoothing.zwidth), '.', 'p'));
    outdir = fullfile(mKDir, 'measurements') ;
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

    % Prepare for plots
    t0 = tubi.t0set() ;
    tunit = [' ' tubi.timeUnits] ;

    % Load the image to put flow on top
    if strcmp(pivimCoords, 'sp_sme')
        im = imread(sprintf(tubi.fullFileBase.im_sp_sme, tp)) ;
        ylims = [0.25 * size(im, 1), 0.75 * size(im, 1)] ;
    else
        error(['Have not coded for this pivimCoords option. Do so here: ' pivimCoords])
    end
    im = cat(3, im, im, im) ;  % convert to rgb for no cmap change


    % Look at smoothed 2d velocity fields
    vn = reshape(vnsm_ii, gridsz) ;
    vx = reshape(v2dsmum_ii(:, 1), gridsz) ;
    vy = reshape(v2dsmum_ii(:, 2), gridsz) ;


    % % DEBUG:
    % tmp = reshape(veln_filt, [tubi.nU, tubi.nV]) ;
    % subplot(1, 2, 1)
    % imagesc(tmp)
    % subplot(1, 2, 2) 
    % imagesc(vn)

    % Generate Plot -----------------------------------------------------------
    % fig = figure('Visible', 'on',  'units', 'centimeters') ;

    plotType = 'tangent' ; %'normal'; % 'gdot' ; % 'tangent';
    convert_to_period = true;

    if strcmpi(plotType, 'tangent')

        vxb = imgaussfilt(vx, 4) ;
        vyb = imgaussfilt(vy, 4) ;
        if invertImage
            imw = (max(im(:))-im) * washout2d + max(im(:)) * (1-washout2d) ;
        else
            imw = im * washout2d + max(im(:)) * (1-washout2d) ;
        end

        close all
        fig = figure('units', 'normalized', 'visible', 'on') ;

        qopts = struct();
        qopts.visibility = 'on';
        qopts.qsubsample = qsubsample ;
        qopts.overlay_quiver = true ;
        qopts.qscale = qscale ;
        qopts.label_interpreter = 'tex';
        qopts.title_font_weight = 'normal';
        qopts.quiver_line_width = 1.5;
        qopts.fig = fig;
        qopts.cbPosition = [0.8 0.25 0.03 0.2]; %[.9 .3 .02 .3] ;
        qopts.title_font_size = 5;
        qopts.label_font_size = 5;
        % qopts.cbLabelPosition = [2.9 1 0];
        qopts.pbPosition = [0.81, 0.6, 0.1, 0.1] ;

        if convert_to_period
            T = 11;
            qopts.label = ['v_{||} [' char(956) 'm/T]'];
            % qopts.label = ['$v_t$ [' tubi.spaceUnits '/T]'] ;
            % qopts.title = ['Smoothed velocity, $v_t$: $t= $ ', ...
            %     num2str((tp - t0)/T), ' T'] ;
            qopts.title = ['tangent velocity, v_{||}'] ;
            vxb = T * vxb;
            vyb = T * vyb;
        else
            qopts.label = ['$v_t$ [' tubi.spaceUnits '/' tubi.timeUnits ']'] ;
            qopts.title = ['Smoothed velocity, $v_{||}$: $t=$ ', ...
                num2str(tp - t0), tunit] ;
        end

        qopts.ylim = ylims ;
        xyf = struct() ;
        xyf.x = xx ;
        xyf.y = yy ;

        % Plot the coarse-grained tangential velocity as heatmap on top of im
        [h1, h2, h3, gax, cax, pbax] = ...
            vectorFieldHeatPhaseOnImage(imw, xyf, vxb, vyb, vtscale, qopts) ;

        gax = h3.Parent;
        set(gax, 'FontSize', 5);

        % pbKids = pbax.Children;
        for i = 1:numel(pbax.Children)
            try
                pbax.Children(i).FontSize = 5;
            catch
                continue;
            end
        end

        axPos = get(gax, 'Position');
        axRatio = axPos(4) / axPos(3);
        axPos(1) = 0.025;
        axPos(2) = 0.15;
        axPos(3) = 0.7;
        axPos(4) = axRatio * axPos(3);
        set(gax, 'Position', axPos);

        set(gcf, 'Color', [1 1 1]);


        % save data that is plotted
        dataPlotted = struct() ;
        dataPlotted.xx = xx';
        dataPlotted.yy = yy ;
        dataPlotted.vx = vxb ;
        dataPlotted.vy = vyb ;
        dataPlotted.climit = vtscale ;
        dataPlotted.ylims = ylims ;

        % export_fig(sprintf('TubULAR_Paper_Figures/Smoothed_Velocity_Figures/Tangent_Velocity_T%06d.png', timePoints(tp)), '-png', '-r200')

    elseif strcmpi(plotType, 'normal')

        close all
        fig = figure('units', 'normalized', ...
            'outerposition', [0 0 1 1], 'visible', 'on') ;
        labelOpts = struct();
        labelOpts.visibility = 'on';

        if convert_to_period
            T = 11;
            labelOpts.label = ['v_n [' char(956) 'm/T]'];
            % labelOpts.label = ['$v_n$ [' tubi.spaceUnits '/T]'] ;
            % labelOpts.title = ['normal velocity, $v_n$: $t=$ ' ...
            %     num2str((tp - t0)/T), ' T'] ;
            labelOpts.title = 'normal velocity, v_n' ;
            vn = T * vn;

        else
            labelOpts.label = ['$v_n$ [' tubi.spaceUnits '/' tubi.timeUnits ']'] ;
            labelOpts.title = ['normal velocity, $v_n$: $t=$ ' num2str(tp - t0) tunit] ;
        end

        % labelOpts.cmap = twilight_shifted_mod(256) ;
        labelOpts.cmap = brewermap(256, '*RdBu');
        if invertImage
            imw = (max(im(:))-im) * washout2d + max(im(:)) * (1-washout2d) ;
        else
            imw = im * washout2d + max(im(:)) * (1-washout2d) ;
        end

        labelOpts.cbPosition = [0.775 0.25 0.045 0.4];
        labelOpts.cbLabelPosition = [2.9 1 0];

        scalarFieldOnImage(imw, [xx', yy], vn, alphaVal, vnscale, ...
            labelOpts, 'Interpreter', 'tex', ...
            'TitleFontWeight', 'normal' ) ;

        ylim(ylims)

        axPos = get(gca, 'Position');
        axRatio = axPos(4) / axPos(3);
        axPos(1) = 0.025;
        axPos(2) = 0.05;
        axPos(3) = 0.7;
        axPos(4) = axRatio * axPos(3);
        set(gca, 'Position', axPos);


        % save data that is plotted
        dataPlotted = struct() ;
        dataPlotted.xx = xx';
        dataPlotted.yy = yy ;
        dataPlotted.vn = vn ;
        dataPlotted.climit = vnscale ;
        dataPlotted.ylims = ylims ;

        set(gcf, 'Color', [1 1 1]);

        % export_fig(sprintf('TubULAR_Paper_Figures/Smoothed_Velocity_Figures/Normal_Velocity_T%06d.png', timePoints(tp)), '-png', '-r200')

    elseif strcmpi(plotType, 'gdot')

        close all
        fig = figure('units', 'normalized', ...
            'outerposition', [0 0 1 1], 'visible', 'on') ;
        labelOpts = struct();
        labelOpts.visibility = 'on';

        if convert_to_period
            T = 11;
            labelOpts.label = ['$\frac{1}{2}Tr[g^{-1}\dot{g}] \,\,[1/T]$ '];
            % labelOpts.label = ['$v_n$ [' tubi.spaceUnits '/T]'] ;
            % labelOpts.title = ['normal velocity, $v_n$: $t=$ ' ...
            %     num2str((tp - t0)/T), ' T'] ;
            labelOpts.title = 'areal strain rate' ;
            gdot = T * gdot_filt' ;

        else
            labelOpts.label = ['areal strain rate [1/' tubi.timeUnits ']'] ;
            labelOpts.title = ['$\frac{1}{2}Tr[g^{-1}\dot{g}]$: $t=$ ' num2str(tp - t0) tunit] ;
        end

        % labelOpts.cmap = twilight_shifted_mod(256) ;
        labelOpts.cmap = brewermap(256, '*RdBu');
        if invertImage
            imw = (max(im(:))-im) * washout2d + max(im(:)) * (1-washout2d) ;
        else
            imw = im * washout2d + max(im(:)) * (1-washout2d) ;
        end

        labelOpts.cbPosition = [0.775 0.25 0.045 0.4];
        labelOpts.cbLabelPosition = [2.9 1 0];

        % scalarFieldOnImage(imw, [xx', yy], [gdot; gdot], alphaVal, gdotscale, ...
        %     labelOpts, 'Interpreter', 'latex', ...
        %     'TitleFontWeight', 'normal' ) ;
        xx = linspace(0, 100, 100) ;
        yy = linspace(0, 100, 100) ;
        imagesc(xx, yy, gdot) ;
        clim(gdotscale * [-1,1])
        daspect([1,2,1])
        title(labelOpts.title)
        cb = colorbar ;
        colormap(labelOpts.cmap)
        ylabel(cb, labelOpts.label, 'Interpreter', 'latex')
        axis off

        % ylim(ylims)

        axPos = get(gca, 'Position');
        axRatio = axPos(4) / axPos(3);
        axPos(1) = 0.025;
        axPos(2) = 0.05;
        axPos(3) = 0.7;
        axPos(4) = axRatio * axPos(3);
        set(gca, 'Position', axPos);

        set(gcf, 'Color', [1 1 1]);

        % save data that is plotted
        dataPlotted = struct() ;
        dataPlotted.xx = xx ;
        dataPlotted.yy = yy ;
        dataPlotted.gdot = gdot ;
        dataPlotted.gdotscale = gdotscale ;
        dataPlotted.ylims = [0,1] ;
    else

        error('Invalid plot type');

    end

    % Resize Figure for Paper -------------------------------------------------
    set(fig, 'Units', 'centimeters');

    ratio = fig.Position(4) ./ fig.Position(3);
    % fig.Position(3) = 4;
    % fig.Position(4) = ratio * fig.Position(3);

    set(gca, 'FontSize', 5);
    set(gca, 'FontWeight', 'normal');

    set(fig, 'PaperPositionMode', 'auto');
    set(fig.Children, 'FontName', 'Helvetica');
    set(fig.Children, 'FontSize', 5);
    % set(fig.Children, 'Interpreter', 'tex');

    if  strcmpi(plotType, 'gdot')
        fn = sprintf('TubULAR_Paper_Figures/Smoothed_Velocity_Figures/gdot_Velocity_T%06d', timePoints(tp)) ;
        export_fig([fn '.png'], '-png', '-r200')
        F = getframe(gca) ;
        imwrite(F.cdata, [fn '_cdata.png'])
        save([fn '.mat'], 'dataPlotted')
    elseif strcmpi(plotType, 'normal')
        fn = sprintf('TubULAR_Paper_Figures/Smoothed_Velocity_Figures/Normal_Velocity_T%06d', timePoints(tp)) ;
        export_fig([fn '.png'], '-png', '-r200')
        F = getframe(gca) ;
        imwrite(F.cdata, [fn '_cdata.png'])
        save([fn '.mat'], 'dataPlotted')
    elseif strcmpi(plotType, 'tangent')
        fn = sprintf('TubULAR_Paper_Figures/Smoothed_Velocity_Figures_Tangent/Tangent_Velocity_T%06d', timePoints(tp)) ;
        export_fig([fn '.png'], '-png', '-r200')
        F = getframe(gax) ;
        imwrite(F.cdata, [fn '_cdata.png'])
        save([fn '.mat'], 'dataPlotted')

        ofn = [fn '_vx.txt'] ;
        write_txt_with_header(ofn, dataPlotted.vx, 'tangential velocity projected onto pullback coordinates, component in s direction, in microns per beat period')
        
        ofn = [fn '_vy.txt'] ;
        write_txt_with_header(ofn, dataPlotted.vy, 'tangential velocity projected onto pullback coordinates, component in phi direction, in microns per beat period')
    else
        error('what field is this?')
    end
end

close all
disp('done')