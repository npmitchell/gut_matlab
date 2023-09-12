%% Start by running example_zebrafish_heart.m until tubi definition.

cd /mnt/data/tubular_test/zebrafish_heart/analysis



%% TubULAR Analysis Pipeline: Zebrafish Heart (3D + time)
% by NPMitchell and Dillon Cislo

% Add necessary directories to the path
tubularDir = '/mnt/data/code/tubular';
addpath(genpath(tubularDir));

% Add optional external code to the path
% rmpath(genpath('mnt/data/code/gptooolbox'));
addpath(genpath('/mnt/data/code/gptoolbox'));
addpath(genpath('/mnt/data/code/mesh2d'));
addpath(genpath('/home/npmitchell/MATLAB Add-Ons/Apps/PIVlab'));

%% TubULAR Pipeline Initialization ========================================

% We start by clearing the memory and closing all figures
clear; close all; clc;

% The directory containing the zebrafish heart data
dataDir = '/mnt/data/tubular_test/zebrafish_heart/';

% The directory where project files will be generated and saved
[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(projectDir);

% Define TubULAR master settings
overwriteSettings = false;
if (~exist(fullfile(projectDir, 'masterSettings.mat'), 'file') || overwriteSettings)

    stackResolution = [.3524 .3524 2];  % resolution in spaceUnits per pixel
    nChannels = 1;                      % how many channels is the data (ex 2 for GFP + RFP)
    channelsUsed = 1;                   % which channels are used for analysis
    timePoints = 1:30; % 1:50;          % timepoints to include in the analysis
    ssfactor = 4;                       % subsampling factor
    flipy = false ;                     % whether the data is stored inverted relative to real position in lab frame
    timeInterval = 1;                   % physical interval between timepoints
    timeUnits = 'min';                  % physical unit of time between timepoints
    spaceUnits = '$\mu$m';              % physical unit of time between timepoints
    fn = 'Stack_Repeat_014_Time_%03d';	% filename string pattern
    set_preilastikaxisorder = 'xyzc';	% data axis order for subsampled h5 data (ilastik input)
    swapZT = true;                      % whether to swap the z and t dimensions
    
    masterSettings = struct( ...
        'stackResolution', stackResolution, ...
        'nChannels', nChannels, ...
        'channelsUsed', channelsUsed, ...
        'timePoints', timePoints, ...
        'ssfactor', ssfactor, ...
        'flipy', flipy, ...
        'timeInterval', timeInterval, ...
        'timeUnits', timeUnits, ...
        'spaceUnits', spaceUnits, ...
        'fn', fn,...
        'swapZT', swapZT, ...
        'set_preilastikaxisorder', set_preilastikaxisorder, ...
        'nU', 100, ...
        'nV', 100 );
    
    disp(['Saving master settings to ' ...
        fullfile(projectDir, 'masterSettings.mat')]);
    
    save(fullfile(projectDir, 'masterSettings.mat'), 'masterSettings');
    
    clear stackResolution nChannels channelsUsed timePoints
    clear ssfactor flipy timeInterval timeUnits spaceUnits fn
    clear set_preilastikaxisorder swapZT
    
    disp('Saving masterSettings to ./masterSettings.mat')
    
else
    
    disp(['Loading master settings from ' ...
        fullfile(projectDir, 'masterSettings.mat')]);
    
   load(fullfile(projectDir, 'masterSettings.mat'), 'masterSettings');
    
    
end

clear overwriteSettings

%% ************************************************************************
% *************************************************************************
%      PART 1: SURFACE DETECTION USING ImSAnE's 'integralDetector'
% *************************************************************************
% *************************************************************************

% Add ImSAnE to the path
addpath(genpath('/mnt/data/code/imsane_for_git/imsane'));
rmpath('/mnt/data/code/imsane_for_git/imsane/external/CGAL_Code/IsotropicRemeshing');
rmpath('/mnt/data/code/imsane_for_git/imsane/external');

%% Initialize ImSAnE Project ==============================================

% We start by creating the experiment object, which holds the project
% metadata and serves as a front-end for a number of tasks, including data
% loading, detection, fitting, etc.  When instantiating the experiment
% object we must indicate the path of the project directory and the data
% directory, either by passing them as arguments to the constructor or by
% picking them manually from a dialog box

xp = project.Experiment(projectDir, dataDir);

% Set File Metadata -------------------------------------------------------
% First we set the metadata pertaining to the raw data files in the
% structure 'fileMeta'.  ImSAnE assumes that individual time points are
% saved as separate image stacks and that filenames in a time series are
% identical up to a integer specifying the timepoint.
%
% The following file metadata information is required:
%
% * 'directory'         , the project directory (full path)
% * 'dataDir'           , the data directory (full path)
% * 'filenameFormat'    , fprintf type format spec of file name
% * 'timePoints'        , list of times available stored as a vector
% * 'stackResolution'   , stack resolution in microns, e.g. [0.25 0.25 1]
%
% The following file metadata information is optional:
%
% * 'imageSpace'        , bit depth of image, such as uint16 etc., defined
%                         in Stack class
% * 'stackSize'         , size of stack in pixels per dimension 
%                         [xSize ySize zSize]
% * 'swapZT'            , set=1 if time is 3rd dimension and z is 4th
% * 'nChannels'         , The number of channels in the raw data

fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [ masterSettings.fn, '.tif' ];
fileMeta.nChannels          = masterSettings.nChannels;
fileMeta.timePoints         = masterSettings.timePoints;
fileMeta.stackResolution    = masterSettings.stackResolution;
fileMeta.swapZT             = masterSettings.swapZT;

xp.setFileMeta(fileMeta);

% Set Experiment Metadata -------------------------------------------------
% Next we set additional information regarding our experiment as fields in
% the 'expMeta' structure.
%
% The following project metadata information is required:
%
% * 'channelsUsed'      , the channels used, e.g. [1 3] for RGB
% * 'channelColor'      , mapping from element in channels used to RGB=123
% * 'dynamicSurface'    , Boolean, false: static surface
% * 'detectorType'      , name of detector class
% * 'fitterType'        , name of fitter class
%
% The following project meta data information is optional:
%
% * 'description'     , string describing the data set
% * 'jitterCorrection', Boolean, false: No fft based jitter correction 

expMeta                     = struct();
expMeta.channelsUsed        = masterSettings.channelsUsed;
expMeta.channelColor        = 1;
expMeta.description         = 'A beating Zebrafish heart';
expMeta.dynamicSurface      = 1;
expMeta.jitterCorrection    = 0;
expMeta.fitTime             = fileMeta.timePoints(1);
expMeta.detectorType        = 'surfaceDetection.morphsnakesDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';

xp.setExpMeta(expMeta);

% Initialize New Experiment -----------------------------------------------
% Finally we call initNew(), which reads the stack size from the first 
% available time point, then initializes fitter and detector and creates 
% fitOptions and detectOptions based on their defaults.

xp.initNew();

clear fileMeta expMeta

%% Set Surface Detection Objects ==========================================
% We now attempt to detect the surface of interest.  In this pipeline, we
% use the 'integralDetector'.  This detector is essentially a wrapper for
% the morphological snakes external code module.  The input to this module
% is a sub-sampled Ilastik pixel probability map.  Given this map, the
% method attempts to segmemt the image volume into distinct regions by
% minimizing the following (schematic) energy functional
%
% Etot = Est + Ep + Ein + Eout
%
% Est amounts to a surface tension, Ep amounts to a pressure, Ein tries to
% make all voxels inside the segmented region have similar intensities, and
% Eout tries to make all of the voxels outside the segmented region have
% similar intensities. It is assumed that the region of interest is a
% closed volume.
%
% See ImSAnE's 'surfaceDetection.integralDetector' for option documentation

msls_detOpts_fn = fullfile(projectDir, 'msls_detectOpts.mat');

if exist(msls_detOpts_fn, 'file')
    
    load(msls_detOpts_fn, 'detectOptions');
    
else
    
    % The name of the meshlab script that generates a mesh from the point
    % cloud
    mlxprogram = 'laplace_surface_rm_resample30k_reconstruct_LS3_1p2pc_ssfactor4.mlx';
    mlxprogram = fullfile( '/mnt/data/code/meshlab_codes/', mlxprogram );
    
    detectOptions = struct();
    detectOptions.channel = masterSettings.channelsUsed(1);
    detectOptions.ssfactor = masterSettings.ssfactor;
    detectOptions.niter = 200;
    detectOptions.niter0 = 400;
    detectOptions.lambda1 = 1;
    detectOptions.lambda2 = 1;
    detectOptions.pressure = 0.1;
    detectOptions.tension = 2;
    detectOptions.pre_pressure = 0;
    detectOptions.pre_tension = 0;
    detectOptions.post_pressure = 0;
    detectOptions.post_tension = 0;
    detectOptions.exit_thres = 1e-6;
    detectOptions.foreGroundChannel = masterSettings.channelsUsed(1);
    detectOptions.fileName = sprintf( masterSettings.fn, xp.expMeta.fitTime );
    detectOptions.mslsDir = fullfile( projectDir, 'MorphSnakesOutput' );
    detectOptions.ofn_ls = 'msls_DP_';
    detectOptions.ofn_ply = 'mesh_DP_ms_%06d';
    detectOptions.ms_scriptDir = '/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper';
    detectOptions.timepoint = xp.expMeta.fitTime;
    detectOptions.zdim = 3;
    detectOptions.ofn_smoothply = 'mesh_DP_%06d';
    detectOptions.mlxprogram = fullfile( projectDir, mlxprogram );
    detectOptions.init_ls_fn = 'empty_string';
    detectOptions.run_full_dataset = false;
    detectOptions.radius_guess = 200;
    detectOptions.dset_name = 'exported_data';
    detectOptions.center_guess = '30,125,125';
    detectOptions.save = true;
    detectOptions.plot_mesh3d = false;
    detectOptions.dtype = 'h5';
    detectOptions.mask = 'none';
    detectOptions.mesh_from_pointcloud = true;
    detectOptions.prob_searchstr = '_Probabilities.h5';
    detectOptions.preilastikaxisorder = masterSettings.set_preilastikaxisorder;
    detectOptions.ilastikaxisorder = 'cxyz'; %'cxyz';
    detectOptions.physicalaxisorder = 'yxzc';
    detectOptions.include_boundary_faces = true;
    detectOptions.smooth_with_matlab = 0.01;
    detectOptions.pythonVersion = '';
    
end

xp.setDetectOptions( detectOptions );

%% ************************************************************************
% *************************************************************************
%               PART 2: TubULAR -- SURFACE PARAMETERIZATION
% *************************************************************************
% *************************************************************************

%% Initialize the TubULAR Object ==========================================
close all; clc;

tubOpts = struct();
tubOpts.meshDir = fullfile(projectDir, 'TubULAR_Results');	% Director where meshes produced in previous part reside. All TubULAR results will be stored relative to this directory
tubOpts.flipy = masterSettings.flipy;               % Set to true if data volume axes are inverted in chirality wrt physical lab coordinates               
tubOpts.timeInterval = masterSettings.timeInterval; % Spacing between adjacent timepoints in units of timeUnits
tubOpts.timeUnits = masterSettings.timeUnits;       % Units of time, so that adjacent timepoints are timeUnits * timeInterval apart
tubOpts.spaceUnits = masterSettings.spaceUnits;     % Units of space in LaTeX, for ex '$mu$m' for micron
tubOpts.nU = masterSettings.nU;         % How many points along the longitudinal axis to sample surface
tubOpts.nV = masterSettings.nV;         % How many points along the circumferential axis to sample surface
tubOpts.t0 = xp.fileMeta.timePoints(1);	% Reference timepoint used to define surface-Lagrangian and Lagrangian measurements
tubOpts.normalShift = 2;        % Additional dilation acting on surface for texture mapping
tubOpts.a_fixed = 2.0;          % Fixed aspect ratio of pullback images. Setting to 1.0 is most conformal mapping option.
tubOpts.adjustlow = 1.00;       % Floor for intensity adjustment
tubOpts.adjusthigh = 99.9;      % ceil for intensity adjustment (clip)
tubOpts.phiMethod = 'curved3d';	% Method for following surface in surface-Lagrangian mapping [(s,phi) coordinates]
tubOpts.lambda_mesh = 0.0;        % Smoothing applied to the mesh before DEC measurements
tubOpts.lambda = 0.02;             % Smoothing applied to computed values on the surface
tubOpts.lambda_err = 0;         % Additional smoothing parameter, optional

disp('defining TubULAR class instance (tubi= tubular instance)')
tubi = TubULAR(xp, tubOpts) ;
disp('done defining TubULAR instance')

clear tubOpts



%%

tubi.smoothing.zwidth = 3 ;
tubi.smoothing.nmodes = 5 ;

%% Load pullback pathline metric Kinematics
t0Pathline = 1; 
sresStr = '' ;
nU = tubi.nU ;
outdir = './TubULAR_Paper_Figures/Kymographs_Pathlines_lambda0p020_nmodes05w03/' ;
mkdir(outdir)

tubi.smoothing.nmodes = 5 ; tubi.smoothing.zwidth = 3 ; 
ntp = length(tubi.xp.fileMeta.timePoints) ;
H2vn_apM = zeros(ntp-1, nU) ;
divv_apM = zeros(ntp-1, nU) ;
radius_apM = zeros(ntp-1, nU) ;
veln_apM = zeros(ntp-1, nU) ;
HH_apM = zeros(ntp-1, nU) ;
gdot_apM = zeros(ntp-1, nU) ;
for tidx = 1:ntp-1
    disp(['loading tidx = '  num2str(tidx)])
    tp = tubi.xp.fileMeta.timePoints(tidx) ;
    mKDir = fullfile(tubi.dir.metricKinematics.root, ...
    strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d'], ...
    tubi.smoothing.lambda, tubi.smoothing.lambda_mesh, ...
    tubi.smoothing.lambda_err, tubi.smoothing.nmodes, ...
    tubi.smoothing.zwidth), '.', 'p'));
    mKPDir = fullfile(mKDir, sprintf('pathline_%04dt0', t0Pathline)) ;
    tmpdir = fullfile(mKPDir, 'measurements') ;
    Hfn = fullfile(tmpdir, sprintf('HH_pathline%04d_%03d.mat', t0Pathline, tp))   ;
    efn = fullfile(tmpdir, sprintf('gdot_pathline%04d_%03d.mat', t0Pathline, tp)) ;
    dfn = fullfile(tmpdir, sprintf('divv_pathline%04d_%03d.mat', t0Pathline, tp)) ;
    nfn = fullfile(tmpdir, sprintf('veln_pathline%04d_%03d.mat', t0Pathline, tp)) ;
    rfn = fullfile(tmpdir, sprintf('radius_pathline%04d_%03d.mat', t0Pathline, tp)) ;
    H2vnfn = fullfile(tmpdir, sprintf('H2vn_pathline%04d_%03d.mat', t0Pathline, tp)) ;

    tmp1 = load(H2vnfn) ;
    H2vn_apM(tidx, :) = tmp1.H2vn_ap ;
    tmp2 = load(dfn) ;
    divv_apM(tidx, :) = tmp2.divv_ap ;
    tmp = load(rfn) ;
    radius_apM(tidx, :) = tmp.radius_ap ;
    tmp = load(nfn) ;
    veln_apM(tidx, :) = tmp.veln_ap ;
    tmp = load(Hfn) ;
    HH_apM(tidx, :) = tmp.HH_ap ;
    tmp = load(efn) ;
    gdot_apM(tidx, :) = tmp.gdot_ap ;
end


%% Some plotting params

climit = 0.075*1.1; % 0.2 ;
climit_err = 0.75*1.1; % 0.2 ;
climit_veln = climit * 10 ;
climit_H = climit * 2 ;
climit_radius = 45;

%% Generate Kymograph Plot ================================================
close all; clc;
% figure; 

% Unit definitions for axis labels
convert_to_period = true;
if convert_to_period
    
    % unitstr = '[1/T]' ;
    % Hunitstr = [ '[1/' tubi.spaceUnits ']' ];
    % vunitstr = [ '[' tubi.spaceUnits '/T]' ];
    
    unitstr = '[1/T]' ;
    Hunitstr = [ '[1/' char(956) 'm]' ];
    vunitstr = [ '[' char(956) 'm/T]' ];
    
else
    
    unitstr = [ '[1/' tubi.timeUnits ']' ];
    Hunitstr = [ '[1/' tubi.spaceUnits ']' ];
    vunitstr = [ '[' tubi.spaceUnits '/' tubi.timeUnits ']' ];
    
end

titleadd = ': circumferentially averaged';

plotTypes = {'gdot' , 'divv', 'H2vn', 'veln', 'HH', 'radius'} ;

for ii = length(plotTypes)
    plotType= plotTypes{ii} ;
    if strcmpi(plotType, 'gdot')
        m2plot = gdot_apM;
        titles = '$\frac{1}{2}\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$';
        labels = ['$\frac{1}{2}\textrm{Tr}[g^{-1}\dot{g}]$ ' unitstr];
        climits = climit;
    elseif strcmpi(plotType, 'HH')
        m2plot = HH_apM;
        % titles = 'mean curvature, $H$';
        % labels = ['mean curvature, $H$ ' Hunitstr];
        titles = 'mean curvature, H';
        labels = ['mean curvature, H ' Hunitstr];
        climits = climit_H;
    elseif strcmpi(plotType, 'divv')
        m2plot = divv_apM;
        % titles = 'divergence of flow, $\nabla \cdot \mathbf{v}$';
        % labels = ['$\nabla \cdot \mathbf{v}$ ' unitstr];
        titles = 'flow divergence';
        labels = ['\nabla \cdot v_{||} ' unitstr];
        climits = climit;
    elseif strcmpi(plotType, 'veln')
        m2plot = veln_apM;
        % titles = 'normal velocity, $v_n$';
        % labels = ['normal velocity, $v_n$ ' vunitstr];
        titles = 'normal velocity, v_n';
        labels = ['normal velocity, v_n ' vunitstr];
        climits = climit_veln;
    elseif strcmpi(plotType, 'H2vn')
        m2plot = H2vn_apM;
        % titles = 'normal motion, $v_n 2 H$';
        % labels = ['normal motion, $v_n 2 H $ ' unitstr];
        titles = 'normal motion';
        labels = ['2 H v_n' unitstr];
        climits = climit;
    elseif strcmpi(plotType, 'radius')
        m2plot = radius_apM;
        % titles = 'radius';
        % labels = ['radius [' tubi.spaceUnits ']'];
        titles = 'radius';
        labels = ['radius [' char(956) ']'];
        climits = climit_radius;
    else
        error('Invalid plot type');
    end

    % We relate the normal velocities to the divergence / 2 * H.
    tps = tubi.xp.fileMeta.timePoints(1:end-1) - tubi.t0;

    fig = figure('Visible', 'on',  'units', 'centimeters') ;

    if convert_to_period

        T = 11; % Length of period
        tps = tps ./ T;

        if ismember(lower(plotType), ...
                {'gdot', 'divv', 'veln', 'h2vn'})
            m2plot = T * m2plot;
            climits = T * climits;
        end

        imagesc((1:nU)/nU, tps, m2plot)

        if strcmpi(plotType, 'radius')
            if climits > 0
                % caxis([min(radius_apM(:)), max(radius_apM(:))]);
                caxis([20 45]);
            end
            cmap = brewermap(512, 'RdBu');
            colormap(cmap((257:end).', :));

        else
            if climits > 0
                caxis([-climits, climits])
            end
            colormap(brewermap(256, '*RdBu'));
        end

        axis square

        xticks(0:0.2:1);
        yticks(0:0.5:2.5);

        % title and save
        % title([titles, titleadd], 'Interpreter', 'Latex')
        % ylabel(['time [t/T]'], 'Interpreter', 'Latex')
        % xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
        % cb = colorbar() ;
        % ylabel(cb, labels, 'Interpreter', 'Latex')
        % set(gcf, 'Color', [1 1 1]);

        % title([titles, titleadd], 'FontWeight', 'normal')
        if contains(titles, '$')
            title(titles, 'FontWeight', 'normal', 'interpreter', 'latex');
        else
            title(['\langle' titles char(9002)], 'FontWeight', 'normal');
        end
        ylabel('time [T]');
        xlabel('position [s/L]');

        cb = colorbar() ;
        if contains(titles, '$')
            cb.Label.Interpreter = 'latex';
        end
        cb.Label.String = labels;
        cb.Position(1) = 0.75;
        cb.Position(2) = 0.25;
        cb.Position(3) = 0.045;
        cb.Position(4) = 0.6;
        cb.Label.Position(1) = 3.5;

        axPos = get(gca, 'Position');
        axRatio = axPos(4) / axPos(3);
        axPos(1) = 0.10;
        axPos(2) = 0.20;
        axPos(3) = 0.675;
        axPos(4) = axRatio * axPos(3);
        set(gca, 'Position', axPos);



    else

        imagesc((1:nU)/nU, tps, m2plot)
        if climits > 0
            caxis([-climits, climits])
        end
        colormap(brewermap(256, '*RdBu'));
        axis equal

        % title and save
        title([titles, titleadd], 'Interpreter', 'Latex')
        ylabel(['time [' tubi.timeUnits ']'], 'Interpreter', 'Latex')
        xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
        cb = colorbar() ;
        ylabel(cb, labels, 'Interpreter', 'Latex')

    end

    set(gcf, 'Color', [1 1 1]);

    % Resize Figure for Paper -------------------------------------------------
    set(fig, 'Units', 'centimeters');

    set(fig, 'Position', [0,0,4,4]) ;
    % ratio = fig.Position(4) ./ fig.Position(3);
    % fig.Position(3) = 3.6;
    % fig.Position(4) = ratio * fig.Position(3);

    set(gca, 'FontSize', 10);
    set(gca, 'FontWeight', 'normal');

    set(fig, 'PaperPositionMode', 'auto');
    set(fig.Children, 'FontName', 'Helvetica');
    
    % save data
    soverL = (1:nU)/nU ; 
    save(fullfile(outdir, [ plotType '.mat']), 'tps', 'm2plot', 'soverL')

    saveas(gcf, fullfile(outdir, [ plotType '.pdf']))

end



%% Generate Cross-Correlation Plot ========================================
close all; clc;

T = 11; % Period size

% A = A(4:4+2*T, :); B = B(4:4+2*T, :);

% Remove some timepoints so that we have 2 periods total (22 tps)
% Also remove two pixels at each vertical (spatial) edge of the kymo 
% A = A(4:4+2*T, 3:end-2); B = B(4:4+2*T, 3:end-2);
% 

% Note that we remove the edge effects with the following reasoning:
% We diffuse across 0.02*1000 steps = 20 timesteps.
% The  length scale of diffusion  is sqrt((#dims)*2*D*t) , with D=1
% and t = 20 timesteps. We live in 2D meshes, so 
% 2*sqrt(D*t) = 2*sqrt(20). 
% So we have 2*sqrt(20) = 8.94427191 um.
% The organ length is 223 um as measured by max(mesh.s). 
% So this is 4.01088 % on each side.
% Add one unit for the zwidth == 3.
% 4% of 100 s steps is 4, add 1 gives 5 columns.


nEdges2consider = 1:15 ;
maxLocs = zeros(size(nEdges2consider)) ;
for nEdge = nEdges2consider
    
    A = H2vn_apM;   
    B = divv_apM;
    
    A = T * A; B = T* B;

    A = A(4:4+2*T, nEdge+1:end-nEdge); B = B(4:4+2*T, nEdge+1:end-nEdge);
    
    % %% Preview the circshift action here
    % for reapeat = 1:5
    %     for i = 0:(size(B,1)-1)
    %         subplot(1,2,1)
    %         imagesc(A); caxis([-climit, climit])
    %         subplot(1,2,2)
    %         imagesc(circshift(conj(B), [i 0]));
    %         caxis([-climit, climit])
    %         sgtitle(['shift ' num2str(i)])
    %         pause(0.5)    
    %     end
    % end
    % 
    % %% 2D cross correlation
    % C2d = zeros(size(B));
    % for i = 0:(size(B,1)-1)
    %     for j = 0:(size(B,2)-1)
    % 
    %         C2d(i+1,j+1) = sum(A .* circshift(conj(B), [i j]), 'all');
    % 
    %     end
    % end
    % imagesc(C2d); colorbar; colormap bwr
    
    % 1D cross correlation
    C = zeros(size(B,1), 1);
    shifts = (-size(B,1)+1)*0.5:(size(B,1)-1)*0.5 ;
    for idx = 1:length(shifts)
        i = shifts(idx) ;
        C(idx) = sum(A .* circshift(conj(B), [i 0]), 'all');
    end
    [~, t0_idx] = min(abs(shifts)) ;
    timeSteps = shifts / T ;
    plot(timeSteps, C)
    
    % NEW VERSION for plotting
    
    % Fit nonlinear model to correlation function
    % fitFunc = @(p, t) p(1) * cos(2 * pi .* t - p(2)) ;
    % param0 = [ max(C) 0.5 ];
    % options = optimoptions( 'lsqcurvefit', ...
    %     'Algorithm', 'levenberg-marquardt', ...
    %     'MaxIterations', 1000, ...
    %     'MaxFunctionEvaluations', 10000, ...
    %     'Display', 'none');
    % [param, resnorm, ~, exitflag, output] = ...
    %     lsqcurvefit(fitFunc, param0, timeSteps, C', [], [], options);
    
    fitOptions = fitoptions('Method', 'NonlinearLeastSquares', ...
        'Lower', [-Inf,2*pi-eps,-Inf], 'Upper', [Inf, 2*pi+eps, Inf], ...
        'StartPoint', [400, 2*pi, 0]) ;
    % modelFit = fit(timeSteps', C,'sin1', fitOptions) ;
    
    
    fitOptions = fitoptions('Method', 'NonlinearLeastSquares', ...
        'Lower', [-Inf,6,-Inf], 'Upper', [Inf, 6.5, Inf], ...
        'StartPoint', [400, 2*pi, 0]) ;
    modelFit = fit([timeSteps'], [C],'sin1', fitOptions) ;
    
    % Generate Visualization --------------------------------------------------
    
    plot_lw = 0.75;
    scatter_size = 7.5;
    
    fig = figure('Visible', 'on',  'units', 'centimeters') ;
    
    % Assumes A and B are the same size
    scatter(timeSteps, C, scatter_size, 'o', ...
        'MarkerEdgeColor', tubi.plotting.colors(1,:), ...
        'LineWidth', plot_lw );
    
    hold on
    
    tplot = linspace(timeSteps(1), timeSteps(end), 300);
    %plot(tplot, fitFunc(param, tplot), '-', 'LineWidth', plot_lw, ...
    %    'Color', tubi.plotting.colors(1,:) );
    plot(tplot, modelFit.a1 * sin(modelFit.b1 * tplot + modelFit.c1), ...
        '-', 'LineWidth', plot_lw, 'Color', tubi.plotting.colors(1,:) )
    
    % Zero point of the sin is = -modelFit.c1, so the MAX is 1/4 * T later.
    % Multiply modelFit.c1 by 1/2pi to convert to period, then add 1/4 T, where
    % T = 1 here because we put timeSteps in units of T. 
    maxLoc = -modelFit.c1 / (2*pi) + 0.25 ;
    ci = confint(modelFit) ;
    ci = ci(:, 3) ;
    maxLoc_ci = -ci / (2*pi) + 0.25 ;
    assert(all(abs(maxLoc - maxLoc_ci) - abs(maxLoc - maxLoc_ci(1))<10^-15))
    maxLoc_unc = (maxLoc - maxLoc_ci) * 0.5 ;
    maxLoc_unc = maxLoc_unc(maxLoc_unc > 0) ;
    scatter(maxLoc, modelFit(maxLoc), 1.5 * scatter_size, 'filled', 'k');
    ylims = ylim ;
    plot([maxLoc, maxLoc], [ylims(1) modelFit(maxLoc)], '--k', ...
        'LineWidth', max(0.75 * plot_lw, 0.5))
    
    % hold off
    
    % Restrict attention to near zero.
    xlim([-1,1])
    xlims = xlim ;
    ylim(ylims);
    
    % xlabel(['time [' tubi.timeUnits ']'], 'Interpreter', 'Latex');
    xlabel(['time shift \Delta [T]']);
    ylabel({'cross correlation', 'C [2Hv_n(t), \nabla\cdotv_{||}(t+\Delta)]'});
    
    try
        tanno = annotation('textarrow', ...
           (maxLoc-xlims(1))/diff(xlims) * [1.15 1.1], ...
           ( modelFit(maxLoc)-ylims(1))/diff(ylims) * [1.025 0.975], ...
           'String', sprintf('\\Delta = %0.3f\\pm%0.3f T', maxLoc, maxLoc_unc), ...
           'LineWidth', plot_lw, ...
           'FontSize', 10, 'HeadLength', 2.5, 'HeadWidth', 2.5);
    catch
        print('annotation failed')
    end

    set(gcf, 'Color', [1 1 1]);
    
    % Resize Figure for Paper -------------------------------------------------
    set(fig, 'Units', 'centimeters');
    
    % ratio = fig.Position(4) ./ fig.Position(3);
    % fig.Position(3) = 3.6;
    % fig.Position(4) = ratio * fig.Position(3);
    set(fig, 'Position', [0, 0, 7, 4])
    
    set(gca, 'FontSize', 10);
    set(gca, 'FontWeight', 'normal');
    
    set(fig, 'PaperPositionMode', 'auto');
    set(fig.Children, 'FontName', 'Helvetica');
    
    maxLocs(nEdge) = maxLoc ;
    maxLoc_uncs(nEdge) = maxLoc_unc ;
    close all
end

%% Look at how result varies with edge truncation
errorbar(maxLocs, maxLoc_uncs, '.')
xlabel('# columns truncated near each edge')
ylabel('phase shift [1/T]')

saveas(gcf, fullfile('./TubULAR_Paper_Figures', 'metricKinematics_phaseRelation_offset_vs_nEdgeMask.pdf'))
fn = fullfile('./TubULAR_Paper_Figures', 'metricKinematics_phaseRelation_offset_vs_nEdgeMask.mat') ;
save(fn, 'maxLocs', 'maxLoc_uncs')

%% Repeat for chosen nEdge

nEdge = 5 ;
A = H2vn_apM;   
B = divv_apM;

A = T * A; B = T* B;

A = A(4:4+2*T, nEdge+1:end-nEdge); B = B(4:4+2*T, nEdge+1:end-nEdge);

% %% Preview the circshift action here
% for reapeat = 1:5
%     for i = 0:(size(B,1)-1)
%         subplot(1,2,1)
%         imagesc(A); caxis([-climit, climit])
%         subplot(1,2,2)
%         imagesc(circshift(conj(B), [i 0]));
%         caxis([-climit, climit])
%         sgtitle(['shift ' num2str(i)])
%         pause(0.5)    
%     end
% end
% 
% %% 2D cross correlation
% C2d = zeros(size(B));
% for i = 0:(size(B,1)-1)
%     for j = 0:(size(B,2)-1)
% 
%         C2d(i+1,j+1) = sum(A .* circshift(conj(B), [i j]), 'all');
% 
%     end
% end
% imagesc(C2d); colorbar; colormap bwr

% 1D cross correlation
C = zeros(size(B,1), 1);
shifts = (-size(B,1)+1)*0.5:(size(B,1)-1)*0.5 ;
for idx = 1:length(shifts)
    i = shifts(idx) ;
    C(idx) = sum(A .* circshift(conj(B), [i 0]), 'all');
end
[~, t0_idx] = min(abs(shifts)) ;
timeSteps = shifts / T ;
plot(timeSteps, C)

% NEW VERSION for plotting

% Fit nonlinear model to correlation function
% fitFunc = @(p, t) p(1) * cos(2 * pi .* t - p(2)) ;
% param0 = [ max(C) 0.5 ];
% options = optimoptions( 'lsqcurvefit', ...
%     'Algorithm', 'levenberg-marquardt', ...
%     'MaxIterations', 1000, ...
%     'MaxFunctionEvaluations', 10000, ...
%     'Display', 'none');
% [param, resnorm, ~, exitflag, output] = ...
%     lsqcurvefit(fitFunc, param0, timeSteps, C', [], [], options);

fitOptions = fitoptions('Method', 'NonlinearLeastSquares', ...
    'Lower', [-Inf,2*pi-eps,-Inf], 'Upper', [Inf, 2*pi+eps, Inf], ...
    'StartPoint', [400, 2*pi, 0]) ;
% modelFit = fit(timeSteps', C,'sin1', fitOptions) ;


fitOptions = fitoptions('Method', 'NonlinearLeastSquares', ...
    'Lower', [-Inf,6,-Inf], 'Upper', [Inf, 6.5, Inf], ...
    'StartPoint', [400, 2*pi, 0]) ;
modelFit = fit([timeSteps'], [C],'sin1', fitOptions) ;

% Generate Visualization --------------------------------------------------

plot_lw = 0.75;
scatter_size = 7.5;

fig = figure('Visible', 'on',  'units', 'centimeters') ;

% Assumes A and B are the same size
scatter(timeSteps, C, scatter_size, 'o', ...
    'MarkerEdgeColor', tubi.plotting.colors(1,:), ...
    'LineWidth', plot_lw );

hold on

tplot = linspace(timeSteps(1), timeSteps(end), 300);
%plot(tplot, fitFunc(param, tplot), '-', 'LineWidth', plot_lw, ...
%    'Color', tubi.plotting.colors(1,:) );
plot(tplot, modelFit.a1 * sin(modelFit.b1 * tplot + modelFit.c1), ...
    '-', 'LineWidth', plot_lw, 'Color', tubi.plotting.colors(1,:) )

% Zero point of the sin is = -modelFit.c1, so the MAX is 1/4 * T later.
% Multiply modelFit.c1 by 1/2pi to convert to period, then add 1/4 T, where
% T = 1 here because we put timeSteps in units of T. 
maxLoc = -modelFit.c1 / (2*pi) + 0.25 ;
ci = confint(modelFit) ;
ci = ci(:, 3) ;
maxLoc_ci = -ci / (2*pi) + 0.25 ;
assert(all(abs(maxLoc - maxLoc_ci) - abs(maxLoc - maxLoc_ci(1))<10^-15))
maxLoc_unc = (maxLoc - maxLoc_ci) * 0.5 ;
maxLoc_unc = maxLoc_unc(maxLoc_unc > 0) ;
scatter(maxLoc, modelFit(maxLoc), 1.5 * scatter_size, 'filled', 'k');
ylims = ylim ;
plot([maxLoc, maxLoc], [ylims(1) modelFit(maxLoc)], '--k', ...
    'LineWidth', max(0.75 * plot_lw, 0.5))

% hold off

% Restrict attention to near zero.
xlim([-1,1])
xlims = xlim ;
ylim(ylims);

% xlabel(['time [' tubi.timeUnits ']'], 'Interpreter', 'Latex');
xlabel(['time shift \Delta [T]']);
ylabel({'cross correlation', 'C [2Hv_n(t), \nabla\cdotv_{||}(t+\Delta)]'});

try
    tanno = annotation('textarrow', ...
       (maxLoc-xlims(1))/diff(xlims) * [1.15 1.1], ...
       ( modelFit(maxLoc)-ylims(1))/diff(ylims) * [1.025 0.975], ...
       'String', sprintf('\\Delta = %0.3f\\pm%0.3f T', maxLoc, maxLoc_unc), ...
       'LineWidth', plot_lw, ...
       'FontSize', 10, 'HeadLength', 2.5, 'HeadWidth', 2.5);
catch
    print('annotation failed')
end

set(gcf, 'Color', [1 1 1]);

% Resize Figure for Paper -------------------------------------------------
set(fig, 'Units', 'centimeters');

% ratio = fig.Position(4) ./ fig.Position(3);
% fig.Position(3) = 3.6;
% fig.Position(4) = ratio * fig.Position(3);
set(fig, 'Position', [0, 0, 7, 4])

set(gca, 'FontSize', 10);
set(gca, 'FontWeight', 'normal');

set(fig, 'PaperPositionMode', 'auto');
set(fig.Children, 'FontName', 'Helvetica');
axis square

%% Save data
fitcurve_x = tplot' ;
fitcurve_y = modelFit.a1 * sin(modelFit.b1 * tplot + modelFit.c1) ;
fitcurve_y = fitcurve_y' ;
x = timeSteps' ;
y = C ;
save(fullfile('./TubULAR_Paper_Figures', 'metricKinematics_phaseRelation.mat'), ...
    'x', 'y', 'fitcurve_x', 'fitcurve_y')
m2save = [x, y] ; % , fitcurve_x, fitcurve_y] ;
fn = fullfile('./TubULAR_Paper_Figures', 'metricKinematics_phaseRelation.txt') ;
header = 'time shift (in units of T), cross correlation' ;
write_txt_with_header(fn, m2save, header)

% save figure
saveas(gcf, fullfile('./TubULAR_Paper_Figures', 'metricKinematics_phaseRelation_20230905.pdf'))


%% OLD VERSION

% % tileT = T*2;
% B = B(4:end, :) ;
% tileA1 = A(4:4+T-1,:) ;
% tileA2 = A(4+T:4+2*T-1, :) ;
% AA1 = cat(1, tileA1, tileA1, tileA1, tileA1, tileA1, tileA1, tileA1, tileA1) ;
% AA2 = cat(1, tileA2, tileA2, tileA2, tileA2, tileA2, tileA2, tileA2, tileA2) ;
% 
% corrAB2 = xcorr2(AA1, B);
% % Take the column that only does cross correlation in TIME, not space
% corrAB_1 = corrAB2(:, size(A,2));
% corrAB2 = xcorr2(AA1, B);
% % Take the column that only does cross correlation in TIME, not space
% corrAB_2 = corrAB2(:, size(A,2));
% 
% % Average the two periods together 
% corrAB = (corrAB_1 + corrAB_2)*0.5 ;
%
% timeSteps = -(size(A,1)-1):(size(A,1)-1);
% timeSteps = timeSteps / T;
%
% Find locations of maximum signal
% [maxCorr, maxLoc] = max(corrAB);
% 
% xLim = [timeSteps(1) timeSteps(end)];
% yLim = [-5, 5];
% yLim = 550 * [-1 1];

% 
% %%
% % Fit nonlinear model to correlation function
% % fitFunc = @(p, t) p(1) * cos(p(2) .* (t-p(4)) - p(3)) .* exp(-(t-p(4)).^2 / p(5).^2);
% % param0 = [ min(corrAB) 2*pi 0 0 1];
% fitFunc = @(p, t) p(1) * cos(2 * pi .* (t-p(2)) - p(3)) .* exp(-(t-p(2)).^2 / p(4).^2);
% fitFunc = @(p, t) p(1) * cos(2 * pi .* (t-p(2)) - p(3)) .* exp(-(t-p(2)).^2 / p(4).^2);
% param0 = [ min(corrAB) 0 0 1];
% options = optimoptions( 'lsqcurvefit', ...
%     'Algorithm', 'levenberg-marquardt', ...
%     'MaxIterations', 1000, ...
%     'MaxFunctionEvaluations', 10000, ...
%     'Display', 'none');
% [param, resnorm, ~, exitflag, output] = ...
%     lsqcurvefit(fitFunc, param0, timeSteps, corrAB.', [], [], options);
% 
% % Generate Visualization --------------------------------------------------
% 
% plot_lw = 0.75;
% scatter_size = 7.5;
% 
% fig = figure('Visible', 'on',  'units', 'centimeters') ;
% 
% % Assumes A and B are the same size
% scatter(timeSteps, corrAB, scatter_size, 'o', ...
%     'MarkerEdgeColor', tubi.plotting.colors(1,:), ...
%     'LineWidth', plot_lw );
% 
% hold on
% 
% tplot = linspace(timeSteps(1), timeSteps(end), 300);
% plot(tplot, fitFunc(param, tplot), '-', 'LineWidth', plot_lw, ...
%     'Color', tubi.plotting.colors(1,:) );
% 
% maxLoc = (2*pi*param(2) + param(3))  +1;
% scatter(maxLoc, fitFunc(param, maxLoc), 1.5 * scatter_size, 'filled', 'k');
% plot([maxLoc, maxLoc], [yLim(1) fitFunc(param, maxLoc)], '--k', ...
%     'LineWidth', max(0.75 * plot_lw, 0.5))
% 
% hold off
% 
% xlim(xLim);
% ylim(yLim);
% 
% % xlabel(['time [' tubi.timeUnits ']'], 'Interpreter', 'Latex');
% xlabel(['time shift ' char(916) ' [T]']);
% ylabel({'cross correlation', 'C [2Hv_n, \nabla\cdotv_{||}]'});
% 
% 
% tanno = annotation('textarrow', ...
%     (maxLoc-xLim(1))/diff(xLim) * [1.15 1.1], ...
%     ( fitFunc(param, maxLoc)-yLim(1))/diff(yLim) * [1.025 0.975], ...
%     'String', sprintf('\\Delta = %0.2f T', maxLoc), ...
%     'LineWidth', plot_lw, ...
%     'FontSize', 5, 'HeadLength', 2.5, 'HeadWidth', 2.5);
% 
% set(gcf, 'Color', [1 1 1]);
% 
% % Resize Figure for Paper -------------------------------------------------
% set(fig, 'Units', 'centimeters');
% 
% ratio = fig.Position(4) ./ fig.Position(3);
% fig.Position(3) = 3.6;
% fig.Position(4) = ratio * fig.Position(3);
% 
% set(gca, 'FontSize', 5);
% set(gca, 'FontWeight', 'normal');
% 
% set(fig, 'PaperPositionMode', 'auto');
% set(fig.Children, 'FontName', 'Helvetica');
% 
% 
% 
% %%
% 
% 
% plot(timeSteps, corrAB, 'LineWidth', 3);
% hold on
% scatter(timeSteps(maxLoc), maxCorr, 100, 'filled', 'k');
% plot([timeSteps(maxLoc), timeSteps(maxLoc)], [yLim(1) maxCorr], '--k', ...
%     'LineWidth', 3)
% hold off
% 
% 
% 
% 
% % tanno = annotation('textarrow', ...
% %     (timeSteps(maxLoc)-xLim(1))/diff(xLim) * [1.05 1], ...
% %     (maxCorr-yLim(1))/diff(yLim) * [1.025 0.975], ...
% %     'String', sprintf('$$\\Delta \\frac{t}{T} = %0.2f$$', timeSteps(maxLoc)), ...
% %     'LineWidth', 2, 'Interpreter', 'Latex' );
% tanno = annotation('textarrow', ...
%     (timeSteps(maxLoc)-xLim(1))/diff(xLim) * [1.075 1], ...
%     (maxCorr-yLim(1))/diff(yLim) * [1.00 0.975], ...
%     'String', sprintf('$$\\Delta = %0.2f$$ T', timeSteps(maxLoc)), ...
%     'LineWidth', 2, 'Interpreter', 'Latex' );
% 
% 
% 
% %% 
% % get rid of the first and last few timepoints and endcaps
% ecap = 3 ;
% H2vnM = H2vn_apM(3:end-3, ecap:end-ecap) ;
% divvM = divv_apM(3:end-3, ecap:end-ecap) ;
% binscatter(H2vn_apM(:), divv_apM(:), 100)

