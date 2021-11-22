%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quick instantiation of a QuapSlap object 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace ========================================================
% We start by clearing the memory and closing all figures
clear; close all; clc;
% change this path, for convenience
% cd /mnt/crunch/48Ygal4-UAShistRFP/201904031830_great/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2/data/
% cd /mnt/crunch/48YGal4UasLifeActRuby/201904021800_great/Time6views_60sec_1p4um_25x_1p0mW_exp0p150_3/data/
% cd /mnt/data/48YGal4UasLifeActRuby/201902201200_unusualfolds/Time6views_60sec_1p4um_25x_obis1_exp0p35_3/data/
%cd /mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data
 cd /home/yuzhenglin/membraneproject/caaxdataset
% .=========.
% |  VIP10  |
% .=========.
% cd /mnt/crunch/gut/48YGal4UasLifeActRuby/201907311600_48YGal4UasLifeActRuby_60s_exp0p150_1p0mW_25x_1p4um
% cd /mnt/crunch/gut/48YGal4klarUASCAAXmChHiFP/202001221000_60sec_1p4um_25x_1mW_2mW_exp0p25_exp0p7/Time3views_1017/data/
% cd /mnt/crunch/gut/Mef2Gal4klarUASCAAXmChHiFP/202003151700_1p4um_0p5ms3msexp/Time3views_1/data/
dataDir = cd ;
meshDir = fullfile(dataDir, 'deconvolved_16bit', 'msls_output') ;

% PATHS ==================================================================
origpath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(origpath))
aux_paths_and_colors
cd(dataDir)

% DEFINE NEW MASTER SETTINGS
overwrite_masterSettings = false ;
if overwrite_masterSettings || ~exist('./masterSettings.mat', 'file')
    % Metadata about the experiment
    stackResolution = [.2619 .2619 .2619] ;
    nChannels = 1 ;
    channelsUsed = 1 ;
    timePoints = 10:263 ;
    ssfactor = 4 ;
    % whether the data is stored inverted relative to real position
    flipy = true ; 
    timeInterval = 1 ;  % physical interval between timepoints
    timeUnits = 'min' ; % physical unit of time between timepoints
    spaceUnits = '$\mu$m' ; % physical unit of time between timepoints
    scale = 0.04 ;      % scale for conversion to 16 bit
    % file32Base = 'TP%d_Ch0_Ill0_Ang0,45,90,135,180,225,270,315.tif'; 
    file32Base = 'TP%d_Ch0_Ill0_Ang0,60,120,180,240,300.tif'; 
    % file32Base = 'TP%d_Ch0_Ill0_Ang0,60,120,180,240,300.tif'; 
    fn = 'Time_%06d_c1_stab';
    fn_prestab = 'Time_%06d_c1.tif';
    set_preilastikaxisorder = 'xyzc' ;
    swapZT = 1 ;
    masterSettings = struct('stackResolution', stackResolution, ...
        'nChannels', nChannels, ...
        'channelsUsed', channelsUsed, ...
        'timePoints', timePoints, ...
        'ssfactor', ssfactor, ...
        'flipy', flipy, ...
        'timeInterval', timeInterval, ...
        'timeUnits', timeUnits, ...
        'spaceUnits', timeUnits, ...
        'scale', scale, ...
        'file32Base', file32Base, ...
        'fn', fn,...
        'fn_prestab', fn_prestab, ...
        'swapZT', swapZT, ...
        'set_preilastikaxisorder', set_preilastikaxisorder, ...
        't0_for_phi0', 110, ... % 40 for mef2 single channel, 110 for CAAX excellent
        'nU', 100, ...  % 150 for mef2 data with posterior midgut loop
        'nV', 100); 
    disp('Saving masterSettings to ./masterSettings.mat')
    if exist('./masterSettings.mat', 'file')
        ui = input('This will overwrite the masterSettings. Proceed (Y/n)?', 's') ;
        if ~isempty(ui) && (strcmp(ui(1), 'Y') || strcmp(ui(1), 'y'))
            save('./masterSettings.mat', 'masterSettings')
            loadMaster = false ;
        else
            disp('Loading masterSettings from disk instead of overwriting')
            loadMaster = true ;
        end
    else
        save('./masterSettings.mat', 'masterSettings')
        loadMaster = false ;
    end
else
    loadMaster = true ;
end


% LOAD EXISTING MASTER SETTINGS
disp('Loading masterSettings from ./masterSettings.mat')
load('./masterSettings.mat', 'masterSettings')
% Unpack existing master settings
stackResolution = masterSettings.stackResolution ;
nChannels = masterSettings.nChannels ;
channelsUsed = masterSettings.channelsUsed ;
timePoints = masterSettings.timePoints ;
ssfactor = masterSettings.ssfactor ;
% whether the data is stored inverted relative to real position
flipy = masterSettings.flipy ; 
timeInterval = masterSettings.timeInterval ;  % physical interval between timepoints
timeUnits = masterSettings.timeUnits ; % physical unit of time between timepoints
spaceUnits = masterSettings.spaceUnits ; % unit of distance of full resolution data pixels ('$\mu$m')
scale = masterSettings.scale ;      % scale for conversion to 16 bit
file32Base = masterSettings.file32Base ; 
fn = masterSettings.fn ;
fn_prestab = masterSettings.fn_prestab ;
set_preilastikaxisorder = masterSettings.set_preilastikaxisorder ;
swapZT = masterSettings.swapZT ;
t0_for_phi0 = masterSettings.t0_for_phi0 ;
nU = masterSettings.nU ;
nV = masterSettings.nV ;


% Fill in
scale = masterSettings.scale ;      % scale for conversion to 16 bit
file32Base = masterSettings.file32Base ; 
fn = masterSettings.fn ;
fn_prestab = masterSettings.fn_prestab ;
set_preilastikaxisorder = masterSettings.set_preilastikaxisorder ;

dir32bit = fullfile(dataDir, 'deconvolved_32bit') ;
dir16bit = fullfile(dataDir, 'deconvolved_16bit') ;
dir16bit_prestab = fullfile(dir16bit, 'data_pre_stabilization') ;

% I. INITIALIZE ImSAnE PROJECT ===========================================
% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.
cd(dir16bit)
dataDir = cd ;
projectDir = dataDir ;
% [ projectDir, ~, ~ ] = fileparts(matlab.desktop.editor.getActiveFilename); 
cd(projectDir);
if projectDir(end) ~= '/'
    projectDir = [projectDir '/'];
end

% Start by creating an experiment object, optionally pass on the project
% directory (otherwise it will ask), and change into the directory of the
% data.  This serves as a front-end for data loading, detection, fitting
% etc.
xp = struct() ;
% A filename base template - to be used throughout this script
fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fn, '.tif'];
fileMeta.nChannels          = nChannels;
fileMeta.timePoints         = timePoints ;
fileMeta.stackResolution    = stackResolution;
fileMeta.swapZT             = 1;

% first_tp is also required, which sets the tp to do individually.
first_tp = 1 ;
expMeta                     = struct();
expMeta.channelsUsed        = channelsUsed ;
expMeta.channelColor        = 1;
expMeta.description         = 'Drosophila gut';
expMeta.dynamicSurface      = 1;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.fitTime             = fileMeta.timePoints(first_tp);
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';

% Now set the meta data in the experiment.
xp.fileMeta = fileMeta;
xp.expMeta = expMeta;

clear fileMeta expMeta
 
% % SET DETECT OPTIONS =====================================================
% % Must run this section for later functionality.
% % Load/define the surface detection parameters
% msls_detOpts_fn = fullfile(projectDir, 'msls_detectOpts.mat') ;
% load(msls_detOpts_fn, 'detectOptions')
% 
% % Overwrite certain parameters for script structure
% detectOptions.fileName = sprintf( fn, xp.currentTime ) ;
% detectOptions.run_full_dataset = false ;
% detectOptions.ms_scriptDir = ms_scriptDir ;
% meshDir = detectOptions.mslsDir ;
% % These are now in QS.
% % meshFileBase = [ofn_smoothply '%06d'] ;
% % alignedMeshBase = [ofn_smoothply '%06d_APDV_um'] ;
% 
% % Set detect options ------------------------------------------------------
% xp.setDetectOptions( detectOptions );
% disp('done')
xp.detectOptions = struct() ;
xp.detectOptions.ssfactor = 4 ;
xp.detector.options.ofn_smoothply = 'mesh_' ;
xp.tIdx = 1:length(timePoints) ;

%% QS DEFINITION
opts.meshDir = meshDir ;
opts.flipy = flipy ;
opts.timeInterval = timeInterval ;
opts.timeUnits = timeUnits ;
opts.spaceUnits = '$\mu$m' ;
opts.nV = nV ;
opts.nU = nU ;
opts.normalShift = 10 ;
opts.a_fixed = 2.0 ;
opts.adjustlow = 1.00 ;                  %  floor for intensity adjustment
opts.adjusthigh = 99.9 ;                 % ceil for intensity adjustment (clip)
% opts.adjustlow = 0 ;                  %  floor for intensity adjustment
% opts.adjusthigh = 0 ;                 % ceil for intensity adjustment (clip)
opts.phiMethod = 'curves3d' ;

opts.lambda_mesh = 0 ;
opts.lambda = 0 ;
opts.lambda_err = 0 ;
%  opts.lambda_mesh = 0.002 ;
%  opts.lambda = 0.01 ;
%  opts.lambda_err = 0.01 ;
 
disp('defining QS')
QS = QuapSlap(xp, opts) ;
disp('done')
%%
% %% Cell Segmentation
% % options = struct() ;
% % options.overwrite = true;
% % options.cellSize = 50 ;
% % options.strelRadius = 0 ;
% % options.gaussKernel = 1 ;
% % options.heightMinimum = 2.2 ;
% % options.coordSys= 'sprsme';
% % options.timePoints = 213:15:258;% 93:15:200 ;
% % QS.generateCellSegmentation2D(options) 
% 
% % Now manual corrections in GIMP
% options = struct() ;
% options.overwrite = true ;
% options.overwriteImages = true;
% options.timePoints = 186;
% QS.processCorrectedCellSegmentation2D(options) 
% 
% %% Generate 3d
% options = struct() ;
% options.overwrite = false;
% options.timePoints =96:10:206;
% options.correctedSegmentation = true;
% %QS.generateCellSegmentation3D(options) 
% QS.plotSegmentationStatisticsLobes(options)
% 
% %% 
% options = struct() ;
% options.overwrite = false;
% options.timePoints = 96:10:206;
% options.correctedSegmentation = true;
% QS.generateCellSegmentationPathlines3D(options)


%% Testing
cd '/home/yuzhenglin/membraneproject'
load("confirmed_t1.mat");
load('not_always_neighbors.mat');
load("t1_overlap.mat");
options = struct() ;
options.overwrite = false;
options.preview = false ;
timePoints = 96:2:206;
options.scaleByMetricComponents = true;
subdir = '/home/yuzhenglin/membraneproject';
labelDir = fullfile(subdir, 'labeled_groundTruth');
% consider each T1 in question
overlap_t1 = suspects(:, [1,3]);
split_events = setdiff(split_events, overlap_t1, 'rows');
theta_merge = zeros(size(merge_events, 1), 1);
t1timepoints_merge = zeros(size(merge_events, 1), 1);
theta_split = zeros(size(split_events, 1), 1);
t1timepoints_split = zeros(size(split_events, 1), 1);

for pp = 1:  size(merge_events, 1)
    cellpair = not_always_neighbors_label(merge_events(pp,1),:);
    load_label = load(fullfile(labelDir, sprintf("tracks_label_%06d.mat", timePoints(merge_events(pp,2)+1)))); % +1 since T1 happens a frame after 
    label = load_label.imlabel;
    cells = zeros(size(label));
    cells(label == cellpair(1)) = 1;
    cells(label == cellpair(2)) = 2;
    regp = regionprops(cells,'Centroid');
    ctd = vertcat(regp.Centroid);
    % load cell pairs and their associated timepoints
    %pts = t1info(pp, 1:2) ;
    %t1tp = t1info(pp, 3) ;

    % Find their local pullback coordinate orientation
    QS.setTime(timePoints(merge_events(pp,2)+1)) ;
    [subm, newpts] = QS.computeLocalSurfacePatch(ctd, options);

    % Get angle
    dY = newpts(2, 2) - newpts(1, 2) ;
    dX = newpts(2, 1) - newpts(1, 1) ;
    theta_merge(pp) = atan2(dY, dX) ;
    t1timepoints_merge(pp) = timePoints(merge_events(pp,2)+1) ;
end

for pp = 1:  size(split_events, 1)
    cellpair = not_always_neighbors_label(split_events(pp,1),:);
    load_label = load(fullfile(labelDir, sprintf("tracks_label_%06d.mat", timePoints(split_events(pp,2)+1)))); % +1 since T1 happens a frame after 
    label = load_label.imlabel;
    cells = zeros(size(label));
    cells(label == cellpair(1)) = 1;
    cells(label == cellpair(2)) = 2;
    regp = regionprops(cells,'Centroid');
    ctd = vertcat(regp.Centroid);
    % load cell pairs and their associated timepoints
    %pts = t1info(pp, 1:2) ;
    %t1tp = t1info(pp, 3) ;

    % Find their local pullback coordinate orientation
    QS.setTime(timePoints(split_events(pp,2)+1)) ;
    [subm, newpts] = QS.computeLocalSurfacePatch(ctd, options);

    % Get angle
    dY = newpts(2, 2) - newpts(1, 2) ;
    dX = newpts(2, 1) - newpts(1, 1) ;
    theta_split(pp) = atan2(dY, dX) ;
    t1timepoints_split(pp) = timePoints(split_events(pp,2)+1) ;
end

theta_merge = mod(theta_merge, pi) ;
theta_split = mod(theta_split, pi) ;

save('confirmed_t1_angles.mat', 'theta_merge', 'theta_split', ...
    't1timepoints_merge', 't1timepoints_split')


%% plot statistics
%load("confirmed_t1_angles.mat");

t_min = 96;
t_max = 206;
a_min = -pi/4;
a_max = 3*pi/4;
t_size = 56;
a_size = 10;

dt = (t_max-t_min)/t_size;
da = (a_max-a_min)/a_size;



kymo_merge = zeros(t_size, a_size);
for tt = 1: t_size
    for aa = 1: a_size
        t_in_range = find(t1timepoints_merge >= t_min+ (tt-1)*dt ...
                            & t1timepoints_merge < t_min+ tt*dt);
        a_in_range = find(theta_merge >= a_min+ (aa-1)*da ...
                            & theta_merge < a_min+aa*da);
        kymo_merge(tt,aa) = numel(intersect(t_in_range, a_in_range));
    end
end

kymo_split = zeros(t_size, a_size);
for tt = 1: t_size
    for aa = 1: a_size
        t_in_range = find(t1timepoints_split >= t_min+ (tt-1)*dt ...
                            & t1timepoints_split < t_min+ tt*dt);
        a_in_range = find(theta_split >= a_min+ (aa-1)*da ...
                            & theta_split < a_min+aa*da);
        kymo_split(tt,aa) = numel(intersect(t_in_range, a_in_range));
    end
end

%%
% subplot(1,2,1)
% imagesc(kymo_merge);
% set(gca,'YDir','normal')
% set(gca,'Xtick',linspace(1,a_size,5),'XTickLabel',{'-\pi/4', '0', '\pi/4', '\pi/2', '3\pi/4'})
% xlabel("merge angle")
% set(gca,'Ytick',linspace(1,t_size,10),'YTickLabel',round(linspace(t_min,t_max,10)))
% 
% colorbar;
% hold on
% subplot(1,2,2)
% imagesc(kymo_split);
% set(gca,'YDir','normal')
% set(gca,'Xtick',linspace(1,a_size,5),'XTickLabel',{'-\pi/4', '0', '\pi/4', '\pi/2', '3\pi/4'})
% xlabel("split angle")
% set(gca,'Ytick',linspace(1,t_size,10),'YTickLabel',round(linspace(t_min,t_max,10)))
% colorbar;


%% first half
merge_ap = sum(kymo_merge(1:36, 1:5), 'all');
merge_dv = sum(kymo_merge(1:36, 6:10), 'all');
split_ap = sum(kymo_split(1:36, 1:5), 'all');
split_dv = sum(kymo_split(1:36, 6:10), 'all');
contributing1 = merge_dv + split_ap;
detracting1 = split_dv + merge_ap;

AP1 = contributing1;
APe1 = sqrt(AP1);
DV1 = detracting1;
DVe1 = sqrt(DV1);
pval1 = normcdf((contributing1 - detracting1)/sqrt(contributing1 + detracting1));



%% second half
merge_ap = sum(kymo_merge(37:56, 1:5), 'all');
merge_dv = sum(kymo_merge(37:56, 6:10), 'all');
split_ap = sum(kymo_split(37:56, 1:5), 'all');
split_dv = sum(kymo_split(37:56, 6:10), 'all');

contributing2 = merge_dv + split_ap;
detracting2 = split_dv + merge_ap;
AP2 = contributing2;
APe2 = sqrt(AP2);
DV2 = detracting2;
DVe2 = sqrt(DV2);
pval2 = normcdf(-(contributing2 - detracting2)/sqrt(contributing2 + detracting2));

%% Significance plot for first and second half of timecourse
hf = figure('Position', [100 100 800 700], 'units', 'centimeters');
clf;
Y = [AP1, DV1, AP2, DV2];
E = cat(3, APe1, DVe1, APe2, DVe2) ;

Colors = [
    0.90    0.55    0.55
    0.62    0.76    0.84
    0.89    0.10    0.11
    0.12    0.47    0.70
    ];
Colors = reshape(Colors, [2 2 3]);

P = [pval1,  pval1, NaN, NaN; pval1, pval1, NaN, NaN; ...
    NaN, NaN, pval2, pval2; NaN, NaN, pval2, pval2;];
% Make P symmetric, by copying the upper triangle onto the lower triangle
% PT = P';
% lidx = tril(true(size(P)), -1);
% P(lidx) = PT(lidx);

superbar([1, 2, 3, 4], Y, 'E', E, 'P', P, 'BarFaceColor', Colors, 'Orientation', 'v', ...
    'ErrorbarStyle', 'I', 'PLineOffset', 0.1, 'PStarShowGT', false);
ylabel('number of T1 events')
xticks([1, 2, 3, 4])
set(gca, 'xticklabels', {'AP', 'DV', 'AP', 'DV', })
sgtitle('T1 orientations', 'interpreter', 'latex')
xlabel('$t < 40$ min \hspace{25pt} $t>40$ min', 'interpreter', 'latex')
saveas(gcf, 'confirmed_t1_contribution_timeCourse.pdf')


%% -Merge + split
%theta_merge = mod(theta_merge, pi) ;
%theta_split = mod(theta_split, pi) ;
% theta_merge(theta_merge > 3*pi/4) = theta_merge(theta_merge > 3*pi/4) - pi ;
% theta_split(theta_split > 3*pi/4) = theta_split(theta_split > 3*pi/4) - pi ;
% % number of T1s from merging oriented such that elongation is along AP
% nAPmerge1 = length( find(theta_merge > pi/4 )) ;
% % number of T1s from splitting oriented such that elongation is along AP
% nAPsplit0 = length( find(theta_split < pi/4 )) ;
% % number of T1s from merging oriented such that elongation is along DV
% nDVmerge0 = length( find(theta_merge < pi/4 )) ;
% % number of T1s from splitting oriented such that elongation is along DV
% nDVsplit1 = length( find(theta_split > pi/4 )) ;
% bar([0, pi/2], [nDVmerge0, nAPmerge1])
% hold on;
% bar([0, pi/2], [nAPsplit0, nDVsplit1], 'facealpha', 0.5)
% 
% % so we measure an AP t1 rate over the whole timecourse:
% AP = nAPmerge1 + nAPsplit0 ;
% DV = nDVmerge0 + nDVsplit1 ;
% % uncertainties are sqrt(N) ?
% APe = sqrt(AP) ;
% DVe = sqrt(DV) ;
% % Convert to p value
% zscore = -(AP - DV) / sqrt(APe.^2 + DVe.^2) ;
% pval = normcdf(zscore) ;


%save('confirmed_t1_angles.mat', 'theta_merge', 'theta_split', ...
%    't1timepoints_merge', 't1timepoints_split')

% save('confirmed_t1_angles_stats.mat', 'AP', 'DV', 'nAPmerge1', 'nAPsplit0', ...
%     'nDVmerge0', 'nDVsplit1')
% 
% % Significance plot
% hf = figure('Position', [100 100 400 300], 'units', 'centimeters');
% clf;
% Y = [AP, DV];
% E = cat(3, APe, DVe) ;
% 
% Colors = [
%     0.90    0.55    0.55
%     0.62    0.76    0.84
%     0.89    0.10    0.11
%     0.12    0.47    0.70
%     ];
% Colors = reshape(Colors, [2 2 3]);
% 
% P = [pval,  pval; pval, pval];
% % Make P symmetric, by copying the upper triangle onto the lower triangle
% % PT = P';
% % lidx = tril(true(size(P)), -1);
% % P(lidx) = PT(lidx);
% 
% superbar(Y, 'E', E, 'P', P, 'BarFaceColor', Colors, 'Orientation', 'v', ...
%     'ErrorbarStyle', 'I', 'PLineOffset', 0.1, 'PStarShowGT', false);
% ylabel('number of T1 events')
% xticks([1, 2])
% set(gca, 'xticklabels', {'contributing', 'detracting'})
% sgtitle('T1 orientations', 'interpreter', 'latex')
% saveas(gcf, 'confirmed_t1_contribution.pdf')
% 
% % for pp = 1: size(split_events, 1)
% %     % load cell pairs and their associated timepoints
% %     %pts = t1info(pp, 1:2) ;
% %     %t1tp = t1info(pp, 3) ;
% % 
% %     % Find their local pullback coordinate orientation
% %     QS.setTime(t1tp) ;
% %     %[subm, newpts] = QS.computeLocalSurfacePatch(pts, options);
% % 
% %     % Get angle
% %     dY = newpts(2, 2) - newpts(1, 2) ;
% %     dX = newpts(2, 1) - newpts(1, 1) ;
% %     theta(pp) = atan2(dY, dX) ;
% %     t1timepoints(pp) = t1tp ;
% % end
% % save results as .mat
% %save(outfn, 'theta', 't1timepoints')
% 
% % Make kymo 
% % imagesc(AngleTimeMat)
% % cb = colorbar() ;
% % ylabel(cb, 'occurences')
% 
% 
% %% Significance plot for first and second half of timecourse
% hf = figure('Position', [100 100 400 300], 'units', 'centimeters');
% clf;
% Y = [AP1, DV1, AP2, DV2];
% E = cat(3, APe1, DVe1, APe2, DVe2) ;
% 
% Colors = [
%     0.90    0.55    0.55
%     0.62    0.76    0.84
%     0.89    0.10    0.11
%     0.12    0.47    0.70
%     ];
% Colors = reshape(Colors, [2 2 3]);
% 
% P = [pval1,  pval1, NaN, NaN; pval1, pval1, NaN, NaN; ...
%     NaN, NaN, pval2, pval2; NaN, NaN, pval2, pval2;];
% % Make P symmetric, by copying the upper triangle onto the lower triangle
% % PT = P';
% % lidx = tril(true(size(P)), -1);
% % P(lidx) = PT(lidx);
% 
% superbar([1, 2, 3, 4], Y, 'E', E, 'P', P, 'BarFaceColor', Colors, 'Orientation', 'v', ...
%     'ErrorbarStyle', 'I', 'PLineOffset', 0.1, 'PStarShowGT', false);
% ylabel('number of T1 events')
% xticks([1, 2, 3, 4])
% set(gca, 'xticklabels', {'AP', 'DV', 'AP', 'DV', })
% sgtitle('T1 orientations', 'interpreter', 'latex')
% xlabel('$t < 40$ min \hspace{25pt} $t>40$ min', 'interpreter', 'latex')
% saveas(gcf, 'confirmed_t1_contribution_timeCourse.pdf')
