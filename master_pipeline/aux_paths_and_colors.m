disp('Adding paths')
ms_scriptDir = '/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/' ;
codepath = '/mnt/data/code/' ;
gutpath = fullfile(codepath, 'gut_matlab') ;
meshlabCodeDir = fullfile(codepath, 'meshlab_codes') ;

% Add the addpath_recurse function
addpath(fullfile(codepath, 'gut_matlab', 'addpath_recurse'))

% Codepath + gutpath folders
addpath('./plotting_functions')
addpath_recurse(fullfile(codepath, 'gptoolbox'))
addpath(fullfile(gutpath, 'master_pipeline'))
addpath_recurse(fullfile(gutpath, 'data_handling'))
addpath_recurse(fullfile(gutpath, 'mesh_handling'))
addpath(fullfile(gutpath, 'basics'))
addpath(fullfile(gutpath, 'distanceTransform'))
addpath_recurse(fullfile(gutpath, 'plotting'))
addpath_recurse(fullfile(gutpath, 'master_pipeline', 'nes_functions'))
addpath(fullfile(gutpath, 'savgol'))
addpath_recurse(fullfile(gutpath, 'toolbox_fast_marching/'));
addpath(genpath(fullfile(gutpath, 'euclidean_orbifolds')));
addpath(genpath(fullfile(gutpath, 'TexturePatch')));
addpath_recurse(fullfile(gutpath, ['axisymmetric_pullbacks' filesep])) ;
addpath_recurse(fullfile(gutpath, 'plotting' )) ;
addpath_recurse(fullfile(gutpath, 'h5_handling' )) ;
addpath_recurse(fullfile(gutpath, 'tiff_handling')) ;
addpath_recurse(fullfile(gutpath, 'geometry')) ;
addpath_recurse(fullfile(gutpath, 'curve_functions')) ;
addpath(fullfile(gutpath, ['ExtPhaseCorrelation' filesep])) ;
addpath(fullfile(gutpath, 'savgol')) ;
addpath(fullfile(codepath, 'DEC')) ;
addpath(fullfile(codepath, 'TexturePatch_for_git', 'TexturePatch')) ;
addpath_recurse('/mnt/crunch/djcislo/MATLAB/CGAL_Code/')
addpath('/mnt/data/code/RicciFlow_MATLAB/')
addpath(fullfile(gutpath, 'tissueAnalysisSuite'))
% addpath(genpath('/mnt/crunch/djcislo/MATLAB/TexturePatch'));
disp('done')

%% Define some colors
[colors, color_names] = define_colors() ;
blue = [0 0.4470 0.7410] ;
orange = [0.8500 0.3250 0.0980] ;
yellow = [0.9290, 0.6940, 0.1250] ;
purple = [0.4940, 0.1840, 0.5560] ;
green = [0.4660, 0.6740, 0.1880] ;
sky = [0.3010, 0.7450, 0.9330] ;
red = [0.6350, 0.0780, 0.1840] ;
brick = [0.800000 0.250000 0.330000] ;
light_green =[0.560000 0.930000 0.560000] ;
light_gray = [0.830000 0.830000 0.830000] ;
bwr = diverging_cmap([0:0.01:1], 1, 2) ;
disp('done adding colors')