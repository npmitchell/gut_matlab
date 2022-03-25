%% eLife reviews for anisotropy 

datdir = '/Users/npmitchell/Dropbox/Soft_Matter/PAPER/gut_paper/figure_drafting/figure2_kinematics/seg3d_corrected' ;
fns = dir(fullfile(datdir, 'Time*.mat')) ;

% Colormap
addpath(genpath('~/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/gut_matlab/plotting/'))
colors = colormap(0.9 * cmap('I1', 'N', length(fns))) ;


%% WARMUP: UV coordinates
% Plot each timepoint
clf ;
for tidx = 1:length(fns)
    fn = fullfile(fns(tidx).folder, fns(tidx).name) ;
    seg = load(fn) ;
    seg = seg.seg3d ;
    stat = seg.statistics ;
    ap = stat.apStats.aspectWeighted.apBins ;
    c2t = stat.apStats.aspectWeighted.apCos2Theta ;
    cs = stat.apStats.aspectWeighted.apCos2ThetaStd ;
    ce = stat.apStats.aspectWeighted.apCos2ThetaSte ;
    
    % Plot the anisotropy
    lineProps = {'-','color', colors(tidx, :)} ;
    % shadedErrorBar(ap, c2t, cs, 'lineProps', lineProps)
    hold on;
    shadedErrorBar(ap, c2t, ce, 'lineProps', lineProps)
    
end

cb = colorbar ;

%% Find where cells go

% build needed files
QS = QuapSlap(xp) ;
cellVertexPathlineFn = fullfile(QS.dir.segmentation, 'pathlines', ...
    sprintf('cellVertexPathlines_%06dt0.mat', QS.t0set())) ;
load(cellVertexPathlineFn, 'segVertexPathlines2D', ...
         'segVertexPathlines3D', 'cellIDs')
segP2 = segVertexPathlines2D ;
segP3 = segVertexPathlines3D ;

% find position where folds occur in pullback space



%% 

    
 