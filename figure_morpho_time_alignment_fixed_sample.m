% figure_morpho_time_alignment_fixed_sample
clearvars
addpath(genpath('/mnt/data/code/gut_matlab/'))

%% Options
sigmaTime = 5 ;         % smoothing window for errorbar plot
ssfactorMedianThres = 8 ;
rsubsampling = 1 ;
ssample_factor = 40 ;
t0ref = 123 ;

%% Directories
refdir = '/mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/' ;
datdir = '/mnt/data/antibodies_lightsheet/48YGAL4CAAXmCh_antp8C11_1t50/202108171735_e2_16a_1p2um_0p3ms0p5ms_3mWGFP_15mWRFP/data/' ;

refAlignedMeshDir = fullfile(refdir, 'msls_output', 'alignedMesh') ;
meshdir = fullfile(datdir, 'msls_output') ;
alignedMeshDir = fullfile(meshdir, 'aligned_meshes') ;

ssrDir = fullfile(meshdir, 'time_alignment') ;

% Load the mesh
% mesh = read_ply_mod(meshFn) ;

%% Build mesh cell array and timestamps
mca = {} ;
timstamps = [] ;
ameshDirs = {refAlignedMeshDir, alignedMeshDir} ;
meshes = dir(fullfile(refAlignedMeshDir, 'mesh_*_APDV_um.ply')) ;
ntps = 0 ;
for kk = 1:length(meshes)
    mca{1, kk} = meshes(kk);
    % Get timestamp of this mesh
    tmp = strsplit(mca{1, kk}.name, '_0') ;
    tmp = strsplit(tmp{2}, '_APDV') ;
    timestamps(1, kk) = str2double(tmp{1});
    ntps(1) = ntps(1)+1 ;
end
meshes = dir(fullfile(alignedMeshDir, 'mesh_*_APDV_um_original.ply')) ;
ntps(2) = 0 ;
for kk = 1:length(meshes)
    mca{2, kk} = meshes(kk);
    % Get timestamp of this mesh
    tmp = strsplit(mca{2, kk}.name, '_0') ;
    tmp = strsplit(tmp{2}, '_APDV') ;
    timestamps(2, kk) = str2double(tmp{1});
    ntps(2) = ntps(2)+1 ;
end


%% Compare mesh to morphological timeline
refID = 1 ;
clist = 1 ; lastTP = 1 ;
rlist = 1:length(mca(1, :)) ;
cExptID = '48YGAL4CAAXmCh_antp8C11_1t50_202108171735_e2_16a' ;
refExptID = 'caax201902072000'; 
options = struct() ;
options.cc = 2 ;
options.refID = refID ;
options.cExptID = cExptID ;
options.refExptID = refExptID ;
options.ssrccDir = ssrDir ;
options.ssfactorMedianThres = ssfactorMedianThres ;
options.rsubsampling = rsubsampling ;
options.ssample_factor = ssample_factor ;
options.t0ref = t0ref ;
[ssrM, minddssr, minname, minerror, ...
    minweights, ssr_minimum ] = ...
    computeSSRMatrixMeshes(mca, clist, lastTP, rlist, options) ;

%% Save result
timestamp_idx = timestamps(refID, minddssr) ;
timestamp = timestamp_idx-t0ref ;
uncertainty = minerror ;
mean_distance_from_reference = sqrt(ssr_minimum) ;
save(fullfile(ssrDir, 'timestamp.mat'), 'timestamp', 'timestamp_idx', 'uncertainty', ...
    'mean_distance_from_reference')
write_txt_with_header(...
    fullfile(ssrDir, 'timestamp.txt'), ...
    [timestamp, uncertainty, mean_distance_from_reference], ...
    ['timestamp, uncertainty, mean_distance_from_reference'])

% Save pdf for figure inset
close all
figure('Units', 'centimeters', 'Position', [0 0 8 6])
plot((timestamps(refID, :)-t0ref)/60, ssrM, '.')
hold on; 
errorbar(timestamp/60, ssr_minimum, NaN, NaN, uncertainty/60, uncertainty/60, 'o')
xlabel('morphological timeline [a.u.]')
ylabel(['geometric mismatch [' char(956) 'm^2]'])
xlim([timestamp-5*uncertainty, timestamp+5*uncertainty]/60)
saveas(gcf, fullfile(ssrDir, ['timestamping_against_' refExptID '_ssr.pdf']))

close all
figure('Units', 'centimeters', 'Position', [0 0 8 6])
plot((timestamps(refID, :)-t0ref)/60, sqrt(ssrM), '.')
hold on; 
errorbar(timestamp/60, sqrt(ssr_minimum), NaN, NaN, uncertainty/60, uncertainty/60, 'o')
xlabel('morphological timeline [a.u.]')
ylabel(['geometric mismatch [' char(956) 'm]'])
xlim([timestamp-5*uncertainty, timestamp+5*uncertainty]/60)
saveas(gcf, fullfile(ssrDir, ['timestamping_against_' refExptID '_sqrtssr.pdf']))


% % CORRESPONDENCE CURVES
% options = struct() ;
% options.overwrite = false ;
% options.cc = 2 ;
% options.refID = 1 ;
% options.ssrDir = ssrDir;
% options.timestamps = timestamps;
% options.ntps = ntps ;
% options.refExptID = refExptID ;
% options.cExptID = cExptID ;
% options.sigmaTime = sigmaTime ;
% options.rsubsampling = rsubsampling ;
% [corrPath, corrRaw, corrError, ssrPath, ssrPathError] = ...
%     computeCorrespondenceCurve(ssrM, minddssr, minerror, options) ;
