% Compute the surface area and volume of a time series of meshes
% NPMitchell 2019

clear; close all; clc;
addpath('/mnt/data/code/gut_matlab/plotting/')
addpath('/mnt/data/code/gut_matlab/mesh_handling/')

% Prepare path for meshes ===========================================
% path = '../../data/48Ygal4-UAShisRFP/2016021015120_objFiles/';
% meshes = dir([path, 'pointCloud_T*.ply']) ;
%path = '../data/48Ygal4-UAShisRFP/20170329_objFiles/';
%meshes = dir([path, '*_ascii.ply']) ;
% path = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/data/48Ygal4-UAShisRFP/2016021015120_objFiles/cleaned/';
% meshes = dir([path, 'cleaned_pointCloud*_mesh.ply']) ;
%%%%%%%%%%%%%%%%%%%%%
% 48Ygal4-UASHistRFP
% path = '/mnt/crunch/48Ygal4-UAShistRFP/201901021550_folded_2part/Time12views_60sec_1.2um_25x_4/data/deconvolved_16bit/msls_output_nu0p10_s1_pn4_ps4_l1_l1/';
% path = '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/';
path = './' ;
meshes = dir(fullfile(path, 'mesh_apical_ms_stab_0*.ply'))  ; 
%%%%%%%%%%%%%%%%%%%%%

thres_fracda = 1.0 ;
outdir = path ;
clims = [-1 1] ;
check = false ;
dx = 0.2619 ;
ii = 1;
xlims = [0, 200] ;
ylims = [-50, 50] ;
zlims = [-50, 50] ;

% Create output dirs
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Prepare for iteration
% chisums is the integral chirality
dt = 1;
todo = 1:dt:min(length(meshes), 111) ;
vvs = zeros(length(todo), 1) ;
aas = zeros(length(todo), 1) ;
dmyk = 1;

 %% Compute surface area and volume for each mesh
for ii = todo
%     try
        % Load the mesh
        meshfn = fullfile(meshes(ii).folder, meshes(ii).name) ;
        disp(['analyzing ', meshfn])
        [tri, pts] = ply_read(meshfn, 'tri') ;
        pts = transpose(pts) * dx;
        pts(:, 3) = - pts(:, 3) + mean(pts(:, 3));
        tri = transpose(tri) ;

        % View Result --------------------------------------------------------
        if check
            trisurf(tri, pts(:,1), pts(:,2), pts(:,3));
            axis equal
        end

        % Compute surface area and volume
        [vv, aa] = meshVolumeArea(pts, tri) ;
        if dmyk > 1
            fracda = abs(abs(aa) - abs(prev_a)) / abs(prev_a) ;
            disp(['fracda = ', num2str(fracda)])
            if fracda > thres_fracda
                vvs(dmyk) = NaN ;
                aas(dmyk) = NaN ;
                disp('Fractional change in aa is too great, skipping')
            else
                vvs(dmyk) = abs(vv) ;
                aas(dmyk) = aa ;
                prev_a = aa ;
            end
        else
            vvs(dmyk) = abs(vv) ;
            aas(dmyk) = aa ;
            prev_a = aa ;
            disp(['First step: aa = ', num2str(aa)])
        end

        dmyk = dmyk + 1;
%     catch
%         disp('could not process this time point')
%         vvs(dmyk) = NaN ;
%         aas(dmyk) = NaN ;
%         dmyk = dmyk + 1 ;
%     end
end

% Save the data
save(fullfile(outdir, 'surfacearea_volume_stab.m'), 'aas', 'vvs', 'dt')

% Save the image
figh = figure();
hold on;
ah = plot(1:dt:dt*length(aas), aas / aas(1)) ;
vh = plot(1:dt:dt*length(aas), vvs / vvs(1)); 
legend('area', 'volume')
title('Surface area and volume')
xlabel('Time [min]')
ylabel('Normalized area or volume')
saveas(figh, fullfile(outdir, 'area_volume_over_time_stab.pdf'))
saveas(figh, fullfile(outdir, 'area_volume_over_time_stab.png'))
