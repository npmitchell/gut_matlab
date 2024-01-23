%% Align meshes in Space and Time to a reference mesh
%
% Isaac Breinyn 2019
%
% This is an edited combination of mesh_temporal_alignment.m and
% plot_aligned meshes.m
% (Master script)
%
% To run before: compute_mesh_surfacearea_volume.m, setup.m
%
% Note: All scaling and alignment is done relative to CAAX Excellent

clear
close all
clc
addpath('/mnt/data/code/gut_matlab/curve_functions/')

%% USER OPTIONS %%
%\/\/\/\/\/\/\/\/%

dim = 'three' ; % 'two' for 2D plots, 'three' for 3D plots
orient = 'lateral' ; % (ONLY FOR 2D PLOTS) Specify 'lateral' for a lateral view, and 'posterior' for... well you get it
slice_AP = 50; % The location of the yz slice in the AP direction


%% Prepare paths, Load Data, Calculate Scaling Factors
outdir = '/mnt/data/analysis/SA_volume/' ;
markers = {'caax', 'hrfp', 'la'} ;
labels = {'Membrane', 'Nuclei', 'Actin'};

% Prepare paths to data
rootdir = '/mnt/crunch/' ;

% membrane
caax_root = [rootdir '48Ygal4UASCAAXmCherry'] ;
caax_paths = {'201902072000_excellent/Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/', ...
    '201903211930_great/Time6views_60sec_1p4um_25x_1p0mW_exp0p150/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/'} ;

% nuclei
hrfp_root = [rootdir '48Ygal4-UAShistRFP/'] ;
hrfp_paths = {'201901021550_folded_2part/Time12views_60sec_1.2um_25x_4/data/deconvolved_16bit/msls_output_prnu0_prs0_nu0p10_s1p00_pn4_ps4_l1_l1/', ...
    '201904031830_great/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/', ...
    '201903312000_closure_folding_errorduringtwist/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2_folding/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/',...
    %'201903312000_closure_folding_errorduringtwist/Time4views_180sec_1p4um_25x_1p0mW_exp0p35_dorsalclosure/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/',...
    } ;

% actin
la_root = [rootdir '48YGal4UasLifeActRuby'] ;
la_paths = {'201904021800_great/Time6views_60sec_1p4um_25x_1p0mW_exp0p150_3/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/'};

roots = {caax_root, hrfp_root, la_root} ;
paths = {{caax_paths}, {hrfp_paths}, {la_paths}};
npaths = [length(caax_paths), length(hrfp_paths), length(la_paths)] ;

mca = {} ;

% Iterate over each marker
for mi = 1:length(markers)
    
    % Obtain the label for this marker
    label = labels{mi} ;
    these_paths = paths{mi} ;
    these_paths = these_paths{1} ;
    
    % Cycle through all datasets of this marker
    for j=1:length(these_paths)
        mpath = these_paths{j} ;
        matdir = fullfile(roots{mi}, mpath) ;
        fn = fullfile(fullfile(roots{mi}, mpath), 'surfacearea_volume_stab.mat') ;
        if exist(fn, 'file')
            load(fn) % Load the surface area and volume from disk
            
            P{mi,j} = (aas)./(vvs.^(2/3)) ; % dedimensionalized value
            
            if mi == 1 && j == 1
                P{mi,j} = P{mi,j}(30:110);
            else
                P{mi,j} = vertcat((P{mi,j}(2)+zeros(length(P{1,1})-length(P{mi,j}),1)),P{mi,j},(P{mi,j}(end-2)+zeros(1,1)));
            end
            time = zeros(length(P{mi,j}),1);
            for i = 1:length(time)
                time(i) = i;
            end
            
            P1{mi,j} = horzcat(time, P{mi,j});
            %options = optimset('PlotFcns',@optimplotfval);
            %toffmin{mi,j} = fminsearch(@(vars)rotTransMesh(vars, P1{mi,j}, P1{1,1}), [0, 1, 0, 1]) ;%, options)
            
            % get col offset
            if strcmp(label, 'Membrane')
                col_off = 0 ;
            elseif strcmp(label, 'Nuclei')
                col_off = 1 ;
            elseif strcmp(label, 'Actin')
                col_off = 3 ;
            end
            
            meshes = dir(fullfile(matdir, 'aligned_meshes',  'mesh_apical_stab_0*_APDV_um.ply')) ;
            
            % Write each mesh to the cell array mca
            for k = 1:length(meshes)
                %toff = round(toffmin{mi,j}(1));
                mca{k, mi-1+j+col_off} = meshes(k);
            end
        else
            disp(['Could not find ' fn])
        end
    end
end

%% Create Output folders

outdir = '/mnt/data/analysis/' ;
cd(outdir)

if ~isfile([outdir, 'aligned_meshes_xy_2D'])
    mkdir /mnt/data/analysis/aligned_meshes_xy_2D/ ;
end

if ~isfile([outdir, 'aligned_meshes_xy_3D'])
    mkdir /mnt/data/analysis/aligned_meshes_xy_3D/ ;
end

if ~isfile([outdir, 'aligned_meshes_yz_2D'])
    mkdir /mnt/data/analysis/aligned_meshes_yz_2D/ ;
end

if ~isfile([outdir, 'aligned_meshes_yz_3D'])
    mkdir /mnt/data/analysis/aligned_meshes_yz_3D/ ;
end

leg = {};

%% Specify Plot Type %%

dimsfa = {} ; % array of scaling factors in each dimension
dimsfa{1,1} = 1 ; %scaling factor for x of CAAX excelent is 1
dimsfa{2,1} = 1 ; %scaling factor for y of CAAX excelent is 1
dimsfa{3,1} = 1 ; %scaling factor for z of CAAX excelent is 1

for c = 1:6
    r = 1;
    scaleFound = false;
    while ~scaleFound
        if ~isfield(mca{r,c}, 'name') == 1 % If TP doesn't exist
            r = r+1; % go on to next row and repeat
        else
            if c == 1
                caaxmesh = ply_read_with_normals(fullfile(mca{r,1}.folder, mca{r,1}.name));
                caax_xrange = max(caaxmesh.vertex.x(:))-min(caaxmesh.vertex.x(:));
                caax_yrange = max(caaxmesh.vertex.y(:))-min(caaxmesh.vertex.y(:));
                caax_zrange = max(caaxmesh.vertex.z(:))-min(caaxmesh.vertex.z(:));
            else
                fmesh = ply_read_with_normals(fullfile(mca{r,c}.folder, mca{r,c}.name));
                xrange = max(fmesh.vertex.x(:))-min(fmesh.vertex.x(:));
                xscfa = caax_xrange/xrange;
                dimsfa{1,c} = xscfa;
                yrange = max(fmesh.vertex.y(:))-min(fmesh.vertex.y(:));
                yscfa = caax_yrange/yrange;
                dimsfa{2,c} = yscfa;
                zrange = max(fmesh.vertex.z(:))-min(fmesh.vertex.z(:));
                zscfa = caax_zrange/zrange;
                dimsfa{3,c} = zscfa;
            end
            scaleFound = true;
        end
    end
end

%% Scale and Plot the meshes

color = {'red', 'k', 'y', 'g', 'b', 'c'} ;
labels = {'CAAX 201902072000', 'CAAX 201903211930', 'RFP 201901021550', 'RFP 201904031830', 'RFP 201903312000', 'LifeAct 201904021800'};

for r = 10:10:100 %length(mca) % Iterate over rows of mca (corresponding to timepoints)
    fig1 = figure('Name', 'Lateral View', 'Visible', 'Off'); hold on
    for c = 1:6 %Iterate over columns of mca (corresponding to experiments)
        if isfield(mca{r,c}, 'name') == 1 % If that TP exists
            cmesh = ply_read_with_normals(fullfile(mca{r,c}.folder, mca{r,c}.name)); % load mesh
            axis equal
            xoff = min(cmesh.vertex.x(:)) ; % y offset
            yoff = mean(cmesh.vertex.y(:)) ; % y offset
            zoff = mean(cmesh.vertex.z(:)) ; % z offset
            if strncmp(dim, 'two', 3)
                epsz = .5 ; % epsilon in z direction
                epsx = .5 ; % epsilon in x direction
                if strncmp(orient, 'lateral', 6)
                    x_zslice = dimsfa{1,c}*(cmesh.vertex.x(cmesh.vertex.z < epsz & cmesh.vertex.z > -epsz)-xoff) ;
                    y_zslice = dimsfa{2,c}*(cmesh.vertex.y(cmesh.vertex.z < epsz & cmesh.vertex.z > -epsz)-yoff) ;
                    scatter(x_zslice, y_zslice, '.', 'MarkerFaceColor', color{c});
                else
                    y_xslice = dimsfa{1,c}*(cmesh.vertex.y(cmesh.vertex.x < (slice_AP+epsx) & cmesh.vertex.x > (slice_AP-epsx))-yoff) ;
                    z_xslice = dimsfa{1,c}*(cmesh.vertex.z(cmesh.vertex.x < (slice_AP+epsx) & cmesh.vertex.x > (slice_AP-epsx))-zoff) ;
                    scatter(y_xslice, z_xslice, '.', 'MarkerFaceColor', color{c});
                end
            else
                x = dimsfa{1,c}*(cmesh.vertex.x(:)-xoff) ; % scale and offset x
                y = dimsfa{2,c}*(cmesh.vertex.y(:)-yoff) ; % scale and offset y
                z = dimsfa{3,c}*(cmesh.vertex.z(:)-zoff) ; % scale and offset z
                trisurf(cell2mat(cmesh.face.vertex_indices)+1, x, y, z,...
                    ones(size(cmesh.vertex.z(:))), 'Facecolor', color{c}, 'FaceAlpha', .1, 'edgecolor', 'none'); % plot mesh
            end
            leg{length(leg)+1} = labels{c};
        else
        end
    end
    
    % Plot Meshes
    
    ax.XDir = 'normal' ;
    ax.YDir = 'normal'  ;
    ax.ZDir = 'normal';
    
    xlabel('AP [\mum]');
    ylabel('Lateral [\mum]');
    zlabel('DV [\mum]');
    lgd = legend(leg);
    lgd.FontSize=5;
    
    if strncmp(dim, 'two', 3)
        if strncmp(orient, 'lateral', 6)
            view(2)
            title('Lateral View');
            saveas(gcf, [outdir, 'aligned_meshes_xy_2D', '/aligned_meshes_' sprintf('%03d', r), '_2D.png']);
            fprintf('Saving: %s \n', [outdir, 'aligned_meshes_xy_2D', '/aligned_meshes_' sprintf('%03d', r), '_2D.png']);
        else
            title('Posterior View');
            saveas(gcf, [outdir, 'aligned_meshes_yz_2D', '/aligned_meshes_' sprintf('%03d', r), '_2D.png']);
            fprintf('Saving: %s \n', [outdir, 'aligned_meshes_yz_2D', '/aligned_meshes_' sprintf('%03d', r), '_2D.png']);
        end
    else
        view(2) % xy view
        title('Lateral View')
        saveas(gcf, [outdir, 'aligned_meshes_xy_3D', '/aligned_meshes_' sprintf('%03d', r), '.png']);
        fprintf('Saving: %s \n', [outdir, 'aligned_meshes_xy_3D', '/aligned_meshes_' sprintf('%03d', r), '.png']);
        
        view(90, 180) % yz view
        title('Posterior View')
        saveas(gcf, [outdir, 'aligned_meshes_yz_3D', '/aligned_meshes_' sprintf('%03d', r), '.png'])
        fprintf('Saving: %s \n', [outdir, 'aligned_meshes_yz_3D', '/aligned_meshes_' sprintf('%03d', r), '.png'])
    end
    leg = {};
end