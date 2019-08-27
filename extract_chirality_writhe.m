%% Extract the centerlines from a series of meshes (PLY files)
% Noah Mitchell 2019
% This version relies on Gabriel Peyre's toolbox called
% toolbox_fast_marching/
% Run from the msls_output directory
clear ;

%% First, compile required c code
% mex ./FastMarching_version3b/shortestpath/rk4
close all ;
odir = pwd ;
codepath = '/mnt/data/code/gut_matlab/' ;
if ~exist(codepath, 'dir')
    codepath = [pwd filesep] ;
end
addpath(codepath)
addpath([codepath 'addpath_recurse' filesep]) ;
addpath([codepath 'mesh_handling' filesep]);
addpath([codepath 'inpolyhedron' filesep]);
addpath([codepath 'savgol' filesep])
addpath_recurse('/mnt/data/code/gptoolbox/')
% addpath_recurse([codepath 'gptoolbox' filesep])

toolbox_path = [codepath 'toolbox_fast_marching/toolbox_fast_marching/'];
dtpath = [codepath 'distanceTransform/'] ;
addpath_recurse(toolbox_path)
addpath(dtpath)
% compile_c_files
cd(odir)

%% Parameters
res = 1 ;
resolution = 0.2619 ;
buffer = 5 ;
plot_buffer = 30;
ssfactor = 4; 
weight = 0.1;
normal_step = 0.5 ; 
preview = false ;
eps = 0.01 ;
meshorder = 'zyx' ;
exponent = 1;
% figure parameters
xwidth = 16 ; % cm
ywidth = 10 ; % cm

% Find all meshes to consider
meshdir = pwd ;
cd ../
rootdir = pwd ;
cd(meshdir)
% rootpath = '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/' ;
% rootpath = [rootpath 'Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/'] ;
% if ~exist(rootpath, 'dir')
%     rootpath = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/' ;
%     rootpath = [rootpath 'data/48Ygal4UasCAAXmCherry/201902072000_excellent/'] ;
% end
% meshdir = [rootpath 'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/'];
fns = dir(fullfile(meshdir, 'mesh_apical_stab_0*.ply')) ;
rotname = fullfile(meshdir, 'rotation_APDV') ;
transname = fullfile(meshdir, 'translation_APDV') ;
xyzlimname = fullfile(meshdir, 'xyzlim_APDV') ;

% Name output directory
outdir = [fullfile(fns(1).folder, 'centerline') filesep ];
if ~exist(outdir, 'dir')
    mkdir(outdir) ;
end
% figure 1
choutdir = [outdir 'chirality' filesep];
if ~exist(choutdir, 'dir')
    mkdir(choutdir) ;
end
figoutdir = [choutdir 'images' filesep];
if ~exist(figoutdir, 'dir')
    mkdir(figoutdir) ;
end

fig_smoothoutdir = [outdir 'images' filesep 'smoothed' filesep];
if ~exist(fig_smoothoutdir, 'dir')
    mkdir(fig_smoothoutdir) ;
end

ii = 1 ;

%% Load transformations
% Load the rotation matrix
rot = importdata([rotname '.txt']) ;
% Load the translation to put anterior to origin
trans = importdata([transname '.txt']) ;
% Load plotting limits
xyzlim = importdata([xyzlimname '.txt']) ;
xmin = xyzlim(1) * resolution * ssfactor ;
ymin = xyzlim(2) * resolution * ssfactor ;
zmin = xyzlim(3) * resolution * ssfactor ;
xmax = xyzlim(4) * resolution * ssfactor ;
ymax = xyzlim(5) * resolution * ssfactor ;
zmax = xyzlim(6) * resolution * ssfactor ;

%% Iterate through each mesh
for ii=1:length(fns)
    %% Name the output centerline
    name_split = strsplit(fns(ii).name, '.ply') ;
    name = name_split{1} ; 
    expstr = strrep(num2str(exponent, '%0.1f'), '.', 'p') ;
    centerlinename = [fullfile(outdir, name) '_centerline_scaled_exp' expstr] ;
    figsmoutname = [fullfile(fig_smoothoutdir, name) '_centerline_smoothed_exp' expstr] ;
    tmp = strsplit(name, '_') ;
    timestr = tmp{length(tmp)} ;
    
    % Load centerline 
    disp(['Loading centerline from txt: ', centerlinename, '.txt'])
    try
        sskelrs = importdata([centerlinename '.txt']) ;
        ss = sskelrs(:, 1) ;
        ds = diff(ss) ; 
        ds = [ds; ds(length(ds))] ;
        skelrs = sskelrs(:, 2:4) ;
    except
        % Save the rotated, translated, scaled curve
        ss_s = ss * resolution * ssfactor ;
        skelr_s = skelr * resolution * ssfactor ;
        disp(['Saving rotated & scaled skeleton to txt: ', skel_rs_outfn, '.txt'])
        dlmwrite([skel_rs_outfn '.txt'], [ss_s, skelr_s])
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute writhe of the centerline
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filter the result
    % Smoothing parameters
    framelen = 37 ;  % must be odd
    polyorder = 2 ;
    xskel = skelrs(:, 1) ;
    yskel = skelrs(:, 2) ;
    zskel = skelrs(:, 3) ;
    
    % Note we ignore the variations in ds to do this fit
    scx = savgol(xskel, polyorder, framelen)' ;
    scy = savgol(yskel, polyorder, framelen)' ;
    scz = savgol(zskel, polyorder, framelen)' ;
    sc = [scx, scy, scz] ;
    
    %% Check the smoothing
    if preview
        close all
        fig = figure;
        set(gcf, 'Visible', 'Off')
        tmp = trisurf(tri + 1, xyzr(:, 1), xyzr(:,2), xyzr(:, 3), ...
            xyzr(:, 1), 'edgecolor', 'none', 'FaceAlpha', 0.1) ;
        hold on
        plot3(xskel, yskel, zskel, '-.')
        hold on 
        plot3(scx, scy, scz, 's')
        title('Smoothing')
        % save the figure
        xlim([xmin xmax])
        ylim([ymin ymax])
        zlim([zmin zmax])
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);
        saveas(fig, [figsmoutname '.png'])
        if preview
            set(gcf, 'Visible', 'On')
        end
        
        close all
        plot(ss, scx - xskel)
        hold on
        plot(ss, scy - yskel)
        plot(ss, scz - zskel)
        title('\Delta = smoothing - original curve')
        ylabel(['\Delta [\mu' 'm]'])
        xlabel(['s [\mu' 'm]'])
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);
        saveas(fig, [figsmoutname '_diff.png'])
    end
    
    %% Get rate of change of curve along curve
    gradc_raw = [gradient(sc(:, 1)), gradient(sc(:, 2)), gradient(sc(:, 3))] ; 
    gradc = bsxfun(@rdivide, gradc_raw, ds(:)) ;
    gradc_ds = vecnorm(gradc, 2, 2) ;
    % Compute the tangent to the curve
    tangent = bsxfun(@rdivide, gradc, gradc_ds(:)) ;
    
%     % Fit the result
%     xcoeffs = polyfit(ss, scx, 7) ;
%     ycoeffs = polyfit(ss, scy, 7) ;
%     zcoeffs = polyfit(ss, scz, 7) ;
%     xp = poly1d(xcoeffs) ;
%     yp = poly1d(ycoeffs) ;
%     zp = poly1d(zcoeffs) ;
% 
%     % Re-parametrize the curve with equidistant ss
%     ss = linspace(min(ss), max(ss), 100)
%     sc = [xp(ss), yp(ss), zp(ss)]
%     % Re-compute rate of change of curve along curve
%     gradc = grad_curve(sc)
%     ds = sqrt(sum(gradc ^ 2, 2))
%     % Compute the tangent to the curve
%     tangent = bsxfun(@rdivide, gradc, gradc_ds(:)) ;

    % Compute normal 
    normal_raw = [gradient(tangent(:, 1)), ...
        gradient(tangent(:, 2)), ...
        gradient(tangent(:, 3))] ; 
    normalc = bsxfun(@rdivide, normal_raw, ds(:)) ;
    normalc_ds = vecnorm(normalc, 2, 2) ;
    normal = bsxfun(@rdivide, normalc, normalc_ds(:)) ;
    
    % Compute binormal from cross product
    binormalc = cross(tangent, normal) ;
    binormalc_ds = vecnorm(binormalc, 2, 2) ;
    binormal = bsxfun(@rdivide, binormalc, binormalc_ds(:)) ;
    
    % Now compute how the normal changes along the AP axis    
    v1 = normal(1:end-1, :) ;
    v2 = normal(2:end, :) ;
    dphi = acos(sum(v1 .* v2, 2)) ;
    signphi = sign(sum(tangent(1:end-1, :) .* cross(v1, v2), 2)) ;
    dphi = dphi .* signphi ;
    chirality = dphi ./ ss(1:end-1) ;

    % Clean up measurement
    % chirs = savgol_filter(chirality, window_length=351, polyorder=polyorder)
    % try:
    %     chirs = moving_average(chirs, n=avg_window)
    %     sss = ss[int(np.floor(avg_window * 0.5)):]
    %     sss = sss[:-int(np.floor(avg_window * 0.5))]
    % except ValueError:
    %     sss = ss
    
    % coeffs = np.polyfit(ss, chirality, 5)
    % polynomial = np.poly1d(coeffs)
    % sss = np.linspace(min(ss), max(ss), 100)
    % chirs = polynomial(sss)

    % Check the angle change along AP axis
    if preview
        plt.plot(np.arange(len(dphi)), dphi)
        plt.plot(np.arange(len(ds)), ds)
        plt.ylabel('$d\phi$, $ds$', 'Interpreter', 'Latex')
        plt.xlabel('index')
        plt.show()
    end
        
    % Save the image of chirality
    close all
    fig = figure;
    set(gcf, 'Visible', 'Off')
    plot(ss(1:end-1), chirality)
    % plot(sss, chirs, '-')
    title('chirality of centerline')
    ylabel('$d \phi / d s$', 'Interpreter', 'Latex')
    xlabel('path length, $s$', 'Interpreter', 'Latex')
    xlim([xmin, xmax])
    ylim([ymin, ymax])
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=16cm
    saveas(fig, fullfile(figoutdir, [name '_chirality.png']))
    close all

    % Save the data
    datfn = fullfile(choutdir, [name '_chirality.txt' ]);
    header = 'Chirality of centerline for timepoint ' ;
    header = [header name ': s, dphi/ds'] ;
    fid = fopen(datfn, 'wt');
    fprintf(fid, header); 
    fclose(fid);
    dlmwrite(datfn, [ss(1:end-1) chirality])
    
end

disp('done')
