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
preview = true ;
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
% xmin = xyzlim(1) * resolution * ssfactor ;
% ymin = xyzlim(2) * resolution * ssfactor ;
% zmin = xyzlim(3) * resolution * ssfactor ;
% xmax = xyzlim(4) * resolution * ssfactor ;
% ymax = xyzlim(5) * resolution * ssfactor ;
% zmax = xyzlim(6) * resolution * ssfactor ;

%% Get xyzlimits of the centerline
xsmin = 0 ;
ysmin = 0 ;
zsmin = 0 ;
xsmax = 0 ;
ysmax = 0 ;
zsmax = 0 ;
smax = 0 ;
for ii=1:length(fns)
    %% Name the output centerline
    name_split = strsplit(fns(ii).name, '.ply') ;
    name = name_split{1} ; 
    expstr = strrep(num2str(exponent, '%0.1f'), '.', 'p') ;
    skel_rs_outfn = [fullfile(outdir, name) '_centerline_scaled_exp' expstr] ;
    skel_outfn = [fullfile(outdir, name) '_centerline_exp' expstr ] ;
    tmp = strsplit(name, '_') ;
    timestr = tmp{length(tmp)} ;
    
    % Load centerline 
    disp(['Loading centerline from txt: ', skel_rs_outfn, '.txt'])
    try
        sskelrs = importdata([skel_rs_outfn '.txt']) ;
        ss = sskelrs(:, 1) ;
        ds = diff(ss) ; 
        ds = [ds; ds(length(ds))] ;
        skelrs = sskelrs(:, 2:4) ;
    catch
        % Save the rotated, translated, scaled curve
        skel = importdata([skel_outfn '.txt']) ;
        skelr = (rot * skel')' + trans ; 
        
        % get distance increment
        ds = vecnorm(diff(skel), 2, 2) ;
        % get pathlength at each skeleton point
        ss = [0; cumsum(ds)] ;
        
        % Create rotated and scaled skeleton
        ss_s = ss * resolution * ssfactor ;
        skelr_s = skelr * resolution * ssfactor ;
        disp(['Saving rotated & scaled skeleton to txt: ', skel_rs_outfn, '.txt'])
        dlmwrite([skel_rs_outfn '.txt'], [ss_s, skelr_s])
        
        % Recompute ds        
        sskelrs = importdata([skel_rs_outfn '.txt']) ;
        ss = sskelrs(:, 1) ;
        ds = diff(ss) ; 
        ds = [ds; ds(length(ds))] ;
        skelrs = sskelrs(:, 2:4) ;
    end
    xsmin = min(xsmin, min(skelrs(:, 1))) ;
    ysmin = min(ysmin, min(skelrs(:, 2))) ;
    zsmin = min(zsmin, min(skelrs(:, 3))) ;
    xsmax = max(xsmax, min(skelrs(:, 1))) ;
    ysmax = max(ysmax, min(skelrs(:, 2))) ;
    zsmax = max(zsmax, min(skelrs(:, 3))) ;
    smax = max(max(ss), smax) ;

end


%% Iterate through each mesh
for ii=1:length(fns)
    %% Name the output centerline
    name_split = strsplit(fns(ii).name, '.ply') ;
    name = name_split{1} ; 
    expstr = strrep(num2str(exponent, '%0.1f'), '.', 'p') ;
    skel_rs_outfn = [fullfile(outdir, name) '_centerline_scaled_exp' expstr] ;
    skel_outfn = [fullfile(outdir, name) '_centerline_exp' expstr ] ;
    figsmoutname = [fullfile(fig_smoothoutdir, name) '_centerline_smoothed_exp' expstr] ;
    tmp = strsplit(name, '_') ;
    timestr = tmp{length(tmp)} ;
    
    % Load centerline 
    disp(['Loading centerline from txt: ', skel_rs_outfn, '.txt'])
    sskelrs = importdata([skel_rs_outfn '.txt']) ;
    ss = sskelrs(:, 1) ;
    ds = diff(ss) ; 
    ds = [ds; ds(length(ds))] ;
    skelrs = sskelrs(:, 2:4) ;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute writhe of the centerline
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filter the result
    % Smoothing parameters
    framelen = 17 ;  % must be odd
    polyorder = 2 ;
    xskel = skelrs(:, 1) ;
    yskel = skelrs(:, 2) ;
    zskel = skelrs(:, 3) ;
    
    % Note we ignore the variations in ds to do this fit
    scx = savgol(xskel, polyorder, framelen)' ;
    scy = savgol(yskel, polyorder, framelen)' ;
    scz = savgol(zskel, polyorder, framelen)' ;
    sc = [scx, scy, scz] ;
    
    % Fit the smoothed curve
    xcoeffs = polyfit(ss, scx, 7) ;
    ycoeffs = polyfit(ss, scy, 7) ;
    zcoeffs = polyfit(ss, scz, 7) ;
    ssx = linspace(min(ss), max(ss), 100) ;
    dsx = gradient(ssx) ;
    xp = polyval(xcoeffs, ssx);
    yp = polyval(ycoeffs, ssx);
    zp = polyval(zcoeffs, ssx);
    
    %% Check the smoothing
    if preview
        close all
        fig = figure;
        set(gcf, 'Visible', 'Off')
        hold on
        plot3(xskel, yskel, zskel, '-.')
        hold on 
        plot3(scx, scy, scz, 's')
        title('Smoothing')
        % save the figure
        xlim([xsmin xsmax])
        ylim([ysmin ysmax])
        zlim([zsmin zsmax])
        axis equal
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);
        saveas(fig, [figsmoutname '.png'])
        if preview
            set(gcf, 'Visible', 'On')
        end
        
        % Save the difference
        close all
        fig = figure;
        set(fig, 'Visible', 'Off')
        plot(ss, scx - xskel)
        hold on
        plot(ss, scy - yskel)
        plot(ss, scz - zskel)
        title('\Delta = smoothing - original curve')
        ylabel(['\Delta [\mu' 'm]'])
        xlabel(['s [\mu' 'm]'])
        set(fig, 'PaperUnits', 'centimeters');
        set(fig, 'PaperPosition', [0 0 xwidth ywidth]);
        saveas(fig, [figsmoutname '_diff.png'])
        
        % Save the polynomial fit
        close all
        fig = figure;
        set(fig, 'Visible', 'Off')
        plot3(xskel, yskel, zskel, '-.')
        hold on 
        plot3(xp, yp, zp, 's')
        title('Polynomial fit to the centerline')
        xlabel('x')
        ylabel('y')
        % save the figure
        xlim([xsmin xsmax])
        ylim([ysmin ysmax])
        zlim([zsmin zsmax])
        axis equal
        set(fig, 'PaperUnits', 'centimeters');
        set(fig, 'PaperPosition', [0 0 xwidth ywidth]);
        saveas(fig, [figsmoutname '_fit.png'])
    end
    
    %% Get rate of change of curve along curve
    gradc_raw = [gradient(xp'), gradient(yp'), gradient(zp')] ; 
    gradc = bsxfun(@rdivide, gradc_raw, dsx(:)) ;
    gradc_ds = vecnorm(gradc, 2, 2) ;
    % Compute the tangent to the curve
    tangent = bsxfun(@rdivide, gradc, gradc_ds(:)) ;

    % Compute normal 
    normal_raw = [gradient(tangent(:, 1)), ...
        gradient(tangent(:, 2)), ...
        gradient(tangent(:, 3))] ; 
    normalc = bsxfun(@rdivide, normal_raw, dsx(:)) ;
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
    chirality = dphi ./ ssx(1:end-1)' ;

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
        close all
        fig = figure ;
        set(fig, 'Visible', 'Off')
        plot(ssx(1:end-1), dphi)
        title('$d\phi/ds$', 'Interpreter', 'Latex')
        ylabel('$d\phi/ds$', 'Interpreter', 'Latex')
        xlabel('$s$ [$\mu$m]', 'Interpreter', 'Latex')
        saveas(fig, [figsmoutname '_dphids.png'])
    end
        
    % Save the image of chirality
    close all
    fig = figure;
    set(gcf, 'Visible', 'Off')
    plot(ssx(1:end-1)', chirality)
    % plot(sss, chirs, '-')
    title('chirality of centerline')
    ylabel('$d \phi / d s$', 'Interpreter', 'Latex')
    xlabel('path length, $s$', 'Interpreter', 'Latex')
    xlim([0, smax])
    ylim([-0.06, 0.06])
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
    dlmwrite(datfn, [ssx(1:end-1)' chirality])
    
    % Compute the writhe
    % Wr = 1/4pi \int \int T(s) xx T(s') \cdot [R(s) - R(s') / |R(s) - R(s')|^3]
    % first compute vec from each point to every other point
    for jj=1:length(ssx)
        
    end
    
    
    
    
end

disp('done')
