% figure_tubular_generate_synthetic_dataset

close all 
clear all

%% Add paths
dataDir = '/mnt/data/tubular_test/single_coil_dynamic3/' ;
mkdir(dataDir)
% dataDir = '/mnt/data/tubular_test/tube_centerlines/single_coil_path_independence/' ;
origpath = '/mnt/data/code/tubular/' ;
cd(origpath)
addpath(origpath) ;
addpath(genpath('/mnt/data/code/gptoolbox'))
addpath(genpath(fullfile(origpath, 'utility')))
addpath('TexturePatch')
addpath('DECLab')
addpath('RicciFlow_MATLAB')
addpath(fullfile('utility','plotting'))
addpath(fullfile('utility','plotting'))
% go back to the data
cd(dataDir)


%%
outdir = dataDir ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the meshes and data
overwrite = true ;
TT = 20 ;
for tp = 1
    
    meshfn = fullfile(outdir, 'mesh_output', sprintf('mesh_ls_%06d.ply', tp)) ;
    meshsmfn = fullfile(outdir, 'mesh_output', sprintf('mesh_%06d.ply', tp)) ;
    
    if ~exist(meshfn, 'file')  || overwrite
        RR = (1 + 1.5*sin(((tp-1)/TT-0.25)*2*pi)); % radius of coil (width of lateral motion)
        
        % if tp == 5
        %     R1 = 0.35 + 0.3*sin(((tp-1)/TT-0.25)*2*pi); % amount of radial undulation / variation along pathlength
        % else
        R1 = 0.15 + 0.02*sin(((tp-1)/TT-0.25)*2*pi); % amount of radial undulation / variation along pathlength
        % end
        debug = false ;

        % make spiral tube
        LL = 10 ;
        zz = linspace(0, LL, 100)  ;
        % width from each point along z to consider inside mesh
        ww = max(0.02, 0.25 + R1 * (cos(0.5*(zz/LL+0.5)*2*pi))+0.1) ;
        
        % if tp < 5
        %     ww = ww + (5-tp) * 0.02 * (1 + sin(((zz-min(zz))/LL-0.25)*2*pi)) ;
        % end
        
        % if tp == 5
        %     thres = 0.25 ;
        %     ww(ww < thres) = ww(ww < thres) - 2 * (thres - ww(ww < thres)) ;
        % end
                
        plot(zz, ww)
        ylim([0,1])
        pause(2)
        noiseX = rand(size(zz)) * max(0,0.5-RR)  ;
        noiseY = rand(size(zz)) * max(0,0.5-RR) ;
        
        noiseX = movmean(noiseX, 15) ;
        noiseY = movmean(noiseY, 15) ;
        
        % points in 3D of the coil
        xcurv0 = RR * cos(zz*2*pi/LL) .* hanning(numel(zz))' + noiseX ; 
        ycurv0 = RR * sin(zz*2*pi/LL) .* hanning(numel(zz))' + noiseY ;
        zcurv0 = LL * 0.5 * ((2*zz/LL-1).^3 + 1)  + 0.5*RR*(sin(zz*2*pi/LL)) ;
        ss = cumsum([0, vecnorm(diff([xcurv0; ycurv0; zcurv0]', 1), 2, 2)']) ;
        xcurv = equidistantSampling1D(xcurv0, ss, numel(zz)) ;
        ycurv = equidistantSampling1D(ycurv0, ss, numel(zz)) ;
        
        if tp == 5
            ycurv = (ycurv)*2 ;
            % ycurv(1:49) = ycurv(1:49) + 0.01*(1:49) ;
            % ycurv(51:100) = ycurv(51:100) - 0.01*(1:50) ;
        end
        
        zcurv = equidistantSampling1D(zcurv0, ss, numel(zz)) ;
        % zcurv = mean(zz, zcurv) ;
        rr = [ xcurv', ycurv' ...
           zcurv' ] ;
        
        plot3(rr(:, 1), rr(:, 2), rr(:, 3), '.-')
        axis equal
        
        pause(0.1)
%     end
% end
% %
% try
%     if true
        % Build 3D volume segmentation
        dx = 0.04 ;
        xlin = -3:dx:3 ;
        ylin = -3:dx:3 ;
        zlin = -2:dx:LL+2 ;
        [segx, segy, segz] = meshgrid(xlin, ylin, zlin) ;
        seg = zeros(size(segx)) ;

        % segment the volume
        disp('segmenting area near curve')
        if tp < TT*0.5+0.01
            pts2do = tp:(length(zz)-tp) ;
        else
            pts2do = max((TT-tp), 1):length(zz)-max((TT-tp), 1) ;
        end
        for ii = pts2do
            if mod(ii, 10) == 0 
                disp(['ii = ' num2str(ii) ' / ' num2str(length(zz))])
            end
            % include voxels near this point ii
            xpt = rr(ii, 1) ;
            ypt = rr(ii, 2) ;
            zpt = rr(ii, 3) ;
            seg = seg + ((segx-xpt).^2 + (segy-ypt).^2 + (segz-zpt).^2 < ww(ii)) ;
        end
        seg = seg>0 ;

        % Inspect seg
        if debug
            for ii = 1:5:max(size(seg))
                if ii<= size(seg, 1) 
                    subplot(1, 3, 1)
                    imagesc(squeeze(seg(ii, :, :)))
                    axis equal
                end
                if ii<= size(seg, 2) 
                    subplot(1, 3, 2)
                    imagesc(squeeze(seg( :, ii, :)))
                    axis equal
                end
                if ii <= size(seg, 3)
                    subplot(1, 3, 3)
                    imagesc(squeeze(seg( :, :, ii)))
                    axis equal
                end
                pause(1e-9)
            end
        end

        % Make mesh 
        disp('making mesh from isosurface')
        mesh = isosurface(seg) ;
        mesh.v = mesh.vertices ;
        mesh.f = mesh.faces ;

        mesh = rmfield(mesh, 'vertices') ;
        mesh = rmfield(mesh, 'faces') ;
        
               
        % SAVE result
        plywrite_with_normals(meshfn, mesh.f, mesh.v) 
        save(fullfile(outdir, 'segs', sprintf('seg_%06d.mat', tp)), 'seg')

        % Plot result
        [colors, names] = define_colors ;
        sky = colors(6, :) ;
        disp('plotting resulting mesh')
        clf
        trisurf(triangulation(mesh.f, mesh.v(:, [3,1,2])), 'edgecolor', 'none',...
            'facecolor', sky, 'facealpha', 0.5)
        axis equal
        lighting gouraud
        material dull
        camlight 
        axis off
        set(gcf, 'color', 'white')
        pause(0.01)
        export_fig(gcf, fullfile(outdir, ...
            sprintf('mesh_and_skeleton_raw_%06d.png', tp)), '-nocrop', '-r100')
    else
        disp('already on disk, loading...')
        mesh = read_ply_mod(meshfn) ;
        load(fullfile(outdir, 'segs', sprintf('seg_%06d.mat', tp)), 'seg')
    end

end

%%
for tp = 1:TT
    
    meshfn = fullfile(outdir, 'mesh_output', sprintf('mesh_ls_%06d.ply', tp)) ;
    meshsmfn = fullfile(outdir, 'mesh_output', sprintf('mesh_%06d.ply', tp)) ;
    
    mesh = read_ply_mod(meshfn) ;
    load(fullfile(outdir, 'segs', sprintf('seg_%06d.mat', tp)), 'seg')

    if ~exist(meshsmfn, 'file') || overwrite
        %
        disp(['tp = ' num2str(tp) ': smoothing vertices'])
        % targetEdgeLen = 0.1 ;
        % numIterations = 100 ;
        % [ff, vv, fn, vn] = isotropic_remeshing(mesh.f, mesh.v, targetEdgeLen, numIterations) ;
        % mesh_sm = reducepatch(mesh,0.4) ;
        % mesh_sm.v = mesh_sm.vertices ;
        % mesh_sm.f = mesh_sm.faces ;
        % mesh_sm = rmfield(mesh_sm, 'vertices') ;
        % mesh_sm = rmfield(mesh_sm, 'faces') ;

        mesh_sm = mesh ;
        mesh_sm.v = laplacian_smooth(mesh.v, mesh.f, 'cotan', ...
            [], 0.01, 'implicit', mesh.v) ;

        disp('reorienting_facets')
        mesh_sm.f = reorient_facets(mesh_sm.v, mesh_sm.f) ;
        mesh_sm.vn = per_vertex_normals(mesh_sm.v, mesh_sm.f) ;
        % mesh_sm = struct() ;
        % mesh_sm.f = ff ;
        % mesh_sm.v = vv ;
        % mesh_sm.vn = vn ; 
        % mesh_sm.fn = fn ;

        % Compute skeleton
        disp('Computing skeleton')
        [sy, sx, sz] = ind2sub(size(seg), find(bwskel(seg))) ;
        skel = [sx, sy, sz] ;
        
        plywrite_with_normals(meshsmfn, mesh_sm.f, mesh_sm.v, mesh_sm.vn) 
        save(fullfile(outdir, ...
            sprintf('bwskel_centerline_%06d.mat', tp)), 'skel', 'mesh', 'mesh_sm', 'seg')
        save(fullfile(outdir, ...
            sprintf('input_centerline_%06d.mat', tp)), 'xlin', 'ylin', 'zlin', 'rr') ;

        % Plot result
        [colors, names] = define_colors ;
        sky = colors(6, :) ;
        disp('plotting resulting mesh')
        clf
        trisurf(triangulation(mesh_sm.f, mesh_sm.v(:, [3,1,2])), 'edgecolor', 'none',...
            'facecolor', sky, 'facealpha', 0.5)
        axis equal
        hold on;
        plot3(sz, sx, sy, 'k-', 'linewidth', 2)
        lighting gouraud
        material dull
        camlight 
        axis off
        set(gcf, 'color', 'white')
        export_fig(gcf, fullfile(outdir, ...
            sprintf('mesh_and_skeleton_%06d.png', tp)), '-nocrop', '-r600')
        
        % Save segmentation as volume
        writeTiff5D(reshape(uint16(seg*65535), ...
            [size(seg, 1), size(seg, 2), 1, size(seg, 3)]), ...
            sprintf('Time_%06d_segmentation.tif', tp), 16)

        % Save intensity data as volume
        intens = cos(segx*10) .* sin(segy*10) + 1 ; 
        % 
        writeTiff5D(reshape(uint16(intens*65535*0.5), ...
             [size(seg, 1), size(seg, 2), 1, size(seg, 3)]), ...
             sprintf('Time_%06d.tif', tp), 16)
    end
end


%% Now run the TubULAR pipeline


if ~exist(fullfile(dataDir, 'xp.mat'), 'file') 
    %% DEFINE NEW MASTER SETTINGS
    if ~exist('./masterSettings.mat', 'file') 
        % Metadata about the experiment
        stackResolution = [1, 1, 1] ;  % resolution in spaceUnits per pixel
        nChannels = 2 ;             % how many channels is the data (ex 2 for GFP + RFP)
        channelsUsed = [1 2] ;          % which channels are used for analysis
        timePoints = 1:15;       % timepoints to include in the analysis
        ssfactor = 4 ;              % subsampling factor
        flipy = false ;             % whether the data is stored inverted relative to real position in lab frame
        timeInterval = 1 ;          % physical interval between timepoints
        timeUnits = 'min' ;         % physical unit of time between timepoints
        spaceUnits = [char(956) 'm'] ;     % physical unit of time between timepoints
        fn = 'Time_%06d';        % filename string pattern
        set_preilastikaxisorder = 'xyzc' ; % data axis order for subsampled h5 data (ilastik input)
        swapZT = 0 ;                % whether to swap the z and t dimensions
        masterSettings = struct('stackResolution', stackResolution, ...
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

    if loadMaster
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
        fn = masterSettings.fn ;
        set_preilastikaxisorder = masterSettings.set_preilastikaxisorder ;
        swapZT = masterSettings.swapZT ;
        nU = masterSettings.nU ;
        nV = masterSettings.nV ;
    end
    dir16bit = fullfile(dataDir) ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PART 1: Define the metadata for the project
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd(dir16bit)
    dataDir = cd ;
    projectDir = dataDir ;

    % A filename base template - to be used throughout this script
    fileMeta                    = struct();
    fileMeta.dataDir            = dataDir;
    fileMeta.filenameFormat     = [fn, '.tif'];
    fileMeta.nChannels          = nChannels;
    fileMeta.timePoints         = timePoints ;
    fileMeta.stackResolution    = stackResolution;
    fileMeta.swapZT             = masterSettings.swapZT;

    % first_tp is also required, which sets the tp to do individually.
    first_tp = 1 ;
    expMeta                     = struct();
    expMeta.channelsUsed        = channelsUsed ;
    expMeta.channelColor        = 1;
    expMeta.description         = 'example tube';
    expMeta.dynamicSurface      = 1;
    expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
    expMeta.fitTime             = fileMeta.timePoints(first_tp);

    %% SET DETECTION OPTIONS ==================================================
    % Load/define the surface detection parameters
    detOpts_fn = fullfile(projectDir, 'detectOpts.mat') ;
    if exist(detOpts_fn, 'file')
        disp('loading detectOptions')
        load(detOpts_fn, 'detectOptions')
    else
        outputfilename_ply='mesh_ls_' ;
        outputfilename_ls='ls_' ;
        outputfilename_smoothply = 'mesh_' ;
        init_ls_fn = 'ls_initguess' ;
        prob_searchstr = '_stab_Probabilities.h5' ;
        preilastikaxisorder = set_preilastikaxisorder; ... % axis order in input to ilastik as h5s. To keep as saved coords use xyzc
        ilastikaxisorder= 'cxyz'; ... % axis order as output by ilastik probabilities h5
        imsaneaxisorder = 'xyzc'; ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc

        % Name the output mesh directory --------------------------------------
        meshDir = [fullfile(projectDir, 'mesh_output') filesep];

        % Surface detection parameters ----------------------------------------
        detectOptions = struct( 'channel', 1, ...
            'ssfactor', 1, ...
            'niter', 35,...
            'niter0', 160, ...
            'pre_pressure', -5, ...
            'pre_tension', 0, ...
            'pressure', 0, ...
            'tension', 0.5, ...
            'post_pressure', 2, ...
            'post_tension', 3, ...
            'exit_thres', 1e-7, ...
            'foreGroundChannel', 1, ...
            'fileName', sprintf( fn, 0 ), ...
            'meshDir', meshDir, ...
            'ofn_ls', outputfilename_ls, ...
            'ofn_ply', outputfilename_ply,...
            'timepoint', 0, ...
            'zdim', 2, ...
            'ofn_smoothply', outputfilename_smoothply, ...
            'init_ls_fn', init_ls_fn, ... % set to none to load prev tp
            'run_full_dataset', projectDir,... % projectDir, ... % set to 'none' for single tp
            'radius_guess', 40, ...
            'dset_name', 'exported_data',...
            'center_guess', '200,75,75',... % xyz of the initial guess sphere ;
            'save', true, ... % whether to save images of debugging output
            'plot_mesh3d', false, ...
            'dtype', 'h5',...
            'mask', 'none',...
            'mesh_from_pointcloud', false, ...
            'prob_searchstr', prob_searchstr, ...
            'preilastikaxisorder', preilastikaxisorder, ... 
            'ilastikaxisorder', ilastikaxisorder, ... 
            'physicalaxisorder', imsaneaxisorder, ... 
            'include_boundary_faces', true, ...
            'smooth_with_matlab', 0.01) ;

        % save options
        if exist(detOpts_fn, 'file')
            disp('Overwriting detectOptions --> renaming existing as backup')
            backupfn1 = [detOpts_fn '_backup1'] ;
            if exist(backupfn1, 'file')
                backupfn2 = [detOpts_fn '_backup2'] ; 
                system(['mv ' backupfn1 ' ' backupfn2])
            end
            system(['mv ' detOpts_fn ' ' backupfn1])
        end
        disp('Saving detect Options to disk')
        save(detOpts_fn, 'detectOptions') ;
    end

    % Overwrite certain parameters for script structure
    meshDir = detectOptions.meshDir ;

    %% Define Experiment as struct
    xp = struct('fileMeta', fileMeta, ...
        'expMeta', expMeta, 'detectOptions', detectOptions) ;
    disp('done')
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PART 2: TubULAR -- surface parameterization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Now we have 3d data volumes and surfaces. Define a TubULAR object. 
    % To visualize data on these surfaces and compute how these surfaces deform
    % we now define TubULAR object.
    nU = masterSettings.nU ;
    nV = masterSettings.nV ;
    opts = struct() ;
    opts.meshDir = meshDir ;        % Directory where meshes reside
    opts.flipy = flipy ;            % Set to true if data volume axes are inverted in chirality wrt physical lab coordinates
    opts.timeInterval = timeInterval ; % Spacing between adjacent timepoints in units of timeUnits 
    opts.timeUnits = timeUnits ;    % units of time, so that adjacent timepoints are timeUnits * timeInterval apart
    opts.spaceUnits = spaceUnits ;  % Units of space in LaTeX, for ex '$mu$m' for micron
    opts.nU = nU ;                  % How many points along the longitudinal axis to sample surface
    opts.nV = nV ;                  % How many points along the circumferential axis to sample surface
    opts.normalShift = 0 ;         % Additional dilation acting on surface for texture mapping
    opts.a_fixed = 1.0 ;            % Fixed aspect ratio of pullback images. Setting to 1.0 is most conformal mapping option.
    opts.adjustlow = 0.00 ;         % floor for intensity adjustment
    opts.adjusthigh = 100. ;        % ceil for intensity adjustment (clip)
    opts.phiMethod = 'curves3d' ;   % Method for following surface in surface-Lagrangian mapping [(s,phi) coordinates]
    opts.lambda_mesh = 0.00 ;       % Smoothing applied to the mesh before DEC measurements
    opts.lambda = 0.0 ;             % Smoothing applied to computed values on the surface
    opts.lambda_err = 0.0 ;         % Additional smoothing parameter, optional
    opts.zwidth = 1 ;
    opts.nmodes = 7 ;
    % opts.t0 = xp.fileMeta.timePoints(1) ;   % reference timepoint used to define surface-Lagrangian and Lagrangian measurements
    % opts.t0 = 123 ;
    % opts.t0 = 37 ;
    % opts.t0 = 1 ;

    disp('saving xp struct and opts to disk')
    save(fullfile(dataDir, 'xp.mat'), 'xp', 'opts')
else
    disp('loading xp struct from disk')
    load(fullfile(dataDir, 'xp.mat'), 'xp', 'opts')
end

%% TubULAR class instance
disp('defining TubULAR class instance (tubi= tubular instance)')
tubi = TubULAR(xp, opts) ;
disp('done defining TubULAR instance')


%% Define global orientation frame (for viewing in canonical frame)
% Compute APDV coordinate system
alignAPDVOpts = struct() ;
alignAPDVOpts.overwrite = false ;
tubi.computeAPDVCoords(alignAPDVOpts) ;

%% Select the endcaps for the centerline computation (A and P) and a point
% along which we will form a branch cut for mapping to the plane (D).
apdvOpts = struct() ;
apdvOpts.overwrite = false ;
apdvOpts.autoAP = true ;
apdvOpts.smwindow = 0 ;
[apts_sm, ppts_sm] = tubi.computeAPDpoints(apdvOpts) ;

% Align the meshes in the APDV global frame & plot them
alignAPDVOpts.overwrite = true ;
alignAPDVOpts.forceEndpointsInside = true ;
alignAPDVOpts.normal_step = 2 ;
tubi.alignMeshesAPDV(alignAPDVOpts) ;

disp('done')

% EXTRACT CENTERLINES
% Note: these just need to be 'reasonable' centerlines for topological
% checks on the orbifold cuts. Therefore, use as large a resolution ('res')
% as possible that still forms a centerline passing through the mesh
% surface, since the centerline computed here is just for constraining the 
% mapping to the plane.
cntrlineOpts.overwrite = false ;         % overwrite previous results
cntrlineOpts.overwrite_ims = false ;     % overwrite previous results
cntrlineOpts.weight = 0.1;               % for speedup of centerline extraction. Larger is less precise
cntrlineOpts.exponent = 1.0 ;            % how heavily to scale distance transform for speed through voxel
cntrlineOpts.res = 1 ;                 % resolution of distance tranform grid in which to compute centerlines
cntrlineOpts.preview = false ;           % preview intermediate results
cntrlineOpts.reorient_faces = false ;    % not needed for our well-constructed meshes
cntrlineOpts.dilation = 0 ;              % how many voxels to dilate the segmentation inside/outside before path computation

% Note: this can take about 400s per timepoint for res=2.0, so use as big a 
%   res value as possible.
%
tubi.generateFastMarchingCenterlines(cntrlineOpts)
disp('done with centerlines')

%% Identify anomalies in centerline data
idOptions.ssr_thres = 155 ;  % distance of sum squared residuals in um as threshold for removing spurious centerlines
idOptions.overwrite = true ;
tubi.cleanFastMarchingCenterlines(idOptions) ;
disp('done with cleaning up centerlines')

% Cylinder cut mesh --> transforms a topological sphere into a topological cylinder
% Look for options on disk. If not saved, define options.
if ~exist(tubi.fileName.endcapOptions, 'file') || overwrite
    endcapOpts = struct( 'adist_thres', 30, ...  % 20, distance threshold for cutting off anterior in pix
                'pdist_thres',30, ...  % 15-20, distance threshold for cutting off posterior in pix
                'tref', tubi.t0) ;  % reference timepoint at which time dorsal-most endcap vertices are defined
    
    endcapOpts.aCapMethod = 'ball' ;
    endcapOpts.pCapMethod = 'ball' ;
    % endcapOpts.aCapMethod = 'zvalue' ;
    % endcapOpts.pCapMethod = 'zvalue' ;
    tubi.setEndcapOptions(endcapOpts) ;
    % Save the options to disk
    disp('saving endcap options to disk')
    tubi.saveEndcapOptions() ;
else
    % load endcapOpts
    tubi.loadEndcapOptions() ;
    endcapOpts = tubi.endcapOptions ;
end

methodOpts = struct() ;
methodOpts.overwrite = false ;
methodOpts.save_figs = true ;   % save images of cutMeshes along the way
methodOpts.preview = false  ;     % display intermediate results
methodOpts.quickScan = false ;
tubi.sliceMeshEndcaps(endcapOpts, methodOpts) ;

% Clean Cylinder Meshes
% This removes "ears" from the endcaps of the tubular meshes (cylindrical
% meshes)
cleanCylOptions = struct() ;
cleanCylOptions.overwrite = true ;
tubi.cleanCylMeshes(cleanCylOptions)
disp('done cleaning cylinder meshes')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORBIFOLD -> begin populating tubi.dir.mesh/gridCoords_nUXXXX_nVXXXX/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overwrite = false ;
% Iterate Through Time Points to Create Pullbacks ========================
t0_for_phi0= 7 ;  % Do this "reference" timepoint first.
tubi.t0set(t0_for_phi0)
tps = tubi.xp.fileMeta.timePoints ;
tp2do = [t0_for_phi0, fliplr(tps(tps < t0_for_phi0)), tps(tps > t0_for_phi0)] ;
for tt = tp2do
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    
    % Load the data for the current time point ------------------------
    tubi.setTime(tt) ;
    
    %----------------------------------------------------------------------
    % Create the Cut Mesh
    %----------------------------------------------------------------------
    cutMeshfn = sprintf(tubi.fullFileBase.cutMesh, tt) ;
    cutPathfn = sprintf(tubi.fullFileBase.cutPath, tt) ;
    if ~exist(cutMeshfn, 'file') || ~exist(cutPathfn, 'file') || overwrite
        if exist(cutMeshfn, 'file')
            disp('Overwriting cutMesh...') ;
        else
            disp('cutMesh not found on disk. Generating cutMesh... ');
        end
        options = struct() ;
        options.maxJitter = 15 ;
        options.maxTwistChange = 1.5 ;
        options.useMaxDeviationPtsForCorrection = false ;
        tubi.generateCurrentCutMesh(options)
        disp('Saving cutP image')
        % Plot the cutPath (cutP) in 3D
        tubi.plotCutPath(tubi.currentMesh.cutMesh, tubi.currentMesh.cutPath)
        compute_pullback = true ;
    else
        compute_pullback = false ;
    end
    
    if overwrite
        options = struct() ;
        options.overwrite = true ;
        tubi.generateCurrentUVCutMesh(options) ;
    else
        uvcutMesh = tubi.getCurrentUVCutMesh() ;
    end
    
    % Compute the pullback if the cutMesh is ok
    pbOptions = struct() ;
    pbOptions.overwrite = false ; 
    pbOptions.generate_uv = true;
    pbOptions.generate_sphi = false;
    tubi.generateCurrentPullbacks([], [], [], pbOptions) ;
    
    options = struct() ;
    options.coordSys = 'uv' ;
    tubi.coordinateSystemDemo(options)
end

%% Now transform from uv to sphi coordinates
overwrite = true ;
for tt = tp2do(4:end)
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)])
    
    % Load the data for the current time point ------------------------
    tubi.setTime(tt) ;
        
    spcutMeshOptions = struct() ;
    spcutMeshOptions.t0_for_phi0 = t0_for_phi0 ;  % which timepoint do we define corners of pullback map
    spcutMeshOptions.save_phi0patch = false ;
    spcutMeshOptions.iterative_phi0 = false ;
    spcutMeshOptions.smoothingMethod = 'none' ;
    spcutMeshOptions.overwrite = overwrite ;
    spcutMeshOptions.smoothingMethod = 'savgol' ;
    spcutMeshOptions.smoothingWidth = 35 ;
    spcutMeshOptions.smoothingOrder = 1 ;
    tubi.plotting.preview = false ;
    tubi.generateCurrentSPCutMesh([], spcutMeshOptions) ;
    
    options = struct() ;
    options.coordSys = 'sp' ;
    tubi.coordinateSystemDemo(options)
    
    % Compute the pullback if the cutMesh is ok
    pbOptions = struct() ;
    pbOptions.overwrite = false ; 
    pbOptions.generate_uv = false;
    tubi.generateCurrentPullbacks([], [], [], pbOptions) ;
        
end
disp('Done with generating spcutMeshes and cutMeshes')


%% CREATE SYNTHETIC DATA
% Evenly pattern dots at timepoint 6
overwrite = true ;
tubi.setTime(9)
mesh = tubi.getCurrentSPCutMesh() ;
zphi = mesh.sphi ;
zphi(:, 1) = zphi(:, 1) / max(zphi(:, 1)) ;
ptfn = sprintf('farthest_points_tp%06d.mat', tubi.currentTime) ;
if ~exist(ptfn, 'file') || overwrite
    nNuc = 120 ;
    % [P3d, PI] = farthest_points(mesh.v, nNuc, 'F', mesh.f, ...
    %    'Distance', 'geodesic') ;
    [p3d, pInd] = farthest_points(mesh.v, nNuc) ;
    p2d = zphi(pInd, :) ;
    close all
    figure()
    plot(p2d(:, 1), p2d(:, 2), '.')
    
    % Tube version
    % tube_zphi = zphi ;
    % tube_zphi(:, 2) = cos(zphi(:, 2)*2*pi) ;
    % tube_zphi(:, 3) = sin(zphi(:, 2)*2*pi) ;
    % [p2d, pInd] = farthest_points(tube_zphi, nNuc) ;
    % p2d(:,2) = atan2(2*pi*p2d(:, 3), 2*pi*p2d(:, 2))/2*pi + 0.5 ;
    % plot(p2d(:, 1), p2d(:, 2), '.')
    
    % Torus version
    % tube_zphi = zphi ;
    % tube_zphi(:, 1) = cos(zphi(:, 1)*2*pi) ;
    % tube_zphi(:, 2) = sin(zphi(:, 1)*2*pi) ;
    % tube_zphi(:, 3) = cos(zphi(:, 2)*2*pi) ;
    % tube_zphi(:, 4) = sin(zphi(:, 2)*2*pi) ;
    % [pts2d, pInd] = farthest_points(tube_zphi, nNuc) ;
    % p2dU = atan2(2*pi*pts2d(:, 2), 2*pi*pts2d(:, 1))/(2*pi) + 0.5 ;
    % p2dV = atan2(2*pi*pts2d(:, 4), 2*pi*pts2d(:, 3))/(2*pi) + 0.5 ;
    % p2d = [p2dU, p2dV] ;
    % clf
    % plot(p2d(:, 1), p2d(:, 2), '.')
    % % select three points in middle top and remove
    % ptid = find(p2d(:, 1) > 0.38 & p2d(:, 1) < 0.55 & p2d(:, 2) > 0.5) ;
    % if length(ptid) > 5
    %     idx = randperm(length(ptid));
    %     idx = idx(1:5);
    %     ptid = ptid(idx) ;
    % end
    % pInd(ptid) = [] ;
    % p2d(ptid, :) = [] ;
    % hold on; plot(p2d(:, 1), p2d(:, 2), 'o')
    % p3d = mesh.v(pInd, :) ;
    
    
    % Do 2D version of farthest point search
    % Nsearch = 400 ;
    % not_enough_points= true ;
    % while not_enough_points
    %     pts2d = rand(10000, 2) * 2 - [0.5, 0.5] ;
    %     pts2d = farthest_points( ...
    %         [pts2d(:, 1), sin(pts2d(:, 2)*2*pi), cos(pts2d(:, 2)*2*pi)], ...
    %         Nsearch) ;
    %     pts2d = pts2d(pts2d(:, 1) < 1 & pts2d(:, 1) > 0 & pts2d(:, 2) < 1 & pts2d(:, 2) > 0, :) ;
    %     plot(pts2d(:, 1), atan2(2*pi*pts2d(:, 2), 2*pi*pts2d(:, 3)), '.')
    %     pause(0.001)
    %     not_enough_points = length(pts2d) < nNuc ;
    %     Nsearch = Nsearch * 1.30 ;
    % end
    % pInd = pointMatch(pts2d, zphi) ;
    % p2d = zphi(pInd, :) ;
    % p3d = mesh.v(pInd, :) ;

    % Save results
    save(ptfn, 'p3d', 'p2d', 'pInd') ;
else
    load(ptfn, 'p3d', 'p2d', 'pInd') ;
end
sizeVariation = 0.5 ;
sizes = (p3d(:, 3)-mean(p3d(:, 3))).^2 ;
sizes = sizes ./ max(sizes) ;
sizes = (1 + sizes*sizeVariation) * 600 / (max(p3d(:,3)) - min(p3d(:,3))) ;

% plot it
[x,y,z] = sphere;
figure
trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none')
hold on;
for ii = 1:length(sizes)
    % scatter3(p3d(:, 1), p3d(:, 2), p3d(:, 3), sizes, 'filled')
    surf(sizes(ii)*x+p3d(ii, 1), ...
        p3d(ii, 2) + sizes(ii)*y, ...
        p3d(ii, 3) + sizes(ii)*z, 'facecolor', 'k')
end
axis equal ;


%% Intensity as a function of distance from each point 
clf
hold on;
res = 2 ;
overwrite= true ;
write_highres = false ;
for tidx = tp2do
    tp = tubi.xp.fileMeta.timePoints(tidx) ;
    disp(['Creating t = ' num2str(tp)])
    dataFn = fullfile(tubi.dir.data, sprintf('Time_%06d.tif', tp)) ;
    dataFnRes = fullfile(tubi.dir.data, sprintf('Time_highres_%06d.tif', tp)) ;
    %memDataFn = fullfile(tubi.dir.data, sprintf('Time_membrane_%06d.tif', tp)) ;
    labelFn = fullfile(tubi.dir.data, sprintf('Label_%06d.tif', tp)) ;
    if ~exist(dataFn, 'file') || ~exist(labelFn, 'file') || overwrite
        tubi.setTime(tp)
        load(fullfile(outdir, 'segs', sprintf('seg_%06d.mat', tp)), 'seg')
        
        % Make slab of segmentation
        ball = strel('sphere', 3) ;
        segSlab = imdilate(seg, ball) ;
        seg_tmp = imerode(seg, ball) ;
        segSlab = double(segSlab & ~seg_tmp) ;
        segSlabT = permute(segSlab, [2, 1, 3]) ;
        
        intens = zeros(size(seg)) ;
        nuclei3d = zeros(size(seg)) ;
        mem = zeros(size(seg)) ;
        
        % high res version
        if write_highres
            disp('preallocating high res')
            intens_res = zeros(res*size(seg)) ;
            nuclei3d_res = zeros(res*size(seg)) ;
            mem_res = zeros(res*size(seg)) ;
            % Make slab of segmentation
            ball = strel('sphere', 3*res) ;
            segSlabRes = imdilate(imresize3(double(seg), res), ball) ;
            seg_tmpRes = imerode(imresize3(double(seg), res), ball) ;
            segSlabRes = segSlabRes & ~seg_tmpRes ;
            [xxres, yyres, zzres] = meshgrid(1:res*size(seg, 1), ...
                1:res*size(seg, 2), 1:res*size(seg, 3)) ;
        end
        
        % preallocate gridded positions
        [xx, yy, zz] = meshgrid(1:size(seg, 1), 1:size(seg, 2), 1:size(seg, 3)) ;
        
        % Load the mesh for pushing forward into lab frame
        mesh = tubi.getCurrentSPCutMesh() ;
            
        p3dfn = fullfile(tubi.dir.data, sprintf('nuclei_positions_%06d.mat', tp)) ;
        if ~exist(p3dfn, 'file') || overwrite
            % Get p3d for this timepoint
            zphi = mesh.sphi ;
            zphi(:, 1) = zphi(:, 1) ./ max(zphi(:, 1)) ;
            p2d_zp = p2d+rand(size(p2d))* 0.005 ;
            p2d_zp(p2d_zp < 0) = eps ;
            p2d_zp(p2d_zp > 1) = 1-eps ;
            
            [p3d, vertexMeshFaces] = ...
                interpolate2Dpts_3Dmesh(mesh.f, zphi, mesh.v, p2d_zp) ;

            assert(~any(isnan(p3d(:))))
            
            % Make dual for membrane 
            p2d_zp1 = p2d_zp ;
            p2d_zp1(:, 2) = p2d_zp1(:, 2) + 1 ; 
            p2d_zp2 = p2d_zp ;
            p2d_zp2(:, 2) = p2d_zp2(:, 2) - 1 ;
            
            p2d_zpX1 = p2d_zp ;
            p2d_zpX1(:, 1) = p2d_zp1(:, 1) - 1 ; 
            p2d_zpX2 = p2d_zp ;
            p2d_zpX2(:, 1) = p2d_zp2(:, 1) + 1 ;
            
            p2d_zpTiled = cat(1, p2d_zp, p2d_zp1, p2d_zp2, p2d_zpX1, p2d_zpX2) ;
            DT = delaunayTriangulation(p2d_zpTiled(:, 1), p2d_zpTiled(:, 2)) ;
            triTiled = DT.ConnectivityList ;
            
            % Compare Wigner-Seitz cells:
            % clf; subplot(1, 2, 1)
            % [dualPts, regs] = voronoiDiagram(DT) ;
            % for ii = 1:length(regs)
            %     plot(dualPts(regs{ii}, 1), dualPts(regs{ii}, 2), 'k-')
            %     hold on;
            % end
            % axis equal
            % xlim([0,1]) ; ylim([0, 1])
            
            % Make a bounding box for edge cells
            Le = linspace(0, 1, 100)' ;
            boxV = [Le, 0*Le; Le, ones(size(Le)); ...
                0*Le, Le; ones(size(Le)), Le] ;
            
            % Compare to centroidal dual cells
            close all
            rescale = 3 ;
            fig = figure('Units', 'centimeters', 'visible', 'off') ;
            set(fig, 'Position', [0, 0, 60, 60])
            p2d_zpTiled(:, 2) = 1/rescale * p2d_zpTiled(:, 2) ;
            bcs = barycenter(p2d_zpTiled, triTiled) ;
            EE = edges(DT) ;
            ID = edgeAttachments(DT,EE) ;
            % get all edge attachments --> faces connected 
            hold on;
            first = true ;
            disp('sampling membrane points...')
            for eid = 1:length(ID)
                eids = ID{eid} ;
                if length(eids) > 1
                    plot([bcs(eids(1), 1), bcs(eids(2), 1)], ...
                        [bcs(eids(1), 2), bcs(eids(2), 2)], 'k-')
                    
                    % get these coordinates directly here
                    tt = linspace(0, 1, 20) ;
                    mem2d_edgeX = tt * bcs(eids(1), 1) + (1-tt)* bcs(eids(2), 1) ;
                    mem2d_edgeY = tt * bcs(eids(1), 2) + (1-tt)* bcs(eids(2), 2) ;
                    mem2d_e = [mem2d_edgeX', mem2d_edgeY'] ;
                    if first
                        memp2d = mem2d_e ;
                        first = false ;
                    else
                        memp2d = [memp2d; mem2d_e] ;
                    end
                    % else
                    % boxId = pointMatch(bcs(eids, :), boxV) ;
                    % plot([bcs(eids, 1), boxV(boxId, 1)], ...
                    %    [bcs(eids, 2), boxV(boxId, 2)], 'k-')
                end
            end
            memp2d(:, 2) = rescale * memp2d(:, 2) ;
            xlim([0,1]) ; ylim([0, 1/rescale])
            axis square
            xlim([0,1]) ; ylim([0, 1/rescale])
            axis off;
            FF = getframe() ;
            close all
            memIm = 255 - double(FF.cdata) ;
            memIm = memIm - min(memIm(:)) ;
            memIm = memIm ./ max(memIm(:)) * 255 ;
            
            memp2d = memp2d(memp2d(:, 1) > 0, :) ;
            memp2d = memp2d(memp2d(:, 2) > 0, :) ;
            memp2d = memp2d(memp2d(:, 1) < 1, :) ;
            memp2d = memp2d(memp2d(:, 2) < 1, :) ;
            
            % Push membrane points to 3d
            [memp3d, vertexMeshFaces] = ...
                interpolate2Dpts_3Dmesh(mesh.f, zphi, mesh.v, memp2d) ;
            
            % inds = pointMatch(dualPts, bcs) ;
            % dualBcs = bcs(inds, :) ;
            % 
            % for ii = 1:length(regs)
            %     plot(dualBcs(regs{ii}, 1), dualBcs(regs{ii}, 2), 'b.-')
            %     hold on;
            % end
            
            save(p3dfn, 'p3d', 'sizes', 'p2d_zp', 'memp2d', 'memp3d', 'p2d_zpTiled', ...
                'bcs', 'memIm', 'triTiled', 'EE', 'ID')
        else
            load(p3dfn, 'p3d', 'sizes', 'p2d_zp',  'memp2d', 'memp3d', 'p2d_zpTiled', ...
                'bcs', 'memIm', 'triTiled', 'EE', 'ID')
        end
        close all
        
        % Compute nuclear channel
        br0 = 0*xx ;
        if write_highres
            br0res = 0*xxres ;
        end
        for ii = 1:numel(sizes)
            if mod(ii, 10) == 1
                disp(['ii = ' num2str(ii)])
            end
            dist_to_mem = min(vecnorm(memp3d - p3d(ii, :), 2, 2)) ;
            sz = max(1, min(sizes(ii), dist_to_mem)) ;
            d2 = (xx- p3d(ii, 2)).^2 + (yy-p3d(ii, 1)).^2 + (zz-p3d(ii, 3)).^2 ;
            % inds = find(d2 < sz.^2) ;
            % d2 = d2(inds) ;
            if any(isnan(d2(:)))
                error('d2 has nan')
            end
            
            % Determine brightness as a function of distance from nucleus
            % central point:
            % br = ((sizes(ii)-min(sizes(:))) * (d2/sizes(ii).^2) + min(sizes(:))/sizes(ii).^2)...
            %     .* max(exp(-d2/(0.7*sizes(ii).^2))-0.3, 0)  .* max(0, min(0.3, (-sqrt(d2) + sizes(ii))/sizes(ii))) ;
            small_num = min(min(sizes(:)), sz) ;
            br = ((sz-small_num) * (d2/sz.^2) + small_num/sz.^2)...
                .* max(exp(-d2/(0.7*sz.^2))-0.3, 0)  .* max(0, min(0.3, (-sqrt(d2) + sz)/sz)) ;
            
            if isnan(max(br(:)))
                error('bad brightness')
            end
            assert(max(br(:)) > 0) 
            
            br = segSlabT .* br ./ max(br(:)) ;
            assert(max(br(:)) > 0) 
            
            if any(isnan(br(:)))
                error('bad division')
            end
            % plot(dd, br)
            intens = intens + br ;
            nuclei3d(br > 0) = ii ; 
            
            if write_highres
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % High res version                              %%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                dist_to_mem = min(vecnorm(memp3d * res - res* p3d(ii, :), 2, 2)) ;
                sz = max(1, min(res * sizes(ii), dist_to_mem)) ;
                d2 = (xxres- res* p3d(ii, 2)).^2 + (yyres-res* p3d(ii, 1)).^2 + (zzres-res* p3d(ii, 3)).^2 ;
                inds = find(d2 < sz.^2) ;
                d2 = d2(inds) ;
                % Determine brightness as a function of distance from nucleus
                % central point:
                % br = ((sizes(ii)-min(sizes(:))) * (d2/sizes(ii).^2) + min(sizes(:))/sizes(ii).^2)...
                %     .* max(exp(-d2/(0.7*sizes(ii).^2))-0.3, 0)  .* max(0, min(0.3, (-sqrt(d2) + sizes(ii))/sizes(ii))) ;
                small_num = min(min(sizes(:)*res), sz) ;
                br = br0res ;
                br(inds) = ((sz-small_num) * (d2/sz.^2) + small_num/sz.^2)...
                    .* max(exp(-d2/(0.7*sz.^2))-0.3, 0)  .* max(0, min(0.3, (-sqrt(d2) + sz)/sz)) ;
                if isnan(max(br(:)))
                    error('bad brightness')
                end
                assert(max(br(:)) > 0) 

                br = segSlabRes .* br ./ max(br(:)) ;

                if any(isnan(br(:)))
                    error('bad division')
                end
                % plot(dd, br)
                intens_res = intens_res + br ;
                nuclei3d_res(br > 0) = ii ; 
            end
        end
        
        % Compute membrane channel
        % [rr, cc] = find(squeeze(memIm(:,:,1))) ;
        % [Lx, Ly, ~] = size(memIm) ;
        % pts2push = [(Ly - cc + 1) ./ Ly, rr ./ Lx] ;
        % [mem3d, memVertexMeshFaces] = ...
        %     interpolate2Dpts_3Dmesh(mesh.f, zphi, mesh.v, pts2push) ;
        % mempix = round(mem3d) ;
        % disp('stuffing membrane channel...')
        % for mid = 1:length(mempix)
        %    mem(mempix(mid, 1), mempix(mid, 2), mempix(mid, 3)) = rand(1) ; 
        % end
        disp('stuffing membrane channel...')
        for qq = 1:length(memp3d)
            mem(round(memp3d(qq, 1)), round(memp3d(qq, 2)), round(memp3d(qq, 3))) = 1 - 0.3*rand(1) ;
            mem_res(round(res*memp3d(qq, 1)), ...
                round(res*memp3d(qq, 2)), ...
                round(res*memp3d(qq, 3))) = 1 - 0.3*rand(1) ;
        end
        
        clearvars xx_res yy_res zz_res mempix
        
        % Save image and volume
        if write_highres
            tmp = squeeze(mean(intens_res, 2)) ;
            tmp2 = squeeze(mean(mem_res, 2)) ;
        else            
            tmp = squeeze(mean(intens, 2)) ;
            tmp2 = squeeze(mean(mem, 2)) ;
        end
        close all
        figure('units', 'centimeters', 'position', [0,0,30,50])
        subplot(3, 1, 1)
        imagesc(tmp) ;
        axis equal ;
        subplot(3, 1, 2)
        imagesc(tmp2) ;
        axis equal
        subplot(3, 1, 3)
        imshow(cat(3, 255 * tmp, 255 * tmp2, 255*tmp))
        pause(0.01)
        
        % Check in 3d
        figure()
        plot3(memp3d(:, 1), memp3d(:, 2), memp3d(:, 3), '.')
        axis equal
        hold on;
        scatter3(p3d(:, 1), p3d(:, 2), p3d(:, 3), 50, 'filled')
        
        if write_highres
            
            imwrite(squeeze(uint16(65535*max(intens_res, [], 2))), ...
                sprintf('snap_highres_max_Time_%06d.png', tp))
            imwrite(squeeze(uint16(65535*mean(intens_res, 2))), ...
                sprintf('snap_highres_mean_Time_%06d.png', tp))
            imwrite(squeeze(uint16(65535*mean(mem_res, 2))), ...
                sprintf('snap_highres_mem_Time_%06d.png', tp))

        end
        imwrite(squeeze(uint16(65535*max(intens, [], 2))), ...
            sprintf('snap_max_Time_%06d.png', tp))
        imwrite(squeeze(uint16(65535*mean(intens, 2))), ...
            sprintf('snap_mean_Time_%06d.png', tp))
        imwrite(squeeze(uint16(65535*mean(mem, 2))), ...
            sprintf('snap_mem_Time_%06d.png', tp))
        
        
        % RGB image from label matrix 3d
        cols = distinguishable_colors(numel(sizes), [0,0,0]) ;
        RR = zeros(size(nuclei3d)) ;
        GG = zeros(size(nuclei3d)) ;
        BB = zeros(size(nuclei3d)) ;
        for colID = 1:numel(sizes)
            assert(~isempty(find(nuclei3d==colID, 1)))
            RR(nuclei3d==colID) = cols(colID, 1) ;
            GG(nuclei3d==colID) = cols(colID, 2) ;
            BB(nuclei3d==colID) = cols(colID, 3) ;
        end
        
        maxR = squeeze(uint16(65535*max(RR, [], 2))) ;
        maxG = squeeze(uint16(65535*max(GG, [], 2))) ;
        maxB = squeeze(uint16(65535*max(BB, [], 2))) ;
        imwrite(cat(3, maxR, maxG, maxB), ...
            sprintf('snap_label_Time_%06d.png', tp))
        
        % Write data volume to disk
        fulldat = reshape(uint16(intens*65535), ...
             [size(seg, 1), size(seg, 2), 1, size(seg, 3)]) ;
        fulldat(:, :, 2, :) = reshape(uint16(mem*65535), ...
             [size(seg, 1), size(seg, 2), 1, size(seg, 3)]) ;
        writeTiff5D(fulldat, dataFn, 16)
        
        
        % Write full volume for nuclei only
        % writeTiff5D(reshape(uint16(intens*65535), ...
        %      [size(seg, 1), size(seg, 2), 1, size(seg, 3)]), ...
        %      dataFn, 16)
         
        % Write full volume data for membrane only
        % writeTiff5D(reshape(uint16(mem*65535), ...
        %      [size(seg, 1), size(seg, 2), 1, size(seg, 3)]), ...
        %      memDataFn, 16)
         
        % Label matrix in 3d (nuclei segmentation)
        writeTiff5D(reshape(uint16(nuclei3d), ...
             [size(seg, 1), size(seg, 2), 1, size(seg, 3)]), ...
             labelFn, 16)
         
         if write_highres

            % Write data volume to disk -- high res
            fulldat_res = reshape(uint16(intens_res*65535), ...
                 [size(intens_res, 1), size(intens_res, 2), 1, size(intens_res, 3)]) ;
            fulldat_res(:, :, 2, :) = reshape(uint16(mem_res*65535), ...
                 [size(intens_res, 1), size(intens_res, 2), 1, size(intens_res, 3)]) ;
            writeTiff5D(fulldat_res, dataFnRes, 16)

            % Label matrix in 3d (nuclei segmentation)
            writeTiff5D(reshape(uint16(nuclei3d_res), ...
                 [size(nuclei3d_res, 1), size(nuclei3d_res, 2), 1, size(nuclei3d_res, 3)]), ...
                 labelFn, 16)
         end
    end
end

disp('done!')

%% Generate synthetic dataset with segmentation as a third channel
for tp = tubi.xp.fileMeta.timePoints
    dataFn = sprintf('Time_%06d.tif', tp) ;
    
    % read in the 2-channel data (membrane + nuclei)
    nColors = 2; 
    fn2c = fullfile('separated_channels', dataFn) ;
    dat2c = readTiff4D(fn2c, nColors) ;
    
    % read in the segmentation
    fnseg = fullfile('separated_channels', ...
        sprintf('Time_%06d_segmentation.tif', tp)) ;
    seg = readTiff4D(fnseg, 1) ;
    seg = double(permute(seg, [2, 1, 3]))/double(max(seg(:))) ;
    % smooth the segmentation with gaussian blur
    seg = imgaussfilt3(seg, 2) ;
    
    % Reshape the channels into a 4D tiff
    dat3c = zeros(size(seg, 1), size(seg, 2), 3, ...
        size(seg, 3), 'uint16') ;
    dat3c(:, :, 1, :) = reshape(dat2c{1}, ...
         [size(seg, 1), size(seg, 2), 1, size(seg, 3)]) ;
    dat3c(:, :, 2, :) = reshape(dat2c{2}, ...
         [size(seg, 1), size(seg, 2), 1, size(seg, 3)]) ;
    dat3c(:, :, 3, :) = reshape(uint16(seg*65535), ...
         [size(seg, 1), size(seg, 2), 1, size(seg, 3)]) ;
    
    writeTiff5D(dat3c, dataFn, 16)
    tmp = readTiff4D(dataFn) ;
    
    % Take a max intensity projection
    mip = reshape(max(dat3c, [], 2), [size(dat3c, 1), 3, size(dat3c, 4)]) ;
    mip = permute(mip, [1,3,2]) ;
    mipFn = sprintf('mip_Time_%06d.png', tp) ;
    imwrite(mip, mipFn) 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 3: Further refinement of dynamic meshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Smooth the sphi grid meshes in time ====================================
options = struct() ;
options.overwrite = true ;
options.width = 1 ;  % width of kernel, in #timepoints, to use in smoothing meshes
tubi.smoothDynamicSPhiMeshes(options) ;

%% Plot the time-smoothed meshes
tubi.plotSPCutMeshSmRS(options) ;

%% Inspect coordinate system charts using smoothed meshes
for tidx = 1:length(tubi.xp.fileMeta.timePoints)
    tp = tubi.xp.fileMeta.timePoints(tidx) ;
    options = struct() ;
    options.coordSys = 'sp' ;
    tubi.setTime(tp) ;
    tubi.coordinateSystemDemo(options)
end

%% Redo Pullbacks with time-smoothed meshes ===============================
disp('Create pullback using S,Phi coords with time-averaged Meshes')
for tt = tubi.xp.fileMeta.timePoints
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    tidx = tubi.xp.tIdx(tt);
    
    % Load the data for the current time point ------------------------
    tubi.setTime(tt) ;
    
    % Establish custom Options for MIP --> choose which pullbacks to use
    pbOptions = struct() ;
    pbOptions.numLayers = [0 0] ; % how many onion layers over which to take MIP
    pbOptions.generate_spsm = true ;
    pbOptions.generate_sp = false ;
    pbOptions.generate_uv = false ;
    pbOptions.overwrite = true ;
    pbOptions.falseColors = [1,0,1; 0,1,0] ;
    tubi.generateCurrentPullbacks([], [], [], pbOptions) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 4: Computation of tissue deformation, with in-plane and out-of-plane flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TILE/EXTEND SMOOTHED IMAGES IN Y AND RESAVE ============================
% Skip if already done
options = struct() ;
options.coordsys = 'spsm' ;
tubi.doubleCoverPullbackImages(options)
options.coordsys = 'sp' ;
tubi.doubleCoverPullbackImages(options) 
disp('done')

%% PERFORM PIV ON PULLBACK MIPS ===========================================
% % Compute PIV either with built-in phase correlation or in PIVLab
options = struct() ;
options.overwrite = true ;
tubi.measurePIV2d(options) ;

%% Measure velocities =====================================================
disp('Making map from pixel to xyz to compute velocities in 3d for smoothed meshes...')
options = struct() ;
options.show_v3d_on_data = false ;
tubi.measurePIV3d(options) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lagrangian dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pullback pathline time averaging of velocities
options = struct() ;
tubi.timeAverageVelocities(options)
%% Velocity plots for pathline time averaging 
options.plot_vxyz = false ;
options.invertImage = true ;
options.averagingStyle = 'Lagrangian'; 
tubi.plotTimeAvgVelocities(options)
% Divergence and Curl (Helmholtz-Hodge) for Lagrangian
options = struct() ;
options.averagingStyle = 'Lagrangian' ;
options.lambda = 0 ;
options.lambda_mesh = 0 ; 
tubi.helmholtzHodge(options) ;

% Compressibility & kinematics for Lagrangian
options = struct() ;
tubi.measureMetricKinematics(options)

%% Metric Kinematics Kymographs & Correlations -- Bandwidth Filtered
options = struct() ;
tubi.plotMetricKinematics(options)

%% Pullback pathlines connecting Lagrangian grids
options = struct() ;
tubi.measurePullbackPathlines(options)

%% Query velocities along pathlines
options = struct() ;
tubi.measurePathlineVelocities(options)
% plot the pathline velocities 
options = struct() ;
options.gridTopology = 'triangulated' ;
tubi.plotPathlineVelocities(options)

% Measure Pathline Kinematics
options = struct() ;
tubi.measurePathlineMetricKinematics(options)

% Plot Pathline Kinematics
options = struct() ;
tubi.plotPathlineMetricKinematics(options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create ricci mesh at t0 to measure Beltrami coefficient in pathlines
options = struct() ;
options.climit = 1 ;
options.coordSys = 'ricci' ;
tubi.measureBeltramiCoefficient(options) ;

%% Strain rate (epsilon = 1/2 (djvi+divj) -vn bij)
options = struct() ;
tubi.measureStrainRate(options) 

%% Plot time-averaged strain rates in 3d on mesh
options = struct() ;
tubi.plotStrainRate3DFiltered(options) 

%% Kymograph strain rates
options = struct() ;
options.clim_trace = 0.05 ;
options.clim_deviatoric = 0.05 ;
tubi.plotStrainRate(options)

% Measure strain rate along pathlines
options = struct() ;
options.overwriteImages = false ;
options.plot_dzdp = false ;
tubi.measurePathlineStrainRate(options)

%% Measure divergence and out-of-plane deformation along pathlines
tubi.measurePathlineMetricKinematics()

% Pathline strain rate plots
options = struct() ;
options.climit = 0.05 ;
options.climitWide = 1.0 ;
tubi.plotPathlineStrainRate(options)

%% Measure strain along pathlines -- note this is from pathlines, not integrating rates
options = struct() ;
options.plot_dzdp = false ;
options.climitInitial = 0.05 ;
options.climitRamp = 0.01 ;
options.climitRatio = 1 ;
tubi.measurePathlineStrain(options)
tubi.plotPathlineStrain(options)



%% PCA decomposition
pcaTypes = {'vnVector', 'v3d', 'vt', 'H2vn', 'vnScalar', 'divv', 'gdot'} ;
% pcaTypes = {'H2vn', 'vnScalar', 'divv', 'gdot'} ;
options = struct('overwrite', true, ...
    'overwriteImages', true) ;
options.pcaTypes = pcaTypes ;
% options.meshStyles = 'sphi' ;
tubi.spaceUnits = [char(181) 'm'] ;
tubi.getPCAoverTime(options)

%% Laplace-Beltrami Spectral (LBS) decomposition
close all; clc;

% lbsTypes = {'vnVector', 'v3d', 'vt', 'H2vn', 'vnScalar', 'divv', 'gdot'} ;
lbsTypes = {'H2vn', 'vnScalar', 'divv', 'gdot'} ;
options = struct('overwrite', true, ...
    'overwriteImages', true) ;
options.lbsTypes = lbsTypes ;
% options.meshStyles = 'sphi' ;
tubi.spaceUnits = [char(181) 'm'] ;
tubi.getLBSoverTime(options)
