function [acom_sm, pcom_sm, dcom] = computeAPDCOMs(QS, opts)
%[acom_sm, pcom_sm, dcom] = COMPUTEAPDCOMS(opts)
% Compute the anterior, posterior, and dorsal centers of mass from training
%
% Parameters
% ----------
% opts : struct with fields
%   - timePoints
%   - dorsal_thres : float between 0 and 1
%   - anteriorChannel : int
%   - posteriorChannel
%   - dorsalChannel
%   - overwrite : bool
%   - apdvoutdir : string
%   - meshDir : string
%   - apdProbFileName : string
%   - preview_com : bool
%   - check_slices : bool
%   - axorder : length 3 int array
%   - smwindow : float or int (optional, default=30)
%       number of timepoints over which we smooth
%   - preview : bool (optional, default=false)
%       
%
% OUTPUTS
% -------
% apdv_coms_from_training
% 
% NPMitchell 2020

timePoints = QS.xp.fileMeta.timePoints ;
apdvoutdir = QS.dir.cntrline ;
meshDir = QS.dir.mesh ;
axorder = QS.data.axisOrder ;
apdProbFileName = QS.fullFileBase.apdProb ;

% Default options
overwrite = false ; 
preview_com = false ;
check_slices = false ;

% Unpack opts
dorsal_thres = opts.dorsal_thres ;
anteriorChannel = opts.anteriorChannel ;
posteriorChannel = opts.posteriorChannel ;
dorsalChannel = opts.dorsalChannel ;
if isfield(opts, 'overwrite')
    overwrite = opts.overwrite ;
end
if isfield(opts, 'preview_com')
    preview_com = opts.preview_com ;
end
if isfield(opts, 'check_slices')
    check_slices = opts.check_slices ;
end

% Default valued options
smwindow = 30 ;
if isfield(opts, 'smwindow')
    smwindow = opts.smwindow ;
end

dcomname = fullfile(meshDir, 'dcom_for_rot.txt') ;
rawapdvname = QS.fileName.apdv ;
rawapdvmatname = fullfile(apdvoutdir, 'apdv_coms_from_training.mat') ;
preview = false ;
if isfield(opts, 'preview')
    preview = opts.preview ;
end

%% Iterate through each mesh to compute acom(t) and pcom(t). Prepare file.
acoms = zeros(length(timePoints), 3) ;
pcoms = zeros(length(timePoints), 3) ;
load_from_disk = false ;
if exist(rawapdvname, 'file') && ~overwrite
    load_from_disk = true ;
    try
        h5create(rawapdvname, '/acom_sm', size(acoms)) ;
        load_from_disk = false ;
    catch
        try
            acom_sm = h5read(rawapdvname, '/acom_sm') ;
            disp('acom_sm already exists')
        catch
            load_from_disk = false;
        end
    end
    try
        h5create(rawapdvname, '/pcom_sm', size(pcoms)) ;
        load_from_disk = false ;
    catch
        try
            pcom_sm = h5read(rawapdvname, '/pcom_sm') ;
            disp('pcom_sm already exists')
        catch
            load_from_disk = false;
        end
    end
end
if ~load_from_disk
    disp('acom and/or pcom not already saved on disk. Compute them')
end

disp(['Load from disk? =>', num2str(load_from_disk)])

%% Compute acom and pcom if not loaded from disk -- RAW XYZ coords
if ~load_from_disk || overwrite
    if ~exist(rawapdvmatname, 'file') || overwrite
        for tidx = 1:length(timePoints)
            tt = timePoints(tidx) ;
            %% Load the AP axis determination
            msg = ['Computing acom, pcom for ' num2str(tt) ] ;
            disp(msg)
            thres = 0.5 ;
            apfn = sprintf(apdProbFileName, tt);
            apdat = h5read(apfn, '/exported_data');

            % rawfn = fullfile(rootdir, ['Time_' timestr '_c1_stab.h5' ]);
            % rawdat = h5read(rawfn, '/inputData');
            adat = squeeze(apdat(anteriorChannel,:,:,:)) ;
            pdat = squeeze(apdat(posteriorChannel,:,:,:)) ;

            % define axis order: 
            % if 1, 2, 3: axes will be yxz
            % if 1, 3, 2: axes will be yzx
            % if 2, 1, 3: axes will be xyz (ie first second third axes, ie --> 
            % so that bright spot at im(1,2,3) gives com=[1,2,3]
            adat = permute(adat, axorder) ;
            pdat = permute(pdat, axorder) ;
            options.check = preview_com ;
            disp('Extracting acom')
            options.color = 'red' ;
            acom = com_region(adat, thres, options) ;
            disp('Extracting pcom')
            options.color = 'blue' ;
            pcom = com_region(pdat, thres, options) ;
            clearvars options
            % [~, acom] = match_training_to_vertex(adat, thres, vertices, options) ;
            % [~, pcom] = match_training_to_vertex(pdat, thres, vertices, options) ;
            acoms(tidx, :) = acom ;
            pcoms(tidx, :) = pcom ; 
        end
        % Save raw data to .mat
        save(rawapdvmatname, 'acoms', 'pcoms')
        clearvars adat pdat
    else
        % load raw data from .mat
        load(rawapdvmatname, 'acoms', 'pcoms')
    end
    disp('done determining acoms, pcoms')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Smooth the acom and pcom data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Smoothing acom and pcom...')
    acom_sm = 0 * acoms ;
    pcom_sm = 0 * acoms ;
    % fraction of data for smoothing window
    smfrac = smwindow / double(length(timePoints)) ;  
    acom_sm(:, 1) = smooth(timePoints, acoms(:, 1), smfrac, 'rloess');
    pcom_sm(:, 1) = smooth(timePoints, pcoms(:, 1), smfrac, 'rloess');
    acom_sm(:, 2) = smooth(timePoints, acoms(:, 2), smfrac, 'rloess');
    pcom_sm(:, 2) = smooth(timePoints, pcoms(:, 2), smfrac, 'rloess');
    acom_sm(:, 3) = smooth(timePoints, acoms(:, 3), smfrac, 'rloess');
    pcom_sm(:, 3) = smooth(timePoints, pcoms(:, 3), smfrac, 'rloess');
    
    if preview
        plot(timePoints, acoms - mean(acoms,1), '.')
        hold on
        plot(timepts, acom_sm - mean(acoms, 1), '-')
        title('Smoothed COMs for AP')
    end
    clear acom pcom
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save smoothed anterior and posterior centers of mass ===============
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        h5create(rawapdvname, '/acom', size(acoms)) ;
    catch
        disp('acom already exists as h5 file. Overwriting.')
    end
    try
        h5create(rawapdvname, '/pcom', size(pcoms)) ;
    catch
        disp('pcom already exists as h5 file. Overwriting.')
    end
    try
        h5create(rawapdvname, '/acom_sm', size(acom_sm)) ;
    catch
        disp('acom_sm already exists as h5 file. Overwriting.')
    end
    try
        h5create(rawapdvname, '/pcom_sm', size(pcom_sm)) ;
    catch
        disp('pcom_sm already exists as h5 file. Overwriting.')
    end
    h5write(rawapdvname, '/acom', acoms) ;
    h5write(rawapdvname, '/pcom', pcoms) ;
    h5write(rawapdvname, '/acom_sm', acom_sm) ;
    h5write(rawapdvname, '/pcom_sm', pcom_sm) ;
    clear acoms pcoms
else
    disp('Skipping, since already loaded acom_sm and pcom_sm')
    if preview
        acom_sm = h5read(rawapdvname, '/acom_sm');
        pcom_sm = h5read(rawapdvname, '/pcom_sm');
        plot3(acom_sm(:, 1), acom_sm(:, 2), acom_sm(:, 3))
        hold on;
        plot3(pcom_sm(:, 1), pcom_sm(:, 2), pcom_sm(:, 3))
        xlabel('x [subsampled pix]')
        ylabel('y [subsampled pix]')
        zlabel('z [subsampled pix]')
        axis equal
        pause(1)
    end
end

disp('done with AP COMs')


%% Dorsal COM for first timepoint
tt = timePoints(1) ;
% load the probabilities for anterior posterior dorsal
apfn = sprintf(apdProbFileName, tt);
disp(['Reading ' apfn])
apdat = h5read(apfn, '/exported_data');
ddat = permute(squeeze(apdat(dorsalChannel, :, :, :)), axorder) ;
options.check = preview_com ; 
options.check_slices = check_slices ; 
options.color = 'green' ;

% Load dcom if already on disk
if exist(dcomname, 'file') && ~overwrite
    disp('Loading dorsal COM from disk')
    dcom = dlmread(dcomname) ;
else
    if exist(dcomname, 'file')
        disp('Overwriting existing dorsal COM on disk')
    else
        disp('Computing dorsal COM for the first time')
    end
    search4com = true ;
    % start with a threshold == dorsal_thres, iteratively lower if
    % necessary
    tmp_dorsal_thres = dorsal_thres ;
    while search4com 
        try
            msg = 'Finding com region of dorsal data: ' ;
            disp([msg 'thres=' num2str(tmp_dorsal_thres)])
            close all
            dcom = com_region(ddat, tmp_dorsal_thres, options) ;
            search4com = false ;
        catch
            disp('no region found, lowering dorsal threshold for prob cloud') ;
            tmp_dorsal_thres = 0.9 * tmp_dorsal_thres ;
        end

        if tmp_dorsal_thres < 1e-9
            msg = 'Could not find any dorsal signal within 1e-9' ;
            disp(msg)
            disp('Showing dorsal signal volume, leaf by leaf')
            clf; set(gcf, 'visible', 'on')
            for qq=1:size(ddat, 1)
                imshow(squeeze(ddat(qq, :, :)))
                title(['dorsal signal, x=' num2str(qq)])
                pause(0.01)                            
            end
            error(msg)
        end
    end

    % SAVE DCOM
    dlmwrite(dcomname, dcom) ;
end
%%%%%%%%%%%%%%%%%%%%%%
if preview
    % % disp('Showing dorsal segmentation...')
    % clf
    % for slice=1:2:size(ddat, 2)
    %     im = squeeze(ddat(:, slice, :)) ;
    %     % im(im < dorsal_thres) = 0 ;
    %     imshow(im)
    %     xlabel('x')
    %     ylabel('z')
    %     hold on
    %     plot(dcom(:, 1), dcom(:, 3), 'o')
    %     title([num2str(slice) '/' num2str(size(apdat, 3))])
    %     pause(0.001)
    % end
    %%%%%%%%%%%%%%%%%%%%%%
    fig = figure ;
    disp('Displaying mesh in figure ...')
    % iso = isosurface(rawdat, 880) ;
    % patch(iso,'facecolor',[1 0 0],'facealpha',0.1,'edgecolor','none');
    % view(3)
    % camlight
    % hold on;
    trimesh(fvsub.faces, ...
        vtx_sub(:, 1), vtx_sub(:,2), vtx_sub(:, 3), ...
        vtx_sub(:, 1), 'edgecolor', 'none', 'FaceAlpha', 0.1) ;
    hold on;
    plot3(acom(1), acom(2), acom(3), 'ro')
    plot3(pcom(1), pcom(2), pcom(3), 'bo')
    plot3(dcom(1), dcom(2), dcom(3), 'go')
    xlabel('x [subsampled pixels]')
    ylabel('y [subsampled pixels]')
    zlabel('z [subsampled pixels]')
    title('Original mesh in subsampled pixels, with APD marked')
    axis equal
    %%%%%%%%%%%%%%%%%%%%%%
    waitfor(fig)
end


disp('done')