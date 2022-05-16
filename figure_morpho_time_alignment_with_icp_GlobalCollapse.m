%% Align meshes in Space and Time to a reference mesh
% for dynamicAtlas paper
% NPM based on Isaac Breinyn 
%
% This is an edited combination of mesh_temporal_alignment.m and
% plot_aligned meshes.m
% (Master script)
%
% To run before: compute_mesh_surfacearea_volume.m, setup.m
%
% Note: All scaling and alignment is done relative to CAAX Excellent
% Note: Previously this code had run over reftime, then ctime. Now ctime
% then reftime, so minddssr and other objects have been effectively
% 'transposed'.
%
% ICPt(cc, ii) is the reference dataset timepoint match of dataset cc's
% iith timepoint.
%
%
% rlist: reference timepoints for which ssr is computed
% minddssr : (int) actual index of reference timpoint that minimizes ssr
%   "minimum index for ssr"
% numpts1 is the number of vertices in the subsampled reference surface
% numpts2 is the number of vertices in the subsampled current surface being
%   compared to the reference
% 

%% Clean MATLAB and add paths
clear
close all
clc
codeDir = '/mnt/data/code/' ;
gutDir = fullfile(codeDir, 'gut_matlab') ;
atlasDir = fullfile(codeDir, 'dynamicAtlas') ;

addpath(codeDir)
addpath(gutDir)
addpath(fullfile(gutDir, 'addpath_recurse'))
addpath(fullfile(gutDir, 'curve_functions'))
addpath_recurse(fullfile(gutDir, 'plotting'))
% addpath(fullfile(imsaneDir, 'generalfunctions'))
addpath(fullfile(gutDir, 'data_handling'))
addpath(fullfile(atlasDir, 'timeline_handling'))
addpath_recurse(fullfile(gutDir, 'toolbox_fast_marching'))
addpath(fullfile(codeDir, 'imsane_for_git', 'imsane', 'generalfunctions'))

outdir = '/mnt/data/analysis/' ;
cd(outdir)
icpDir = fullfile(outdir, 'ICP_Plots_2022/') ;
if ~exist(icpDir, 'dir')
    mkdir(icpDir) ;
end
%% Global options
sigmaTime = 5 ;         % smoothing window for errorbar plot
ssample_factor = 40 ;   % subsampling of point clouds for ICP
rsubsampling = 1 ;     % subsampling of timepoints for reference dataset
csubsampling = 1 ;     % subsampling of timepoints for matched dataset
t0ref = 113 ;           % reference time index
overwrite_figures = false ;  % overwrite current results on disk, if any
ignoreLast = 10 ;
ignoreFirst = 30 ;
ssfactorMedianThres = 8 ;
fs = 10 ;
extn = [sprintf('_ss%02d', rsubsampling) '_ssR'] ;
fitrange = max(3, round(20/rsubsampling)) ;         
overwrite = false ;      % overwrite previous results
if overwrite == 1
    disp('WARNING: OVERWRITE IS TURNED ON');
end
% Plotting
figWidth = 7 ;
figHeight = 6 ;

%% Create Cell Array with Mesh Metadata
dmap = buildLookupMap('/mnt/data/analysis/lookupMeta_only6.txt'); % build the map that contains all data and its metadata
mca = {} ; % initiate mesh cell array

%% Iterate over each marker
timestamps = [] ;
for mi = 1:length(dmap.keys)
    labels = dmap.keys ;
    label = labels{mi} ;
    % Cycle through all datasets of this marker
    for jj=1:length(dmap(label).folders)
        disp([label ' jj = ' num2str(jj)])
        % get col
        if strcmp(label, 'caax')
            col = 1 ;
        elseif strcmp(label, 'hrfp')
            col = length(dmap('caax').folders) + 1 ;
        elseif strcmp(label, 'la')
            col = length(dmap('caax').folders) + length(dmap('hrfp').folders) + 1 ;
        else
            error(['did  not recognize label: ' label])
        end
        
        meshes = dir(fullfile(dmap(label).folders{jj}, 'alignedMesh',  'mesh_apical_stab_0*_APDV_um.ply')) ;
        
        if isempty(meshes)
            dmap(label).folders{jj}
            disp(' --> was not found!')
        end
        
        % Write each mesh to the "mesh cell array" mca
        for kk = 1:length(meshes)
            mca{col+jj-1, kk} = meshes(kk);
            
            % Get timestamp of this mesh
            tmp = strsplit(mca{col+jj-1, kk}.name, '_0') ;
            tmp = strsplit(tmp{2}, '_APDV') ;
            timestamps(col+jj-1, kk) = str2double(tmp{1});
        end
    end
end
disp('done building mca & timestamps')
clearvars kk jj

% Define maximum number of timepoints
max_ntp = size(mca, 2) ;
ndatasets = size(mca, 1) ;

% Build labels cell from dmap, also define number of timepoints ntps
labels = cell(1, ndatasets) ;
areas = cell(1, ndatasets) ;
ntps = zeros(1, ndatasets) ;
dmyk = 1 ;
dkeys = dmap.keys ;
for kk = 1:length(dkeys)
    key = dkeys{kk} ;
    for jj = 1:length(dmap(key).folders)
        datestamp = strsplit(dmap(key).folders{jj}, '/') ;
        datestamp = datestamp{5} ;
        % check that this looks like a timestamp
        assert(strcmp(datestamp(1:2), '20'))
        datestamp = strsplit(datestamp, '_') ;
        datestamp = datestamp{1} ;
        datestamp = datestamp(3:end-4) ;
        if strcmp(key, 'la')
            jjkey = 'LifeAct' ;
        else
            jjkey = upper(key) ;
        end
        
        labels{dmyk} = [ jjkey ' ' datestamp] ;
        ntps(dmyk) = length(dmap(key).area{jj}) ;
        areas{dmyk} = dmap(key).area{jj} ;
        fluors{dmyk} = key ;
        datestamps{dmyk} = datestamp ; 
        
        % update the labels index
        dmyk = dmyk + 1 ;
    end
end
% cleanup from labels generation
clearvars dmyk jjkey key datestamp kk
disp('done with labels generation')

%% Declare the range of morphological time and datasets
rlist = 1:rsubsampling:max_ntp;
clist = 1:csubsampling:max_ntp;

%% Use ICP to align meshes in time
forceTrue = false ;
minddssr = zeros(ndatasets, max_ntp); % initiate matrix that will store SSR data for each dataset
ssr_minimum = 0*minddssr ;
numpts1 = zeros(ndatasets, max_ntp) ;
numpts2 = zeros(ndatasets, max_ntp) ;


for refID = 1:ndatasets
    % rename experiment ID from referenceID and list of labels
    refExptID = strrep(labels{refID}, ' ', '') ;

    mindFn = fullfile(icpDir, sprintf('minddssr_ref%s_rsub%03d.mat', refExptID,rsubsampling)) ;
    if exist(mindFn, 'file') && ~overwrite
        answer = questdlg('Min matrix already exists? Overwrite?') ;
    else
        answer = 'Yes' ;
    end
    % Note that c=1 is the reference dataset, so cycle through all others
    % starting with index=2
    corrPaths = {} ;
    corrErrors = {} ;
    ssrPaths = {} ;
    ssrDir = fullfile(icpDir, 'ssr') ;
    if strcmp(answer, 'Yes')
        % Ensure the output directory for the plots
        if ~exist(ssrDir, 'dir')
            mkdir(ssrDir)
        end
        for cc = 1:ndatasets
            cExptID = strrep(labels{cc}, ' ', '') ;
            
            % define the SSR directory for this dataset
            ssrccDir = fullfile(ssrDir, sprintf(['dataset_' fluors{cc} '_' datestamps{cc}])) ;
            if ~exist(ssrccDir, 'dir')
                mkdir(ssrccDir)
            end

            % Load the areas if we use them to do SSR normalization. Otherwise,
            % use length of point cloud.
            area = areas{cc} ;
            lastTP = ntps(cc) ;
            ntimepoints_ref = ntps(refID) ;
            assert(length(area) == lastTP)
            ssrM = [] ;

            % Consider each dataset TP (t_c) and match to reference mesh
            for ii = clist(clist < (lastTP + 1)) % index of dataset thats being parsed against reference
                
                tic 
                
                % Consider the mesh only if mca{i,c} is populated with
                % struct that has field 'name'.
                if isfield(mca{cc, ii}, 'name') == 1

                    % Check if SSR(cc, ii, rlist) has already been saved
                    ssrii_fn = fullfile(ssrccDir, ...
                        sprintf('ssr_ref%s_tp%04d_rsub%03d_ptsub%03d.mat',...
                        refExptID, ii, rsubsampling, ssample_factor)) ;

                    % Assume that we must redo calculation unless proven
                    % that otherwise if file exists and overwrite == false
                    % ---------------------------------------------------
                    redoii = true ;
                    if exist(ssrii_fn, 'file') && ~overwrite
                        disp(['Loading ssr from: ' ssrii_fn])
                        tmp = load(ssrii_fn, 'rlist', 'ssr', 'numpts1_cc', 'numpts2_cc') ;
                        % Check that the rlist is indeed the same
                        if length(tmp.rlist) == length(rlist)
                            if all(tmp.rlist == rlist)
                                % Calculation exists for the same rlist for this ii
                                redoii = false ;
                                % assign ssr for this cc, ii
                                ssr = tmp.ssr ;
                                try
                                    numpts1(cc,:) = tmp.numpts1_cc ;
                                    numpts2(cc,:) = tmp.numpts2_cc ;
                                catch
                                    disp('no numpts on disk!')
                                end
                            end
                        end
                    end

                    % (Re)Compute the SSR(cc, ii, rr) if we could not load ssr
                    if redoii
                        % prepare the dataset pointcloud
                        cCloud = pcread(fullfile(mca{cc, ii}.folder, mca{cc, ii}.name)); % the mesh in question
                        cxyz = cCloud.Location(1:ssample_factor:end,:) ;

                        % preallocate ssr
                        ssr = zeros(length(rlist), 1) ;
                        tforms = zeros(length(rlist), 4, 4) ;

                        % Consider each reference mesh to find match
                        % Use dummy index qq to allow subsampled rlist
                        for qq = 1:length(rlist)
                            rr = rlist(qq) ;
                            if rr < ntimepoints_ref + 1
                                if mod(qq, 50) == 0 || qq == 1
                                    disp(['dataset_c=',sprintf('%d',cc), ...
                                        ' timepoint_i=',sprintf('%d/%d',ii,lastTP), ...
                                        ' against refTP_r=',sprintf('%d',rr), ...
                                        ' of reference dataset=', sprintf('%s', refExptID)])
                                end

                                % Load the CAAX mesh used for comparison
                                refCloud = pcread(fullfile(mca{refID, rr}.folder, mca{refID, rr}.name));
                                refxyz = refCloud.Location;
                                refxyz = refxyz(1:ssample_factor:end,:);
                                tform = pcregistericp(cCloud,refCloud,'Extrapolate',true); % extracting the 4x4 transform from cmesh to refmesh
                                cxyzpad = horzcat(cxyz, (zeros(length(cxyz),1)+1)); % pad the locations with a ones vector
                                cxyztrans = cxyzpad*tform.T;    % apply the transform
                                cxyztrans = cxyztrans(:,1:3);   % remove padding

                                numpts1(cc,rr) = length(refxyz) ;
                                numpts2(cc,rr) = length(cxyztrans) ;

                                % find sum of squared residuals using two different ref orders
                                % and finding geometric mean to penalize
                                % dissimilarities in both directions of time
                                indx1 = pointMatch(refxyz, cxyztrans);
                                indx2 = pointMatch(cxyztrans, refxyz);

                                % Normalize by the number of points considered in each
                                % mesh separately to get true mean distance between
                                % correspondence points
                                ssr1 = sum(sum((cxyztrans(indx1,:)-refxyz).^2,2)) / length(indx1) ;
                                ssr2 = sum(sum((refxyz(indx2,:)-cxyztrans).^2,2)) / length(indx2) ;
                                % Geometric mean of the two SSRs
                                ssr(qq) = sqrt(ssr1*ssr2);
                                tforms(qq, :, :) = tform.T ;
                            end

                        end
                        % save matching curve for this cc tp == ii
                        disp(['saving tforms: ' ssrii_fn])
                        numpts1_cc = numpts1(cc, :) ;
                        numpts2_cc = numpts2(cc, :) ;
                        save(ssrii_fn, 'rlist', 'ssr', 'tforms', ...
                            'numpts1_cc', 'numpts2_cc')
                    end

                    ssrM(ii, :) = ssr ;                

                    % -------------------------------------------------
                    % 1D FITTING TO SSR CURVE AS FUNCTION OF REFTIME
                    % -------------------------------------------------
                    % Apply median filter if we have a fine enough sampling
                    if rsubsampling < ssfactorMedianThres
                        ssr = movmedian(ssr, 11);
                    end
                    [minssr_value,ix] = min(ssr);

                    %find polyfit of min to find error of match
                    if ix + fitrange > length(ssr)
                        endmax = length(ssr) ;
                    else
                        endmax = ix + fitrange ;
                    end
                    if ix > fitrange && (ix + fitrange) <= length(rlist)
                        xpoly = rlist(ix-fitrange:endmax) ;
                        ypoly = ssr(ix-fitrange:endmax) ;
                    elseif ix < fitrange
                        xpoly = rlist(1:2*fitrange);
                        ypoly = ssr(1:2*fitrange);
                    elseif ix + fitrange > length(rlist)
                        xpoly = rlist(length(rlist)-2*fitrange:end) ;
                        ypoly = ssr(length(rlist)-2*fitrange:end) ;
                    elseif ix == fitrange && (ix + fitrange) <= length(rlist)
                        xpoly = rlist(ix-fitrange+1:endmax) ;
                        ypoly = ssr(ix-fitrange+1:endmax) ;
                    else
                        error('could not assign xpoly & ypoly')
                    end

                    % Fit to parabola
                    [f, S] = polyfix(xpoly', ypoly, 2, rlist(ix), min(ssr), rlist(ix), 0) ;
                    [poly, delta] = polyval(f, xpoly, S) ;
                    Rinv = inv(S.R) ;
                    covm = (Rinv*Rinv')*S.normr^2/S.df ;
                    vara = covm(1,1) ;
                    varb = covm(2,2) ;
                    varc = covm(3,3) ;
                    dix = ( ix* 0.5 ) * sqrt(vara/f(1)^2 + varb/f(2)^2) ;
                    disp(['Uncertainty = ' num2str(dix)])
                    % store minimum index as minddssr
                    minddssr(cc, ii) = rlist(ix);
                    % Place reference timepoint that matches this (cc, ii)
                    % tuple into minname
                    minname{cc, ii} = mca{1, rlist(ix)}.name;
                    minerror(cc, ii) = dix ; % store error to use later in plotting
                    minweights(cc, ii) = 1/dix^2 ;
                    ssr_minimum(cc,ii) = min(ssr) ;
                    assert(abs(polyval(f, rlist(ix)) - min(ssr))<1e-7)

                    % make figure invisible for speedup
                    figDirOut = fullfile(ssrccDir, ['figs_' sprintf('rsub%03d',rsubsampling)]) ;
                    if ~exist(figDirOut, 'dir')
                        mkdir(figDirOut)
                    end
                    figoutfn = fullfile(figDirOut, ['ICPPlot2_Dataset_', ...
                        sprintf('%s',cExptID),'_', ...
                        sprintf('%s',refExptID),'_', ...
                        sprintf('%03d',ii),'.png']) ;    
                    if redoii || ~exist(figoutfn, 'file') || overwrite_figures
                        fig = figure('visible', 'off') ;
                        plot(rlist - t0ref, ssr, '.');
                        hold on
                        plot(xpoly - t0ref, ypoly, '-') ;
                        plot(rlist(ix) - t0ref, minssr_value, 'o')
                        plot(xpoly - t0ref,poly, 'color', 'magenta', 'linewidth', 1) ;
                        ylim([min(ssr)-10 max(ssr)+10])
                        title(['dataset ', cExptID, ' vs ', refExptID,...
                            ': $t=$', sprintf('%03d', ii), '$\rightarrow$', ...
                            sprintf('%03.1f', rlist(ix) - t0ref), '$\pm$', ...
                            sprintf('%0.1f', dix), ' min'], ...
                            'Interpreter', 'Latex')
                        xlabel('Reference Timepoint [min]', 'Interpreter', 'Latex');
                        labstr = 'Mismatch, $\langle \Sigma |\mathrm{min}(\vec{x}_c - \vec{x}_{\mathrm{ref}})|^2 \rangle$ [$\mu$m]$^2$' ;
                        ylabel(labstr, 'Interpreter', 'Latex');
                        disp(['Saving ' figoutfn]);
                        saveas(fig, figoutfn);
                        close all
                    end
                end
                toc
            end

%             % Now ssrM is fully filled in
%             corrPathFn = fullfile(ssrDir, ...
%                 [sprintf('correspondencePath_r%s_c%s', refExptID, cExptID), ...
%                 extn '.mat']) ;
% 
%             if ~exist(corrPathFn, 'file') || overwrite || forceTrue
%                 clf
%                 imagesc(timestamps(cc, 1:ntps(cc))/60, timestamps(1, :)/60, ssrM') ;
%                 axis equal
%                 axis tight
%                 cb = colorbar() ;
%                 ylabel(cb, '$\sqrt{\langle \sigma_{ij} \rangle \langle \sigma_{ji} \rangle}$', ...
%                     'interpreter', 'latex')
%                 xlabel(['time ' cExptID ' [hr]'], 'interpreter', 'latex')
%                 ylabel(['time ' refExptID ' [hr]'], 'interpreter', 'latex')
%                 saveas(gcf, fullfile( ssrDir, ...
%                     sprintf(['ssr_heatmap_c%s_r%s' extn '.png'], cExptID, refExptID)))
%                 save(fullfile( ssrDir, ...
%                     sprintf(['ssr_c%s_r%s' extn '.mat'], cExptID, refExptID)), ...
%                     'ssrM')
% 
% 
%                 % Get shortest path
%                 % Define correspondence pairs like in dynamicAtlas
%                 ssr4path = imgaussfilt(-ssrM, sigmaTime / rsubsampling) ;
%                 ssr4path = ssr4path - min(ssr4path(:)) ;
%                 pathOpts = struct('exponent', 1.0) ;
%                 corrPath = shortestPathInImage(ssr4path, pathOpts) ;
% 
%                 % Interpolate to get measure of the ssR along the path
%                 % imagesc(timestamps(cc, 1:ntps(cc)), timestamps(1,:), ssrM')
%                 [tc, tr] = meshgrid(timestamps(cc, 1:ntps(cc))-timestamps(cc,1), ...
%                     timestamps(1,:)- timestamps(1, 1)) ;
%                 ssrPath = interp2(tc, tr, ssrM', corrPath(:, 1), corrPath(:, 2), 'linear') ;
% 
% 
%                 % Save corrPaths and save image
%                 corrRaw = minddssr(cc, 1:ntps(cc)) ;
%                 corrError = movmean(minerror(cc, 1:ntps(cc)), 5) ;
%                 % Check if we need to truncate the correpondence path at start
%                 % or end
%                 if length(corrPath) < length(corrRaw)
%                     if corrPath(1, 1) == 1
%                         rawID = 1:length(corrPath) ;
%                         corrRaw = corrRaw(rawID) ;
%                         corrError = corrError(rawID) ;
%                     else
%                         rawID = (length(corrRaw)-length(corrPath)):length(corrRaw) ;
%                         corrRaw = corrRaw(rawID) ;
%                         corrError = corrError(rawID) ;
%                     end
%                 else
%                     rawID = 1:length(corrRaw) ;
%                 end
% 
%                 % Estimate uncertainty by difference between ssR at raw and
%                 % smoothed paths
%                 rawcorrSSRs = [] ;
%                 for ptID = 1:length(corrRaw)
%                     rrowP = corrRaw(ptID) ;
%                     ccolP = rawID(ptID) ;
%                     rawcorrSSRs(ptID) = ssrM(ccolP,rrowP) ;
%                 end
%                 ssrPathError = abs(ssrPath - rawcorrSSRs') ;
% 
%                 % Visualization
%                 clf
%                 imagesc(ssrM)
%                 hold on;
%                 plot(corrPath(:, 2), corrPath(:, 1), 'o') ;
%                 plot(corrRaw, rawID, 'k.') 
%                 plot(corrPath(:, 2) -corrError', corrPath(:, 1), 'k-')
%                 plot(corrPath(:, 2) +corrError', corrPath(:, 1), 'k-')
%                 cb = colorbar() ;
%                 ylabel(cb, '$\langle \sigma_{ij} \rangle \langle \sigma_{ji} \rangle$', ...
%                     'interpreter', 'latex')  
%                 colormap(viridis_r)
%                 axis equal 
%                 axis tight
%                 xlabel(['time ' refExptID ' [min]'], 'interpreter', 'latex')
%                 ylabel(['time ' cExptID ' [min]'], 'interpreter', 'latex')
%                 corrPathFigFn = fullfile(ssrDir, ...
%                     sprintf(['correspondencePath_r%s_c%s' extn '.pdf'], ...
%                     refExptID, cExptID)) ;
%                 saveas(gcf, corrPathFigFn)
% 
%                 % SAVE DATA RESULT
%                 save(corrPathFn, 'corrPath', 'corrRaw', 'corrError', 'ssrPath', 'ssrPathError') ;
%             else
%                 load(corrPathFn, 'corrPath', 'corrRaw', 'corrError', 'ssrPath', 'ssrPathError') ;
%             end
%             corrPaths{cc} = corrPath ;
%             corrErrors{cc} = corrError ;
%             ssrPaths{cc} = ssrPath ; 
%             ssrPathErrors{cc} = ssrPathError ;

        end

        disp('Saving minddssr, minname, minweights, minerror, numpts')
        ssStr = sprintf('_ref%s_rsub%03d', refExptID, rsubsampling) ;
%         save(fullfile(icpDir, ['corrPaths' ssStr '.mat']), 'corrPaths', ...
%             'corrErrors', 'ssrPaths', 'ssrPathErrors') ;
        save(fullfile(icpDir, ['ssr_minimum' ssStr '.mat']), 'ssr_minimum');
        save(fullfile(icpDir, ['minddssr' ssStr '.mat']), 'minddssr');
        save(fullfile(icpDir, ['minname' ssStr '.mat']), 'minname');
        save(fullfile(icpDir, ['minweights' ssStr '.mat']), 'minweights');
        save(fullfile(icpDir, ['minerror' ssStr '.mat']), 'minerror');
        save(fullfile(icpDir, ['numpts1' ssStr '.mat']), 'numpts1');
        save(fullfile(icpDir, ['numpts2' ssStr '.mat']), 'numpts2');
    else
        % load results
        disp('Loading minssr results from disk')
        ssStr = sprintf('_ref%s_rsub%03d', refExptID, rsubsampling) ;
        load(fullfile(icpDir, ['corrPaths' ssStr '.mat']), 'corrPaths', ...
            'corrErrors', 'ssrPaths', 'ssrPathErrors') ;
        load(fullfile(icpDir, ['ssr_minimum' ssStr '.mat']), 'ssr_minimum');
        load(fullfile(icpDir, ['minddssr' ssStr '.mat']), 'minddssr');
        load(fullfile(icpDir, ['minname' ssStr '.mat']), 'minname');
        load(fullfile(icpDir, ['minweights' ssStr '.mat']), 'minweights');
        load(fullfile(icpDir, ['minerror' ssStr '.mat']), 'minerror');
        load(fullfile(icpDir, ['numpts1' ssStr '.mat']), 'numpts1');
        load(fullfile(icpDir, ['numpts2' ssStr '.mat']), 'numpts2');

    end
    
    % collect all the info for this reference dataset refID
    corrPathCollection{refID} = corrPaths ;
    corrErrorsCollection{refID} = corrErrors ;
    minddssrCollection{refID} = minddssr ;
    numpts1Collection{refID} = numpts1 ;
    numpts2Collection{refID} = numpts2 ;
end

% What is the sampling density for each surface?
clf
for cc = 1:ndatasets
    length(find(numpts2(cc, :)))
    % plot(numpts1(cc, :))
    try
        plot(areas{cc} ./ numpts2(cc, find(numpts2(cc, :)))')
    catch
        disp(['could not plot this area for cc=' num2str(cc)])
    end
    hold on;
end
ylabel('area per sample point')
xlabel('time')

%% Plot an example for a figure -- dynamicAtlas main text figure
colors = define_colors ;
cmp = twilight_shifted(256) ;
cmp = cmp(length(cmp)*0.5+1:end, :) ;
for cc = 1:ndatasets
    % Load residuals
    ssrccDir = fullfile(ssrDir, ...
        sprintf(['dataset_' fluors{cc} '_' datestamps{cc}])) ;
    if ~exist(ssrccDir, 'dir')
        mkdir(ssrccDir)
    end

    % Load the areas if we use them to do SSR normalization. Otherwise,
    % use length of point cloud.
    area = areas{cc} ;
    lastTP = ntps(cc) ;
    assert(length(area) == lastTP)

    % Consider each dataset TP (t_c) and match to reference mesh
    ssrMatrix = zeros(lastTP, length(clist)) ;
    for ii = clist(clist < (lastTP + 1)) % index of dataset thats being parsed against CAAX
        % Consider the mesh only if mca{i,c} is populated with
        % struct that has field 'name'.
        if isfield(mca{cc, ii}, 'name') == 1

            % Check if SSR(cc, ii, rlist) has already been saved
            ssrii_fn = fullfile(ssrccDir, ...
                sprintf('ssr_tp%04d_rsub%03d_ptsub%03d.mat', ii, ...
                rsubsampling, ssample_factor)) ;

            tmp = load(ssrii_fn) ;
            ssrMatrix(ii, :) = tmp.ssr ;
        end
    end

    % Load path
    corrPathFn = fullfile(ssrDir, ...
        [sprintf('correspondencePath_r%02dc%02d', refID, cc), ...
        extn '.mat']) ;
    tmp = load(corrPathFn) ;
    corrPath = tmp.corrPath ; % path of correspondence
    corrError = tmp.corrError ; % unc in time of path
    ssrPath = tmp.ssrPath ; % ssr along the path against reference
    ssrError = tmp.ssrPathError ; % uncertainty in ssr along the path against reference
    
    % Plot it
    close all
    figure('Position', [0 0 160 140], 'Units', 'pixels')
    imagesc((1:size(ssrMatrix, 1))/60, ...
        (1:size(ssrMatrix,2))/60, ssrMatrix')
    axis equal
    axis tight
    caxis([0, Inf]) ;
    colormap(cmp)
    hold on;
    plot(corrPath(:, 1)/60, corrPath(:, 2)/60, '-', 'linewidth', 2, 'color', colors(3, :))
    lineProps = {'-','color', colors(3, :)} ;
    shadedErrorBar(corrPath(:, 1)/60, corrPath(:, 2)/60, corrError/60, 'lineprops', lineProps)
    cb = colorbar ;
    labelstr = ['$\sqrt{d^2_{' num2str(cc) ',' num2str(refID) '}'] ;
    labelstr = [labelstr 'd^2_{' num2str(refID) ',' num2str(cc) '}}$'] ;
    ylabel(cb, labelstr, 'interpreter', 'latex')
    xlabel(['emyro ' num2str(cc) ' timeline [hr]'])
    ylabel(['reference embryo ' num2str(refID) ' timeline [hr]'])
    title(mca{cc, refID}.folder)
    fn = fullfile(icpDir, ...
        [sprintf('correspondencePath_r%02dc%02d', refID, cc), ...
        extn '.pdf']) ;
    disp(['saving ' fn])
    saveas(gcf, fn)
    
    colormap(flipud(cmp))
    
    fn = fullfile(icpDir, ...
        [sprintf('correspondencePath_r%02dc%02d', refID, cc), ...
        extn '_flipcolormap.pdf']) ;
    disp(['saving ' fn])
    saveas(gcf, fn)
end

%% Align to global timeline
[expts, exptIDs] = unpackDataMap(dmap) ;

% % KEEP ONLY 1-6 right now...
% expts2 = {} ; exptIDs2 = {} ;
% for pp = 1:6
%     expts2{pp} = expts{pp} ;
%     exptIDs2{pp} = exptIDs{pp} ;
% end
% expts= expts2 ;
% exptIDs = exptIDs2 ;

hard = 1 ;
timelineDir = fullfile(icpDir, 'timeline') ;
if ~exist(timelineDir, 'dir')
    mkdir(timelineDir)
end
cpathFnBase = fullfile(icpDir, 'ssr', ['correspondencePath_r%02dc%02d' extn '.mat']) ;
[ttc, expts, exptIDs] = ...
    buildTimeTimeCorrespondences(expts, exptIDs, hard, ...
    timelineDir, extn, cpathFnBase, overwrite) ;

% Report average rate across the datasets
for pp = 1:length(expts)
    rates(pp) = (ttc{1}{pp}(end, 1) - ttc{1}{pp}(1, 1)) / ...
        (ttc{1}{pp}(end, 2) - ttc{1}{pp}(1, 2)) ;
end
disp(['rates = ', num2str(rates)])

relaxPairwiseCorrespondenceNetwork(ttc, hard, ...
    expts, exptIDs, timelineDir, 'ssR', overwrite)


%% Plotting "metabolism" (relative development rates) using ICP
[color, colornames] = define_colors() ;
% note: previously defined colors by hand
% color = {[217 54 104]/255, [144 115 50]/255, [76 133 50]/255,...
%     [54 131 123]/255, [59 125 171]/255, ...
%     [186 64 204]/255, [217 54 104]/255};
close all
figure('visible', 'off')
shape = define_markers(ndatasets) ;
% note: previously defined markers by hand
% shape = {'o', 'square', '^', 'V', 'd', 'p'} ;

% Plot the curves
hold on
leg = {};
% Blanket zero-valued minddssr as NaN
minddssr(minddssr<1e-5) = nan;
pss = 1 ; % plot subsampling
% Plot each dataset's correspondence times
ICPt = zeros(ndatasets, length(mca)) ;
slopes = zeros(ndatasets, 1) ;
for cc = 1:ndatasets
    % reference times -- note most ctimes are NaN if subsampled
    % Here this is really just a long list of indices
    rr = 1:max_ntp ;
    ctime = minddssr(cc, :);
    % now, rr becomes just the current dataset's time (cc's timepoints)
    rr = rr(~isnan(ctime));
    y = ctime(~isnan(ctime));
    yerror = minerror(cc, :);
    yerror = yerror(~isnan(ctime));
    
    if rsubsampling < ssfactorMedianThres
        y = movmedian(y,3) ;
    end
    
    % Truncate first few and last few timepoints if any timematches 
    % are before t=ignoreFirst or after t=ignoreLast in units of reference
    % timepoints.
    % Get LAST timepoint that hits min for domain of fit
    [~, fitmin] = min(abs(ignoreFirst - flipud(y(:)))) ;
    fitmin = numel(y) - fitmin + 1 ;
    % Get FIRST timepoint that hits max for domain of fit
    [~, fitmax] = min(abs(max(rlist)-ignoreLast-y(:))) ;
    
    %p = polyfit(rr(fitmin:fitmax), y(fitmin:fitmax), 1) ;
    w = minweights(cc, :) ;
    w = w(~isnan(ctime)) ;
    linfun = @(a,rr) a(1).*rr + a(2) ;
    p = nlinfit(rr(fitmin:fitmax) - t0ref, ...
        y(fitmin:fitmax) - t0ref, linfun, [0 0],...
        'Weights', w(fitmin:fitmax)) ;
    pfit = p(1)*(rr-t0ref)+p(2);
    xshift = p(2)/p(1) ;
    rr = rr + xshift ;
    %plot(rr(1:pss:end), y(1:pss:end), shape{cc}, 'Color', color(cc, :));
    sh = shadedErrorBar(rr(1:pss:end)- t0ref, y(1:pss:end)- t0ref, yerror(1:pss:end),...
        'lineprops', '-', 'lineprops', {'color', color(cc,:)}) ;
    plot(rr - t0ref, pfit, '--', 'Color', color(cc,:));
    leg{length(leg)+1} = labels{cc} ;
    % Fit to a quadratic profile
    leg{length(leg)+1} = ['Slope = ', num2str(p(1),'%10.2f')];
    lgd = legend(leg);
    for tt = 1:length(rlist)
        if any(tt>pfit) && any(tt<pfit)
            [~, ix] = min(abs(tt-pfit)) ;
            ICPt(cc,tt) = ix ;
            % plot(1:length(rlist), pfit)
        end
    end
    xlim([min(rlist)-ignoreFirst, max(rlist)-t0ref-ignoreLast]) ;
    ylim([min(rlist)-ignoreFirst, max(rlist)-t0ref-ignoreLast]) ;
    slopes(cc) = p(1) ;
end
clearvars y r p yy ii cc

lgd = legend(leg, 'Location', 'NorthEastOutside');
ax.XDir = 'normal' ;
ax.YDir = 'normal' ;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 20]);
title('Relative Rate of Development', 'fontsize', fs, 'Interpreter', 'Latex')
xlabel('Time [min]', 'fontsize', fs, 'Interpreter', 'Latex');
ylabel('Morphological time [min]', 'fontsize', fs, 'Interpreter', 'Latex');
pbaspect([1,1,1])
axis equal
lgd.FontSize=fs;
% Save to disk
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 figWidth figHeight]); 
figfn = fullfile(icpDir, 'rate_of_development') ;
disp(['Saving figure: ' figfn])
saveas(gcf, [figfn '.png']) ;
saveas(gcf, [figfn '.pdf']) ;
clearvars lgd ax


%% RMS of slopes
std(slopes)

%% Plotting "metabolism" (relative development rates) using ICP corrPaths
[color, colornames] = define_colors() ;
% note: previously defined colors by hand
% color = {[217 54 104]/255, [144 115 50]/255, [76 133 50]/255,...
%     [54 131 123]/255, [59 125 171]/255, ...
%     [186 64 204]/255, [217 54 104]/255};
close all
figure('visible', 'off')
ssStr = sprintf('_rsub%03d', rsubsampling) ;
% save(fullfile(icpDir, ['corrPaths' ssStr]), 'corrPaths', 'corrErrors', 'ssrPaths') ;
shape = define_markers(ndatasets) ;
% note: previously defined markers by hand
% shape = {'o', 'square', '^', 'V', 'd', 'p'} ;

% Plot the curves
hold on
leg = {};
% Blanket zero-valued minddssr as NaN
minddssr(minddssr<1e-5) = nan;
pss = 1 ; % plot subsampling
% Plot each dataset's correspondence times
ICPt = zeros(length(ndatasets), length(mca)) ;
for cc = 1:ndatasets
    % Load correspondences
    ctime = corrPaths{cc} ;
    % now, rr becomes just the current dataset's time (cc's timepoints)
    tt = ctime(:, 1) ;
    y = ctime(:, 2);
    yerror = corrErrors{cc} ;
    
    % if rsubsampling < ssfactorMedianThres
    %     y = movmedian(y,3) ;
    % end
    
    % Truncate first few and last few timepoints if any timematches 
    % are before t=ignoreFirst or after t=ignoreLast in units of reference
    % timepoints.
    % Get LAST timepoint that hits min for domain of fit
    [~, fitmin] = min(abs(ignoreFirst - flipud(y(:)))) ;
    fitmin = numel(y) - fitmin + 1 ;
    % Get FIRST timepoint that hits max for domain of fit
    [~, fitmax] = min(abs(max(rlist)-ignoreLast-y(:))) ;
    
    % Plot weighted fit
    w = minweights(cc, :) ;
    w = w(~isnan(ctime(:, 1))) ;
    w = w(:) ;
    linfun = @(a,xxx) a(1).*xxx + a(2) ;
    p = nlinfit(tt(fitmin:fitmax) - t0ref, ...
        y(fitmin:fitmax) - t0ref, linfun, [0 0],...
        'Weights', w(fitmin:fitmax)) ;
    pfit = p(1)*(tt-t0ref)+p(2);
    xshift = p(2)/p(1) ;
    tt = tt + xshift ;
    
    % Option 1: use spline for intercept offset
    % tref_c = interp1(y-t0ref, rr, 0, 'linear', 'extrap');
    % sh = shadedErrorBar((rr(1:pss:end)-tref_c) / 60,...
    %     (y(1:pss:end)- t0ref)/60, ...
    %     yerror(1:pss:end) / 60,...
    %     'lineprops', '-', 'lineprops', {'color', color(cc,:)}) ;
    
    % Option 2: use fit for intercept offset
    tref_c = interp1(pfit, tt, 0, 'linear', 'extrap');
    sh = shadedErrorBar((tt(1:pss:end)-tref_c) / 60,...
        (y(1:pss:end)- t0ref)/60, ...
        yerror(1:pss:end) / 60,...
        'lineprops', '-', 'lineprops', {'color', color(cc,:)}) ;
    
    % Plot fit
    plot((tt-t0ref)/60, pfit / 60, '--', 'Color', color(cc,:));
    
    
    leg{length(leg)+1} = labels{cc} ;
    % Fit to a quadratic profile
    leg{length(leg)+1} = ['Slope = ', num2str(p(1),'%10.2f')];
    for tt = 1:length(rlist)
        if any(tt>pfit) && any(tt<pfit)
            [~, ix] = min(abs(tt-pfit)) ;
            ICPt(cc,tt) = ix ;
            % plot(1:length(rlist), pfit)
        end
    end
    xlim([(min(rlist)-ignoreFirst)/60, (max(rlist)-t0ref-ignoreLast)/60]) ;
    ylim([(min(rlist)-ignoreFirst)/60, (max(rlist)-t0ref-ignoreLast)/60]) ;
end
clearvars y r p yy ii cc

% lgd = legend(leg, 'Location', 'NorthEastOutside');
ax.XDir = 'normal' ;
ax.YDir = 'normal' ;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 20]);
title('Relative Rate of Development', 'fontsize', fs, 'Interpreter', 'Latex')
xlabel('time [hr]', 'fontsize', fs, 'Interpreter', 'Latex');
ylabel('morphological time [hr]', 'fontsize', fs, 'Interpreter', 'Latex');
pbaspect([1,1,1])
xlim([-0.1, 1.8])
ylim([-0.1, 1.8])
xticks( [0, 0.5, 1, 1.5])
yticks( [0, 0.5, 1, 1.5])
axis equal
lgd.FontSize=fs;
% Save to disk
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figWidth figHeight]); 
figfn = fullfile(icpDir, 'rate_of_development_corrPaths') ;
disp(['Saving figure: ' figfn])
saveas(gcf, [figfn '.png']) ;
saveas(gcf, [figfn '.pdf']) ;
clearvars lgd ax


%% Plot residual
close all
hold on
leg = {};
% Plot each dataset's ssr over morphological time
first = true ;
for cc = setdiff(1:ndatasets, [refID])
    % Load correspondences
    ctime = corrPaths{cc} ;
    % now, rr becomes just the current dataset's time (cc's timepoints)
    tt = ctime(:, 1) ;
    tau = ctime(:, 2);
    tau_error = corrErrors{cc} ;
    ssrPath = sqrt(ssrPaths{cc}) ;
    ssrPathSm = movmedian(ssrPath, 7, 'omitnan')
    % apply chain rule for uncertainty: y=sqrt(x) --> dy = dx / [2sqrt(x)]
    ssrError = max(0.1, ssrPathErrors{cc} ./ (2 * sqrt(ssrPaths{cc}))) ;
    ssrErrorSm = movmedian(ssrError, 7, 'omitnan') ;
    
    
    
    if first
        vbinSSR = [(tau-t0ref)/60, ssrPathSm, ssrErrorSm] ;
        first = false ;
    else
        dat2add = [(tau-t0ref)/60, ssrPathSm, ssrErrorSm] ;
        vbinSSR = [vbinSSR; dat2add] ;
    end
    
    sh = shadedErrorBar((tau - t0ref)/ 60,...
        ssrPathSm, ssrErrorSm,...
        'lineprops', '-', 'lineprops', {'color', color(cc,:)}) ;
    
    leg{length(leg)+1} = labels{cc} ;
    xlim([(min(rlist)-ignoreFirst)/60, (max(rlist)-t0ref-ignoreLast)/60]) ;
    % ylim([(min(rlist)-ignoreFirst)/60, (max(rlist)-t0ref-ignoreLast)/60]) ;
end
clearvars y r p yy ii cc

lgd = legend(leg, 'Location', 'NorthEastOutside');
ax.XDir = 'normal' ;
ax.YDir = 'normal' ;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 20]);
title('Relative Rate of Development', 'fontsize', fs, 'Interpreter', 'Latex')
xlabel('time [hr]', 'fontsize', fs, 'Interpreter', 'Latex');
ylabel('difference between organ shapes [$\mu$m]', ...
    'fontsize', fs, 'Interpreter', 'Latex');
pbaspect([1,1,1])
axis equal
lgd.FontSize=fs;
% Save to disk
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figWidth figHeight]); 
figfn = fullfile(icpDir, 'ssr_between_corrPaths') ;
disp(['Saving figure: ' figfn])
saveas(gcf, [figfn '.png']) ;
saveas(gcf, [figfn '.pdf']) ;
clearvars lgd ax


%% Now plot the binstats
[midx, meany, stdy, ny, stey] = ...
    binDataMeanStdWeighted(vbinSSR(:, 1), vbinSSR(:, 2), ...
    [-0.5:2/60:2], 1./(vbinSSR(:, 3).^2), 100, vbinSSR(:, 3))

close all
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'position', [0, 0, 250, 200]) ;
set(gcf, 'PaperPosition', [0 0 figWidth*0.5 figHeight]); 
sh = shadedErrorBar(midx,...
    meany, stdy,...
    'lineprops', '-') ;
xlim([-0.1, 1.8])
xlabel('morphological time [a.u.]', 'interpreter', 'latex')
ylabel('difference between organ shapes [$\mu$m]', ...
    'fontsize', fs, 'Interpreter', 'Latex');
figfn = fullfile(icpDir, 'ssr_between_corrPathsBinStats') ;
disp(['Saving figure: ' figfn])
saveas(gcf, [figfn '.png']) ;
saveas(gcf, [figfn '.pdf']) ;



%% Plotting Overlays using Morphological Time Plot with corrPaths
close all

% TODO: make it clear which ttype is which (ie what is ttype == 3)
mintransDir = fullfile(icpDir, 'MinTransformCorrPath') ;
tp1transDir = fullfile(icpDir, 'FirstTPTransformCorrPath') ;
tavgtransDir = fullfile(icpDir, 'TimeAvgTransformCorrPath') ;
vartransDir = fullfile(icpDir, 'VariableTransformCorrPath') ;
quatDir = fullfile(icpDir, 'QuaternionTransformCorrPath') ;
% make all directories if they don't exist
dirs = {mintransDir, tp1transDir, tavgtransDir, vartransDir, quatDir} ;
for ii = 1:length(dirs)
    if ~exist(dirs{ii}, 'dir')
        mkdir(dirs{ii})
    end
end

transform_meshes = true ;
ttypes = {'pcregistericp', '', '', '', 'quaternion'} ;
for ttype = 1
    if ttype == 5 %quaternion method
        clear tform
        Q = [] ; % initiate matrix of quaternions
        % Load timevals
        load(fullfile(icpDir, 'tvals.mat')) ;
        for cc = 1:ndatasets
            for rr = 1:max_ntp
                clear x; clear y; clear z; clear w;
                syms x y z w ;
                % linearly index each tform value into an equation to solve
                eqn1 = 1 - 2*y^2 - 2*z^2 == tvals(cc,rr,1,1) ;
                eqn2 = 2*x*y + 2*z*w == tvals(cc,rr,2,1) ;
                eqn3 = 2*x*z - 2*y*w == tvals(cc,rr,3,1) ;
                eqn4 = 2*x*y-2*z*w == tvals(cc,rr,1,2) ;
                % eqn5 = 1 - 2*x^2 - 2*z^2 == tvals(cc,rr,2,2) ;
                % eqn6 = 2*y*z + 2*x*w == tvals(cc,rr,3,2) ;
                % eqn7 = 2*z*x + 2*y*w == tvals(cc,rr,1,3) ;
                % eqn8 = 2*y*z - 2*x*w == tvals(cc,rr,2,3) ;
                % eqn9 = 1 - 2*x^2 - 2*y^2 == tvals(cc,rr,3,3) ;
                % solve for components of the quaternion for c and r
                sols = solve(eqn1, eqn2, eqn3, eqn4) ;
                % build matrix of quaternions
                Q(1,rr) = sols.w(1) ;
                Q(2,rr) = sols.x(1) ;
                Q(3,rr) = sols.y(1) ;
                Q(4,rr) = sols.z(1) ;
            end
            Q2 = Q*Q' ; % multiply Q by its transpose
            [V,D] = eig(Q2) ; % find all eigenvalues of Q2 and their corresponding eigenvectors
            linind = find(D == max(D,[],'all')) ; % find linear index of max eigenvalue
            [m,n] = ind2sub([size(D,1) size(D,2)], linind) ; % turn that linear index into a 2d index
            maxev = V(:,n) ; % use column of max eigenvalue to find corresponding eigenvector
            q = maxev/norm(maxev) ; % normalize the averaged quaternion
            % assign component values and rebuild rotation matrix
            tformq(1,1,cc) = 1 - 2*q(3)^2 - 2*q(4)^2 ;
            tformq(1,2,cc) = 2*q(2)*q(3) + 2*q(4)*q(1) ;
            tformq(1,3,cc) = 2*q(2)*q(4) - 2*q(3)*q(4) ;
            tformq(2,1,cc) = 2*q(2)*q(3)-2*q(4)*q(1) ;
            tformq(2,2,cc) = 1 - 2*q(2)^2 - 2*q(4)^2 ;
            tformq(2,3,cc) = 2*q(3)*q(4) + 2*q(2)*q(1) ;
            tformq(3,1,cc) = 2*q(4)*q(2) + 2*q(3)*q(1) ;
            tformq(3,2,cc) = 2*q(3)*q(4) - 2*q(2)*q(1) ;
            tformq(3,3,cc) = 1 - 2*q(2)^2 - 2*q(3)^2 ;
            % finally, take mean of translation vector!
            tformq(4,1,cc) = mean(tvals(cc,:,4,1),2) ;
            tformq(4,2,cc) = mean(tvals(cc,:,4,2),2) ;
            tformq(4,3,cc) = mean(tvals(cc,:,4,3),2) ;
            tformq(:,4,cc) = 0 ;
            tformq(4,4,cc) = 1 ;
        end
    end
    
    % Consider all reference timepoints that have been matched. Plot this
    % rr reference mesh and all other meshes matched to this reference time
    rr2do = 113:30:size(mca, 2) ;
    rr2do = [rr2do, setdiff(1:size(mca, 2), rr2do)] ;
    for rr = rr2do
        disp(['rr = ' num2str(rr)])
        leg = {};
        leg{1} = labels{1};
        close all
        fig = figure('visible', 'off') ;
        
        hold on
        % First plot the reference mesh (c==1)
        refCloud = pcread(fullfile(mca{1, rr}.folder, mca{1, rr}.name));
        refxyz = refCloud.Location;
        refxoff = min(refxyz(:,1)) ;
        refxyz(:,1) = refxyz(:,1) - refxoff ;
        %refMesh = read_ply_mod(fullfile(mca{1, rr}.folder, mca{1, rr}.name)) ;
        
        drefxy = double(refxyz(:,1:2));
        %aShape = alphaShape(drefxy,5);
        %bnd = boundaryFacets(aShape) ;
        %bnd = [bnd(:, 1); bnd(1,1)] ;
        
        % Note: the third argument is the tightness parameter (S=1 is 100%
        % tight, while S=0 is convex hull)
        bnd = boundary(drefxy(:,1),drefxy(:,2), 1) ;
        plot(drefxy(bnd, 1), drefxy(bnd,2)) ;
        
        %trisurf(refMesh.f, refMesh.v(:, 1), refMesh.v(:, 2), refMesh.v(:, 3), ...
        %   'FaceColor', color{c}, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
        
        % Now consider all datasets (c~=1) and find matches to this
        % reference time rr
        found_any_match = false ;
        for cc = 2:ndatasets
            
            % Load correspondences
            ctime = corrPaths{cc} ;
            % ctime(:, 1) is just the current dataset's time (cc's timepoints)
            yy = ctime(:, 2);
            [minDT, matchID] = min(abs(yy - rr)) ;
            yerror = corrErrors{cc} ;

            % If closest match is within uncertainty, plot it!
            if minDT < max(1, yerror(matchID))
                found_any_match = true ;
                matchnum = round(matchID) ;
                cCloud = pcread(fullfile(mca{cc, matchnum}.folder, mca{cc, matchnum}.name));
                cxyz = cCloud.Location ;
                
                % CHANGE THE WAY YOU TRANSFORM GIVEN 'ttype'
                % ------------------------------------------
                if ttype == 1
                    clear tform
                    tform1 = pcregistericp(cCloud,refCloud,'Extrapolate',true); % extracting the 4x4 transform from cmesh to refmesh
                    tvals(cc,rr,:,:) = tform1.T;
                    tform = tform1.T ;
                    pn = fullfile(mintransDir, ['Compare_TP',sprintf('%d',rr),'_MinTransform.png']) ;
                    pn_pdf = fullfile(mintransDir, ['Compare_TP',sprintf('%d',rr),'_MinTransform.pdf']) ;
                    
                    % Save the time matching here
                    save(fullfile(icpDir, 'tvals'), 'tvals');
                    
                elseif ttype == 2
                    clear tform
                    load '/mnt/data/analysis/ICP_Plots/tvals.mat';
                    tform = squeeze(tvals(cc,1,:,:)); % 4x4 transform corresponding to first TP
                    pn = fullfile(tp1transDir, ['Compare_TP',sprintf('%d',rr),'_FirstTPTransform.png']) ;
                elseif ttype == 3
                    clear tform
                    load '/mnt/data/analysis/ICP_Plots/tvals.mat';
                    for m = 1:4
                        for n = 1:4
                            tform(m,n) = squeeze(mean(tvals(cc,:,m,n), 2)); % 4x4 transform corresponding to (elementwise) average of all time for a certain c
                        end
                    end
                    pn = fullfile(tavgtransDir, ['Compare_TP',sprintf('%d',rr),'_TimeAvgTransform.png']) ;
                elseif ttype == 4
                    clear tform
                    load '/mnt/data/analysis/ICP_Plots/tvals.mat';
                    TP = 60;
                    tform = squeeze(tvals(cc,TP,:,:)); % 4x4 transform corresponding to TP (change TP above to alter)
                    pn = fullfile(vartransDir, ['Compare_TP',sprintf('%d',rr),'_TP', sprintf('%d',TP), 'Transform.png']) ;
                elseif ttype == 5
                    tform = tformq(:,:,cc) ;
                    pn = fullfile(quatDir, ['Compare_TP',sprintf('%d',rr),'_QuaternionTransform.png']) ;
                end
                %--------------------------------------------------------------------------------------------------------------------
                
                if transform_meshes
                    cxyz = horzcat(cxyz, ones(length(cxyz),1)); % pad the 4th dim with ones to apply tform
                    cxyz = cxyz*tform; % apply the transform
                    cxyz = cxyz(:,1:3);
                end
                % AP offset
                % cxoff = min(cxyz(:,1)) ;
                % cxyz(:,1) = (cxyz(:,1)-cxoff) ;
                
                
                % Plot the transformed mesh surface
                dcxy = double(cxyz(:,1:2));
                %aShape = alphaShape(dcxy,5);
                %bnd = boundaryFacets(aShape) ;
                %bnd = [bnd(:, 1); bnd(1,1)] ;
                bnd = boundary(dcxy(:,1),dcxy(:,2), 1) ;
                plot(dcxy(bnd, 1), dcxy(bnd,2))
                %trisurf(cMesh.f, cxyz(:, 1), cxyz(:, 2), cxyz(:, 3), ...
                %  'FaceColor', color{c}, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
                
                leg{length(leg)+1} = labels{cc};
            end
        end

        % If there was at least one other mesh that matched in time, save
        % plot
        if found_any_match
            xlim([-30, 280])
            xlabel('AP position [$\mu$m]', 'Interpreter', 'latex')
            ylabel('lateral position [$\mu$m]', 'Interpreter', 'latex')
            axis equal
            ylim([-80, 80])
            lgd = legend(leg, 'Location', 'NorthEastOutside');
            set(gca, 'Position', [0.12 0.21 0.6 0.6]) 
            lgd.FontSize=8;
            title(['$t_\textrm{ref}=$', sprintf('%d',rr-t0ref), ' min'], 'Interpreter', 'Latex');
            disp(['Saving overlay figure to ', pn])
            saveas(fig,pn);
            saveas(fig,pn_pdf);
            clf
            hold off
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error('all done')



%% Using Residuals to Compute Chi Squared (Goodness of Fit)

load(fullfile(icpDir, 'numpts1')) ;
load(fullfile(icpDir, 'numpts2')) ;
load(fullfile(icpDir, 'minddssr')) ;
load(fullfile(icpDir, 'minname')) ;
load(fullfile(icpDir, 'minerror')) ;
ssrDir = fullfile(icpDir, 'ssr') ;

if exist(fullfile(icpDir, 'ChiSqs.mat'), 'file') && ~overwrite
    answer = questdlg('Chi Squareds already exists? Overwrite?') ;
else
    answer = 'Yes' ;
end

% Note that c=1 is the reference dataset, so cycle through all others
% starting with index=2
if strcmp(answer, 'Yes')
    % Calculate var of each TP combination
    for cc = 2:size(mca, 1)
        for rr = 1:size(mca, 2)
            if ~isempty(mca{cc,rr}) % check is cell array is populated
                cCloud = pcread(fullfile(mca{cc, rr}.folder, mca{cc, rr}.name));
                cxyz = cCloud.Location(1:ssample_factor:end,:) ;
                tree = KDTreeSearcher(cxyz);
                for ii = 1:size(cxyz,1)
                    [~,D] = knnsearch(tree, tree.X(ii,:), 'K', 6);
                    avg = mean(D(2:end));
                    avgs(ii) = avg;
                end
                std = mean(avgs) ;
                variance(cc,rr) = std^2 ;
                fprintf('c= %02d r= %03d\n', cc, rr);
            end
        end
    end
    save(fullfile(icpDir, 'variance'), 'variance');
    
    % Calculating Chi Squared for each considered timepoint as a function of
    % reference timepoints
    for cc = 2:ndatasets
        
        % define the SSR directory for this dataset
        ssrccDir = fullfile(ssrDir, ...
            sprintf(['dataset_' fluors{cc} '_' datestamps{cc}])) ;
        if ~exist(ssrccDir, 'dir')
            mkdir(ssrccDir)
        end
        
        % Load the areas if we use them to do SSR normalization. Otherwise,
        % use length of point cloud.
        area = areas{cc} ;
        lastTP = ntps(cc) ;
        assert(length(area) == lastTP)
        
        % Consider each dataset TP (t_c) and match to reference mesh
        for ii = clist(clist < (lastTP + 1)) % index of dataset thats being parsed against CAAX
            % Consider the mesh only if mca{i,c} is populated with
            % struct that has field 'name'.
            if isfield(mca{cc, ii}, 'name') == 1
                
                % Check if SSR(cc, ii, rlist) has already been saved
                ssrii_fn = fullfile(ssrccDir, sprintf('ssr_tp%04d_rsub%03d_ptsub%03d.mat', ii, rsubsampling, ssample_factor)) ;
                
                % Assume that we must redo calculation unless proven
                % that otherwise if file exists and overwrite == false
                redoii = true ;
                if exist(ssrii_fn, 'file') && overwrite
                    tmp = load(ssrii_fn, 'rlist', 'ssr') ;
                    % Check that the rlist is indeed the same
                    if all(tmp.rlist == rlist)
                        % Calculation exists for the same rlist for this ii
                        redoii = false ;
                        % assign ssr for this cc, ii
                        ssr = tmp.ssr ;
                    end
                end
                chi2 = ssr.^2/variance(cc,ii) ; % calculate Chi Squared
                ChiSqs{cc,ii} = chi2*sqrt(numpts1(cc,ii))*sqrt(numpts2(cc,ii)) ; % Normalize and store values in a cell array
            end
        end
        
    end
    save(fullfile(icpDir, 'ChiSqs'), 'ChiSqs');
end
disp('done with chisquared calculation')

%% Plotting and extracting match info from Chi Squareds

load(fullfile(icpDir, 'ChiSqs')) ;
load(fullfile(icpDir, 'numpts1')) ;

rlist = 1:max_ntp;
chisqDir = fullfile(icpDir, 'Residuals_chisq') ;
if ~exist(chisqDir, 'dir')
    mkdir(chisqDir)
end

% Cycle through ChiSqs Array
for cc = 2:size(ChiSqs, 1)
    for rr = 1:size(ChiSqs,2)
        if ~isempty(ChiSqs{cc,rr})
            clearvars xpoly ypoly
            Chis2 = ChiSqs{cc,rr} ; % extract Chis2 info
            if rsubsampling < ssfactorMedianThres
                Chis2 = movmedian(Chis2, 5); % smooth Chis2
            end
            [~,ix] = min(Chis2); % extract min of Chis2 = best fit
            %find polyfit of min to find error of match
            rlist = 1:max_ntp;
            if ix+fitrange > size(Chis2,2)
                maxend = size(Chis2,2);
            else
                maxend = ix + fitrange ;
            end
            if ix <= fitrange
                xpoly = rlist(1:maxend) ;
                ypoly = Chis2(1:maxend) ;
            elseif ix > fitrange
                xpoly = rlist(ix-fitrange:maxend) ;
                ypoly = Chis2(ix-fitrange:maxend) ;
            elseif ix + fitrange > length(rlist)
                xpoly = rlist(length(rlist)-2*fitrange:end) ;
                ypoly = Chis2(length(rlist)-2*fitrange:end) ;
            end
            %p = polyfit(xpoly,ypoly',2) ;
            pfix = polyfix(xpoly, ypoly', 2, ix, min(Chis2),ix, 0) ;
            p = polyval(pfix, round(linspace(ix-15, ix+15, 30)));
            % store minimum index as ChiSqsmin
            ChiSqsMin(cc, rr) = rlist(ix);
            
            % find error using this fit
            [~, xer] = min(abs(polyval(pfix,rlist)-Chis2(ix)-1)) ;
            ChiSqsError(cc,rr) = abs(ix-xer) ;
            
            fig = figure('visible', 'off','Position', [10 10 900 600]) ;
            plot(Chis2, '.--', 'Color', 'blue')
            hold on
            plot(round(linspace(ix-15, ix+15, 30)),p, 'color', 'magenta', 'linewidth', 1.5) ;
            title(['Goodness of Fit', sprintf(' Dataset %01d TP= %03d', cc, rr )]);
            xlabel('Reference Time [min]');
            ylabel('$$\chi^2$$', 'Interpreter', 'Latex');
            % ylim([0 max(Chis2)]) ;
            % xlim([0 length(rlist)]) ;
            outfn = fullfile(chisqDir, ...
                ['Chi2_Dataset_', num2str(cc),'_','TP_', num2str(rr), '.png']) ;
            disp(['Saving ' outfn]);
            saveas(fig, outfn);
            close all
        end
    end
end
ChiSqsError(isnan(ChiSqsError)) = 0 ;
save(fullfile(icpDir, 'ChiSqsMin'), 'ChiSqsMin');
save(fullfile(icpDir, 'ChiSqsError'), 'ChiSqsError');
%% Plotting "metabolism" (relative development rates) using Chi Squared
[color, ~] = define_colors() ;
% note: previously defined colors by hand
% color = {[217 54 104]/255, [144 115 50]/255, [76 133 50]/255,...
%     [54 131 123]/255, [59 125 171]/255, ...
%     [186 64 204]/255, [217 54 104]/255};
close all
load(fullfile(icpDir, 'ChiSqsMin')) ;
load(fullfile(icpDir, 'ChiSqsError')) ;
shape = define_markers(ndatasets) ;
% note: previously defined markers by hand
% shape = {'o', 'square', '^', 'V', 'd', 'p'} ;
ChiSqsMin(ChiSqsMin == 0) = nan ;
ChiSqsError(ChiSqsError == 0) = nan ;
% Plot the curves
hold on
leg = {};
rr = linspace(1, max_ntp, max_ntp);
plot(rr, rr, 'k-');
leg{1} = labels{1};
lgd = legend(leg, 'Location', 'NorthEastOutside');
% Blanket zero-valued minddssr as NaN
pss = 5 ; % plot subsampling
% Plot each dataset's correspondence times
for cc = 2:ndatasets
    ctime = ChiSqsMin(cc,:);
    rr = rr(~isnan(ctime));
    y = ctime(~isnan(ctime));
    yerror = ChiSqsError(cc, :);
    yerror = yerror(~isnan(ctime));
    y = movmedian(y,3) ;
    p = polyfit(rr, y, 1) ;
    %fit = p(1)*rr;
    y = y - p(2);
    plot(rr(1:pss:end), y(1:pss:end), shape{cc}, 'Color', color(cc, :));
    shadedErrorBar(rr(1:pss:end), y(1:pss:end), yerror(1:pss:end), 'lineprops', '-') ;
    % Fit to a quadratic profile
    %plot(rr, fit, '--', 'Color', color(cc,:));
    leg{length(leg)+1} = labels{cc} ;
    leg{length(leg)+1} = [labels{cc}, ' Error'] ;
    %leg{length(leg)+1} = ['Slope = ', num2str(p(1),'%10.2f')];
    lgd = legend(leg);
    rr = linspace(1, max_ntp, max_ntp);
end
clearvars y r p yy ii cc

ax.XDir = 'normal' ;
ax.YDir = 'normal' ;
title('Relative Rate of Development Chi2')
xlabel('Considered Mesh Timepoint [min]');
ylabel('Reference Mesh Timepoint [min]');
axis equal
pbaspect([1,1,1])
lgd.FontSize=10;
%saveas(gcf, fullfile(icpDir, 'rate_of_development.png')) ;
clearvars lgd ax



%% Plotting Overlays using Morphological Time Plot -- outlines via TSP travelling salesman
close all

% TODO: make it clear which ttype is which (ie what is ttype == 3)
mintransDir = fullfile(icpDir, 'MinTransform') ;
tp1transDir = fullfile(icpDir, 'FirstTPTransform') ;
tavgtransDir = fullfile(icpDir, 'TimeAvgTransform') ;
vartransDir = fullfile(icpDir, 'VariableTransform') ;
quatDir = fullfile(icpDir, 'QuaternionTransform') ;
% make all directories if they don't exist
dirs = {mintransDir, tp1transDir, tavgtransDir, vartransDir, quatDir} ;
for ii = 1:length(dirs)
    if ~exist(dirs{ii}, 'dir')
        mkdir(dirs{ii})
    end
end

ttypes = {'pcregistericp', '', '', '', 'quaternion'} ;
for ttype = 1
    if ttype == 5 %quaternion method
        clear tform
        Q = [] ; % initiate matrix of quaternions
        % Load timevals
        load(fullfile(icpDir, 'tvals.mat')) ;
        for cc = 1:ndatasets
            for rr = 1:max_ntp
                clear x; clear y; clear z; clear w;
                syms x y z w ;
                % linearly index each tform value into an equation to solve
                eqn1 = 1 - 2*y^2 - 2*z^2 == tvals(cc,rr,1,1) ;
                eqn2 = 2*x*y + 2*z*w == tvals(cc,rr,2,1) ;
                eqn3 = 2*x*z - 2*y*w == tvals(cc,rr,3,1) ;
                eqn4 = 2*x*y-2*z*w == tvals(cc,rr,1,2) ;
                eqn5 = 1 - 2*x^2 - 2*z^2 == tvals(cc,rr,2,2) ;
                eqn6 = 2*y*z + 2*x*w == tvals(cc,rr,3,2) ;
                eqn7 = 2*z*x + 2*y*w == tvals(cc,rr,1,3) ;
                eqn8 = 2*y*z - 2*x*w == tvals(cc,rr,2,3) ;
                eqn9 = 1 - 2*x^2 - 2*y^2 == tvals(cc,rr,3,3) ;
                % solve for components of the quaternion for c and r
                sols = solve(eqn1, eqn2, eqn3, eqn4) ;
                % build matrix of quaternions
                Q(1,rr) = sols.w(1) ;
                Q(2,rr) = sols.x(1) ;
                Q(3,rr) = sols.y(1) ;
                Q(4,rr) = sols.z(1) ;
            end
            Q2 = Q*Q' ; % multiply Q by its transpose
            [V,D] = eig(Q2) ; % find all eigenvalues of Q2 and their corresponding eigenvectors
            linind = find(D == max(D,[],'all')) ; % find linear index of max eigenvalue
            [m,n] = ind2sub([size(D,1) size(D,2)], linind) ; % turn that linear index into a 2d index
            maxev = V(:,n) ; % use column of max eigenvalue to find corresponding eigenvector
            q = maxev/norm(maxev) ; % normalize the averaged quaternion
            % assign component values and rebuild rotation matrix
            tformq(1,1,cc) = 1 - 2*q(3)^2 - 2*q(4)^2 ;
            tformq(1,2,cc) = 2*q(2)*q(3) + 2*q(4)*q(1) ;
            tformq(1,3,cc) = 2*q(2)*q(4) - 2*q(3)*q(4) ;
            tformq(2,1,cc) = 2*q(2)*q(3)-2*q(4)*q(1) ;
            tformq(2,2,cc) = 1 - 2*q(2)^2 - 2*q(4)^2 ;
            tformq(2,3,cc) = 2*q(3)*q(4) + 2*q(2)*q(1) ;
            tformq(3,1,cc) = 2*q(4)*q(2) + 2*q(3)*q(1) ;
            tformq(3,2,cc) = 2*q(3)*q(4) - 2*q(2)*q(1) ;
            tformq(3,3,cc) = 1 - 2*q(2)^2 - 2*q(3)^2 ;
            % finally, take mean of translation vector!
            tformq(4,1,cc) = mean(tvals(cc,:,4,1),2) ;
            tformq(4,2,cc) = mean(tvals(cc,:,4,2),2) ;
            tformq(4,3,cc) = mean(tvals(cc,:,4,3),2) ;
            tformq(:,4,cc) = 0 ;
            tformq(4,4,cc) = 1 ;
        end
    end
    
    % Consider all reference timepoints that have been matched. Plot this
    % rr reference mesh and all other meshes matched to this reference time
    rr2do = 1:20:size(mca, 2) ;
    rr2do = [rr2do, setdiff(1:size(mca, 2), rr2do)] ;
    for rr = rr2do
        disp(['rr = ' num2str(rr)])
        leg = {};
        leg{1} = labels{1};
        close all
        fig = figure('visible', 'off') ;
        
        hold on
        % First plot the reference mesh (c==1)
        refCloud = pcread(fullfile(mca{1, rr}.folder, mca{1, rr}.name));
        refxyz = refCloud.Location;
        refxoff = min(refxyz(:,1)) ;
        refxyz(:,1) = refxyz(:,1) - refxoff ;
        %refMesh = read_ply_mod(fullfile(mca{1, rr}.folder, mca{1, rr}.name)) ;
        
        drefxy = double(refxyz(:,1:2));
        %aShape = alphaShape(drefxy,5);
        %bnd = boundaryFacets(aShape) ;
        %bnd = [bnd(:, 1); bnd(1,1)] ;
        
        % Note: the third argument is the tightness parameter (S=1 is 100%
        % tight, while S=0 is convex hull)
        bnd = boundary(drefxy(:,1),drefxy(:,2), 1) ;
        plot(drefxy(bnd, 1), drefxy(bnd,2)) ;
        
        %trisurf(refMesh.f, refMesh.v(:, 1), refMesh.v(:, 2), refMesh.v(:, 3), ...
        %   'FaceColor', color{c}, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
        
        % Now consider all other datasets (c~=1) and find matches to this
        % reference time rr
        for cc = 2:ndatasets
            if ICPt(cc,rr) ~= 0
                matchnum = ICPt(cc,rr) ;
                cCloud = pcread(fullfile(mca{cc, matchnum}.folder, mca{cc, matchnum}.name));
                cxyz = cCloud.Location ;
                
                % CHANGE THE WAY YOU TRANSFORM GIVEN 'ttype'
                % ------------------------------------------
                if ttype == 1
                    clear tform
                    tform1 = pcregistericp(cCloud,refCloud,'Extrapolate',true); % extracting the 4x4 transform from cmesh to refmesh
                    for m = 1:4
                        for n = 1:4
                            tvals(cc,rr,m,n) = tform1.T(m,n);
                        end
                    end
                    tform = tform1.T ;
                    pn = fullfile(mintransDir, ['Compare_TP',sprintf('%d',rr),'_MinTransform.png']) ;
                    
                    % Save the time matching here
                    save(fullfile(icpDir, 'tvals'), 'tvals');
                    
                elseif ttype == 2
                    clear tform
                    load '/mnt/data/analysis/ICP_Plots/tvals.mat';
                    tform = squeeze(tvals(cc,1,:,:)); % 4x4 transform corresponding to first TP
                    pn = fullfile(tp1transDir, ['Compare_TP',sprintf('%d',rr),'_FirstTPTransform.png']) ;
                elseif ttype == 3
                    clear tform
                    load '/mnt/data/analysis/ICP_Plots/tvals.mat';
                    for m = 1:4
                        for n = 1:4
                            tform(m,n) = squeeze(mean(tvals(cc,:,m,n), 2)); % 4x4 transform corresponding to (elementwise) average of all time for a certain c
                        end
                    end
                    pn = fullfile(tavgtransDir, ['Compare_TP',sprintf('%d',rr),'_TimeAvgTransform.png']) ;
                elseif ttype == 4
                    clear tform
                    load '/mnt/data/analysis/ICP_Plots/tvals.mat';
                    TP = 60;
                    tform = squeeze(tvals(cc,TP,:,:)); % 4x4 transform corresponding to TP (change TP above to alter)
                    pn = fullfile(vartransDir, ['Compare_TP',sprintf('%d',rr),'_TP', sprintf('%d',TP), 'Transform.png']) ;
                elseif ttype == 5
                    tform = tformq(:,:,cc) ;
                    pn = fullfile(quatDir, ['Compare_TP',sprintf('%d',rr),'_QuaternionTransform.png']) ;
                end
                %--------------------------------------------------------------------------------------------------------------------
                
                cxyz = horzcat(cxyz, ones(length(cxyz),1)); % pad the 4th dim with ones to apply tform
                cxyz = cxyz*tform; % apply the transform
                cxyz = cxyz(:,1:3);
                cxoff = min(cxyz(:,1)) ;
                cxyz(:,1) = (cxyz(:,1)-cxoff) ;
                % Plot the transformed mesh surface
                dcxy = double(cxyz(:,1:2));
                %aShape = alphaShape(dcxy,5);
                %bnd = boundaryFacets(aShape) ;
                %bnd = [bnd(:, 1); bnd(1,1)] ;
                bnd = boundary(dcxy(:,1),dcxy(:,2), 1) ;
                plot(dcxy(bnd, 1), dcxy(bnd,2))
                %trisurf(cMesh.f, cxyz(:, 1), cxyz(:, 2), cxyz(:, 3), ...
                %  'FaceColor', color{c}, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
                
                leg{length(leg)+1} = labels{cc};
            end
        end

        xlim([-30, 280])
        xlabel('AP position [$\mu$m]', 'Interpreter', 'latex')
        ylabel('lateral position [$\mu$m]', 'Interpreter', 'latex')
        axis equal
        ylim([-80, 80])
        lgd = legend(leg, 'Location', 'NorthEastOutside');
        set(gca, 'Position', [0.12 0.21 0.6 0.6]) 
        lgd.FontSize=8;
        title(['$t_\textrm{ref}=$', sprintf('%d',rr), ' min'], 'Interpreter', 'Latex');
        disp(['Saving overlay figure to ', pn])
        saveas(fig,pn);
        clf
        hold off
    end
end