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
% Note: Previously this code had run over reftime, then ctime. Now ctime
% then reftime, so minddssr and other objects have been effectively
% 'transposed'.

%% Clean MATLAB and add paths
clear
close all
clc
addpath('/mnt/data/code/')
addpath('/mnt/data/code/gut_matlab/')
addpath('/mnt/data/code/gut_matlab/curve_functions/')
addpath('/mnt/data/code/gut_matlab/plotting')
addpath('/mnt/data/code/imsaneV1.2.3/generalfunctions')

outdir = '/mnt/data/analysis/' ;
cd(outdir)
icpDir = fullfile(outdir, 'ICP_Plots/') ;
if ~exist(icpDir, 'dir')
    mkdir(icpDir) ;
end
%% Global options
ssample_factor = 40 ;   % subsampling of point clouds for ICP
rsubsampling = 1 ;     % subsampling of timepoints for reference dataset
csubsampling = 1 ;     % subsampling of timepoints for matched dataset
overwrite = false ;      % overwrite previous results
if overwrite == 1
    disp('WARNING: OVERWRITE IS TURNED ON');
end
% Plotting
figWidth = 12 ;
figHeight = 9 ;
%% Create Cell Array with Mesh Metadata

dmap = buildLookupMap('/mnt/data/analysis/lookupMeta.txt'); % build the map that contains all data and its metadata
mca = {} ; % initiate mesh cell array

% Iterate over each marker
for mi = 1:length(keys(dmap))
    labels = keys(dmap) ;
    label = labels{mi} ;
    % Cycle through all datasets of this marker
    for jj=1:length(dmap(label).folders)
        
        % get col
        if strcmp(label, 'caax')
            col = 1 ;
        elseif strcmp(label, 'hrfp')
            col = length(dmap('caax').folders) + 1 ;
        elseif strcmp(label, 'la')
            col = length(dmap('caax').folders) + length(dmap('hrfp').folders) + 1 ;
        end
        
        meshes = dir(fullfile(dmap(label).folders{jj}, 'aligned_meshes',  'mesh_apical_stab_0*_APDV_um.ply')) ;
        
        % Write each mesh to the cell array mca
        for kk = 1:length(meshes)
            mca{col+jj-1, kk} = meshes(kk);
        end
    end
end
disp('done building mca')
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
        % update the labels index
        dmyk = dmyk + 1 ;
    end
end
% cleanup from labels generation
clearvars dmyk jjkey key datestamp kk

%% Declare the range of morphological time and datasets
rlist = 1:rsubsampling:max_ntp;
clist = 1:csubsampling:max_ntp;

%% Use ICP to align meshes in time
minddssr = zeros(ndatasets, max_ntp); % initiate matrix that will store SSR data for each dataset

if exist(fullfile(icpDir, 'minname.mat'), 'file') && ~overwrite
    answer = questdlg('Min matrix already exists? Overwrite?') ;
else
    answer = 'Yes' ;
end
% Note that c=1 is the reference dataset, so cycle through all others
% starting with index=2
if strcmp(answer, 'Yes')
    % Ensure the output directory for the plots
    ssrDir = fullfile(icpDir, 'ssr') ;
    if ~exist(ssrDir, 'dir')
        mkdir(ssrDir)
    end
    for cc = 2:ndatasets
        
        % define the SSR directory for this dataset
        ssrccDir = fullfile(ssrDir, sprintf('dataset_%03d', cc)) ;
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
                ssrii_fn = fullfile(ssrccDir, sprintf('ssr_tp%04d_rsub%03d.mat', ii, rsubsampling)) ;
                
                % Assume that we must redo calculation unless proven
                % that otherwise if file exists and overwrite == false
                % ---------------------------------------------------
                redoii = true ;
                if exist(ssrii_fn, 'file') && ~overwrite
                    tmp = load(ssrii_fn, 'rlist', 'ssr') ;
                    % Check that the rlist is indeed the same
                    if all(tmp.rlist == rlist)
                        % Calculation exists for the same rlist for this ii
                        redoii = false ;
                        % assign ssr for this cc, ii
                        ssr = tmp.ssr ;
                    end
                end
                
                % (Re)Compute the SSR(cc, ii, rr)
                if redoii
                    % prepare the dataset pointcloud
                    cCloud = pcread(fullfile(mca{cc, ii}.folder, mca{cc, ii}.name)); % the mesh in question
                    cxyz = cCloud.Location(1:ssample_factor:end,:) ;
                    
                    % preallocate ssr
                    ssr = zeros(length(rlist), 1) ;
                    
                    % Consider each reference mesh to find match
                    % Use dummy index qq to allow subsampled rlist
                    for qq = 1:length(rlist)
                        rr = rlist(qq) ;
                        disp(['c=',sprintf('%d',cc),' r=',sprintf('%d',rr),' i=',sprintf('%d',ii)])
                        
                        % Load the CAAX mesh used for comparison
                        refCloud = pcread(fullfile(mca{1, rr}.folder, mca{1, rr}.name));
                        refxyz = refCloud.Location;
                        refxyz = refxyz(1:ssample_factor:end,:);
                        tform = pcregistericp(cCloud,refCloud,'Extrapolate',true); % extracting the 4x4 transform from cmesh to refmesh
                        cxyzpad = horzcat(cxyz, (zeros(length(cxyz),1)+1)); % pad the locations with a ones vector
                        cxyztrans = cxyzpad*tform.T; % apply the transform
                        cxyztrans = cxyztrans(:,1:3); % remove padding
                        
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
                        ssr(qq) = sqrt(ssr1*ssr2);
                    end
                    
                    % save matching curve for this cc tp == ii
                    save(ssrii_fn, 'rlist', 'ssr')
                end
                % -------------------------------------------------
                
                ssr = movmedian(ssr, 11);
                [~,ix] = min(ssr);
                
                %find polyfit of min to find error of match
                range = 20 ;
                rlist = 1:max_ntp;
                if ix + range > length(ssr)
                    endmax = length(ssr) ;
                else
                    endmax = ix + range ;
                end
                if ix > range
                    ypoly = ssr(ix-range:endmax) ;
                    xpoly = rlist(ix-range:endmax) ;
                elseif ix < range
                    xpoly = rlist(1:2*range);
                    ypoly = ssr(1:2*range);
                elseif ix + range > length(rlist)
                    xpoly = rlist(length(rlist)-2*range:end) ;
                    ypoly = ssr(length(rlist)-2*range:end) ;
                end
                
                [f, S] = polyfix(xpoly', ypoly, 2, ix, min(ssr), ix, 0) ;
                [poly, delta] = polyval(f, round(linspace(ix-15, ix+15, 30)), S) ;
                Rinv = inv(S.R) ;
                covm = (Rinv*Rinv')*S.normr^2/S.df ;
                vara = covm(1,1) ;
                varb = covm(2,2) ;
                varc = covm(3,3) ;
                dix = (ix/2)*sqrt(vara/f(1)^2 + varb/f(2)^2) ;
                dix
                % store minimum index as minddssr
                minddssr(cc, ii) = rlist(ix);
                % Place reference timepoint that matches this (cc, ii)
                % tuple into minname
                minname{cc, ii} = mca{1, rlist(ix)}.name;
                minerror(cc, ii) = dix ; % store error to use later in plotting
                minweights(cc, ii) = 1/dix^2 ;
                
                % make figure invisible for speedup
                fig = figure('visible', 'off') ;
                plot(rlist, ssr, '.');
                hold on
                plot(round(linspace(ix-15, ix+15, 30)),poly, 'color', 'magenta', 'linewidth', 1.5) ;
                ylim([min(ssr)-10 max(ssr)+10])
                title(['Matching Dataset ' num2str(cc),' TP=' sprintf('%03d', ii)])
                xlabel('Reference Timepoint [min]', 'Interpreter', 'Latex');
                labstr = 'Mismatch, $\langle \Sigma |\mathrm{min}(\vec{x}_c - \vec{x}_{\mathrm{ref}})|^2 \rangle$ [$\mu$m]$^2$' ;
                ylabel(labstr, 'Interpreter', 'Latex');
                outfn = fullfile(ssrccDir, ['ICPPlot2_Dataset_',sprintf('%03d',cc),'_',sprintf('%03d',ii),'.png']) ;
                disp(['Saving ' outfn]);
                saveas(fig, outfn);
                close all
            end
        end
    end
    
    save(fullfile(icpDir, 'minddssr'), 'minddssr');
    save(fullfile(icpDir, 'minname'), 'minname');
    save(fullfile(icpDir, 'minweights'), 'minweights');
    save(fullfile(icpDir, 'minerror'), 'minerror');
    save(fullfile(icpDir, 'numpts1'), 'numpts1');
    save(fullfile(icpDir, 'numpts2'), 'numpts2');
end

%% Plotting "metabolism" (relative development rates) using ICP
[color, colornames] = define_colors() ;
% note: previously defined colors by hand
% color = {[217 54 104]/255, [144 115 50]/255, [76 133 50]/255,...
%     [54 131 123]/255, [59 125 171]/255, ...
%     [186 64 204]/255, [217 54 104]/255};
close all
figure('visible', 'off')
load(fullfile(icpDir, 'minddssr')) ;
load(fullfile(icpDir, 'minerror')) ;
load(fullfile(icpDir, 'minname')) ;
load(fullfile(icpDir, 'minweights')) ;
shape = define_markers(ndatasets) ;
% note: previously defined markers by hand
% shape = {'o', 'square', '^', 'V', 'd', 'p'} ;

% Plot the curves
hold on
leg = {};
rr = linspace(1, max_ntp, max_ntp);
plot(rr, rr, 'k-');
leg{1} = labels{1};
lgd = legend(leg, 'Location', 'NorthEastOutside');
% Blanket zero-valued minddssr as NaN
minddssr(minddssr<1e-5) = nan;
pss = 1 ; % plot subsampling
% Plot each dataset's correspondence times
ICPt = zeros(length(ndatasets), length(mca)) ;
for cc = 2:ndatasets
    ctime = minddssr(cc, :);
    rr = rr(~isnan(ctime));
    y = ctime(~isnan(ctime));
    yerror = minerror(cc, :);
    yerror = yerror(~isnan(ctime));
    y = movmedian(y,3) ;
    [~, fitmin] = min(abs(30-y)) ;
    [~, fitmax] = min(abs(length(rlist)-20-y)) ;
    %p = polyfit(rr(fitmin:fitmax), y(fitmin:fitmax), 1) ;
    w = minweights(cc, :) ;
    w = w(~isnan(ctime)) ;
    linfun = @(a,rr) a(1).*rr + a(2) ;
    p = nlinfit(rr(fitmin:fitmax), y(fitmin:fitmax), linfun, [0 0],...
        'Weights', w(fitmin:fitmax)) ;
    pfit = p(1)*rr+p(2);
    xshift = p(2)/p(1) ;
    rr = rr + xshift ;
    %plot(rr(1:pss:end), y(1:pss:end), shape{cc}, 'Color', color(cc, :));
    sh = shadedErrorBar(rr(1:pss:end), y(1:pss:end), yerror(1:pss:end),...
        'lineprops', '-', 'lineprops', {'color', color(cc,:)}) ;
    leg{length(leg)+1} = labels{cc} ;
    % Fit to a quadratic profile
    plot(rr, pfit, '--', 'Color', color(cc,:));
    leg{length(leg)+1} = ['Slope = ', num2str(p(1),'%10.2f')];
    lgd = legend(leg);
    rr = linspace(1, max_ntp, max_ntp);
    for tt = 1:length(rlist)
        if any(tt>pfit) && any(tt<pfit)
            [~, ix] = min(abs(tt-pfit)) ;
            ICPt(cc,tt) = ix ;
        end
    end
    xlim([0 154]) ;
    ylim([0 154]) ;
end
clearvars y r p yy ii cc

ax.XDir = 'normal' ;
ax.YDir = 'normal' ;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 20]);
title('Relative Rate of Development', 'fontsize', fs)
xlabel('Time [min]', 'fontsize', fs);
ylabel('Morphological time [min]', 'fontsize', fs);
pbaspect([1,1,1])
lgd.FontSize=fs;
% Save to disk
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figWidth figHeight]);  
saveas(gcf, fullfile(icpDir, 'rate_of_development.png')) ;
saveas(gcf, fullfile(icpDir, 'rate_of_development.pdf')) ;
clearvars lgd ax

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
        ssrccDir = fullfile(ssrDir, sprintf('dataset_%03d', cc)) ;
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
                ssrii_fn = fullfile(ssrccDir, sprintf('ssr_tp%04d_rsub%03d.mat', ii, rsubsampling)) ;
                
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

%% Plotting and extracting match info from Chi Squareds

load(fullfile(icpDir, 'ChiSqs')) ;
load(fullfile(icpDir, 'numpts1')) ;

% Cycle through ChiSqs Array
for cc = 2:size(ChiSqs, 1)
    for rr = 1:size(ChiSqs,2)
        if ~isempty(ChiSqs{cc,rr})
            clearvars xpoly ypoly
            Chis2 = ChiSqs{cc,rr} ; % extract Chis2 info
            Chis2 = movmedian(Chis2, 5); % smooth Chis2
            [~,ix] = min(Chis2); % extract min of Chis2 = best fit
            %find polyfit of min to find error of match
            range =  20;
            rlist = 1:max_ntp;
            if ix+range > size(ChiSqs,2)
                maxend = size(ChiSqs,2);
            else
                maxend = ix + range ;
            end
            if ix <= range
                xpoly = rlist(1:maxend) ;
                ypoly = Chis2(1:maxend) ;
            elseif ix > range
                xpoly = rlist(ix-range:maxend) ;
                ypoly = Chis2(ix-range:maxend) ;
            elseif ix + range > length(rlist)
                xpoly = rlist(length(rlist)-2*range:end) ;
                ypoly = Chis2(length(rlist)-2*range:end) ;
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
            ylabel('$$X^2$$', 'Interpreter', 'Latex');
            ylim([0 max(Chis2)]) ;
            xlim([0 length(rlist)]) ;
            outfn = fullfile('/mnt/data/analysis/ICP_Plots/Residuals', ['Chi2_Dataset_', num2str(cc),'_','TP_', num2str(rr), '.png']) ;
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


%% Plotting Overlays using Morphological Time Plot
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
    for rr = 1:120
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
