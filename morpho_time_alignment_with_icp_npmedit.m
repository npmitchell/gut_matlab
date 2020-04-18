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
ssample_factor = 30 ;   % subsampling of point clouds for ICP 
rsubsampling = 1 ;     % subsampling of timepoints for reference dataset
csubsampling = 1 ;     % subsampling of timepoints for matched dataset
overwrite = false ;      % overwrite previous results

%% Create Cell Array with Meshe Metadata

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

%% Use ICP to align meshes in time

rlist = 1:rsubsampling:max_ntp;
clist = 1:csubsampling:max_ntp;
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
                    cxyz = cCloud.Location ;

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
                        refxyz = refxyz(1:30:end,:);

                        tform = pcregistericp(cCloud,refCloud,'Extrapolate',true); % extracting the 4x4 transform from cmesh to refmesh
                        cxyzpad = horzcat(cxyz, (zeros(length(cxyz),1)+1)); % pad the locations with a ones vector
                        cxyztrans = cxyzpad*tform.T; % apply the transform
                        cxyztrans = cxyztrans(1:ssample_factor:end,1:3); % remove padding

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
                        ssr(qq) = (sqrt(ssr1) + sqrt(ssr2)) / 2;

                    end
                    
                    % save matching curve for this cc tp == ii
                    save(ssrii_fn, 'rlist', 'ssr')
                end

                [~,ix] = min(smooth(ssr));
                % store minimum index as minddssr
                minddssr(cc, ii) = rlist(ix);
                % Place reference timepoint that matches this (cc, ii) 
                % tuple into minname
                minname{cc, ii} = mca{1, rlist(ix)}.name;

                %Comparing min mesh to ref mesh------------
                %             cCloud = pcread(fullfile(mca{c, ix}.folder, mca{c, ix}.name)); % the mesh in question
                %             cxyz = cCloud.Location ;
                %             refCloud = pcread(fullfile(mca{1, r}.folder, mca{1, r}.name));
                %             refxyz = refCloud.Location;
                %             refxyz = refxyz(1:30:end, :);
                %             tform = pcregistericp(cCloud,refCloud,'Extrapolate',true); % extracting the 4x4 transform from cmesh to refmesh
                %             cxyz = horzcat(cxyz, (zeros(length(cxyz),1)+1));
                %             cxyz = cxyz*tform.T; % apply the transform
                %             cxyz = cxyz(1:end,1:3);
                %             refCloud = pcread(fullfile(mca{1, r}.folder, mca{1, r}.name));
                %             refxyz = refCloud.Location;
                %             clf
                %             hold on
                %             pcshow(cxyz, 'g');
                %             pcshow(refxyz, 'r');
                %             view(2);
                %             title(['Reference r=', sprintf('%d',r),' Dataset c=',sprintf('%d',c), ' minTP=',sprintf('%d',ix),' Lateral View']);
                %             saveas(gcf,[icpDir,'MinCompare2_dataset_',sprintf('%d',c),'_minTP_',sprintf('%d',ix),'_reftp_',sprintf('%d',r),'_LateralView.png']);
                %             hold off
                %             clf
                %             hold on
                %             pcshow(cxyz, 'g');
                %             pcshow(refxyz, 'r');
                %             title(['Reference r=', sprintf('%d',r),' Dataset c=',sprintf('%d',c), ' minTP=',sprintf('%d',ix),' 3D View']);
                %             saveas(gcf,[icpDir,'MinCompare2_dataset_',sprintf('%d',c),'_minTP_',sprintf('%d',ix),'_reftp_',sprintf('%d',r),'_3DView.png']);
                %             hold off
                %------------------------------

                % make figure invisible for speedup
                fig = figure('visible', 'off') ;
                plot(rlist, ssr', '.');
                hold on
                plot(rlist, smooth(ssr), '--') ;
                title(['Matching reference time to dataset ' num2str(cc) ', tp=' sprintf('%03d', ii)], 'Interpreter', 'Latex')
                xlabel('reference timepoint [min]', 'Interpreter', 'Latex');
                labstr = 'mismatch, $\langle \Sigma |\mathrm{min}(\vec{x}_c - \vec{x}_{\mathrm{ref}})| \rangle$ [$\mu$m]' ;
                ylabel(labstr, 'Interpreter', 'Latex');
                outfn = fullfile(ssrccDir, ['ICPPlot2_Dataset_',sprintf('%03d',cc),'_',sprintf('%03d',ii),'.png']) ;
                disp(['Saving ' outfn]);
                saveas(fig, outfn);
                close all
            end
        end
    end
    
    % Save global results for current csubsampling and rsubsampling
    disp('Saving minddssr and minname for all datasets')
    save(fullfile(icpDir, 'minddssr'), 'minddssr');
    save(fullfile(icpDir, 'minname'), 'minname');
end

%% Plotting "metabolism" (relative development rates)
[color, colornames] = define_colors() ;
% note: previously defined colors by hand
% color = {[217 54 104]/255, [144 115 50]/255, [76 133 50]/255,...
%     [54 131 123]/255, [59 125 171]/255, ...
%     [186 64 204]/255, [217 54 104]/255};
close all
load(fullfile(icpDir, 'minddssr'))
load(fullfile(icpDir, 'minname'))
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
ICPt = {};

% Plot each dataset's correspondence times
for cc = 2:ndatasets
    ctime = minddssr(cc, :);
    rr = rr(~isnan(ctime));
    y = ctime(~isnan(ctime));

    plot(rr, y, shape{cc}, 'Color', color(cc, :));
    leg{length(leg)+1} = labels{cc};
    
    % Fit to a quadratic profile
    p = polyfit(rr, y, 2);
    % fit = p(1)*r.^2 + p(2)*r + p(3);
    % plot(r, fit, ['--', color{cc}]);
    % leg{length(leg)+1} = [labels{cc}, ' Best Fit'];  
    
    lgd = legend(leg);
    rr = linspace(1, max_ntp, max_ntp);
    fit = p(1)*rr.^2 + p(2)*rr + p(3);
    for ii =1:154
        [diff, ix] = min(abs(fit - ii));
        ICPt{ii, cc} = ix;
    end
end
clearvars y r p yy ii cc 

ax.XDir = 'normal' ;
ax.YDir = 'normal' ;
title('Rate of development')
xlabel('Considered Mesh Timepoint [min]');
ylabel('Reference Mesh Timepoint [min]');
pbaspect([1,1,1])
lgd.FontSize=10;
saveas(gcf, fullfile(icpDir, 'rate_of_development.png')) ;
clearvars lgd ax

%% Plotting meshes in 3D using ICP
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

        
for ttype = 1:5
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
    for rr = rlist
        disp(['rr = ' num2str(rr)])
        leg = {};
        leg{1} = labels{1};
        close all
        fig = figure('visible', 'off') ;
        
        hold on
        % First plot the reference mesh (c==1)
        refCloud = pcread(fullfile(mca{1, rr}.folder, mca{1, rr}.name));
        refxyz = refCloud.Location;
        refxrange = max(refxyz(:,1)) - min(refxyz(:,1)) ;
        refyrange = max(refxyz(:,2)) - min(refxyz(:,2)) ;
        refzrange = max(refxyz(:,3)) - min(refxyz(:,3)) ;
        refxoff = min(refxyz(:,1)) ;
        refxyz(:,1) = refxyz(:,1) - refxoff ;
        refMesh = read_ply_mod(fullfile(mca{1, rr}.folder, mca{1, rr}.name)) ;
        
        drefxy = double(refxyz(:,1:2));
        %aShape = alphaShape(drefxy,5);
        %bnd = boundaryFacets(aShape) ;
        %bnd = [bnd(:, 1); bnd(1,1)] ;
        
        % Note: the third argument is the tightness parameter (S=1 is 100%
        % tight, while S=0 is convex hull)
        bnd = boundary(drefxy(:,1),drefxy(:,2), 0.95) ;
        plot(drefxy(bnd, 1), drefxy(bnd,2)) ;
        
        %trisurf(refMesh.f, refMesh.v(:, 1), refMesh.v(:, 2), refMesh.v(:, 3), ...
        %   'FaceColor', color{c}, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
        
        % Now consider all other datasets (c~=1) and find matches to this
        % reference time rr
        for cc = 2:ndatasets
            % Is there a match for this dataset at this time in c?
            
            % Find the timepoint in dataset c that matches this t_ref
            for ii =1:111
                if ~isempty(mca{cc, ii}) 
                    % Note: we search for a match to mca{rr, 1}.name
                    % in minname{cc, ii}
                    if strcmp(mca{1, rr}.name, minname{cc, ii}) 
                        disp(['Found match: d=' num2str(cc) ', ii=' num2str(ii)])
                        cCloud = pcread(fullfile(mca{cc, ii}.folder, mca{cc, ii}.name));
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
                        cMesh = read_ply_mod(fullfile(mca{cc, ii}.folder, mca{cc, ii}.name)) ;
                        dcxy = double(cxyz(:,1:2));
                        %aShape = alphaShape(dcxy,5);
                        %bnd = boundaryFacets(aShape) ;
                        %bnd = [bnd(:, 1); bnd(1,1)] ;
                        bnd = boundary(dcxy(:,1),dcxy(:,2)) ;
                        plot(dcxy(bnd, 1), dcxy(bnd,2))
                        %trisurf(cMesh.f, cxyz(:, 1), cxyz(:, 2), cxyz(:, 3), ...
                        %  'FaceColor', color{c}, 'FaceAlpha', 0.2, 'EdgeColor', 'none')

                        leg{length(leg)+1} = labels{cc};
                    end
                end

            end
        end
        xlim([-30, 280])
        ylim([-80, 80])
        xlabel('AP position [$\mu$m]', 'Interpreter', 'latex')
        ylabel('lateral position [$\mu$m]', 'Interpreter', 'latex')
        axis equal
        lgd = legend(leg, 'Location', 'NorthEastOutside');
        lgd.FontSize=8;
        title(['$t_\textrm{ref}=$', sprintf('%d',rr), ' min'], 'Interpreter', 'Latex');
        disp(['Saving overlay figure to ', pn])
        saveas(fig,pn);
        clf
        hold off
    end
end

%% Plotting tform values over time
load(fullfile(icpDir, 'tvals.mat'))
tvals(tvals == 0) = nan;
for cc = 2%:6
    clf
    figure
    hold on
    leg = {};
    for n = 1:3
        for m = 1:3
            y = tvals(cc,:,m,n);
            y = y(~isnan(y));
            %             for ii =2:length(y)
            %                 y(i) = abs(y(i)/y(1));
            %             end
            plot(y);
            leg{length(leg)+1} = ['n=',sprintf('%0d', n), ' m=', sprintf('%0d',m)];
            lgd = legend(leg);
        end
    end
    title(['tform rot entries for dataset c=', sprintf('%0d',cc)])
    xlabel('Entry');
    ylabel('Value');
    lgd.FontSize = 10;
    %saveas(gcf,[icpDir,'tform_Dataset_',sprintf('%d',c),'.png']);
end
for cc = 2%:6
    figure
    hold on
    leg = {};
    for n = 1:3
        for m = 4
            y = tvals(cc,:,m,n);
            y = y(~isnan(y));
            %             for ii =1:length(y)
            %                 y(i) = abs(y(i)/y(1));
            %             end
            plot(y);
            leg{length(leg)+1} = ['n=',sprintf('%0d', n), ' m=', sprintf('%0d',m)];
            lgd = legend(leg);
        end
    end
    title(['tform transl entries for dataset c=', sprintf('%0d',cc)])
    xlabel('Entry');
    ylabel('Value');
    lgd.FontSize=10;
    %saveas(gcf,[icpDir,'tform_Dataset_',sprintf('%d',c),'.png']);
end