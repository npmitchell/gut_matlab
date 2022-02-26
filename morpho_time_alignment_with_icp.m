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
addpath('/mnt/data/code/gut_matlab/')
addpath('/mnt/data/code/gut_matlab/curve_functions/')
addpath('/mnt/data/code/imsaneV1.2.3/generalfunctions')

%% Create Cell Array with Meshe Metadata

dmap = buildLookupMap('/mnt/data/analysis/lookupMeta.txt'); % build the map that contains all data and its metadata
mca = {} ; % initiate mesh cell array

% Iterate over each marker
for mi = 1:length(keys(dmap))
    labels = keys(dmap) ;
    label = labels{mi} ;
    % Cycle through all datasets of this marker
    for j=1:length(dmap(label).folders)
        
        % get col
        if strcmp(label, 'caax')
            col = 1 ;
        elseif strcmp(label, 'hrfp')
            col = length(dmap('caax').folders) + 1 ;
        elseif strcmp(label, 'la')
            col = length(dmap('caax').folders) + length(dmap('hrfp').folders) + 1 ;
        end
        
        meshes = dir(fullfile(dmap(label).folders{j}, 'aligned_meshes',  'mesh_apical_stab_0*_APDV_um.ply')) ;
        
        % Write each mesh to the cell array mca
        for k = 1:length(meshes)
            mca{k, col+j-1} = meshes(k);
        end
    end
end

%% Use ICP to align meshes in time

outdir = '/mnt/data/analysis/' ;
cd(outdir)

if ~exist([outdir, 'ICP_Plots'], 'dir')
    mkdir /mnt/data/analysis/ICP_Plots/ ;
end

rlist = 1:size(mca,1);
minssr = zeros(size(mca,1), size(mca,2)); % initiate matrix that will store SSR data for each dataset

if exist('/mnt/data/analysis/ICP_Plots/minname.mat')
    answer = questdlg('Min matrix already exists? Overwrite?') ;
end
% Note that c=1 is the reference dataset, so cycle through all others
% starting with index=2
if strcmp(answer, 'Yes')
    for c = 2:size(mca,2)
        if c == 2
            area = dmap('caax').area{2};
        elseif c == 3
            area = dmap('hrfp').area{1};
        elseif c == 4
            area = dmap('hrfp').area{2};
        elseif c == 5
            area = dmap('hrfp').area{3};
        elseif c == 6
            area = dmap('la').area{1};
        end
        lastTP = length(area) ;
        for r = rlist
            %the CAAX mesh used for comparison
            clear ssr
            refCloud = pcread(fullfile(mca{r,1}.folder, mca{r,1}.name));
            refxyz = refCloud.Location;
            refxyz = refxyz(1:30:end,:);
            if r<=lastTP
                for i = 1:lastTP % index of dataset thats being parsed against CAAX
                    if isfield(mca{i,c}, 'name') == 1
                        
                        disp(['c=',sprintf('%d',c),' r=',sprintf('%d',r),' i=',sprintf('%d',i)])
                        cCloud = pcread(fullfile(mca{i,c}.folder, mca{i,c}.name)); % the mesh in question
                        cxyz = cCloud.Location ;
                        tform = pcregistericp(cCloud,refCloud,'Extrapolate',true); % extracting the 4x4 transform from cmesh to refmesh
                        cxyz = horzcat(cxyz, (zeros(length(cxyz),1)+1)); % pad the locations with a ones vector
                        cxyz = cxyz*tform.T; % apply the transform
                        cxyz = cxyz(1:30:end,1:3); % remove padding
                        
                        % find sum of squared residuals using two different ref orders
                        % and finding geometric mean to penalize
                        % dissimilarities in both directions of time
                        indx1 = pointMatch(refxyz, cxyz);
                        indx2 = pointMatch(cxyz, refxyz);
                        ssr1 = sum(sum((cxyz(indx1,:)-refxyz).^2,2));
                        ssr2 = sum(sum((refxyz(indx2,:)-cxyz).^2,2));
                        ssr(i) = (ssr1+ssr2)/2;
                    end
                end
                
                if length(area) > length(ssr)
                    area = area(ssr' ~= 0) ;
                elseif length(area) < length(ssr)
                    ssr = ssr(area' ~= 0);
                end
                
                ddssr = ssr'./area ;
                ddssr = ddssr(ddssr>0);
                [~,ix] = min(smooth(ddssr));
                minddssr(r,c) = ix;
                minname{r,c} = mca{ix,c}.name;
                
                %Comparing min mesh to ref mesh------------
                %             cCloud = pcread(fullfile(mca{ix,c}.folder, mca{ix,c}.name)); % the mesh in question
                %             cxyz = cCloud.Location ;
                %             refCloud = pcread(fullfile(mca{r,1}.folder, mca{r,1}.name));
                %             refxyz = refCloud.Location;
                %             refxyz = refxyz(1:30:end, :);
                %             tform = pcregistericp(cCloud,refCloud,'Extrapolate',true); % extracting the 4x4 transform from cmesh to refmesh
                %             cxyz = horzcat(cxyz, (zeros(length(cxyz),1)+1));
                %             cxyz = cxyz*tform.T; % apply the transform
                %             cxyz = cxyz(1:end,1:3);
                %             refCloud = pcread(fullfile(mca{r,1}.folder, mca{r,1}.name));
                %             refxyz = refCloud.Location;
                %             clf
                %             hold on
                %             pcshow(cxyz, 'g');
                %             pcshow(refxyz, 'r');
                %             view(2);
                %             title(['Reference r=', sprintf('%d',r),' Dataset c=',sprintf('%d',c), ' minTP=',sprintf('%d',ix),' Lateral View']);
                %             saveas(gcf,[outdir,'ICP_Plots/','MinCompare2_dataset_',sprintf('%d',c),'_minTP_',sprintf('%d',ix),'_reftp_',sprintf('%d',r),'_LateralView.png']);
                %             hold off
                %             clf
                %             hold on
                %             pcshow(cxyz, 'g');
                %             pcshow(refxyz, 'r');
                %             title(['Reference r=', sprintf('%d',r),' Dataset c=',sprintf('%d',c), ' minTP=',sprintf('%d',ix),' 3D View']);
                %             saveas(gcf,[outdir,'ICP_Plots/','MinCompare2_dataset_',sprintf('%d',c),'_minTP_',sprintf('%d',ix),'_reftp_',sprintf('%d',r),'_3DView.png']);
                %             hold off
                %------------------------------
                
                plot(ddssr');
                title(['Similarity of Dataset ',sprintf('%d',c), ' to reference TP=',sprintf('%d',r)])
                xlabel('TP of considered Dataset [min]');
                ylabel('SSR/SA');
                saveas(gcf,[outdir,'ICP_Plots/','ICPPlot2_Dataset_',sprintf('%d',c),'_RefTP_',sprintf('%d',r),'.png']);
                disp(['Saving ',outdir,'ICP_Plots/','ICPPlot2_Dataset_',sprintf('%d',c),'_RefTP_',sprintf('%d',r),'.png']);
                close
            end
        end
    end
    save('/mnt/data/analysis/ICP_Plots/minddssr', 'minddssr');
    save('/mnt/data/analysis/ICP_Plots/minname', 'minname');
end
%% Plotting "metabolism" (relative development rates)
color = {[217 54 104]/255, [144 115 50]/255, [76 133 50]/255, [54 131 123]/255, [59 125 171]/255, [186 64 204]/255, [217 54 104]/255};
close all
load /mnt/data/analysis/ICP_Plots/minddssr.mat
load /mnt/data/analysis/ICP_Plots/minname.mat
minddssr = minddssr';
minname = minname' ;
shape = {'o', 'square', '^', 'V', 'd', 'p'} ;
labels = {'CAAX 201902072000', 'CAAX 201903211930', 'RFP 201901021550', 'RFP 201904031830', 'RFP 201903312000', 'LifeAct 201904021800'};
hold on
leg = {};
r = linspace(1,154,154);
z = r;
plot(r, z, 'r-');
leg{1} = labels{1};
lgd = legend(leg);
minddssr(minddssr<5) = nan;
ICPt = {};

for c = 2:6
    y = minddssr(:,c);
    r = r(~isnan(y));
    y = y(~isnan(y));
    p = polyfit(r', y, 2);
    %fit = p(1)*r.^2 + p(2)*r + p(3);
    plot(r,y, [color{c},shape{c}]);
    leg{length(leg)+1} = labels{c};
    leg{length(leg)+1} = [labels{c}, ' Best Fit'];
    %plot(r, fit, ['--', color{c}]);
    lgd = legend(leg);
    r = linspace(1,154,154);
    fit = p(1)*r.^2 + p(2)*r + p(3);
    for i = 1:154
        [diff, ix] = min(abs(fit-i));
        ICPt{i,c} = ix;
    end
end

ax.XDir = 'normal' ;
ax.YDir = 'normal' ;
title('Development of Various Drosophila Embryo')
xlabel('Considered Timepoints [min]');
ylabel('Corresponding Timepoints [min]');
lgd.FontSize=10;
%% Plotting meshes in 3D using ICP
outdir = '/mnt/data/analysis/' ;
labels = {'CAAX 201902072000', 'CAAX 201903211930', 'RFP 201901021550', 'RFP 201904031830', 'RFP 201903312000', 'LifeAct 201904021800'};;
for ttype = 1%1:5
    if ttype == 5 %quaternion method
        clear tform
        Q = [] ; % initiate matrix of quaternions
        load '/mnt/data/analysis/ICP_Plots/tvals.mat';
        for c = 1:size(minddssr,2)
            for r = 1:size(tvals,2)
                clear x; clear y; clear z; clear w;  
                syms x y z w ;
                % linearly index each tform value into an equation to solve
                eqn1 = 1 - 2*y^2 - 2*z^2 == tvals(c,r,1,1) ;
                eqn2 = 2*x*y + 2*z*w == tvals(c,r,2,1) ;
                eqn3 = 2*x*z - 2*y*w == tvals(c,r,3,1) ;
                eqn4 = 2*x*y-2*z*w == tvals(c,r,1,2) ;
                eqn5 = 1 - 2*x^2 - 2*z^2 == tvals(c,r,2,2) ;
                eqn6 = 2*y*z + 2*x*w == tvals(c,r,3,2) ;
                eqn7 = 2*z*x + 2*y*w == tvals(c,r,1,3) ;
                eqn8 = 2*y*z - 2*x*w == tvals(c,r,2,3) ;
                eqn9 = 1 - 2*x^2 - 2*y^2 == tvals(c,r,3,3) ;
                % solve for components of the quaternion for c and r
                sols = solve(eqn1, eqn2, eqn3, eqn4) ;
                % build matrix of quaternions
                Q(1,r) = sols.w(1) ;
                Q(2,r) = sols.x(1) ;
                Q(3,r) = sols.y(1) ;
                Q(4,r) = sols.z(1) ;
            end
            Q2 = Q*Q' ; % multiply Q by its transpose 
            [V,D] = eig(Q2) ; % find all eigenvalues of Q2 and their corresponding eigenvectors
            linind = find(D == max(D,[],'all')) ; % find linear index of max eigenvalue
            [m,n] = ind2sub([size(D,1) size(D,2)], linind) ; % turn that linear index into a 2d index
            maxev = V(:,n) ; % use column of max eigenvalue to find corresponding eigenvector
            q = maxev/norm(maxev) ; % normalize the averaged quaternion
            % assign component values and rebuild rotation matrix
            tformq(1,1,c) = 1 - 2*q(3)^2 - 2*q(4)^2 ;
           	tformq(1,2,c) = 2*q(2)*q(3) + 2*q(4)*q(1) ; 
            tformq(1,3,c) = 2*q(2)*q(4) - 2*q(3)*q(4) ;
           	tformq(2,1,c) = 2*q(2)*q(3)-2*q(4)*q(1) ;
            tformq(2,2,c) = 1 - 2*q(2)^2 - 2*q(4)^2 ;
            tformq(2,3,c) = 2*q(3)*q(4) + 2*q(2)*q(1) ;
            tformq(3,1,c) = 2*q(4)*q(2) + 2*q(3)*q(1) ;
            tformq(3,2,c) = 2*q(3)*q(4) - 2*q(2)*q(1) ;
            tformq(3,3,c) = 1 - 2*q(2)^2 - 2*q(3)^2 ;
            % finally, take mean of translation vector!
            tformq(4,1,c) = mean(tvals(c,:,4,1),2) ;
            tformq(4,2,c) = mean(tvals(c,:,4,2),2) ;
            tformq(4,3,c) = mean(tvals(c,:,4,3),2) ;
            tformq(:,4,c) = 0 ;
            tformq(4,4,c) = 1 ;
        end
    end
    for r = 1:size(minname, 1)
        leg = {};
        leg{1} = labels{1};
        close all
        fig = figure('visible', 'off') ;
        
        hold on
        % First plot the reference mesh (c==1)
        refCloud = pcread(fullfile(mca{r,1}.folder, mca{r,1}.name));
        refxyz = refCloud.Location;
        refxrange = max(refxyz(:,1)) - min(refxyz(:,1)) ;
        refyrange = max(refxyz(:,2)) - min(refxyz(:,2)) ;
        refzrange = max(refxyz(:,3)) - min(refxyz(:,3)) ;
        refxoff = min(refxyz(:,1)) ;
        refxyz(:,1) = refxyz(:,1) - refxoff ;
        refMesh = read_ply_mod(fullfile(mca{r,1}.folder, mca{r,1}.name)) ;
        
        drefxy = double(refxyz(:,1:2));
        %aShape = alphaShape(drefxy,5);
        %bnd = boundaryFacets(aShape) ;
        %bnd = [bnd(:, 1); bnd(1,1)] ;
        bnd = boundary(drefxy(:,1),drefxy(:,2)) ;
        plot(drefxy(bnd, 1), drefxy(bnd,2)) ;
        
        %trisurf(refMesh.f, refMesh.v(:, 1), refMesh.v(:, 2), refMesh.v(:, 3), ...
        %   'FaceColor', color{c}, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
        
        % Now consider all other datasets (c~=1)
        for c = 2:size(minname,2)
            if minname{r,c} ~= 0
                for i = 1:size(mca,1)
                    if ~isempty(mca{i,c}) == 1
                        if strcmp(mca{i,c}.name, minname{r,c}) == 1
                            cCloud = pcread(fullfile(mca{i,c}.folder, mca{i,c}.name));
                            cxyz = cCloud.Location ;
                            
                            % CHANGE THE WAY YOU TRANSFORM GIVEN 'ttype'
                            % ------------------------------------------
                            if ttype == 1
                                clear tform
                                tform1 = pcregistericp(cCloud,refCloud,'Extrapolate',true); % extracting the 4x4 transform from cmesh to refmesh
                                for m = 1:4
                                    for n = 1:4
                                        tvals(c,r,m,n) = tform1.T(m,n);
                                    end
                                end
                                tform = tform1.T ;
                                pn = [outdir,'ICP_Plots/','Compare_TP',sprintf('%d',r),'_MinTransform.png'] ;
                                save('/mnt/data/analysis/ICP_Plots/tvals', 'tvals');
                            elseif ttype == 2
                                clear tform
                                load '/mnt/data/analysis/ICP_Plots/tvals.mat';
                                tform = squeeze(tvals(c,1,:,:)); % 4x4 transform corresponding to first TP
                                pn = [outdir,'ICP_Plots/','Compare_TP',sprintf('%d',r),'_FirstTPTransform.png'] ;
                            elseif ttype == 3
                                clear tform
                                load '/mnt/data/analysis/ICP_Plots/tvals.mat';
                                for m = 1:4
                                    for n = 1:4
                                        tform(m,n) = squeeze(mean(tvals(c,:,m,n), 2)); % 4x4 transform corresponding to (elementwise) average of all time for a certain c
                                    end
                                end
                                pn = [outdir,'ICP_Plots/','Compare_TP',sprintf('%d',r),'_TimeAvgTransform.png'] ;
                            elseif ttype == 4
                                clear tform
                                load '/mnt/data/analysis/ICP_Plots/tvals.mat';
                                TP = 60;
                                tform = squeeze(tvals(c,TP,:,:)); % 4x4 transform corresponding to TP (change TP above to alter)
                                pn = [outdir,'ICP_Plots/','Compare_TP',sprintf('%d',r),'_TP', sprintf('%d',TP), 'Transform.png'] ;
                            elseif ttype == 5
                                tform = tformq(:,:,c) ;
                                pn = [outdir,'ICP_Plots/','Compare_TP',sprintf('%d',r),'_QuaternionTransform.png'] ;
                            end
                            %--------------------------------------------------------------------------------------------------------------------
                            
                            cxyz = horzcat(cxyz, ones(length(cxyz),1)); % pad the 4th dim with ones to apply tform
                            cxyz = cxyz*tform; % apply the transform
                            cxyz = cxyz(:,1:3);
                            cxoff = min(cxyz(:,1)) ;
                            cxyz(:,1) = (cxyz(:,1)-cxoff) ;
                            % Plot the transformed mesh surface
                            cMesh = read_ply_mod(fullfile(mca{i,c}.folder, mca{i,c}.name)) ;
                            dcxy = double(cxyz(:,1:2));
                            %aShape = alphaShape(dcxy,5);
                            %bnd = boundaryFacets(aShape) ;
                            %bnd = [bnd(:, 1); bnd(1,1)] ;
                            bnd = boundary(dcxy(:,1),dcxy(:,2)) ;
                            plot(dcxy(bnd, 1), dcxy(bnd,2))
                            %trisurf(cMesh.f, cxyz(:, 1), cxyz(:, 2), cxyz(:, 3), ...
                            %  'FaceColor', color{c}, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
                            
                            leg{length(leg)+1} = labels{c};
                        end
                    end
                end
            end
        end
        axis equal
        xlim([-30 300]);
        ylim([-100 100]);
        xlabel('AP Position [$\mu$m]', 'Interpreter', 'Latex');
        ylabel('Lateral Position [$\mu$m]', 'Interpreter', 'Latex');
        lgd = legend(leg);
        lgd.FontSize=8;
        title(['Overlayed Midguts TP=', sprintf('%d',r)]);
        saveas(fig,pn);
        clf
        hold off
    end
end
%% Plotting tform values over time
load /mnt/data/analysis/ICP_Plots/tvals.mat
tvals(tvals == 0) = nan;
for c = 2%:6
    clf
    figure
    hold on
    leg = {};
    for n = 1:3
        for m = 1:3
            y = tvals(c,:,m,n);
            y = y(~isnan(y));
            %             for i = 2:length(y)
            %                 y(i) = abs(y(i)/y(1));
            %             end
            plot(y);
            leg{length(leg)+1} = ['n=',sprintf('%0d', n), ' m=', sprintf('%0d',m)];
            lgd = legend(leg);
        end
    end
    title(['tform rot entries for dataset c=', sprintf('%0d',c)])
    xlabel('Entry');
    ylabel('Value');
    lgd.FontSize=10;
    %saveas(gcf,[outdir,'ICP_Plots/','tform_Dataset_',sprintf('%d',c),'.png']);
end
for c = 2%:6
    figure
    hold on
    leg = {};
    for n = 1:3
        for m = 4
            y = tvals(c,:,m,n);
            y = y(~isnan(y));
            %             for i = 1:length(y)
            %                 y(i) = abs(y(i)/y(1));
            %             end
            plot(y);
            leg{length(leg)+1} = ['n=',sprintf('%0d', n), ' m=', sprintf('%0d',m)];
            lgd = legend(leg);
        end
    end
    title(['tform transl entries for dataset c=', sprintf('%0d',c)])
    xlabel('Entry');
    ylabel('Value');
    lgd.FontSize=10;
    %saveas(gcf,[outdir,'ICP_Plots/','tform_Dataset_',sprintf('%d',c),'.png']);
end