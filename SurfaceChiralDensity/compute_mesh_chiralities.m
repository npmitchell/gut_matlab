%% SURFACE CHIRAL DENSITY =========================================
% An example detailing the functionality of the discrete surface chiral
% density pseudotensor.
% By Noah P Mitchell 02/2019
% with input from Dillon Cislo
%==========================================================================

clear; close all; clc;
addpath('../mesh_handling/')
addpath('../plotting/')
coolwarm = coolwarm();

% Prepare path for meshes ===========================================
path = '../../data/48Ygal4-UAShistRFP/2016021015120_objFiles/';
meshes = dir([path, 'pointCloud_T*_ascii.ply']) ;
% path = '../../data/48Ygal4-UAShisRFP/20170329_objFiles_test_chirality_and_curvature/';
% meshes = dir([path, '*test_mesh_ascii.ply']) ;
% path = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/data/48Ygal4-UAShisRFP/2016021015120_objFiles/';
% meshes = dir([path, 'cleaned/cleaned_pointCloud*_mesh.ply']) ;

outdir = fullfile(path, 'surface_chirality/') ;
imdir = fullfile(path, 'surface_chirality_images/') ;
clims = [-1 1] ;
check = false ;
dx = 0.2619 * 2 ;
ii = 1;
xlims = [0, 200] ;
ylims = [-50, 50] ;
zlims = [-50, 50] ;

% Create output dirs
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
if ~exist(imdir, 'dir')
    mkdir(imdir)
end

% Prepare for iteration
% chisums is the integral chirality
dmyk = 1
dt = 1;
todo = 1:dt:min(length(meshes), 555) ;
chisums = zeros(length(todo), 1) ;


%% Compute chirality for each mesh
for ii = todo
    try
        % Load the mesh
        meshfn = fullfile(meshes(ii).folder, meshes(ii).name) ;
        disp(['analyzing ', meshfn])
        [tri, pts] = ply_read(meshfn, 'tri') ;
        pts = transpose(pts) * dx;
        pts(:, 3) = pts(:, 3) - mean(pts(:, 3));
        % Some samples should be flipped
        % pts(:, 3) = -pts(:, 3) + mean(pts(:, 3));
        tri = transpose(tri) ;

        % View Result --------------------------------------------------------
        if check
            trisurf(tri, pts(:,1), pts(:,2), pts(:,3));
            axis equal
        end

        % Check if any faces have normals on both sides (ie are their own front
        % and back)
        trisort = sort(tri, 2) ;
        trisort = sortrows(trisort) ;
        [trisu, ia, ib] = unique(trisort, 'rows', 'stable') ;
        keep_face = false(size(trisort, 1), 1);
        keep_face(ia) = true;
        faces_discard = find(~keep_face) ;
        % Get the rows that were duplicates
        offending = trisort(~keep_face, :);
        tri2rm = [];
        for jj=1:size(offending, 1)
            newrm = find(all(trisort == offending(jj, :), 2)) ;
            tri2rm = [tri2rm; newrm(:)] ;
        end
        tri2rm = unique(tri2rm) ;
        trikeep = transpose(setdiff(1:size(trisort, 1), tri2rm)) ;
        face = trisort(trikeep, :) ;
        vertex = pts;

        % face = trisort(keep, :) ;
        % 
        % trisurf(face, pts(:,1), pts(:,2), pts(:,3));
        % 
        % test = 0 * trisort;
        % for i =1:length(trisort)-1
        %     test(i, :) = trisort(i + 1,:) - trisort(i,:) ;
        % end
        % subtr0 = find(all(test == 0, 2)) + 1 ;

        % % Remove the edges that are shared only by these faces
        % surfTri = triangulation(face, pts);
        % % Vertex IDs defining each edge
        % eIDx = surfTri.edges;
        % edgeFace_cell = surfTri.edgeAttachments( eIDx );
        % % If any edge is shared only by these faces, remove that edge
        % for jj=1:length(edgeFace_cell)
        %     if isempty(setdiff(edgeFace_cell{jj}, faces_discard))
        %         killface = faces_discard(jj) ;
        %         find(all(edgeFace_cell == killface, 2))
        %     end
        % end

        % % Remove points associated with these faces in 'vertex'
        % toremove = trisort(discard, :) ;
        % toremove = unique(toremove(:), 'sorted') ;
        % keep = true(size(pts, 1), 1) ;
        % keep(toremove) = false ;
        % vertex = pts(keep, :) ;

        % % Remove points associated with these faces in 'face'
        % for jj=1:length(toremove)
        %     rmvpt = toremove(jj) ;
        %     face(face > rmvpt) = face(face > rmvpt) - 1 ;
        % end

        % trisurf(tri, pts(:,1), pts(:,2), pts(:,3));
        % hold on;
        % trisurf(face, vertex(:,1), vertex(:,2), vertex(:,3), -vertex(:, 3));

        %% CONSTRUCT CHIRAL DENSITY PSEUDOTENSOR ==============================
        chiral = calculate_chiral_density(face, vertex);

        %% VIEW RESULTS =======================================================

        % Create surface vector field based on faces --------------------------
        fvf = vertex(face(:, 2), :) - vertex(face(:, 1), :);
        e1 = vertex(face(:,3),:) - vertex(face(:,2),:);
        e2 = vertex(face(:,1),:) - vertex(face(:,3),:);
        fvf = fvf ./ sqrt( sum( fvf.^2, 2 ) );

        % View surface vector field -----------------------------------------------
        % CoM = cat( 3, vertex(face(:,1),:), ...
        %     vertex(face(:,2),:), vertex(face(:,3),:) );
        % CoM = mean( CoM, 3 ); % Face centroids
        % 
        % trisurf(tri, pts(:,1), pts(:,2), pts(:,3),'FaceColor',[0.8 0.8 1.0]);
        % axis equal
        % hold on
        % quiver3( CoM(:,1), CoM(:,2), CoM(:,3), ...
        %     fvf(:,1), fvf(:,2), fvf(:,3), ...
        %     0.5, 'Color', 'b' );
        % 
        % clear CoM

        %%
        % Act the chirality operators on each face vector ---------------------
        % chi = zeros( size(face,1), 1 );
        chisum = zeros(3, 3) ;
        for i = 1:size(face,1)
            % chi(i) = fvf(i,:) * chiral{i} * fvf(i,:)';
            % chi(i) = chiral{i}(1) ;
            c = chiral{i} ;
            % chi(i) = 0.5 * (trace(c)^2 - trace(c*c)) ;
            chisum = 
            
            % [eigvects, eigvals] = eig(c);
            % or equivalently, 
            % chi(i) = c(1,1)*c(2,2)+ c(2,2)*c(3,3)+c(1,1)*c(3,3) - c(2,1)*c(1,2)- c(3,2)*c(2,3)-c(1,3)*c(3,1);
        end
        chisums(dmyk) = nansum(chi) / length(vertex);

        % View the chirality --------------------------------------------------
        figh = figure('visible','off');
        colormap(coolwarm)
        patch('Faces', face, 'Vertices', vertex, ...
            'FaceVertexCData', chi, 'FaceColor', 'flat', ...
            'edgecolor', [0.5 0.5 0.5], 'linewidth', 0.1);

        axis equal
        caxis(clims)
        cbar = colorbar ;
        cbar.Label.String = 'chirality, \chi' ;
        xlabel('AP position [\mum]')
        ylabel('DV position [\mum]')
        xlim(xlims) ;
        zlim(zlims) ;   
        view(0,0) 
        % Save the figure
        namestr = split(meshes(ii).name, '_mesh') ;
        namestr = split(namestr(1), 'pointCloud_') ;
        namestr = namestr(2) ;
        name = join(['surfchirality_', namestr, '.png'], '');
        filename = fullfile(imdir, name) ;
        saveas(figh, filename{1})
        close(figh)
        
        disp(['Analyzed ii = ', num2str(ii)])
        dmyk = dmyk + 1 ;
    catch
        disp('could not process this time point')
        chisums(dmyk) = NaN;
        dmyk = dmyk + 1 ;
    end
end
disp('done analyzing PLYs')

%% Plot the results
figh = figure('visible','off');
plot(dt * (0:length(todo)-1), chisums, '.-')
ylabel('Chirality, ||\chi||')
xlabel('time [min]')
title('Surface chirality over time')
filename = fullfile(imdir, 'summed_chirality.png') ;
saveas(figh, filename)

