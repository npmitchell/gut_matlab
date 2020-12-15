function [length_lobes, area_lobes, volume_lobes] = ...
    aux_plot_lobe_dynamics_3d(folds, ssfold, ssmax, lobeDir, timePoints, ...
    spcutMeshBase, nV, nU, rot, trans, resolution, flipy, xyzlim, colors, save_ims, overwrite_lobeims)
%AUX_PLOT_LOBE_DYNAMICS_3D auxiliary function for Generate_Axisymmetric_Pullbacks_Orbifold.m
%   Plot the lobe dynamics for all timepoints without computing
%   area/volume/etc
%   todo: generalize to #lobes = other than 4
% 
% Parameters
% ----------
% folds : 
% ssfold :
% ssmax : 
% lobeDir :  str
%   the path to where lobe information/data/images are stored
% timePoints : #timepoints x 1 float array
%   the time in minutes of each frame
% spcutMeshBase : str
%   full file path to the (s,phi) coord sys cylinder cut mesh
% nV : int
%   sampling number along DV axis
% nU : int 
%   sampling number along AP axis
% rot : 3x3 float array
%   rotation matrix transforming data into APDV coord sys
% trans : 1x3 float array
%   translation to perform after rotation to move APDV coord sys to origin
% resolution : float 
%   um / pixel conversion for data 
% colors : (>=#lobes)x3 float array
%   colors to use to color each lobe
% save_ims : bool
%   save an image of the mesh in 3d, with each lobe colored differently
% 
%
%
% NPMitchell 2019

disp('Plotting lobes in 3d...')
% Now compute Surface area and Volume in one pass.
lobeImDir = fullfile(lobeDir, 'images_lobes3D') ;
if ~exist(lobeImDir, 'dir')
    mkdir(lobeImDir)
end

disp('Computing surface area and volume of each lobe...')
for tidx = 1:length(timePoints)
    % Translate to which timestamp
    tp = timePoints(tidx) ;
    timestr = sprintf('%04d', tp) ;
    load(sprintf(spcutMeshBase, tp), 'spcutMesh') ;

    % Load the centerline too
    % fn = sprintf(clineDVhoopBase, t) ;
    % load(fn, 'avgpts')
    avgpts = spcutMesh.avgpts ;

    % rename the fold indices (in U)
    f1 = folds(tidx, 1) ;
    f2 = folds(tidx, 2) ;
    f3 = folds(tidx, 3) ;

    % Create figure for plotting lobes
    lobeimfn = fullfile(lobeImDir, ['lobes_' timestr '.png']) ;
    redo_lobeims = save_ims && (~exist(lobeimfn, 'file') || overwrite_lobeims) ;
    
    
    % Save plot of lobes, each a different color
    if redo_lobeims
        close all
        fig = figure('visible', 'off') ;
        
        % Announce which timestamp we consider
        if mod(tp, 40) == 0
            disp([' ... saving each figure to ' lobeimfn])
        end

        % Define rotated, translated mesh
        vrs = ((rot * spcutMesh.v')' + trans) * resolution ;
        % use QS to allow flipy
        if flipy
            vrs(:, 2) = -vrs(:, 2) ;
        end

        for lobe = 1:4
            % Note that the vertices are ordered in AP strips, with
            % nU elements for each value of nV. 

            % This is for transposed ordering
            % if lobe == 1
            %     rmvtx = ((f1+1)*nV):size(vrs, 1) ;
            %     rear = ((f1-1) * nV + 1):(((f1 + 1) * nV) - 1) ;
            %     front = 1:nV ;
            % elseif lobe == 2
            %     rmvtx = [1:(f1*nV), ((f2 + 1) * nV):size(vrs, 1) ] ;
            %     rear = ((f2-1) * nV + 1):(((f2 + 1) * nV) - 1) ;
            %     front = ((f1-1) * nV + 1):(((f1 + 1) * nV) - 1) ;
            % elseif lobe == 3
            %     rmvtx = [1:(f2*nV), ((f3 + 1) * nV):size(vrs, 1) ] ;
            %     rear = ((f3-1) * nV + 1):(((f3 + 1) * nV) - 1) ;
            %     front = ((f2-1) * nV + 1):(((f2 + 1) * nV) - 1) ;
            % elseif lobe == 4
            %     rmvtx = 1:(f3*nV) ;
            %     rear = ((f4-1) * nV + 1):(((f4 + 1) * nV) - 1) ;
            %     front = ((f3-1) * nV + 1):(((f3 + 1) * nV) - 1) ;
            % end

            % Create matrix of indices, each row identical (y location)
            allvs = (0:(nV-1)) * nU ;
            if lobe == 1
                urm = (f1+1):nU ;
                rear_id = f1 ;
                front_id = 1 ;
            elseif lobe == 2
                urm = [1:(f1-1), (f2+1):nU] ;
                rear_id = f2 ;
                front_id = f1 ;
            elseif lobe == 3
                urm = [1:(f2-1), (f3+1):nU] ;
                rear_id = f3 ;
                front_id = f2 ;
            elseif lobe == 4
                urm = 1:(f3-1) ;
                rear_id = nU ;
                front_id = f3 ;
            end
            rear_id = double(rear_id) ;
            front_id = double(front_id) ;
            urm = double(urm) ;

            % Find all indices on the rear hoop and front hoop
            rear = (rear_id - front_id + 1) + ...
                (0:(nV-1)) * (rear_id - front_id + 1) ;
            front = 1 + (0:(nV-1)) * (rear_id - front_id + 1) ;
            % Create brick of indices, each column identical (x location)
            rmvtx = (urm' .* ones(length(urm), nV))' ;
            % Add y location and unravel
            rmvtx = rmvtx + allvs' .* ones(nV, length(urm)) ;
            rmvtx = rmvtx(:) ;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute area of mesh up to first fold
            % extract faces of first lobe
            % remove the faces that do not include 1:fold1
            [ newF, newV, ~ ] = remove_vertex_from_mesh(spcutMesh.f, ...
                vrs, rmvtx) ;

            % Close the lobe volume with additional triangles
            % first close the rear
            addpt = avgpts(rear_id, :) ;
            newV = cat(1, newV, addpt) ;
            addptidx = length(newV);
            backf = zeros(length(rear), 3) ;
            for qq = 1:length(rear)
                nid = mod(qq + 1, length(rear)) ;
                if nid == 0
                    nid = length(rear) ;
                end
                backf(qq, :) = [rear(qq), rear(nid), addptidx];
            end
            newF = cat(1, newF, backf) ;

            % Close the front face
            addpt = avgpts(front_id, :) ;
            newV = cat(1, newV, addpt) ;
            addptidx = length(newV);
            frontf = zeros(length(front), 3) ;
            for qq = 1:length(front)
                nid = mod(qq + 1, length(front)) ;
                if nid == 0
                    nid = length(front) ;
                end
                frontf(qq, :) = [front(qq), addptidx, front(nid)];
            end
            newF = cat(1, newF, frontf) ;
            
            % Plot lobes in 3D
            trisurf(newF, newV(:, 1), newV(:, 2), newV(:, 3), ...
                'EdgeColor', colors(lobe, :) * 0.5, 'FaceColor', colors(lobe, :));
            hold on;
        end

        axis equal
        title(['Lobes, t = ' sprintf('%03d', tp)])
        xlabel('x [\mum]')
        ylabel('y [\mum]')
        zlabel('z [\mum]')
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        disp(['Saving figure: ' lobeimfn])
        saveas(fig, lobeimfn)
        close all
    end
end
