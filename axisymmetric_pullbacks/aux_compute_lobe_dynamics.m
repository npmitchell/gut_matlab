function aux_compute_lobe_dynamics(ssfold, ssmax, lobeDir, timePoints, ...
    sphiBase, rot, trans, xyzlim, colors)
%AUX_COMPUTE_LOBE_DYNAMICS auxiliary function for Generate_Axisymmetric_Pullbacks_Orbifold.m
%   Compute the lobe dynamics for all timepoints
% 
% NPMitchell 2019

disp('Computing length_lobes...')
% Length is given by ssfold.
ss_lobes = cat(2, ssfold, ssmax) ;
length_lobes = ss_lobes ;
for qq=2:4
    length_lobes(:, qq) = ss_lobes(:, qq) - ss_lobes(:, qq-1) ;
end

disp('Computing area_lobes and volume_lobes...')
% Now compute Surface area and Volume in one pass.
lobeImDir = fullfile(lobeDir, 'images_lobes3D') ;
if ~exist(lobeImDir, 'dir')
    mkdir(lobeImDir)
end

% Preallocate arrays
area_lobes = zeros(length(timePoints), 4) ;
volume_lobes = zeros(length(timePoints), 4) ;

disp('Computing surface area and volume of each lobe...')
for kk = 1:length(timePoints)
    % Translate to which timestamp
    t = timePoints(kk) ;
    timestr = sprintf('%04d', t) ;
    load(sprintf(sphiBase, t), 'spcutMesh') ;

    % Load the centerline too
    % fn = sprintf(clineDVhoopBase, t) ;
    % load(fn, 'avgpts')
    avgpts = spcutMesh.avgpts ;

    % rename the fold indices (in U)
    f1 = folds(kk, 1) ;
    f2 = folds(kk, 2) ;
    f3 = folds(kk, 3) ;

    % preallocate
    vol_kk = [0, 0, 0, 0] ;
    area_kk = [0, 0, 0, 0] ;

    % Create figure for plotting lobes
    lobeimfn = fullfile(lobeImDir, ['lobes_' timestr '.png']) ;
    redo_lobeims = save_ims && (~exist(lobeimfn, 'file') || overwrite_lobeims) ;
    if redo_lobeims
        close all
        fig = figure('visible', 'off') ;
    end

    % Announce which timestamp we consider
    if mod(t, 40) == 0
        disp(['SA and Vol: t = ', num2str(t)])
        if redo_lobeims
            disp([' ... and saving each figure to ' lobeimfn])
        end
    end

    % Define rotated, translated mesh
    vrs = ((rot * spcutMesh.v')' + trans) * resolution ;
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

        % Find all indices on the rear hoop and front hoop
        rear = (rear_id - front_id + 1) + (0:(nV-1)) * (rear_id - front_id + 1) ;
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
        [ newF, newV, oldVertexIDx ] = remove_vertex_from_mesh(spcutMesh.f, ...
            vrs, rmvtx) ;
        % Compute area without closing the volume
        area_kk(lobe) = meshArea(newV, newF) ;

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
                nid = length(rear) ;
            end
            frontf(qq, :) = [front(qq), addptidx, front(nid)];
        end
        newF = cat(1, newF, frontf) ;

        % Finally, compute volume
        [vol_kk(lobe), ~] = meshVolumeArea(newV, newF) ;

        % Plot lobes in 3D
        if redo_lobeims
            trisurf(newF, newV(:, 1), newV(:, 2), newV(:, 3), ...
                'EdgeColor', colors(lobe, :) * 0.5, 'FaceColor', colors(lobe, :));
            hold on;
        end
    end

    % Add area_kk and vol_kk to lists
    area_lobes(kk, :) = area_kk ;
    volume_lobes(kk, :) = abs(vol_kk) ;   

    % Save plot of lobes, each a different color
    if redo_lobeims
        axis equal
        title(['Lobes, t = ' sprintf('%03d', t)])
        xlabel('x [\mum]')
        ylabel('y [\mum]')
        zlabel('z [\mum]')
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        saveas(fig, lobeimfn)
        close all
    end
end
