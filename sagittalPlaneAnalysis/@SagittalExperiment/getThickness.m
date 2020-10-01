function thickOut = getThickness(SE, overwrite)
%
%
%

if nargin < 2
    overwrite = false ;
end

% Check that all fields are present in SE.thickness
allFieldsPresent = true ;
try
    SE.thickness.thickness_eval ;
    SE.thickness.thickness_medial ;
    SE.thickness.avg_thickness ;
    SE.thickness.unc_thickness ;
    SE.thickness.avg_thickness_um ;
    SE.thickness.unc_thickness_um ;
    SE.thickness.x0 ;
    SE.thickness.y0 ;
    SE.thickness.landmarkIds ;
    SE.thickness.landmarkNames ;
catch
    allFieldsPresent = false ;
end


if allFieldsPresent && ~any( structfun(@isempty, SE.thickness) )
    if nargout > 0
        thickOut = SE.thickness ;
    end
else
    SE.getImL() ;
    SE.getBW() ;
    % Landmarks for curve extraction
    SE.getLandmarks() ;
    % Find curve for each pair of landmarks
    SE.getSkel() ;
    skel_ss = SE.skel.skel_ss ;

    % Load thickness DTs
    SE.getThicknessMasks() ;

    % Extract thickness from distance transform (DT) values at medial
    % curve points
    thicknessFn = fullfile(SE.zDir.thickness, ...
        sprintf('thickness_%04d.mat', SE.currentTime)) ;
    tFigFn = fullfile(SE.zDir.thickness, sprintf('thickness_%04d', SE.currentTime)) ;
    recompute = true ;
    if exist(thicknessFn, 'file') && ~overwrite 
        thickOut = load(thicknessFn) ;
        recompute = false ;
        try            
            % Check that all fields are present
            thickness_eval = thickOut.thickness_eval ;
            thickness_medial = thickOut.thickness_medial ;
            avg_thickness = thickOut.avg_thickness ;
            unc_thickness = thickOut.unc_thickness ;
            avg_thickness_um = thickOut.avg_thickness_um ;
            unc_thickness_um = thickOut.unc_thickness_um ;
            xo = thickOut.xo ;
            yo = thickOut.yo ;
            landmarkIds = thickOut.landmarkIds ;
            landmarkNames = thickOut.landmarkNames ;
            
            % Attribute to self
            SE.thickness = thickOut ;
        catch
            recompute = true ;
        end
    end
    
    if recompute
        % For each cycle of landmarks
        avg_thickness_um = cell(0) ;
        unc_thickness_um = cell(0) ;
        avg_thickness = cell(0) ;
        unc_thickness = cell(0) ;
        thickness_medial = cell(0) ;
        thickness_eval = cell(0) ;
        xos = cell(0) ;
        yos = cell(0) ;
        for qq = 1:length(SE.landmarks)
            skel_qq = skel_ss{qq} ;
            lmnameCell = strsplit(SE.landmarks{qq}.id, '/') ;

            % For each linesegment connecting landmark pairs
            for pp = 1:size(SE.landmarks{qq}.v, 1)

                % Grab skel_ss indices between each pair of landmarks
                pairId = [pp mod(pp+1, size(SE.landmarks{qq}.v, 1))] ;
                if pairId(2) == 0
                    pairId(2) = size(SE.landmarks{qq}.v, 1) ;
                end
                lmPair = SE.landmarks{qq}.v(pairId, :) ;
                pair = pointMatch(lmPair, skel_ss{qq}) ;
                if pair(1) == SE.npts_skel 
                    pair(1) = 1 ;
                end
                if pair(2) == 1 
                    pair(2) = SE.npts_skel ;
                end

                % Inspect pair
                % Comput DT for this leg of the landmark cycle
                disp('Computing DT') ;
                if ~isempty(SE.thicknessMasks{qq}.mask{pp})
                    outside = 1 - SE.bw .* logical(SE.thicknessMasks{qq}.mask{pp})' ;
                    medial = bwskel(logical(SE.bw .* SE.thicknessMasks{qq}.mask{pp}')) ;
                else
                    outside = 1 - SE.bw ;
                    medial = bwskel(logical(SE.bw)) ;
                end
                DT = bwdistsc(outside) ;
                [xo, yo] = find(medial) ;
                if SE.plotting.preview
                    figure('visible', 'off')
                    SE.getImL() ;
                    imshow(SE.imL); hold on
                    h1 = imagesc(DT') ; 
                    colormap(brewermap(256, '*OrRd'))
                    dtim = DT ;
                    dtim(dtim > 0) = dtim(dtim > 0) + max(dtim(:)) ;
                    set(h1, 'AlphaData', double(dtim'/max(dtim(:))));
                    plot(xo, yo, '.');
                    plot(skel_qq(pair(1):pair(2), 1), ...
                        skel_qq(pair(1):pair(2), 2), '.', 'color', [0.99, 0.99, 0.99])
                    title(['medial curve & skeleton, t= ' num2str(SE.currentTime)])
                    scatter(skel_qq(pair, 1), skel_qq(pair, 2), 100, 'filled', 'w')
                    axis equal
                    saveas(gcf, [tFigFn sprintf('_cycle%02d_seg%02d.png', qq, pp)])
                    close all
                end
                
                xos{qq}{pp} = xo ;
                yos{qq}{pp} = yo ;

                % OPTION 1 -- from medial line's DT
                % Point-match each skel pt to medial curve for thickness(s)
                idx = pointMatch(skel_qq(pair(1):pair(2), :), ...
                    [xo, yo]) ;
                imIdx = sub2ind(size(medial), xo(idx), yo(idx)) ;
                thickness_medial{qq}(pair(1):pair(2)) = DT(imIdx) ;

                % OPTION 2 -- directly from DT eval at skel_ss
                imIdx = sub2ind(size(DT), ...
                    round(skel_qq(pair(1):pair(2), 1)), ...
                    round(skel_qq(pair(1):pair(2), 2))) ;
                thickness_eval{qq}(pair(1):pair(2)) = DT(imIdx) ;
                
                % Save landmark IDs from pointmatching
                landmarkIds{qq}(pp) = pair(1) ;
                landmarkIds{qq}(pp+1) = pair(2) ;
                landmarkNames{qq}{pp} = lmnameCell{pp} ;
                landmarkNames{qq}{pp+1} = lmnameCell{pp+1} ;
            end
            assert(pair(2) == SE.npts_skel ) ;

            % Set uncertainty = difference between two measurements
            avg_thickness{qq} = 0.5 * (thickness_medial{qq} + thickness_eval{qq}) ;
            unc_thickness{qq} = 0.5 * abs(thickness_medial{qq} - thickness_eval{qq}) ;

        end

        % Save thickness measurement
        for qq = 1:length(avg_thickness)
            avg_thickness_um{qq} = avg_thickness{qq} * SE.resolution ;
            unc_thickness_um{qq} = unc_thickness{qq} * SE.resolution ;
        end
        xo = xos ;
        yo = yos ;
        save(thicknessFn, 'thickness_eval', 'thickness_medial', ...
            'avg_thickness', 'unc_thickness', ...
            'avg_thickness_um', 'unc_thickness_um', ...
            'xo', 'yo', 'landmarkIds', 'landmarkNames')

        thickOut = struct('thickness_eval', thickness_eval, ...
            'thickness_medial', thickness_medial, ...
            'avg_thickness', avg_thickness, ...
            'unc_thickness', unc_thickness, ...
            'avg_thickness_um', avg_thickness_um, ...
            'unc_thickness_um', unc_thickness_um, ...
            'xo', xos, ...
            'yo', yos, ...
            'landmarkIds', landmarkIds, ...
            'landmarkNames', landmarkNames) ;
        
        % Attribute thickness
        SE.thickness = thickOut ;
    end
end

% Save image
thickImFn = fullfile(SE.zDir.thickness_images, sprintf('thickness_%04d.png', SE.currentTime)) ;
if ~exist(thickImFn, 'file') || overwrite 

    close all
    figure('visible', 'off')
    imagesc(1:size(SE.imL, 2), 1:size(SE.imL,1), SE.imL); hold on;
    axis equal
    title(['skeleton, $t=$', num2str(SE.currentTime)], 'interpreter', 'latex')
    for qq = 1:length(skel_ss)
        scatter(skel_ss{qq}(:, 1), skel_ss{qq}(:, 2), ...
            2, avg_thickness_um{qq}, 'filled')
        caxis([2, 7])
    end
    colormap(hot(256))
    cbar = colorbar() ;
    ylabel(cbar, ['thickness [' SE.spaceUnits ']'], ...
        'interpreter', 'latex')
    axis off
    pos = get(gca, 'position') ;
    cb_pos = get(cbar, 'position') ;
    % set(cbar, 'position', [cb_pos(1) pos(2) cb_pos(3) pos(4)])
    saveas(gcf, thickImFn)
end