function kymoData = getKymographs(SE, overwriteData, overwriteImages)
% Compute kymograph intensity data for a given zplane (current zplane)
%
%

if nargin < 2
    overwriteData = false ;
    overwriteImages = false ;
elseif nargin < 3
    overwriteImages = false ;
end
originalTime = SE.currentTime ;

kDataFn = sprintf(SE.filename.kymo, SE.zplane) ;
if ~exist(SE.zDir.kymographs, 'dir')
    mkdir(SE.zDir.kymographs)
end
vkfigFn = fullfile(SE.zDir.kymographs, 'sagittal_kymograph_ventral.png') ;
dkfigFn = fullfile(SE.zDir.kymographs, 'sagittal_kymograph_dorsal.png') ;
vkfigFnSwap = fullfile(SE.zDir.kymographs, 'sagittal_kymograph_ventral_swap.png') ;
dkfigFnSwap = fullfile(SE.zDir.kymographs, 'sagittal_kymograph_dorsal_swap.png') ;

vkResampleFigFn = fullfile(SE.zDir.kymographs, ...
    'sagittal_kymograph_resampled_ventral.png') ;
dkResampleFigFn = fullfile(SE.zDir.kymographs, ...
    'sagittal_kymograph_resampled_dorsal.png') ;
vkResampleFigFnSwap = fullfile(SE.zDir.kymographs, ...
    'sagittal_kymograph_resampled_ventral_swap.png') ;
dkResampleFigFnSwap = fullfile(SE.zDir.kymographs, ...
    'sagittal_kymograph_resampled_dorsal_swap.png') ;

if ~exist(vkfigFn, 'file') || ~exist(dkfigFn, 'file') || ...
        ~exist(kDataFn, 'file') || overwriteData || overwriteImages
    if exist(kDataFn, 'file') && ~overwriteData
        kymoData = load(kDataFn, 'landmarks', 'VSS', 'VTT', 'Vkymo', ...
            'DSS', 'DTT', 'Dkymo', 'resampled', 'LR_or_DV') ;
        LR_or_DV = kymoData.LR_or_DV ;
        VSS = kymoData.VSS ;
        VTT = kymoData.VTT ;
        DSS = kymoData.DSS ;
        DTT = kymoData.DTT ;
        Vkymo = kymoData.Vkymo ;
        Dkymo = kymoData.Dkymo ;
        try
            v1 = kymoData.landmarks.v1 ;
            v2 = kymoData.landmarks.v2 ;
            v3 = kymoData.landmarks.v3 ;
            d1 = kymoData.landmarks.d1 ;
            d2 = kymoData.landmarks.d2 ;
            d3 = kymoData.landmarks.d3 ;
            idorsal = kymoData.resampled.idorsal ;
            iventral = kymoData.resampled.iventral ;
            sdorsal = kymoData.resampled.sdorsal ;
            sventral = kymoData.resampled.sventral ;
        catch
            v1 = kymoData.landmarks.r1 ;
            v2 = kymoData.landmarks.r2 ;
            v3 = kymoData.landmarks.r3 ;
            d1 = kymoData.landmarks.l1 ;
            d2 = kymoData.landmarks.l2 ;
            d3 = kymoData.landmarks.l3 ;
            idorsal = kymoData.resampled.ileft ;
            iventral = kymoData.resampled.iright ;
            sdorsal = kymoData.resampled.sleft ;
            sventral = kymoData.resampled.sright ;
        end
    else
        for tidx = 1:length(SE.timePoints)
            tp = SE.timePoints(tidx) ;
            disp(['kymos, t= ', num2str(tp)])
            SE.setZPlane(SE.zplane)
            SE.setTime(tp)

            % Load skel
            skel = SE.getSkel() ;
            skel_ss = skel.skel_ss ;
            
            % Load skelIntensity for kymograph
            SE.getSkelIntensity() ;
            s1 = SE.skelIntensity.s1 ;
            s2 = SE.skelIntensity.s2 ;
            assert(any(s2{1}))
            assert(~isempty(s2{1}))
            pathlength_um = SE.skelIntensity.pathlength_um ;
            landmarkIds = SE.skelIntensity.landmarkIds ;
            landmarkNames = SE.skelIntensity.landmarkNames ;

            % Preallocate kymograph gridded data
            if tidx == 1
                VXkymo = nan * zeros(length(SE.timePoints), 2*SE.npts_skel, length(SE.channels)) ;
                VYkymo = nan * zeros(length(SE.timePoints), 2*SE.npts_skel, length(SE.channels)) ;
                DXkymo = nan * zeros(length(SE.timePoints), 2*SE.npts_skel, length(SE.channels)) ;
                DYkymo = nan * zeros(length(SE.timePoints), 2*SE.npts_skel, length(SE.channels)) ;
                
                Vkymo = nan * zeros(length(SE.timePoints), 2*SE.npts_skel, length(SE.channels)) ;
                Dkymo = nan * zeros(length(SE.timePoints), 2*SE.npts_skel, length(SE.channels)) ;
                VSS = nan * zeros(length(SE.timePoints), 2*SE.npts_skel) ;
                VTT = nan * zeros(length(SE.timePoints), 2*SE.npts_skel) ;
                DSS = nan * zeros(length(SE.timePoints), 2*SE.npts_skel) ;
                DTT = nan * zeros(length(SE.timePoints), 2*SE.npts_skel) ;
            end

            clearvars Vcolors Dcolors
            % Build curve from anterior to posterior
            sId = 1 ;
            pathL0 = 0 ;
            for qq = 1:length(s1)
                aId = landmarkIds{qq}(find(contains(landmarkNames{1}, 'a'), 1)) ;
                pId = landmarkIds{qq}(find(contains(landmarkNames{1}, 'p'), 1)) ;
                eId = sId + pId - aId ;

                Vkymo(tidx*2-1, sId:eId, 1) = min(double(s1{qq}(aId:pId)) / SE.norms(1), 1) ;
                Vkymo(tidx*2-1, sId:eId, 2) = min(double(s2{qq}(aId:pId)) / SE.norms(2), 1) ;
                Vkymo(tidx*2, sId:eId, 1) = Vkymo(tidx*2-1, sId:eId, 1) ;
                Vkymo(tidx*2, sId:eId, 2) = Vkymo(tidx*2-1, sId:eId, 2) ;
                VSS(tidx*2-1, sId:eId) = pathL0 + pathlength_um{qq}(aId:pId) ;
                VTT(tidx*2-1, sId:eId) = (tp - 0.5) * ones(length(pathlength_um{qq}(aId:pId)),1) * SE.dt ;
                VSS(tidx*2, sId:eId) = VSS(tidx*2-1, sId:eId) ;
                VTT(tidx*2, sId:eId) = VTT(tidx*2-1, sId:eId) + SE.dt ;

                % Positional kymos
                VXkymo(tidx*2-1, sId:eId) = skel_ss{qq}(aId:pId, 1) ;
                VYkymo(tidx*2-1, sId:eId) = skel_ss{qq}(aId:pId, 2) ;
                VXkymo(tidx*2, sId:eId) = skel_ss{qq}(aId:pId, 1) ;
                VYkymo(tidx*2, sId:eId) = skel_ss{qq}(aId:pId, 2) ;
                
                % Store Offsets
                plOffV(tidx, qq) = pathL0 ;
                
                % Update offsets and indices for next skeleton segment
                sId = pId + 1;
                pathL0 = max(pathlength_um{qq}(aId:pId)) ;
                
            end

            % Build curve from posterior to anterior along dorsal side
            % Determine what pathlength value to begin the posterior point
            % This will be sum of ventral pathlengths, but the regions are
            % reversed for dorsal side.
            pathL0 = nanmax(VSS(tidx*2-1,:)) ;
            sId = 1 ;
            for qqReverse = 1:length(s1)
                % Query the lobes in reverse order: moving from P -> A
                qq = length(s1) - qqReverse + 1 ;

                pId = landmarkIds{qq}(find(contains(landmarkNames{1}, 'p'), 1)) ;
                aId = landmarkIds{qq}(find(contains(landmarkNames{1}, 'a'), 2)) ;
                aId = aId(2) ;
                eId = sId + aId - pId ;

                Dkymo(tidx*2-1, sId:eId, 1) = min(double(s1{qq}(pId:aId)) / SE.norms(1), 1) ;
                Dkymo(tidx*2-1, sId:eId, 2) = min(double(s2{qq}(pId:aId)) / SE.norms(2), 1) ;
                Dkymo(tidx*2, sId:eId, 1) = Dkymo(tidx*2-1, sId:eId, 1) ;
                Dkymo(tidx*2, sId:eId, 2) = Dkymo(tidx*2-1, sId:eId, 2) ;
                start_um = min(pathlength_um{qq}(pId:aId)) ;
                DSS(tidx*2-1, sId:eId) = pathL0 - start_um + pathlength_um{qq}(pId:aId) ;
                DTT(tidx*2-1, sId:eId) = (tp - 0.5) * ones(length(pathlength_um{qq}(pId:aId)),1) * SE.dt ;
                DSS(tidx*2, sId:eId) = DSS(tidx*2-1, sId:eId) ;
                DTT(tidx*2, sId:eId) = DTT(tidx*2-1, sId:eId) + SE.dt ;

                % Positional kymos
                DXkymo(tidx*2-1, sId:eId) = skel_ss{qq}(pId:aId, 1) ;
                DYkymo(tidx*2-1, sId:eId) = skel_ss{qq}(pId:aId, 2) ;
                DXkymo(tidx*2, sId:eId) = skel_ss{qq}(pId:aId, 1) ;
                DYkymo(tidx*2, sId:eId) = skel_ss{qq}(pId:aId, 2) ;
                
                % Store offsets
                plOffD(tidx, qq) = pathL0 - start_um ;

                % Update offsets and indices for next skeleton segment 
                pathL0 = nanmax(DSS(tidx*2-1, sId:eId)) ;
                sId = sId + length(sId:eId) ;
            end
        end

        % Mask out zeros --> note that time TT is never zero since
        % we +/-(dt * 0.5)
        VTT(VTT == 0 ) = NaN ;
        VSS(VTT == 0 ) = NaN ;
        DTT(DTT == 0 ) = NaN ;
        DSS(DTT == 0 ) = NaN ;

        % Plot landmarks
        v1 = nan * zeros(length(SE.timePoints), 1) ;
        v2 = nan * zeros(length(SE.timePoints), 1) ;
        v3 = nan * zeros(length(SE.timePoints), 1) ;
        d1 = nan * zeros(length(SE.timePoints), 1) ;
        d2 = nan * zeros(length(SE.timePoints), 1) ;
        d3 = nan * zeros(length(SE.timePoints), 1) ;
        LR_or_DV = '' ;  % decide if left-right sagittal or dorsal-ventral
        for tidx = 1:length(SE.timePoints)
            tp = SE.timePoints(tidx) ;
            disp(['kymos, t= ', num2str(tp)])
            SE.setZPlane(SE.zplane)
            SE.setTime(tp)
            SE.getSkelIntensity()
            
            for skelId = 1:length(SE.skelIntensity.landmarkNames)
                lmNames = SE.skelIntensity.landmarkNames{skelId} ;
                pum = SE.skelIntensity.pathlength_um{skelId} ;
                lmids = SE.skelIntensity.landmarkIds{skelId} ;
                for lm = 1:length(lmNames)
                    
                    if isempty(LR_or_DV)
                        if contains(lmNames{lm}, 'r') ||  contains(lmNames{lm}, 'l')
                            LR_or_DV = 'LR' ;
                        elseif contains(lmNames{lm}, 'r') ||  contains(lmNames{lm}, 'l')
                            LR_or_DV = 'DV' ;
                        end
                    end
                        
                    if contains(lmNames{lm}, 'v1') || contains(lmNames{lm}, 'r1')
                        v1(tidx) = pum(lmids(lm)) + plOffV(tidx, skelId) ;
                    elseif contains(lmNames{lm}, 'v2') || contains(lmNames{lm}, 'r2')
                        v2(tidx) = pum(lmids(lm)) + plOffV(tidx, skelId) ;
                    elseif contains(lmNames{lm}, 'v3') || contains(lmNames{lm}, 'r3')
                        v3(tidx) = pum(lmids(lm)) + plOffV(tidx, skelId) ;
                    elseif contains(lmNames{lm}, 'p2')
                        v2(tidx) = pum(lmids(lm)) + plOffV(tidx, skelId) ;
                    elseif contains(lmNames{lm}, 'p3')
                        v3(tidx) = pum(lmids(lm)) + plOffV(tidx, skelId) ;
                    end
                end

                for lm = 1:length(lmNames)
                    if contains(lmNames{lm}, 'd1') || contains(lmNames{lm}, 'l1')
                        d1(tidx) = pum(lmids(lm)) + plOffD(tidx, skelId) ;
                    elseif contains(lmNames{lm}, 'd2') || contains(lmNames{lm}, 'l2')
                        d2(tidx) = pum(lmids(lm)) + plOffD(tidx, skelId) ;
                    elseif contains(lmNames{lm}, 'd3') || contains(lmNames{lm}, 'l3')
                        d3(tidx) = pum(lmids(lm)) + plOffD(tidx, skelId);
                    elseif contains(lmNames{lm}, 'p1')
                        d1(tidx) = pum(lmids(lm)) + plOffD(tidx, skelId) ;
                    elseif contains(lmNames{lm}, 'p2')
                        d2(tidx) = pum(lmids(lm)) + plOffD(tidx, skelId);
                    elseif contains(lmNames{lm}, 'p3')
                        d3(tidx) = pum(lmids(lm)) + plOffD(tidx, skelId);
                    elseif contains(lmNames{lm}, 'a2')
                        d1(tidx) = pum(lmids(lm)) + plOffD(tidx, skelId) ;
                    elseif contains(lmNames{lm}, 'a3')
                        d2(tidx) = pum(lmids(lm)) + plOffD(tidx, skelId);
                    elseif contains(lmNames{lm}, 'a4')
                        d3(tidx) = pum(lmids(lm)) + plOffD(tidx, skelId);
                    end
                end
            end
        end
        if strcmpi(LR_or_DV, 'lr')
            landmarks = struct('r1', v1, 'r2', v2, 'r3', v3, ...
                'l1', d1, 'l2', d2, 'l3', d3) ;
        else
            landmarks = struct('v1', v1, 'v2', v2, 'v3', v3, ...
                'd1', d1, 'd2', d2, 'd3', d3) ;
        end
        
        % Resample tKymo data 
        iventral = zeros(length(SE.timePoints), 1000, length(SE.channels)) ;
        sventral = zeros(length(SE.timePoints), 1000) ;
        x_ventral = zeros(length(SE.timePoints), 1000) ;
        y_ventral = zeros(length(SE.timePoints), 1000) ;
        idorsal = zeros(length(SE.timePoints), 1000, length(SE.channels)) ;
        sdorsal = zeros(length(SE.timePoints), 1000) ;
        x_dorsal = zeros(length(SE.timePoints), 1000) ;
        y_dorsal = zeros(length(SE.timePoints), 1000) ;
        for tidx = 1:length(SE.timePoints)
            tmod = tidx*2-1 ;
            
            % Resample ventral side
            [Sventral, inds] = unique(VSS(tmod, :)) ;
            Iventral = squeeze(Vkymo(tmod, inds, :)) ;
            Xventral = VXkymo(tmod, inds) ;
            Yventral = VYkymo(tmod, inds) ;
            keep = find(~isnan(Sventral)) ;
            sV = linspace(0, nanmax(Sventral, [], 2), 1000) ;
            sventral(tidx,:) = sV ;
            iventral(tidx,:,1) = interp1(Sventral(keep), Iventral(keep, 1), sV, 'linear') ;
            iventral(tidx,:,2) = interp1(Sventral(keep), Iventral(keep, 2), sV, 'linear') ;
            % Resample ventral positions
            x_ventral(tidx,:) = interp1(Sventral(keep), Xventral(keep), sV, 'linear') ;
            y_ventral(tidx,:) = interp1(Sventral(keep), Yventral(keep), sV, 'linear') ;
            

            % Resample dorsal side
            keep = find(~isnan(DSS(tmod, :)) & DSS(tmod, :) > 0) ;
            Sdorsal = DSS(tmod, keep) ;
            Idorsal = squeeze(Dkymo(tmod, keep, :)) ;
            Xdorsal = DXkymo(tmod, keep) ;
            Ydorsal = DYkymo(tmod, keep) ;
            [Sdorsal, inds] = unique(Sdorsal) ;
            Idorsal = Idorsal(inds, :) ;
            Xdorsal = Xdorsal(inds) ;
            Ydorsal = Ydorsal(inds) ;
            sD = linspace(nanmin(DSS(tmod, keep), [], 2), ...
                nanmax(DSS(tmod, keep), [], 2), 1000) ;
            sdorsal(tidx,:) = sD ;
            idorsal(tidx,:,1) = interp1(Sdorsal, Idorsal(:, 1), sD, 'linear') ;
            idorsal(tidx,:,2) = interp1(Sdorsal, Idorsal(:, 2), sD, 'linear') ;
            % Resample ventral positions
            x_dorsal(tidx,:) = interp1(Sdorsal, Xdorsal, sD, 'linear') ;
            y_dorsal(tidx,:) = interp1(Sdorsal, Ydorsal, sD, 'linear') ;

        end
        
        % Save tKymo data
        if strcmpi(LR_or_DV, 'lr')
            resampled = struct('sleft', sdorsal, 'ileft', idorsal, ...
                'sright', sventral, 'iright', iventral, ...
                'x_right', x_ventral, 'y_right', y_ventral, ...
                'x_left', x_dorsal, 'y_left', y_dorsal) ;
            save(kDataFn, 'landmarks', 'VSS', 'VTT', 'Vkymo', ...
                'DSS', 'DTT', 'Dkymo', 'resampled', 'LR_or_DV') ;
        else
            resampled = struct('sdorsal', sdorsal, 'idorsal', idorsal, ...
                'sventral', sventral, 'iventral', iventral, ...
                'x_ventral', x_ventral, 'y_ventral', y_ventral, ...
                'x_dorsal', x_dorsal, 'y_dorsal', y_dorsal) ;
            save(kDataFn, 'landmarks', 'VSS', 'VTT', 'Vkymo', ...
                'DSS', 'DTT', 'Dkymo', 'resampled', 'LR_or_DV') ;
        end
    end
    
    if ~exist(vkfigFn, 'file') || ~exist(dkfigFn, 'file') ||...
            ~exist(vkResampleFigFn, 'file') || ...
            ~exist(dkResampleFigFn, 'file') ||...
            overwriteImages
        for swapColor = 1:2
            for ii = 1:2
                if ii == 1
                    %% KYMOGRAPHS 
                    kymos = {Vkymo, Dkymo} ;
                    if swapColor == 1
                        kfns = {vkfigFn, dkfigFn} ;
                    else
                        kfns = {vkfigFnSwap, dkfigFnSwap} ;
                    end
                    ps = max(DSS, [], 2) ;
                    ps2 = ps(1:2:end) ;
                    XXs = {VSS, ps - DSS} ;
                    YYs = {VTT, DTT} ;
                    lms = {{v1, v2, v3}, {ps2 - d1, ps2 - d2, ps2 - d3}} ;
                elseif ii == 2
                    %% Resampled kymographs
                    kymos = {iventral, idorsal} ;
                    if swapColor == 1
                        kfns = {vkResampleFigFn, dkResampleFigFn} ;
                    else
                        kfns = {vkResampleFigFnSwap, dkResampleFigFnSwap} ;
                    end
                    ps = max(sdorsal, [], 2) ;  % posterior pathlength (s) coordinates
                    XXs = {sventral, ps - sdorsal} ;
                    YYs = {SE.dt * (SE.timePoints .* ones(1000, length(SE.timePoints)))', ...
                        SE.dt * (SE.timePoints .* ones(1000, length(SE.timePoints)))'} ;
                    lms = {{v1, v2, v3}, {ps - d1, ps - d2, ps - d3}} ;
                end

                % lms = {{v1, v2, v3}, {d1, d2, d3}} ;
                if strcmpi(LR_or_DV, 'dv')
                    ktitleAdd = {'ventral side', 'dorsal side'} ;
                else
                    ktitleAdd = {'right side', 'left side'} ;
                end
                for kId = 1:2
                    kymo = kymos{kId} ;

                    clf
                    % Background
                    % set(gca,'Color','k'); hold on;

                    colors = {[1,1,0], [0,0,0] } ;    
                    % Blue layer: 
                    p1 = surf(XXs{kId},YYs{kId}, zeros(size(XXs{kId})),'AlphaData', kymo(:, :, 1),...
                        'FaceAlpha','interp',...
                        'FaceColor', colors{mod(swapColor, 2)+1},...
                        'edgecolor','none');
                    % Red layer:  
                    hold on
                    p2 = surf(XXs{kId},YYs{kId},zeros(size(XXs{kId})),'AlphaData', kymo(:, :, 2),...
                        'FaceAlpha','interp',...
                        'FaceColor', colors{mod(swapColor+1, 2)+1},...
                        'edgecolor','none');
                    view(2)
                    grid off
                    set(gca,'Ydir','reverse')
                    ylabel(['time, [' SE.timeUnits ']'], 'interpreter', 'latex')
                    xlabel(['pathlength, $s$ [' SE.spaceUnits ']'], 'interpreter', 'latex')
                    title(['Sagittal kymograph: ' ktitleAdd{kId}])
                    % Note: https://www.reddit.com/r/matlab/comments/469qp1/overlaying_two_pcolor_plots/

                    % Plot landmarks
                    lm = lms{kId} ;
                    for lid = 1:length(lm) 
                        plot(lm{lid}, SE.timePoints * SE.dt, '-', ...
                            'color', SE.plotting.colors(1, :))
                    end

                    ylim([min(SE.timePoints-0.5)*SE.dt, max(SE.timePoints+0.5)*SE.dt])
                    saveas(gcf, kfns{kId})  
                end
            end
        end
    end
end

% Reset time
if ~isempty(originalTime)
    SE.setTime(originalTime) 
end

if nargout > 0
    kymoData = load(sprintf(SE.filename.kymo, SE.zplane), 'landmarks', 'VSS', 'VTT', 'Vkymo', ...
        'DSS', 'DTT', 'Dkymo', 'resampled') ;
end
