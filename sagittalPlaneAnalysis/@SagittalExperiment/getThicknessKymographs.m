function tKymoData = getThicknessKymographs(SE, overwriteData, overwriteImages)
% Compute kymograph thickness data for a given zplane (current zplane)
%
%

if nargin < 2
    overwriteData = false ;
    overwriteImages = false ;
elseif nargin < 3
    overwriteImages = false ;
end
originalTime = SE.currentTime ;

kThicknessFn = sprintf(SE.filename.tKymo, SE.zplane) ;
if ~exist(SE.zDir.thickness_kymographs, 'dir')
    mkdir(SE.zDir.thickness_kymographs)
end
vtfigFn = fullfile(SE.zDir.thickness_kymographs, ...
    'thickness_kymograph_ventral.png') ;
dtfigFn = fullfile(SE.zDir.thickness_kymographs, ...
    'thickness_kymograph_dorsal.png') ;
vtfigFn_zoom = fullfile(SE.zDir.thickness_kymographs, ...
    'thickness_kymograph_ventral_zoom.png') ;
dtfigFn_zoom = fullfile(SE.zDir.thickness_kymographs, ...
    'thickness_kymograph_dorsal_zoom.png') ;
vufigFn = fullfile(SE.zDir.thickness_kymographs, ...
    'thickness_kymograph_ventral_unc.png') ;
dufigFn = fullfile(SE.zDir.thickness_kymographs, ...
    'thickness_kymograph_dorsal_unc.png') ;
vufigFn_zoom = fullfile(SE.zDir.thickness_kymographs, ...
    'thickness_kymograph_ventral_unc_zoom.png') ;
dufigFn_zoom = fullfile(SE.zDir.thickness_kymographs, ...
    'thickness_kymograph_dorsal_unc_zoom.png') ;


vtResampleFigFn = fullfile(sprintf(SE.dir.thickness_kymographs, SE.zplane), ...
    'thickness_kymograph_resample_ventral.png') ;
dtResampleFigFn = fullfile(sprintf(SE.dir.thickness_kymographs, SE.zplane), ...
    'thickness_kymograph_resample_dorsal.png') ;
vtResampleFigFn_zoom = fullfile(sprintf(SE.dir.thickness_kymographs, SE.zplane), ...
    'thickness_kymograph_resample_ventral_zoom.png') ;
dtResampleFigFn_zoom = fullfile(sprintf(SE.dir.thickness_kymographs, SE.zplane), ...
    'thickness_kymograph_resample_dorsal_zoom.png') ;
vuResampleFigFn = fullfile(sprintf(SE.dir.thickness_kymographs, SE.zplane), ...
    'thickness_kymograph_resample_ventral_unc.png') ;
duResampleFigFn = fullfile(sprintf(SE.dir.thickness_kymographs, SE.zplane), ...
    'thickness_kymograph_resample_dorsal_unc.png') ;
vuResampleFigFn_zoom = fullfile(sprintf(SE.dir.thickness_kymographs, SE.zplane), ...
    'thickness_kymograph_resample_ventral_unc_zoom.png') ;
duResampleFigFn_zoom = fullfile(sprintf(SE.dir.thickness_kymographs, SE.zplane), ...
    'thickness_kymograph_resample_dorsal_unc_zoom.png') ;


TdtVFigFn = fullfile(sprintf(SE.dir.thickness_kymographs, SE.zplane), ...
    'thicknessDynamics_kymograph_resample_ventral.png') ;
TdtDFigFn = fullfile(sprintf(SE.dir.thickness_kymographs, SE.zplane), ...
    'thicknessDynamics_kymograph_resample_dorsal.png') ;

if ~exist(vtfigFn, 'file') || ~exist(dtfigFn, 'file') || ...
        ~exist(kThicknessFn, 'file') || overwriteData || overwriteImages
    if exist(kThicknessFn, 'file') && ~overwriteData
        tKymoData = load(kThicknessFn, 'landmarks', ...
            'VSS', 'VTT', 'VtKymo', 'VuKymo', ...
            'DSS', 'DTT', 'DtKymo', 'DuKymo', 'resampled') ;
        resampled = tKymoData.resampled ;
        tventral = resampled.tventral ;
        tdorsal = resampled.tdorsal ;
        unc_ventral = resampled.unc_ventral ;
        unc_dorsal = resampled.unc_dorsal ;
        sventral = resampled.sventral ;
        sdorsal = resampled.sdorsal ;
        v1 = tKymoData.landmarks.v1 ;
        v2 = tKymoData.landmarks.v2 ;
        v3 = tKymoData.landmarks.v3 ;
        d1 = tKymoData.landmarks.d1 ;
        d2 = tKymoData.landmarks.d2 ;
        d3 = tKymoData.landmarks.d3 ;
        VSS = tKymoData.VSS ;
        VTT = tKymoData.VTT ;
        DSS = tKymoData.DSS ;
        DTT = tKymoData.DTT ;
        VtKymo = tKymoData.VtKymo ;
        DtKymo = tKymoData.DtKymo ;
        VuKymo = tKymoData.VuKymo ;
        DuKymo = tKymoData.DuKymo ;
    else
        for tidx = 1:length(SE.timePoints)
            tp = SE.timePoints(tidx) ;
            disp(['tKymos, t= ', num2str(tp)])
            SE.setZPlane(SE.zplane)
            SE.setTime(tp)

            % Load skelIntensity for tKymograph
            SE.getSkel() ;
            pathlength_um = SE.skel.pathlength ;
            for qq = 1:length(pathlength_um)
                pathlength_um{qq} = pathlength_um{qq} * SE.resolution ;
            end
            SE.getThickness() ;
            th = SE.thickness.avg_thickness_um ;
            unc = SE.thickness.unc_thickness_um ;
            landmarkIds = SE.thickness.landmarkIds ;
            landmarkNames = SE.thickness.landmarkNames ;

            % Preallocate tKymograph gridded data
            if tidx == 1
                VtKymo = nan * zeros(length(SE.timePoints), 2*SE.npts_skel) ;
                DtKymo = nan * zeros(length(SE.timePoints), 2*SE.npts_skel) ;
                VuKymo = nan * zeros(length(SE.timePoints), 2*SE.npts_skel) ;
                DuKymo = nan * zeros(length(SE.timePoints), 2*SE.npts_skel) ;
                VSS = nan * zeros(length(SE.timePoints), 2*SE.npts_skel) ;
                VTT = nan * zeros(length(SE.timePoints), 2*SE.npts_skel) ;
                DSS = nan * zeros(length(SE.timePoints), 2*SE.npts_skel) ;
                DTT = nan * zeros(length(SE.timePoints), 2*SE.npts_skel) ;
            end

            clearvars Vcolors Dcolors
            % Build curve from anterior to posterior
            sId = 1 ;
            pathL0 = 0 ;
            for qq = 1:length(th)
                aId = landmarkIds{qq}(find(contains(landmarkNames{1}, 'a'), 1)) ;
                pId = landmarkIds{qq}(find(contains(landmarkNames{1}, 'p'), 1)) ;
                eId = sId + pId - aId ;

                VtKymo(tidx*2-1, sId:eId) = double(th{qq}(aId:pId)) ;
                VtKymo(tidx*2, sId:eId) = VtKymo(tidx*2-1, sId:eId) ;
                VuKymo(tidx*2-1, sId:eId) = double(unc{qq}(aId:pId)) ;
                VuKymo(tidx*2, sId:eId) = VuKymo(tidx*2-1, sId:eId) ;
                VSS(tidx*2-1, sId:eId) = pathL0 + pathlength_um{qq}(aId:pId) ;
                VTT(tidx*2-1, sId:eId) = (tp - 0.5) * ones(length(pathlength_um{qq}(aId:pId)),1) * SE.dt ;
                VSS(tidx*2, sId:eId) = VSS(tidx*2-1, sId:eId) ;
                VTT(tidx*2, sId:eId) = VTT(tidx*2-1, sId:eId) + SE.dt ;

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
            for qqReverse = 1:length(th)
                % Query the lobes in reverse order: moving from P -> A
                qq = length(th) - qqReverse + 1 ;

                pId = landmarkIds{qq}(find(contains(landmarkNames{1}, 'p'), 1)) ;
                aId = landmarkIds{qq}(find(contains(landmarkNames{1}, 'a'), 2)) ;
                aId = aId(2) ;
                eId = sId + aId - pId ;

                DtKymo(tidx*2-1, sId:eId) = double(th{qq}(pId:aId))  ;
                DtKymo(tidx*2, sId:eId) = DtKymo(tidx*2-1, sId:eId) ;
                DuKymo(tidx*2-1, sId:eId) = double(unc{qq}(pId:aId))  ;
                DuKymo(tidx*2, sId:eId) = DuKymo(tidx*2-1, sId:eId) ;
                start_um = min(pathlength_um{qq}(pId:aId)) ;
                DSS(tidx*2-1, sId:eId) = pathL0 - start_um + pathlength_um{qq}(pId:aId) ;
                DTT(tidx*2-1, sId:eId) = (tp - 0.5) * ones(length(pathlength_um{qq}(pId:aId)),1) * SE.dt ;
                DSS(tidx*2, sId:eId) = DSS(tidx*2-1, sId:eId) ;
                DTT(tidx*2, sId:eId) = DTT(tidx*2-1, sId:eId) + SE.dt ;

                % Store offsets
                plOffD(tidx, qq) = pathL0 - start_um ;

                % Update offsets and indices for next skeleton segment 
                pathL0 = nanmax(DSS(tidx*2-1, sId:eId)) ;
                sId = sId + length(sId:eId) ;
            end
        end

        % Mask out zeros
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
        for tidx = 1:length(SE.timePoints)
            tp = SE.timePoints(tidx) ;
            disp(['tKymos, t= ', num2str(tp)])
            SE.setZPlane(SE.zplane)
            SE.setTime(tp)
            SE.getSkelIntensity()

            for skelId = 1:length(SE.skelIntensity.landmarkNames)
                lmNames = SE.skelIntensity.landmarkNames{skelId} ;
                pum = SE.skelIntensity.pathlength_um{skelId} ;
                lmids = SE.skelIntensity.landmarkIds{skelId} ;
                for lm = 1:length(lmNames)
                    if contains(lmNames{lm}, 'v1') 
                        v1(tidx) = pum(lmids(lm)) + plOffV(tidx, skelId) ;
                    elseif contains(lmNames{lm}, 'v2')
                        v2(tidx) = pum(lmids(lm)) + plOffV(tidx, skelId) ;
                    elseif contains(lmNames{lm}, 'v3')
                        v3(tidx) = pum(lmids(lm)) + plOffV(tidx, skelId) ;
                    elseif contains(lmNames{lm}, 'p2')
                        v2(tidx) = pum(lmids(lm)) + plOffV(tidx, skelId) ;
                    elseif contains(lmNames{lm}, 'p3')
                        v3(tidx) = pum(lmids(lm)) + plOffV(tidx, skelId) ;
                    end
                end

                for lm = 1:length(lmNames)
                    if contains(lmNames{lm}, 'd1')
                        d1(tidx) = pum(lmids(lm)) + plOffD(tidx, skelId) ;
                    elseif contains(lmNames{lm}, 'd2')
                        d2(tidx) = pum(lmids(lm)) + plOffD(tidx, skelId) ;
                    elseif contains(lmNames{lm}, 'd3')
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
        landmarks = struct('v1', v1, 'v2', v2, 'v3', v3, ...
            'd1', d1, 'd2', d2, 'd3', d3) ;

        % Resample tKymo data 
        tventral = zeros(length(SE.timePoints), 1000) ;
        unc_ventral = zeros(length(SE.timePoints), 1000) ;
        sventral = zeros(length(SE.timePoints), 1000) ;
        tdorsal = zeros(length(SE.timePoints), 1000) ;
        unc_dorsal = zeros(length(SE.timePoints), 1000) ;
        sdorsal = zeros(length(SE.timePoints), 1000) ;
        for tidx = 1:length(SE.timePoints)
            tmod = tidx*2-1 ;
            
            % Resample ventral side
            [Xventral, inds] = unique(VSS(tmod, :)) ;
            Tventral = VtKymo(tmod, inds) ;
            Uventral = VuKymo(tmod, inds) ;
            keep = find(~isnan(Xventral)) ;
            sV = linspace(0, nanmax(Xventral, [], 2), 1000) ;
            sventral(tidx,:) = sV ;
            tventral(tidx,:) = interp1(Xventral(keep), Tventral(keep), sV, 'linear') ;
            unc_ventral(tidx,:) = interp1(Xventral(keep), Uventral(keep), sV, 'linear') ;

            % Resample dorsal side
            keep = find(~isnan(DSS(tmod, :)) & DSS(tmod, :) > 0) ;
            Xdorsal = DSS(tmod, keep) ;
            Tdorsal = DtKymo(tmod, keep) ;
            Udorsal = DuKymo(tmod, keep) ;
            [Xdorsal, inds] = unique(Xdorsal) ;
            Tdorsal = Tdorsal(inds) ;
            Udorsal = Udorsal(inds) ;
            sD = linspace(nanmin(DSS(tmod, keep), [], 2), ...
                nanmax(DSS(tmod, keep), [], 2), 1000) ;
            sdorsal(tidx,:) = sD ;
            tdorsal(tidx,:) = interp1(Xdorsal, Tdorsal, sD, 'linear') ;
            unc_dorsal(tidx,:) = interp1(Xdorsal, Udorsal, sD, 'linear') ;

            assert(~any(isnan(tventral(tidx, :))))
            assert(~any(isnan(tdorsal(tidx, :))))
        end
        
        % Save tKymo data
        resampled = struct('sdorsal', sdorsal, 'tdorsal', tdorsal, ...
            'sventral', sventral, 'tventral', tventral, ...
            'unc_ventral', unc_ventral, 'unc_dorsal', unc_dorsal) ;
        save(kThicknessFn, 'landmarks', 'VSS', 'VTT', ...
            'VtKymo', 'VuKymo', ...
            'DSS', 'DTT', 'DtKymo', 'DuKymo', 'resampled') ;
    end
    
    for ii = 1:2
        if ii == 1
            %% KYMOGRAPHS 
            tKymos = {VtKymo, DtKymo} ;
            uKymos = {VuKymo, DuKymo} ;
            tfns = {vtfigFn, dtfigFn} ;
            tfns_zoom = {vtfigFn_zoom, dtfigFn_zoom} ;
            ufns = {vufigFn, dufigFn} ;
            ufns_zoom = {vufigFn_zoom, dufigFn_zoom} ;
            ps = max(DSS, [], 2) ;
            ps2 = ps(1:2:end) ;
            XXs = {VSS, ps - DSS} ;
            YYs = {VTT, DTT} ;
            lms = {{v1, v2, v3}, {ps2 - d1, ps2 - d2, ps2 - d3}} ;
        elseif ii == 2
            %% Resampled kymographs
            tKymos = {tventral, tdorsal} ;
            uKymos = {unc_ventral, unc_dorsal} ;
            tfns = {vtResampleFigFn, dtResampleFigFn} ;
            tfns_zoom = {vtResampleFigFn_zoom, dtResampleFigFn_zoom} ;
            ufns = {vuResampleFigFn, duResampleFigFn} ;
            ufns_zoom = {vuResampleFigFn_zoom, duResampleFigFn_zoom} ;
            ps = max(sdorsal, [], 2) ;  % posterior pathlength (s) coordinates
            XXs = {sventral, ps - sdorsal} ;
            YYs = {SE.dt * (SE.timePoints .* ones(1000, length(SE.timePoints)))', ...
                SE.dt * (SE.timePoints .* ones(1000, length(SE.timePoints)))'} ;
            lms = {{v1, v2, v3}, {ps - d1, ps - d2, ps - d3}} ;
        end
        
        % Now plot it
        ktitleAdd = {'ventral side', 'dorsal side'} ;
        for kId = 1:2
            tKymo = tKymos{kId} ;
            uKymo = uKymos{kId} ;
            clf

            % thickness: 
            p1 = surf(XXs{kId}, YYs{kId}, tKymo,... % 'AlphaData', tKymo(:, :, 1),...
                'FaceColor','interp',...
                'edgecolor','none');
            cb = colorbar() ;
            colormap(parula)
            ylabel(cb, ['thickness [' SE.spaceUnits ']'], 'Interpreter', 'Latex') ;
            hold on ;
            view(2)
            grid off
            set(gca,'Ydir','reverse')
            ylabel(['time, [' SE.timeUnits ']'], 'interpreter', 'latex')
            xlabel(['pathlength, $s$ [' SE.spaceUnits ']'], 'interpreter', 'latex')
            title(['Thickness: ' ktitleAdd{kId}])
            % Note: https://www.reddit.com/r/matlab/comments/469qp1/overlaying_two_pcolor_plots/

            % Plot landmarks
            lm = lms{kId} ;
            for lid = 1:length(lm) 
                plot3(lm{lid}, SE.timePoints * SE.dt, ...
                    (max(tKymo(:))+1) * ones(size(SE.timePoints)), '-', ...
                    'color', SE.plotting.colors(2, :)); 
                hold on;
            end
            ylim([min(SE.timePoints-0.5)*SE.dt, max(SE.timePoints+0.5)*SE.dt])
            saveas(gcf, tfns{kId})  
            caxis([0, 6])
            saveas(gcf, tfns_zoom{kId})  

            % Plot uncertainty
            clf        
            p1 = surf(XXs{kId},YYs{kId}, uKymo,... % 'AlphaData', tKymo(:, :, 1),...
                'FaceColor','interp',...
                'edgecolor','none');
            cb = colorbar() ;
            colormap(parula)
            ylabel(cb, ['uncertainty in thickness [' SE.spaceUnits ']'], 'Interpreter', 'Latex') ;
            hold on ;
            view(2)
            grid off
            set(gca,'Ydir','reverse')
            ylabel(['time, [' SE.timeUnits ']'], 'interpreter', 'latex')
            xlabel(['pathlength, $s$ [' SE.spaceUnits ']'], 'interpreter', 'latex')
            title(['Thickness uncertainty: ' ktitleAdd{kId}])

            % Plot landmarks
            lm = lms{kId} ;
            for lid = 1:length(lm) 
                plot3(lm{lid}, SE.timePoints * SE.dt, ...
                    (max(tKymo(:))+1) * ones(size(SE.timePoints)), '-', ...
                    'color', SE.plotting.colors(2, :)); 
                hold on;
            end
            ylim([min(SE.timePoints-0.5)*SE.dt, max(SE.timePoints+0.5)*SE.dt])
            saveas(gcf, ufns{kId})
            caxis([0, 1])
            saveas(gcf, ufns_zoom{kId})  
        end
    end
    
    % Plot gradients in thickness
    close all
    thick = {tventral, tdorsal} ;
    ps = max(sdorsal, [], 2) ;  % posterior pathlength (s) coordinates
    ss = {sventral, ps - sdorsal} ;
    fns = {TdtVFigFn, TdtDFigFn} ;
    for kId = 1:2
        Tblur = imgaussfilt(thick{kId}, 2);
        [~,gradT] = gradient(Tblur, 1, SE.dt) ;
        Tdt = imgaussfilt(gradT) / SE.dt ;
        % convolve with square pulse function
        A = 1:10 ;
        A = [A fliplr(A)] ;
        A = A ./ sum(A(:)) ;
        Tdt2 = conv2(Tdt, A, 'same');
        
        %% Save as square array
        % imagesc(Tdt2)
        % caxis([-max(abs(Tdt(:))), max(abs(Tdt(:)))])
        % bwr256 = bluewhitered(256) ;
        % colormap(bwr256)
        % cb = colorbar() ;
        % title('Change in thickness over time')
        % ylabel(cb, ['$\partial_t h$ [' SE.spaceUnits '/' SE.timeUnits ']'], ...
        %     'Interpreter', 'latex')
        % clf

        timeM = SE.dt * (SE.timePoints .* ones(1000, length(SE.timePoints)))' ;
        p1 = surf(ss{kId},timeM, Tdt2,... % 'AlphaData', tKymo(:, :, 1),...
            'FaceColor','interp',...
            'edgecolor','none');
        hold on;
        caxis([-max(abs(Tdt(:))), max(abs(Tdt(:)))])
        bwr256 = bluewhitered(256) ;
        colormap(bwr256)
        cb = colorbar() ;
        title('Change in thickness over time')
        ylabel(cb, ['$\partial_t h$ [' SE.spaceUnits '/' SE.timeUnits ']'], ...
            'Interpreter', 'latex')
        view(2)
        grid off
        set(gca,'Ydir','reverse')
        ylabel(['time, [' SE.timeUnits ']'], 'interpreter', 'latex')
        xlabel(['pathlength, $s$ [' SE.spaceUnits ']'], 'interpreter', 'latex')
        title(['Thickness dynamics: ' ktitleAdd{kId}])

        % Plot landmarks
        lm = lms{kId} ;
        for lid = 1:length(lm) 
            plot3(lm{lid}, SE.timePoints * SE.dt, ...
                (max(tKymo(:))+1) * ones(size(SE.timePoints)), '-', ...
                'color', SE.plotting.colors(2, :)); 
            hold on;
        end
        ylim([min(SE.timePoints)*SE.dt, max(SE.timePoints)*SE.dt])
        
        saveas(gcf, fns{kId})
        close all
    end
    close all
end



% Reset time
if ~isempty(originalTime)
    SE.setTime(originalTime) 
end
