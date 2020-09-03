function featureIDs = measurePathlineFeatureIDs(QS, pathlineType, options)
%measurePathlineFeatureIDs(QS, options)
% Interactively identify one longitudinal (zeta) position per feature for 
% Lagrangian pathlines using field measurements along Lagrangian pathlines. 
%
% One can use, for example, radius of the pathlines on the mesh from the 
% centerline, or the normal velocity, or dz -- the ratio of projected 
% (embedding space) proper distance of a unit vector along longitudinal 
% mesh direction to the pullback distance along the longitudinal direction,
% or dp -- similar for circumferential direction, divv -- the divergence of
% the velocity field, etc. Uses two fields to identify the feature
% positions along the longitudinal dimension.
%
% todo: generalize beyond vP pathline array
% 
% Parameters
% ----------
% QS : QuapSlap class object instance
% pathlineType : str specifier ('vertices', 'piv', 'faces')
%   whether pathlines in question thread through mesh vertices, piv
%   evaluation coordinates, or face barycenters at t=t0Pathlines
% options : optional struct with fields
%   overwrite : bool
%   field1 : str specifier ('radius', 'dz', 'dp', 'divv', 'veln')
%   field2 : str specifier ('radius', 'dz', 'dp', 'divv', 'veln')
% 
% Returns
% -------
% featureIDs : #features x 1 int array
%   longitudinal pullback coordinate for features, as an index into the
%   pullback pathlines (which are grid-like).
% 
% NPMitchell 2020

% Default options
overwrite = false ;
field1 = 'radius' ;
field2 = 'veln' ;
if nargin < 2
    pathlineType = 'vertices' ;
end
if isempty(QS.pathlines.t0)
    QS.pathlines.t0 = QS.t0set() ;
end
t0Pathline = QS.pathlines.t0 ;

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end

%% Identify pathline coodinates using desired scalar fields for inspection
if strcmpi(pathlineType, 'vertices')
    %% Obtain mat filenames for loading fields
    loadDir = QS.dir.metricKinematics.pathline.measurements ;
    apKymoMetricKinFn = fullfile(loadDir, 'apKymographMetricKinematics.mat') ;

    %% To grab fold location in Lagrangian coords robustly, find minima of 
    % radii from ap average as kymograph and grab folds and lobes indexed in 
    % Lagrangian coords
    % fIDfn is the feature ID filename for these Lagrangian data
    fIDfn = sprintf(QS.fileName.pathlines.featureIDs, t0Pathline) ;
    if exist(fIDfn, 'file') || ~overwrite
        load(fIDfn, 'featureIDs')
    else
        %% Interactively choose the feature locations

        % Load scalar field data to choose featureIDs (ex, folds)
        if strcmpi(field1, 'radius') || strcmpi(field2, 'radius')
            try
                tmp = load(apKymoRadiusFn, 'radius_apM') ;
                if strcmpi(field1, 'radius') 
                    sf1 = tmp.radius_apM ;
                else
                    sf2 = tmp.radius_apM ;
                end    
            catch
                error('Run QS.plotPathlineMetricKinematics() before QS.plotPathlineStrainRate()')
            end
        elseif strcmpi(field1, 'divv') || strcmpi(field2, 'divv')
            try
                tmp = load(apKymoMetricKinFn, 'divv_apM') ;
                if strcmpi(field1, 'divv') 
                    sf1 = tmp.divv_apM ;
                else
                    sf2 = tmp.divv_apM ;
                end
            catch
                error('Run QS.plotPathlineMetricKinematics() before QS.plotPathlineStrainRate()')
            end
        elseif strcmpi(field1, 'veln') || strcmpi(field2, 'veln')
            try
                tmp = load(apKymoMetricKinFn, 'veln_apM') ;
                if strcmpi(field1, 'veln') 
                    sf1 = tmp.veln_apM ;
                else
                    sf2 = tmp.veln_apM ;
                end
            catch
                error('Run QS.plotPathlineMetricKinematics() before QS.plotPathlineStrainRate()')
            end
        end

        %% Interactively adjust feature locations
        nfeatureIDs = input('How many featureIDs (ex folds) to identify in Lagrangian data? [Default=3]') ;
        if isempty(nfeatureIDs)
            nfeatureIDs = 3 ;
        end

        % Make a guess as to the features using minima of field1
        div1d = mean(sf1(tps > max(20, min(tps)) & ...
            tps < min(max(tps), 60), :), 1) ;
        div1dsm = savgol(div1d, 2, 11) ;
        [~, featureIDs] = maxk(-islocalmin(div1dsm) .* div1dsm, nfeatureIDs) ;
        featureIDs = sort(featureIDs) ;

        % Show guess overlaying field1 data
        figure ;
        set(gcf, 'visible', 'on')
        subplot(1, 2, 1)
        imagesc((1:nU)/nU, tps, sf1)
        colormap(bwr256)
        caxis([-climit, climit])
        hold on;
        for qq = 1:nfeatureIDs
            plot(featureIDs(qq)/nU * ones(size(tps)), tps) ;
        end
        title('Guess for featureIDs on divergence(v)')
        subplot(1, 2, 2) ;
        imagesc((1:nU)/nU, tps, sf2)
        colormap(bwr256)
        caxis([-max(abs(sf2(:))), max(abs(sf2(:)))])
        hold on;
        for qq = 1:nfeatureIDs
            plot(featureIDs(qq)/nU * ones(size(tps)), tps) ;
        end
        title('Guess for featureIDs on normal velocity')
        disp(['Guessed automatic featureIDs to be: [' num2str(featureIDs) ']'])

        % Update the guess
        for qq = 1:nfeatureIDs
            qok = false ;
            while ~qok
                msg = 'What is the Lagrangian zeta ID of feature ' ;
                msg = [msg num2str(qq) '? '] ;
                msg = [msg '[Default=' num2str(featureIDs(qq)) ']'] ;
                newvalley = input(msg) ;
                if isa(newvalley, 'double')
                    if ~isempty(newvalley)
                        featureIDs(qq) = newvalley ;
                    end
                end

                % Show guess overlaying div(v) data
                clf;
                set(gcf, 'visible', 'on')
                ax1 = subplot(1, 2, 1) ;
                imagesc((1:nU)/nU, tps, divv_apM)
                colormap(bwr256)
                caxis([-climit, climit])
                hold on;
                title('Guess for featureIDs on divergence(v)')
                ax2 = subplot(1, 2, 2) ;
                imagesc((1:nU)/nU, tps, veln_apM)
                colormap(bwr256)
                caxis([-max(abs(veln_apM(:))), max(abs(veln_apM(:)))])
                hold on;
                for pp = 1:nfeatureIDs
                    axes(ax1)
                    plot(featureIDs(pp)/nU * ones(size(tps)), tps) ;
                    axes(ax2)
                    plot(featureIDs(pp)/nU * ones(size(tps)), tps) ;
                end
                axes(ax2)
                title('Guess for featureIDs on normal velocity')
                disp(['Guessed automatic featureIDs to be: [' num2str(featureIDs) ']'])

                % Check it -- is the new feature look good?
                qYN = input(['does feature ' num2str(qq) ' look ok? [y/n]'], 's') ;
                if ~contains(lower(qYN), 'n')
                    qok = true ;
                end
            end
        end
        %% Save featureIDs (valleys of div(v))
        save(fIDfn, 'featureIDs')
    end
    
    % Store in QS
    QS.pathlines.featureIDs.vertices = featureIDs ;
else
    error('Code for this pathlineType here')
end





