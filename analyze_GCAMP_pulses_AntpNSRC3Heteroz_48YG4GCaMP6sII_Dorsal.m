% Analyze the GCAMP pulses for width of anterior fold using 48YG4 x 
% UAS-GCaMP6s(II) in Antp^NS+RC3 heterozygotes

% Use posterior fold as timing and anterior cap as spatial reference
% Alt: use anterior fold as timing and as spatial reference
% NPMitchell 2021


preview = false ;
previewMethod = false ;
overwrite = false ;
resDirFn = 'anteriorMutantHeterozygotes' ;
pausetime = 3 ;

gutMatlabDir = '/mnt/data/code/gut_matlab/' ;
% gutMatlabDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/gut_matlab/' ;

addpath(fullfile(gutMatlabDir, 'addpath_recurse')) ;
addpath_recurse(gutMatlabDir) 

% datdir = '/mnt/data/confocal_data/gut/Mef2GAL4klarUASGCAMP6sIII/analysis_mef2G4kGCaMP6sIII/' ;
% datdir = '/Users/npmitchell/Desktop/gut/GCaMP6sIII_analysis/Mef2GAL4klarGCaMP6sIII_analysis' ;
datdir='/mnt/data/confocal_data/gut/48YGAL4GCaMP6sIIAntpNSRC3/analysis_AntpNSRC3_48YGAL4GCaMP6sII/';

expts = {'202107151200_AntpNSRC3Heteroz_e2',...
    '202107151739_AntpNSRC3Heteroz_e1',...
    '202107161150_AntpNSRC3Heteroz_e3',...
    '202107161818_AntpNSRC3Heteroz_e2', ...
    '202107162200_AntpNSRC3Heteroz_e3', ...
    '202107171200_AntpNSRC3Heteroz_e3', ...
    '202107172106_AntpNSRC3Heteroz_e2', ...
    '202107172106_AntpNSRC3Heteroz_e3', ...
    '202107172300_AntpNSRC3Heteroz_e1', ...
    '202107172300_AntpNSRC3Heteroz_e2', ...
    '202107181040_AntpNSRC3Heteroz_e1', ...
    ... %'202107181556_AntpNSRC3Heteroz_e1', ...  exclude since rocks early
    '202107181556_AntpNSRC3Heteroz_e2', ...
    '202107182209_AntpNSRC3Heteroz_e1', ...
    '202107261411_AntpNSRC3Heteroz_e1', ...     % LONG EXPOSURE 202107191411 e1
    '202107261411_AntpNSRC3Heteroz_e1', ...     % LONG EXPOSURE 202107191411 e2
    '202107261723_AntpNSRC3Heteroz_e1', ...     % LONG EXPOSURE 202107191723 e1
    '202107261723_AntpNSRC3Heteroz_e2', ...     % LONG EXPOSURE 202107191723 e2
    '202107261104_AntpNSRC3Heteroz_e1', ...     % LONG EXPOSURE EXPTS
    '202107261104_AntpNSRC3Heteroz_e2', ...
    } ;

% poster frames: 202106211253_e3, 202106211440_e1, 202106211440_e2

clipY0s = {[70, 70+245], ...    % 202107151200 e2
    [76,76+275], ...              % 202107151739 e1 
    [68,68+255], ...             % 202107161150 e3
    [78,78+236], ...             % 202107161818 e2
    [50,50+280], ...              % 202107162200 e3
    [111,111+240], ...               % 202107171200 e3
    [96,96+242], ...            % 202107172106 e2
    [71,71+250], ...            % 202107172106 e3
    [64,64+268], ...              % 202107172300 e1
    [89,89+250], ...              % 202107172300 e2
    [57,57+245], ...              % 202107181040 e2
    ... %[116, 116+234], ...       % 202107181556 e1
    [99,99+242], ...       % 202107181556 e2
    [90,90+240], ...       % 202107182209 e1
    [], ...                % LONG EXPOSURE 202107191411 e1
    [], ...                % LONG EXPOSURE 202107191411 e2
    [], ...                % LONG EXPOSURE 202107191723 e1
    [], ...                % LONG EXPOSURE 202107191723 e2
    [], ...                % LONG EXPOSURE 202107191104 e1
    [], ...                % LONG EXPOSURE 202107191104 e2
    } ;                
clipYAdjs = {[0, 0], [50,-50]}; 
clipXs = {[1, Inf], [1, Inf], [1, Inf], ...
    [1, Inf], [1, Inf], [1, Inf], [1, Inf], ...
    [1, Inf], [1, Inf], [1, Inf], [1, Inf], ...
    [1, Inf], [1, Inf], [1, Inf], [1, Inf]} ;
% Timestep in MINUTES
dts = 1/60 * [91, ...      % 202107151200 e2
    92, ...         % 202107151739 e1
    92, ...         % 202107161150 e3
    93,  ...        % 202107161818 e2
    93, ...         % 202107162200 e3
    93,  ...        % 202107171200 e3
    93.5, 93.5, ... % 202107172106 e2, e3
    93, 93, ...     % 202107172300 e1, e2
    93, ...         % 202107181040 e2
    ... % 90,       202107181556 e1, 
    90, ...       % 202107181556 e2
    93, ...       % 202107182209 e1, e3
    120, ...                % LONG EXPOSURE 202107191411 e1
    120, ...                % LONG EXPOSURE 202107191411 e2
    120, ...                % LONG EXPOSURE 202107191411 e2
    120, ...                % LONG EXPOSURE 202107191723 e1
    120, ...                % LONG EXPOSURE 202107191723 e2
    120, ...                % LONG EXPOSURE 202107191104 e1
    120, ...                % LONG EXPOSURE 202107191104 e2
    ];

pcPower = [12, 12, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, ...
    20,20,20,20,20,20];
dz = [2.5, 2.5, 3,3,3,3,3,3,3,3,3,3,3, ...
    3,3,3,3,3,3] ;
dates = [15,15,16,16,16,17,17,17,17,17, 18,18,18, ...
    19,19,19,19,19,26,26,] ;

% Conversion in micron / XY pixels
pix2um = [0.331610137, ...  % 202107151200 e2
    0.331611 , ...      % 202107151739 e1 
    0.32999875, ...     % 202107161150 e3
    0.32999787, ...       % 202107161818 e2
    0.3294632, ...      % 202107162200 e3
    0.330669, ...       % 202107171200 e3
    0.331610638, ...    % 202107172106 e2
    0.331610638, ...    % 202107172106 e3
    0.3305353, ...       %202107172300 e1
    0.3305353, ...       %202107172300 e2
    0.3294627, ...       % 202107181040 e2
    ... % 0.33107218, ...     % 202107181556 e1
    0.33107218, ...     % 202107181556 e1
    0.329999, ...       % 202107182209 e1
    [], ...                % LONG EXPOSURE 202107191411 e1
    [], ...                % LONG EXPOSURE 202107191411 e2
    [], ...                % LONG EXPOSURE 202107191411 e3
    [], ...                % LONG EXPOSURE 202107191723 e1
    [], ...                % LONG EXPOSURE 202107191723 e2
    0.30277471839, ...                % LONG EXPOSURE 202107191104 e1
    0.30277471839, ...                % LONG EXPOSURE 202107191104 e2
     ] ;

% FOR POSTERIOR FOLD LATERAL VIEWS
% DEFINE FOLD TIME AS time of 10 micron indentation on ventral side at
% 16*2.5um = 40 um into the tissue -- this is basically the saggital plane
% fold X position at frame foldT, in pixels

% FOR POSTERIOR FOLD DORSAL VIEWS
% position of anterior face when posterior fold self-contact visible
anteriorXs = [147, ...   % 202107151200 e2
    132, ...      % 202107151739 e1 
    171, ...     % 202107161150 e3
    126, ...       % 202107161818 e2
    181, ...       % 202107162200 e3
    175, ...           % 202107171200 e3
    197, ...    % 202107172106 e2
    160, ...    % 202107172106 e3
    159, ...       % 202107172300 e1
    212, ...       % 202107172300 e2
    194, ...       % 202107181040 e2
   ... % 184, ...       % 202107181556 e1
    180, ...       % 202107181556 e2
    164, ...       % 202107182209 e1
    ] ;          
    
% FOR POSTERIOR FOLD DORSAL VIEWS
% time of self-contact visible
% fold frame # at which anterior fold begins & is measured
posteriorTs = [22, ...   % 202107151200 e2
    9, ...       % 202107151739 e1 
    19, ...       % 202107161150 e3
    46, ...       % 202107161818 e2
    31,...          % 202107162200 e3
    13, ...       % 202107171200 e3
    37, ...    % 202107172106 e2
    18, ...    % 202107172106 e3
    5,...      %202107172300 e1
    31, ...      %202107172300 e2
    23, ...       % 202107181040 e2
    ... %17, ...       % 202107181556 e1
    46, ...       % 202107181556 e2
    39, ...       % 202107182209 e1
    ] ;          

% x positions of folds --- just for reference
anteriorFoldXs = [270, ... % 202107151200 e2
    254, ...      % 202107151739 e1 
    286, ...     % 202107161150 e3
    246, ...       % 202107161818 e2
    302, ...       % 202107162200 e3
    309, ...         % 202107171200 e3
    329, ... % 202107172106 e2
    273, ...    % 202107172106 e3
    273, ... % 202107172300 e1 -- previously 274
    354,... % 202107172300 e2
    302, ...       % 202107181040 e2 -- previously 306, rounded muscle
    ... %325, ...       % 202107181556 e1
    308, ...       % 202107181556 e2
    270, ...       % 202107182209 e1
    ] ;  

% fold frame # at which anterior fold begins & is measured
anteriorFoldTs = [24, ...   % 202107151200 e2
    13, ...       % 202107151739 e1 
    25, ...       % 202107161150 e3
    54, ...       % 202107161818 e2
    34, ... % 202107162200 e3
    22, ...       % 202107171200 e3
    41, ...    % 202107172106 e2
    23, ...    % 202107172106 e3
    14,... % 202107172300 e1
    37,...  % 202107172300 e2
    28,...    % 202107181040 e2
    ... %21, ...       % 202107181556 e1
    55, ...       % 202107181556 e2
    47, ...       % 202107182209 e1
    ] ;        

% tOffset is the offset time between morphological stage of posterior
% folding and anterior fold onset in heterozygotes (controls)
toffs = (anteriorFoldTs - posteriorTs) .* dts ;
tOffset = mean(toffs) ;
tOffset_std = std(toffs) ;
tOffset_unc = std(toffs) / sqrt(length(toffs)) ;
% xOffset is the offset distance between anterior fold and anterior face 
% of the midgut
xoffs = (anteriorFoldXs - anteriorXs) .* pix2um ;
xOffset = mean(xoffs) ;
xOffset_std = std(xoffs) ;
xOffset_unc = std(xoffs) / sqrt(length(xoffs)) ;
save(fullfile(datdir, resDirFn, 'spacetime_offsets_controls.mat'), ...
    'xoffs', 'xOffset', 'xOffset_std', 'xOffset_unc', ...
    'toffs', 'tOffset', 'tOffset_std', 'tOffset_unc') ;

xfixed_nonFoldRef = linspace(0, 100, 100/0.33) ;
xfixed_antFoldRef = linspace(-50, 50, 100/0.33) ;
save(fullfile(datdir, resDirFn, 'AntpNSRC3HeterozSettings.mat'), ...
    'clipY0s', 'dts', 'pix2um', 'anteriorXs', 'posteriorTs', ...
    'anteriorFoldXs', 'anteriorFoldTs', 'pcPower', 'dz', 'dates', ...
    'expts', 'xfixed_antFoldRef', 'xfixed_nonFoldRef') ;

nExpts = length(expts) ;
assert(length(anteriorFoldXs) == nExpts)
assert(length(anteriorFoldTs) == nExpts)
assert(length(anteriorXs) == nExpts)
assert(length(posteriorTs) == nExpts)
assert(length(pix2um) == nExpts)
assert(length(clipY0s) == nExpts)

for refID = [2,1]
    
    if refID == 1
        % onset of POSTERIOR folding in minutes
        t0 = dts .* posteriorTs ;
        refStr = '_nonFoldRef' ;
        xfixed = xfixed_nonFoldRef ;
        fixXticks = [0, 25, 50, 75, 100] ;
        xlimFix = [0, 100] ;
    else
        % onset of ANTERIOR folding in minutes
        t0 = dts .* anteriorFoldTs ;
        refStr = '_antFoldRef' ;
        xfixed = xfixed_antFoldRef ;
        fixXticks = [-40, -20, 0, 20, 40] ;
        xlimFix = [-50,50] ;
    end

    % Fixed xlimits in microns
    ylimFix = [] ;
    ylimFixAll = [-15, 90] ;
    cLimits = [0, 10] ;

    % Filtering of drift deltaX options
    capDeltaX = [0, 0, 0, ...
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ;
    medFiltDeltaXW = [0, 0, 0, ...
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ;

    % Antp domain is about 40 microns, broadens on dorsal side.
    close all

    % If drift in X, define it here as [frame#, offsetX] 
    for ee = 1:length(expts)
        fns1 = dir(fullfile(datdir, expts{ee}, 'images', '*c001.png')) ;
        fns2 = dir(fullfile(datdir, expts{ee}, 'images', '*c002.png')) ;
        fns3 = dir(fullfile(datdir, expts{ee}, 'images', '*c003.png')) ;

        % For each yrange
        for clipyPairIdx = 1:length(clipYAdjs)
            if refID == 1
                outfn = sprintf([expts{ee} '_results_Yrange%d' refStr '.mat'], clipyPairIdx) ;
            else
                outfn = sprintf([expts{ee} '_results_Yrange%d' refStr '.mat'], clipyPairIdx) ;
            end
            if ~exist(fullfile(datdir, resDirFn, outfn), 'file') || overwrite
                if ~exist(fullfile(datdir, resDirFn, outfn), 'file')  
                    disp(['results not on disk: ee = ' num2str(ee)])
                else
                    disp(['overwriting: ee = ' num2str(ee)])
                end
                % this y range is defined
                clipY = clipY0s{ee} + clipYAdjs{clipyPairIdx} ;
                clipX = clipXs{ee} ;
                if refID == 1
                    foldX = anteriorXs(ee) ;
                    foldT = posteriorTs(ee) ;
                else
                    foldX = anteriorFoldXs(ee) ;
                    foldT = anteriorFoldTs(ee) ;
                end

                % Look at differences along y for each timepoint
                ntps = length(fns1) ;
                cmap = viridis(ntps) ;
                timestamps = dts(ee) * (0:ntps-1) ;
                timestamps = timestamps - t0(ee) ;

                % Hack
                % if clipyPairIdx == 1
                %     deltaX = zeros(ntps, 1) ;
                % end

                clf
                for ii = 1:ntps
                    im1 = double(imread(fullfile(fns1(ii).folder, fns1(ii).name))) ;
                    im2 = double(imread(fullfile(fns3(ii).folder, fns2(ii).name))) ;
                    im3 = double(imread(fullfile(fns3(ii).folder, fns3(ii).name))) ;

                    if isfinite(clipX(2))
                        c1 = im1(clipY(1):clipY(2), clipX(1):clipX(2)) ;
                        c2 = im2(clipY(1):clipY(2), clipX(1):clipX(2)) ;
                        c3 = im3(clipY(1):clipY(2), clipX(1):clipX(2)) ;
                    else
                        c1 = im1(clipY(1):clipY(2), clipX(1):end) ;
                        c2 = im2(clipY(1):clipY(2), clipX(1):end) ;
                        c3 = im3(clipY(1):clipY(2), clipX(1):end) ;
                    end

                    % mean image
                    curr = (c1 + c2 + c3) / 3 ; 

                    if ii == 50
                        pause(1)
                    end

                    % top-hat TRANSFORM 
                    % First do top-hat FILTERING, then subtract
                    % -- make an elliptical strel
                    % diffA = imgaussfilt(c2, 0.1 / pix2um(ee))
                    % diffB = 
                    % diffC = 
                    % c2fs = gaussian_spectral_filter(c2, 50)
                    
                    % c1fs = gaussianEllipse_spectral_filter(c1, 100, 10) ;
                    % c2fs = gaussianEllipse_spectral_filter(c2, 100, 10) ;
                    % c3fs = gaussianEllipse_spectral_filter(c3, 100, 10) ;
                    % s12 = abs(c1fs - c2fs) ; 
                    % s13 = abs(c1fs - c3fs) ;
                    % s23 = abs(c2fs - c3fs) ;

                    c1fg = imgaussfilt(c1, 1 / pix2um(ee)) ;
                    c2fg = imgaussfilt(c2, 1 / pix2um(ee)) ;
                    c3fg = imgaussfilt(c3, 1 / pix2um(ee)) ;
                    d12 = abs(c1fg - c2fg) ; 
                    d13 = abs(c1fg - c3fg) ;
                    d23 = abs(c2fg - c3fg) ;
                    
                    if ii == 24 && false                
                        figure(1)
                        npanel = 5 ;
                        subplot(npanel, 1, 1)
                        imshow(mat2gray(c2, [1, 40]))
                        subplot(npanel, 1, 2)
                        imshow(mat2gray(c2fs, [0, 40]))
                        subplot(npanel, 1, 3)
                        imshow(mat2gray(c2fg, [0, 40]))
                        subplot(npanel, 1, 4)
                        imshow(mat2gray(d12+d23+d13, [0, 20]))
                        subplot(npanel, 1, 5)
                        imshow(mat2gray(g12+g23+g13, [0, 20]))
                        plot(sum(d12+d23+d13, 1)) ; hold on;
                        plot(sum(g12+g23+g13, 1))
                    end
                    
                    % KEEP BIGGER FEATURES
                    % lOpening = 2 ;  % DV extent of opening
                    % wOpening = 2 ;  % AP extent of opening
                    % se_nlong = round(lOpening / pix2um(ee)) ;
                    % se_nwide = round(wOpening / pix2um(ee)) ;
                    % % se0 = strel('disk', se_nlong) ;
                    % % se_nwide = round(wObjects / pix2um(ee)) ; % 2um diameter / pix2um
                    % % se_sampl = round(size(se0.Neighborhood, 2) / se_nwide) ;
                    % % seN = se0.Neighborhood(:, 1:se_sampl:end) ;
                    % se = strel('arbitrary', ones([se_nlong, se_nwide])) ;
                    % hatt1 = imtophat(d12, se) ;
                    % d12new = d12 - hatt1 ;
                    % figure(1); 
                    % subplot(2, 2, 1); imagesc(d12); axis equal; colorbar
                    % subplot(2, 2, 2); imagesc(hatt1); axis equal; colorbar
                    % subplot(2, 2, 3); imagesc(d12-hatt1); axis equal; colorbar
                    % subplot(2, 2, 4); imagesc(hatt1-d12); axis equal; colorbar

                    % BLUR to KEEP BIGGER FEATURES
                    d12new = imgaussfilt(d12, 0.5 / pix2um(ee)) ;
                    d13new = imgaussfilt(d13, 0.5 / pix2um(ee)) ;
                    d23new = imgaussfilt(d23, 0.5 / pix2um(ee)) ;
                    % d12blur = (d12blur - mean(d12blur) > 0) .* d12blur ;

                    % KEEP SMALLER FEATURES
                    lOpening = 20 ;  % DV extent of opening
                    wOpening = 20 ;  % AP extent of opening
                    se_nlong = round(lOpening / pix2um(ee)) ;
                    se_nwide = round(wOpening / pix2um(ee)) ;
                    % se0 = strel('disk', se_nlong) ;
                    % se_nwide = round(wObjects / pix2um(ee)) ; % 2um diameter / pix2um
                    % se_sampl = round(size(se0.Neighborhood, 2) / se_nwide) ;
                    % seN = se0.Neighborhood(:, 1:se_sampl:end) ;
                    se = strel('arbitrary', ones([se_nlong, se_nwide])) ;
                    hatt12 = imtophat(d12new, se) ;
                    hatt23 = imtophat(d23new, se) ;
                    hatt13 = imtophat(d13new, se) ;


                    if previewMethod && mod(ii, 25) == 0
                        figure(2); 
                        subplot(2, 2, 1); imagesc(d12new); axis equal; colorbar
                        subplot(2, 2, 2); imagesc(hatt12); axis equal; colorbar
                        subplot(2, 2, 3); imagesc(d12new-hatt12); axis equal; colorbar
                        subplot(2, 2, 4); imagesc(hatt12-d12new); axis equal; colorbar
                        sgtitle('d12 and hatt12')
                        pause(pausetime)
                        subplot(2, 2, 1); imagesc(d23new); axis equal; colorbar
                        subplot(2, 2, 2); imagesc(hatt23); axis equal; colorbar
                        subplot(2, 2, 3); imagesc(d23new-hatt23); axis equal; colorbar
                        subplot(2, 2, 4); imagesc(hatt23-d23new); axis equal; colorbar
                        sgtitle('d23 and hatt23')
                        pause(pausetime)
                        subplot(2, 2, 1); imagesc(d13new); axis equal; colorbar
                        subplot(2, 2, 2); imagesc(hatt13); axis equal; colorbar
                        subplot(2, 2, 3); imagesc(d13new-hatt13); axis equal; colorbar
                        subplot(2, 2, 4); imagesc(hatt13-d13new); axis equal; colorbar
                        sgtitle('d13 and hatt13')
                        pause(pausetime)
                        clf
                        plot(sum(hatt12, 1)); hold on; 
                        plot(sum(hatt23, 1)); hold on; 
                        plot(sum(hatt13, 1)); hold on; 
                        sgtitle('hatt12, hatt23, hatt13')
                        pause(pausetime)
                        
                        figure(1); 
                    end

                    if ii == 1 
                        if isfinite(clipX(2))
                            wX = clipX(2)-clipX(1) + 1 ;
                        else
                            wX = size(c1, 2)-clipX(1) + 1 ;
                        end
                        dd = ones(length(fns1), wX) ;
                        ss = ones(length(fns1), wX) ;
                        tophats = ones(length(fns1), wX) ;
                        xshifted = zeros(length(fns1), wX) ;
                        xx = pix2um(ee) * (0:wX-1) ;
                        dinterp = zeros(ntps, size(xfixed, 2)) ; 
                        sinterp = zeros(ntps, size(xfixed, 2)) ;
                        deltaX = zeros(ntps, 1) ;
                    elseif clipyPairIdx < 3
                        % Find offset in x wrt previous image
                        
                        try
                            assert(any(any(abs(curr-prev))))
                        catch
                            error('current and previous frames are identical')
                        end
                        try
                            [deltaX(ii), ~] = ExtPhaseCorrelation(curr, prev) ;
                            LX = clipXs{ee}(2) - clipXs{ee}(1) ;
                            if abs(deltaX(ii)) > capDeltaX(ee)
                                deltaX(ii) = 0 ;
                            end
                            if deltaX(ii) > LX * 0.5 
                                deltaX(ii) = mod(deltaX(ii), LX) ;
                            elseif deltaX(ii) < -LX* 0.5 
                                deltaX(ii) = -mod(abs(deltaX(ii)), LX) ;                    
                            end
                        catch
                            disp('WARNING: could not compute phase corr!')
                            deltaX(ii) = 0 ;
                        end
                    else
                        deltaXfn = sprintf([expts{ee} '_results_Yrange%d.mat'], 2) ;
                        load(fullfile(datdir, resDirFn, deltaXfn), 'deltaX')
                    end

                    % s1 = sum(c1, 1) ;
                    % s2 = sum(c2, 1) ;
                    % s3 = sum(c3, 1) ;                 
                    % d12 = abs(s1 - s2); 
                    % d13 = abs(s1 - s3) ;
                    % d23 = abs(s2 - s3) ;
                    % d12 = abs(sum(c1 - c2, 1)) ; 
                    % d13 = abs(sum(c1 - c3, 1)) ;
                    % d23 = abs(sum(c2 - c3, 1)) ;

                    % Save as image
                    assert(max(hatt12(:)) < 255)
                    diffDir = sprintf(fullfile(fns1(ii).folder, 'diffs_clipY%02d'), clipyPairIdx) ;
                    if ~exist(diffDir, 'dir')
                         mkdir(diffDir)
                     end
                     fn12 = fullfile(diffDir, [fns1(ii).name(1:end-4) '_d12' refStr '.png']) ;
                     fn23 = fullfile(diffDir, [fns1(ii).name(1:end-4) '_d23' refStr '.png']) ;
                     fn13 = fullfile(diffDir, [fns1(ii).name(1:end-4) '_d13' refStr '.png']) ;
                     if ~exist(fn12, 'file') 
                         disp('writing hatt12/23/13')
                         imwrite(uint16(hatt12 * 2^8), fn12) ;
                         imwrite(uint16(hatt23 * 2^8), fn23) ;
                         imwrite(uint16(hatt13 * 2^8), fn13) ;
                     end
                     
                    % Save transient signal in matrix
                    sh12 = sum(hatt12, 1) ;
                    sh23 = sum(hatt23, 1) ; 
                    sh13 = sum(hatt13, 1) ;
                    ddii = mean([sh12; sh23; sh13], 1) ;
                    tophats(ii, :) = ddii ;

                    % Background estimation
                    sii = min([sh12; sh23; sh13], [], 1) ;
                    ss(ii, :) = sii ;
                    dd(ii, :) = ddii - sii ;

                    % Save image for phase correlation with next image
                    prev = curr ;
                end

                % Median filter deltaX
                if medFiltDeltaXW(ee) > 0
                    deltaX = movmedian(deltaX, medFiltDeltaXW(ee)) ;
                end

                % Register the position so that (1) stabilized and (2) 0um is anterior fold
                for ii = 1:ntps
                    % recall dd and ss for this tp
                    ddii = dd(ii, :) ;
                    sii = ss(ii, :) ;

                    % Save shifted ap positions in matrix
                    xshiftii = xx+sum(deltaX(1:ii)) ;
                    xshifted(ii, :) = xshiftii ;
                end
                % Now overall offset
                xshiftPix = xshifted(foldT, foldX) ;
                xshifted = xshifted - xshiftPix ;

                for ii = 1:ntps
                    % recall dd and ss for this tp
                    ddii = dd(ii, :) ;
                    sii = ss(ii, :) ;

                    % Save shifted ap positions in matrix
                    xshiftii = xshifted(ii, :) ;

                    % Interpolate shifted results
                    dinterp(ii, :) = interp1(xshiftii, ddii, xfixed) ;
                    sinterp(ii, :) = interp1(xshiftii, sii, xfixed) ;

                    % Factor in translation of the embryo over the timecourse
                    plot(xx+sum(deltaX(1:ii)), dd(ii, :), 'color', cmap(ii, :)) ; hold on;

                end

                % First plot -- raw curves
                set(gcf, 'Visible', 'off')
                xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
                ylabel('GCAMP transient signal [a.u.]', 'interpreter', 'latex')    
                title('Transient GCAMP activity')
                cb = colorbar ;
                caxis([min(timestamps), max(timestamps)])
                ylabel(cb, 'time [min]')
                outfigfn = sprintf(fullfile(datdir, expts{ee}, ['00_difference_Yrange%d' refStr '.png']), clipyPairIdx) ;
                saveas(gcf, outfigfn) 

                if preview
                    set(gcf, 'visible', 'on')
                    pause(1)
                end
                
                % Kymograph -- pre shift
                clf
                imagesc(xx, timestamps, dd)
                xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
                ylabel('time [min]', 'interpreter', 'latex')
                title('Transient GCAMP activity')
                cb = colorbar ;
                ylabel(cb, '$ \delta I $ [a.u.]', 'interpreter', 'latex')
                outfigfn = sprintf(fullfile(datdir, expts{ee}, ['01_difference_kymo_Yrange%d' refStr '.png']), clipyPairIdx) ;
                saveas(gcf, outfigfn)

                if preview
                    set(gcf, 'visible', 'on')
                    pause(1)
                end
                
                % Shifted kymo
                clf
                for row = 1:ntps
                    imagesc(xshifted(row, :), timestamps(row), dd(row, :));
                    hold on;
                end
                xlabel('stabilized ap position [$\mu$m]', 'interpreter', 'latex')
                ylabel('time [min]', 'interpreter', 'latex')
                title('Transient GCAMP signal', ...
                    'interpreter', 'latex')
                xlim(xlimFix)
                if ~isempty(ylimFix)
                    ylim(ylimFix)
                else
                    ylim([min(timestamps), max(timestamps)])
                end
                cb = colorbar ;
                ylabel(cb, '$ \delta I $ [a.u.]', 'interpreter', 'latex')
                outfigfn = sprintf(fullfile(datdir, expts{ee}, ['02_difference_kymoShift_Yrange%d' refStr '.png']), clipyPairIdx) ;
                saveas(gcf, outfigfn)

                if preview
                    set(gcf, 'visible', 'on')
                    pause(1)
                end

                % PLot deltaX
                clf
                subplot(2, 1, 1)
                plot(timestamps, deltaX, '.-') ;
                ylabel('$\delta x$ [$\mu$m]', 'interpreter', 'latex')
                xlabel('time [min]', 'interpreter', 'latex')
                subplot(2, 1, 2)
                plot(timestamps, cumsum(deltaX), '.-') ;
                ylabel('$\int \delta x$ [$\mu$m]', 'interpreter', 'latex')
                xlabel('time [min]', 'interpreter', 'latex')
                outfigfn = sprintf(fullfile(datdir, expts{ee}, ['02b_shift_deltaX_Yrange%d' refStr '.png']), clipyPairIdx) ;
                saveas(gcf, outfigfn)

                % Plot interpolated results
                clf
                timeGrid = timestamps' .* ones(size(dinterp)) ;
                apGrid = xfixed .* ones(size(dinterp)) ;
                imagesc(xfixed, timestamps, dinterp) ;
                xlim(xlimFix)
                if ~isempty(ylimFix)
                    ylim(ylimFix)
                else
                    ylim([min(timestamps), max(timestamps)])
                end
                xlabel('stabilized, re-gridded ap position [$\mu$m]', 'interpreter', 'latex')
                ylabel('time [min]', 'interpreter', 'latex')
                title('Transient GCAMP signal [registered frame]', ...
                    'interpreter', 'latex')
                cb = colorbar ;
                ylabel(cb, '$ \delta I $ [a.u.]', 'interpreter', 'latex')
                outfigfn = sprintf(fullfile(datdir, expts{ee}, ...
                    ['03_difference_kymoShiftInterp_Yrange%d' refStr '.png']), clipyPairIdx) ;
                saveas(gcf, outfigfn)
                xticks(fixXticks)
                ylim(ylimFixAll)
                saveas(gcf, ...
                    fullfile(datdir, resDirFn, sprintf(['expt%02d_kymo_clipY%d' refStr '.png'], ...
                    ee, clipyPairIdx)))

                if preview
                    set(gcf, 'visible', 'on')
                    pause(1)
                end
                
                % Average over time
                % Get t=0 idx
                [~, tidx0] = min(abs(timestamps)) ;
                end5 = find(timestamps - 5 > 0, 1);
                end10 = find(timestamps - 10 > 0, 1);
                end15 = find(timestamps - 15 > 0, 1);
                end20 = find(timestamps - 20 > 0, 1);
                end25 = find(timestamps - 25 > 0, 1);
                end30 = find(timestamps - 30 > 0, 1);
                end45 = find(timestamps - 45 > 0, 1);
                end60 = find(timestamps - 60 > 0, 1);
                activity00 = mean(dinterp(1:tidx0, :))' ;
                activity05 = mean(dinterp(tidx0:end5, :))' ;
                activity10 = mean(dinterp(tidx0:end10, :))' ;
                activity15 = mean(dinterp(tidx0:end15, :))' ;
                activity20 = mean(dinterp(tidx0:end20, :))' ;
                activity25 = mean(dinterp(tidx0:end25, :))' ;
                activity30 = mean(dinterp(tidx0:end30, :))' ;
                activity45 = mean(dinterp(tidx0:end45, :))' ;
                activity60 = mean(dinterp(tidx0:end60, :))' ;
                legends = {'$\langle t<$0 minutes$\rangle$', ...
                    '$\langle 0<t<$5 minutes$\rangle$', ...
                    '$\langle 0<t<$10 minutes$\rangle$', ...
                    '$\langle 0<t<$15 minutes$\rangle$', ...
                    '$\langle 0<t<$20 minutes$\rangle$', ...
                    '$\langle 0<t<$25 minutes$\rangle$', ...
                    '$\langle 0<t<$30 minutes$\rangle$'} ;
                colors = flipud(viridis(7)) ;

                clf
                plot(xfixed, activity00, 'color', colors(1, :)); hold on;
                plot(xfixed, activity05, 'color', colors(2, :))
                plot(xfixed, activity10, 'color', colors(3, :))
                plot(xfixed, activity15, 'color', colors(4, :))
                plot(xfixed, activity20, 'color', colors(5, :))
                plot(xfixed, activity25, 'color', colors(6, :))
                plot(xfixed, activity30, 'color', colors(7, :))
                legend(legends, 'interpreter', 'latex')
                xlabel('stabilized ap position [$\mu$m]', 'interpreter', 'latex')
                ylabel('GCAMP transient activity [a.u.]', 'interpreter', 'latex')    
                title('Transient GCAMP activity', ...
                    'interpreter', 'latex')
                xlim(xlimFix)
                outfigfn = sprintf(fullfile(datdir, expts{ee}, ['07_difference_meanBgSub_Yrange%d' refStr '.png']), clipyPairIdx) ;
                saveas(gcf, outfigfn) 

                if preview
                    set(gcf, 'visible', 'on')
                    pause(1)
                end
                
                % Filter in space
                close all
                w1um = round(1 / pix2um(ee)) ;
                plot(xfixed, movmedian(activity00, w1um), 'color', colors(1, :)); 
                hold on;
                plot(xfixed, movmedian(activity05, w1um), 'color', colors(2, :)); 
                plot(xfixed, movmedian(activity10, w1um), 'color', colors(3, :)); 
                plot(xfixed, movmedian(activity15, w1um), 'color', colors(4, :)); 
                plot(xfixed, movmedian(activity20, w1um), 'color', colors(5, :)); 
                plot(xfixed, movmedian(activity25, w1um), 'color', colors(6, :)); 
                plot(xfixed, movmedian(activity30, w1um), 'color', colors(7, :)); 
                legend(legends, 'interpreter', 'latex')
                xlabel('stabilized ap position [$\mu$m]', 'interpreter', 'latex')
                ylabel('GCAMP transient activity [a.u.]', 'interpreter', 'latex')    
                title('Transient GCAMP activity, filtered in space 1 $\mu$m', ...
                    'interpreter', 'latex')
                xlim(xlimFix)
                outfigfn = sprintf(fullfile(datdir, ...
                    expts{ee}, ...
                    ['08_difference_meanBgSubSm_Yrange%d' refStr '.png']), clipyPairIdx) ;
                saveas(gcf, outfigfn) 
                saveas(gcf, fullfile(datdir, resDirFn, ...
                    sprintf(['expt%02d_means_clipY%d_%s' refStr '.png'], ...
                    ee, clipyPairIdx, expts{ee})))

                if preview
                    set(gcf, 'visible', 'on')
                    pause(1)
                end
                
                % Save result
                save(fullfile(datdir, resDirFn, outfn), 'xx', 'xshifted', 'dd', 'deltaX', ...
                    'xfixed', 'dinterp', 'sinterp', 'timeGrid', 'apGrid', ...
                    'activity00', 'activity05', 'activity10', ...
                    'activity15', 'activity20', ...
                    'activity25', 'activity30', ...
                    'activity45', 'activity60')
            else
                load(fullfile(datdir, resDirFn, outfn), 'xx', 'xshifted', 'dd', 'deltaX', ...
                    'xfixed', 'dinterp', 'sinterp', 'timeGrid', 'apGrid')
                timestamps = unique(timeGrid)' ;
            end
        end

    end


    %% Average all experiments together
    expts2include = 1:length(expts);
    for clipyPairIdx = 1:2
        colors = define_colors(nExpts) ;
        fixTimeStamps = -20:1.5:90.5 ;  % minutes
        kymoM = zeros(length(xfixed), length(fixTimeStamps)) ;
        nsamples = kymoM ;

        % For different averaging
        avgMin = [15, 20, 25, 30, 45, 60] ;
        allres = [] ;
        for avgID = 1:length(avgMin)
            edmy = 1 ;
            disp(['Performing ' num2str(avgMin(avgID)) ' min average'])
            close all
            statAll = zeros(length(xfixed), length(expts2include)) ;
            for ee = expts2include
                disp(['loading ee = ' num2str(ee)])
                outfn = sprintf([expts{ee} '_results_Yrange%d' refStr '.mat'], clipyPairIdx) ;
                resfn = fullfile(datdir, resDirFn, outfn) ;
                load(resfn, 'xfixed', 'activity15', ...
                    'activity20', 'activity25', 'activity30', ...
                    'activity45', 'activity60', ...
                    'timeGrid', 'dinterp')
                timestamps = unique(timeGrid) ;
                acts = {activity15, activity20, activity25, ...
                    activity30, activity45, activity60} ;

                w1um = round(1 / pix2um(ee)) ;
                activity = movmedian(acts{avgID}, w1um) ;
                padd = round(5/pix2um(ee)) ;

                % no need for this anymore -- all x are xfixed already
                % keep = xfixed > xlimFix(1) & xfixed < xlimFix(2) ;
                % idx2keep = find(keep) ;
                % normIdx = (idx2keep(1) : idx2keep(1)+padd) ;
                % normIdx = [normIdx, idx2keep(end)-padd:idx2keep] ;

                normVal = nanmedian(activity) ;
                bgVal = mean(mink([activity(1:padd); activity(end-padd:end)], 10)) ;
                maxVal = maxk(activity(1+padd:end-padd), 10) ;
                maxVal = mean(maxVal) ;
                anorm = (activity - bgVal) / (maxVal - bgVal); 

                % AVERAGE KYMO
                if avgID == 1
                    dnorm = (dinterp - normVal) / (maxVal - normVal) ;
                    intrp = griddedInterpolant({timestamps, xfixed'}, dnorm, 'nearest', 'none') ;
                    [tt, xx] = ndgrid(fixTimeStamps, xfixed) ;
                    newKymo = intrp(tt, xx) ;
                    nsamples = nsamples + ~isnan(newKymo)' ;
                    newKymo(isnan(newKymo)) = 0 ;
                    kymoM = kymoM + newKymo' ;
                end

                if avgID == 4
                    pause(0.001)
                end
                hold on; 
                plot(xfixed, anorm, 'color', colors(edmy, :)); 
                statAll(:, edmy) = anorm ;
                hold on;

                if ee == expts2include(1)
                    xall = xfixed ;
                end
                % min(xfixed(keep))
                % max(xfixed(keep))

                % Collate all results into table
                % keep4all = ismember(xall, xfixed(keep)) ;
                allres(:, edmy) = anorm ;
                edmy = edmy + 1;
            end
            xlabel('ap position from anterior fold [$\mu$m]', 'interpreter', 'latex')
            ylabel('normalized transient GCAMP activity [a.u.]', ...
                'interpreter', 'latex')
            title(sprintf('transient calcium activity, %d min average', avgMin(avgID)), ...
                'interpreter', 'latex')
            resfn = fullfile(datdir, resDirFn, 'gcamp_activity_results') ;
            saveas(gcf, [resfn sprintf(['_clipY%d_avg%d' refStr '.png'], clipyPairIdx, avgMin(avgID))])
            saveas(gcf, [resfn sprintf(['_clipY%d_avg%d' refStr '.pdf'], clipyPairIdx, avgMin(avgID))])


            % Take stats
            close all
            fig = figure('units','centimeters','position',[0,0,6,10]);

            % kymo panel
            subplot(2, 1, 1)
            kymoM(~isfinite(kymoM)) = NaN;
            kymoMean = (kymoM - min(kymoM(:))) ./ nsamples  ;
            imagesc(xfixed, fixTimeStamps, kymoMean' ); 
            caxis(cLimits)
            % xlabel('ap position from anterior fold [$\mu$m]', 'interpreter', 'latex')
            ylabel('time from fold onset [min]', 'interpreter', 'latex')

            hold on;
            redcolor = colors(2, :) ;
            plot(xlimFix, 0*xlimFix + 5, '--', 'color', redcolor)
            plot(xlimFix, 0*xlimFix - 5, '--', 'color', redcolor)
            plot(xlimFix, 0*xlimFix, '-', 'color', redcolor)
            xticks(fixXticks)
            set(gca, 'xticklabels', [])
            ylim(ylimFixAll)
            colormap(viridis)

            % curve panel
            subplot(2, 1, 2)
            avgact = mean(statAll, 2, 'omitnan') ;
            stdact = std(statAll, 0, 2, 'omitnan') ;
            lineProps = {'-','color', colors(1, :)} ;
            factor = 1.0 / max(avgact(:)) ;
            avgact = movmean(avgact, 3) ;
            stdact = movmedian(stdact, 15) ;
            h1=shadedErrorBar(xfixed, avgact*factor, stdact*factor, 'lineProps', lineProps) ;
            % ylim([min(avgact*factor - stdact*factor), max(avgact*factor + stdact*factor)])
            ylim([-0.25, 1.25])
            xticks(fixXticks)
            xlabel('ap position from anterior fold [$\mu$m]', 'interpreter', 'latex')
            ylabel('GCaMP activity [a.u.]', ...
                'interpreter', 'latex')
            title(sprintf('%d min average', avgMin(avgID)), ...
                'interpreter', 'latex')
            resfn = fullfile(datdir, resDirFn, 'gcamp_mean_results') ;
            saveas(gcf, [resfn sprintf(['_clipY%d_avg%d' refStr '.png'], clipyPairIdx, avgMin(avgID))])
            saveas(gcf, [resfn sprintf(['_clipY%d_avg%d' refStr '.pdf'], clipyPairIdx, avgMin(avgID))])

            % Plot mean kymo for first pass
            if avgID == 1
                close all; 
                chelixMap = cubehelix(128,0.43,-0.68,1.3,0.4,[0,1.0],[0,1.0]) ;
                kymoMean = (kymoM - min(kymoM(:))) ./ nsamples  ;
                imagesc(xfixed, fixTimeStamps, kymoMean' ); 
                cb = colorbar;
                caxis(cLimits)
                xlabel('ap position from anterior fold [$\mu$m]', 'interpreter', 'latex')
                ylabel('time from fold onset [min]', 'interpreter', 'latex')
                ylabel(cb,'transient GCAMP activity [a.u.]')
                xticks(fixXticks)
                resfn = fullfile(datdir, resDirFn, 'gcamp_kymo_results') ;
                colormap(viridis)
                saveas(gcf, [resfn sprintf(['_clipY%d' refStr '.png'], clipyPairIdx)])
                saveas(gcf, [resfn sprintf(['_clipY%d' refStr '.pdf'], clipyPairIdx)])
                colormap(chelixMap)
                saveas(gcf, [resfn sprintf(['_clipY%d_cubeHelix' refStr '.png'], clipyPairIdx)])
                saveas(gcf, [resfn sprintf(['_clipY%d_cubeHelix' refStr '.pdf'], clipyPairIdx)])

            end
        end
    end

    disp('done')
end
