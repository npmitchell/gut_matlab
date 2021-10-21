function outstats = aux_computeQstats(QQ, ang1, MOIeigvalRatio, ...
    keep, areas, strength, eccentricity, ...
    foldt, nLobes, ap_pos, xedges, maxCellSize)
% Compute exhaustive statistics on Qtensors, weighted by weights (areas) as
% well as scaled by nematic strength (strength) or alternative strength
% measure (eccentricity). Stat all cells, but also divide into nLobes lobes
% between foldt =[0, ... 1] and in binned ap positions ap_pos, in bins
% given by xedges. Bound the largest weight by areas = maxCellSize 
%
%
% See also
% --------
% QtensorStats.m
% QtensorAspectTheta.m
%
% Npmitchell 2021



%% prepare weights (typically passed as areas)
allweights = zeros(size(areas)) ;
weight = areas(keep) ;
weights = weight ./ nansum(weight) ;
allweights(keep) = weights ;
    
%% NO WEIGHTS
noweight = 1/length(keep) * ones(length(keep), 1) ;

% Take mean shape from nematic tensor scaled by ar-1
[meanQA, stdQA, steQA] = QtensorStats(strength(keep) .* QQ(keep, :, :), noweight) ;
[meanQAAspect, meanQATheta, ...
    meanQAAspectStd, meanQAThetaStd] = ...
    QtensorAspectTheta(meanQA, stdQA) ;
[~,~, meanQAAspectSte, meanQAThetaSte] = ...
    QtensorAspectTheta(meanQA, steQA) ;

% Take mean shape from nematic tensor scaled by e
[meanQE, stdQE, steQE] = QtensorStats(eccentricity(keep) .* QQ(keep, :, :), noweight) ;
[meanQEAspect, meanQETheta, meanQEAspectStd, meanQEThetaStd] = ...
    QtensorAspectTheta(meanQE, stdQE) ;
[~,~, meanQEAspectSte, meanQEThetaSte] = QtensorAspectTheta(meanQE, steQE) ;


%% Weight by areas
disp('meanQ weighted by area...')
% Take mean shape from nematic tensor scaled by ar-1, weighted
[meanQAW, stdQAW, steQAW] = QtensorStats(strength(keep) .* QQ(keep, :, :), weights) ;
[meanQAAspectW, meanQAThetaW, ...
    meanQAAspectStdW, meanQAThetaStdW] = ...
    QtensorAspectTheta(meanQAW, stdQAW) ;
[~,~, meanQAAspectSteW, meanQAThetaSteW] = ...
    QtensorAspectTheta(meanQAW, steQAW) ;

% Take mean shape from nematic tensor scaled by e, weighted
[meanQEW, stdQEW, steQEW] = QtensorStats(eccentricity(keep) .* QQ(keep, :, :), weights) ;
[meanQEAspectW, meanQEThetaW, meanQEAspectStdW, meanQEThetaStdW] = ...
    QtensorAspectTheta(meanQEW, stdQEW) ;
[~,~, meanQEAspectSteW, meanQEThetaSteW] = ...
    QtensorAspectTheta(meanQEW, steQEW) ;

%% area-weighted, roll off weights for largest cells
disp('meanQ weighted by bounded area...')
weight(weight > 0.5 * maxCellSize) = maxCellSize - weight(weight > 0.5 * maxCellSize) ;
weights = weight ./ nansum(weight) ;
allweightsBounded = zeros(size(areas)) ;
allweightsBounded(keep) = weights ;

% Take mean shape from nematic tensor scaled by ar-1, weighted
[meanQAWB, stdQAWB, steQAWB] = QtensorStats(strength(keep) .* QQ(keep, :, :), weights) ;
[meanQAAspectWB, meanQAThetaWB, ...
    meanQAAspectStdWB, meanQAThetaStdWB] = ...
    QtensorAspectTheta(meanQAWB, stdQAWB) ;
[~,~, meanQAAspectSteWB, meanQAThetaSteWB] = ...
    QtensorAspectTheta(meanQAWB, steQAWB) ;

% Take mean shape from nematic tensor scaled by e, weighted
[meanQEWB, stdQEWB, steQEWB] = QtensorStats(eccentricity(keep) .* QQ(keep, :, :), weights) ;
[meanQEAspectWB, meanQEThetaWB, ...
    meanQEAspectStdWB, meanQEThetaStdWB] = ...
    QtensorAspectTheta(meanQEWB, stdQEWB) ;
[~,~, meanQEAspectSteWB, meanQEThetaSteWB] = ...
    QtensorAspectTheta(meanQEWB, steQEWB) ;

% Other measures
disp('other measures like cos2theta and sin2theta...')
% c1 = sqrt(mratio_principal(:)) .* cos(2 * ang1(keep)) ;
% s1 = sqrt(mratio_principal(:)) .* sin(2 * ang1(keep)) ;
% ar = vecnorm([mean(c1), mean(s1)]) ;
% theta = 0.5 * atan2(nanmean(s1), nanmean(c1)) ;
mratio_principal = MOIeigvalRatio(keep) ;
ars = sqrt(mratio_principal(:)) ;
[aspectMeanWeighted, aspectStdWeighted, aspectSteWeighted] = ...
    weightedMeanStdSte(ars, weights) ;
thetas = ang1(keep) ;
cos2thetas = cos(2 * ang1(keep)) ;
sin2thetas = sin(2 * ang1(keep)) ;
[mt,st,et1,et2 ] = angleStats(2*thetas, weights) ;
meanThetaWeighted =mt*0.5 ;
stdThetaWeighted = st*0.5 ;
steThetaWeighted1 = et1*0.5 ;
steThetaWeighted2 = et2*0.5 ;
[cos2thetaMeanWeighted, cos2thetaStdWeighted, cos2thetaSteWeighted] = ...
    weightedMeanStdSte(cos2thetas, weights) ;
[sin2thetaMeanWeighted, sin2thetaStdWeighted, sin2thetaSteWeighted] = ...
    weightedMeanStdSte(sin2thetas, weights) ;

% The aspect ratio is related to the meanQ 
% For perfectly aligned cells, meanQ would have a norm(Q) = 0.5, so
% abs(norm(meanQ)) * 2 = |ar| - 1
try
    assert(abs(sqrt(abs(4 * det(meanQAWB))) - 2 * norm(meanQAWB)) < 1e-7)
    assert(abs(sqrt(abs(4 * det(meanQEWB))) - 2 * norm(meanQEWB)) < 1e-7)
catch
    error('meanQ does not satisfy traceless norm condition')
end

%% Statistics by AP position
disp('statistics by AP position...')
[mid_ap, mean_QAc2t_ap, std_QAc2t_ap, ~, ste_QAc2t_ap] = ...
    binDataMeanStdWeighted(ap_pos, strength(keep) .* cos2thetas, ...
        xedges, weights) ;
[mid_ap, mean_QAs2t_ap, std_QAs2t_ap, ~, ste_QAs2t_ap] = ...
    binDataMeanStdWeighted(ap_pos, strength(keep) .* sin2thetas, ...
        xedges, weights) ;

% eccentricity-weighted
[mid_ap, mean_QEc2t_ap, std_QEc2t_ap, ~, ste_QEc2t_ap] = ...
    binDataMeanStdWeighted(ap_pos, eccentricity(keep) .* cos2thetas, ...
        xedges, weights) ;
[mid_ap, mean_QEs2t_ap, std_QEs2t_ap, ~, ste_QEs2t_ap] = ...
    binDataMeanStdWeighted(ap_pos, eccentricity(keep) .* sin2thetas, ...
        xedges, weights) ;

%% Statistics by Lobe (between features.folds)    
disp('statistics by lobe...')
% Statistics by Lobe -- ars.*cos(2theta), ars.*sin(2theta)
[~, mean_QAc2t_lobes, std_QAc2t_lobes, ~, ste_QAc2t_lobes] = binDataMeanStdWeighted(ap_pos, ...
    strength(keep) .* cos(2 * ang1(keep)), foldt, weights) ;
[~, mean_QAs2t_lobes, std_QAs2t_lobes, ste_QAs2t_lobes] = binDataMeanStdWeighted(ap_pos, ...
    strength(keep) .* sin(2 * ang1(keep)), foldt, weights) ;
[~, mean_QEc2t_lobes, std_QEc2t_lobes, ~, ste_QEc2t_lobes] = binDataMeanStdWeighted(ap_pos, ...
    strength(keep) .* cos(2 * ang1(keep)), foldt, weights) ;
[~, mean_QEs2t_lobes, std_QEs2t_lobes, ste_QEs2t_lobes] = binDataMeanStdWeighted(ap_pos, ...
    strength(keep) .* sin(2 * ang1(keep)), foldt, weights) ;

% Other measures
[~, lobes_Q11A, lobes_std_Q11A, ~, lobes_ste_Q11A] = ...
    binDataMeanStdWeighted(ap_pos, ...
    strength(keep) .* squeeze(QQ(keep, 1, 1)), foldt, weights) ;
[~, lobes_Q12A, lobes_std_Q12A, ~, lobes_ste_Q12A] = binDataMeanStdWeighted(ap_pos, ...
    strength(keep) .* squeeze(QQ(keep, 1, 2)), foldt, weights) ;
[~, lobes_Q21A, lobes_std_Q21A, ~, lobes_ste_Q21A] = binDataMeanStdWeighted(ap_pos, ...
    strength(keep) .* squeeze(QQ(keep, 2, 1)), foldt, weights) ;
[~, lobes_Q22A, lobes_std_Q22A, ~, lobes_ste_Q22A] = binDataMeanStdWeighted(ap_pos, ...
    strength(keep) .* squeeze(QQ(keep, 2, 2)), foldt, weights) ;

% Eccentricity-weighted
[~, lobes_Q11E, lobes_std_Q11E, ~, lobes_ste_Q11E] = ...
    binDataMeanStdWeighted(ap_pos, ...
    eccentricity(keep) .* squeeze(QQ(keep, 1, 1)), foldt, weights) ;
[~, lobes_Q12E, lobes_std_Q12E, ~, lobes_ste_Q12E] = ...
    binDataMeanStdWeighted(ap_pos, ...
    eccentricity(keep) .* squeeze(QQ(keep, 1, 2)), foldt, weights) ;
[~, lobes_Q21E, lobes_std_Q21E, ~, lobes_ste_Q21E] = ...
    binDataMeanStdWeighted(ap_pos, ...
    eccentricity(keep) .* squeeze(QQ(keep, 2, 1)), foldt, weights) ;
[~, lobes_Q22E, lobes_std_Q22E, ~, lobes_ste_Q22E] = ...
    binDataMeanStdWeighted(ap_pos, ...
    eccentricity(keep) .* squeeze(QQ(keep, 2, 2)), foldt, weights) ;

% Check that result is still traceless and symmetric
assert(all(abs(lobes_Q11A + lobes_Q22A) < 1e-7))
assert(all(abs(lobes_Q12A - lobes_Q12A) < 1e-7))
assert(all(abs(lobes_Q11E + lobes_Q22E) < 1e-7))
assert(all(abs(lobes_Q12E - lobes_Q12E) < 1e-7))

% Collate lobe/chamber information
for AorE = 1:2
    meanQLobeAspect = zeros(nLobes, 1) ;
    meanQLobeAspectStd = zeros(nLobes, 1) ;
    meanQLobeAspectSte = zeros(nLobes, 1) ;
    meanQLobeTheta = zeros(nLobes, 1) ;
    meanQLobeThetaStd = zeros(nLobes, 1) ;
    meanQLobeThetaSte = zeros(nLobes, 1) ;
    for lobe = 1:nLobes
        if AorE == 1
            % aspect-weighted
            lobes_Q11 = lobes_Q11A ;
            lobes_Q12 = lobes_Q12A ;
            lobes_Q21 = lobes_Q21A ;
            lobes_Q22 = lobes_Q22A ;
            lobes_std_Q11 = lobes_std_Q11A ;
            lobes_std_Q12 = lobes_std_Q12A ;
            lobes_std_Q21 = lobes_std_Q21A ;
            lobes_std_Q22 = lobes_std_Q22A ;
            lobes_ste_Q11 = lobes_ste_Q11A ;
            lobes_ste_Q12 = lobes_ste_Q12A ;
            lobes_ste_Q21 = lobes_ste_Q21A ;
            lobes_ste_Q22 = lobes_ste_Q22A ;
        else
            % eccentricity-weighted
            lobes_Q11 = lobes_Q11E ;
            lobes_Q12 = lobes_Q12E ;
            lobes_Q21 = lobes_Q21E ;
            lobes_Q22 = lobes_Q22E ;
            lobes_std_Q11 = lobes_std_Q11E ;
            lobes_std_Q12 = lobes_std_Q12E ;
            lobes_std_Q21 = lobes_std_Q21E ;
            lobes_std_Q22 = lobes_std_Q22E ;
            lobes_ste_Q11 = lobes_ste_Q11E ;
            lobes_ste_Q12 = lobes_ste_Q12E ;
            lobes_ste_Q21 = lobes_ste_Q21E ;
            lobes_ste_Q22 = lobes_ste_Q22E ;
        end
        meanQ_lobes{lobe} = [lobes_Q11(lobe), lobes_Q12(lobe); ...
            lobes_Q21(lobe), lobes_Q22(lobe)] ;
        stdQ_lobes{lobe} = [lobes_std_Q11(lobe), lobes_std_Q12(lobe); ...
            lobes_std_Q21(lobe), lobes_std_Q22(lobe)] ;
        steQ_lobes{lobe} = [lobes_ste_Q11(lobe), lobes_ste_Q12(lobe); ...
            lobes_ste_Q21(lobe), lobes_ste_Q22(lobe)] ;

        % diagonalize this lobeQ
        [meanQLobeAspect(lobe), meanQLobeTheta(lobe), ...
            meanQLobeAspectStd(lobe), meanQLobeThetaStd(lobe)] ...
            = QtensorAspectTheta(meanQ_lobes{lobe}, stdQ_lobes{lobe}) ;
        [ ~, ~, meanQLobeAspectSte(lobe), meanQLobeThetaSte(lobe)] ...
            = QtensorAspectTheta(meanQ_lobes{lobe}, steQ_lobes{lobe}) ;
    end

    if AorE == 1
        % aspect-weighted result
        meanQALobeAspect = meanQLobeAspect ;
        meanQALobeAspectStd = meanQLobeAspectStd ;
        meanQALobeAspectSte = meanQLobeAspectSte ;
        meanQALobeTheta = meanQLobeTheta ;
        meanQALobeThetaStd = meanQLobeThetaStd ;
        meanQALobeThetaSte = meanQLobeThetaSte ;
        meanQA_lobes = meanQ_lobes ;
        stdQA_lobes = stdQ_lobes ;
        steQA_lobes = steQ_lobes ;
    else
        % eccentricity-weighted result
        meanQELobeAspect = meanQLobeAspect ;
        meanQELobeAspectStd = meanQLobeAspectStd ;
        meanQELobeAspectSte = meanQLobeAspectSte ;
        meanQELobeTheta = meanQLobeTheta ;
        meanQELobeThetaStd = meanQLobeThetaStd ;
        meanQELobeThetaSte = meanQLobeThetaSte ;
        meanQE_lobes = meanQ_lobes ;
        stdQE_lobes = stdQ_lobes ;
        steQE_lobes = steQ_lobes ;
    end
end


%% Save  in struct
% STORE NEMATIC INFO IN QUALITIES
seg3d.qualities.nematicTensor = QQ ;
seg3d.qualities.nematicStrength = strength ; % aspect ratio - 1.
seg3d.qualities.eccentricity = eccentricity ;


outstats = struct() ;
outstats.keep = keep ;
outstats.maxCellSize = maxCellSize ;
outstats.weights = allweights ;
outstats.weightsBounded = allweightsBounded ;

% Raw distributions    
outstats.aspectMean = aspectMeanWeighted ;
outstats.aspectMedian = median(ars) ;
outstats.aspectStd = aspectStdWeighted ;
outstats.aspectSte = aspectSteWeighted ;
outstats.aspect25 = prctile(ars, 25.0) ;
outstats.aspect75 = prctile(ars, 75.0) ;
outstats.thetaMean = meanThetaWeighted ;    
outstats.thetaStd = stdThetaWeighted ;    
outstats.thetaSte = max(steThetaWeighted1, steThetaWeighted2) ;

% Note: Multiply by two for aspect since norm gives 1/2*(strength_Q)
% mean_aspectRatio_A = [mean_aspectRatio_A, norm(meanQAW) * 2 + 1 ] ;    
% mean_aspectRatio_E = [mean_aspectRatio_E, 1/(sqrt(1 - (norm(meanQEW)*2)^2))  ] ;

outstats.cos2thetaMean = cos2thetaMeanWeighted ;
outstats.sin2thetaMean = sin2thetaMeanWeighted ;
outstats.cos2thetaStd = cos2thetaStdWeighted ;
outstats.cos2thetaSte = cos2thetaSteWeighted ;
outstats.sin2thetaStd = sin2thetaStdWeighted ;
outstats.sin2thetaSte = sin2thetaSteWeighted ;
outstats.cos2theta25 = prctile(cos2thetas, 25.0) ;
outstats.cos2theta75 = prctile(cos2thetas, 75.0) ;
outstats.sin2theta25 = prctile(sin2thetas, 25.0) ;
outstats.sin2theta75 = prctile(sin2thetas, 75.0) ;

%% Mean tensor stats

outstats.meanQ = struct() ;
% weighted by (aspect ratio - 1) "cell anisotropy tensor"
QAStats = struct() ;
QAStats.meanQ = meanQA ;
QAStats.meanQWeighted = meanQAW ;
QAStats.meanQWeightBounded = meanQAWB ;
QAStats.meanQAspect = meanQAAspect ;
QAStats.meanQAspectStd = meanQAAspectStd ;
QAStats.meanQAspectSte = meanQAAspectSte ;
QAStats.meanQAspectWeighted = meanQAAspectW ;
QAStats.meanQAspectStdWeighted = meanQAAspectStdW ;
QAStats.meanQAspectSteWeighted = meanQAAspectSteW ;
QAStats.meanQAspectWeightBounded = meanQAAspectWB ;
QAStats.meanQAspectStdWeightBounded = meanQAAspectStdWB ;
QAStats.meanQAspectSteWeightBounded = meanQAAspectSteWB ;
QAStats.meanQTheta = meanQATheta ;
QAStats.meanQThetaStd = meanQAThetaStd ;
QAStats.meanQThetaSte = meanQAThetaSte ;
QAStats.meanQThetaWeighted = meanQAThetaW ;
QAStats.meanQThetaStdWeighted = meanQAThetaStdW ;
QAStats.meanQThetaSteWeighted = meanQAThetaSteW ;
QAStats.meanQThetaWeightBounded = meanQAThetaWB ;
QAStats.meanQThetaStdWeightBounded = meanQAThetaStdWB ;
QAStats.meanQThetaSteWeightBounded = meanQAThetaSteWB ;
outstats.meanQ.aspectWeighted = QAStats ;

% weighted by eccentricity "cell eccentricity tensor"
QEStats = struct() ;
QEStats.meanQ = meanQE ;
QEStats.meanQWeighted = meanQEW ;
QEStats.meanQWeightBounded = meanQEWB ;
QEStats.meanQAspect = meanQEAspect ;
QEStats.meanQAspectStd = meanQEAspectStd ;
QEStats.meanQAspectSte = meanQEAspectSte ;
QEStats.meanQAspectWeighted = meanQEAspectW ;
QEStats.meanQAspectStdWeighted = meanQEAspectStdW ;
QEStats.meanQAspectSteWeighted = meanQEAspectSteW ;
QEStats.meanQAspectWeightBounded = meanQEAspectWB ;
QEStats.meanQAspectStdWeightBounded = meanQEAspectStdWB ;
QEStats.meanQAspectSteWeightBounded = meanQEAspectSteWB ;
QEStats.meanQTheta = meanQETheta ;
QEStats.meanQThetaStd = meanQEThetaStd ;
QEStats.meanQThetaSte = meanQEThetaSte ;
QEStats.meanQThetaWeighted = meanQEThetaW ;
QEStats.meanQThetaStdWeighted = meanQEThetaStdW ;
QEStats.meanQThetaSteWeighted = meanQEThetaSteW ;
QEStats.meanQThetaWeightBounded = meanQEThetaWB ;
QEStats.meanQThetaStdWeightBounded = meanQEThetaStdWB ;
QEStats.meanQThetaSteWeightBounded = meanQEThetaSteWB ;
outstats.meanQ.eccentricityWeighted = QEStats ;

% Mean Q tensor (weightedBounded) for each lobe
outstats.lobes = struct() ;
% aspect-weighted Q
lobeQAStats = struct() ;
lobeQAStats.meanQLobes = meanQA_lobes ;
lobeQAStats.stdQLobes = stdQA_lobes ;
lobeQAStats.steQLobes = steQA_lobes ;
lobeQAStats.meanQLobeAspect = meanQALobeAspect ;
lobeQAStats.meanQLobeAspectStd = meanQALobeAspectStd ;
lobeQAStats.meanQLobeAspectSte = meanQALobeAspectSte ;
lobeQAStats.meanQLobeTheta = meanQALobeTheta ;
lobeQAStats.meanQLobeThetaStd = meanQALobeThetaStd ;
lobeQAStats.meanQLobeThetaSte = meanQALobeThetaSte ;
lobeQAStats.meanQCos2Theta = mean_QAc2t_lobes ;
lobeQAStats.meanQSin2Theta = mean_QAs2t_lobes ;
lobeQAStats.stdQCos2Theta = std_QAc2t_lobes ;
lobeQAStats.stdQSin2Theta = std_QAs2t_lobes ;
lobeQAStats.steQCos2Theta = ste_QAc2t_lobes ;
lobeQAStats.steQSin2Theta = ste_QAs2t_lobes ;
outstats.lobes.aspectWeighted = lobeQAStats ;

% eccentricity-weighted Q
lobeQEStats = struct() ;
lobeQEStats.meanQLobes = meanQE_lobes ;
lobeQEStats.stdQLobes = stdQE_lobes ;
lobeQEStats.steQLobes = steQE_lobes ;
lobeQEStats.meanQLobeAspect = meanQELobeAspect ;
lobeQEStats.meanQLobeAspectStd = meanQELobeAspectStd ;
lobeQEStats.meanQLobeAspectSte = meanQELobeAspectSte ;
lobeQEStats.meanQLobeTheta = meanQELobeTheta ;
lobeQEStats.meanQLobeThetaStd = meanQELobeThetaStd ;
lobeQEStats.meanQLobeThetaSte = meanQELobeThetaSte ;
lobeQEStats.meanQCos2Theta = mean_QEc2t_lobes ;
lobeQEStats.meanQSin2Theta = mean_QEs2t_lobes ;
lobeQEStats.stdQCos2Theta = std_QEc2t_lobes ;
lobeQEStats.stdQSin2Theta = std_QEs2t_lobes ;
lobeQEStats.steQCos2Theta = ste_QEc2t_lobes ;
lobeQEStats.steQSin2Theta = ste_QEs2t_lobes ;
outstats.lobes.eccentricityWeighted = lobeQEStats ;

%% AP averaging
outstats.apStats = struct() ;

% anisotropy-weighted
QAapStats = struct() ;
QAapStats.apBins = mid_ap ;
QAapStats.apCos2Theta = mean_QAc2t_ap ;
QAapStats.apSin2Theta = mean_QAs2t_ap ;
QAapStats.apCos2ThetaStd = std_QAc2t_ap ;
QAapStats.apSin2ThetaStd = std_QAs2t_ap ;
QAapStats.apCos2ThetaSte = ste_QAc2t_ap ;
QAapStats.apSin2ThetaSte = ste_QAs2t_ap ;
outstats.apStats.aspectWeighted = QAapStats ;

% eccentricity-weighted
QEapStats = struct() ;
QEapStats.apBins = mid_ap ;
QEapStats.apCos2Theta = mean_QEc2t_ap ;
QEapStats.apSin2Theta = mean_QEs2t_ap ;
QEapStats.apCos2ThetaStd = std_QEc2t_ap ;
QEapStats.apSin2ThetaStd = std_QEs2t_ap ;
QEapStats.apCos2ThetaSte = ste_QEc2t_ap ;
QEapStats.apSin2ThetaSte = ste_QEs2t_ap ;
outstats.apStats.eccentricityWeighted = QEapStats ;

    