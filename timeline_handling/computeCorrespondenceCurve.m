function [corrPath, corrRaw, corrError, ssrPath, ssrPathError] = ...
    computeCorrespondenceCurve(ssrM, minddssr, minerror, options)
%[corrPath, corrRaw, corrError, ssrPath, ssrPathError] = ...
%   computeCorrespondenceCurve(ssrM, minddssr, minerror, options)
%
% Parameters
% ----------
% ssrM 
% minddssr 
% minerror  
% options : struct with fields
%   overwrite : bool (default=false)
%   ssrDir : string path to ssr output dir
%   timestamps : maxNtp x MaxNtp numeric
%       timestamps for comparison timelines, in collated array, such that
%       timestamps(cc, 1:ntps(cc)) gives the comparison timeline cc, while
%       timestamps(refID, :) gives the reference timeline
%   refExptID : string
%   cExptID : string
%   sigmaTime : float
%   rsubsampling : int
%       subsampling factor of the reference timeline 
%   
%
% Returns
% -------
% corrPath :
% corrRaw :
% corrError :
% ssrPath :
% ssrPathError :
%
% NPMitchell 2020-2022

% Default options
overwrite = false ;

% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
cc = options.cc ;
refID = options.refID ;
ssrDir = options.ssrDir ;
timestamps = options.timestamps ;
ntps = options.ntps ;
refExptID = options.refExptID ;
cExptID = options.cExptID ;
sigmaTime = options.sigmaTime ;
rsubsampling = options.rsubsampling ;

%% process options
extn = [sprintf('_ss%02d', rsubsampling) '_ssR'] ;


% Now ssrM is fully filled in
corrPathFn = fullfile(ssrDir, ...
    [sprintf('correspondencePath_r%s_c%s', refExptID, cExptID), ...
    extn '.mat']) ;

if ~exist(corrPathFn, 'file') || overwrite
    clf
    imagesc(timestamps(cc, 1:ntps(cc))/60, timestamps(refID, :)/60, ssrM') ;
    axis equal
    axis tight
    cb = colorbar() ;
    ylabel(cb, '$\sqrt{\langle \sigma_{ij} \rangle \langle \sigma_{ji} \rangle}$', ...
        'interpreter', 'latex')
    xlabel(['time ' cExptID ' [hr]'], 'interpreter', 'latex')
    ylabel(['time ' refExptID ' [hr]'], 'interpreter', 'latex')
    saveas(gcf, fullfile( ssrDir, ...
        sprintf(['ssr_heatmap_c%s_r%s' extn '.png'], cExptID, refExptID)))
    save(fullfile( ssrDir, ...
        sprintf(['ssr_c%s_r%s' extn '.mat'], cExptID, refExptID)), ...
        'ssrM')

    % Get shortest path
    % Define correspondence pairs like in dynamicAtlas
    ssr4path = imgaussfilt(-ssrM, sigmaTime / rsubsampling) ;
    ssr4path = ssr4path - min(ssr4path(:)) ;
    pathOpts = struct('exponent', 1.0) ;
    corrPath = shortestPathInImage(ssr4path, pathOpts) ;

    % Interpolate to get measure of the ssR along the path
    % imagesc(timestamps(cc, 1:ntps(cc)), timestamps(1,:), ssrM')
    [tc, tr] = meshgrid(timestamps(cc, 1:ntps(cc))-timestamps(cc,1), ...
        timestamps(1,:)- timestamps(1, 1)) ;
    ssrPath = interp2(tc, tr, ssrM', corrPath(:, 1), corrPath(:, 2), 'linear') ;

    % Save corrPaths and save image
    corrRaw = minddssr(cc, 1:ntps(cc)) ;
    corrError = movmean(minerror(cc, 1:ntps(cc)), 5) ;
    % Check if we need to truncate the correpondence path at start
    % or end
    if length(corrPath) < length(corrRaw)
        if corrPath(1, 1) == 1
            rawID = 1:length(corrPath) ;
            corrRaw = corrRaw(rawID) ;
            corrError = corrError(rawID) ;
        else
            rawID = (length(corrRaw)-length(corrPath)):length(corrRaw) ;
            corrRaw = corrRaw(rawID) ;
            corrError = corrError(rawID) ;
        end
    else
        rawID = 1:length(corrRaw) ;
    end

    % Estimate uncertainty by difference between ssR at raw and
    % smoothed paths
    rawcorrSSRs = [] ;
    for ptID = 1:length(corrRaw)
        rrowP = corrRaw(ptID) ;
        ccolP = rawID(ptID) ;
        rawcorrSSRs(ptID) = ssrM(ccolP,rrowP) ;
    end
    ssrPathError = abs(ssrPath - rawcorrSSRs') ;

    % Visualization
    clf
    imagesc(ssrM)
    hold on;
    plot(corrPath(:, 2), corrPath(:, 1), 'o') ;
    plot(corrRaw, rawID, 'k.') 
    plot(corrPath(:, 2) -corrError', corrPath(:, 1), 'k-')
    plot(corrPath(:, 2) +corrError', corrPath(:, 1), 'k-')
    cb = colorbar() ;
    ylabel(cb, '$\langle \sigma_{ij} \rangle \langle \sigma_{ji} \rangle$', ...
        'interpreter', 'latex')  
    colormap(viridis_r)
    axis equal 
    axis tight
    xlabel(['time ' refExptID ' [min]'], 'interpreter', 'latex')
    ylabel(['time ' cExptID ' [min]'], 'interpreter', 'latex')
    corrPathFigFn = fullfile(ssrDir, ...
        sprintf(['correspondencePath_r%s_c%s' extn '.pdf'], ...
        refExptID, cExptID)) ;
    saveas(gcf, corrPathFigFn)

    % SAVE DATA RESULT
    save(corrPathFn, 'corrPath', 'corrRaw', 'corrError', 'ssrPath', 'ssrPathError') ;
else
    load(corrPathFn, 'corrPath', 'corrRaw', 'corrError', 'ssrPath', 'ssrPathError') ;
end