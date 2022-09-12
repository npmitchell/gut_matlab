function [h1, h2, q1, q2, cb] = plotTensorFieldCrossesOnImage(im, xx, yy, QQ, options)
% [h1, h2, q1, q2, cb] = plotTensorFieldCrossesOnImage(im, xx, yy, QQ, options)
%   Plot a tensor field (Qxx,Qxy;Qyx,Qyy) evaluated at grid[xx, yy] on an 
%   image im as quiverplot (subsampled quiver). 
%   Plots a heatmap as plotStyle. 
%   Note that you can incorporate a metric in options and can specify the 
%   resolution of the subsampled quiver plot. 
%   
%
% Parameters
% ----------
% im : PxQx1 or PxQx3 numeric
%   image over which to plot tensor crosses
% xx : N x 1 float array
%   x values of tensor evaluation points
% yy : M x 1 float array
%   y values of tensor evaluation points
% QQ : NxM x 2x2 float array
%   tensor field defined across xx and yy in x direction
% qopts : struct with fields
%   overlayHeatmap : bool
%       overlay a heatmap on top of the image based on the tensors:
%       dev/mag/theta
%   label : str
%       colorbar label. Default is '$||||$' 
%   qsubsample : 2x1 int
%       final size to which to subsample the supplied tensor field for
%       crosses
%   qscale : float
%       overall scale of the quivers
%   outfn : str
%       output filename for figure as png 
%   colorTensors : bool
%       color the vectors by qopts.colors (default parula colors)
%   climv : length 2 float
%       magnitude associated with min/maximum color/intensity in image
%   colors : Rx3 rgb value array for tensors
%       if colorTensors==true, use these colors for the tensors
%   eig1color : color specifier (eg, 'r' or [1,0,0])
%       if colorTensors==false, use this color to color the first bar
%   eig2color : color specifier (eg, 'r' or [1,0,0])
%       if colorTensors==false, use this color to color the second bar
%   metric : NxM x 2x2
%       metric tensor for each element of the tensor field, for taking
%       norm
%
%
% Returns
% -------
% h1 : handle for imshow
% h2 : handle for quiverplot long axis of tensors
% h3 : handle for quiverplot short axis of tensors
%
% Example Usage 
% -------------
%
% See also
% --------
% PLOTNEMATICFIELD()
% vectorFieldQuiverOnImage()
%
% 
%
%
% NPMitchell 2022

% Default options
fig = 'A' ;  % initialize as something not a handle
overlayHeatmap = false ;
plotStyle = 'mag' ;
qsubsample = 10 ;
qscale = 10 ;
transposeVectorShape = false ;  % Consider using this if there is some meshgrid/ndgrid monkeybusiness
lw = 0.5 ;
colorTensors = false ;
makeCbar = true ;
axisOff = true ;
eig1color = [ 0.929, 0.694, 0.125] ;
eig2color = [ 0.929, 0.694, 0.125] ;
clim = [] ;

% Unpack options
if isfield(options, 'plotStyle')
    plotStyle = options.plotStyle ;
end
if isfield(options, 'overlayHeatmap')
    overlayHeatmap = options.overlayHeatmap ;
end
if isfield(options, 'qsubsample')
    qsubsample = options.qsubsample ;
end
if isfield(options, 'qscale') 
    qscale = options.qscale ;
end
if isfield(options, 'eig1color')
    eig1color = options.eig1color ;
end
if isfield(options, 'eig2color')
    eig1color = options.eig1color ;
end
if isfield(options, 'fig') 
    fig = options.fig ;  % figure on which to plot, if supplied
end
if isfield(options, 'colorVectors') 
    colorTensors = options.colorVectors ;
end
if isfield(options, 'lw') 
    lw = options.lw ;
end
if isfield(options, 'fig')
    fig = options.fig ;
end
if isfield(options, 'clim')
    clim = options.clim ;  % for heatmapQ, not for image. For image LUT, use RGB image.
end
if isfield(options, 'makeCbar')
    makeCbar = options.makeCbar ;
end
if isfield(options, 'label')
    labelstr = options.label ;
elseif strcmpi(plotStyle, 'mag')
    labelstr = '||Q||' ;
elseif strcmpi(plotStyle, 'dev')
    labelstr = '||Dev(Q)||' ;    
elseif strcmpi(plotStyle, 'tr')
    labelstr = '||Tr(Q)||' ;    
elseif strcmpi(plotStyle, 'theta')
    labelstr = '||arg[Dev(Q)]||' ;    
end
if isfield(options, 'axisOff')
    axisOff = options.axisOff ;
elseif isfield(options, 'axisOn')
    axisOff = ~options.axisOn ;
end

if colorTensors
    if isfield(options, 'climv')
        climv = options.climv ;
    else
        maxv = max(mags(:)) ;
        climv = [0,maxv] ;
    end
    if isfield(options, 'colors')
        colors = options.colors ;
    else
        colors = parula(100) ;
    end
    
end

if ~ishandle(fig)
    close all;
    fig = figure('Units', 'centimeters', 'Position', [0,0,4.5,4.5]) ;
end

% Figure out how large the array of crosses is, in #crossesW x #crossesH
nDV = round(size(QQ, 2) / qsubsample) ;
nAP = round(size(QQ, 1) / qsubsample) ;
if qsubsample ~= 1
    % Downsample the image entry-by-entry
    QQDS = zeros(nDV, nAP, 2,2) ;
    for pp = 1:2
        for qq = 1:2
            Qij = squeeze(QQ(:,:,pp, qq)) ;
            QQDS( :, :, pp, qq) = imresize(Qij, [nDV, nAP]) ;
        end
    end
else
    QQDS = QQ ;
end


% Fine tr/dev/theta if that's what we are plotting
if overlayHeatmap 
    if any(strcmpi({'dev','tr','theta', 'mag'}, plotStyle))
        heatmapQ = zeros(size(QQ, 1), size(QQ, 2)) ;
        for dv = 1:size(QQ, 1)
            if mod(dv, 200) == 0
                disp(['computing heatmapQ row=' num2str(dv) '/' num2str(size(QQ,1))])
            end
            for ap = 1:size(QQ, 2)
                mij = squeeze(QQ(dv, ap, :, :)) ;

                if any(strcmpi({'dev','tr','theta'}, plotStyle))
                    [tr, dev, theta] = ...
                        trace_deviator(mij, [1.0,0;0,1.0]) ;
                    if strcmpi('dev', plotStyle)
                        heatmapQ(dv,ap) = dev ;
                    elseif strcmpi('tr', plotStyle)
                        heatmapQ(dv,ap) = tr ;
                    elseif strcmpi('theta', plotStyle)
                        heatmapQ(dv,ap) = theta ;
                    end    
                else
                    heatmapQ(dv,ap) = sqrt(trace((mij * mij))) ;
                end
                % % Full diagonalization -- checking that these give the
                % % same results, also storing both eigs
                % [evec, evals] = eig( mij ) ;
                % [~, idx] = sort(diag(evals)) ;
                % evec = evec(:, idx) ;
                % pevec = evec(:, end) ;
                % theta2(ap, dv) = atan2(pevec(2), pevec(1)) ;
                % eval1(ap, dv) = evals(2, 2)  ;
                % eval2(ap, dv) = evals(1, 1)  ;            
            end
        end
    else
        error('Could not recognize plotStyle: dev, tr, theta, or mag')
    end
end

% trace deviator decomposition -- DownSampled
theta2 = zeros(nDV, nAP) ;
eval1 = zeros(nDV, nAP) ;
eval2 = zeros(nDV, nAP) ;
for dv = 1:nDV
    for ap = 1:nAP
        mij = squeeze(QQDS(dv, ap, :, :)) ;

        % Full diagonalization -- checking that these give the
        % same results, also storing both eigs
        [evec, evals] = eig( mij ) ;
        [~, idx] = sort(diag(evals)) ;
        evec = evec(:, idx) ;
        pevec = evec(:, end) ;
        theta2(dv, ap) = atan2(pevec(2), pevec(1)) ;
        eval1(dv, ap) = evals(2, 2)  ;
        eval2(dv, ap) = evals(1, 1)  ;            
    end
end

% Plot as crosses
% Note: x is the second dimension (horizontal) in the image
% and y is the first (vertical, going down). Matlab will be Matlab, alas...

if isempty(im)
    WW = size(QQ, 2) ;
    HH = size(QQ, 1) ;
else
    WW = size(im, 2) ;
    HH = size(im, 1) ;
end

if overlayHeatmap
    xx_fine = 1:WW ;
    yy_fine = 1:HH ;
end

xi = linspace(1, nAP, nAP) * WW / double(nAP) ;
yi = linspace(1, nDV, nDV) * HH / double(nDV) ;
[xx, yy] = meshgrid(xi, yi) ;
v1x = eval1 .* cos(theta2) * qscale ;
v1y = eval1 .* sin(theta2) * qscale ;
v2x = eval2 .* cos(theta2+pi/2) * qscale ;
v2y = eval2 .* sin(theta2+pi/2) * qscale ;

%% Plot it
if ~isempty(im)
    h1 = imshow(im) ;
    axis equal
    hold on;
else
    h1 = [] ;
end
if overlayHeatmap
    h2 = imagesc(xx_fine, yy_fine, heatmapQ) ;
    axis equal
    hold on;
else
    h2 = [] ;
end

if colorTensors
    q1 = quiverColorVectors2D(xx'-v1x, yy'-v1y,...
        2*v1x, 2*v1y, 0, colors, climv, lw) ;
    q1.ShowArrowHead = 'off';
    q2 = quiverColorVectors2D(xx'-v2x, yy'-v2y,...
        2*v2x, 2*v2y, 0, colors, climv, lw) ;
    q2.ShowArrowHead = 'off';
else
    q1 = quiver(xx-v1x, yy-v1y, 2*v1x, 2*v1y, 0, 'linewidth', lw, 'color', eig1color) ;
    q1.ShowArrowHead = 'off';
    q2 = quiver(xx-v2x, yy-v2y, 2*v2x, 2*v2y, 0, 'linewidth', lw, 'color', eig2color) ;
    q2.ShowArrowHead = 'off';
end

axis equal ;
if makeCbar
    cb = colorbar ;
    ylabel(cb, labelstr) ;
else
    cb = [] ;
end
if axisOff
    axis off
end
colormap gray

if ~isempty(clim)
    caxis(clim)
end

