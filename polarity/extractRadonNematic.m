function [angle, magnitude, results] = extractRadonNematic(im, options)
% EXTRACTRADONNEMATIC Convert an image patch to a nematic director
%   Find peaks in radon transform, compute nematic from each peak (cos2t)
%   and average vectors in 2*theta space. Then divide angle by 2 and return
%   the vector magnitude.
% 
% Todo: hmaximum transform (a morphological operation)
%
% Requires PeakFinding/ to be in current path
% 
% Parameters
% ----------
% d     The 2D data raw image - assumes a Double\Single-precision
%       floating-point, uint8 or unit16 array. Please note that the code
%       casts the raw image to uint16 if needed.  If the image dynamic range
%       is between 0 and 1, I multiplied to fit uint16. This might not be
%       optimal for generic use, so modify according to your needs.
% options   struct with optional fields: thres, filt, edg, res, method
%         thres     A number between 0 and 1 denoting the CDF below which to remove as background
%         filt      A filter matrix used to smooth the image. The filter size
%                   should correspond the characteristic size of the peaks
%         edg       A number>1 for skipping the first few and the last few 'edge' pixels
%         res       A handle that switches between two peak finding methods:
%                   1 - the local maxima method (default).
%                   2 - the weighted centroid sub-pixel resolution method.
%                   Note that res=2 method takes ~20% more time on average.
%         theta     array of angles onto which to project
%         savefn    full path of where to save figure summarizing result
%         show_fig  bool, show figure of the peaks detected
%         title     Title of the figure to save if savefn == true
%         Additional options for res == 1
%         peak_ar   aspect ratio in pix/Dela(deg) of neighborhood
%         peakw     width of area to sample around peak to measure weight 
%                   (ie robustness)
%         num2keep  int, number of peaks to use in finding orientation
%
% Returns
% -------
% angle : float in (0, pi)
%   The angle of the nematic alignment in half of S1
% magnitude : float
%   The magnitude of the nematic "vector" relating to the strength of the
%   peak relative to the rest of the image
% results : struct with fields
%       R : radon transform
%       xp : x values of radon
%       
% NPMitchell
 
%% Default options
num2keep = Inf ;
thres = 0.97 ;
show_fig = false ;
theta = -40:220 ;
filt = fspecial('gaussian', 7, 2) ;
edg = 2 ;
options.res = 1 ;
peakw = 6.5 ;
peak_ar = 3. ; % aspect ratio in pix/Dela(deg) of neighborhood

%% Find peaks in the transform

% thresholding on CDF to find peaks
if isfield(options, 'thres')
    thres = options.thres ;
end

if isfield(options, 'num2keep')
    num2keep = options.num2keep ;
end

% show the figure or not
if isfield(options, 'show_fig')
    show_fig = options.show_fig ;
end

% Angles to sample -- note must extend below 0 and above 180 
if isfield(options, 'theta')
    theta = options.theta ;
end

% filter to use on radon transform
if isfield(options, 'filt')
    filt = options.filt ; 
end

% width of boundary pixels to ignore
if isfield(options, 'edg')
    edg = options.edg ;
end

% Which method to use: local peak finding or weighted centroid
if ~isfield(options, 'res')
    options.res = 1 ;
end

% For method 1, additional parameters
% peak width for region around peak to measure average intensity
if isfield(options, 'peakw')
    peakw = options.peakw ;
end
if isfield(options, 'peak_ar')
    peak_ar = options.peak_ar ; 
end


% Take the Radon transform
[R,xp] = radon(im, theta);

% Check it
% fig = figure;
% imagesc(R)
% waitfor(fig)
% error('break')

if options.res == 1
    [cent d, ~, cm] = RadonPeakFind(R, thres, filt, edg, 1) ;
    xyc = [min(theta) + cent(1:2:end) - 1, min(xp) + cent(2:2:end) - 1] ;

    % Filter to (0, 180) deg
    keep = xyc(:, 1) < 180 & xyc(:, 1) > 0 ;
    xyc = xyc(keep, :) ;

    % Build measure of order from peaks
    [tg, pg] = meshgrid(theta, xp) ;
    weights = zeros(size(xyc, 1), 1) ;
    for i = 1:size(xyc, 1)
        % measure local intensity
        x2 = (tg - xyc(i, 1)) .^2 ;
        y2 = (peak_ar * (pg - xyc(i, 2))).^2 ;  % use aspect ratio
        near = (x2 + y2) < peakw ^2 ;
        % npix = length(nonzeros(near)) ;
        % Compare local brightness to average at this value of xp
        % neary = y2 < peakw ^ 2 ;
        % Rband = R .* neary ;
        % Rband_median = median(Rband(:)) ;
        peak_mean = sum(R .* near, 'all') / sum(near(:)) ;
        weights(i) = peak_mean / mean(R(:)) - 1 ;
        
        % check
        % fig = figure
        % imagesc(R .* near)
        % waitfor(fig)
    end

else
    %% Other option is hard threshold and centroids
    [cent, d, weights, cm] = RadonPeakFind(R, thres, filt, 2, 2) ;
    xyc = [min(theta) + cent(1:2:end) - 1, min(xp) + cent(2:2:end) - 1] ;

    % Filter to (0, 180) deg
    keep = xyc(:, 1) < 180 & xyc(:, 1) > 0 ;
    weights = weights(keep) ;
    xyc = xyc(keep, :) ;
end

% Filter by weight : ensure that this peak is significant in a global way
keep = weights > 0 ;
xyc = xyc(keep, :) ;
weights = weights(weights > 0) ;

% Keep k peaks
if num2keep > 1
    if isfinite(num2keep)
        [~, idx] = maxk(weights, num2keep) ;
        xyc = xyc(idx, :) ;
        weights = weights(idx, :) ;
    end
    
    % Take average vector to the various positions [cos(2t), sin(2t)]
    vecs = [cos(2 * xyc(:, 1) * pi / 180), sin(2 * xyc(:, 1) * pi / 180) ] ;
    % Weight vector magnitudes by weights
    vecw = [weights .* vecs(:, 1), weights .* vecs(:, 2)] ;
    vec = sum(vecw, 1) ;
    % Note: divide by two since we used cos(2t), sin(2t)
    angle = mod(atan2(vec(2), vec(1)) * 0.5, pi) ;
    magnitude = norm(vec) ;
else
    
    [magnitude, idx] = max(weights) ;
    angle = xyc(idx, 1) * pi / 180 ;
    
end

% Visualize the result if we don't return output
save_fig = isfield(options, 'savefn') ;
if save_fig
    save_fig = ~isempty(options.savefn)  ;
end
if show_fig || save_fig
    close all
    fig = figure ;
    cmap = colormap ;
    imagesc(theta,xp,d); hold on;
    tmp = max(1, uint8(length(colormap) * (weights - min(weights)) / (max(weights) - min(weights)))) ;
    scatter(xyc(:, 1), xyc(:, 2), 30, cmap(tmp, :), 'filled', 'markeredgecolor', 'k') ;
    plot(angle * 180 / pi * [1, 1], [min(xp), max(xp)], 'r--') ;
    % plot(xyc(weights < 0, 1), xyc(weights < 0, 2), 'rx')
    xlabel('\theta [deg]')
    ylabel('Projection axis')
    if isfield(options, 'title')
        title(options.title)
    else
        title('Filtered projected intensity')
    end
    cb = colorbar ;
    ylabel(cb, 'weights')
    caxis([min(weights), max(weights)])
    ylim([min(xp), max(xp)])
    if isfield(options, 'savefn') && options.savefn
        saveas(gcf, options.savefn)
    else
        waitfor(fig) 
    end
    
end

% Collate the results of the radon transform into a struct
results.R = R;
results.xp = xp;
results.cm = cm;
results.xyc = xyc ;
results.theta = theta ;

return 

