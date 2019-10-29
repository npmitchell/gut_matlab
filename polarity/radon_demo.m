%% Demo for radon transform
close all
addpath('/mnt/data/code/gut_matlab/PeakFinding/')

% Prepare the outdir
outdir = './radon_demo/';
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

%% Unit square
mask = ones(40, 40) ;
[R,xp] = radon(mask, theta);
imagesc(theta,xp,R);
xlabel('\theta [degrees]')
ylabel('Projected intensity')
title('Radon transform of the unit square')
colorbar
saveas(gcf, fullfile(outdir, 'unit_square_radon.png'))
close all
imshow(mask) ;
print(fullfile(outdir, 'unit_square.png'))
close all

%% Unit disk
[xsz, ysz] = size(mask) ;
xcenter = xsz * 0.5 ;
ycenter = ysz * 0.5 ;
% Mask out a circle from the patch
% create a xygrid
[xc, yc] = meshgrid(1:xsz, 1:ysz) ;
dist = (xc - xcenter) .^2 + (yc - ycenter) .^2 ;
mask = dist' < w^2 ;
imagesc(mask)
axis equal
print(fullfile(outdir, 'unit_disk.png'))
close all
[R,xp] = radon(mask, theta);
imagesc(theta,xp,R);
xlabel('\theta [degrees]')
ylabel('Projected intensity')
title('Radon transform of the unit disk')
colorbar
saveas(gcf, fullfile(outdir, 'unit_disk_radon.png'))
close all

%% Test radon on one line to find PSF
theta = 0:180 ;
chunk = 0 * mask ;
for q = 1:length(chunk)
    if abs(q - 0.5 * length(chunk)) < 10
        chunk(q, q) = 1 ;
    end
end
imagesc(chunk)
axis equal
print(fullfile(outdir, 'lineseg.png'))
close all
[R,xp] = radon(chunk, theta);
imagesc(theta,xp,R);
xlabel('\theta [degrees]')
ylabel('Projected intensity')
title('Radon transform of linesegment')
colorbar
saveas(gcf, fullfile(outdir, 'lineseg_radon.png'))
close all
xsel = find(abs(xp) < 20);
ysel = 30:61 ;
psfi = R(xsel, ysel) ;
imagesc(xp(xsel), ysel, psfi) ;
hold on;
plot(mean(xp(xsel)), mean(ysel), 'o')
saveas(gcf, fullfile(outdir, 'lineseg_psf.png'))
close all

% Deconvolve
J = deconvreg(R, psfi) ;
imagesc(theta, xp, J) ;

xlabel('\theta [degrees]')
ylabel('Projected intensity')
title('Deconvolved Radon transform of linesegment')
saveas(gcf, fullfile(outdir, 'lineseg_deconvolved.png'))
close all

%% Test radon on two lines
chunk = 0 * mask ;
for q = 1:length(chunk)
    ri = mod(q - 5, length(chunk)) ;
    si = mod(q + 5, length(chunk)) ;
    if abs(ri - q) < 6 && ri > 0
        chunk(q, ri) = 1 ;
    end
    if abs(si - q) < 6 && si > 0
        chunk(q, si) = 1 ;
    end
end
imagesc(chunk)
axis equal
print(fullfile(outdir, 'twolines.png'))
close all
[R,xp] = radon(chunk, theta);
imagesc(theta,xp,R);
xlabel('\theta [degrees]')
ylabel('Projected intensity')
title('Radon transform of two lines')
colorbar
saveas(gcf, fullfile(outdir, 'twolines_radon.png'))
close all

%% Test radon on three lines
for q = 1:length(chunk)
    ri = mod(length(chunk) - q , length(chunk)) ;
    if ri > 0
        chunk(q, ri) = 1 ;
    end
end
imagesc(chunk)
axis equal
print(fullfile(outdir, 'threelines.png'))
close all
[R,xp] = radon(chunk, theta);
imagesc(theta,xp,R);
xlabel('\theta [degrees]')
ylabel('Projected intensity')
title('Radon transform of three lines')
colorbar
saveas(gcf, fullfile(outdir, 'threelines_radon.png'))
close all

% Build measure of nematic ordering
rvar = var(R, 1, 1) ;
plot(theta, rvar)
xlabel('\theta [deg]')
ylabel('variance')
title('Radon transform variance for 3 lines')
saveas(gcf, fullfile(outdir, 'threelines_variance.png'))

% Build measure of nematic ordering
rvar = sum(R.^2, 1) ;
plot(theta, rvar)
xlabel('\theta [deg]')
ylabel('\Sigma_i R(d_i, \theta)^2')
title('Radon transform squared for 3 lines')
saveas(gcf, fullfile(outdir, 'threelines_radonsquared.png'))

% Build measure of nematic ordering
rvar = sum(R.^3, 1) ;
plot(theta, rvar)
xlabel('\theta [deg]')
ylabel('\Sigma_i R(d_i, \theta)^3')
title('Radon transform cubed for 3 lines')
saveas(gcf, fullfile(outdir, 'threelines_radoncubed.png'))
close all

% Build measure of nematic ordering
c = mean(R(:)) + std(R(:)) ;
R = R - c;
R(R<0) = 0;
rvar = sum(R, 1) ;
plot(theta, rvar)
xlabel('\theta [deg]')
ylabel('\Sigma_{i for R(d_i,\theta) > c} R(d_i, \theta) ')
title('Radon transform thresholded for 3 lines')
saveas(gcf, fullfile(outdir, 'threelines_radonthreshold.png'))
close all

%% Five lines
% Test radon on five lines
chunk = 0 * mask ;
for q = 1:length(chunk)
    ri = mod(q - 5, length(chunk)) ;
    si = mod(q + 5, length(chunk)) ;
    if abs(ri - q) < 6 && ri > 0
        chunk(q, ri) = 1 ;
    end
    if abs(si - q) < 6 && si > 0
        chunk(q, si) = 1 ;
    end
    
    % Add other angle
    si = 2*q ;
    if si > 0 && si < length(chunk)
        chunk(q, si) = 2 ;
    end
    
    % Add backward diagonal
    ri = mod(length(chunk) - q , length(chunk)) ;
    if ri > 0
        chunk(q, ri) = 1 ;
    end
    
end
% image the chunk
imagesc(chunk)
axis equal
saveas(gcf, fullfile(outdir, 'fourlines.png'))
close all

% Radon and find peaks
[R,xp] = radon(chunk, theta);
imagesc(theta,xp,R);
xlabel('\theta [degrees]')
ylabel('Projected intensity')
title('Radon transform of five lines')
colorbar
saveas(gcf, fullfile(outdir, 'fourlines_radon.png'))
close all
filt = fspecial('gaussian', 7, 3) ;

% Deconvolve
R2=conv2(single(R),filt,'same') ;
J = deconvreg(R2, psfi) ;
imagesc(theta,xp,J);
xlabel('\theta [deg]')
ylabel('Projected intensity')
title('Deconvolved Radon transform of four lines')
colorbar
saveas(gcf, fullfile(outdir, 'fourlines_radondeconvolved.png'))
close all

cent = RadonPeakFind(R2, 0.97, filt, 2, 1) ;
xyc = [cent(1:2:end), min(xp) + cent(2:2:end)] ;
imagesc(theta,xp,R2);
hold on;
scatter(xyc(:, 1), xyc(:, 2), 40, 'k', 'filled') ;
xlabel('\theta [deg]')
ylabel('Projected intensity')
title('Radon transform of four lines')
colorbar
saveas(gcf, fullfile(outdir, 'fourlines_radonpeaks.png'))


%% Compute via automated process
options.savefn = fullfile(outdir, 'fourlines_radonpeaks_algorithm1.png') ;
options.res = 1 ;
options.title = 'Filtered Radon, res=1' ;
[xyc, weights, result] = extractRadonNematic(chunk, options) ;
weights

options.savefn = fullfile(outdir, 'fourlines_radonpeaks_algorithm2.png') ;
options.res = 2 ;
options.thres = 0.99 ;
options.title = 'Filtered Radon, res=2' ;
[xyc, weights, result] = extractRadonNematic(chunk, options) ;
weights


%% Patch of image (chunk)
chunk = [ 0     0     0     0     0     0     0   178   178   165   161   157   153     0     0     0     0     0     0     0 ; ...
     0     0     0     0     0   214   210   198   198   174   174   165   161   178   210     0     0     0     0     0 ; ...
     0     0     0   117   133   161   153   170   170   145   157   149   186   178   149   178   178     0     0     0 ; ...
     0     0   129   101   113   113   109    97    97   105   113   129   145   145   129   121   113   113     0     0 ; ...  
     0     0   133   133   129   129   113   113   125   125   125   137   137   137   121   113   101   101     0     0 ; ... 
     0   117   117   117   121   121   125   133   133   133   141   157   157   157   129   129   109   125   125     0 ; ...
     0   117   121   121   137   133   133   133   133   133   149   157   149   153   153   153   141   141   153     0 ; ...
   133   137   137   137   137   137   149   149   145   145   121   133   137   149   133   133   149   170   182   182 ; ...
   141   141   149   149   178   178   186   170   133   133   121   121   149   149   153   178   178   186   182   165 ; ...
   153   153   161   161   174   182   178   178   129   125   129   141   157   157   182   182   165   165   153   157 ; ...
   149   149   153   153   170   174   149   137   137   141   149   198   182   182   157   157   157   157   157   137 ; ...
   117   129   161   182   182   153   149   149   157   165   186   186   186   153   145   149   149   165   170   170 ; ...
   133   141   165   157   157   137   137   149   161   161   198   178   125   125   149   153   174   174   161   161 ; ...
     0   153   153   125   117   125   137   137   157   153   153   141   129   129   153   153   161   161   161     0 ; ...
     0   129   133   157   214   202   186   186   182   178   174   161   157   157   157   157   170   170   161     0 ; ...
     0     0   153   190   206   194   174   170   170   161   157   157   161   161   170   170   178   182     0     0 ; ...
     0     0   178   178   161   170   178   182   174   165   165   161   157   157   157   161   145   145     0     0 ; ...
     0     0     0   137   137   165   190   190   174   174   174   174   145   105   105   109   157     0     0     0 ; ...
     0     0     0     0     0   121   174   145   145   153   145   145   157   149   149     0     0     0     0     0 ; ... 
     0     0     0     0     0     0     0   141   157   149   149   174   145     0     0     0     0     0     0     0 ; ...
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 ]

% Check with image patch
% chunk = double(chunk.^2) / mean(double(chunk(:))) ;
[R,xp] = radon(chunk, theta);
imagesc(theta,xp,R);
xlabel('\theta [deg]')
ylabel('Projected intensity')
title('Radon transform of tissue')
colorbar
saveas(gcf, fullfile(outdir, 'tissue_radon.png'))
close all
imagesc(chunk)
axis equal
saveas(gcf, fullfile(outdir, 'tissue.png'))
close all

%% Automated scripts for tissue chunk
options.savefn = fullfile(outdir, 'tissue_radonpeaks_algorithm1.png') ;
options.res = 1 ;
options.title = 'Filtered Radon, res=1' ;
[xyc, weights, result] = extractRadonNematic(chunk, options) ;
weights

options.savefn = fullfile(outdir, 'tissue_radonpeaks_algorithm2.png') ;
options.res = 2 ;
options.thres = 0.97 ;
options.title = 'Filtered Radon, res=2' ;
[xyc, weights, result] = extractRadonNematic(chunk, options) ;
weights

%% Find peaks in the transform -- playing here
theta = -30:210 ;
[R,xp] = radon(chunk, theta);
imagesc(theta,xp,R);
% Inputs:
% d     The 2D data raw image - assumes a Double\Single-precision
%       floating-point, uint8 or unit16 array. Please note that the code
%       casts the raw image to uint16 if needed.  If the image dynamic range
%       is between 0 and 1, I multiplied to fit uint16. This might not be
%       optimal for generic use, so modify according to your needs.
% thres A number between 0 and max(raw_image(:)) to remove  background
% filt  A filter matrix used to smooth the image. The filter size
%       should correspond the characteristic size of the peaks
% edg   A number>1 for skipping the first few and the last few 'edge' pixels
% res   A handle that switches between two peak finding methods:
%       1 - the local maxima method (default).
%       2 - the weighted centroid sub-pixel resolution method.
%       Note that the latter method takes ~20% more time on average.
% fid   In case the user would like to save the peak positions to a file,
%       the code assumes a "fid = fopen([filename], 'w+');" line in the
%       script that uses this function.
filt = fspecial('gaussian', 7, 2) ;
cent = RadonPeakFind(R, 0.95, filt, 2, 1) ;
xyc = [min(theta) + cent(1:2:end), min(xp) + cent(2:2:end)] ;
imagesc(theta,xp,R);
hold on;
scatter(xyc(:, 1), xyc(:, 2), 40, 'k', 'filled') ;

% Filter to (0, 180) deg
keep = xyc(:, 1) < 180 & xyc(:, 1) > 0 ;
xyc = xyc(keep, :) ;

% Build measure of order from peaks
peakw = 6.5 ;
peak_ar = 3. ;  % aspect ratio in pix/Dela(deg) of neighborhood
[tg, pg] = meshgrid(theta, xp) ;
weights = zeros(size(xyc, 1), 1) ;
for i = 1:size(xyc, 1)
    % measure local intensity
    x2 = (tg - xyc(i, 1)) .^2 ;
    y2 = (peak_ar * (pg - xyc(i, 2))).^2 ;  % use aspect ratio
    near = (x2 + y2) < peakw ^2 ;
    npix = length(nonzeros(near)) ;
    % Compare local brightness to average at this value of xp
    neary = y2 < peakw ^ 2 ;
    Rband = R .* neary ;
    Rband_mean = sum(Rband(:)) / length(find(neary)) ;
    peak_mean = sum(R .* near, 'all') / sum(near(:)) ;
    weights(i) = peak_mean / Rband_mean - 1 ;
end
cmap = colormap ;
close all
imagesc(theta,xp,R); hold on;
tmp = max(1, uint8(length(colormap) * (weights - min(weights)) / (max(weights) - min(weights)))) ;
scatter(xyc(:, 1), xyc(:, 2), 30, cmap(tmp, :), 'filled', 'markeredgecolor', 'k') ;
plot(xyc(weights < 0, 1), xyc(weights < 0, 2), 'rx')
xlabel('\theta [deg]')
ylabel('Projected intensity')
title('Radon transform of tissue with peaks from gradients')
colorbar
saveas(gcf, fullfile(outdir, 'patch_radonpeakeval.png'))

%% Compare against hard threshold and centroids
[cent, cm, weights] = RadonPeakFind(R, 0.98, filt, 2, 2) ;
xyc = [min(theta) + cent(1:2:end), min(xp) + cent(2:2:end)] ;
imagesc(theta,xp,R);
hold on;
scatter(xyc(:, 1), xyc(:, 2), 40, 'k', 'filled') ;

% Filter to (0, 180) deg
keep = xyc(:, 1) < 180 & xyc(:, 1) > 0 ;
weights = weights(keep) ;
xyc = xyc(keep, :) ;

xlabel('\theta [deg]')
ylabel('Projected intensity')
title('Radon transform of tissue')
colorbar
saveas(gcf, fullfile(outdir, 'patch_radonpeaks.png'))

% Build measure of order from peaks
peakw = 6.5 ;
peak_ar = 3. ;  % aspect ratio in pix/Dela(deg) of neighborhood
[tg, pg] = meshgrid(theta, xp) ;
weight = zeros(size(xyc, 1), 1) ;
for i = 1:size(xyc, 1)
    % measure local intensity
    x2 = (tg - xyc(i, 1)) .^2 ;
    y2 = (peak_ar * (pg - xyc(i, 2))).^2 ;  % use aspect ratio
    near = (x2 + y2) < peakw ^2 ;
    npix = length(nonzeros(near)) ;
    % Compare local brightness to average at this value of xp
    neary = y2 < peakw ^ 2 ;
    Rband = R .* neary ;
    Rband_mean = sum(Rband(:)) / length(find(neary)) ;
    peak_mean = sum(R .* near, 'all') / sum(near(:)) ;
    weight(i) = peak_mean / Rband_mean - 1 ;
end
cmap = colormap ;
close all
imagesc(theta,xp,R); hold on;
tmp = max(1, uint8(length(colormap) * (weights - min(weights)) / (max(weights) - min(weights)))) ;
scatter(xyc(:, 1), xyc(:, 2), 30, cmap(tmp, :), 'filled', 'markeredgecolor', 'k') ;
plot(xyc(weights < 0, 1), xyc(weights < 0, 2), 'rx')
xlabel('\theta [deg]')
ylabel('Projected intensity')
title('Radon transform of tissue with peaks from gradients')
colorbar
saveas(gcf, fullfile(outdir, 'patch_radonpeakeval.png'))


% Take average vector to the various positions
vecs = [cos(xyc(:, 1) * pi / 180), sin(xyc(:, 1) * pi / 180) ] ;
% Weight vector magnitudes by weights
vecw = [weights .* vecs(:, 1), weights .* vecs(:, 2)] ;
vec = sum(vecw, 1) ;
angle = atan2(vec(2), vec(1)) ;
magnitude = norm(vec) ;
plot(angle * 180 / pi * [1, 1], [-15, 15], 'r--')
saveas(gcf, fullfile(outdir, 'patch_radonpeakevalavg.png'))

%% Inspect finite size of circles by tissue mask
mask = chunk > 0;
[R,xp] = radon(mask, theta);
imagesc(theta,xp,R);
xlabel('\theta [deg]')
ylabel('Projected intensity')
title('Radon transform of tissue')
colorbar
saveas(gcf, fullfile(outdir, 'mask_radon.png'))
close all
imagesc(mask)
axis equal
saveas(gcf, fullfile(outdir, 'mask.png'))
close all

%% 


