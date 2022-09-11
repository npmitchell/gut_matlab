%% Figure for radon demo
% figure_tubular_polarity_radon.m

im = imread(fullfile('/mnt/data/mef2GAL4klarUASCAAXmChHiFP/', ...
    '202003151700_1p4um_0p5ms3msexp/data/deconvolved_16bit/', ...
    'msls_output/gridCoords_nU0100_nV0100/', ...
    'PullbackImages_010step_uv/channel2/', ...
    'Time_000000_stab_pbuv.tif')) ;

im = imresize(im, [max(size(im)), max(size(im))]) ;

imshow(im)

[xe0, ye0] = meshgrid(200:20:400, 1050:20:1250) ;
xe = xe0(:) ;
ye = ye0(:) ;
angles = zeros(size(xe)) ;
mags = zeros(size(xe)) ;
results = cell(size(xe)) ;
for ideval = 1:length(xe)
    disp(['ii = ' num2str(ideval) '/' num2str(length(xe))])
    WW = 100 ;
    xcrop = xe(ideval) + [0, WW] ;
    ycrop = ye(ideval) + [0, WW] ;

    % Crop image to be a circle
    imc = im(ycrop(1):ycrop(2), xcrop(1):xcrop(2)) ;
    szx = diff(xcrop) ;
    szy = diff(ycrop) ;
    [xx, yy] = meshgrid(1:szx, 1:szy) ;
    mid = WW * 0.5 ;
    rr = sqrt((xx-mid).^2 + (yy-mid).^2) ;
    [rr, cc] = find(rr > WW*0.5) ;
    for ii = 1:length(rr)
        ri = rr(ii); ci = cc(ii) ;
        imc(ri, ci) = 0 ;
    end
    % imshow(imc)

    opts = struct('res', 1, 'show_fig', false, 'num2keep', 3) ;
    [angle, magnitude, result] = extractRadonNematic(imc, opts) ;

    angles(ideval) = angle ;
    mags(ideval) = magnitude ;
    results{ideval} = result ;
    % opts = struct('res', 2) ;
    % [angle2, magnitude2, results2] = extractRadonNematic(imc, opts) ;
end

%%
close all
angles = reshape(angles, size(xe0)) ;
mags = reshape(mags, size(xe0)) ;
opts = struct('colormap', phasemap, 'addQuiver', true, 'qsub', 2, 'makeCbar', false) ;
plotNematicField(mags, angles+ pi/2, opts)
view([0,90])
axis off
saveas(gcf, '/mnt/data/analysis/figure_tubular_polarity_radon_qsub2.png')
saveas(gcf, '/mnt/data/analysis/figure_tubular_polarity_radon_qsub2.pdf')

%% Raw image no quiver
theta = angles + pi/2 ;
cmap = phasemap ;
clim_mag = max(mags(:))
indx = max(1, round(mod(2*theta(:), 2*pi)*size(cmap, 1)/(2 * pi))) ;
colors = cmap(indx, :) ;
colors = min(mags(:) / clim_mag, 1) .* colors ;
% colors = reshape(colors, [size(angles,1), size(angles, 2), 3]) ;
ri = scatteredInterpolant(xe, ye, colors(:, 1), 'natural') ;
gi = scatteredInterpolant(xe, ye, colors(:, 2), 'natural') ;
bi = scatteredInterpolant(xe, ye, colors(:, 3), 'natural') ;
[xi, yi] = meshgrid(min(xe):max(xe), min(ye):max(ye)) ;
RR = ri(xi, yi) ;
GG = gi(xi, yi) ;
BB = bi(xi, yi) ;
imCC = reshape([RR,GG,BB], [size(xi, 1), size(xi, 2), 3]) ;
imagesc(imCC) ;
imwrite(imCC, '/mnt/data/analysis/figure_tubular_polarity_radon_Interp.png')


%% Corresponding image
imC = im(min(ye):max(ye)+WW, min(xe):max(xe)+WW) ;
imshow(imC)
imwrite(imC, '/mnt/data/analysis/figure_tubular_polarity_radon_Image.png')

%%
scatter(xe, ye,round(100*mags(:)), angles(:), 'filled')
colormap(isolum)

% 
imagesc(angles) ;
