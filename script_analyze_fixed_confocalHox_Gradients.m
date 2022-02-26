%% K99 figure for gradients in hox genes
addpath('/mnt/data/code/gut_matlab/addpath_recurse/')
addpath_recurse('/mnt/data/code/gut_matlab/')

%%
rootdir = '/mnt/data/antibodies/antp/UbxMutants/20210808200_UbxMutant1854_fixAntp8C11_1t50_ms488_1t500_1pc405flat_5pc488Flat/';
rootdir = fullfile(rootdir, 'e17_15a_RightDorsLatBand') ;
imfn = fullfile(rootdir, 'e17_15a_2p25x0p5um_RightDorsLatBand_C01_MAX') ;
im = imread([imfn '.tif']) ;
prob = h5read([imfn '_Probabilities.h5'], '/exported_data') ;
prob = squeeze(prob(1, :, :)) ;
prob = imgaussfilt(prob, 1) ;


%% 
I = uint8(prob*255)' ;
gmag = imgradient(I);
figure(2)
imagesc(im)
figure(1)
imagesc(gmag)
[locx, locy] = getpts() ;
xy = [locx, locy] ;
save([imfn '_nucleiLocations.mat'], 'xy')

%% Review nuclei
clf
imagesc(im) ;
for ptId = 1:length(locx)
    hold on; 
    plot(locx(ptId), locy(ptId), 'ro')
    title('Press enter to accept, backspace to correct pt')
    was_a_key = waitforbuttonpress;
    
    if was_a_key && strcmp(get(gcf, 'CurrentKey'), 'backspace')
        [newx, newy] = getpts() ;
        plot(locx(ptId), locy(ptId), 's')
        locx(ptId) = newx ;
        locy(ptId) = newy ;
    end
end

%%
radthres = 2 ;
[xx, yy] = meshgrid(1:size(im, 1), 1:size(im, 2)) ;
minima = zeros(size(im)) ;
for ptId = 1:length(locx)
    dist = (yy- locx(ptId)).^2 + (xx- locy(ptId)).^2 ;
    minima = minima | (dist' < radthres.^2) ;
end
gmag2 = imimposemin(gmag, minima);
LL = watershed(gmag2);
labels = imdilate(LL==0,ones(3,3)) ;
I4 = labeloverlay(I,labels);
imagesc(I4)
title('Markers and Object Boundaries Superimposed on Original Image')
Lrgb = label2rgb(LL,'jet','w','shuffle');
imshow(Lrgb)
title('Colored Watershed Label Matrix')
figure
imagesc(I)
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title('Colored Labels Superimposed Transparently on Original Image')
%% Remove largest or border
prop = regionprops(LL) ;
[maxArea, ind] = max([prop.Area]) ;
nuclei2measure = setdiff(1:max(LL(:)), [ind]) ;
imf = medfilt2(imgaussfilt(double(im), 1), [4,4]) ;

centroids = zeros(length(nuclei2measure), 2) ;
dmy = 1;
Imed = [] ; Imax = [] ; Imin = [] ; Imean = [] ; areas = [] ;
for qq = nuclei2measure
    centroids(dmy, :) = prop(qq).Centroid ;
    mask = (LL == qq) ;
    Imed(dmy) = median(imf(mask)) ;
    Imean(dmy) = mean(imf(mask)) ;
    Imax(dmy) = max(imf(mask)) ;
    Imin(dmy) = min(imf(mask)) ;
    areas(dmy) = prop(qq).Area ;
    dmy = dmy + 1 ;
end

%% Get point where inflection of nuclei is
close all
imagesc(im)
inflection = getpts() ;
flipx = true ;
pix2um = 0.12610821114 ;
if flipx
    positionAP = inflection - centroids(:, 1) ;
else
    positionAP = centroids(:, 1) - inflection ;
end
positionAP = positionAP * pix2um ;

%% Plot it
figure()
colors = define_colors ;
lineProps = {'-','color', colors(1, :)} ;
% h2=shadedErrorBar(positionAP, Imed, [Imax-Imed; Imed-Imin], 'lineProps', lineProps) ;
errorbar(positionAP, Imed, Imed-Imin, Imax-Imed, 'o') ;
xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
ylabel('Antp hox expression [a.u.]', 'interpreter', 'latex')
saveas(gcf, [imfn, '_nucleiAntpExpression.pdf'])
saveas(gcf, [imfn, '_nucleiAntpExpression.fig'])




%% K99 plot
close all
xedges = -35:4:35 ;
% binDataMeanStdWeighted(positionAP, Imax, xedges, areas)
[apbins, Ivals, Istd, ~, Iste] = binDataMeanStd(positionAP, Imax / mean(Imax(:)), xedges) ;
Istd2plot = movmax(Istd,4, 'omitnan') ;
% Istd2plot(Istd2plot==0) = NaN;
h2=shadedErrorBar(apbins, Ivals, [Istd2plot, Istd2plot]', 'lineProps', lineProps) ;
xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
ylabel('Antp hox expression [a.u.]', 'interpreter', 'latex')
saveas(gcf, [imfn, '_nucleiBinAntpExpression.pdf'])
saveas(gcf, [imfn, '_nucleiBinAntpExpression.fig'])


%% Save all data
save([imfn, '_nucleiMeasurements.mat'], 'positionAP', ...
    'centroids', 'Imed', 'Imax', 'Imin', 'Imean', 'LL', 'areas', 'pix2um', ...
    'apbins', 'Ivals', 'Istd', 'Iste')

% Save masked image
Lim = LL ;
Lim(LL == ind) = 0 ;
se = strel('disk', 10);
Lim = imdilate(Lim, se)
maskBlur = imgaussfilt(double(Lim > 0), 4) ;
outim = uint8(maskBlur .* medfilt2(imgaussfilt(double(im), 1), [4,4])) ;
imwrite(outim, [imfn, '_nucleiMaskBlur.png'])








%%

%%
% se = strel('disk',10);
% I = uint8(prob*255)' ;
% imshow(Io)
% title('Opening')
% Ie = imerode(im,se);
% Iobr = imreconstruct(Ie,I);
% imagesc(Iobr)
% title('Opening-by-Reconstruction')
% Ioc = imclose(Io,se);
% imagesc(Ioc)
% title('Opening-Closing')
% Iobrd = imdilate(Iobr,se);
% Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
% Iobrcbr = imcomplement(Iobrcbr);
% imagesc(Iobrcbr); axis equal
% title('Opening-Closing by Reconstruction')
% fgm = imregionalmax(Iobrcbr);
% imagesc(fgm)
% title('Regional Maxima of Opening-Closing by Reconstruction')
% I2 = labeloverlay(I,fgm);
% imagesc(I2)
% title('Regional Maxima Superimposed on Original Image')
% se2 = strel(ones(5,5));
% fgm2 = imclose(fgm,se2);
% fgm3 = imerode(fgm2,se2);
% fgm4 = bwareaopen(fgm3,20);
% I3 = labeloverlay(I,fgm4);
% imagesc(I3)
% title('Modified Regional Maxima Superimposed on Original Image')
% bw = imbinarize(Iobrcbr);
% imshow(bw)
% title('Thresholded Opening-Closing by Reconstruction')
% D = bwdist(bw);
% DL = watershed(D);
% bgm = DL == 0;
% imagesc(bgm)
% title('Watershed Ridge Lines')
% gmag = imgradient(I);
% gmag2 = imimposemin(gmag, bgm | fgm4);
% L = watershed(gmag2);
% labels = imdilate(L==0,ones(3,3)) + 2*bgm + 3*fgm4;
% I4 = labeloverlay(I,labels);
% imagesc(I4)
% title('Markers and Object Boundaries Superimposed on Original Image')
% Lrgb = label2rgb(L,'jet','w','shuffle');
% imshow(Lrgb)
% title('Colored Watershed Label Matrix')
% figure
% imagesc(I)
% hold on
% himage = imshow(Lrgb);
% himage.AlphaData = 0.3;
% title('Colored Labels Superimposed Transparently on Original Image')

%% 
% adaphisteqClip = 0 ;
% strelRadius = 0 ;
% seg = segmentImageWatershedSimple(prob, adaphisteqClip, strelRadius) ;

%%
% mask = prob > 0.5 ;
% bw = prob > 0.7 ;
% dd = bwdist(bw);
% ws = watershed(dd) ;
% ws(~mask) = 0 ;
% rgb = label2rgb(ws,'jet',[.5 .5 .5]);
% imshow(rgb)
% title('Watershed Transform')
% % imagesc(dd)
% xlim(xlims)
% ylim(ylims) ;