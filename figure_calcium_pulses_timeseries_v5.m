%% Make first GCaMP panel

clear all
close all
clc

rootdir = '/mnt/data/confocal_data/gut/Mef2GAL4klarUASGCAMP6sIII/' ;
rootdir = fullfile(rootdir, 'analysis_mef2G4kGCaMP6sIII') ;
resDirFn = 'anteriorFoldResults' ;

% 4 8 or 12 for exptID=3, with d=8 and spacing=18 looks good for nregions=3
exptID = 3 ;
dThres = 8 ; % distance in um
nregions = 3 ;
spacing = 18 ;
posterFrameOffsets = [5,10,4,10,10 ] ;
lnstyle1= 's-' ;
lnstyle2= '^-' ;
lnstyle3= 'o-' ;
sps = 9/60 ;
 wframe = 5 ;
 forceXoff = 0 ;

load(fullfile(rootdir, resDirFn, 'Mef2G4kAnteriorFoldDorsalSettings.mat'), ...
    'clipY0s', 'dts', 'pix2um', 'foldXs', 'foldTs', 't0', ...
    'pcPower', 'dz', 'dates', 'expts', 'xfixed') ;

datdir = fullfile(rootdir, expts{exptID}) ;

% fns12 = dir(fullfile(datdir, 'images', 'diffs_clipY02', '*d12*.png')) ;
% fns23 = dir(fullfile(datdir, 'images', 'diffs_clipY02', '*d23*.png')) ;
% fns13 = dir(fullfile(datdir, 'images', 'diffs_clipY02', '*d13*.png')) ;

raw1 = dir(fullfile(datdir, 'images', '*c001.png')) ;
raw2 = dir(fullfile(datdir, 'images', '*c002.png')) ;
raw3 = dir(fullfile(datdir, 'images', '*c003.png')) ;

% preallocate
fns = raw1; 
dThres = dThres / pix2um(exptID) ;
spacing = spacing / pix2um(exptID) ;
reg00a = zeros(1, length(fns)) ;
regp1a = zeros(1, length(fns)) ;
regp2a = zeros(1, length(fns)) ;
regn1a = zeros(1, length(fns)) ;
regn2a = zeros(1, length(fns)) ;
reg00b = zeros(1, length(fns)) ;
regp1b = zeros(1, length(fns)) ;
regp2b = zeros(1, length(fns)) ;
regn1b = zeros(1, length(fns)) ;
regn2b = zeros(1, length(fns)) ;
reg00c = zeros(1, length(fns)) ;
regp1c = zeros(1, length(fns)) ;
regp2c = zeros(1, length(fns)) ;
regn1c = zeros(1, length(fns)) ;
regn2c = zeros(1, length(fns)) ;
for ii = 1:length(fns)
    im1 = imread(fullfile(raw1(ii).folder, raw1(ii).name)) ;
    im2 = imread(fullfile(raw2(ii).folder, raw2(ii).name)) ;
    im3 = imread(fullfile(raw3(ii).folder, raw3(ii).name)) ;
    
    % hatt12 = imread(fullfile(fns12(ii).folder, fns12(ii).name)) ;
    % hatt13 = imread(fullfile(fns13(ii).folder, fns13(ii).name)) ;
    % hatt23 = imread(fullfile(fns23(ii).folder, fns23(ii).name)) ;
    
    % Save transient signal in matrix
    sa = sum(im1, 1) ;
    sb = sum(im2, 1) ; 
    sc = sum(im3, 1) ;
    dd = mean([sa; sb; sc], 1) ;

    % Background estimation
    % sii = min([sh12; sh23; sh13], [], 1) ;
    % dd = dd - sii ;
    
    % Collect activity in different regions
    if ii == 1
        xx = 1:length(dd) ;
        ind00 = find(abs(xx-foldXs(exptID) - forceXoff) < dThres ) ;
        beg00 = min(ind00) ;
        end00 = max(ind00) ;
        indp1 = find(abs(xx-foldXs(exptID) - spacing- forceXoff) < dThres ) ;
        begp1 = min(indp1) ;
        endp1 = max(indp1) ;
        indn1 = find(abs(xx-foldXs(exptID) +  spacing- forceXoff) < dThres ) ;
        begn1 = min(indn1) ;
        endn1 = max(indn1) ;
    end
    
    reg00a(ii) = mean_or_max(sa(beg00:end00)) ;
    reg00b(ii) = mean_or_max(sb(beg00:end00)) ;
    reg00c(ii) = mean_or_max(sc(beg00:end00)) ;
    regp1a(ii) = mean_or_max(sa(begp1:endp1)) ;
    regp1b(ii) = mean_or_max(sb(begp1:endp1)) ;
    regp1c(ii) = mean_or_max(sc(begp1:endp1)) ;
    regn1a(ii) = mean_or_max(sa(begn1:endn1)) ;
    regn1b(ii) = mean_or_max(sb(begn1:endn1)) ;
    regn1c(ii) = mean_or_max(sc(begn1:endn1)) ;
    
    if nregions > 3
        indp2 = find(abs(xx-foldXs(exptID) - ( spacing) * 2) < dThres/pix2um(exptID) ) ;
        begp2 = min(indp2) ;
        endp2 = max(indp2) ;
        regp2(ii) = mean_or_max(sa(begp2:endp2)) ;
        indn2 = find(abs(xx-foldXs(exptID) + (spacing) * 2) < dThres/pix2um(exptID) ) ;
        begn2 = min(indn2) ;
        endn2 = max(indn2) ;
        regn2(ii) = mean_or_max(sa(begn2:endn2)) ;
    end
end

timestamps = dts(exptID) * (1:length(raw1)) - t0(exptID) ;

% Get colors
caxis([-1, 1])
colors = blueblackred(nregions) ;
colors = define_colors(3) ;

r00 = [reg00a; reg00b; reg00c]; %; NaN*reg00a] ;
rn1 = [regp1a; regp1b; regp1c]; %; NaN*reg00a] ;
rp1 = [regn1a; regn1b; regn1c] ; %; NaN*reg00a] ;
timestamps = [timestamps; timestamps + sps; timestamps + 2*sps]; %
% If you wish to disconnect each triplet, add the following line to the
% previous:
%; NaN*reg00a] ;

r00 = r00(:) ;
rn1 = rn1(:) ;
rp1 = rp1(:) ;
timestamps = timestamps(:) ;
% med00 = movmedian(r00, 20) ;
% medp1 = movmedian(rp1, 20) ;
% medn1 = movmedian(rn1, 20) ;
med00 = movmin(r00, 15) ;
medp1 = movmin(rp1, 15) ;
medn1 = movmin(rn1, 15) ;

m00 = movmax(r00- med00, 10) ;
mp1 = movmax(rp1- medp1, 10) ;
mn1 = movmax(rn1- medn1, 10) ;

close all
fig = figure('units','centimeters','position',[0,0,10,4]);
subplot(1, 2, 2)
hold on;
if nregions == 3
    % plot(timestamps, rn1 - medn1, lnstyle2, 'color', colors(2, :)) ;
    % plot(timestamps, rp1 - medp1, lnstyle3, 'color', colors(3, :)) ;
    plot(timestamps, r00 - med00, lnstyle1, 'color', colors(1, :)) ;
    %plot(timestamps, r00 - med00, '.-', 'color', colors(1, :)) ;
    %plot(timestamps, m00, '--', 'color', colors(1, :)) ;
    %plot(timestamps, mn1, '--', 'color', colors(2, :)) ;
    %plot(timestamps, mp1, '--', 'color', colors(3, :)) ;
    
elseif nregions == 5
    plot(timestamps, rn2, lnstyle1, 'color', colors(1, :)) ;
    plot(timestamps, rn1, lnstyle1, 'color', colors(2, :)) ;
    plot(timestamps, reg00, lnstyle1, 'color', colors(3, :)) ;
    plot(timestamps, regp1, lnstyle1, 'color', colors(4, :)) ;
    plot(timestamps, regp2, lnstyle1, 'color', colors(5, :)) ;
end
xlim([-15, 30])
yticklabels([])
ylims = ylim ;

xlabel('time [min]', 'interpreter', 'latex')
ylabel('fluctuating GCaMP intensity, $\delta I$ [a.u.]', 'interpreter', 'latex')
title('Max of dv-integrated signal within each region')

% Save what was plotted
outrootdir = '/mnt/data/analysis/eLife_gut2022_minimalData/';
fn = fullfile(outrootdir, 'figure5_calcium_example_embryo_dorsalAnterior/timetrace_anteriorConstriction.mat') ;
save(fn, 'r00', 'med00', 'timestamps')

%% read poster frame
% for posterFrameOffsets = -1:25 ;
poster = t0(exptID) / dts(exptID) + posterFrameOffsets(exptID) ;
im1 = imread(fullfile(raw1(poster).folder, raw1(poster).name)) ;
im2 = imread(fullfile(raw2(poster).folder, raw2(poster).name)) ;
im3 = imread(fullfile(raw3(poster).folder, raw3(poster).name)) ;


clipI = 50 ;
imR = mat2gray(im1, [0, clipI]) ;
imG = mat2gray(im2, [0, clipI]) ;
imB = mat2gray(im3, [0, clipI]) ;
rgb = cat(3, imR, imG, imB) ;
[szy, szx, szz] = size(rgb) ;
xoffbar = 170 ;
yoffbar = 50 ;
wbar = 15 ;
% Mark scalebar
rgb(yoffbar:yoffbar+wbar, xoffbar:xoffbar + 25/pix2um(exptID), :) = 1 ;

% Mark regions
rgb(szy-wframe:szy, beg00:end00, 1) = colors(1, 1) ;
rgb(szy-wframe:szy, beg00:end00, 2) = colors(1, 2) ;
rgb(szy-wframe:szy, beg00:end00, 3) = colors(1, 3) ;
rgb(szy-wframe:szy, begn1:endn1, 1) = colors(2, 1) ;
rgb(szy-wframe:szy, begn1:endn1, 2) = colors(2, 2) ;
rgb(szy-wframe:szy, begn1:endn1, 3) = colors(2, 3) ;
rgb(szy-wframe:szy, begp1:endp1, 1) = colors(3, 1) ;
rgb(szy-wframe:szy, begp1:endp1, 2) = colors(3, 2) ;
rgb(szy-wframe:szy, begp1:endp1, 3) = colors(3, 3) ;
subplot(1, 2, 1)
imshow(rgb)
% title(num2str(posterFrameOffsets))
pause(1)
%end
exten = sprintf('_expt%d', exptID) ; 
outFigFn = fullfile(rootdir, resDirFn, ['posterFigure' exten '.pdf']) ;
if exist(outFigFn, 'file')
    exten = sprintf('%d_v2', exptID) ; 
end
saveas(gcf, fullfile(rootdir, resDirFn, ['posterFigure' exten '_v5.pdf'])) ;
imwrite(rgb, fullfile(rootdir, resDirFn, ['posterImage' exten '_v5.png'])) ;


%% Second panel -- used to filter the curves, now we sum them
close all
fig = figure('units','centimeters','position',[0,0,10,4]) ;
subplot(1, 2, 2)
cla
hold on;

% % filter the curves
% medianKernel = 10 ; 
% meanKernel = 25 ;
% minKernel = 15;
% maxKernel = 15 ;
% filtn1 = movmedian(rn1 - medn1, medianKernel) ;
% filtp1 = movmedian(rp1 - medp1, medianKernel) ;
% filt00 = movmedian(r00 - med00, medianKernel) ;
% filtn1 = movmean(filtn1, meanKernel) ;
% filtp1 = movmean(filtp1, meanKernel) ;
% filt00 = movmean(filt00, meanKernel) ;
% lowern1 = movmin(rn1 - medn1, minKernel) ;
% lowerp1 = movmin(rp1 - medp1, minKernel) ;
% lower00 = movmin(r00 - med00, minKernel) ;
% uppern1 = movmean(rn1 - medn1, maxKernel) ;
% upperp1 = movmean(rp1 - medp1, maxKernel) ;
% upper00 = movmean(r00 - med00, maxKernel) ;
% m00 = movmax(filt00, 1) ;
% mp1 = movmax(filtp1, 1) ;
% mn1 = movmax(filtn1, 1) ;

% % Sum the curves
m00 = cumsum(r00 - med00) ;
mp1 = cumsum(rp1 - medp1) ;
mn1 = cumsum(rn1 - medn1) ;

plot(timestamps, m00, '-', 'color', colors(1, :)) ;
plot(timestamps, mn1, '-', 'color', colors(2, :)) ;
plot(timestamps, mp1, '-', 'color', colors(3, :)) ;

xlim([-15, 30])
yticklabels([])
ylimits = [min([min(m00(timestamps >= -15)), min(mn1(timestamps >= -15)), min(mp1(timestamps >= -15))]), ...
    max([max(m00(timestamps <=30)), max(mn1(timestamps <= 30)), max(mp1(timestamps <=30))])]; 
xlabel('time [min]', 'interpreter', 'latex')
ylim(ylimits) ;
ylabel('integrated fluctuating GCaMP intensity, $\int_0^t \delta I \, dt$ [a.u.]', 'interpreter', 'latex')
saveas(gcf, fullfile(rootdir, resDirFn, ['posterFigure' exten '_part2_v5.pdf'])) ;


%% Save stuff
fluctuatingData = r00 - med00 ;
fluctuatingAhead = rp1 - medp1 ;
fluctuatingBehind = rn1 - medn1 ;
save(fullfile(rootdir, resDirFn, 'figureData.mat'), ...
    'm00', 'mp1', 'mn1', 'fluctuatingData', ...
    'fluctuatingAhead', 'fluctuatingBehind')
%% Function definitions
function output = mean_or_max(input)
    output = max(input) ;
end