%% Make first GCaMP panel

rootdir = '/mnt/data/confocal_data/gut/Mef2GAL4klarUASGCAMP6sIII/' ;
rootdir = fullfile(rootdir, 'analysis_mef2G4kGCaMP6sIII') ;
resDirFn = 'anteriorFoldResults' ;

exptID = 2 ;
dThres = 10 ; % distance in um
nregions = 3 ;
spacing = 30 ;
posterFrameOffsets = [10 ] ;
lnstyle= '-' ;

load(fullfile(rootdir, resDirFn, 'Mef2G4kAnteriorFoldDorsalSettings.mat'), ...
    'clipY0s', 'dts', 'pix2um', 'foldXs', 'foldTs', 't0', ...
    'pcPower', 'dz', 'dates', 'expts', 'xfixed') ;

datdir = fullfile(rootdir, expts{exptID}) ;
fns12 = dir(fullfile(datdir, 'images', 'diffs_clipY02', '*d12*.png')) ;
fns23 = dir(fullfile(datdir, 'images', 'diffs_clipY02', '*d23*.png')) ;
fns13 = dir(fullfile(datdir, 'images', 'diffs_clipY02', '*d13*.png')) ;

% preallocate
dThres = dThres / pix2um(exptID) ;
spacing = spacing / pix2um(exptID) ;
reg00 = zeros(1, length(fns12)) ;
regp1 = zeros(1, length(fns12)) ;
regp2 = zeros(1, length(fns12)) ;
regn1 = zeros(1, length(fns12)) ;
regn2 = zeros(1, length(fns12)) ;
for ii = 1:length(fns12)
    % im1 = imread(fullfile(raw1(ii).folder, raw1(ii).name)) ;
    % im2 = imread(fullfile(raw2(ii).folder, raw2(ii).name)) ;
    % im3 = imread(fullfile(raw3(ii).folder, raw3(ii).name)) ;
    
    hatt12 = imread(fullfile(fns12(ii).folder, fns12(ii).name)) ;
    hatt13 = imread(fullfile(fns13(ii).folder, fns13(ii).name)) ;
    hatt23 = imread(fullfile(fns23(ii).folder, fns23(ii).name)) ;
    
    % Save transient signal in matrix
    sh12 = sum(hatt12, 1) ;
    sh23 = sum(hatt23, 1) ; 
    sh13 = sum(hatt13, 1) ;
    dd = mean([sh12; sh23; sh13], 1) ;

    % Background estimation
    % sii = min([sh12; sh23; sh13], [], 1) ;
    % dd = dd - sii ;
    
    % Collect activity in different regions
    xx = 1:length(dd) ;
    ind00 = find(abs(xx-foldXs(exptID)) < dThres ) ;
    beg00 = min(ind00) ;
    end00 = max(ind00) ;
    indp1 = find(abs(xx-foldXs(exptID) - spacing) < dThres ) ;
    begp1 = min(indp1) ;
    endp1 = max(indp1) ;
    indn1 = find(abs(xx-foldXs(exptID) +  spacing) < dThres ) ;
    begn1 = min(indn1) ;
    endn1 = max(indn1) ;
    
    reg00(ii) = mean_or_max(dd(beg00:end00)) ;
    regp1(ii) = mean_or_max(dd(begp1:endp1)) ;
    regn1(ii) = mean_or_max(dd(begn1:endn1)) ;
    
    if nregions > 3
        indp2 = find(abs(xx-foldXs(exptID) - ( spacing) * 2) < dThres/pix2um(exptID) ) ;
        begp2 = min(indp2) ;
        endp2 = max(indp2) ;
        regp2(ii) = mean_or_max(dd(begp2:endp2)) ;
        indn2 = find(abs(xx-foldXs(exptID) + (spacing) * 2) < dThres/pix2um(exptID) ) ;
        begn2 = min(indn2) ;
        endn2 = max(indn2) ;
        regn2(ii) = mean_or_max(dd(begn2:endn2)) ;
    end
end

close all
timestamps = dts(exptID) * ((1:length(fns12)) - t0(exptID)) ;
caxis([-1, 1])
colors = blueblackred(nregions) ;
hold on;
if nregions == 3
    plot(timestamps, regn1, lnstyle, 'color', colors(1, :)) ;
    plot(timestamps, reg00, lnstyle, 'color', colors(2, :)) ;
    plot(timestamps, regp1, lnstyle, 'color', colors(3, :)) ;
elseif nregions == 5
    plot(timestamps, regn2, lnstyle, 'color', colors(1, :)) ;
    plot(timestamps, regn1, lnstyle, 'color', colors(2, :)) ;
    plot(timestamps, reg00, lnstyle, 'color', colors(3, :)) ;
    plot(timestamps, regp1, lnstyle, 'color', colors(4, :)) ;
    plot(timestamps, regp2, lnstyle, 'color', colors(5, :)) ;
end
xlim([-60, 60])

% read poster frame
% imread(fullfile(rootdir, ))


function output = mean_or_max(input)
    output = max(input) ;
end

