%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TURN COLLECTION OF TIFF STACKS INTO A MOSOAIC (STITCHING)
%
% convention : str
%   across, up
%
% see also
% --------
% stiching_tiffs.m
%
%
%%%% Isaac Breinyn 2019, NPMitchell 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

clear; clc ;
close all

% rootdir = '/mnt/data/RNAi/Control_mef2Gal4klar/' ;
% rootdir = fullfile(rootdir, '202004071914_control_mef2GAL4klar_hs60min32C_1730to1830_inHCOil27') ;
rootdir = cd ;
oper = 11/100; % overlap percentage (i.e. 10/100 --> 10%)
ncol = 7 ;
nrow = 6 ;
clipI = 0.9 ;

window = 200 ;
embryoDir = fullfile(rootdir, 'embryos') ;
snapDir = fullfile(embryoDir, 'snaps') ;
seqDir = fullfile(rootdir, 'imagesequence') ;

% Make dirs
dirs2make = {embryoDir, snapDir, seqDir} ;
for ii = 1:length(dirs2make)
    dir2make = dirs2make{ii} ;
    if ~exist(dir2make, 'dir')
        mkdir(dir2make)
    end
end

% Move directories & grab filenames
cd(rootdir) % Directory with image folders
addpath('/mnt/data/code/imsane_for_git/imsane/external/bfmatlab/')
addpath('/mnt/data/code/gut_matlab/tiff_handling/')
addpath('/mnt/data/code/gut_matlab/')
addpath('/mnt/data/code/gut_matlab/basics/')

%% Unpack into h5
fns = dir('*2017*.tif') ;
datfn = fullfile(rootdir, 'mosaic.h5') ;
if ~exist(datfn, 'file')
    if ncol > 1 || nrow > 0
        [data, frame1] = stitchTiffMosaic(fns, nrow, ncol, datfn, ...
            oper, 'uint8', clipI) ;
    else
        data = loadtiff(fullfile(fns(1).folder, fns(1).name)) ;
    end
else
    data = h5read(datfn, '/tileScan') ;
end
frame1 = squeeze(data(:, :, 1)) ;

%% Add aftermath images
% Unpack aftermath into h5
datfn = fullfile(rootdir, 'mosaic_aftermath.h5') ;
fns = dir('*aftermath*.tif') ;
if ~exist(datfn, 'file')
    [data2, frame2] = stitchTiffMosaic(fns, nrow, ncol, datfn, ...
        oper, 'uint8', clipI) ;
else
    data2 = h5read(datfn, '/tileScan') ;
end

%% Combine data + aftermath images
data = cat(3, data, data2) ;

%% Save data as image sequence
outfn = fullfile(seqDir, 'fullFOV_%06d.png') ;
for ii = 1:size(data, 3)
    disp(['saving image tp=' num2str(ii)])
    fn = sprintf(outfn, ii) ;
    imwrite(squeeze(data(:, :, ii)), fn) ;
end

%% Identify initial positions
ptfn = fullfile(rootdir, 'embryo_indices.mat') ;
if ~exist(ptfn, 'file')
    % Click on initial positions
    disp('Click each embryo, press Return when finished')
    pts = readPoints(frame1) ;

    % sort points by Y values
    [~, idx] = sort(pts(:, 2)) ;
    pts = pts(idx, :) ;

    % save points 
    save(ptfn, 'pts') ;
else
    load(ptfn, 'pts') ;
end
close all

%% Crop image for each embryo
fnbase = 'embryo%04d.tif' ;
for ii = 1:size(pts, 1)
    disp(num2str(ii))
    fn = fullfile(embryoDir, sprintf(fnbase, ii)) ;
    if ~exist(fn, 'file')
        pt = pts(ii, :) ;
        xmin = max(1, round(pt(1) - window)) ;
        xmax = min(size(data, 2), round(pt(1) + window)) ;
        ymin = max(1, round(pt(2) - window)) ;
        ymax = min(size(data, 2), round(pt(2) + window)) ;
        edat = data(ymin:ymax, xmin:xmax, :) ;

        % Save the data
        disp(['Saving stack to ' fn])
        for qq = 1:size(edat, 3)
            % write the first page as overwrite mode, then append
            if qq == 1
                imwrite(edat(:,:,qq), fn, 'tiff', 'Compression','none');
            else
                imwrite(edat(:,:,qq), fn, 'tiff', 'Compression','none','WriteMode','append');    
            end
        end
    end
end
disp('done')

%% Classification -- initialize fate log
classfn = fullfile(embryoDir, 'embryo_fate.txt') ;
fileID = fopen(classfn, 'a');
fprintf(fileID, ['EmbryoID, stageInit (hs init), stageFinal (hs end),' ...
    'fateID [0=WT,1=disrupted,2=death,3=indeterminate], missing folds, notes\n']);
fclose(fileID) ;

%% Classify each embryo
% Load initial and final conditions of heatshock
initfn = dir(fullfile(rootdir, '*stageInit.jpg')) ;
if length(initfn) > 0
    initim = imread(fullfile(initfn.folder, initfn.name)) ;
end
finalfn = dir(fullfile(rootdir, '*stageFinal.jpg')) ;
if length(finalfn) > 0
    finalim = imread(fullfile(finalfn.folder, finalfn.name)) ;
end

% rotate images (optional)
% initim = rot90(initim) ;
% finalim = rot90(finalim) ;


fileID = fopen(classfn, 'r');
nlines = linecount(fileID) ;
fclose(fileID) ;

% Go through each embryo to classify
for ii = nlines:size(pts, 1)
    pt = pts(ii, :) ;
    continue2class = false ;
    while ~continue2class
        close all    
        % Identify where we are in the FOV
        if exist('initim', 'var') && exist('finalim', 'var')
            subplot(3,3,[1 4 7])
            imshow(initim) ; 
            title('before heatshock') ;
            subplot(3,3,[2 5 8])
            imshow(finalim); 
            title('after heatshock') ;
            subplot(3,3, [3 6 9])
        end 
        imshow(frame1) ; hold on
        xmin = max(1, round(pt(1) - window)) ;
        xmax = min(size(data, 2), round(pt(1) + window)) ;
        ymin = max(1, round(pt(2) - window)) ;
        ymax = min(size(data, 2), round(pt(2) + window)) ;
        edat = data(ymin:ymax, xmin:xmax, :) ;
        xs = [xmin, xmax, xmax, xmin, xmin] ;
        ys = [ymin, ymin, ymax, ymax, ymin] ;
        plot(xs, ys, 'r-')
        plot(pt(1), pt(2), 'go')
        title(['ID of embryo' num2str(ii) '/' num2str(size(pts, 1))])
        disp(['EmbryoID=' num2str(ii) '/' num2str(size(pts, 1))])
        stageInit = input('stageInit = ', 's') ;
        if ~strcmp(stageInit, 'redo') && ~strcmp(stageInit, '')
            continue2class = true ;
        end
        stageFinal = input('stageFinal = ', 's') ;   
        if strcmp(stageFinal, 'redo') || strcmp(stageFinal, '')
            continue2class = false ;
        end
    end
    
    % Load each and toggle
    close all
    fn = fullfile(embryoDir, sprintf(fnbase, ii)) ;
    snapfn = fullfile(snapDir, sprintf(fnbase, ii)) ;
    edat = loadtiff(fn) ;
    titlestr = ['Classify embryo ' num2str(ii) '/' num2str(size(pts,1))] ;
    k = flipThroughStackFindLayer(edat, titlestr, 3, 10) ;
    snap = squeeze(edat(:, :, k)) ;
    imwrite(snap, snapfn) ;
    fateID = input('fateID [0=WT(default),1=disrupted,2=death,3=?] = ') ;
    if isempty(fateID)
        fateID = 0 ;
        disp('classified as WT')
    end
    if fateID > 0
        missingfolds = input('missing folds [0(default),1,2,3] = ') ;
        if isempty(missingfolds)
            missingfolds = 0 ;
            disp('no missing folds declared')
        end
    else
        missingfolds = 0 ;
    end
    notes = input('notes = ', 's') ;
    if isempty(notes)
        notes = 'none' ;
    end
    
    % clear figure handles
    clf
    close all
    
    % Save the resulting classification
    fileID = fopen(classfn, 'a');
    fprintf(fileID, '%d %s %s %d %d %s\n', ii, stageInit, ...
        stageFinal, fateID, missingfolds, notes);
    fclose(fileID) ;
    pause(0.001)
end

