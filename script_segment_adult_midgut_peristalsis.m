%% Segment the cells during peristalsis of the gut

clear all
close all
rootdir = '/Users/npmitchell/Dropbox/Soft_Matter/PAPER/tubular_paper/gut_peristalsis/';
datdir = fullfile(rootdir, 'image_sequence/raw_data/') ;
fn = fullfile(datdir, 'gut_peristalsis_motility_%04d.png') ;
filebase = fullfile(datdir, 'gut_peristalsis_motility_%04d_Probabilities.h5') ;
segImFn = fullfile(rootdir, 'seg2d', 'seg_%04d.png') ;
trackImFn = fullfile(rootdir, 'track2d', 'track', 'track_%04d.png') ;
segFn = fullfile(rootdir, 'seg2d', 'seg_%04d.mat') ;
trackFn = fullfile(rootdir, 'track2d', 'track_%04d.mat') ;
outdir = rootdir ;

trackdir = fullfile(rootdir, 'track2d') ;
if ~exist(trackdir, 'dir')
    mkdir(trackdir)
    mkdir(fullfile(trackdir, 'track')) ;
end

if ~exist(fullfile(outdir, 'seg2d'), 'dir')
    mkdir(fullfile(outdir, 'seg2d'))
end

%%
addpath(genpath('/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/gut_matlab'))
addpath(genpath('/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/tissueAnalysisSuite'))

%%
iLastikVersion = 0  ;
cellSize = 100 ;
strelRadius = -5 ; % negative values close cells
gaussKernel = 1 ;
heightMinimum = 3.5 ;
% how far in pixels is too far to connect two cell vertices
very_far = 550 ;
overwrite = true ;
overwriteImages = true ;
timeUnits = 's' ;
areaOpen = 5 ;
skipPolygons = true ;

%% Run through all timepoints

timePoints = 0:29 ;

for tp = timePoints
    
    outfn = sprintf(segFn, tp) ;
    if ~exist(outfn, 'file') || overwrite

        % Define path to this timePoint's hdf5 probabilities file
        h5fn = sprintf(filebase, tp) ;
        [ mem ] = load.ilastikh5Single( h5fn, iLastikVersion );

        %% Segment the membrane.
        % cellSize : int (default=200)
        %   kernel size for laplacian of Gaussian : set to scale of curvature
        %   picking out, around a cell size or higher, in units of area (pix^2)
        % strelRadius : int (default=1)
        %   strel disk radius for dilation of segmented image
        % gaussKernel : float (default=2)
        %   kernel size for Gaussian filter. Set to a couple pixels. 
        %   Has units of length (pix).
        % heightMinimum : float (default=3.5)
        %   height of any local minima to merge, to reduce noise at rugged
        %   minima
        L = seg.memWS(mem, cellSize, strelRadius, gaussKernel, heightMinimum, 200);
        
        % imagesc(L)
        % waitfor(gcf)
        
        rp = regionprops(L, 'centroid', 'area', 'perimeter') ;
                
        % Kill the largest region (background)
        [~, maxID] = max([rp.Area]) ;
        L(L==maxID) = 0 ;
        L(L>maxID) = L(L>maxID) - 1 ;
        
        % Remove all but the desired cell to track
        if tp == timePoints(1)
            
            rp = regionprops(L, 'centroid', 'area', 'perimeter') ;
            imagesc(L); axis equal; colorbar
            N2keep = 2 ;
            % select the cells to keep
            title('click on cells to keep')
            xy = ginput(N2keep) ;
            
            % collate centroids of all regions/cells
            xyc = zeros(length(rp), 2) ;
            for cid = 1:length(rp)
                xyc(cid, :) = rp(cid).Centroid ;
            end
            
            % find the IDs of the cells to keep
            keep = [] ;
            for tmp = 1:N2keep
                [dist, minID] = min(sqrt(sum((xyc- xy(tmp, :)).^2, 2))) ;
                keep = [keep, minID] ;
            end
            
            % toss all other cells
            reg2rm = setdiff(1:length(rp), keep) ;
            for rr = reg2rm
                L(L == rr) = 0 ;
            end
            
            % relabel keeping only the cells to keep
            for tmp = 1:N2keep
                L(L == keep(tmp)) = tmp ;
            end
            
            imshow(L)
            pause(1)
        end
        rp = regionprops(L, 'centroid', 'area', 'perimeter') ;
        
        % Now link to previous timepoint
        xys = zeros(length(rp), 5) ;
        xys(:, 6) = tp ;
        xys(:, 5) = 1:length(rp) ;
        for cid = 1:length(rp)
            xys(cid, 1:2) = rp(cid).Centroid ;
            xys(cid, 3) = rp(cid).Area ;
            xys(cid, 4) = rp(cid).Perimeter ;
        end
        if tp > timePoints(1)
            xyap = [xyap; xys] ;
             
            % % Tracking
            % param = struct() ;
            % param.dim = 2 ;
            % xyapSorted = track(xyap, 100, param) ;
            % 
            % % Use sorting to relabel
            % % map is from column 5 to column 7
            % tId = xyapSorted(:, 6) ;
            % map1L = xyapSorted(tId==tp, [5,7]) ;
            % L2 = L ;
            % for cidx = 1:size(map1L, 1)
            %     c2replace = map1L(cidx, 1) ;
            %     L2(L == c2replace) = map1L(cidx, 2) ;
            % end
            % L = L2 ;
            
        else
            xyap = xys ;
        end    
        
        [seg2d, segIm] = constructSegAndSave(L, very_far, outfn, skipPolygons) ;
        
        
    else
        disp(['already on disk: ' outfn])
        
        % %% Convert to simpler format
        load(outfn, 'seg2d', 'segIm')
        
    end
    
    close all
    outim = saveSegImage(segIm, segImFn, fn, tp, overwrite || overwriteImages) ;
    imshow(outim) 
    set(gcf, 'visible', 'on')
    pause(1e-6)
end

% Save xyap
disp('saving xyap to disk...')
readme = 'x,y,area,perimeter, timepoint' ;
xyapfn = fullfile(rootdir, 'xyapt.mat') ;
save(xyapfn, 'xyap', 'readme')

clearvars outfn seg2d segIm tp L cid xys mem rp

%% Manual tracking by clicking
% [tracks, trackGraph] = manualTrack2D({}, fileBase, timePoints, trackOutfn, tracks2Add, tidx0)


%% Get scale
pix2um = 20.0 / 87.0 ; % um per pixel

%% Go back and redo earlier timepoints images with tracked results

% Tracking
param = struct() ;
param.dim = 2 ;
xyapSorted = trackNearest(xyap, 60, param) ;
disp('saving xyapSorted to disk...')
readme = 'x, y, area, perimeter, timepoint, trackID' ;
xyapfn = fullfile(rootdir, 'xyapSorted.mat') ;
save(xyapfn, 'xyapSorted', 'readme')

for tidx = 1:length(timePoints)
    tp = timePoints(tidx) ;
    
    infn = sprintf(segFn, tp) ;
    outfn = sprintf(trackFn, tp) ;
    load(infn, 'segIm')   
    tId = xyapSorted(:, 6) ;
    mapL = xyapSorted(tId==tp, [5,7]) ;
    mapL(:, 1) = mapL(:, 1) + 1 ;
    
    % Use sorting to relabel
    L2 = 0*segIm ;
    for cidx = 1:size(mapL, 1)
        c2replace = mapL(cidx, 1) ;
        L2(segIm == c2replace) = mapL(cidx, 2) ;
    end       
    % [seg2d, segIm] = constructSegAndSave(L2, very_far, outfn, skipPolygons) ; 
    save(outfn, 'L2')    
    
    
    % Areas of tracked regions
    ucells = unique(xyapSorted(:, end)) ;
    nCells = length(ucells) ;
    areas = cell(nCells, 1) ;
    peris = cell(nCells, 1) ;
    shapes = cell(nCells, 1) ;
    for cidx = 1:nCells
        cid = ucells(cidx) ;
        areas{cid} = xyapSorted(xyapSorted(:, end) == cid, 3) ;
        areas{cid} = areas{cid} * pix2um.^2 ;
        peris{cid} = xyapSorted(xyapSorted(:, end) == cid, 4) ;
        peris{cid} = peris{cid} * pix2um ;
        shapes{cid} = sqrt(areas{cid}) ./ peris{cid} ;
    end
    
    % Plot area over time
    if nCells < 3
        cmpQ = distinguishable_colors(nCells, [0,0,0;1,1,1;0,0,1;1,0,0;0,1,0;1,0.10345,0.72414]);
    else
        cmpQ = isolum(nCells+1) ; 
        cmpQ = cmpQ(1:nCells, :) ;
    end
    
    outim = saveSegImage(L2, trackImFn, fn, tp, true, cmpQ) ; 
    % for some reason we flip the color order here to match the image --
    % not sure why!!!
    cmpQ = flipud(0.95 * cmpQ ./ vecnorm(cmpQ, 2, 2)) ;
    
    clf
    % area
    subplot(2, 2, 1)
    for cidx = 1:nCells
        cid = ucells(cidx) ;
        plot(timePoints, areas{cid}, 'color', cmpQ(cid, :))
        hold on;
        plot(tp, areas{cid}(tidx), 'o', 'color', cmpQ(cid, :))
    end
    xlabel('time [s]', 'interpreter', 'latex')
    ylabel('area, $A$ [$\mu$m$^2$]', 'interpreter', 'latex')
    ylim([0, Inf])
    
    % perimeter
    subplot(2, 2, 2)
    for cidx = 1:nCells
        cid = ucells(cidx) ;
        plot(timePoints, peris{cid}, 'color', cmpQ(cid, :))
        hold on;
        plot(tp, peris{cid}(tidx), 'o', 'color', cmpQ(cid, :))
    end 
    xlabel('time [s]', 'interpreter', 'latex')
    ylabel('perimeter, $P$ [$\mu$m]', 'interpreter', 'latex')
    ylim([0, Inf])
    
    % shape index
    subplot(2, 2, 3)
    for cidx = 1:nCells
        cid = ucells(cidx) ;
        plot(timePoints, shapes{cid}, 'color', cmpQ(cid, :))
        hold on;
        plot(tp, shapes{cid}(tidx), 'o', 'color', cmpQ(cid, :))    
    end
    xlabel('time [s]', 'interpreter', 'latex')
    ylabel('shape index, $P/\sqrt{A}$', 'interpreter', 'latex')
    ylim([0, Inf])
        
    % cell identifier in image
    subplot(2, 2, 4)
    imshow(outim)
    hold on; 
    xy = xyapSorted(tId==tp, 1:2) ;    
    % cmpQ = distinguishable_colors(max(segIm(:)), [0,0,0;1,1,1;0,0,1;1,0,0;0,1,0;1,0.10345,0.72414]);
    scatter(xy(:, 1), xy(:, 2), 50, 1:length(xy(:, 1)), 'filled', 'markerfacecolor', 'none')
    sgtitle(['$t=' sprintf('%03d$', tp) ' ' timeUnits], 'interpreter', 'latex')
    pause(0.0001)
    saveas(gcf, fullfile(trackdir, 'metrics', sprintf('segL_%03d.png', tp)))
end


%% AUX FUNCTIONS
function [seg2d, segIm] = constructSegAndSave(L, very_far, outfn, skipPolygons)

        % Set bond=0 and clear_border = 1
        [L, Struct] = seg.generate_structs(L, 0, 1, 0, very_far);
        % Bad cells are bubble cells, which is a segmentation that forked and
        % reconnected.
        % ToDo: Should we do this? Can we skip it or does that lead to issues?
        % disp('removing bad cells')
        % L = seg.removeBadCells(Struct, L);
        % Now change label matrix after removing bad cells
        % disp('relabelling cells')
        % L = seg.relabelL(L);
        % Now also synchronize Struct after removing bad cells
        [segIm, seg2d] = seg.generate_structs(L, 0, 0, 0, very_far);
        disp('done with segmentation')

        %% Prepare data structure for inverse (optional? Does this improve segmentation?)
        % % put a parameter in the cdat of Struct, a boolean of whether every vertex
        % % is 3-fold.
        % Struct = seg.threefold_cell(Struct);
        % % generate the bdat structure in Struct
        % Struct = seg.recordBonds(Struct, L);
        % disp('generated the bond structure')
        % % Segment the curvature of each bond
        % Struct = seg.curvature(Struct, size(L));
        % disp('segmented the curvature of each bond')
        % % Remove all fourfold vertices, recursively if there are z>4
        % Struct = seg.removeFourFold(Struct);
        % disp('removed fourfold vertices')
        % % The inverse is ill-posed if we have convex cells, so hack those to be
        % % convex
        % Struct = seg.makeConvexArray(Struct);
        % disp('done with data preparation')

        %% Convert to simpler format
        disp('Constructing vdat')
        vdat = struct() ;
        vdat.v = zeros(length(seg2d.Vdat), 2) ;
        vdat.NL = zeros(length(seg2d.Vdat), 4) ;
        vdat.fourfold = false(length(seg2d.Vdat), 1) ;
        for qq = 1:length(seg2d.Vdat)
            vdat.v(qq, :) = [seg2d.Vdat(qq).vertxcoord, seg2d.Vdat(qq).vertycoord] ;
            nv = length(seg2d.Vdat(qq).nverts) ;
            try
                vdat.NL(qq, 1:nv) = seg2d.Vdat(qq).nverts ;
            catch
                % Increase the size of NL to accomodate more neighbors
                disp('Increasing NL size (dim 2)')
                swap = vdat.NL ;
                vdat.NL = zeros(length(seg2d.Vdat), nv) ;
                vdat.NL(1:qq, 1:size(swap, 2)) = swap(1:qq, :) ;
                vdat.NL(qq, 1:nv) = seg2d.Vdat(qq).nverts ;
            end
            if isfield(seg2d.Vdat(qq), 'fourfold')
                vdat.fourfold(qq) = ~isempty(seg2d.Vdat(qq).fourfold) ;
            else
                vdat.fourfold(qq) = false ;
            end
        end
        
        disp('generating bond list')
        BL = Vdat2BL(seg2d.Vdat) ;
        vdat.BL = BL ;
        
        % Assign polygons to vdat
        NL = vdat.NL ;
        
        seg2d.vdat = vdat ;
        seg2d.cdat = struct() ;

        % make polygons
        if ~skipPolygons
            polygons = Cdat2polygons(seg2d.Cdat, vdat.v, BL, NL, opts) ;
            seg2d.cdat.polygons = polygons ;
        end
        
        seg2d.cdat.centroid = zeros(length(seg2d.Cdat), 2) ;
        for cid = 1:length(seg2d.Cdat)
            seg2d.cdat.centroid(cid, :) = seg2d.Cdat(cid).centroid.coord ;
        end
        
        %% Save the segmentation to disk
        save(outfn, 'seg2d', 'segIm')
end


function outim = saveSegImage(segIm, segImFn, fn, tp, overwrite, cmpQ)

    %% Save image of the segmentation
    outfn = sprintf(segImFn, tp) ;
    if ~exist(segImFn, 'file') || overwrite
        
        imageFn = sprintf(fn, tp) ;
        im = imread(imageFn) ;
        
        if nargin < 6
            cmpQ = distinguishable_colors(max(segIm(:)), [0,0,0;1,1,1;0,0,1;1,0,0;0,1,0;1,0.10345,0.72414]);
        end
        coloredLabelsImage = label2rgb(segIm, cmpQ, 'k', 'shuffle');
        % Display the pseudo-colored image.
        coloredLabelsImage = (coloredLabelsImage) .* uint8(segIm > 0) ;
        outim = uint8(0.5 * coloredLabelsImage + im) ;
        % outim = imfuse(im, coloredLabelsImage,'blend','Scaling','joint');
        
        disp(['saving image: ' outfn])
        imwrite(outim, outfn) 
        % 
        % clf
        % hold on;
        % axis equal
        % title(['t = ' sprintf('%03d', tp) ' ' timeUnits])
        % saveas(gcf, imfn)
    end
    
end
