function GG = unpackManualIlastikGroundTruthH5(h5fileBase, timePoints, ...
    imOutDir, maxN_for_plot, rawImFileBase)
%
%
%
%
%
%

if ~strcmp(h5fileBase(end-2:end), '.h5')
    h5fileBase = [h5fileBase '.h5'] ;
end

if nargin < 2
    timePoints = 1:length(dir(h5fileBase)) ;
    try
        assert(length(timePoints) > 0)
    catch
        error('No matching files found')
    end
end
if nargin < 3
    imOutDir = [] ;
end
if nargin < 4
    maxN_for_plot = 5000 ;
end
if ~isempty(imOutDir)
    jets = jet(maxN_for_plot) ;
    jetshuffle = jets(randperm(maxN_for_plot), :) ;
end

mfmoveList = [] ;
for tidx = 1:length(timePoints)
    tp = timePoints(tidx) ;
    disp(['condidering t = ' num2str(tp)])
    tracksfn_tt = sprintf(h5fileBase, tidx-1) ;
    
    im = h5read(tracksfn_tt, '/segmentation/labels') ;
    im = squeeze(im(1, :, :)) ;
    
    % Query all tracked objects among those segmented
    if tidx > 1
        moves = h5read(tracksfn_tt, '/tracking/Moves')' ;
        try
            appr = h5read(tracksfn_tt, '/tracking/Appearances')' ;
        catch
            appr = [] ;
        end
        try
            disa = h5read(tracksfn_tt, '/tracking/Disappearances')' ;
        catch
            disa = [] ;
        end
        try
            splt = h5read(tracksfn_tt, '/tracking/Splits')' ;
        catch
            splt = [] ;
        end
        try
            merg = h5read(tracksfn_tt, '/tracking/Mergers')' ;
        catch
            merg = [] ;
        end
        
        try 
            mfmoves = h5read(tracksfn_tt, '/tracking/MultiFrameMoves') ;
            disp('Found multi-frame moves')
            if ~isempty(mfmoveList)
                mfmoveList = [mfmoveList; cat(2, mfmoves', tp * ones(size(mfmoves', 1), 1))];
            else
                mfmoveList = cat(2, mfmoves', tp * ones(size(mfmoves', 1), 1)) ;
            end
        catch
            mfmoves = [] ;
        end
        
        regp = regionprops(im, 'Centroid') ;
        com1 = vertcat(regp.Centroid) ;
        
        % Thread tracks
        if tidx == 2
            % Assign first tp track endpts
            spt = com0(moves(:, 1), :) ;
            allCellCenters = spt ;
            numCellsT = size(spt, 1) ;
            allCellTimes = repmat(timePoints(1), numCellsT, 1) ;
            allCellTimes0 = allCellTimes ;
            prevEndptSegIdx = moves(:, 1) ;
            prevSegmentLabels = (1:length(prevEndptSegIdx))' ;
            allSegmentLabelsArr = prevSegmentLabels ;
        end
        sptIdx = moves(:, 1) ;
        
        % All current COM from moves, appearances, etc
        eptIdx = [moves(:, 2), appr] ;
        ept = com1(eptIdx, :) ;
        numCellsT = size(ept, 1) ;
        count0 = length(allCellTimes0) ;
        
        % Create edges of graph
        globalIdxStart = zeros(size(sptIdx)) ;
        donePrevEndptIdx = [] ;
        newTrackId = max(allSegmentLabelsArr(:)) + 1 ;
        newSegmentLabels = zeros(size(sptIdx)) ;
        for qq = 1:length(sptIdx)
            sqq = sptIdx(qq) ;
            % Find spt Index in previous endpoints
            match = find(prevEndptSegIdx == sqq) ;
            if length(match) == 1
                % There is a single match, so connect the current idx to
                % this single previous endpoint (now start point)
                
                % This is the node index: number of nodes before previous
                % timepoint (count0) + index of node in previous 
                % timepoint's endpoint node list.
                globalIdxStart(qq) = match + count0 ;
                try
                    newSegmentLabels(qq) = prevSegmentLabels(match) ;
                catch
                    error('Could not index into prevSegmentLabels')
                end
            elseif ~isempty(match)
                % There are multiple matches --> the current idx split
                match2keep = ~ismember(match, donePrevEndptIdx) ;
                match = match(match2keep) ;
                if ~isempty(match)
                    globalIdxStart(qq) = match(1) + count0 ;
                    % Mark that we have used this previous endpoint already
                    donePrevEndptIdx = [ donePrevEndptIdx, match(1) ] ;
                    newSegmentLabels(qq) = prevSegmentLabels(match(1)) ;
                else
                    disp('new track??')
                    % assign new track label to segment label array
                    newSegmentLabels(qq) = newTrackId ;
                    newTrackId = newTrackId + 1 ;
                end
            else
                disp('could not find start point in previous endpoint list -- likely multi-Frame Move')
                newSegmentLabels(qq) = newTrackId ;
                newTrackId = newTrackId + 1 ;
            end
        end
        globalIdxEnd = ((1:size(moves, 1)) + size(sptIdx, 1) + count0)' ;
        newEdges = [globalIdxStart, globalIdxEnd] ;
        
        % Remove any rows that have null matches (not sure why these exist)
        [row, col] = find(newEdges == 0) ;
        newEdges(row, :) = [] ;
        
        if tidx == 2
            edgesT = newEdges ;
        else
            edgesT = [edgesT; newEdges] ;
        end
        try
            assert(all(edgesT(:) > 0))
        catch
            error('Some edges are assigned to be zero nodes!')
        end
        
        % Update the global cell center list
        allCellCenters = [allCellCenters; ept];
        
        % Update the times in the global list
        allCellTimes0 = allCellTimes ;
        allCellTimes =  [allCellTimes; repmat(tp, numCellsT, 1)];
    
        % Update the segment labels (track IDs)
        allSegmentLabelsArr = [allSegmentLabelsArr; newSegmentLabels ] ;
        
        % assign current centers of mass to previous array (for next tp)
        com0 = com1 ;
        % current timepoint's endpoint segmentation indices are now old
        prevEndptSegIdx = eptIdx ;
        % current timepoint's segment labels (track IDs) are now old
        prevSegmentLabels = newSegmentLabels ;
        
        assert(numel(eptIdx) == numel(newSegmentLabels)) ;
        
    else
        regp = regionprops(im, 'Centroid') ;
        com0 = vertcat(regp.Centroid) ;
        disp('Waiting for second tp to assign object tracks in first tp')
        im0 = im ;
    end
    
    % Output image of tracks if output image directory specified
    if ~isempty(imOutDir)
        if tidx > 1
            if tidx == 2
                raw0 = imread(sprintf(rawImFileBase, timePoints(1))) ;
                im2 = 0*im0 ;
                for qq = 1:length(sptIdx)
                    im2(im0 == sptIdx(qq)) = newSegmentLabels(qq) ;
                end
                if any(newSegmentLabels > maxN_for_plot)
                    im2 = mod(im2, maxN_for_plot) ;
                end
                rgb0 = label2rgb(im2, jetshuffle, [0,0,0]);
                zeroIdx = all(rgb0 == 0, 3) ;
                red = squeeze(rgb0(:, :, 1)) ;
                grn = squeeze(rgb0(:, :, 2)) ;
                blu = squeeze(rgb0(:, :, 3)) ;
                red(zeroIdx) = raw0(zeroIdx) ;
                grn(zeroIdx) = raw0(zeroIdx) ;
                blu(zeroIdx) = raw0(zeroIdx) ;
                rgb0 = cat(3, red, grn, blu) ;

                if ~exist(imOutDir, 'dir')
                    mkdir(imOutDir)
                end
                imwrite(rgb0, fullfile(imOutDir, ...
                    sprintf('tracks_color_%06d.png', timePoints(1))))
            end

            raw = imread(sprintf(rawImFileBase, tp)) ;
            im2 = 0*im ;
            for qq = 1:length(eptIdx)
                im2(im == eptIdx(qq)) = newSegmentLabels(qq) ;
            end
            if any(newSegmentLabels > maxN_for_plot)
                im2 = mod(im2, maxN_for_plot) ;
            end
            
            rgb = label2rgb(im2, jetshuffle, [0,0,0]);
            zeroIdx = all(rgb == 0, 3) ;
            red = squeeze(rgb(:, :, 1)) ;
            grn = squeeze(rgb(:, :, 2)) ;
            blu = squeeze(rgb(:, :, 3)) ;
            red(zeroIdx) = raw(zeroIdx) ;
            grn(zeroIdx) = raw(zeroIdx) ;
            blu(zeroIdx) = raw(zeroIdx) ;
            rgb = cat(3, red, grn, blu) ;

            imwrite(rgb, fullfile(imOutDir, ...
                sprintf('tracks_color_%06d.png', tp)))
        end
    end
end

% Convert segment Label array into cell
nCells = size(allCellTimes, 1) ;
allSegmentLabels = cell(nCells, 1) ;
ndigits = num2str(length(num2str(max(allSegmentLabelsArr)))) ;
for qq = 1:length(allSegmentLabelsArr)
    allSegmentLabels{qq} = sprintf(...
        ['%0' ndigits 'd'], allSegmentLabelsArr(qq)) ;
end
allGenerations = ones([nCells, 1]);

% A table containing node properties with which to construct the digraph
nodeTable = table( allCellCenters, allCellTimes, ...
    allSegmentLabels, allGenerations, ...
    'VariableNames', {'UPix', 'T', 'Segment', 'Generation'} );


% Construct the graph and update node properties
GG = digraph;
GG = GG.addnode(nodeTable);

for eID = 1:size(edgesT, 1)
    
   % Update the edge list in the digraph
   GG = addedge(GG, edgesT(eID, 1), edgesT(eID, 2));
    
end
