function processCorrectedCellSegmentation2D(QS, options)
% Unfinished code -- Lin working on it.
% ToDo: 
%
%
% tissueAnalysisSuite fields are different than QuapSlap's. 
% For tissueAnalysisSuite, we have:
% Vdat: 
%     nverts : neighbor list for vertices
%     ncells : 
%     vertxcoord : column of data in which vertex lives
%     vertycoord : row of data in which each vertex lives
% Cdat : cell data    
%     ncells : indices of neighboring cells
%     nverts : vertices that define the cell
%     centroid.coord(1) x position of cell centroid
%     centroid.coord(2) y position of cell centroid
% Bdat : 
%     nverts
%     ncells
%     pix : linear indices of the pixels associated with that bond
%
% Note on exterior calculus objects
% d0 is an e x c matrix of exterior derivatives with +1 and -1s
% at the endpts of each bond
% d1 is a v x e matrix of exterior derivatives. Upstream is +1,
% downstream is -1 when moving counterclockwise around a
% tension plaquette.
% d0 and d1 are matrices that take derivatives 

%% Parameters
overwrite = false ;
overwriteImages = false ;
% timepoints to process
timePoints = QS.xp.fileMeta.timePoints ;
% how far in pixels is too far to connect two cell vertices
very_far = 250 ;
% which coordinate system to use for segmentation
coordSys = QS.currentSegmentation.coordSys ; 
% Toggle for iLastik version control -- zero for newer version
iLastikVersion = 0;
cellSize = 50 ;
strelRadius = 0 ;
gaussKernel = 1 ;
heightMinimum = 3.5 ;


%% unpack options
if nargin < 2 
    options = struct() ;
end

if isfield(options, 'overwrite') 
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages') 
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'timePoints') 
    timePoints = options.timePoints ;
end
if isfield(options, 'coordSys') 
    coordSys = options.coordSys ;
end

%% Load in h5 from ilastik.
if strcmpi(erase(coordSys, '_'), 'spsme') 
    Folder = [QS.dir.im_sp_sme, '_pixelClassification'] ;
    if ~exist(Folder, 'dir')
        mkdir(Folder)
        error(['Populate ' Folder ' with pixelClassification on pullbacks with coordSys ' coordSys])
    end
    coordSys = 'spsme' ;
elseif strcmpi(erase(coordSys, '_'), 'sprsme') || ...
        strcmpi(erase(coordSys, '_'), 'rspsme') || ...
        strcmpi(erase(coordSys, '_'), 'rsme')
    Folder = [QS.dir.im_r_sme, '_pixelClassification'] ;
    if ~exist(Folder, 'dir')
        mkdir(Folder)
        error(['Populate ' Folder ' with pixelClassification on pullbacks with coordSys ' coordSys])
    end
    coordSys = 'sprsme' ;
else
    error('Have not coded for this coordinate system yet. Do so here')
end

for tp = timePoints
    
    % The results file with segmentation and polygons is called outfn
    outfn = sprintf(QS.fullFileBase.segmentation2dCorrected, coordSys, tp)  ;
    if strcmpi(erase(coordSys, '_'), 'spsme') 
        imageFn = sprintf(QS.fullFileBase.im_sp_sme, tp) ; 
    elseif strcmpi(erase(coordSys, '_'), 'sprsme') 
        imageFn = sprintf(QS.fullFileBase.im_r_sme, tp) ;
    end
    
    if ~exist(outfn, 'file') || overwrite

        bwfn = sprintf(QS.fullFileBase.segmentation2dCorrectedBinary, coordSys, tp) ;
        bw = imread(bwfn) ;
        
        cc = bwconncomp(~bw, 4) ;
        segIm = labelmatrix(cc) ;

        % Check it -- Apply a variety of pseudo-colors to the regions.
        coloredLabelsImage = label2rgb(segIm, 'hsv', 'k', 'shuffle');
        % Display the pseudo-colored image.
        % [row,col] = find(segIm < 2) ;
        % coloredLabelsImage(row, col, :) = 0 ;
        imshow(coloredLabelsImage);
        title('Pseudocolored Labeled Image', 'Interpreter', 'None');
        % impixelinfo;
        
        coordSys = lower(erase(coordSys, '_')) ;
        if strcmpi(coordSys, 'spsme') || ...
                strcmpi(coordSys, 'sprsme') || ...
                strcmpi(coordSys, 'rspsme') || ...
                strcmpi(coordSys, 'rsme')
            % provide ROI in [minx, maxx; miny, maxy]
            ROI = [-eps, size(segIm, 2) + eps; round([0.25, 0.75] * size(segIm, 1))] ;
            adjustY = 0.5 * size(segIm, 1) ;
        else
            error('Handle this coordSys here')
        end
        props = regionprops(segIm, 'centroid') ;
        
        seg2d.roi = ROI ;
        % Could find vertices by seeking each point in all vertices, bonds
        % by seeking starting and ending points of shared linesegments
        % between pairs of cells. This seems unnecessary for now.
        % seg2d.vdat = vdat ;
        seg2d.cdat = struct() ;
        seg2d.cdat.centroid = zeros(max(segIm(:)), 2) ;
        for cid = 1:max(segIm(:))
            seg2d.cdat.centroid(cid, :) = props(cid).Centroid ;
        end
        
        xmin = seg2d.roi(1, 1) ;
        xmax = seg2d.roi(1, 2) ;
        ymin = seg2d.roi(2, 1) ;
        ymax = seg2d.roi(2, 2) ;
        insideROI = inpolygon(seg2d.cdat.centroid(:, 1), ...
            seg2d.cdat.centroid(:, 2), ...
            [xmin, xmax, xmax, xmin], [ymin, ymin, ymax, ymax]) ;
        if any(~insideROI)
            c2move = find(~insideROI) ;
            segIm2 = segIm ;
            imagesc(segIm) 
            hold on;
            plot([xmin, xmax, xmax, xmin], [ymin, ymin, ymax, ymax], '-') ;
            pause(1) 
            clf
            for ccId = 1:length(c2move)
                qq = c2move(ccId) ;
                % Check if we need to push it up in Y or down
                insideUp = inpolygon(seg2d.cdat.centroid(qq, 1), ...
                    seg2d.cdat.centroid(qq, 2) + adjustY, ...
                    [xmin, xmax, xmax, xmin], [ymin, ymin, ymax, ymax]) ;
                insideDown = inpolygon(seg2d.cdat.centroid(qq, 1), ...
                    seg2d.cdat.centroid(qq, 2) - adjustY, ...
                    [xmin, xmax, xmax, xmin], [ymin, ymin, ymax, ymax]) ;
                if insideUp || insideDown
                    % get all indices where segIm == c2move
                    [row,col] = find(segIm == qq) ;
                        
                    if insideUp && ~insideDown
                        for id = 1:length(row)
                            % If this isn't merging this cell with another
                            if segIm(row(id) + adjustY, col(id)) < 2
                                segIm(row(id) + adjustY, col(id)) = qq ;
                            end
                            segIm(row(id), col(id)) = 1 ;
                        end
                    elseif insideDown && ~insideUp
                        for id = 1:length(row)
                            % If this isn't merging this cell with another
                            if segIm(row(id) - adjustY, col(id)) < 2                            
                                segIm(row - adjustY, col) = qq ;
                            end
                            segIm(row(id), col(id)) = 1 ;
                        end
                    else 
                        error('Cannot place cell into the ROI')
                    end
                    disp('Check that this motion is correct')
                    imagesc(segIm)
                    pause(0.0001)
                    % title('Press any button to continue')
                    % waitforbuttonpress
                else 
                    error('Cannot place cell into the ROI')
                end
                
            end
            
            % Recompute centroids
            seg2d.cdat = struct() ;
            seg2d.cdat.centroid = zeros(max(segIm(:)), 2) ;
            for cid = 1:max(segIm(:))
                seg2d.cdat.centroid(cid, :) = props(cid).Centroid ;
            end
        end
        
        % Now convert to polygons
        polygons = polygonsFromSegmentation(segIm) ;        
        seg2d.cdat.polygons = polygons ;
        
        %% Save the segmentation to disk
        save(outfn, 'seg2d', 'segIm', 'coordSys')
    else
        disp(['already on disk: ' outfn])
        
        % %% Convert to simpler format
        load(outfn, 'seg2d', 'segIm', 'coordSys')
    end
    
    %% Save image of the segmentation
    imfn = [outfn(1:end-3) 'png'] ;
    if ~exist(imfn, 'file') || overwrite || overwriteImages
        imageFn = sprintf(QS.fullFileBase.im_sp_sme, tp) ;
        im = imread(imageFn) ;
        colors = define_colors() ;
        im2 = imoverlay(im, segIm==0, colors(3, :)) ;
        imwrite(im2, imfn)
    end
    
    
    %% Save image of the polygons
    imfn = [outfn(1:end-4) '_polygons.png'] ;
    if ~exist(imfn, 'file') || overwrite || overwriteImages || true
        im = imread(imageFn) ;
        
        % Check it -- Apply a variety of pseudo-colors to the regions.
        % iters = ceil(max(segIm(:)) / 256) ;
        % cmpQ = tab10 ;
        % for ii = 1:iters
        %     cmpQ = [cmpQ; tab10] ;
        % end
        cmpQ = brewermap(max(segIm(:)), 'Paired') ;
        coloredLabelsImage = label2rgb(segIm, cmpQ, 'k', 'shuffle');
        % Display the pseudo-colored image.
        coloredLabelsImage = (coloredLabelsImage) .* uint8(segIm > 1) ;
        outim = uint8(0.5 * coloredLabelsImage + im) ;
        % h = imshow(outim); 
        
        % saveas(gcf, imfn)
        disp(['saving image: ' imfn])
        % current_ax = getframe(gca) ;
        % imwrite(current_ax.cdata, imfn)
        imwrite(outim, imfn) 
    end
end
