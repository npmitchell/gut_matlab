function extractMedianSurface()
function medianSurfaceMIP()
% Isaac B 2020

%% Extract Median Surface using Masked Data
fnp = [fn, '_Probabilities'] ; % Outputted probabilities file name
ssfArray = ones(ssfactor, ssfactor, ssfactor) ; % create cubic 3D array the size of your ssfactor in each dim to perform kronecker tensor product
% initiate a matrix that will store surface data for each tp
surfaceArray = zeros(ltp, 512, 512) ;
% Create coord list for delaunay
[X,Y] = meshgrid(1:512,1:512) ;
xy = [X(:), Y(:)] ;
for tt = 1:ltp % iterate over all tps in the dataset
    close all
    probField = h5read([sprintf(fnp,tt), '.h5'], '/exported_data') ; % outputed prob from Ilastik corresponding to this timepoint
    probField = squeeze(probField(1,:,:,:) ) ;
    probField = permute(probField, [2 3 1]) ; % permute the prob array
    probField = superkron(probField, ssfArray) ; % take kronicker tensor product to "undownsample" the probability field
    % Implement conncomping to clean up the prob field
    probField(probField<0.4) = 0 ; % threshold the probability field
    cPF = bwconncomp(probField) ; % get the concomp struct of the prob field
    seg = zeros(size(probField), 'logical') ; % build a bw array the same size of the prob field (segmentation)
    seg(cPF.PixelIdxList{1}) = true ;
    stats = regionprops3(cPF) ; % extract geometric stats of each volume in the segmentation
    threshold = 1e5 ; % assign the upper threshold of all volumes to be removed
    removeMask = [stats.Volume]<threshold ;
    seg(cat(1,cPF.PixelIdxList{removeMask})) = false ; % remove all bodies with volume < threshold
    [~,~, zidx] = meshgrid(1:size(seg,1), 1:size(seg,2), 1:size(seg,3)) ; % build an array the same size of segmentation but all values in a slice are equal to that slice's z coord
    bw = double(seg) .* zidx ; % mask the array by the segmentation
    bw(bw == 0) = NaN ; % NaN any zeros
    medianz = nanmedian(bw, 3) ; % take the nanmedian of the array (this is median height(x,y) of the segmentaion)
    medianz = inpaint_nans(medianz, 2) ; % interpolate over NaNs
    % smooth and interpolate data
    %medianz = imgaussfilt(medianz) ; % smooth the data using a gaussian
    % store the surface in an array that will later be used to make MIPs
    surfaceArray(tt, :, :) = medianz ;
    tri = delaunay(xy) ;  % create delauney triangulation of xy coords of centdata
    trisurf(tri, xy(:, 1), xy(:,2), medianz(:), 'EdgeColor', 'none')
    caxis([0 105]) ;
    title('Visualization of Muscle Surface for e4') ;
    view(2)
    fns = [fn '_trisurf'] ;
    saveas(gcf, sprintf(fns,tt), 'png') ;
end
save(fullfile(cd, 'surfaceArray.mat'), 'surfaceArray') ;

%% Make Median MIPs
% This section uses the extracted median surface to make MIPs of the region
% of data that lives nearby the median surface.
if ~exist('surfaceArray') % if the var surfaceArray does not exist, load it in from hardrive
    load(fullfile(cd, 'surfaceArray.mat'), 'surfaceArray') ;
end
% Now make overlayed MIP
for tt = 1:size(surfaceArray, 1)
    mask = round(squeeze(surfaceArray(tt, :, :))) ; % extract z(x,y) position from surfaceArray
    xp.loadTime(tt) ; % load current TP
    xp.rescaleStackToUnitAspect() ; % rescale to unit aspect ratio
    data = xp.stack.image.apply() ; % save to a variable
    for dataChannel = 1:2
        rawData = data{dataChannel} ; % load the appropriate channel
        [xx,yy] = meshgrid(1:size(rawData, 1), 1:size(rawData, 2));
        mask3d = false(size(rawData)) ; % initiate a boolean 3D array that will mask the rawData
        xmax = size(rawData, 1) ;
        ymax = size(rawData, 2) ;
        zmax = size(rawData, 3) ;
        if dataChannel == 1
            lowbnd = -15 ;
            upbnd = 20 ;
        elseif dataChannel == 2
            lowbnd = 23 ;
            upbnd = zmax ;
        end
        for chadd = lowbnd:upbnd
            % what page to set to true?
            zind = mask + chadd ;
            % convert to linear index
            lin = (zind - 1) * xmax * ymax + (xx-1)*ymax + yy;
            % clip linear indices to be inside volume
            lin = lin(lin > 0) ;
            lin = lin(lin < xmax * ymax * zmax) ;
            % Set this layer to be true
            mask3d(lin) = true ;
        end
        maskedData= uint8(rawData).*uint8(mask3d) ; % mask the raw data using our boolean 3d mask
        maskedData = imadjust(uint8(max(maskedData, [], 3))) ; % take the max value along the z-dim and adjust to new map
        if dataChannel == 1
            overlay = cat(3, zeros(size(maskedData)) ,maskedData, maskedData) ;
        elseif dataChannel == 2
            overlay = cat(3, maskedData, zeros(size(maskedData)), maskedData) ;
        end
        imshow(overlay) ;
        fnMasked = 'midfold_e4_T%02d_Ch%01d' ;
        imwrite(maskedData, fullfile(dataDir, 'overlays', ['masked_ch',num2str(dataChannel)], [sprintf(fnMasked, tt, dataChannel), '_overlayMIP_masked.tiff']))  ;
        close all
    end
end