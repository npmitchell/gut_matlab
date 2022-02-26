function com = com_region2d(probability_grid, thres, xyzgrid) 

if nargin < 2
    thres = 0.5 ;
end
if nargin < 3
    xyzgrid = [] ;
end

% Find big region, decreasing thres if needed
bwccOk = false ;
thresIter = thres ;
while ~bwccOk
    try
        bwcc = bwconncomp(probability_grid > thresIter) ; 
        npix = cellfun(@numel,bwcc.PixelIdxList);
        [~, indexOfMax] = max(npix); 
        % isolate the largest connected component
        biggest = zeros(size(probability_grid));
        biggest(bwcc.PixelIdxList{indexOfMax}) = 1;
        bwccOk = true ;
    catch
        thresIter = thresIter * 0.9 ;
        disp(['Could not segment big region, decreasing threshold to ' num2str(thresIter)])
        if thresIter < 0.01
            error('Threshold of 0.01 still did not yield a connected component')
        end
    end
end
% Now multiply with probability
mass = probability_grid .* biggest ;
% Get center of mass. There are two ways
if ~isempty(xyzgrid)
    % Method 1, using mean
    mass = probability_grid ;
    mean_mass = mean(mass(:)) ;
    comX = mean(mass(:) .* xyzgrid(1, :, :, :)) / mean_mass ;
    comY = mean(mass(:) .* xyzgrid(2, :, :, :)) / mean_mass ;
     com = [ comX comY ] ; 
else 
    props = regionprops(true(size(mass)), mass, 'WeightedCentroid');
    com_tmp = props.WeightedCentroid ;
    % Note that we flip XY here since WeightedCentroid treats the 1st
    % dimension as the 2nd (and 2nd as 1st)
    com = [ com_tmp(2) com_tmp(1) ] ;
end
    
end