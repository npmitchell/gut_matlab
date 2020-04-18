function generateMaskedDataSingleMask(QS, dir16bit)
% This is a standalone script that works fine, but doesn't actually seem to
% give good results for the gut. We instead want a time-varying mask, as is
% done in generateMaskedData()

% First create an average stab image to use for training
avgImageFn = fullfile(dir16bit, 'avgImage.h5') ;
maskFn = fullfile(dir16bit, 'avgImage_Probabilities.h5') ;
maskedDir = fullfile(dir16bit, 'masked_data') ;
step = 3 ;
preview = QS.plotting.preview ;
ssfactorMask = 4 ;
last_tp = xp.fileMeta.timePoints(end) ;
tptodo = xp.fileMeta.timePoints([1:10:100, 105:10:last_tp]) ;
if exist(avgImageFn, 'file') && ~overwrite_avgImage
    disp('avgImageFn already on disk -- train in ilastik')
else
    for tt = tptodo
        disp(['Loading t = ' num2str(tt)])
        xp.loadTime(tt);
        xp.rescaleStackToUnitAspect();
        IV = xp.stack.image.apply() ;
        if tt > xp.fileMeta.timePoints(1)
            avgIm = avgIm + double(IV{1}) ;
        else
            avgIm = double(IV{1}) ;
        end
        clf
        imagesc(squeeze(avgIm(:, :, 300)))
        title(['t= ' num2str(tt)])
    end
    avgImage = avgIm / length(tptodo) ;
    vlo = double(prctile(avgImage(:) , 1.00 )) / double(max(avgImage(:))) ;
    vhi = double(prctile(avgImage(:) , 99.99 )) / double(max(avgImage(:))) ;
    disp(['--> ' num2str(vlo) ', ' num2str(vhi)])
    avgImage2 = imadjustn(mat2gray(avgImage), [double(vlo); double(vhi)]) * 2^16 ;
    avgImage2 = uint16(avgImage2(1:ssfactorMask:end, 1:ssfactorMask:end,...
        1:ssfactorMask:end)) ;
    % Write it
    h5create(avgImageFn, '/input_data', size(avgImage2), ...
        'Datatype','uint16');
    h5write(avgImageFn, '/input_data', avgImage2) ;
    disp('done writing average Image')
end

% Mask data, rotate to lateral and ventral views
maskprob = h5read(maskFn, '/exported_data') ;
mask = squeeze(maskprob(1, :, :, :)) > 0.5 ;
% Take largest connected component, erode and dilate.
se3 = strel('sphere', 3);
se4 = strel('sphere', 4);
mask0 = imerode(mask, se3) ;
mask0 = imdilate(mask0, se4) ;
if preview
    for qq = 1:step:size(mask0, 3)
        imshow(squeeze(mask0(:, :, qq)))
        title('eroded & dilated volume')
    end
end
mask0 = imerode(mask0, se3) ;
mask0 = imdilate(mask0, se4) ;
if preview
    for qq = 1:step:size(mask0, 3)
        imshow(squeeze(mask0(:, :, qq)))
        title('eroded & dilated volume -- pass 2')
    end
end
CC = bwconncomp(mask0) ;
numPixels = cellfun(@numel,CC.PixelIdxList) ;
[vol, bigID] = max(numPixels) ;
biggest = zeros(size(mask0), 'logical') ;
biggest(CC.PixelIdxList{bigID}) = 1 ;
if preview
    for qq = 1:step:size(biggest, 3)
        imshow(squeeze(biggest(:, :, qq)))
        title('biggest volume')
    end
end
% Now fill holes by keeping only the biggest outside region as zero
mask1 = (1-biggest) ;
CC = bwconncomp(mask1) ;
numPixels = cellfun(@numel,CC.PixelIdxList) ;
[vol, bigID] = max(numPixels) ;
bwsolid = ones(size(mask1), 'logical') ;
bwsolid(CC.PixelIdxList{bigID}) = 0 ;
if preview
    for qq = 1:step:size(bwsolid, 3)
        title('filled solid')
        imshow(squeeze(bwsolid(:, :, qq)))
    end
end

% blur it & threshold to smoothen
blurred = imgaussfilt3(uint16(bwsolid), 20) ;
bwsolid = blurred > 0.5 ;

% active contour to smoothen more
bwsolid = activecontour(squeeze(maskprob(1, :, :, :)), bwsolid, ...
    10, 'Chan-Vese',...
    'ContractionBias', -0.1, 'SmoothFactor', 5);
mask_final = superkron(bwsolid, ones(ssfactorMask, ssfactorMask, ssfactorMask)) ;

% blur the mask to soften edges
mask_final = imdilate(mask_final, se3) ;
mask_final = imgaussfilt3(uint16(mask_final), 20) ;

if preview
    for qq = 1:step:size(bwsolid, 3)
        title('final mask')
        imshow(squeeze(mask_final(:, :, qq)))
        title('final mask')
    end
end


%% Apply the mask to the data
old_adjusthigh = QS.data.adjusthigh ;
QS.data.adjusthigh = 99.999 ;
for tt = QS.xp.fileMeta.timePoints(30:30:end)
    disp(['Loading t = ' num2str(tt)])
    QS.setTime(tt)
    QS.getCurrentData() ;
    IV = QS.currentData.IV ;
    IVmasked = IV ;
    
    % clip the mask if larger than data
    if tt == QS.xp.fileMeta.timePoints(1) 
        mask_final = mask_final(1:size(IV{1}, 1), ...
            1:size(IV{1}, 2), 1:size(IV{1}, 3)) ; 
    end
    
    % inspect it
    for ii = 1:length(IV)
        if nChannels > 1
            out16name = fullfile(maskedDir, sprintf([fn '_masked.tif'], tt, ii)) ;
        else
            out16name = fullfile(maskedDir, sprintf([fn '_masked.tif'], tt)) ;
        end
        IVii = IV{ii} ;
        IVmasked{ii} = IVii .* uint16(mask_final) ;
        
        % check it
        if preview
            clf
            for qq = 1:step:size(IVmasked{ii}, 3)
                page = squeeze(IVmasked{ii}(:, :, qq)) ;
                bwpage = squeeze(uint16(mask_final(:, :, qq)*2^16)) ;
                imshow(cat(3, 0.5*bwpage, page, page + 0.5*bwpage))
                title(['t= ' num2str(tt)])
                pause(0.1)
                % imshow(bwpage)
            end
        end
        
        % Save each image
        disp(['writing data to ' out16name])
        out = uint16(IVmasked{ii}) ;
        for z = 1:size(out, 3)
            if mod(z, 100) == 0
                disp(['page ' num2str(z) ' / ' num2str(size(out,3))])
            end
            % write the first page as overwrite mode, then append
            if z == 1
                imwrite(out(:,:,z), out16name, 'tiff',...
                    'Compression','none');
            else
                imwrite(out(:,:,z), out16name, 'tiff',...
                    'Compression','none','WriteMode','append');    
            end
        end
    end
    
end