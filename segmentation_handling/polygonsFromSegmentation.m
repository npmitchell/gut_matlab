function polygonsFromSegmentation(labelMat) 


for obj = 1:max(labelMat(:))
    bw = labelMat == obj ;
    bnd = bwboundaries(imdilate(bw), 4);
    boundaries{obj} = bnd ;
end

% Check it -- Apply a variety of pseudo-colors to the regions.
coloredLabelsImage = label2rgb(labelMat, 'hsv', 'k', 'shuffle');
% Display the pseudo-colored image.
subplot(2, 2, 4);
imshow(coloredLabelsImage);
title('Pseudocolored Labeled Image', 'FontSize', fontSize, 'Interpreter', 'None');
impixelinfo;

%--------------------------------------------------------------------------------------------------------
% MEASUREMENT OF BOUNDARIES
boundaries = bwboundaries(binaryImage, 4);