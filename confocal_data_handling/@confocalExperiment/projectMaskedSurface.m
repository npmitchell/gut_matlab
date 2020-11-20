function [maskedData, medianValues] = projectMaskedSurface(rawData, probField , ssfactor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Isaac B Breinyn 2020
%
% Args :
% rawData : the data to be masked and projected
% probField : the field to use as a mask
% ssfactor : the sub sampling factor between the Data and the Mask
%
% Outputs :
% maskedData : the rawData masked using the probability field (same size as
% rawData and probField)
% medianValues : the height of the median of the surface of maskedData (N-1
% dims than rawData and probField)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Making the masked Data

ssfArray = ones(ssfactor, ssfactor, ssfactor) ;

probField = superkron(probField, ssfArray) ; % take kronicker tensor product to "upsample" the probability field

rawData = superkron(rawData, ssfArray) ; % take kronicker tensor product to "upsample" the probability field


% remove all disconnected objects from the probField so that the only thing
% left is a single, continuous surface
% 
% cPF = bwconncomp(probField) ;
% pSize = length(cPF.PixelIdxList) ;
% probabilityField = bwareaopen(probField, pSize) ;

% Threshold the probability field and then make it binary

% probabilityField(probField > 0.7) = 1 ;
% probabilityField(probField < 0.7) = 0 ;

% Now mask the raw data with the cleaned prob field

maskedData = uint16(rawData.*probField) ; % and there is the masked Data

%% extracting median values

[X,Y] = meshgrid(1:size(rawData, 1), 1:size(rawData, 2)) ;
xy = [X(:), Y(:)] ;

[~,~, zidx] = meshgrid(1:size(probField,2), 1:size(probField,1), 1:size(probField,3)) ; % build an array the same size of segmentation but all values in a slice are equal to that slice's z coord
bw = double(probField) .* zidx ; % mask the array by the segmentation
bw(bw == 0) = NaN ; % NaN any zeros

medianz = nanmedian(bw, 3) ; % take the nanmedian of the array (this is median height(x,y) of the segmentaion)
medianValues = inpaint_nans(medianz, 2) ; % interpolate over NaNs

end
