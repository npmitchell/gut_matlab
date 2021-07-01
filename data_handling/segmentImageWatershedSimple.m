function segmentImageWatershedSimple(prob, cellSize, strelRadius, gaussKernel)
% 
% Simple method for segmenting cells, following Nick Noll's procedure
%
% Nick Noll & NPMitchell 2021

%% WATERSHED to segment cells
% Filter output: laplacian of Gaussian sharpens, then Gaussian smooths        
h1 = fspecial('log', cellSize);
seD1 = strel('disk', strelRadius);

% Pack label image L with watershed results
mem = mat2gray(prob, [0, double(max(prob(:)))]) ;
% L = zeros(size(mem));
disp(['memWS: segmenting timepoint ', num2str(tp)]) 
if gaussKernel > 0
    gk = fspecial('gaussian', gaussKernel);
    mem = imfilter(mem,gk);
end
cyto = 1 - mem;

% mem = imfilter(mem,h1);
% mem(:,:,t) = imclose(mem(:,:,t),strel('disk',3));

lev = graythresh(cyto);
seed = im2bw(cyto, lev);  % consider swapping to imbinarize
seed = imdilate(seed, seD1);
seed = bwareaopen(seed, 25);

pre_water = imhmin(mem, heightMiminum);
pre_water = imimposemin(pre_water, seed);
LL = watershed(pre_water);
skel = LL == 0 ;
