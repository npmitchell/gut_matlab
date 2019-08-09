%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Go through many images that are different sizes in pwd, plot them on a 
% canvas that is as large as the biggest one so that one can make a single 
% tiff of all the images.
% SJS and NPM 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tInit = 110 ; 
tEnd  = 242 ;
step  = 1;
chan  = 1;

nameStruct = ['cmp_',num2str(chan),'_1_T%04d.tif'];

% determine the maximum with and height in the image series; 
w = zeros(1,tEnd-tInit+1);
h = zeros(1,tEnd-tInit+1);
s = struct('inf',[]);
tInd = tInit : step: tEnd;
for time = 1:length(tInd)
    name = sprintf(nameStruct,tInd(time));
    %name = ['Cyl2_p_',num2str(time),'.tif'];
    inf = imfinfo(name);
    s(time).inf = inf; 
    
    w(time) = s(time).inf.Width;
    h(time) = s(time).inf.Height;
end

% Now take the max of both width and height
mw = max(w);
mh = max(h);

% adjust the canvas size by initializing large empty image
im = zeros(mh,mw,tEnd-tInit+1,'uint16');
for time = 1:length(tInd)
    name = sprintf(nameStruct,tInd(time));
    %name = ['Cyl2_p_',num2str(time),'.tif'];
    temp = imread(name);
    
    st = size(temp);
    stD = floor(([mh,mw]-st)/2)+1;
    
    im( stD(1):(stD(1)+st(1)-1) , stD(2):(stD(2)+st(2)-1) ,time) = temp;
end

% Write the result as a tiff 
for k = 1 : size(im,3)    
    imwrite(im(:,:,k),['adjust',num2str(chan),'.tif'],'tiff','Compression','none','WriteMode','append');
end
