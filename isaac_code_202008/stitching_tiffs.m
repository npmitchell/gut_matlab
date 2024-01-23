%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TURN COLLECTION OF TIFFS INTO A MOSOAIC (STITCHING)
%%%%
%%%% Isaac Breinyn 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

clear all, clc

root_dir = '/mnt/data/development_videos/201911051100_byncheck/byn_grid_40xbf_orcaflash4.0_650nm_10po_1/' ;
cd(root_dir) % Directory with image folders

oper = 78/100; % overlap percentage (i.e. 10/100 --> 10%)

% iterate over each tif slice
for t = 1:141 
    % Stitch together all rows (hotdog --> hamburger)
    for r = 0:24 
        % Stitch together all images in one row (image --> hotdog)
        for c = 0:25  
            cd(fullfile(root_dir, 'Tiles'))
            fnr = ['Pos_', sprintf('%03d', r), '_', sprintf('%03d', c)];
            if c == 0
                R{c+1} = imread(fnr,t) ;  %save image data to array R
            else
                R{c+1} = imread(fnr,t);
                R{c+1}(:,(512-round(oper*512)):512) = [] ; % clip the overlap from the correct side of the image
                R{c+1} = cat(2, R{c+1}, R{c}) ;
            end
        end
        stitched_row = R{c+1};
        imwrite(stitched_row, ['Stitched_row_', sprintf('%02d',r),'_Time_', sprintf('%03d', t), '.tif']);
        fnc = ['Stitched_row_', sprintf('%02d',r),'_Time_', sprintf('%03d', t), '.tif'] ;
        if r == 0
            C{r+1} = imread(fnc) ;  %save image data to array C
        else
            C{r+1} = imread(fnc);
            C{r+1}((1:round(oper*512)),:) = [] ; % clip the overlap from the correct side of the image
            C{r+1} = cat(1, C{r}, C{r+1}) ;
        end
        stitched_col = C{r+1};
        cd(fullfile(root_dir, 'Stitched Tiffs'))
        imwrite(stitched_col, ['StitchedBynCheck_Time', sprintf('%03d%', t), '.tif']);
    end
end

%% Turn the outputted tiffs into one jumbo tiff
% Split this up so as not to reach 4GB max for tiffs
for t = 1:40
    fn = ['StitchedBynCheck_Time', sprintf('%03d%', t), '.tif']
    if t == 1
        imwrite(fn, 'StitchedBynCheck_Final_1_40.tif')
    else
        imwrite(fn, 'StitchedBynCheck_Final_1_40.tif', 'WriteMode', 'append')
    end
end

for t = 41:80
    fn = ['StitchedBynCheck_Time', sprintf('%03d%', t), '.tif']
    if t == 41
        imwrite(fn, 'StitchedBynCheck_Final_41_80.tif')
    else
        imwrite(fn, 'StitchedBynCheck_Final_41_80.tif', 'WriteMode', 'append')
    end
end

for t = 81:120
    fn = ['StitchedBynCheck_Time', sprintf('%03d%', t), '.tif']
    if t == 81
        imwrite(fn, 'StitchedBynCheck_Final_81_120.tif')
    else
        imwrite(fn, 'StitchedBynCheck_Final_81_120.tif', 'WriteMode', 'append')
    end
end

for t = 121:160
    fn = ['StitchedBynCheck_Time', sprintf('%03d%', t), '.tif']
    if t == 121
        imwrite(fn, 'StitchedBynCheck_Final_121_160.tif')
    else
        imwrite(fn, 'StitchedBynCheck_Final_121_160.tif', 'WriteMode', 'append')
    end
end

for t = 161:200
    fn = ['StitchedBynCheck_Time', sprintf('%03d%', t), '.tif']
    if t == 161
        imwrite(fn, 'StitchedBynCheck_Final_161_200.tif')
    else
        imwrite(fn, 'StitchedBynCheck_Final_161_200.tif', 'WriteMode', 'append')
    end
end

for t = 201:240
    fn = ['StitchedBynCheck_Time', sprintf('%03d%', t), '.tif']
    if t == 201
        imwrite(fn, 'StitchedBynCheck_Final_201_240.tif')
    else
        imwrite(fn, 'StitchedBynCheck_Final_201_240.tif', 'WriteMode', 'append')
    end
end





