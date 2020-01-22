function extendImages(directory, direc_e, fileNameBase, a_fixed, ntiles, overwrite)
% EXTENDIMAGES(directory, direc_e, fileNameBase) Repeat an image above and
% below
%
% directory : str
%   path to the existing images
% direc_e : str
%   path to the place where extended images are to be saved
% fileNameBase : str
%   The file name of the images to load and save
% a_fixed : float
%   The aspect ratio of the pullback image: Lx / Ly
% ntiles : int 
%   The number of bins in each dimension for histogram equilization for a
%   square original image. That is, the extended image will have (a_fixed *
%   ntiles, 2 * ntiles) bins in (x,y).
% overwrite : bool
%   overwrite the existing extended image on disk
%
% NPMitchell 2019 

fns = dir(strrep(fullfile([directory, '/', fileNameBase, '.tif']), '%06d', '*')) ;
% Get original image size
im = imread(fullfile(fns(1).folder, fns(1).name)) ;
halfsize = round(0.5 * size(im, 1)) ;
% osize = size(im) ;

% Extend each timepoint's image
for i=1:length(fns)
    if ~exist(fullfile(direc_e, fns(i).name), 'file') || overwrite 
        % Declare if we are overwriting the file
        if exist(fullfile(direc_e, fns(i).name), 'file') 
            disp(['Overwriting ' fullfile(direc_e, fns(i).name)])
        end
        disp(['Reading ' fns(i).name])
        % fileName = split(fns(i).name, '.tif') ;
        % fileName = fileName{1} ;
        im = imread(fullfile(fns(i).folder, fns(i).name)) ;

        % im2 is as follows:
        % [ im(end-halfsize) ]
        % [     ...          ]
        % [    im(end)       ]
        % [     im(1)        ]
        % [     ...          ]
        % [    im(end)       ]
        % [     im(1)        ]
        % [     ...          ]
        % [  im(halfsize)    ]
        im2 = uint8(zeros(size(im, 1) + 2 * halfsize, size(im, 2))) ;
        im2(1:halfsize, :) = im(end-halfsize + 1:end, :);
        im2(halfsize + 1:halfsize + size(im, 1), :) = im ;
        im2(halfsize + size(im, 1) + 1:end, :) = im(1:halfsize, :);
        im2 = adapthisteq(im2, 'NumTiles', [round(a_fixed * ntiles), round(2 * ntiles)]) ;
        imwrite( im2, fullfile(direc_e, fns(i).name), 'TIFF' );
    else
        disp('already exists')
    end

end