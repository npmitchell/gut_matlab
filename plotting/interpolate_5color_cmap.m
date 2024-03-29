function[map] = interpolate_5color_cmap(s,rgb1, rgb2, rgb3, rgb4, rgb5, normalize, gammaCorrection)
%This function is based on Kenneth Moreland's code for creating Diverging
% Colormaps, modified by Andy Stein. Created by Noah P Mitchell 2019
% template accessed from: https://www.kennethmoreland.com/color-maps/
%s is a vector that goes between zero and one 
% Common rgb1, rgb2, ... values are defined below, so that rgb1 could be
% set to 1 or 2 or 3 etc instead of [x,y,z].
% 
% Parameters
% ----------
% s : array with values spanning [0, 1]
%   The values to be mapped to colors, spanning [0, 1]
% rgb1 : length(3) array or int
%   An RGB triplet for value 0 or an integer indexing common rgb triplets
% rgb2 : length(3) array or int
%   An RGB triplet for value 0.25 or an integer indexing common rgb triplets
% rgb3 : length(3) array or int
%   An RGB triplet for value 0.5 or an integer indexing common rgb triplets
% rgb4 : length(3) array or int
%   An RGB triplet for value 0.75 or an integer indexing common rgb triplets
% rgb5 : length(3) array or int
%   An RGB triplet for value 1 or an integer indexing common rgb triplets
% normalize : bool (default=false)
%   normalize each row of the RGb colormap
% gammaCorrection : bool (default = false)
%   perform gamma correction to make colors more perceptually uniform
%
% Output
% ------
% map : colormap
%   A mapping from [0, 1] to colors spanning rgb1 and rgb2, with length(s)
%   possible values.
% 
% Example Usage
% -------------
% % Purple to black to green:
% cbkry = interpolate_5color_cmap(linspace(0, 1, 256), ...
%           [0,1,1],[0,0,1],[0,0,0],[1,0,0],[1,1,0]))
%

    if nargin < 7
        normalize = false ;
    end

    if nargin < 8
        gammaCorrection = false ;
    end

    % If the arguments rgb1 and rgb2 are ints, then use these arguments 
    % to index from a list of commonly used rgb values
    rgb_1 = [0.230, 0.299, 0.754] ;  % blue
    rgb_2 = [0.706, 0.016, 0.150] ;  % red
    rgb_3 = [0.436, 0.308, 0.631] ;  % dark purple
    rgb_4 = [0.759, 0.334, 0.046] ;  % off red
    rgb_5 = [0.217, 0.525, 0.910] ;  % light blue
    rgbs = [rgb_1; rgb_2; rgb_3; rgb_4; rgb_5] ;
    if all(size(rgb1) == 1)
        rgb1 = rgbs(rgb1, :) ;
    end
    if all(size(rgb2) == 1)
        rgb2 = rgbs(rgb2, :) ;
    end
    if all(size(rgb3) == 1)
        rgb3 = rgbs(rgb3, :) ;
    end
    if all(size(rgb4) == 1)
        rgb4 = rgbs(rgb4, :) ;
    end
    if all(size(rgb5) == 1)
        rgb5 = rgbs(rgb5, :) ;
    end

    map = zeros(length(s),3);
    for i=1:length(s)
        map(i,:) = fivecolor_map_1val(s(i),rgb1,rgb2, rgb3, rgb4, rgb5, normalize, gammaCorrection);
    end
end

% Interpolate a diverging color map.
function[result] = fivecolor_map_1val(s, rgb1, rgb2, rgb3, rgb4, rgb5, normalize, gammaCorrection)
    % s1 is a number between 0 and 1 tuning through the colormap
    % rgb1,2,3 are RGB triplets

    % First color goes from 1 to zero in (0, 0.25)
    s1 = max(0, 1 - 4 * s) ;
    % Second color goes from 0 to 1 to 0 in (0, 0.5)
    s2 = max(0, 1 - abs(4 * (s - 0.25))) ;
    % Third color goes from 0 to 1 to 0 in (0.25, 0.75)
    s3 = max(0, 1 - abs(4 * (s - 0.5))) ;
    % Fourth color goes from 0 to 1 to 0 in (0.5, 1.)
    s4 = max(0, 1 - abs(4 * (s - 0.75))) ;
    % Fifth color goes from 0 to 1 (0.75, 1.)
    s5 = max(4 * s - 3, 0) ;

    result(1) = s1*rgb1(1) + s2*rgb2(1) + s3*rgb3(1)+ ...
        s4*rgb4(1) + s5*rgb5(1);
    result(2) = s1*rgb1(2) + s2*rgb2(2) + s3*rgb3(2)+ ...
        s4*rgb4(2) + s5*rgb5(2);
    result(3) = s1*rgb1(3) + s2*rgb2(3) + s3*rgb3(3) + ...
        s4*rgb4(3) + s5*rgb5(3);

    % Renormalize the result
    if normalize
        result = result / vecnorm(result) ;
    end
    % OPTIONAL: gammacorrection
    if gammaCorrection
        xyz = RGBToXYZ(result) ;    
        result = XYZToRGB(xyz);
    end
end
