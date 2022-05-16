function pkg = pkg(n)
% purple black green colormap diverging
%
% NPMitchell 2022

if nargin < 1
    n = 256 ;
end
pkg = interpolate_3color_cmap(linspace(0, 1, n), ...
    [185., 110, 255.]/255, [0,0,0], [185., 255,110]/255) ;