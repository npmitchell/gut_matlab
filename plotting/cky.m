function cky = cky(n)
% cyan black yellow colormap diverging
%
% NPMitchell 2022

if nargin < 1
    n = 256 ;
end
cky = interpolate_3color_cmap(linspace(0, 1, n), ...
    [0,1,1], [0,0,0], [1., 1, 0.]) ;