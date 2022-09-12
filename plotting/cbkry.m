function cbkry = cbkry(n)
% cyan blue black red yellow colormap diverging, with uniform V shape in
% lightness.
% 
% Note:
% tmp = rgb2hsv(cbkry(256))
% plot(tmp(:, 3)) % shows V shape in brightness.
%
% NPMitchell 2022

if nargin < 1
    n = 256 ;
end
cbkry = interpolate_5color_cmap(linspace(0, 1, n), ...
    [0.5,1,1], [0.25,0.,0.5], [0,0,0], [0.5,0,0.25], [1., 1, 0.5]) ;