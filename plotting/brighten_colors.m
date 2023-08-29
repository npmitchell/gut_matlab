function colors_out = brighten_colors(colors, brighten, permit_clipping)
% brighten a list of colors uniformly without oversaturation
%
% Parameters
% ----------
% colors : Nx3 float
%   rgb values between 0-1
% brighten : float between 0-1 
%   amount to brighten all colors
%
% Returns
% -------
% colors_out : Nx3 float
%   brightened rgb values between 0-1
%
% NPMitchell 2022

if brighten > 0
    hsv = rgb2hsv(colors) ;
    if any(hsv(:, 1) == 1) && ~permit_clipping
        error('at least one color is already saturatated')
    end
    
    hsv(:,1) = hsv(:,1) + brighten ;
    
    if permit_clipping
        hsv(hsv(:,1)>1, 1) = 1 ; 
    end
    
    colors_out = hsv2rgb(hsv) ;
else
    disp('no brightening, since brighten == 0')
end