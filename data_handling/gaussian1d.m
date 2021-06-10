function yy = gaussian1d(xx, mu, sigma)
% normal distribution (gaussian) functional centered around mu, std sigma
%
% Parameters
% ----------
% xx : numeric array
%   evaluation coordinates 
% mu : float
%   mean of gaussian
% sigma : float
%   standard deviation of gaussian
%   
% Returns
% -------
% yy : size of xx, numeric
%   evaluated gaussian function f(xx)

yy = exp(-(xx-mu).^2 ./ (2 * sigma.^2)) ;