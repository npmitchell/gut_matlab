function outsignal = gaussFilter1D(signal, sigma, windowSize)
%gaussFilter1D(signal, sigma, windowSize)
% Filter a signal along the first dimension with normalized Gaussian. 
% Input signal can be a list of 2D signals of the same length.
%
% Parameters
% ----------
% signal : NxD numeric array
%   The signal(s) to filter along the first dimension (or first of its 
%   dimensions that has more than one element). 
% sigma : int (default=5)
%   stdev of gaussian convolution train
% filterWindowSize : int (default=10)
%   length of gaussFilter vector
%
% Returns
% -------
% outsignal : NxD numeric array
%   filtered signal(s)
%
% NPMitchell 2021

% Default options
if nargin < 2
    sigma = 5 ;
end
if nargin < 3
    windowSize = 10 ;
end

% Make sure signal is the right shape if 1d or 3d
if any(size(signal)) == 1
    if size(signal) > 2
        % If there are only two dims of length > 1, squeeze signal and use
        % as 2d array.
        if sum(size(signal) > 1) == 2
            signal = squeeze(signal) ;
        else
            error('signal must be passed as 1d or 2d array')
        end
    else
        if size(signal, 1) == 1
            signal = reshape(signal, [numel(signal), 1]) ;        
        end
    end
end

% Build the filter to convolve
x = linspace(-windowSize / 2, windowSize / 2, windowSize);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize

% Pad for endcap problems with convolution
fF = [signal(1, :) .* ones(windowSize, size(signal, 2)); ...
    signal; ...
    signal(end, :) .* ones(windowSize, size(signal, 2))];

% Apply filter to each column of signal
for dim = 1:size(signal, 2)
    fF(:, dim) = conv(fF(:, dim), gaussFilter, 'same');
end
outsignal = fF(windowSize+1:end-windowSize, :) ;

end