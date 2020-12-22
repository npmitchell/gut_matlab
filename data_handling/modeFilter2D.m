function modeFilterQuasi1D(gridData, nmodesX, nmodesY)
%[ynew, fft_output] = modeFilter(yy, nmodes, preview)
% FFT reconstruction of lowest nmodes of 1d periodic signal with evenly
% spaced sampling in "time" or other presumably cyclic coordinate.
% Brute force low-pass filter without any fancy Nyquist handling, passband
% ripples, etc. Just select lowest n modes and rebuilt signal from those 
% alone. 
%
% For reference, FFT scalings are:
%     Scale by dt for the FFT, and by Fs for the IFFT
%     Scale by 1/M for the FFT, and by M for the IFFT
%     Scale by 1 for the FFT, and by 1 for the IFFT
%     Scale by 1/sqrt(M) for the FFT, and by sqrt(M) for the IFFT.
%
% Parameters
% ----------
% yy : numeric 1d array
%   input signal to filter
% nmodes : int (default=5)
%   number of modes (including DC offset/average) to use in reconstruction 
% preview : bool (default=false)
%   preview the results in figure form
%
% Returns
% -------
% ynew : numeric 1d array
%   lowpass-filtered output signal reconstructed from lowest nmodes of DFT
%   of input data yy
% fft_output : optional output struct with fields
%   Yin : 
%       two-sided DFT of input data, with phases
%   Yout : 
%       two-sided DFT of filtered output data, with phases
%   SSBand : 
%       single-sided DFT band of filtered output data
%   magnitudes :
%       magnitudes of DFT of filtered output data
%   readme : struct with same fields as above
%       descriptions of output
%
% See also
% --------
% https://medium.com/analytics-vidhya/breaking-down-confusions-over-fast-fourier-transform-fft-1561a029b1ab
%
%
% Example usage
% -------------
% nV = 100 ; % note that period T=1
% tt = linspace(0, (nV-1) / nV, nV-1) ;
% yy = sin(6 * pi * tt) + 0.5 * cos(2 * pi * tt + pi/4) + 0.3 * rand(1, nV-1) ;
% modeFilter(yy, 5, true)
%
% NPMitchell 2020

error('have not written this function')
% Cycle over each y strip
for qq = 1:nU 
    [ynew, data] = modeFilter(ystrip, nmodesY) ;
    dat(qq, :) = ynew ;
    amps(qq, :) = data.amplitudes ;
    thetas(qq, :) = data.theta ;
end

imagesc(dat)


