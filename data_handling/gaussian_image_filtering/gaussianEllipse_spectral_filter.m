function [imf, handles] = gaussian_spectral_filter(im, fsz1, fsz2, angle, hi_lo, preview, cmap)
% Gaussian Filter Response Calculation
%
% Parameters
% ----------
% im : grayscale image/matrix
% fsz : int
%   filter size in frequency space (sigma of gaussian)
% hi_lo : string 'hi' or 'lo'
% preview : bool
% cmap : colormap to use for preview of results
%
%
% NPMitchell 2021, based on Gaussian_image_filtering.m from 
% https://www.mathworks.com/matlabcentral/fileexchange/46812-two-dimensional-gaussian-hi-pass-and-low-pass-image-filter?s_tid=mwa_osa_a
% but with their typo corrected in the equation for a gaussian, input
% handling, rescaling, formatting, and wrapping into function

% Input handling
if nargin < 2
    fsz1 = 10; % filter size parameter 
elseif isempty(fsz1)
    fsz1 = 10 ;
end
if nargin < 3
    fsz2 = 10; % filter size parameter 
elseif isempty(fsz2)
    fsz2 = 10 ;
end
if nargin < 4
    angle = 0; % filter size parameter 
elseif isempty(angle)
    angle = 0 ;
end
if nargin < 5
    hi_lo = 'lo' ;
elseif isempty(hi_lo)
    hi_lo = 'lo' ;
else
    hi_lo =lower(hi_lo) ;
end
if nargin < 6
    preview = false ;
end
if nargin < 7
    if preview
        cmap = viridis ;
    end
elseif isempty(cmap)
    cmap = viridis ;
end


% Take FFT
A = fft2(double(im)); % compute FFT of the grey image
A1 = fftshift(A); % frequency scaling

% Gaussian Filter Response Calculation
[M, N] = size(A); % image size
X = 0:N-1;
Y = 0:M-1;
[X, Y] = meshgrid(X,Y);
Cx = 0.5 * N;
Cy = 0.5 * M;
if angle == 0
    Lo = exp(-((X-Cx).^2 / (2*fsz1.^2))-((Y-Cy).^2 ./ (2*fsz2.^2)));
else
    error('handle tilted gaussian here')
end

% Filtered image=ifft(filter response*fft(original image))
if contains(hi_lo, 'lo')
    J = A1.*Lo;
    J1 = ifftshift(J);
    imf = ifft2(J1);
    % rescale to original amplitude
    % imf = imf ;
elseif contains(hi_lo, 'hi')
    Hi = 1-Lo; % High pass filter=1-low pass filter
    K = A1.*Hi;
    K1 = ifftshift(K);
    imf = ifft2(K1);
    % rescale to original amplitude
    % imf = imf ;
else
    error('hi_lo must contain hi or lo in string specification')
end

%----visualizing the results----------------------------------------------
if preview
    clf
    subplot(2, 2, 1)
    imagesc(im');
    colorbar
    h1 = gca ;
    title('original image','fontsize',14)

    subplot(2, 2, 2)
    climlo = prctile(abs(A1(:)), 0.1) ;
    climhi = prctile(abs(A1(:)), 99) ;
    imagesc(abs(A1)', [climlo climhi])
    colorbar
    h2 = gca ;
    title('fft of original image','fontsize',14)

    if contains(hi_lo, 'lo')
        subplot(2, 2, 3)
          imagesc(abs(imf)')
          colorbar
          h3 = gca; 
          title('low passed','fontsize',14)
        
        subplot(2, 2, 4)
          mesh(X,Y,Lo)
          axis([ 0 N 0 M 0 1])
          h4 = gca; 
          title('Gaussian LPF H(f)','fontsize',14)
          colormap viridis

    elseif contains(hi_lo, 'hi')
        subplot(2, 2, 3)
          imshow(abs(imf)')
          colorbar
          h3 = gca; 
          title('High passed','fontsize',14)

        subplot(2, 2, 4)
          mesh(X,Y,Hi)
          axis([ 0 N 0 M 0 1])
          h4 = gca; 
          title('Gaussian HPF H(f)','fontsize',14)
          colormap(cmap)
    end
    handles = [h1 h2 h3 h4] ;
else
    handles = [] ;    
end
