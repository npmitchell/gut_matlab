function [dxint, dyint] = PhaseCorrelationPotentialBound(img1, img2, potential, dim, varargin)
%Description:
%   Find the translation shift between two image using Kuglin & Hines 1975
%
% Parameters
% ----------
% img1 : N x M float or int array
%   image 1
% img2 : N x M float or int array
%   image 2
% potential : str ('linear', 'quadratic', 'linearY', 'quadraticY')
%   string specifier for potential to apply to translation
% dim : int
%   0 = isotropic potential, 1 = along x only, 2 = along y only
% varargin : optional input arguments
%   LowerBoundX : 
%   UpperBoundX :
%   LowerBoundY :
%   UpperBoundY :
%   Verbose : bool
%
% Returns
% -------
% x : int
%   single-pixel resolution translation between img1 and img2
% y : int
%   single-pixel resolution translation between img1 and img2
%
% NPMitchell 2020

% Default value for number of iterations until give up search in bounds
niter = 100 ;
lowerx = -Inf ;
lowery = -Inf ;
upperx = Inf ;
uppery = Inf ;
verbose = false ;
sigma_filt = 3 ; 
potential_sigmay = 0 ;
if ~isempty(varargin)
    check_bounds = true ;
    for i = 1:length(varargin)
        if isa(varargin{i},'double') 
            continue;
        end
        
        if isa(varargin{i},'logical')
            continue;
        end

        if ~isempty(regexp(varargin{i},'^[Ll]ower[Bb]ound[Xx]','match'))
            lowerx = varargin{i+1} ;
        end
        if ~isempty(regexp(varargin{i},'^[Uu]pper[Bb]ound[Xx]','match'))
            upperx = varargin{i+1} ;
        end
        if ~isempty(regexp(varargin{i},'^[Ll]ower[Bb]ound[Yy]','match'))
            lowery = varargin{i+1} ;
        end
        if ~isempty(regexp(varargin{i},'^[Uu]pper[Bb]ound[Yy]','match'))
            uppery = varargin{i+1} ;
        end
        
        % Options for filtering the inverse FFT of phase correlation
        if ~isempty(regexp(varargin{i},'^[Ss]igma[Ff]ilter','match'))
            sigma_filt = varargin{i+1} ;
        end
        if ~isempty(regexp(varargin{i},'^[Pp]otential[Ss]igma[Yy]','match'))
            potential_sigmay = varargin{i+1} ;
        end
        
                
        if ~isempty(regexp(varargin{i},'^[Vv]erbose','match'))
            verbose = varargin{i+1} ;
        end
    end
else
    check_bounds = false ;
end

if verbose
    disp('Bounds are:')
    disp(['x=(' num2str(lowerx) ', ' num2str(upperx) '), '])
    disp(['y=(' num2str(lowery) ', ' num2str(uppery) ')'])
end

% Phase correlation (Kuglin & Hines, 1975)
%============================================================
af = fftn(double(img1));
bf = fftn(double(img2));
denom = abs(af .* conj(bf)) ;
denom(denom == 0) = 1 ;
cp = af .* conj(bf) ./ denom ;
icp = (ifft2(cp));

% Bound the translation with a potential
if ~isempty(regexp(potential, '^[Ll]inear','match'))
    boundMat = ones(size(icp)) ;
    % Bound the Y dimension with linear potential
    if dim ~= 1
        lin = fliplr(0:ceil(size(icp, 1)*0.5)) ;
        pot = [lin fliplr(lin)] ;
        pot = pot(1:size(icp, 1)) ;
        if potential_sigmay > 0 
            pot = pot - potential_sigmay ;
            pot(pot < 0) = 0 ;
        end
        pot = pot / max(pot) ;
        boundMat = pot' .* boundMat ;
    end
    % Bound the X dimension with linear potential
    if dim ~= 2
        error('do this dim')
    end
elseif ~isempty(regexp(potential, '^[Qq]uadratic','match'))
    boundMat = ones(size(icp)) ;
    % Bound the Y dimension with linear potential
    if dim ~= 1
        lin = 0:ceil(size(icp, 1)*0.5) ;
        lin = max(lin(:)).^2 - lin.^2 ;
        pot = [lin fliplr(lin)] ;
        pot = pot(1:size(icp, 1)) ;
        if potential_sigmay > 0 && any(pot > potential_sigmay.^2)
            pot = pot - potential_sigmay.^2 ;
            pot(pot < 0) = 0 ;
        end
        pot = pot / max(pot) ;
        boundMat = pot' .* boundMat ;
    end
    % Bound the X dimension with linear potential
    if dim ~= 2
        error('do this dim')
    end
elseif ~isempty(regexp(potential, '^[Ss]qrt','match'))
    error('here') 
elseif ~isempty(regexp(potential, '^[Nn]one','match'))
else
    error(['Could not recognize potential ' potential ': must be linear, quadratic, sqrt, or none'])
end

% So long as the bounding potential is not none, apply it
if isempty(regexp(potential, '^[Nn]one','match'))
    icpBounded = icp .* boundMat ;
else
    icpBounded = icp ;
end

% Smooth the correlation peaks if requested
if sigma_filt > 0
    if verbose
        disp('Filtering icpBounded with Gaussian convolution')
    end
    icpBounded = imgaussfilt(icpBounded, sigma_filt, 'Padding', 'circular') ; 
end

% Obtain the maximum 
mmax = max(max(icpBounded));
%Find the main peak
[x, y, ~] = find(mmax == icpBounded) ;

% Inspect the result
if verbose
    figure; 
    % Show the bounding potential if not none
    if isempty(regexp(potential, '^[Nn]one','match'))
        imagesc(boundMat) ; colorbar() ;
        title(['max = ', num2str(max(boundMat(:)))])
        pause(3) 
    end
    imagesc(icpBounded) ; colorbar() ; % ,[],'notruesize');
    pause(3)
    close all
end

%Note: could extend to Subpixel Registration [Foroosh, Zerubia & Berthod, 2002]
%============================================================
[M, N] = size(img1);

% Also translate the single pixel resolution answer
if x < (N/2)
    dyint = x - 1;
else
    dyint = x - M - 1;
end

if y < (M/2)
    dxint = y - 1;
else
    dxint = y - N - 1;
end


% Check that the peaks are inside the supplied bounds
if check_bounds
    kii = 2 ;
    force_exit = false ;
    while (dxint < lowerx || dxint > upperx || dyint < lowery || dyint > uppery) && ~force_exit
        if verbose
            disp(['Current peak is ' num2str(x) ', ' num2str(y)] )
            disp('ExttPhaseCorrelationPotentialBound: peak out of bounds, searching next peak...')
        end
        [~, xyind] = maxk(icpBounded(:), kii) ;
        [x, y] = ind2sub(size(icpBounded), xyind) ;
        x = x(kii) ;
        y = y(kii) ;
        
        % Also translate the single pixel resolution answer
        if x < (N/2)
            dyint = x;
        else
            dyint = x - M;
        end

        if y < (M/2)
            dxint = y;
        else
            dxint = y - N;
        end
        
        % Update the iteration index and try again
        kii = kii + 1 ;
        if kii > niter
            disp(['Current peak is ' num2str(x) ', ' num2str(y)] )
            msg = ['exiting: could not find peak within bounds after niter=' niter ' iterations'] ;
            msg = [msg '. Bounds are:' ] ;
            msg = [msg 'x=(' num2str(lowerx) ', ' num2str(upperx) '), '] ;
            msg = [msg 'y=(' num2str(lowery) ', ' num2str(uppery) ')'] ;
            disp(msg)
            disp('Exiting with NaNs')
            dxint = NaN ;
            dyint = NaN ;
            force_exit = true ;
        end
    end        
end
