%Implementation of:
%   Extension of Phase Correlation to Subpixel Registration
%   (Hassan Foroosh, Josiane B. Zerubia, Marc Berthod)
%Implemented by:
%   Lulu Ardiansyah (halluvme@gmail.com)
%   Noah P Mitchell added bounds, bounding functions, cleanup & arguments
%
%TODO:
%   - Find out whether this implementation is correct :)
%   - Combine the result, overlay the images based on the result
%
%                                               eL-ardiansyah
%                                               January, 2010
%                                                       CMIIW
% extensive edits and additions by NPMitchell 2020
%============================================================
function [deltaX , deltaY] = ExtPhaseCorrelation(img1, img2, varargin)
%Description:
%   Find the translation shift between two image
%
%Parameter:
%   img1 = image 1
%   img2 = image 2
%   image 1 and image 2 , must in the same size

% Default value for number of iterations until give up search in bounds
niter = 100 ;
lowerx = -Inf ;
lowery = -Inf ;
upperx = Inf ;
uppery = Inf ;
verbose = false ;
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
mmax = max(max(icp));

%Find the main peak
[x, y, ~] = find(mmax == icp) ;

% tmp = figure; imagesc(icp) % ,[],'notruesize');
% error('here')

%Extension to Subpixel Registration [Foroosh, Zerubia & Berthod, 2002]
%============================================================
[M, N] = size(img1);

%two side-peaks
xsp = x + 1;
ysp = y + 1;

%if the peaks outsize the image, then use xs-1 and/or ys-1 for the two
%side-peaks
if xsp > M
    xsp = M-1;
end
if ysp > N
    ysp = N-1;
end

%Calculate the translational shift
deltaX1 = ((icp(xsp,y) * xsp + icp(x,y) * x) / (icp(xsp,y) + icp(x,y)))-1;
deltaY1 = ((icp(x,ysp) * ysp + icp(x,y) * y) / (icp(x,ysp) + icp(x,y)))-1;

% Note: The result of deltaX and delta Y is inverted.
% Validate if translation shift is negative
if deltaX1 < (N/2)
    deltaY = deltaX1;
else
    deltaY = deltaX1 - M;
end

if deltaY1 < (M/2)
    deltaX = deltaY1;
else
    deltaX = deltaY1 - N;
end


% Check that the peaks are inside the supplied bounds
if check_bounds
    kii = 2 ;
    while deltaX < lowerx || deltaX > upperx || deltaY < lowery || deltaY > uppery
        disp('ExtPhaseCorrelation: peak out of bounds, consider next peak...')
        [~, xyind] = maxk(icp(:), kii) ;
        [x, y] = ind2sub(size(icp), xyind) ;
        x = x(kii) ;
        y = y(kii) ;
        
        %two side-peaks
        xsp = x + 1;
        ysp = y + 1;

        %if the peaks outsize the image, then use xs-1 and/or ys-1 for the two
        %side-peaks
        if xsp > M
            xsp = M-1;
        end
        if ysp > N
            ysp = N-1;
        end

        %Calculate the translational shift
        deltaX1 = ((icp(xsp,y) * xsp + icp(x,y) * x) / (icp(xsp,y) + icp(x,y)))-1;
        deltaY1 = ((icp(x,ysp) * ysp + icp(x,y) * y) / (icp(x,ysp) + icp(x,y)))-1;

        % Note: The result of deltaX and delta Y is inverted.
        % Validate if translation shift is negative
        if deltaX1 < (N/2)
            deltaY = deltaX1;
        else
            deltaY = deltaX1 - M;
        end

        if deltaY1 < (M/2)
            deltaX = deltaY1;
        else
            deltaX = deltaY1 - N;
        end
        
        % Update the iteration index and try again
        kii = kii + 1 ;
        if kii > niter
            error(['exiting: could not find peak within bounds after niter=' niter ' iterations'])
        end
    end        
end
