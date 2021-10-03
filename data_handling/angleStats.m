function [meanA, stdA, steA, ste2] = angleStats(thetas, weights) 
% Compute mean, std, ste of angles on the unit circle as sqrt(-log (R^2))
% with optional weights associated with each measurement. 
%%
%
% Parameters
% ----------
% thetas : Nx1 float
%   angles to compute stats for
% weights : Nx1 float
%   optional weights associated with each observation
%
% Returns
% -------
% stdA : float
%   stdev of (weighted) angle measurement
% meanA : float
%   (weighted) mean of observations on unit circle
% steA : float
%   standard error on the mean of angles. If weights are provided, this is
%   obtained via bootstrapping.
% steA : float
%   If weights are provided, this standard error on the mean of angles is
%   obtained via bootstrapping in a more transparent and rigorous fashion,
%   where we vary the subsampling number (n) from the total N samples, fit
%   the variance to m*(1/n)+b, and plug in n=N to the fit. 
% 
%
% NPMitchell 2021

if nargin > 1
    ww = weights / nansum(weights) ;
    ct = nansum(ww .* cos(thetas)) ;
    st = nansum(ww .* sin(thetas)) ;
else
    ct = mean(cos(thetas)) ;
    st = mean(sin(thetas)) ;
end
% Mean angle is given by arctan(st/ct)
meanA = atan2(st, ct) ;

if nargout > 1
    % return standard deviation
    stdA = sqrt(-log(ct^2 + st^2)) ;
    
    % return standard error
    if nargout > 2
        % unequal weights
        bootstatC = bootstrp(1000,@mean,cos(thetas),'weights', ww);
        steC = std(bootstatC) ;
        bootstatS = bootstrp(1000,@mean,sin(thetas),'weights', ww);
        steS = std(bootstatS) ;
        
        % Error propagation from errors on the mean of ct and st:
        % \sqrt{\frac{\text{$\delta $y}^2 x^2+\text{$\delta $x}^2 y^2}{\left(x^2+y^2\right)^2}}
        steA = sqrt(steC^2*st^2 + ct^2*steS^2) / ...
            sqrt( (ct^2 + st^2)^2 ) ;
        
        if nargout > 3
            steC = bootstrapErrorWithWeights(cos(thetas), ww) ;
            steS = bootstrapErrorWithWeights(sin(thetas), ww) ;
            
            % Error propagation from errors on the mean of ct and st:
            % \sqrt{\frac{\text{$\delta $y}^2 x^2+\text{$\delta $x}^2 y^2}{\left(x^2+y^2\right)^2}}
            ste2 = sqrt(steC^2*st^2 + ct^2*steS^2) / ...
                sqrt( (ct^2 + st^2)^2 ) ;
            
            % % Test against single size
            % % get distribution of means 
            % meany = zeros(1000, 1) ;
            % for pp = 1:1000
            %     % subsample from xx
            %     [yy, idx] = datasample(cos(thetas), length(cos(thetas)), 'weights', ww) ;
            %     meany(pp) = mean(yy) ;
            % 
            %     % Alternatively, weight in sum:
            %     % renormalize weights
            %     % wk = weights(idx) / nansum(weights(idx)) ;
            %     % meany(pp) = nansum(wk .* yy) ;
            % end
            % % Look at variance of means
            % test = var(meany) ;
        end
    else
        % equal weights
        steA = stdA ./ sqrt(size(thetas, 1)) ;
    end
    
end