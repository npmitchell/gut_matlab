function [meanx, stdx, ste] = weightedMeanStdSte(xx, weights, options)
%[meanx, stdx, ste] = weighedMeanStdSte(xx, weights, options)
%
% Parameters
% ----------
% xx : the observations to stat
% weights : one weight per observation, will be normalized here
% options : struct passed to bootstrapErrorWithWeights
%
% Returns
% -------
% meanx : weighted mean
% stdx : weighted standard deviation 
% ste : standard error on the mean via bootstrapping
%
% See also
% --------
% binDataMeanStdWeighted(x, y, xedges, weights)
%
% NPMitchell 2021

if nargin < 3 && nargout > 2
    options = struct() ;
end

% normalize weights
nx = size(xx, 1) ;
if nx == 1
    xx = xx';
    nx = size(xx, 1) ;
end

ww = weights / nansum(weights) ;
meanx = nansum(ww .* xx) ;


if nargout > 1
    % The var function divides by N instead of (N-1), but we don't want this.
    % stdy = sqrt(var(x, ww)) ;
    
    % We choose to divide by (N-1)
    num = nansum(ww .* (xx-meanx).^2) ;
    % denom is very close to 1, slightly smaller because of N-1
    denom = (nx -1) / nx * nansum(ww) ;
    stdx = sqrt(num / denom) ;

    if isnan(stdx)
        assert(nx == 1)
        ste = NaN ;
    else
        % bootstrap for sterror if desired as output var
        if nargout > 2
            ste = bootstrapErrorWithWeights(xx, ww, options) ;
        end
    end
end