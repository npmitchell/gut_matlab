function [meanQ, stdQ, steQ] = QtensorStats(QQ, weights, options)
%[meanQ, stdQ, steQ] = QtensorStats(QQ, weights)
% todo: extend to rank-3+ tensors
%
% Parameters
% ----------
% QQ : #tensors x N x M
%   the observation tensors/matrices (each NxM) to stat
% weights : # tensors x 1 float
%   one weight per observation, will be normalized here
% options : struct with fields
%   struct passed to bootstrapErrorWithWeights
%
% Returns
% -------
% meanQ : NxM float
%   weighted mean, component-by-component
% stdQ : NxM float
%   weighted standard deviation, component-by-component
% steQ : NxM float
%   standard error on the mean, component-by-component. Found via
%   bootstrapping if weights are variable.
%
% See also
% --------
% binDataMeanStdWeighted(x, y, xedges, weights)
% QtensorAspectTheta()
%
% NPMitchell 2021

if nargin < 3 && nargout > 2
    options = struct() ;
end

% normalize weights
nx = size(QQ, 1) ;
if nx == 1
    error('only one tensor to stat. Is this what you want?')
end

ww = weights / nansum(weights) ;
meanQ = squeeze(nansum(ww .* QQ, 1)) ;


if nargout > 1
    % The var function divides by N instead of (N-1), but we don't want this.
    % stdy = sqrt(var(x, ww)) ;
    
    % We choose to divide by (N-1)
    num = nansum(ww .* (QQ-reshape(meanQ, [1, size(meanQ)])).^2) ;
    % denom is very close to 1, slightly smaller because of N-1
    denom = ((nx -1) ./ nx) * nansum(ww) ;
    stdQ = sqrt(num ./ denom) ;

    if isnan(stdQ)
        assert(nx == 1)
        steQ = NaN ;
    else
        % bootstrap for sterror if desired as output var
        if nargout > 2
            steQ = nan(size(meanQ)) ;
            for ii = 1:numel(meanQ) 
                steQ(ii) = bootstrapErrorWithWeights(QQ(:, ii), ww, options) ;
            end
        end
    end
end

