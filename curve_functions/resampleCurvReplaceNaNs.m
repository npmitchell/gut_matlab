function [curvout] = resampleCurvReplaceNaNs(curv, N, closed)
%resampleCurvReplaceNaNs Resample equally a D-dim curve at N points
%   The endpoints are each part of the sampling, so if this is a closed
%   curve, remove the final endpoint (or startpoint) for an equal sampling.
%   Note: convention is that a closed curve will have equal starting and
%   ending points (so the position appears twice).
%
% Parameters
% ----------
% curv : M x D float array
%   The curve to resample such that each segment is equal Euclidean
%   length
% N : int
%   the number of points desired in the sampling
% closed : bool
%   whether the curve is closed (final point == starting point)
% 
% Returns
% -------
% curvout : N x D float array
%   the resampled curve with equally spaced segment sampling
%
% NPMitchell 2019

if any(isnan(curv))
    startbad = any(isnan(curv(1,:))) ;
    endbad = any(isnan(curv(end,:))) ;
    if startbad
        if ~endbad 
            if closed
                curv(1, :) = curv(end, :) ;
            else
                error('Interpolate startpoint here.')
            end
        else
            error('start and end are bad. Handle this case here.')
        end
    elseif endbad
        if ~startbad
            if closed
                curv(end, :) = curv(1, :) ;
            else
                error('Interpolate final endpoint here.')
            end
        end
    else
        error('NaNs in the middle of curve. Handle this case by interpolation')
    end
end
curvout = curvspace(curv, N) ;

end

