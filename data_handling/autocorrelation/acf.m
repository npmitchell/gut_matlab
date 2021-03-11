function ta = acf(y,p, normalize, meanSubtraction, plotResult)
% ACF - Compute Autocorrelation/Autocovariance Through p Lags
% If normalize == true, then normalize by variance of timeseries so that
% the result is a pearson correlation coefficient. 
% >> myacf = acf(y, lags, normalize, plotResult) 
%
% correlation = Expectation[ X(t) conj(X(t+tau)) ]
% covariance = Expectation[ (X(t)-mu) conj(X(t+tau)-mu)) ]
%
% Inputs:
% y - series to compute acf for, nx1 column vector
%   Note: can be complex to accomodate 2d data. There is a conjugation
%   added in the definition for complex numbers:
%   
% p - total number of lags, 1x1 integer
% normalize - boolean (default=true)
%   compute covariance rather than correlation (divide by variance)
% meanSubtraction - boolean (default=true)
%   subtract mean in expression of covariance/correlation
% plotResult : boolean (default=true)
%   Plot the results 
%
% Output:
% myacf - px1 vector containing autocorrelations
%        (First lag computed is lag 1. Lag 0 not computed)
%
%
% A bar graph of the autocorrelations is also produced, with
% rejection region bands for testing individual autocorrelations = 0.
%
% Note that lag 0 autocorelation is not computed, 
% and is not shown on this graph.
%
% Example:
% >> acf(randn(100,1), 10)
%


% --------------------------
% USER INPUT CHECKS
% --------------------------

[n1, n2] = size(y) ;
if n2 ~=1
    error('Input series y must be an nx1 column vector')
end

[a1, a2] = size(p) ;
if ~((a1==1 && a2==1) && (p<n1))
    error('Input number of lags p must be a 1x1 scalar, and must be less than length of series y')
end

if nargin < 5
    plotResult = true ;
end
if nargin < 4
    meanSubtraction = true ;
end
if nargin < 3
    normalize = true ;
end

% -------------
% BEGIN CODE
% -------------

ta = zeros(p,1) ;
global N 
N = max(size(y)) ;
global ybar 
ybar = mean(y); 

% Collect ACFs at each lag i
for i = 1:p
    ta(i) = acf_k(y,i, normalize, meanSubtraction) ;
end

%% Plot ACF
if plotResult
    % Plot rejection region lines for test of individual autocorrelations
    % H_0: rho(tau) = 0 at alpha=.05
    bar(ta)
    line([0 p+.5], (1.96)*(1/sqrt(N))*ones(1,2))
    line([0 p+.5], (-1.96)*(1/sqrt(N))*ones(1,2))

    % Some figure properties
    line_hi = (1.96)*(1/sqrt(N))*1.05;
    line_lo = -(1.96)*(1/sqrt(N))*1.05;
    bar_hi = max(ta)*1.05 ;
    bar_lo = -max(ta)*1.05 ;

    xlim([0 p+.60])
    if (abs(line_hi) > abs(bar_hi)) % if rejection lines might not appear on graph
        ylim([line_lo line_hi])
    else
        ylim([bar_lo bar_hi])
    end
    title({' ','Sample Autocorrelations',' '})
    xlabel('Lag Length')
    % set number of lag labels shown
    if (p<28 && p>4)
        set(gca,'XTick',floor(linspace(1,p,4)))
    elseif (p>=28)
        set(gca,'XTick',floor(linspace(1,p,8)))
    end
    set(gca,'TickLength',[0 0])
end

% ---------------
% SUB FUNCTION
% ---------------
function ta2 = acf_k(y,k, normalize, meanSubtraction)
% ACF_K - Autocorrelation at Lag k
% acf(y,k)
%
% Inputs:
% y - series to compute acf for
% k - which lag to compute acf
% 
global ybar
global N
cross_sum = zeros(N-k,1) ;

% Numerator, unscaled covariance
if meanSubtraction
    for i = (k+1):N
        cross_sum(i) = conj(y(i)-ybar) * (y(i-k)-ybar) ;
    end
else
    for i = (k+1):N
        cross_sum(i) = conj(y(i)) * (y(i-k)) ;
    end
end

if normalize
    % Denominator, unscaled variance
    yvar = (y-ybar)'*(y-ybar) ;
    ta2 = sum(cross_sum) / yvar ;
else
    ta2 = sum(cross_sum)  ;
end



