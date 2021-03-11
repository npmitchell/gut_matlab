function ta = acfVector(y,p, plotResult)
% ACF - Compute Autocorrelations Through p Lags
% >> myacf = acf(y,p) 
%
% Inputs:
% y - series to compute acf for, nxd column vector
% p - total number of lags, 1x1 integer
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

if nargin < 3
    plotResult = true ;
end

[n1, n2] = size(y) ;
if n2 ==1
    disp('WARNING: input is 1d, using acf()')
    acf(y, p, plotResult)
else

    [a1, a2] = size(p) ;
    if ~((a1==1 && a2==1) && (p<n1))
        error('Input number of lags p must be a 1x1 scalar, and must be less than length of series y')
    end

    % -------------
    % BEGIN CODE
    % -------------

    ta = zeros(p,1) ;
    Ny = size(y,1) ;
    ybar = mean(y, 1); 

    % Collect ACFs at each lag i
    for i = 1:p
       ta(i) = acf_Vk(y, i, ybar, Ny) ; 
    end

    %% Plot ACF
    if plotResult
        % Plot rejection region lines for test of individual autocorrelations
        % H_0: rho(tau) = 0 at alpha=.05
        bar(ta)
        line([0 p+.5], (1.96)*(1/sqrt(Ny))*ones(1,2))
        line([0 p+.5], (-1.96)*(1/sqrt(Ny))*ones(1,2))

        % Some figure properties
        line_hi = (1.96)*(1/sqrt(Ny))+.05;
        line_lo = -(1.96)*(1/sqrt(Ny))-.05;
        bar_hi = max(ta)+.05 ;
        bar_lo = -max(ta)-.05 ;

        if (abs(line_hi) > abs(bar_hi)) % if rejection lines might not appear on graph
            axis([0 p+.60 line_lo line_hi])
        else
            axis([0 p+.60 bar_lo bar_hi])
        end
        title({' ','Sample Autocorrelations',' '})
        xlabel('Lag Length')
        set(gca,'YTick',[-1:.20:1])
        % set number of lag labels shown
        if (p<28 && p>4)
            set(gca,'XTick',floor(linspace(1,p,4)))
        elseif (p>=28)
            set(gca,'XTick',floor(linspace(1,p,8)))
        end
        set(gca,'TickLength',[0 0])
    end
end

% ---------------
% SUB FUNCTION
% ---------------
function ta2 = acf_Vk(y,k,ybar, Ny)
% ACF_K - Autocorrelation at Lag k
% acf(y,k)
%
% Inputs:
% y - series to compute acf for
% k - which lag to compute acf
% 
cross_sum = zeros(Ny-k,1) ;


% Numerator, unscaled covariance
for i = (k+1):Ny
    error('do something like this here:')
    sum(sum(M1(i,:) - ybar).*(M1(i-k, :) - ybar))
    cross_sum(i) = (y(i)-ybar)*(y(i-k)-ybar) ;
end

% Denominator, unscaled variance
yvar = (y-ybar)'*(y-ybar) ;

ta2 = sum(cross_sum) / yvar ;

