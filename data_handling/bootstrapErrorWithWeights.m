function [errMean, fitResult ] = bootstrapErrorWithWeights(xx, weights, options) 
%bootstrapErrorWithWeights(xx, weights, options) 
% Estimate standard error of the mean of observations xx via bootstrapping. 
% In default case ('full'), compute variance of mean values of 
%   submsampled batches of xx. Fit variance(n) for several values
%   of n<N, where N = length(xx). Then plug in N to the fit to learn error
%   on the mean. 
% In simple case ('simple'), subsample nsamples times at n=N and report ste
%
% Parameters 
% ----------
% xx : NxD floats (todo: check that this works for D>1)
%   observation data whose standard error we wish to compute
% weights : N floats
%   weights to ascribe to each observation, could be inverse uncertainty,
%   for ex
% options : optional struct with fields
%   ssf : M x 1 floats between (0,1] 
%       subsampling fractions to query for the fit of linear polynomial
%       variance(((1/n)-mu(1))/mu(2), b) = m*xhat + b, where
%       XHAT = (X-MU(1))/MU(2) where MU(1) = MEAN(X) and MU(2) = STD(X)
%   nsamples : int
%       how many times to subsample with each ssf(i)
%   preview : bool
%       preview the fit at the end
%   weightMeanOrSampling : str ('mean' or 'sampling')
%       whether to utilize weights in biasing the choosing of n 
%       observations (sampling) or whether to choose n samples with uniform
%       probability, then take weighted mean using weights (mean).
%   simpleOrFull : str ('simple' or 'full')
%       'simple' : subsample at N only, report std(means) for samples
%       'full'   : subsample at ssf*N, fit variance to 1/n and evaluate fit
%
% Returns
% -------
% errMean : float
%   estimate for the error on the mean of the weighted observations
% fitResult : struct with fields
%   p : 2 floats
%       flope and intercept of fit of variance = m*mu + b ;
%   S : fit quality
%   mu : 2 floats
%       XHAT = (X-MU(1))/MU(2) where MU(1) = MEAN(X) and MU(2) = STD(X)
%   readme : str
%       description of results
%
% NPMitchell 2021

ssf = linspace(0.3, 1, 20) ;
nsamples = 500 ;
preview = false ;
weightMeanOrSampling = 'mean' ;
simpleOrFull = 'full' ;

if nargin < 3
    options = struct() ;
else
    if isfield(options, 'subsamplingFractions')
        ssf = options.subsamplingFractions ;
    end
    if isfield(options, 'nsamples')
        nsamples = options.nsamples ;
    end
    if isfield(options, 'weightMeanOrSampling')
        weightMeanOrSampling = options.weightMeanOrSampling ;
    end
    if isfield(options, 'preview')
        preview = options.preview ;
    end
end

% renormalize weights
ww = weights / nansum(weights) ;

if strcmpi(simpleOrFull, 'full')
    nx = length(xx) ;
    variances = zeros(length(ssf),1 ) ;
    kks = zeros(length(ssf),1) ;
    for qq = 1:length(ssf)
        % subsample how many items?
        kk = round(nx* ssf(qq)) ;
        kks(qq) = kk ;

        % get distribution of means 
        meany = zeros(nsamples, 1) ;
        for pp = 1:nsamples
            if strcmpi(weightMeanOrSampling, 'sampling')
                % subsample from xx
                yy = datasample(xx, kk, 'weights', ww) ;
                meany(pp) = mean(yy) ;
            else
                % Alternatively, weight in sum:
                [yy, idx] = datasample(xx, kk) ;
                % renormalize weights
                wk = weights(idx) / nansum(weights(idx)) ;
                meany(pp) = nansum(wk .* yy) ;
            end
        end

        % Look at variance of means
        variances(qq) = var(meany) ;
    end

    % Fit to a line
    % [P,S,MU] = POLYFIT(X,Y,N) finds the coefficients of a polynomial in
    %   XHAT = (X-MU(1))/MU(2) where MU(1) = MEAN(X) and MU(2) = STD(X)
    [pp, SS, mu] = polyfit(kks.^(-1), variances, 1) ;
    errMean = sqrt(polyval(pp,1/nx,SS,mu)) ;

    if nargout > 1
        fitResult = struct('p', pp, 'S', SS, 'mu', mu, ...
            'readme', ['fit to subsamples: variance_mu = p(1)/k + p(2) and '...
            'errMean=sqrt(polyval(pp,1/nx,SS,mu)), with linear fit to '...
            'XHAT=((1/n)-MU(1))/MU(2) where MU(1) = MEAN(X) and MU(2) = STD(X)']) ;
    end

    % plot the result if preview==true
    if preview
        clf
        varvals = polyval(pp,[kks.^(-1); 1/nx],SS,mu) ;
        plot(kks.^(-1), variances, '.'); hold on;
        plot([kks.^(-1); 1/nx], varvals, '--') ;
        scatter(1/nx, errMean, 100, 'filled')
        xlabel('inverse subsampling size, $1/n$', 'interpreter', 'latex')
        ylabel('variance of subsample mean, $\sigma^2_{\mu}$', 'interpreter', 'latex')
        title('Bootstrapping fit', 'interpreter', 'latex')
    end
elseif strcmpi(simpleOrFull, 'simple')
    % simply subsample at n=N nsamples times
    bootstatS = bootstrp(nsamples,@mean,xx,'weights', ww);
    errMean = std(bootstatS) ;
    
    if nargout > 1
        error(['Too many output args: no fitting performed for ', ...
            'simple method of bootstrap ste'])
    end
else
    error('could not interpret options.simpleOrFull')
    
end



