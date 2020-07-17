function measureStressPattern(QS, options)
%measureStressPattern(QS, options)
%   Measure required stress pattern from Helm-Hodge decomp potential fields
%
% Parameters
% ----------
%
%
% Returns
% -------
%
% NPMitchell 2020

%% IN PLANE STRESS PATTERN
% sigma = (d^2 phi / da db)  
%          + 1/2 [ (d^2 psi / da* db) + (d^2 psi / da* db)] 
%          + harmonic derivs


