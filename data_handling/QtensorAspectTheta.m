function [meanQLobeAspect, meanQLobeTheta, meanQLobeAspectStd, ...
     meanQLobeThetaStd] = QtensorAspectTheta(meanQ, uncQ, options)
% [meanQLobeAspect, meanQLobeTheta, meanQLobeAspectStd, ...
%     meanQLobeThetaStd] = QtensorAspectTheta(meanQ, uncQ, options)
%     
%
% Parameters
% ----------
% meanQ : 2x2 float
%   nematic tensor whose aspect ratio (strength) and angle to compute
% uncQ : 2x2 float
%   uncertainties in nematic tensor components whose aspect ratio 
%   (strength) and angle to compute, for computing uncertainty of strength
%   and angle.
% options : struct with fields
%
% Returns
% -------
% meanQLobeAspect : float
%   nematic tensor long axis strength (for unweighted Q, this is 0.5)
% meanQLobeTheta : float
%   nematic tensor long axis angle
% meanQLobeAspectStd : float
%   uncertainty in the nematic tensor strength
% meanQLobeThetaStd : float
%   uncertainty in the nematic tensor angle
%
% NPMitchell 2021


% Option parsing
traceless = true ;
if nargin > 2
    if isfield(options, 'traceless')
        traceless = options.traceless ;
    end
end

[ eig_vec, eig_val ] = eig(meanQ);
try
    assert(abs(abs(eig_val(2,2)) - abs(eig_val(1,1))) < 1e-7)
catch
    error('Something is wrong with traceless or symmetry')
end
meanQLobeAspect = norm(meanQ) * 2 + 1 ;
meanQLobeTheta = atan2( eig_vec(2,2), eig_vec(1,2) );


if nargout > 1
    Q11 = meanQ(1, 1) ;
    Q12 = meanQ(1, 2) ;
    Q21 = meanQ(2, 1) ;
    Q22 = meanQ(2, 2) ;
    uncQ11 = uncQ(1, 1) ;
    uncQ12 = uncQ(1, 2) ;
    % stdQ21 = stdQ(2, 1) ;  % = stdQ12, so redundant
    % stdQ22 = stdQ(2, 2) ;  % = stdQ11
    % steQ11 = steQ(1, 1) ;
    % steQ12 = steQ(1, 2) ;
    % steQ21 = steQ(2, 1) ;
    % steQ22 = steQ(2, 2) ;

    % Uncertainty in average is given by error propagation
    % lambda = 0.5 * [trace +/- sqrt(tr^2 - 4*det)] 
    % Now, the trace is guaranteed to be zero, but not sure that means
    % unc_trace = 0. If not, then we would propagate errors to be
    if traceless
        unc_tr = 0 ; 
    else
        unc_tr = sqrt(2 * uncQ11.^2) ;
    end
    unc_det = sqrt(2 * (Q11 * uncQ11).^2 ...
        + 2 * (Q12 * uncQ12).^2) ;
    determ = abs(det(meanQ)) ;
    unc_lambda = 0.5 * sqrt(unc_tr.^2 + unc_det.^2 / (determ)) ; 
    % NOTE: |eigenvalue| of symm traceless matrix == norm(matrix)
    meanQLobeAspectStd = 2 * unc_lambda ;


    % % Similar error propagation for the ste ---------------------------
    % unce_tr = sqrt(2 * steQ11.^2) ;
    % unce_det = sqrt(2 * (Q11 * steQ11).^2 ...
    %     + 2 * (Q12 * steQ12).^2) ;
    % determ = abs(det(meanQ)) ;
    % unce_lambda = 0.5 * sqrt(unce_tr.^2 + unce_det.^2 / (determ)) ;
    % % NOTE: |eigenvalue| of symm traceless matrix == norm(matrix)
    % meanQLobeAspectSte = 2 * unce_lambda ; 

    % For angle uncertainty, note that Q is traceless symmetric so
    % we can say Q = [A,B;B,-A]. 
    % Then the eigenvalues are lambda = +/- sqrt(A^2+B^2)
    % Plugging back in allows us to find the eigvects as theta(A,B)
    % so that we can get dtheta(A,B,dA,dB).
    % dtheta = Sqrt[D[theta,A]^2 dA^2 + D[theta,B]^2 dB^2]
    %        = 0.5 * sqrt[ (B^2 dA^2 + A^2 dB^2) / (A^2+B^2)^2 ]
    assert(abs(Q11 + Q22)<1e-8)
    assert(abs(Q21 - Q12)<1e-8)
    numerator1 = Q11^2 * (uncQ12)^2 ;
    numerator2 = Q12^2 * (uncQ11)^2 ;
    numerator = numerator1+numerator2 ;
    denominator = (Q11^2 + Q12^2)^2 ;
    meanQLobeThetaStd = 0.5 * sqrt(numerator / denominator) ;

    % % Similar error propagation for standard error
    % numerator1 = Q11^2 * (steQ12)^2 ;
    % numerator2 = Q12^2 * (steQ11)^2 ;
    % numerator = numerator1+numerator2 ;
    % denominator = (Q11^2 + Q12^2)^2 ;
    % meanQLobeThetaSte = 0.5 * sqrt(numerator / denominator) ;
end

