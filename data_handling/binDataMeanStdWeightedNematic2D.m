function [meanQ_bins, stdQ_bins, meanQBinStrength, meanQBinTheta, ...
    meanQBinStrengthStd] = ...
    binDataMeanStdWeightedNematic2D(nn, strengths, xx, xedges, weights) 
%
%
% Parameters
% ----------
% nn: 
% strengths : 
% xedges :
% weights 
%
% Returns
% -------
% meanQ_bins 
% stdQ_bins 
% meanQBinAspect 
% meanQBinTheta
% meanQBinAspectStd
% 

nBins = numel(xedges) - 1 ;
nCells = size(nn, 1) ;
if length(size(nn)) == 3 && all(size(nn) > 1)
    QQ = nn ;
else
    QQ = zeros(nCells, 2, 2) ;
    for cc = 1:nCells
        try
            QQ(cc, :, :) = (nn(cc, :)' * nn(cc, :) - 0.5 * [1, 0; 0, 1]) ;
        catch
            error(['could not create Q tensor: n=' num2str(nn)])
        end
    end
end
[~, bins_Q11, bins_std_Q11] = binDataMeanStdWeighted(xx, ...
    strengths .* squeeze(QQ(:, 1, 1)), xedges, weights) ;
[~, bins_Q12, bins_std_Q12] = binDataMeanStdWeighted(xx, ...
    strengths .* squeeze(QQ(:, 1, 2)), xedges, weights) ;
[~, bins_Q21, bins_std_Q21] = binDataMeanStdWeighted(xx, ...
    strengths .* squeeze(QQ(:, 2, 1)), xedges, weights) ;
[~, bins_Q22, bins_std_Q22] = binDataMeanStdWeighted(xx, ...
    strengths .* squeeze(QQ(:, 2, 2)), xedges, weights) ;

% Check that result is still traceless and symmetric
assert(all(abs(bins_Q11 + bins_Q22) < 1e-7))
assert(all(abs(bins_Q12 - bins_Q12) < 1e-7))

% Collate bin information into eigenvalues & directions
if nargout > 2
    meanQBinStrength = zeros(nBins, 1) ;
    meanQBinTheta = zeros(nBins, 1) ;
    meanQBinStrengthStd = zeros(nBins, 1) ;
    for bin = 1:nBins
        meanQ_bins{bin} = [bins_Q11(bin), bins_Q12(bin); ...
            bins_Q21(bin), bins_Q22(bin)] ;
        stdQ_bins{bin} = [bins_std_Q11(bin), bins_std_Q12(bin); ...
            bins_std_Q21(bin), bins_std_Q22(bin)] ;

        % diagonalize this binQ
        [ eig_vec, eig_val ] = eig(meanQ_bins{bin});
        try
            assert(abs(abs(eig_val(2,2)) - abs(eig_val(1,1))) < 1e-7)
        catch
            error('Something is wrong with traceless or symmetry')
        end
        meanQBinStrength(bin) = norm(meanQ_bins{bin}) * 2 ;
        meanQBinTheta(bin) = atan2( eig_vec(2,2), eig_vec(1,2) );

        % Uncertainty in average is given by error propagation
        % lambda = 0.5 * [trace +/- sqrt(tr^2 - 4*det)] 
        % Now, the trace is guaranteed to be zero, but not sure that means
        % unc_trace = 0. If not, then we would propagate errors to be

        unc_tr = sqrt(2 * bins_std_Q11(bin).^2) ;
        unc_det = sqrt(2 * (bins_Q11(bin) * bins_std_Q11(bin)).^2 ...
            + 2 * (bins_Q12(bin) * bins_std_Q12(bin)).^2) ;
        determ = abs(det(meanQ_bins{bin})) ;
        unc_lambda = 0.5 * sqrt(unc_tr.^2 + unc_det.^2 / (determ)) ; 

        % NOTE: |eigenvalue| of symm traceless matrix == norm(matrix)
        meanQBinStrengthStd(bin) = 2 * unc_lambda ;
    end
end