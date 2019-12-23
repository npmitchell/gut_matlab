function [phi0_fit, phi0s] = fitPhiOffsetsFromPrevMesh(TF, TV2D, TV3D, uspace, vspace, prev3d_sphi, lowerbound, upperbound, save_im, plotfn)
%FITPHIOFFSETSFROMPREVMESH(TF, TV2D, TV3D, uspace, vspace, prev3d_sphi, lowerbound, upperbound, save_im, plotfn) 
%   Fit the offset phi values to add to V in UV coords to minimize
%   difference in 3D between current embedding mesh and previous one. This
%   rotates the hoops of the sphicutMesh.
%
% NPMitchell 2019

                       
% Consider each value of u in turn
% Fit for phi0 such that v = phi + phi0, meaning new
% phi = v - phi0.
disp('Minimizing phi0s...')

% Using simple offset
phi0s = phiOffsetsFromPrevMesh(TF, TV2D, TV3D, ...
    uspace, vspace, prev3d_sphi, lowerbound, upperbound, {false}) ;

% Using dilation and offset
% [phi0s, ccoeffs] = phiOffsetsFromPrevMeshWithDilation(TF, TV2D, TV3D, ...
%      uspace, vspace, prev3d_sphi, preview) ;
% phi0s = mod(phi0s, 2 * pi) ;

% Convert phi0s to a smooth polynomial phi(u)
% Smoothing parameters
framelen = 11 ;  % must be odd
polyorder = 2 ;
% Low pass filter (Savitsky-Golay)
% Note we ignore the variations in ds -->
% (instead use du=constant) to do this fit
phi0_fit = savgol(phi0s, polyorder, framelen)' ;

% Fit the smoothed curve
% phicoeffs = polyfit(uspace, phi0_fit, 14) ;
% phi0_fit = polyval(phicoeffs, uspace);

% Plot the fit for reference
if save_im
    close all
    fig = figure('visible', 'off') ;
    plot(uspace, phi0s, '.'); hold on;
    % plot(uspace, phix, '--')
    plot(uspace, phi0_fit, '-')
    legend({'measured shifts', 'SG filter'})
    xlabel('u')
    ylabel('\phi_0')
    title('Shift \phi(u)')
    saveas(fig, plotfn)
    close all
end

end

