function [phi0_fit, phi0s] = fitPhiOffsetsFromPrevMesh(TF, TV2D, TV3D,...
    uspace, vvals, prev3d_sphi, lowerbound_phi0, upperbound_phi0, ...
    save_fit, plotfn, save_phi0patch, varargin)
%FITPHIOFFSETSFROMPREVMESH(TF, TV2D, TV3D, uspace, vspace, prev3d_sphi, lowerbound, upperbound, save_im, plotfn) 
%   Fit the offset phi values to add to V in UV coords to minimize
%   difference in 3D between current embedding mesh and previous one. This
%   rotates the hoops of the sphicutMesh.
%
% Parameters
% ----------
% TF : #faces x 3 int array
%   connectivity list of the mesh vertices into triangulated faces
% TV2D : #vertices x 2 float array
%   the UV mapped coordinates of mesh vertices
% TV3D : #vertices x 3 float array
%   the embedding coordinates of mesh vertices
% uspace : nU float array
%   The values of u for each line of constant v in pullback space
% vvals : nV float array OR nU x nV float array as grid
%   If nV x 1 float array, the values of v for each line of constant u in 
%   pullback space, otherwise the values for the whole grid. This allows
%   you to pass either a linspace for v (independent of u) or a series of
%   v values, one array for each u value
% prev3d_sphi : nU x nV x 3 float array
%   The 3D coordinates of the embedding for the reference timepoint
%   (previous timepoint, for ex) at the 2D locations given by uspace and 
%   vspace. Note that uspace is not used explicitly, only nU is used to 
%   extract the strips over which we iterate, minimizing for phi0 for each
%   strip.
% lowerbound_phi0 : float 
%   lower bound for the fit of phi (offset to v). Must be > -1 and < 1
% upperbound_phi0 : float
%   upper bound for the fit of phi (offset to v). Must be > -1 and < 1.
%   Also must be <= lowerbound_phi0 + 1 to avoid phase ambiguity
% save_im : bool
% plotfn : str
% save_phi0patch : bool
%   save a patch colored by the phi0 motion deduced from the difference
%   between previous and current DVhoop coordinates
% varargin : optional positional arguments
%   preview   : bool (optional) visualize the progress of phi0
%   patchImFn : str  (optional) save the progress of phi0 as texture 
%                                 image saved to this path
%   imfn_sp_prev : str (optional) path to previous timepoint's sp pullback
%   IV        : MxNxP array (optional) intensity data to use for patch 
%
%
% Returns
% -------
% phi0_fit
% phi0s : nU x 1 float array 
%   the additional rotation angles, bounded by (lowerbound, upperbound) 
% 
% NPMitchell 2019

% If a preview input boolean is passed, interpret it
if ~isempty(varargin)
    preview = varargin{1} ;
    if length(varargin) > 1 && save_phi0patch
        patchOpts = varargin{2} ;
        patchImFn = patchOpts.patchImFn ;
        imfn_sp_prev = patchOpts.imfn_sp_prev ;
        IV = patchOpts.IV ;
        ringpath_ss = patchOpts.ringpath_ss ;
        v3d = patchOpts.v3d ;
        Options = patchOpts.Options ;
    else
        error('save_phi0patch is true, but patchOpts are missing!')
    end
else
    preview = false ;
    if save_phi0patch
        error('save_phi0patch is true, but patchOpts are missing!')
    end
end

% if ~isempty(varargin)
%     if isa(varargin{i},'logical')
%         continue;
%     end
% 
%     if ~isempty(regexp(varargin{i},'^[]ace[Nn]ormals','match'))
%         facenormals = varargin{i+1} ;
%     end
% end


% Consider each value of u in turn
% Fit for phi0 such that v = phi + phi0, meaning new
% phi = v - phi0.
disp('Minimizing phi0s...')

% Using simple offset
phi0s = phiOffsetsFromPrevMesh(TF, TV2D, TV3D, ...
    uspace, vvals, prev3d_sphi, lowerbound_phi0, upperbound_phi0, {preview}) ;

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
if save_fit
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


%% OPTIONAL: PLOT CHANGE IN PHI OVER TEXTURE PATCH IMAGES OF PREVIOUS AND CURRENT IMAGE
if save_phi0patch
    % First check that we have vgrid, not vspace 
    nU = length(uspace) ;
    if any(size(vvals) == 1)
        nV = length(vvals) ;
        vgrid = (vvals .* ones(nV, nU))' ;
    else
        vgrid = vvals ;
        nV = size(vgrid, 1) ;
    end
    onesUV = ones(nU, nV) ;
    uu = uspace .* onesUV ; % only used for faces definition
    
    % Load previous sp pullback image
    if exist(imfn_sp_prev, 'file')
        im0 = double(imread(imfn_sp_prev)) / 255.0 ;
    else
        error('Previous timepoint not available')
    end
    tmp = ringpath_ss .* onesUV ;
    svcutMesh.u(:, 1) = tmp(:) ;
    svcutMesh.u(:, 2) = vgrid(:) ;
    svcutMesh.v = v3d ;
    
    % Generate cutMesh face triangulation
    vvtmp = linspace(0, 1, nV) .* ones(nV, nU) ; % only used for faces definition
    svcutMesh.f = defineFacesRectilinearGrid([uu(:), vvtmp(:)], nU, nV) ;
    
    % Generate cutpath pairs
    svcutP1 = 1:nU ;
    svcutP2 = nU*nV - fliplr(0:(nU-1)) ;
    svcutMesh.pathPairs = [ svcutP1', svcutP2' ] ;

    % [phi0_fit, phi0s] = fitPhiOffsetsFromPrevPullback() ;

    %% Generate Tiled Orbifold Triangulation -------------------
    disp('fitPhiOffsetsFromPrevMesh: generating temporary pullback')
    tileCount = [1 1];  % how many above, how many below
    [ TF, TV2D, TV3D ] = tileAnnularCutMesh( svcutMesh, tileCount );

    % Create texture image
    if any(isnan(TV2D))
        error('here -- check for NaNs case')
    end
    patchIm = texture_patch_to_image( TF, TV2D, TF, TV3D(:, [2 1 3]), ...
        IV, Options );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save a preview of the intermediate frame with phi0pix as heatmap,
    % overlaying current and previous patchIm
    close all
    figure('visible', 'off')
    disp(['Saving patchIm to ' patchImFn])
    % convert (ss, phi) to (x, y)
    xx = ringpath_ss * size(patchIm, 2) / max(ringpath_ss) ;
    yy = 1:100:size(patchIm, 1) ;
    phi0grid = (phi0_fit .* ones(length(uspace), length(yy)))' ;
    tmp = cat(3, patchIm, patchIm, im0) ;
    opts.label = '$\phi_0$' ;
    opts.qsubsample = 2 ;
    opts.qscale = 1000 ;
    [~, ~, ~, ax, ~] = vectorFieldHeatPhaseOnImage(tmp, xx, yy', ...
        0*phi0grid, phi0grid, max(abs(phi0_fit))*2, opts) ;
    title(ax, 'Blue is prev timept, yellow is current')
    % F = getframe(gca);
    % Image = frame2im(F);
    % imwrite(Image, patchImFn)
    saveas(gcf, patchImFn) 
    close all
end

end

