function fitPhiOffsetsViaTexture(QS)
%
%
%

imfn_sp_prev = sprintf( QS.fullFileBase.im_sp, ...
                QS.xp.fileMeta.timePoints(tidx-1) ) ;

% Load 3D data for coloring mesh pullback
if isempty(QS.currentData.IV)
    QS.xp.loadTime(tt);
    QS.xp.rescaleStackToUnitAspect();

    % Raw stack data
    IV = QS.xp.stack.image.apply();
    IV = imadjustn(IV{1});         
    QS.currentData.IV = IV ;
else
    % Raw stack data
    IV = QS.currentData.IV ;
end

% Texture patch options
Options.PSize = 5;
Options.EdgeColor = 'none';
% Texture image options
Options.imSize = ceil( 1000 .* [ 1 a_fixed ] );
Options.yLim = [0 1];

% fit the shifts in the y direction
% todo: save uncorrected patchIms,
% could try tiling twice...
dmyk = 0 ;
phi0_fit = zeros(size(uspace)) ;
phi0s = zeros(size(uspace)) ;
phi0_fit_kk = 1 ; % for first pass                
phiv_kk = (vspace .* ones(nU, nV))' ;
ensureDir([sphiDir, '/phi0_correction/'])
while any(phi0_fit_kk > 0.002) && dmyk < 6
    disp(['Iteration ' num2str(dmyk)])
    plotfn = sprintf(phi0fitBase, tt, dmyk);
    patchImFn = sprintf( ...
        fullfile([sphiDir, '/phi0_correction/', fileNameBase, '_prephi0_' num2str(dmyk) '.tif']), ...
        QS.xp.fileMeta.timePoints(tidx-1) )  ;
    [phi0_fit_kk, phi0s_kk] = fitPhiOffsetsFromPrevPullback(IV, ...
        new3d, cutMesh.umax, uspace * cutMesh.umax, phiv_kk, ...
        ringpath_ss, imfn_sp_prev, lowerboundy, upperboundy, ...
        save_ims, plotfn, Options, ...
        step_phi0tile, width_phi0tile, potential_sigmay, 'integer', ...
        patchImFn) ;

    % Update the result: phi0_fit is what MATTERS
    dmyk = dmyk + 1;
    phi0_fit = phi0_fit + phi0_fit_kk ;
    phi0s = phi0s + phi0s_kk ;
    % phiv_kk is this iteration's correction only
    phiv_kk = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;
end