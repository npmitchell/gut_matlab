function plotHelmHodgeDECPullback(im, cutMesh, divs, rots, Options)
% PLOTHELMHODGEDECPULLBACK(im, Options)
%
% Parameters
% ----------
%
% 
% Returns
% -------
%
% NPMitchell 2020

% Unpack the cutMesh
FF = cutMesh ;
V2D = cutMesh.u ;
v3drs = cutMesh.v ;

% compute COM for each triangle in 2d and 3d --> note that we do this
% before gluing so that the 2d case is not messed up. The 3d case is
% the same either way.
bc = cat( 3, v3drs(FF(:,1), :), v3drs(FF(:,2), :), v3drs(FF(:,3), :) );
bc = mean( bc, 3 ) ;
bc2d = cat( 3, V2D(FF(:,1), :), V2D(FF(:,2), :), V2D(FF(:,3), :) );
bc2d = mean( bc2d, 3 ) ;

% Unpack divs and rots
divv = divs.divv ;
divU = divs.divU ;
divU2d = divs.divU2d ;
rotv = rots.rotv ;
rotU = rots.rotU ;
rotU2d = rots.rotU2d ;

addTitleStr = '' ;
labelUnit = '[1/min]' ;
harmLabelUnit = '[$\mu$m/min]' ;
xyzlim = [] ;
% Unpack options from struct
if isfield(Options, 'addTitleStr')
    addTitleStr = Options.addTitleStr ;
end
if isfield(Options, 'labelUnit')
    labelUnit = Options.labelUnit ;
end
if isfield(Options, 'harmLabelUnit')
    harmLabelUnit = Options.harmLabelUnit ;
end


% Save divergence and curl images
titlestrs = {[ 'dilatational flow, $\nabla \cdot v_t$:  ' addTitleStr], ...
    [ 'rotational flow, $\star \mathrm{d} v_t^\flat$:  ' addTitleStr], ...
    [ 'harmonic component of $v_t$:  ' addTitleStr]} ;
labelstrs = {['$\nabla \cdot v_t$, ' labelUnit], ...
    ['$\star$d$v_t^\flat$, ' labelUnit], ...
    ['harm$(v_t)$ ' harmLabelUnit]} ;
fnstrs2d = {div2dfn, curl2dfn} ;
fnstrs3d = {div3dfn, curl3dfn} ;
sfs = {divv, rotv} ;
vf3ds = {divU, rotU} ;
vf2ds = {divU2d, rotU2d} ;
for divcurl = 1:2 
    opts.style = 'diverging' ;
    opts.xlabel = 'AP position [\mum]' ;
    opts.ylabel = 'lateral position [\mum]' ;
    opts.zlabel = 'DV position [\mum]' ;
    if ~isempty(xyzlim)
        opts.xlim = xyzlim(1, :) ;
        opts.ylim = xyzlim(2, :) ;
        opts.zlim = xyzlim(3, :) ;
    end
    opts.titlestr = titlestrs{divcurl} ;
    opts.view = [0, 0] ;
    opts.label = labelstrs{divcurl} ;
    opts.outfn = fnstrs3d{divcurl} ;
    opts.sscale = 0.5 ;
    inds = (0:qsubU:(nU-1))' * 2 * (nU-1) * ones(length(1:qsubV:(nV-1)), 1)' +...
         ones(length(1:qsubU:(nU-1)), 1) * 2 * (1:qsubV:(nV-1)) ;
    vf = vf3ds{divcurl} ;
    scalarVectorFieldsOnSurface(FF, v3drs, sfs{divcurl}, ...
                bc(inds,1), bc(inds,2), bc(inds,3), ...
                vf(inds, 1), vf(inds, 2), vf(inds, 3), opts) ;
    close all    

    % Plot the 2d flow field divergence part
    opts2d.outfn = fnstrs2d{divcurl} ;
    opts2d.qsubsample = 1 ;
    opts2d.faces = FF ;
    opts2d.label = labelstrs{divcurl}  ;
    opts2d.sscale = 0.5 ;
    opts2d.style = 'diverging' ;
    opts2d.title = titlestrs{divcurl} ;
    vf = vf2ds{divcurl} ;
    if divcurl == 1
        % DIVERGENCE COMPONENT
        sf = reshape(sfs{divcurl}, [nU, nV])' ;
        xxs = V2D(1:nU, 1) ;
        yys = V2D(1:nU:nU*nV, 2) ;
    else
        % ROTATIONAL COMPONENT
        sf = sfs{divcurl} ;
        xxs = V2D(:, 1) ;
        yys = V2D(:, 2) ;
    end
    scalarVectorFieldsOnImage(im, xxs, yys, sf,...
        bc2d(inds, 1), bc2d(inds, 2), vf(inds, 1), vf(inds, 2), opts2d) ;
end

% % Plot the harmonic portion of the Helmholtz-Hodge decomposition
% % 3D harmonic field
% opts.style = 'phase' ;
% opts.xlabel = 'AP position [\mum]' ;
% opts.ylabel = 'lateral position [\mum]' ;
% opts.zlabel = 'DV position [\mum]' ;
% opts.xlim = xyzlim(1, :) ;
% opts.ylim = xyzlim(2, :) ;
% opts.zlim = xyzlim(3, :) ;
% opts.titlestr = titlestrs{3} ;
% opts.view = [0, 0] ;
% opts.label = labelstrs{3} ;
% opts.outfn = harm3dfn ;
% opts.sscale = 5 ;
% inds = (0:qsubU:(nU-1))' * 2 * (nU-1) * ones(length(1:qsubV:(nV-1)), 1)' +...
%      ones(length(1:qsubU:(nU-1)), 1) * 2 * (1:qsubV:(nV-1)) ;
% harmvphase = atan2(harmU2d(:, 2), harmU2d(:, 1)) ;
% scalarVectorFieldsOnSurface(FF, v3drs, harmvphase, ...
%             bc(inds,1), bc(inds,2), bc(inds,3), ...
%             harmU(inds, 1), harmU(inds, 2), harmU(inds, 3), opts) ;
% close all    
% 
% % 2D harmonic field
% options.outfn = harm2dfn ;
% vectorFieldHeatPhaseOnImage(im, V2D(:, 1), V2D(:, 2),...
%     harmU2d(:, 1), harmU2d(:, 2), 5, options)


% GPtoolbox for laplacian smoothing on VERTICES for div
%  laplacian_smooth(V,F,L_method,b,lambda,method,S,max_iter)