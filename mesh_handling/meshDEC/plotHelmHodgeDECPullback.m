function plotHelmHodgeDECPullback(im, cutMesh, vf, xy, vf2d, divs, rots, Options, opts2d)
% PLOTHELMHODGEDECPULLBACK(im, Options)
%
% Parameters
% ----------
% im : 2D image pullback
% cutMesh : struct with fields 
%   u : #vertices x 2 float array
%       2D vertices of mesh pullback
%   v : #vertices x 3 float array
%       3D vertices of mesh embedding
%   f : #faces x 3 int
%       face connectivity list
% vf : #faces x 3 float array 
%   original vector field before decomposition, in 3d
% vf2d : #faces x 3 float array 
%   original tangential vector field before decomposition, in 2d, properly
%   scaled by dilation (or not) for visualization
% divs : struct with fields
%   divv 
%       scalar divergence field values (div) defined on vertices
%   divU
%       vector dilatational field flow (curl-free) defined on faces
% rots : struct with fields 
%   rotv 
%       scalar rotational field values (curl) defined on faces
%   rotU
%       vector rotational field flow (div-free) defined on faces
% Options : struct with fields
%   div2dfn : str
%       output filename for dilatational flow in 2d
%   div3dfn : str
%       output filename for dilatational flow in 3d
%   rot2dfn : str
%       output filename for rotational flow in 2d
%   rot3dfn : str
%       output filename for rotational flow in 3d
%   xyzlim : 3 x 2 float array
%       limits in each dimension for 3d plot
%   qsubU : int (optional)
%       subsampling factor for quiver overlay in pullback x coord
%   qsubV : int (optional)
%       subsampling factor for quiver overlay in pullback y coord
%   sscaleDiv : float (optional)
%       scalar field scale, so color limit is (-sscaleDiv, +sscaleDiv)
%   sscaleRot : float (optional)
%       scalar field scale, so color limit is (-sscaleRot, +sscaleRot)
%   addTitleStr : str (optional)
%       addition to title (for ex, timestamp specifier)
%   
% 
% Returns
% -------
%
% NPMitchell 2020

% Unpack the cutMesh
FF = cutMesh.f ;
V2D = cutMesh.u ;
v3drs = cutMesh.v ;
nU = cutMesh.nU ;
nV = cutMesh.nV ;

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

% Unpack required Options
div2dfn = Options.div2dfn ;
div3dfn = Options.div3dfn ;
rot2dfn = Options.rot2dfn ;
rot3dfn = Options.rot3dfn ;

% Unpack other options
qsubU = 5 ;
qsubV = 10 ;
sscaleDiv = 0.5 ;
sscaleRot = 0.2 ;
qscaleDiv = 10 ;
qscaleRot = 25 ;
addTitleStr = '' ;
labelUnit = '[1/min]' ;
harmLabelUnit = '[$\mu$m/min]' ;
xyzlim = [] ;
% Unpack options from struct
if isfield(Options, 'sscaleDiv')
    sscaleDiv = Options.sscaleDiv ;
end
if isfield(Options, 'sscaleRot')
    sscaleRot = Options.sscaleRot ;
end
if isfield(Options, 'qscaleDiv')
    qscaleDiv = Options.qscaleDiv ;
end
if isfield(Options, 'qscaleRot')
    qscaleRot = Options.qscaleRot ;
end
if isfield(Options, 'addTitleStr')
    addTitleStr = Options.addTitleStr ;
end
if isfield(Options, 'labelUnit')
    labelUnit = Options.labelUnit ;
end
if isfield(Options, 'harmLabelUnit')
    harmLabelUnit = Options.harmLabelUnit ;
end
if isfield(Options, 'qsubU')
    qsubU = Options.qsubU ;
end
if isfield(Options, 'qsubV')
    qsubV = Options.qsubV ;
end
if isfield(Options, 'xyzlim')
    xyzlim = Options.xyzlim ;
end

% Pack scales for scalar magnitude and quiver length
sscales = [sscaleDiv, sscaleRot ] ;
qscales = [qscaleDiv, qscaleRot ] ;


% Save divergence and curl images
titlestrs = {[ 'dilatational flow, $\nabla \cdot v_t$:  ' addTitleStr], ...
    [ 'rotational flow, $\star \mathrm{d} v_t^\flat$:  ' addTitleStr], ...
    [ 'harmonic component of $v_t$:  ' addTitleStr]} ;
labelstrs = {['$\nabla \cdot v_t$, ' labelUnit], ...
    ['$\star$d$v_t^\flat$, ' labelUnit], ...
    ['harm$(v_t)$ ' harmLabelUnit]} ;
fnstrs2d = {div2dfn, rot2dfn} ;
fnstrs3d = {div3dfn, rot3dfn} ;
sfs = {divv, rotv} ;
% vf3ds = {divU, rotU} ;
% vf2ds = {divU2d, rotU2d} ;
vf2ds = {vf2d, vf2d } ;
for divcurl = 1:2 
    % 3D mesh plot with scalar and vector fields
    opts.style = 'diverging' ;
    opts.xlabel = 'AP position [\mum]' ;
    opts.ylabel = 'lateral position [\mum]' ;
    opts.zlabel = 'DV position [\mum]' ;
    if ~isempty(xyzlim)
        opts.xlim = xyzlim(1, :) ;
        opts.ylim = xyzlim(2, :) ;
        opts.zlim = xyzlim(3, :) ;
    end
    opts.view = [0, 0] ;
    opts.title = titlestrs{divcurl} ;
    opts.label = labelstrs{divcurl} ;
    opts.outfn = fnstrs3d{divcurl} ;
    opts.sscale = sscales(divcurl) ;
    opts.qscale = 3 ; % qscales(divcurl) ;
    opts.linewidth = 0.8 ;
    inds = (0:qsubU:(nU-1))' * 2 * (nU-1) * ones(length(1:qsubV:(nV-1)), 1)' +...
         ones(length(1:qsubU:(nU-1)), 1) * 2 * (1:qsubV:(nV-1)) ;
    % vfq = vf3ds{divcurl} ;
    % Plot the fields in 3d
    scalarVectorFieldsOnSurface(FF, v3drs, sfs{divcurl}, ...
                bc(inds,1), bc(inds,2), bc(inds,3), ...
                vf(inds, 1), vf(inds, 2), vf(inds, 3), opts) ;
    close all    

    % Plot the 2d flow field divergence part
    opts2d.outfn = fnstrs2d{divcurl} ;
    opts2d.qsubsample = 1 ;
    opts2d.faces = FF ;
    opts2d.title = titlestrs{divcurl} ;
    opts2d.label = labelstrs{divcurl}  ;
    opts2d.sscale = sscales(divcurl) ;
    opts2d.qscale = 20 ; % qscales(divcurl) ;
    opts2d.style = 'diverging' ;
    vfq = vf2ds{divcurl} ;
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
    
    % Option 1
    % scalarVectorFieldsOnImage(im, xxs, yys, sf,...
    %     bc2d(inds, 1), bc2d(inds, 2), vfq(inds, 1), vfq(inds, 2), opts2d) ;
    
    % Option 2
    % qsub = 20 ; % max(qsubU, qsubV) ;
    % xx = xy{1} ;
    % yy = xy{2} ;
    % xarr = xx(1, :)' ;
    % yarr = yy(:, 1) ;
    % vfqx = reshape(vfq(:, 1), length(xarr), length(yarr)) ;
    % vfqy = reshape(vfq(:, 2), length(xarr), length(yarr)) ;
    % % Now coarse-grain
    % xs = xarr(1:qsub:end, 1:qsub:end) ;
    % ys = yarr(1:qsub:end, 1:qsub:end) ;
    % xq = imresize(xx, [length(xs), length(ys)]) ;
    % yq = imresize(yy, [length(xs), length(ys)]) ;
    % vfqx = imresize(vfqx, [length(xs), length(ys)]) ;
    % vfqy = imresize(vfqy, [length(xs), length(ys)]) ;
    % scalarVectorFieldsOnImage(im, xxs, yys, sf,...
    %     xq, yq, vfqx, vfqy, opts2d) ;
    
    % Option 3
    scalarFieldOnImage(im, xxs, yys, sf, 0.6, sscales(divcurl), ...
        labelstrs{divcurl})
    title(labelstrs{divcurl})
    ylim([0.25 * size(im, 2), 0.75 * size(im, 2)])
    saveas(gcf, fnstrs2d{divcurl})
    
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