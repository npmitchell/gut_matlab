function aux_nes_divcurl(QS, FF, VV, V1, climit_div, divcurlfn, tt)
% aux_nes_divcurl(QS, FF, VV, V1, climit_div, divcurlfn)
% 
% 
% NPMitchell 2021 

[~,~,~,xyzlim] = QS.getXYZLims() ;

disp('Creating figure with divergence and curl')
close all      
% Extract face barycenters for this timepoint in 3d
bc0 = barycenter(V1, FF) ;
% Extract for prev timpepoint in 3d
bc = barycenter(VV, FF) ;
v0 = bc - bc0 ;
v0_vtx = VV - V1 ;
% Now resolve the vector field for decomposition
[v0n, v0t ] = ...
    resolveTangentNormalVelocities(FF, V1, v0, 1:length(FF)) ;

% Construct DEC class instance
DEC = DiscreteExteriorCalculus( FF, V1 ) ;
% Compute divergence and rotational flow
divv = DEC.divergence(v0t) ;
rotv = DEC.curl(v0t) ;
rotv(isnan(rotv)) = 0 ;

% Plot div
close all
opts = struct() ;
opts.clims = {climit_div, climit_div} ;
opts.alpha = 1.0 ;
opts.labels = {'divergence', 'curl'} ;
opts.cbarlabels = {'$\star$d$\star\left(v_t^\flat\right)$', ...
    '$\star$d$v_t^\flat$'} ;
opts.xzylims = xyzlim ;
opts.view = [0,0] ;
opts.axisOff = true ;
opts.cmap = bwr ;
nFieldsOnSurface({FF, 0.5*(V1+VV)}, {divv, rotv}, opts) ;

sgtitle(['$t=$' num2str(tt * QS.timeInterval) ' ' QS.timeUnits], ...
    'Interpreter', 'latex')

saveas(gcf, divcurlfn)
close all    