function aux_nes_beltrami(QS, FF, VV, refMesh, capID, nU, nV, tt, bfn, bifn, xyzlim)
% aux_nes_beltrami(QS, FF, VV, capID, nU, nV, tt, bfn, bifn)
% Save and plot beltrami coefficients for NES simulation timestep 
%
% Parameters
% ----------
% 
%
% NPMitchell 2021

% [~,~,~,xyzlim] = QS.getXYZLims() ;

[ff, vv] = remove_vertex_from_mesh(FF, VV, capID) ;
mesh_ii = struct('f', ff, 'v', vv, 'nU', nU, 'nV', nV) ;
cutOpts.ignoreRectilinearConstraint = true ;
meshCut = cutRectilinearCylMesh(mesh_ii, cutOpts) ;
mu = bc_metric(meshCut.f, refMesh.u_ricci, meshCut.v, 3) ;
m2d = meshCut ;

m2d.v = refMesh.u ;
m2d.v(:, 1) = m2d.v(:, 1) / max(m2d.v(:, 1)) ;

%% Save Beltrami
save(bfn, 'mu') ;

%% Plot Beltrami
% Move mu to vertices as a means of smoothing over adjacent triangles
[~, F2V] = meshAveragingOperators(mesh_ii) ;
mu3d = F2V * mu ;
mu2d = mu3d ;
mu2d(nU*(nV-1)+1:nU*nV) = mu2d(1:nU) ;
close all
opts = struct() ;
opts.clim = [-0.1, 0.1] ;
opts.alpha = 1.0 ;
opts.makeCbar = {false, false, true, true} ;
opts.labels = {'$\Re\mu$', '$\Im\mu$', '', ''} ;
opts.cbarlabels = {'', '', '$\Re\mu$', '$\Im\mu$'} ;
opts.xzylims = {xyzlim, xyzlim, [], []} ;
opts.views = {[0,0], [0,0], [0,90], [0,90]} ;
opts.axisOff = true ;
opts.cmap = bwr ;
opts.visible = 'off' ;
[axs, cbs, meshHandles] = nFieldsOnSurface({mesh_ii,mesh_ii,m2d,m2d}, ...
    {real(mu3d), imag(mu3d),real(mu2d), imag(mu2d)}, opts) ;
expandSecondAxesRow(axs, -0.1)
try
    sgtitle(['$t=$' num2str(tt * QS.timeInterval) ' ' QS.timeUnits], ...
        'Interpreter', 'latex')
catch
    title(['$t=$' num2str(tt * QS.timeInterval) ' ' QS.timeUnits], ...
        'Interpreter', 'latex')
end
        
saveas(gcf, bifn)
close all   