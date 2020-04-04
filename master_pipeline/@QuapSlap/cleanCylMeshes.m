function cleanCylMeshes(QS, cleanCylOptions) 
% cleanCylMeshes(QS, cleanCylOptions)
%
% Parameters
% ----------
%
% Returns
% -------
%
%
% NPMitchell 2020

% Unpack options
overwrite = false ;
save_ims = true ;
if isfield(cleanCylOptions, 'overwrite')
    overwrite = cleanCylOptions.overwrite ;
end
if isfield(cleanCylOptions, 'save_ims')
    save_ims = cleanCylOptions.save_ims ;
end
cylinderMeshCleanBase = QS.fullFileBase.cylinderMeshClean ;
cylinderMeshBase = QS.fullFileBase.cylinderMesh ;
dpFile = QS.fileName.apBoundaryDorsalPts ;
outadIDxfn = QS.fileName.aBoundaryDorsalPtsClean ;
outpdIDxfn = QS.fileName.pBoundaryDorsalPtsClean ;
cylinderMeshCleanFigBase = fullfile( QS.dir.cylinderMeshClean, 'images', ...
    [ QS.fileBase.mesh 'cylindercut_clean.png' ]) ;

[~, ~, xyzlim_um] = QS.getXYZLims() ;

% Load or compute clean cylindrical mesh for each timepoint
for t = QS.xp.fileMeta.timePoints 
    disp(['t = ' num2str(t)])
    mesh3dfn =  sprintf( cylinderMeshCleanBase, t ) ;
    if ~exist(mesh3dfn, 'file') || overwrite
        disp('Overwriting/Computing clean cylinderMesh')
        % Load the cylinder mesh
        cylmeshfn = sprintf( cylinderMeshBase, t ) ;
        mesh = read_ply_mod( cylmeshfn );
        % Clean the mesh, including averaging normals onto vertices using
        % angle-weighted scheme
        mesh = cleanCylMesh(mesh) ;   

        % The dataset name base for the AD and PD points
        ADBase = ['/' QS.fileBase.mesh '/adorsal'];
        PDBase = ['/' QS.fileBase.mesh '/pdorsal'];
        [adIDx, pdIDx] = aux_adjust_dIDx(mesh, t, dpFile, ...
            ADBase, PDBase, cylinderMeshCleanBase, ...
            outadIDxfn, outpdIDxfn, QS.xp.fileMeta.timePoints) ;

        %% Save the 3d cut mesh with new indices
        % This is saving the cylinder meshes with no ears. Also adIDx
        % is saved in h5file.
        plywrite_with_normals(mesh3dfn, mesh.f, mesh.v, mesh.vn)
        % Save adIDx with new indices
        save_to_h5(outadIDxfn, ['/' sprintf('%06d', t) ], adIDx, ['adIDx for t=' num2str(t) ' already exists'])
        % Save pdIDx with new indices
        save_to_h5(outpdIDxfn, ['/' sprintf('%06d', t) ], pdIDx, ['pdIDx for t=' num2str(t) ' already exists'])
        disp('done with cylindermesh cleaning')

        % View results --------------------------------------------------------
        mesh3dfigfn = sprintf( cylinderMeshCleanFigBase, t ) ;
        if save_ims
            close all
            fig = figure('visible', 'off') ;
            xyzrs = QS.xyz2APDV(mesh.v) ;
            trisurf( triangulation( mesh.f, xyzrs ), xyzrs(:, 2), ...
                'EdgeColor', 'none', 'FaceAlpha', 0.5);
            hold on
            scatter3( xyzrs(adIDx,1), xyzrs(adIDx,2), xyzrs(adIDx,3), ...
                'filled', 'r' );
            scatter3( xyzrs(pdIDx,1), xyzrs(pdIDx,2), xyzrs(pdIDx,3), ...
                'filled', 'c' );
            hold off
            axis equal
            xlabel('x [\mum]')
            ylabel('y [\mum]')
            zlabel('z [\mum]')
            xlim(xyzlim_um(1, :)); 
            ylim(xyzlim_um(2, :)); 
            zlim(xyzlim_um(3, :));
            title(['cleaned cylinder mesh, t = ' num2str(t)])
            saveas(fig, mesh3dfigfn)
            close all
        end
    else
        mesh = read_ply_mod(mesh3dfn) ;
        adIDx = h5read(outadIDxfn, ['/' sprintf('%06d', t)]) ;
        pdIDx = h5read(outpdIDxfn, ['/' sprintf('%06d', t)]) ;

        % View results --------------------------------------------------------
        mesh3dfigfn = sprintf( cylinderMeshCleanFigBase, t ) ;
        if (~exist(mesh3dfigfn, 'file') || overwrite) && save_ims
            aux_plot_cleanCylMesh
        end
    end
end