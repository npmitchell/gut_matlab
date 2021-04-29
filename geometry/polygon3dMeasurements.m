function [c3d, cellCntrd3d, areas, perim, moment1, ang1, ...
    moment2, ang2, moinertia, cellQ2d, cellMeshFaces, vertexMeshFaces] = ...
    polygon3dMeasurements(faces, v3D, v2D, polygonsXY, cntrds)
% [c3d, areas, perim, moment1, ang1, ...
%     moment2, ang2, moinertia, cellQ2d, vertexMeshFaces] = ...
%     polygon3dMeasurements(faces, v3D, v2D, cellVtx2D, polygons)
% 
% Parameters
% ----------
% faces : #faces x 3 int array
%   indices into vertex arrays of each triangular face
% v3D : #vertices x 3 float array
%   pushforward vertices in which polygons/segmentation are mapped
% v2D : #vertices x 2 float array
%   pullback vertices in which polygons/segmentation are known
% polygonsXY : #polygons x 1 cell array of numeric 2d arrays
%   polygonsXY{i} are the vertices of polygon i in 2d (XY)
% cntrds : #polygons x 2 numeric array
%   positions where cells are considered "centered" in 2D mesh. May be
%   centroids of polygonsXY
%
% Returns
% -------
% centroids2d = '2d evaluation coordinates for tangent plane on embedding -- could/should be the 2d polygon centroids' ;
% centroids3d = 'centroids2d in push-forward / embedding space' ;
% areas = 'area of each polygon in tangent plane';
% perim = 'perimeter of each polygon in tangent plane';
% moment1 = 'maximum moment of inertia (extremal #1)' ;
% ang1 = 'angle of eigenvector of moment of inertia tensor of polygon in 3d in a frame rotated to the 2d pullback coordinate frame' ;
% moment2 = 'minimum moment of inertia (extremal #2)';
% moinertia = 'embeddingspace moment of inertia of pushed forward polygon in basis aligned with pullback' ;
% cellQ2d = 'polygon as quasi-2d' ;
% fieldfaces = 'indices of faces on 3d mesh of evaluation points (centroids)';
%
% See also
% --------
% polygonNetwork3dMeasurements()
%
% Example Usage
% --------------
% polygon3dMeasurements(faces, v3D, v2D, polygons2D) ;


% polygonsXY should be an Nx1 or 1xN cell array
assert(any(size(polygonsXY) == 1) && numel(size(polygonsXY))==2) 
nCells = length(polygonsXY) ;

% cell vertices in 3d
c3d = cell(nCells, 1) ;
vertexMeshFaces = cell(nCells, 1) ;
for pp = 1:nCells
    if ~isempty(polygonsXY{pp})
        [c3d{pp}, vertexMeshFaces{pp}] = ...
            interpolate2Dpts_3Dmesh(faces, v2D, v3D, polygonsXY{pp}) ;
        if any(isnan(c3d{pp}))
            % fix any nans with simple scattered interpolation
            Xai = scatteredInterpolant(v2D(:, 1), v2D(:, 2), v3D(:, 1), 'linear', 'nearest') ;
            Yai = scatteredInterpolant(v2D(:, 1), v2D(:, 2), v3D(:, 2), 'linear', 'nearest') ;
            Zai = scatteredInterpolant(v2D(:, 1), v2D(:, 2), v3D(:, 3), 'linear', 'nearest') ;
            xinterp = Xai(polygonsXY{pp}(:, 1), polygonsXY{pp}(:, 2)) ;
            yinterp = Yai(polygonsXY{pp}(:, 1), polygonsXY{pp}(:, 2)) ;
            zinterp = Zai(polygonsXY{pp}(:, 1), polygonsXY{pp}(:, 2)) ;
            c3d{pp} = [ xinterp, yinterp, zinterp ];
        end
    end
end

% Get fieldfaces for centroids
[cellCntrd3d, cellMeshFaces] = interpolate2Dpts_3Dmesh(faces, v2D, v3D, cntrds) ;
fN = faceNormal(triangulation(faces, v3D)) ;
cellNormals = nan(nCells, 3) ;
cellNormals(~isnan(cellMeshFaces), :) = fN(cellMeshFaces(~isnan(cellMeshFaces)), :) ; 
% Not needed: already normalized by construction
% cellNormals = cellNormals ./ vecnorm(cellNormals, 2, 2) ;

% Also get vectors point along zeta from cell centroid to lie along y
jac2d3d = jacobian2Dto3DMesh(v2D, v3D, faces) ;
% jac3d2d = jacobian3Dto2DMesh(v2D, v3D, faces) ;

%% Get aspect ratios of polygons 
% rotate to 2d plane, with x along zeta, and y along phi
areas = nan(nCells, 1) ;
perim = nan(nCells, 1) ;
moment1 = nan(nCells, 1) ;
ang1 = nan(nCells, 1) ;
moment2 = nan(nCells, 1) ;
ang2 = nan(nCells, 1) ;
moinertia = nan(nCells, 2, 2) ;
cellQ2d = {} ;

for cid = 1:nCells
    if ~isempty( c3d{cid} )
        cellVtx0 = c3d{cid} ;
        try
            cellVtx = cellVtx0 - cellCntrd3d(cid, :) ;
        catch
            error('what is the issue?')
        end
        if ~isempty(cellVtx)

            % Note: this approach is tricky since we have to map to
            % 3d and then back to 2d.
            % To figure out which direction to take to z, map vec to 3d

            % dzeta3d points towards the mapped z axis but in original 3d
            % embedding space
            dzeta3d = (jac2d3d{cellMeshFaces(cid)} * [0, 1]')' ;
            dzeta3d = dzeta3d / vecnorm(dzeta3d, 2, 2) ;

            % rotate cell to nearly yz plane
            rot2xy = rotate3dToAlignAxis(cellNormals(cid, :), dzeta3d) ;
            % Note: this one doesn't do the secondary rotation
            % rot = createRotationVector3d(cellNormals(cid, :), [0, 0, 1]) ;
            cell_quasi2d = (rot2xy * cellVtx')' ;

            % % Alternative method based on jacobian --> actually this 
            % % is a bit limited since it assumes small cells 
            % % (same size or smaller than the face) to get the 
            % % stretching in each dimension correct
            % jac = jac3d2d{cellMeshFaces(cid)} ;
            % dilation = jacobian2dilation(jac) ;
            % cell2d = (jac * cellVtx0')' ;
            % mapped centroid should be the same as pre-mapped centroid
            % up to rescaling in each dimension separately
            % cntrd2d = jac * cellCntrd(cid, :)' ;
            % cntrds(cid, :)


            % Check 2d cell polygon
            % if debug
            %     % dchi3d points towards the mapped x axis in embedding space
            %     dchi3d = (jac2d3d{cellMeshFaces(cid)} * [1, 0]')' ;
            %     dchi3d = dchi3d / vecnorm(dchi3d, 2, 2) ;
            % 
            %     % Plot the cell in 3d and 2d
            %     subplot(2, 2, 1)
            %     plot(cell2d0(:, 1), cell2d0(:, 2), '.-');
            %     axis equal
            %     subplot(2, 2, 2)
            % 
            %     % Vector to transform = dzeta since this emanates from
            %     % centroid just as the normal does.
            %     zeta2d = (rot2xy * dzeta3d')' ;
            %     try
            %         assert(all(abs(zeta2d - [0, 0, 1]) < 1e-7))
            %     catch
            %         error('Rotation was not successful!')
            %     end
            % 
            %     % plot it
            %     plot3(cell_quasi2d(:, 1), cell_quasi2d(:, 2), ...
            %         cell_quasi2d(:, 3), '.-'); 
            %     axis equal
            %     hold on; 
            % 
            %     subplot(2, 2, 3)
            %     plot3(cellVtx0(:, 1), cellVtx0(:, 2), cellVtx0(:, 3), ...
            %         '.-')
            %     hold on;
            %     zplus = cellCntrd(cid, :) + dzeta3d * mean(var(cellVtx0)); 
            %     xplus = cellCntrd(cid, :) + dchi3d * mean(var(cellVtx0));
            %     plot3dpts([cellCntrd(cid, :); zplus])
            %     plot3dpts([cellCntrd(cid, :); xplus])
            %     axis equal
            % end



            % Look at yz plane in which cell lives (very nearly, tangent plane)
            %   POLYGEOM( X, Y ) returns area, X centroid,
            %   Y centroid and perimeter for the planar polygon
            %   specified by vertices in vectors X and Y.
            %
            %   [ GEOM, INER, CPMO ] = POLYGEOM( X, Y ) returns
            %   area, centroid, perimeter and area moments of 
            %   inertia for the polygon.
            %   GEOM = [ area   X_cen  Y_cen  perimeter ]
            %   INER = [ Ixx    Iyy    Ixy    Iuu    Ivv    Iuv ]
            %     u,v are centroidal axes parallel to x,y axes.
            %   CPMO = [ I1     ang1   I2     ang2   J ]
            %     I1,I2 are centroidal principal moments about axes
            %         at angles ang1,ang2.
            %     ang1 and ang2 are in radians.
            %     J is centroidal polar moment.  J = I1 + I2 = Iuu + Ivv

            % Discard 3d info (non-tangent plane info) and 
            % compute polygon geometry
            if size(cellVtx0, 1) > 2
                try
                    [ geom, iner, cpmo ] = polygeom( cell_quasi2d(:, 2), ...
                        cell_quasi2d(:, 3) ) ;
                catch
                    error('here')
                end
                areas(cid) = geom(1) ;
                perim(cid) = geom(4) ;
                cellQ2d{cid} = cell_quasi2d ;
                moment1(cid) = cpmo(1) ;
                ang1(cid) = cpmo(2) ;
                moment2(cid) = cpmo(3) ;
                ang2(cid) = cpmo(4) ;
                moinertia(cid, :, :) = [iner(4) -iner(6); -iner(6) iner(5)] ;

            end
        else
            disp(['bad cell: ' num2str(cid)])
        end
    end
end