function aux_generate_orbifold(cutMesh, a, IV, imfn)

% Texture patch options
Options.PSize = 5;
Options.EdgeColor = 'none';

% Generate Tiled Orbifold Triangulation ------------------------------
tileCount = [1 1];  % how many above, how many below
[ TF, TV2D, TV3D ] = tileAnnularCutMesh( cutMesh, tileCount );

% View Results -------------------------------------------------------
% patch( 'Faces', TF, 'Vertices', TV2D, 'FaceVertexCData', ...
%     TV3D(:,3), 'FaceColor', 'interp', 'EdgeColor', 'k' );
% axis equal

% Texture image options
Options.imSize = ceil( 1000 .* [ 1 a ] );
Options.yLim = [0 1];

% profile on
% Create texture image
if any(isnan(TV2D))
    error('here -- check for NaNs case')
end
patchIm = texture_patch_to_image( TF, TV2D, TF, TV3D(:, [2 1 3]), ...
    IV, Options );
% profile viewer

fprintf('Done\n');

% View results --------------------------------------------
% imshow( patchIm );
% set( gca, 'YDir', 'Normal' );

% Write figure to file
disp(['Writing ' imfn]) 
imwrite( patchIm, imfn, 'TIFF' );            

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Save extended relaxed image
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% disp('Generating relaxed, extended image...')          
% % Format axes
% xlim([0 ar]); ylim([-0.5 1.5]);
% 
% % Extract image from figure axes
% patchIm_e = getframe(gca);
% patchIm_e = rgb2gray(patchIm_e.cdata);
% 
% % Write figure to file
% imwrite( patchIm_e, ...
%     sprintf( fullfile([imFolder_re, '/', fileNameBase, '.tif']), t ), ...
%     'TIFF' );

% Close open figures
close all