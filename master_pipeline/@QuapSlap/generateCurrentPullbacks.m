function generateCurrentPullbacks(QS, cutMesh, spcutMesh)
%
%
%
%

% Unpack options
if nargin < 2 || isempty(cutMesh)
    if isempty(QS.currentMesh.cutMesh)
        QS.loadCurrentCutMesh()
    end
    cutMesh = QS.currentMesh.cutMesh ;
end

if nargin < 3 || isempty(spcutMesh)
    if isempty(QS.currentMesh.spcutMesh)
        QS.loadCurrentSPCutMesh()
    end
    spcutMesh = QS.currentMesh.spcutMesh ;
end

tt = QS.currentTime ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Create pullback using S,Phi coords \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------
% Generate Output Image Files
%--------------------------------------------------------------
imfn = sprintf( fullfile([imFolder, '/', fileNameBase, '.tif']), tt); 
imfn_r = sprintf( fullfile([imFolder_r, '/', fileNameBase, '.tif']), tt) ;
imfn_sp = sprintf( fullfile([imFolder_sp, '/', fileNameBase, '.tif']), tt) ;
imfn_up = sprintf( fullfile([imFolder_up, '/', fileNameBase, '.tif']), tt) ;
pullbacks_exist1 = exist(imfn, 'file') && exist(imfn_r, 'file') ;
pullbacks_exist2 = exist(imfn_sp, 'file') && (exist(imfn_up, 'file') || ~generate_uphi_coord) ;
if (~pullbacks_exist1 || ~pullbacks_exist2 || overwrite_pullbacks) && ~IVloaded
    % Load 3D data for coloring mesh pullback
    xp.loadTime(t);
    xp.rescaleStackToUnitAspect();

    % Raw stack data
    IV = xp.stack.image.apply();
    IV = imadjustn(IV{1});
end

if ~exist(imfn_sp, 'file') || overwrite_pullbacks
    fprintf(['Generating SP output image: ' imfn_sp]);
    % Assigning field spcutMesh.u to be [s, phi] (ringpath
    % and azimuthal angle)
    spcutMesh.u = spcutMesh.sphi ;
    aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_sp)
    spcutMesh = rmfield(spcutMesh, 'u') ;
end

if (~exist(imfn_up, 'file') || overwrite_pullbacks) && generate_uphi_coord
    fprintf(['Generating uphi output image: ' imfn_up]);
    % Assigning field spcutMesh.u to be [s, phi] (ringpath
    % and azimuthal angle)
    spcutMesh.u = spcutMesh.uphi ;
    aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_up)
    spcutMesh = rmfield(spcutMesh, 'u') ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Output Image File -- regular UV coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(imfn, 'file') || overwrite_pullbacks
    % Generate output image in uv
    fprintf(['Generating output image: ' imfn]);
    aux_generate_orbifold(cutMesh, a_fixed, IV, imfn)
else
    disp('Skipping pullback image generation since exists')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save relaxed image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if ~exist(imfn_r, 'file')
    disp('Generating relaxed image for sphi coords...')
    spcutMesh.u = spcutMesh.sphi ;
    aux_generate_orbifold(spcutMesh, spcutMesh.ar, IV, imfn_r)
    spcutMesh = rmfield(spcutMesh, 'u') ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save submesh array. Each cell element contains all the 
% submeshes for that TP, which in this case is just one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meshStack{tidx} = cutMesh ;
% if generate_sphi_coord
%     spmeshStack{tidx} = spcutMesh ;
% end
fprintf('Done\n');