function generateCurrentPullbacks(QS, cutMesh, spcutMesh, pbOptions)
%
%
% Parameters
% ----------
% pbOptions : struct
%   overwrite : bool (default=false)
%       overwrite existing images on disk
%   generate_uphi_coord : bool (default=false)
%       create pullbacks in uphi coords
%   PSize : int (default=5)
%       how many interpolation points along each side of a patch
%   EdgeColor : colorspec (default='none')
%       color of mesh edges in pullback image
%   YLim : 1x2 float (default=[0 1])
%       Y extent of pullback image in v or phi coords
%   Additional options as fields are
%       - Options.imSize:       The size of the output image
%       - Options.baseSize:     The side length in pixels of the smallest
%                               side of the output image when using the
%                               tight mesh bounding box
%       - Options.xLim:         The x-bounds of the output image
%       - Options.yLim:         The y-bounds of the output image
%       - Options.pixelSearch:  The method of searching for the faces
%                               containing pixel centers
%                                   - 'AABB' (requires GPToolBox)
%                                   - 'Default' (MATLAB built-ins, faster than AABB)
%       - Options.numLayers:    The number of onion layers to create
%                               Format is [ (num +), (num -) ]
%       - Options.layerSpacing: The spacing between adjacent onion layers
%                               in units of pixels
%       - Options.smoothIter:   Number of iterations of Laplacian mesh
%                               smoothing to run on the mesh prior to
%                               vertex normal displacement (requires
%                               GPToolBox) (Default is 0)
%       - Options.vertexNormal: User supplied vertex unit normals to the
%                               texture triangulation
%       - Options.Interpolant:  A pre-made texture image volume interpolant
%
% NPMitchell 2020

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

% Unpack pbOptions
overwrite = false ;
generate_uphi = false ;
generate_uv = true ;
generate_relaxed = true ;
if nargin > 3
    disp('Unpacking options')
    if isfield(pbOptions, 'overwrite')
        overwrite = pbOptions.overwrite ;
        pbOptions = rmfield(pbOptions, 'overwrite') ;
    end
    if isfield(pbOptions, 'generate_uv')
        generate_uv = pbOptions.generate_uv ;
        pbOptions = rmfield(pbOptions, 'generate_uv') ;
    end
    if isfield(pbOptions, 'generate_uphi_coord')
        generate_uphi = pbOptions.generate_uphi_coord ;
        pbOptions = rmfield(pbOptions, 'generate_uphi_coord') ;
    end
    if isfield(pbOptions, 'generate_relaxed')
        generate_relaxed = pbOptions.generate_relaxed ;
        pbOptions = rmfield(pbOptions, 'generate_relaxed') ;
    end
end

% Unpack QS
tt = QS.currentTime ;
a_fixed = QS.a_fixed ;
fileNameBase = QS.fileBase.name ;
imFolder = QS.dir.im ;
imFolder_r = QS.dir.im_r ;
imFolder_sp = QS.dir.im_sp ;
imFolder_up = QS.dir.im_up ;
axisorder = QS.data.axisOrder ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Create pullback using S,Phi coords \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------
% Generate Output Image Files
%--------------------------------------------------------------
imfn = sprintf( fullfile([imFolder, '/', fileNameBase, '_pb.tif']), tt); 
imfn_r = sprintf( fullfile([imFolder_r, '/', fileNameBase, '_pr.tif']), tt) ;
imfn_sp = sprintf( fullfile([imFolder_sp, '/', fileNameBase, '_pbsp.tif']), tt) ;
imfn_up = sprintf( fullfile([imFolder_up, '/', fileNameBase, '_pbup.tif']), tt) ;
pullbacks_exist1 = exist(imfn, 'file') && exist(imfn_r, 'file') ;
pullbacks_exist2 = exist(imfn_sp, 'file') && (exist(imfn_up, 'file') || ~generate_uphi) ;

do_pullbacks = (~pullbacks_exist1 || ~pullbacks_exist2 || overwrite) ;

if do_pullbacks
    % Load 3D data for coloring mesh pullback
    QS.getCurrentData()
    % grab raw stack data
    IV = QS.currentData.IV ;
end

if ~exist(imfn_sp, 'file') || overwrite
    fprintf(['Generating SP output image: ' imfn_sp]);
    % Assigning field spcutMesh.u to be [s, phi] (ringpath
    % and azimuthal angle)
    spcutMesh.u = spcutMesh.sphi ;
    aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_sp,...
        pbOptions, axisorder)
    spcutMesh = rmfield(spcutMesh, 'u') ;
end

if (~exist(imfn_up, 'file') || overwrite) && generate_uphi
    fprintf(['Generating uphi output image: ' imfn_up]);
    % Assigning field spcutMesh.u to be [s, phi] (ringpath
    % and azimuthal angle)
    spcutMesh.u = spcutMesh.uphi ;
    aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_up, ...
        pbOptions, axisorder)
    spcutMesh = rmfield(spcutMesh, 'u') ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Output Image File -- regular UV coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist(imfn, 'file') || overwrite) && generate_uv
    % Generate output image in uv
    fprintf(['Generating output image: ' imfn]);
    aux_generate_orbifold(cutMesh, a_fixed, IV, imfn, ...
        pbOptions, axisorder)
else
    disp('Skipping pullback image generation since exists')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save relaxed image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if (~exist(imfn_r, 'file') || overwrite) && generate_relaxed
    disp('Generating relaxed image for sphi coords...')
    spcutMesh.u = spcutMesh.sphi ;
    aux_generate_orbifold(spcutMesh, spcutMesh.ar, IV, imfn_r, ...
        pbOptions, axisorder)
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