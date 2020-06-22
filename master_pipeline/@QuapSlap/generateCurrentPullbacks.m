function generateCurrentPullbacks(QS, cutMesh, spcutMesh, ...
    spcutMeshSm, pbOptions)
%GENERATECURRENTPULLBACKS(QS, cutMesh, spcutMesh, spcutMeshSm, pbOptions)
%   Generate 2D images of tissue mapped from 3D
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

%% Unpack pbOptions
overwrite = false ;         % generate pullbacks even if they exist on disk
generate_sphi    = true  ;  % generate an (s, phi) coord system pullback
generate_relaxed = false ;  % generate a relaxed (s,phi) coord system pullback
generate_uv      = false ;  % generate a (u,v) coord system pullback
generate_uphi    = false ;  % generate a (u, phi) coord system pullback
generate_spsm    = false ;  % generate an (s, phi) coord system smoothed mesh pullback
generate_rsm     = false ;  % generate a relaxed (s, phi) coord system smoothed mesh pullback
% Other options
save_as_stack    = false ;

% Replace defaults
if nargin > 4
    disp('Unpacking options for which pullbacks to generate/overwrite')
    if isfield(pbOptions, 'overwrite')
        overwrite = pbOptions.overwrite ;
        pbOptions = rmfield(pbOptions, 'overwrite') ;
    end
    if isfield(pbOptions, 'generate_sphi')
        generate_sphi = pbOptions.generate_sphi ;
        pbOptions = rmfield(pbOptions, 'generate_sphi') ;
    end
    if isfield(pbOptions, 'generate_relaxed')
        generate_relaxed = pbOptions.generate_relaxed ;
        pbOptions = rmfield(pbOptions, 'generate_relaxed') ;
    end
    if isfield(pbOptions, 'generate_uv')
        generate_uv = pbOptions.generate_uv ;
        pbOptions = rmfield(pbOptions, 'generate_uv') ;
    end
    if isfield(pbOptions, 'generate_uphi_coord')
        generate_uphi = pbOptions.generate_uphi_coord ;
        pbOptions = rmfield(pbOptions, 'generate_uphi_coord') ;
    end
    if isfield(pbOptions, 'generate_spsm')
        generate_spsm = pbOptions.generate_spsm ;
        pbOptions = rmfield(pbOptions, 'generate_spsm') ;
    end
    if isfield(pbOptions, 'generate_rsm')
        generate_rsm = pbOptions.generate_rsm ;
        pbOptions = rmfield(pbOptions, 'generate_rsm') ;
    end
end

%% Unpack options
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

if nargin < 4 || isempty(spcutMeshSm)
    if isempty(QS.currentMesh.spcutMeshSm)
        QS.loadCurrentSPCutMeshSm()
    end
    spcutMeshSm = QS.currentMesh.spcutMeshSm ;
end

%% Unpack QS
tt = QS.currentTime ;
a_fixed = QS.a_fixed ;
% fileNameBase = QS.fileBase.name ;
% imFolder = QS.dir.im ;
% imFolder_r = QS.dir.im_r ;
% imFolder_sp = QS.dir.im_sp ;
% imFolder_up = QS.dir.im_up ;
% imFolder_spsm = QS.dir.im_sp_sm ;
% imFolder_rsm = QS.dir.im_r_sm ;
axisorder = QS.data.axisOrder ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Checking whether to create pullback \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------
% Generate Output Image Files
%--------------------------------------------------------------
imfn_uv = sprintf( QS.fullFileBase.im_uv, tt); 
imfn_r = sprintf( QS.fullFileBase.im_r, tt) ;
imfn_sp = sprintf( QS.fullFileBase.im_sp, tt) ;
imfn_up = sprintf( QS.fullFileBase.im_up, tt) ;
imfn_spsm = sprintf( QS.fullFileBase.im_sp_sm, tt) ;
imfn_rsm = sprintf( QS.fullFileBase.im_r_sm, tt) ;
do_pb1 = ~exist(imfn_uv, 'file') && generate_uv ;
do_pb2 = ~exist(imfn_r, 'file') && generate_relaxed ;
do_pb3 = ~exist(imfn_sp, 'file') && generate_sphi ;
do_pb4 = ~exist(imfn_up, 'file') && generate_uphi ;
do_pb5 = ~exist(imfn_spsm, 'file') && generate_spsm ;
do_pb6 = ~exist(imfn_rsm, 'file') && generate_rsm ;

do_pb = [do_pb1, do_pb2, do_pb3, do_pb4, do_pb5, do_pb6] ;
do_pullbacks = (any(do_pb) || overwrite) ;

if do_pullbacks
    % Declare what needs to be redone
    if overwrite
        disp('All pullback images will be recomputed & saved')
    else
        if do_pb1
            disp(['(u,v) PB will be generated: ', imfn_uv])
        end 
        if do_pb2
            disp(['Relaxed (s,phi) PB will be generated: ', imfn_r])
        end
        if do_pb3
            disp(['(s,phi) PB will be generated: ', imfn_sp])
        end
        if do_pb4
            disp(['(u,phi) PB will be generated: ', imfn_up])
        end
        if do_pb5
            disp(['Smooth (s,phi) PB will be generated: ', imfn_spsm])
        end
        if do_pb6
            disp(['Smooth relaxed (s,phi) PB will be generated: ', imfn_rsm])
        end
    end     
    
    % Load 3D data for coloring mesh pullback
    QS.getCurrentData()
    % grab raw stack data
    IV = QS.currentData.IV ;
end

if (~exist(imfn_sp, 'file') || overwrite) && generate_sphi
    fprintf(['Generating SP output image: ' imfn_sp]);
    % Assigning field spcutMesh.u to be [s, phi] (ringpath
    % and azimuthal angle)
    spcutMesh.u = spcutMesh.sphi ;
    aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_sp,...
        pbOptions, axisorder, save_as_stack)
    spcutMesh = rmfield(spcutMesh, 'u') ;    
else
    disp('Skipping SP pullback image generation ')
end

if (~exist(imfn_up, 'file') || overwrite) && generate_uphi
    fprintf(['Generating uphi output image: ' imfn_up]);
    % Assigning field spcutMesh.u to be [s, phi] (ringpath
    % and azimuthal angle)
    spcutMesh.u = spcutMesh.uphi ;
    aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_up, ...
        pbOptions, axisorder, save_as_stack)
    spcutMesh = rmfield(spcutMesh, 'u') ;
else
    disp('Skipping UP pullback image generation ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Output Image File -- regular UV coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist(imfn_uv, 'file') || overwrite) && generate_uv
    % Generate output image in uv
    fprintf(['Generating UV output image: ' imfn_uv]);
    aux_generate_orbifold(cutMesh, a_fixed, IV, imfn_uv, ...
        pbOptions, axisorder, save_as_stack)
else
    disp('Skipping UV pullback image generation ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save relaxed image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if (~exist(imfn_r, 'file') || overwrite) && generate_relaxed
    disp('Generating relaxed image for sphi coords...')
    spcutMesh.u = spcutMesh.sphi ;
    aux_generate_orbifold(spcutMesh, spcutMesh.ar, IV, imfn_r, ...
        pbOptions, axisorder, save_as_stack)
    spcutMesh = rmfield(spcutMesh, 'u') ;
else
    disp('Skipping relaxed SP pullback image generation ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save smoothed sp image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if (~exist(imfn_spsm, 'file') || overwrite) && generate_spsm
    disp('Generating image for smoothed sphi coords...')
    aux_generate_orbifold(spcutMeshSm, a_fixed, IV, imfn_spsm, ...
        pbOptions, axisorder, save_as_stack)
else
    disp('Skipping SPSm pullback image generation ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save smoothed relaxed image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if (~exist(imfn_rsm, 'file') || overwrite) && generate_rsm
    disp('Generating relaxed image for sphi coords...')
    if ~isfield(spcutMeshSm, 'ar')
        % Compute relaxed aspect ratio
        tmp = spcutMeshSm.u ;
        tmp(:, 1) = tmp(:, 1) / max(tmp(:, 1)) ;
        arspsm = minimizeIsoarealAffineEnergy( spcutMeshSm.f, spcutMeshSm.v, tmp );
        spcutMeshSm.ar = arspsm ;
    end
    aux_generate_orbifold(spcutMeshSm, spcutMeshSm.ar, IV, imfn_rsm, ...
        pbOptions, axisorder, save_as_stack)
else
    disp('Skipping relaxed SPSm pullback image generation ')
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