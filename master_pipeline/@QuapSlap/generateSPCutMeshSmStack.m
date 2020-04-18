function generateSPCutMeshSmStack(QS, spcutMeshSmStackOptions)
% generateSPCutMeshSmStack(QS, spcutMeshSmStackOptions)
%
%
% Parameters
% ----------
% QS : QuapSlap class instance
% spcutMeshSmStackOptions: struct with fields
%   - n_outward : int 
%       number of steps in positive normal direction to sample
%   - n_inward : int
%       number of steps in negative normal direction to sample
%   - overwrite : bool
%       whether to overwrite spcutMeshSmStack on disk
%
% Returns 
% -------
%
% NPMitchell 2020

% Unpack options
n_outward = 10 ;
n_inward = 10 ;
overwrite = false ;
if isfield(spcutMeshSmStackOptions, 'overwrite')
    overwrite = spcutMeshSmStackOptions.overwrite ;
end
if isfield(spcutMeshSmStackOptions, 'n_outward')
    n_outward = spcutMeshSmStackOptions.n_outward ;
end
if isfield(spcutMeshSmStackOptions, 'n_inward')
    n_inward = spcutMeshSmStackOptions.n_inward ;
end


% Unpack QS
spcutMeshSmBase = QS.fullFileBase.spcutMeshSm ;
fileNameBase = QS.fileBase.name ;

for qq = 1:length(QS.xp.fileMeta.timePoints)
    tt = QS.xp.fileMeta.timePoints(qq) ;
    disp(['t = ' num2str(tt)])
    QS.setTime(tt)
    
    % Load time-smoothed mesh
    load(sprintf(spcutMeshSmBase, tt), 'spcutMeshSm') ;
    
    %--------------------------------------------------------------
    % Generate Output Image File
    %--------------------------------------------------------------
    spacingstr = strrep(sprintf('%0.2fum', layer_spacing * resolution), '.', 'p') ;
    imfn_spsm = sprintf( fullfile( QS.dir.im_spsm_e2, ...
        [fileNameBase, '_%02d_%02d_' spacingstr '.tif']), tt, n_outward, n_inward ) ;
    
    if ~exist(imfn_spsm, 'file') || overwrite
        % Load 3D data for coloring mesh pullback
        xp.getCurrentData(QS)
        IV = adjustIV();
        
        fprintf(['Generating SP output image for sm mesh: ' imfn_spsm]);
        % Assigning field spcutMesh.u to be [s, phi] (ringpath
        % and azimuthal angle)
        Options.numLayers = [n_outward, n_inward] ;
        Options.layerSpacing = layer_spacing ;
        Options.smoothIter = 1 ;
        Options.yLim = [-0.5, 1.5] ;
        % Note that we pass a_fixed * 0.5 since the image is extended by a
        % factor of two
        aux_generate_orbifold( spcutMeshSm, a_fixed * 0.5, IV, imfn_spsm, Options)
    end
    clear Options
end
