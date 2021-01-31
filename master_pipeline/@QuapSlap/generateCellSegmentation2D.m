function generateCellSegmentation2D(QS, options)
% Unfinished code -- Lin working on it.
% ToDo: 
%
%
% tissueAnalysisSuite fields are different than QuapSlap's. 
% Tor tissueAnalysisSuite, we have:
% vdat: 
%     nverts : neighbor list for vertices
%     ncells : 
%     vertxcoord : column of data in which vertex lives
%     vertycoord : row of data in which each vertex lives
% cdat : cell data    
%     ncells : indices of neighboring cells
%     nverts : vertices that define the cell
%     centroid.coord(1) x position of cell centroid
%     centroid.coord(2) y position of cell centroid
% bdat : 
%     nverts
%     ncells
%     pix : linear indices of the pixels associated with that bond
%
% Note on exterior calculus objects
% d0 is an e x c matrix of exterior derivatives with +1 and -1s
% at the endpts of each bond
% d1 is a v x e matrix of exterior derivatives. Upstream is +1,
% downstream is -1 when moving counterclockwise around a
% tension plaquette.
% d0 and d1 are matrices that take derivatives 

%% Parameters
% timepoints to process
timePoints = QS.xp.fileMeta.timePoints ;
% how far in pixels is too far to connect two cell vertices
very_far = 150 ;
% which coordinate system to use for segmentation
coordSys = QS.currentSegmentation.coordSys ; 
% Toggle for iLastik version control -- zero for newer version
iLastikVersion = 0;

%% unpack options
if isfield(options, 'very_far') 
    very_far = options.very_far ;
end
if isfield(options, 'coordSys') 
    coordSys = options.coordSys ;
end
if isfield(options, 'iLastikVersion') 
    iLastikVersion = options.iLastikVersion ;
end

%% Load in h5 from ilastik.
if strcmpi(coordSys, 'spsme')
    Folder = [QS.dir.im_r_sme, '_pixelClassification'] ;
    if ~exist(Folder, 'dir')
        mkdir(Folder)
        error(['Populate ' Folder ' with pixelClassification on pullbacks with coordSys ' coordSys])
    end
    filebase = [QS.fileBase.im_sp_sme(1:end-4) '_Probabilities.h5'] ;
else
    error('Have not coded for this coordinate system yet. Do so here')
end

for tp = timePoints
    
    % Define path to this timePoint's hdf5 probabilities file
    h5fn = fullfile(Folder, sprintf(filebase, tp)) ;
    [ mem ] = load.ilastikh5Single( h5fn, iLastikVersion );

    %% Segment the membrane.
    L = seg.memWS(mem, 50, 0, 1, 3.5);
    % Set bond=0 and clear_border = 1
    [L, Struct] = seg.generate_structs(L, 0, 1, 0, very_far);
    % Bad cells are bubble cells, which is a segmentation that forked and
    % reconnected.
    % ToDo: Should we do this? Can we skip it or does that lead to issues?
    L = seg.removeBadCells(Struct, L);
    disp('done removing bad cells')
    % Now change label matrix after removing bad cells
    L = seg.relabelL(L);
    % Now also synchronize Struct after removing bad cells
    [L,Struct] = seg.generate_structs(L, 0, 0, 0, very_far);
    disp('done with segmentation')

    %% Prepare data structure for inverse (optional? Does this improve segmentation?)
    % % put a parameter in the cdat of Struct, a boolean of whether every vertex
    % % is 3-fold.
    % Struct = seg.threefold_cell(Struct);
    % % generate the bdat structure in Struct
    % Struct = seg.recordBonds(Struct, L);
    % disp('generated the bond structure')
    % % Segment the curvature of each bond
    % Struct = seg.curvature(Struct, size(L));
    % disp('segmented the curvature of each bond')
    % % Remove all fourfold vertices, recursively if there are z>4
    % Struct = seg.removeFourFold(Struct);
    % disp('removed fourfold vertices')
    % % The inverse is ill-posed if we have convex cells, so hack those to be
    % % convex
    % Struct = seg.makeConvexArray(Struct);
    % disp('done with data preparation')
    
    %% Consider converting segmentation to simpler struct
    seg2d.v2d = [] ;
    
    %% Save the segmentation to disk
    outfn = sprintf(QS.fullFileBase.segmentation2d, tp) ;
    save(outfn, seg2d)
end
