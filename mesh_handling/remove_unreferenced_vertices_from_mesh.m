function [ face, vertex, unreferenced, oldVertexIDx ] = ...
    remove_unreferenced_vertices_from_mesh( face, vertex )
%REMOVE_UNREFERENCED_VERTICES_FROM_MESH 
%   Removes a subset of the vertices which have no associated faces from 
%   a mesh triangulation.
%
%   INPUT PARAMETERS:
%       - face:         #Fx3 face connectivity list
%       - vertex:       #VxD vertex coordinate list
%
%   OUTPUT PARAMETERS:
%       - newFace:      #F'x3 updated face connectivity list
%       - newVertex:    #V'xD updated vertex coordinate list
%       - unreferenced: #N x1 list of vertex IDs to remove
%       - oldVertexIDx: #V'x1 list of the vertex IDs of the updated
%                       vertices in the old vertex list, so that, for ex,
%                       vn_new = old_mesh.vn(oldVertexIdx, :)
%
% NPMitchell 2020

unreferenced = [] ;
oldVertexIDx = 1:size(vertex, 1) ;
if ~isempty(setdiff(1:length(vertex), face(:)))
    unreferenced = setdiff(1:size(vertex, 1), face(:)) ;
    [ face, vertex, oldVertexIDx] = remove_vertex_from_mesh( fv.f, fv.v, unreferenced ) ;
end
        