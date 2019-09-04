% wtmi boston june 2011
% commented Noah Mitchell 2015
%
% TRI is a (#tris x 3) array with one row for each triangle in the
% tessalation. The elements of the rows are the inds of the particles in
% that tri, in ascending order. If particle N is in 6 tris, it will most
% likely appear in each column twice, once as the lowest index, once as the
% middle, and once as the highest.
%
% TRISm is a (#bonds x 2) array listing the bond pairs, uniquely in such a
% way that the particle with the smaller number index is in the first 
% column and the particle with the larger number index is in the second.
%

function [TRISm] = TRI2TRISm(TRI)

TRISm = TRI(:,[1,2]);
TRISm = cat(1,TRISm,TRI(:,[1,3]));
TRISm = cat(1,TRISm,TRI(:,[2,3]));


TRISm = sort(TRISm,2);
TRISm = unique(TRISm,'rows');