% wtmi boston June 2011
% commented by Noah Mitchell 2015
% 
% TRISm is a (#bonds x 2) array listing the bond pairs, uniquely in such a
% way that the particle with the smaller number index is in the first 
% column and the particle with the larger number index is in the second.
%
% TRI is a (#tris x 3) array with one row for each triangle in the
% tessalation. The elements of the rows are the inds of the particles in
% that tri, in ascending order. If particle N is in 6 tris, it will most
% likely appear in each column twice, twice as the lowest index, twice as the
% middle, and twice as the highest.
%



function [TRI] = TRISm2TRI(TRISm)


nenemat=TRISm2NL(TRISm,9);
TRI = [];

% for each row in NL, check if both elements of TRISm (checked over all rows) are neighbors, and
% if so, add lines to TRI with [ TRISm col 1, TRISm col 2, row #]
for kk=1:size(nenemat,1)
    
    idxttemp = and(ismember(TRISm(:,1),nenemat(kk,:)),ismember(TRISm(:,2),nenemat(kk,:)));
    TRIStemp = TRISm(idxttemp,:);
    
    
    TRI = cat(1,TRI,cat(2,TRIStemp,kk*ones(size(TRIStemp,1),1)));
    
    
end

TRI = sort(TRI,2);
TRI = unique(TRI,'rows');