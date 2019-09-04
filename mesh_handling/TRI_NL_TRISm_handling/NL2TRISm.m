% wtmi 2010 modified:
% - aug 2010 (wtmi)
%commented Noah Mitchell 2015 
%
% NL is a neighbor list (NL), in which, for example,
%if particle 1 is bonded to particle 2 3 5 and 9, then the first row of NL
%is
%2 3 5 9 0 0 0 0,
%where the trailing zero columns are buffers for the possibility that
%some particles have many neighbors (ex in disordered lattice).
%
% TRISm is a (#bonds x 2) array listing the bond pairs, uniquely in such a
% way that the particle with the smaller number index is in the first 
% column and the particle with the larger number index is in the second.


function [TRISm] = NL2TRISm(NL)
TRISm = [];
%x = max(unique(TRIs));
for kk=1:length(NL)
    %Concatenate to TRISm as many rows as there are cols in the kkth row 
    % of NL. Each row has the neighbors listed in NL for that particle.
    TRISm = cat(1,TRISm,[kk*ones(size(NL,2),1),NL(kk,:)']);
end

tempidx = find(TRISm(:,1).*TRISm(:,2));
tempidx = setdiff(1:size(TRISm,1),tempidx);
TRISm(tempidx,:) = [];
TRISm = sort(TRISm,2);
TRISm = unique(TRISm,'rows');


