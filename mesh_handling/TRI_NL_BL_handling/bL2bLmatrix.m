function bLmatrix = bL2bLmatrix(NL,BL,bL)
% bLmatrix = bL2bLmatrix(NL,BL,bL)
% Convert rest bond length list into matrix matching NL with rest bond
% lengths in it for entries.
%
% Parameters
% ----------
% NL : #vertices x max(#NN) int
%   neighbor list
% BL : #bonds x 2 int
%   bond list
% bL : #bonds x 1 float 
%   rest bond lengths 
%
% Returns
% -------
% bLmatrix : #vertices x max(#NN) float
%   rest bond lengths matching NL format
%
% NPMitchell 2023
%
    bLmatrix = zeros(size(NL)) ;
    for id = 1:length(NL)
        neighbors = NL(id, :)  ;
        for nidx = 1:length(neighbors)
            neigh = neighbors(nidx) ;
            if neigh > 0
                bid = find(BL(:, 1) == id & BL(:, 2) == neigh) ;
                if isempty(bid)
                    bid = find(BL(:, 2) == id & BL(:, 1) == neigh) ;
                end
                bLmatrix(id, nidx) = bL(bid) ;
            end
        end
    end
end