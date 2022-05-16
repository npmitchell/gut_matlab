function [expts, exptIDs] = unpackDataMap(dmap)
%[expts, exptIDs] = unpackDataMap(dmap)
% Utility function for a containers map
%
% Parameters
% ----------
% dmap : containers map with keys dmap.keys
%   the container map to unpack
% 
% Returns
% -------
% expts : Nx1 cell of dmap(key).folders for all keys
% exptIDs : Nx1 cell of dmap(key).ids for all keys
%

dmyk = 1;
keys = dmap.keys ;
expts = {} ;
exptIDs = {} ;
for keyID =  1:length(keys)
    key = keys{keyID} ;
    for qq = 1:length(dmap(key).folders)
        expts{dmyk} = dmap(key).folders{qq} ;
        exptIDs{dmyk} = dmap(key).ids(qq) ;
        dmyk = dmyk + 1;
    end
end