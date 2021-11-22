function [expts, exptIDs] = unpackDataMap(dmap)
%
%
%
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