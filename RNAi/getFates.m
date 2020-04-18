function fates = getFates(wcell)
%
%
% NPMitchell 2020

dmyk = 1 ;
for ii = 1:length(wcell)
    for jj = 1:length(wcell{ii}.stageInit)
        % get fate
        fates(dmyk) = uint8(wcell{ii}.fateID(jj)) ;
        dmyk = dmyk + 1 ;
    end
end