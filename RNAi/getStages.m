function stages = getStages(wcell)
%
%
% NPMitchell 2020

dmyk = 1 ;
stages = [0] ;
for ii = 1:length(wcell)
    for jj = 1:length(wcell{ii}.stageInit)
        disp([num2str(ii) ': ' num2str(jj) '/' num2str(length(wcell{ii}.stageInit))])
        % get stage
        if ~isnan(str2double(wcell{ii}.stageInit{jj}))
            stages(dmyk) = uint8(str2double(wcell{ii}.stageInit{jj})) ;
        else
            if strcmp(wcell{ii}.stageInit{jj}, '?')
                stages(dmyk) = 0 ;
            elseif ~isnan(str2double(wcell{ii}.stageInit{jj}(1:2)))
                stages(dmyk) = uint8(str2double(wcell{ii}.stageInit{jj}(1:2))) ;
            else
                stages(dmyk) = uint8(str2double(wcell{ii}.stageInit{jj}(1))) ;
            end
        end
        dmyk = dmyk + 1 ;
                
    end
end