function wcell = buildExpCell(dirs)
% buildExpCell(dirs)
% build a cell array with all experiments in dirs for RNAi experiments
% 
% 
% NPMitchell 2020

wcell = {} ; 
for ii = 1:length(dirs)
    wdir = dirs{ii} ;
    % load fates for this experiment
    fatefn = fullfile(wdir, 'embryos/embryo_fate.txt') ;
    disp(['Opening ' fatefn])
    fid = fopen(fatefn, 'rt') ;
    jj = 1 ;
    while true
        thisline = fgetl(fid) ;
        if ~ischar(thisline)
            % this is the end of the txt file
            break
        elseif strcmp(thisline(1:8), 'EmbryoID')
            disp('read header')
        else
            % split the line
            sline = strsplit(thisline, ' ') ;

            % add to struct
            wcell{ii}.embryoID(jj) = uint8(str2double(sline{1})) ;
            wcell{ii}.stageInit{jj} = sline{2} ;
            wcell{ii}.stageFinal{jj} = sline{3} ;
            wcell{ii}.fateID(jj) = uint8(str2double(sline{4})) ;
            wcell{ii}.missingfolds(jj) = uint8(str2double(sline{5})) ;
            wcell{ii}.notes{jj} = sline{6} ;
            for qq = 7:length(sline)
                wcell{ii}.notes{jj} = [wcell{ii}.notes{jj} ' ' sline{qq}] ;
            end
            
            % prepare for next line
            jj = jj + 1 ;
        end
    end
    fclose(fid);
end